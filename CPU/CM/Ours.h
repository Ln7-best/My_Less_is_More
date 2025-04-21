#ifndef OURS_H
#define OURS_H
#include "heap.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include "readerwriterqueue.h"
#include "config.h"
#include "CM.h"
#include <barrier>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <ratio>
#include <vector>

#define MEASURETIME
// #include "Abstract.h"
const uint64_t PROCESSGAP = 100000;
const uint64_t COUTERGAP = PROCESSGAP / 100;

struct alignas(64) global_sketch_sub_section
{
  std::atomic<uint64_t> counter[HASH_NUM * sub_sketch_length];
  char padding[64 - (sizeof(counter) % 64)];
};

template <typename Key>
struct HH_Bucket
{
  uint64_t vote;
  Key ID[COUNTER_PER_BUCKET];
  uint64_t count[COUNTER_PER_BUCKET];
};

template <typename Key>
struct Stash_Bucket
{
  uint64_t vote;
  Key ID[COUNTER_PER_BUCKET_FILTER];
  uint64_t count[COUNTER_PER_BUCKET_FILTER];
  uint16_t pos[COUNTER_PER_BUCKET_FILTER];
};

template <typename Key>
struct alignas(64) global_buckets_sub_section
{
  HH_Bucket<Key> buckets[BUCKET_LENGTH];
  char padding[64 - (sizeof(buckets) % 64)];
};

template <typename Key>
struct alignas(64) child_buckets_sub_section
{
  Stash_Bucket<Key> buckets[HASH_NUM][FILTER_BUCKET_LENGTH];
  char padding[64 - sizeof(buckets) % 64];
};
template <typename Key, uint32_t thread_num>
class Ours
{
public:
  void update(void *start, uint64_t size, HashMap mp,
              double *throughput = nullptr)
  {
    size = size;
    std::thread parent;
    parent = std::thread(&Ours::ParentThread, this, &parent, start, size, &mp,
                         throughput);
    parent.join();
  }

private:
  struct global_sketch_sub_section global_sketch[thread_num];
  double running_time[thread_num];
  uint32_t dataset_size[thread_num];
  std::atomic<uint32_t> partition_num;
  Paddedatomic process_counter[thread_num];
  Paddedint operation_cnt[thread_num];
  global_buckets_sub_section<Key> global_buckets[thread_num];
  child_buckets_sub_section<Key> child_filters[thread_num];
  typedef ReaderWriterQueue<CM_Entry<Key>> myQueue;
  struct Paddedint global_cnt[thread_num];
  struct Paddedint round[thread_num];
  myQueue que[thread_num][thread_num];
  MyChild_CM<Key> *initialize_child() { return new MyChild_CM<Key>(); }

  static void HHCompare(HashMap test, HashMap real, int32_t threshold, std::ofstream *outputFile = nullptr)
  {
    double estHH = 0, HH = 0, both = 0;
    double CR = 0, PR = 0, AAE = 0, ARE = 0;
    std::cout << "real size" << real.size() << std::endl;
    std::cout << "test size" << test.size() << std::endl;
    uint64_t cnt = 0;
    for (auto it = test.begin(); it != test.end(); ++it)
    {
      if (it->second > threshold)
      {
        estHH += 1;
        if (real[it->first] > threshold)
        {
          both += 1;
          AAE += abs(real[it->first] - it->second);
          ARE += abs(real[it->first] - it->second) / (double)real[it->first];
        }
      }
      else if ((int64_t)threshold - (int64_t)it->second < 10000)
      {
        cnt++;
      }
    }
    uint64_t test_hit = 0;
    for (auto it = real.begin(); it != real.end(); ++it)
    {
      if (it->second > threshold)
      {
        HH += 1;
        if (test.find(it->first) != test.end())
        {
          std::cout << (int64_t)threshold - (int64_t)test[it->first] << std::endl;
          test_hit++;
        }
      }
    }
    if (!outputFile)
      std::cout << "real HH:" << HH << std::endl
                << "correct HH: " << both << std::endl
                << "test hit: " << test_hit << std::endl
                << "num: " << cnt << std::endl
                << "CR: " << both / HH << std::endl
                << "PR: " << both / estHH << std::endl
                << "AAE: " << AAE / both << std::endl
                << "ARE: " << ARE / both << std::endl
                << std::endl;
    else
      *outputFile << "CR: " << std::to_string(both / HH) << std::endl
                  << "PR: " << std::to_string(both / estHH) << std::endl
                  << "AAE: " << std::to_string(AAE / both) << std::endl
                  << "ARE: " << std::to_string(ARE / both) << std::endl;
  }

  static void Partition(Key *start, uint64_t size, uint32_t id, std::vector<Key> &vec)
  {
    uint64_t interval = size / thread_num;
    for (uint64_t i = interval * id; i < interval * (id + 1); i++)
    {
      vec.push_back(start[i]);
    }
  }
  HashMap query_all()
  {
    HashMap ret;
    for (uint32_t i = 0; i < thread_num; i++)
    {
      for (uint32_t j = 0; j < BUCKET_LENGTH; j++)
      {
        for (uint32_t k = 0; k < COUNTER_PER_BUCKET; k++)
        {
          if (global_buckets[i].buckets[j].count[k] == 0)
            continue;
          ret[global_buckets[i].buckets[j].ID[k]] = global_buckets[i].buckets[j].count[k];
        }
      }
    }
    return ret;
  }

  void ParentThread(std::thread *thisThd, void *start, uint64_t size,
                    HashMap *mp, double *throughput)
  {
#ifdef __linux__
    if (!setaffinity(thisThd, thread_num))
      return;
#endif
    std::atomic<int32_t> finish(0);
    std::thread thd[thread_num];
    partition_num.store(0);
    for (uint64_t i = 0; i < thread_num; i++)
    {
      for (uint64_t j = 0; j < HASH_NUM * sub_sketch_length; j++)
      {
        global_sketch[i].counter[j] = 0;
      }
    }
    for (size_t t = 0; t < thread_num; ++t)
    {
      for (size_t h = 0; h < HASH_NUM; ++h)
      {
        for (size_t b = 0; b < FILTER_BUCKET_LENGTH; ++b)
        {
          child_filters[t].buckets[h][b].vote = 0;
          for (uint64_t i = 0; i < COUNTER_PER_BUCKET_FILTER; i++)
          {
            child_filters[t].buckets[h][b].ID[i] = 0;
            child_filters[t].buckets[h][b].count[i] = 0;
            child_filters[t].buckets[h][b].pos[i] = 0;
          }
        }
      }
    }
    for (uint64_t i = 0; i < thread_num; i++)
    {
      for (uint64_t j = 0; j < BUCKET_LENGTH; j++)
      {
        global_buckets[i].buckets[j].vote = 0;
        for (uint64_t k = 0; k < COUNTER_PER_BUCKET; k++)
        {
          global_buckets[i].buckets[j].ID[k] = 0;
          global_buckets[i].buckets[j].count[k] = 0;
        }
      }
    }

    for (uint32_t i = 0; i < thread_num; ++i)
    {
      operation_cnt[i].value = 0;
      global_cnt[i].value = 0;
      round[i].value = 0;
      thd[i] = std::thread(&Ours::ChildThread, this, &(thd[i]), i, start, size,
                           &finish, global_sketch + i);
    }
    while (partition_num < thread_num)
    {
    }
    for (uint32_t i = 0; i < thread_num; ++i)
    {
      thd[i].join();
    }
    // collect(heap, finish);
    while (finish.load(std::memory_order_seq_cst) < thread_num)
    {
    }
    double avg_time = 0;
    double max_time = 0;
    double avg_numa1 = 0;
    uint64_t tot_size = 0;
    double avg_numa2 = 0;
    for (uint32_t i = 0; i < thread_num; i++)
    {
      avg_time += running_time[i];
      max_time = std::max(max_time, running_time[i]);
      tot_size += dataset_size[i];
      std::cout << "thread id: " << i << " dataset size: " << dataset_size[i]
                << " running time: " << running_time[i] << std::endl;
    }
    avg_time /= thread_num;
    *throughput = tot_size / max_time * 1000;
    double avg_throughput = tot_size / avg_time * 1000;
    std::cout << "tot_process" << tot_size << std::endl;
    std::cout << "avg: " << avg_time << std::endl;
    std::cout << "max: " << max_time << std::endl;
    std::cout << "avg throughput: " << avg_throughput << std::endl;
    std::cout << "throughput: "
              << *throughput << std::endl;
    HashMap ret = query_all();
    HHCompare(ret, (*mp), size / sizeof(Key) * ALPHA);
  }

  /**
   * The thread of each worker
   */
  void ChildThread(std::thread *thisThd, uint32_t thread_id, void *start,
                   uint64_t size, std::atomic<int32_t> *finish,
                   struct global_sketch_sub_section *sketch_sub_section)
  {
#ifdef __linux__
    if (!setaffinity(thisThd, thread_id))
      return;
#endif
    MyChild_CM<Key> *sketch = initialize_child();
    std::vector<Key> dataset;
    // uint64_t round = 0;
    Partition((Key *)start, size / sizeof(Key), thread_id,
              dataset);
    partition_num.fetch_add(1);
    while (partition_num < thread_num)
    {
    }
#ifdef MEASURETIME
    auto start_time = std::chrono::high_resolution_clock::now();
#endif
    insert(dataset, sketch, que, thread_id, sketch_sub_section);

#ifdef MEASURETIME
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end_time - start_time;
    running_time[thread_id] = duration.count();
#endif
    dataset_size[thread_id] = dataset.size();
    (*finish).fetch_add(1, std::memory_order_seq_cst);
    while ((*finish).load(std::memory_order_seq_cst) < thread_num)
    {
      process_queue(que, global_sketch, thread_id);
    }
    delete sketch;
  }

  void insert(const std::vector<Key> &dataset, MyChild_CM<Key> *sketch,
              myQueue queue_group[][thread_num], uint32_t thread_id,
              struct global_sketch_sub_section *sketch_sub_section)
  {
    uint32_t length = dataset.size();
    for (uint32_t i = 0; i < length; ++i)
    {
      insert_child(sketch, queue_group, dataset[i], thread_id, global_sketch);
      if (i % COUTERGAP == 0)
        process_counter[thread_id].value.store(i + 1, std::memory_order_relaxed);
      if (i % PROCESSGAP == 0)
        process_queue(queue_group, global_sketch, thread_id);
    }
    process_counter[thread_id].value.store(length);
  }

  void insert_child(MyChild_CM<Key> *p, myQueue q_group[][thread_num],
                    const Key &key,
                    uint64_t thread_id,
                    struct global_sketch_sub_section *sketch_sub_section)
  {
    double hit_rate;
    uint64_t local_min;
    auto sketch = p->sketch;
    uint32_t pos[HASH_NUM];
    for (uint32_t hashPos = 0; hashPos < HASH_NUM; ++hashPos)
    {
      pos[hashPos] = hash(key, hashPos) % LENGTH;
    }

    for (uint32_t hashPos = 0; hashPos < HASH_NUM; ++hashPos)
    {
      sketch[hashPos][pos[hashPos]] += 1;
      if (__builtin_expect(sketch[hashPos][pos[hashPos]] >= PROMASK, 0))
      {
        bool update_flag = false;
        uint64_t minVal = 0x7fffffff;
        uint32_t minPos = 0;
        // uint64_t temp_val = 1;
        sketch[hashPos][pos[hashPos]] = 0;
        InsertStash(q_group, key, thread_id, sketch_sub_section, hashPos, pos[hashPos]);
      }
    }
  }

  inline void InsertStash(myQueue q_group[][thread_num],
                          const Key &key,
                          uint64_t thread_id,
                          struct global_sketch_sub_section *sketch_sub_section, uint16_t hashPos, uint16_t pos)
  {
    bool update_flag = false;
    uint64_t minVal = 0x7fffffff;
    uint32_t minPos = 0;
    // uint64_t temp_val = 1;
    uint64_t idx = pos % FILTER_BUCKET_LENGTH;
    for (uint32_t j = 0; j < COUNTER_PER_BUCKET_FILTER; j++)
    {
      if (this->child_filters[thread_id].buckets[hashPos][idx].ID[j] ==
          key)
      {
        this->child_filters[thread_id].buckets[hashPos][idx].count[j]++;
        if (this->child_filters[thread_id].buckets[hashPos][idx].count[j] == 9)
        {
          uint64_t push_sec_id = pos / sub_sketch_length;
          uint16_t sec_inner_index = pos % sub_sketch_length;
          while (__builtin_expect(!q_group[push_sec_id][thread_id].enqueue_fast(
                                      CM_Entry<Key>(key, hashPos, sec_inner_index,
                                                    9 * PROMASK)),
                                  0))
          {
            process_queue(q_group, sketch_sub_section, thread_id);
          }
          this->child_filters[thread_id].buckets[hashPos][idx].count[j] = 0;
          this->child_filters[thread_id].buckets[hashPos][idx].ID[j] = 0;
          this->child_filters[thread_id].buckets[hashPos][idx].pos[j] = 0;
        }
        update_flag = true;
      }
      if (this->child_filters[thread_id].buckets[hashPos][idx].count[j] <
          minVal)
      {
        minPos = j;
        minVal = this->child_filters[thread_id].buckets[hashPos][idx].count[j];
      }
    }
    if (update_flag)
      return;
    uint64_t insert_value = 1;
    uint64_t insert_pos = pos;
    uint64_t insert_key = key;
    if ((this->child_filters[thread_id].buckets[hashPos][idx].vote + 1) >= minVal * 8 || insert_pos == this->child_filters[thread_id].buckets[hashPos][idx].pos[minPos])
    {
      insert_value = minVal;
      insert_key = this->child_filters[thread_id].buckets[hashPos][idx].ID[minPos];
      insert_pos = this->child_filters[thread_id].buckets[hashPos][idx].pos[minPos];
      this->child_filters[thread_id].buckets[hashPos][idx].vote = 0;
      this->child_filters[thread_id].buckets[hashPos][idx].ID[minPos] =
          key;
      this->child_filters[thread_id].buckets[hashPos][idx].count[minPos] =
          1;
      this->child_filters[thread_id].buckets[hashPos][idx].pos[minPos] =
          pos;
    }
    else
    {
      this->child_filters[thread_id].buckets[hashPos][idx].vote++;
    }
    if (insert_value == 0)
      return;
    uint64_t push_sec_id = insert_pos / sub_sketch_length;
    uint16_t sec_inner_index = insert_pos % sub_sketch_length;
    while (__builtin_expect(!q_group[push_sec_id][thread_id].enqueue_fast(CM_Entry<Key>(
                                insert_key, hashPos, sec_inner_index, insert_value * PROMASK)),
                            0))
    {
      process_queue(q_group, sketch_sub_section, thread_id);
    }
  }

  inline void process_queue(myQueue q_group[][thread_num],
                            struct global_sketch_sub_section *sketch_sub_section,
                            uint64_t thread_id = 0)
  {
    CM_Entry<Key> temp;
    uint64_t cnt = 0;
    for (uint64_t i = 0; i < thread_num; i++)
    {
      while (q_group[thread_id][i].try_dequeue(temp))
      {
        const uint16_t &hash_pos = temp.hashPos;
        const uint16_t &pos = temp.pos;
        uint64_t counter_pos = hash_pos * sub_sketch_length + pos;
        sketch_sub_section[thread_id].counter[counter_pos].fetch_add(
            temp.value);
        uint64_t counter_val =
            sketch_sub_section[thread_id].counter[counter_pos].load();
        if ((double)counter_val >
            this->global_cnt[thread_id].value * (ALPHA / 2))
        {
          if (hash_pos != 1)
            continue;
          if (cnt == 0 && this->round[thread_id].value % 9 == 0)
          {
            this->global_cnt[thread_id].value = 0;
            for (uint64_t k = 0; k < thread_num; k++)
              this->global_cnt[thread_id].value +=
                  this->process_counter[k].value.load();
          }
          int32_t minimum = counter_val;
          int32_t origin_min = minimum;
          for (uint32_t tempHash = 0; tempHash < HASH_NUM; ++tempHash)
          {
            uint32_t tempPos = hash(temp.key, tempHash) % LENGTH;
            uint64_t push_sec_id = tempPos / sub_sketch_length;
            uint16_t sec_inner_index = tempPos % sub_sketch_length;
            minimum =
                MIN(minimum,
                    sketch_sub_section[push_sec_id]
                        .counter[tempHash * sub_sketch_length + sec_inner_index]
                        .load());
          }
          int32_t minVal = 0x7fffffff;
          uint32_t minPos = 0;
          uint64_t idx = hash(temp.key, 101) % BUCKET_LENGTH;
          bool update_flag = 0;
          for (uint32_t j = 0; j < COUNTER_PER_BUCKET; j++)
          {
            if (this->global_buckets[thread_id].buckets[idx].ID[j] ==
                temp.key)
            {
              update_flag = 1;
              this->global_buckets[thread_id].buckets[idx].count[j] = minimum;
              break;
            }
            if (this->global_buckets[thread_id].buckets[idx].count[j] <
                minVal)
            {
              minPos = j;
              minVal = this->global_buckets[thread_id].buckets[idx].count[j];
            }
          }
          if (!update_flag)
          {
            if ((this->global_buckets[thread_id].buckets[idx].vote + minimum) >=
                minVal * 8)
            {
              this->global_buckets[thread_id].buckets[idx].vote = 0;
              this->global_buckets[thread_id].buckets[idx].ID[minPos] =
                  temp.key;
              this->global_buckets[thread_id].buckets[idx].count[minPos] =
                  minimum;
            }
            else
            {
              this->global_buckets[thread_id].buckets[idx].vote += minimum;
            }
          }
          cnt++;
        }
      }
    }
    this->round[thread_id].value++;
  }

};
#endif