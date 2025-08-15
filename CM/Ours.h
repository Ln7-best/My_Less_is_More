#ifndef OURS_H
#define OURS_H
#include <algorithm>
#include <atomic>
#include <chrono>
#include "readerwriterqueue.h"
#include "config.h"
#include "CM.h"
#include "hash.h"
#include <barrier>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <ratio>
#include <vector>
#include <unordered_set>
#include <unordered_map>
// #define ONLINEQUERY
const uint64_t PROCESSGAP = 100000;
const uint64_t COUTERGAP = PROCESSGAP / 100;
struct alignas(64) GlobalSketchSubSection
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
  uint16_t count[COUNTER_PER_BUCKET_FILTER];
  uint16_t pos[COUNTER_PER_BUCKET_FILTER];
};

template <typename Key>
struct alignas(64) GlobalBucketsSubSection
{
  HH_Bucket<Key> buckets[BUCKET_LENGTH];
  char padding[64 - (sizeof(buckets) % 64)];
};

template <typename Key>
struct alignas(64) ChildBucketsSubSection
{
  Stash_Bucket<Key> buckets[HASH_NUM][FILTER_BUCKET_LENGTH];
  // char padding[64 - sizeof(buckets) % 64];
};

#ifdef ONLINEQUERY
template <typename T>
struct alignas(64) CacheArray
{
  T array[ARRAY_SIZE];
  char padding0[64 - ((sizeof(array)) % 64)];
};

template <typename T>
struct alignas(64) DoubleCacheArray
{
  CacheArray<T> double_outcome_cache[2];
};
template <typename T>
struct alignas(64) QueryOutcome
{
  DoubleCacheArray<T> outcome_cache;
  uint64_t process_snapshot[NUM_OUTCOME];
  char padding0[64 - ((sizeof(process_snapshot) + sizeof(outcome_cache)) % 64)];
  T *outcome[NUM_OUTCOME];
  char padding[64 - (sizeof(outcome)) % 64];

  QueryOutcome()
  {
    for (int i = 0; i < NUM_OUTCOME; i++)
    {
      if (ARRAY_SIZE < BUCKET_LENGTH * COUNTER_PER_BUCKET)
      {
        std::cerr << "ARRAY_SIZE is too small, adjust it in config.h" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      outcome[i] = static_cast<T *>(
          aligned_alloc(64, sizeof(T) * ARRAY_SIZE));
      if (!outcome[i])
      {
        std::cerr << "Memory allocation failed for outcome[" << i << "]\n";
        std::exit(EXIT_FAILURE);
      }
      volatile char *touch = reinterpret_cast<volatile char *>(outcome[i]);
      for (size_t j = 0; j < sizeof(T) * ARRAY_SIZE; j += 4096)
      {
        touch[j] = 0;
      }
    }
  }

  ~QueryOutcome()
  {
    for (int i = 0; i < NUM_OUTCOME; i++)
    {
      std::free(outcome[i]);
    }
  }
};
#endif
/**
 * A parallel framework for Count Min Sketch
 *
 * @param Key Type of the key
 *
 * @param thread_num Number of the threads
 */
template <typename Key, uint32_t thread_num>
class Ours
{
public:
  void Update(void *start, uint64_t size, std::unordered_map<Key, int32_t> mp)
  {
    size = size;
    std::thread parent;
    parent = std::thread(&Ours::ParentThread, this, &parent, start, size, &mp);
    parent.join();
  }

private:
#ifdef ONLINEQUERY
  Paddedint finish_cnt;
  Paddedint issue_cnt;
  Paddedatomicint64 query_finish;
  Paddedatomicint64 real_value_for_query[thread_num];
  Paddedint thread_snapshot_round[thread_num];
  Paddedatomicint32 query_flag[thread_num];
  typedef std::pair<Key, uint64_t> KV_Pair;
  QueryOutcome<KV_Pair> threads_outcome[thread_num];
#endif
  double running_time[thread_num];
  typedef std::unordered_map<Key, int32_t> HashMap;
  GlobalSketchSubSection global_sketch[thread_num];
  uint32_t dataset_size[thread_num];
  std::atomic<uint32_t> partition_num;
  Paddedatomic process_counter[thread_num];
  Paddedint operation_cnt[thread_num];
  GlobalBucketsSubSection<Key> global_buckets[thread_num];
  ChildBucketsSubSection<Key> stash[thread_num];
  typedef ReaderWriterQueue<CM_Entry<Key>> myQueue;
  struct Paddedint global_cnt[thread_num];
  struct Paddedint round[thread_num];
  myQueue que[thread_num][thread_num];
  MyChild_CM<Key> *initialize_child() { return new MyChild_CM<Key>(); }
  static void Partition(Key *start, uint64_t size, uint32_t id, std::vector<Key> &vec)
  {
    uint64_t interval = size / thread_num;
    for (uint64_t i = interval * id; i < interval * (id + 1); i++)
    {
      vec.push_back(start[i]);
    }
  }

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

  HashMap GetHHCandidates()
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
  /**
   * @brief function of the parent thread, used to initialize the variables and start the worker threads
   * @param thisThd pointer to the parent thread
   * @param start pointer to the beginning of the dataset
   * @param size size of the dataset in bytes
   */
  void ParentThread(std::thread *thisThd, void *start, uint64_t size, std::unordered_map<Key, int32_t> *mp)
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
          stash[t].buckets[h][b].vote = 0;
          for (uint64_t i = 0; i < COUNTER_PER_BUCKET_FILTER; i++)
          {
            stash[t].buckets[h][b].ID[i] = 0;
            stash[t].buckets[h][b].count[i] = 0;
            stash[t].buckets[h][b].pos[i] = 0;
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
#ifdef ONLINEQUERY
      real_value_for_query[i].value.store(0);
      thread_snapshot_round[i].value = 0;
      query_flag[i].value.store(0);
#endif
      process_counter[i].value.store(0);
      operation_cnt[i].value = 0;
      global_cnt[i].value = 0;
      round[i].value = 0;
    }

    for (uint32_t i = 0; i < thread_num; ++i)
    {
      thd[i] = std::thread(&Ours::ChildThread, this, &(thd[i]), i, start, size,
                           &finish);
    }

    while (partition_num < thread_num)
    {
    }
#ifdef ONLINEQUERY
    Query(finish);
#endif
    uint64_t tot_keys = 0;
    for (uint32_t i = 0; i < thread_num; ++i)
    {
      thd[i].join();
      tot_keys += dataset_size[i];
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

    double min_throughput = tot_size / max_time * 1000;
    double avg_throughput = tot_size / avg_time * 1000;
    std::cout << "tot_process" << tot_size << std::endl;
    std::cout << "avg: " << avg_time << std::endl;
    std::cout << "max: " << max_time << std::endl;
    std::cout << "avg throughput: " << avg_throughput << std::endl;
    std::cout << "min throughput: "
              // << std::fixed << std::setprecision(2)
              << min_throughput << std::endl;
    HashMap ret = GetHHCandidates();
    HHCompare(ret,(*mp), size / sizeof(Key) * ALPHA);
    std::cout << sizeof(stash[0].buckets[0][0])<<" "<<alignof(uint16_t)<<std::endl;
    std::cout << &stash[0].buckets[0][0].vote<<" "<<&stash[0].buckets[0][0].ID[0]<<" "<<&stash[0].buckets[0][0].count[0]<<" "<<&stash[0].buckets[0][0].pos[0]<<" "<<&stash[0].buckets[0][1].vote<<std::endl;

    // std::cout << "Insert " << tot_keys << " keys." << std::endl;
  }
  /**
   * @brief function of the worker thread, used to insert the elements in dataset into the sketch
   * @param thisThd pointer to the worker thread
   * @param thread_id id of the thread
   * @param start pointer to the beginning of the dataset
   * @param size size of the dataset in bytes
   * @param finish counter to record the number of finished worker
   */
  void ChildThread(std::thread *thisThd, uint32_t thread_id, void *start,
                   uint64_t size, std::atomic<int32_t> *finish)
  {
#ifdef __linux__
    if (!setaffinity(thisThd, thread_id))
      return;
#endif
    MyChild_CM<Key> *sketch = initialize_child();
    std::vector<Key> dataset;
    Partition((Key *)start, size / sizeof(Key), thread_id,
              dataset);
    partition_num.fetch_add(1);
    while (partition_num < thread_num)
    {
    }
    auto start_time = std::chrono::high_resolution_clock::now();
    Insert(dataset, sketch, que, thread_id);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end_time - start_time;
    running_time[thread_id] = duration.count();
    dataset_size[thread_id] = dataset.size();
    (*finish).fetch_add(1, std::memory_order_seq_cst);
    while ((*finish).load(std::memory_order_seq_cst) < thread_num)
    {
      ProcessQueue(que, global_sketch, thread_id);
    }
#ifdef ONLINEQUERY
    finish_cnt.value = thread_num;
#endif
    delete sketch;
  }

#ifdef ONLINEQUERY
  void IssueQuery()
  {
    bool finish_flag[thread_num];
    finish_cnt.value = 0;
    for (uint64_t i = 0; i < thread_num; i++)
    {
      finish_flag[i] = false;
    }
    query_finish.value.store(0);
    uint64_t hit_cnt = 0;
    uint64_t round = issue_cnt.value;
    for (uint64_t i = 0; i < thread_num; i++)
    {
      threads_outcome[i].process_snapshot[round] = real_value_for_query[i].value.load();
    }
    // a lightweight RCU based strategy for single reader to get the Heavy Hitter candidates from the HH keeper
    while (finish_cnt.value < thread_num)
    {
      for (uint64_t i = 0; i < thread_num; i++)
      {
        if (finish_flag[i])
          continue;
        // retrieve the snapshot buffer to read according to the writer flag
        // occupy the free snapshot buffer by setting the reader counter to 1
        uint32_t expected[2] = {0x00000000, 0x01000000};
        uint32_t swap[2] = {0x00010000, 0x01000100};
        bool flag0 = query_flag[i].value.compare_exchange_strong(expected[0], swap[0]);
        bool flag1 = query_flag[i].value.compare_exchange_strong(expected[1], swap[1]);
        uint64_t idx = 0;
        if (flag0)
        {
          idx = 1;
        }
        else if (flag1)
        {
          idx = 0;
        }
        else
        {
          continue;
        }
        // read the snapshot buffer
        for (uint64_t j = 0; j < BUCKET_LENGTH * COUNTER_PER_BUCKET; j++)
          threads_outcome[i].outcome[round][j] = this->threads_outcome[i].outcome_cache.double_outcome_cache[idx].array[j];
        // reset the counter to 0 to release the snapshot buffer back to the writer
        this->query_flag[i].value.fetch_and(0x11000000);
        finish_flag[i] = true;
        finish_cnt.value++;
        if (finish_cnt.value >= thread_num)
          break;
      }
    }
  }
  /**
   * @brief function of the query thread, used to issue the query
   * @param finish counter to record the number of finished worker
   */
  void Query(std::atomic<int32_t> &finish)
  {
    uint64_t cnt = 0;
    issue_cnt.value = 0;
    while (finish < thread_num)
    {
      issue_cnt.value++;
      std::this_thread::sleep_for(std::chrono::microseconds(1));
      IssueQuery();
    }
    std::cout << "Issue count:" << issue_cnt.value << std::endl;
  }
#endif

  void Insert(const std::vector<Key> &dataset, MyChild_CM<Key> *sketch,
              myQueue queue_group[][thread_num], uint32_t thread_id)
  {
    uint32_t length = dataset.size();
    for (uint32_t i = 0; i < length; ++i)
    {
      InsertChild(sketch, queue_group, dataset[i], thread_id, global_sketch);
      if (i % COUTERGAP == 0)
        process_counter[thread_id].value.store(i + 1, std::memory_order_relaxed);
      if (i % PROCESSGAP == 0)
        ProcessQueue(queue_group, global_sketch, thread_id);
#ifdef ONLINEQUERY
      if (i % 500 == 0)
        real_value_for_query[thread_id].value.store(i + 1, std::memory_order_relaxed);
#endif
    }
    process_counter[thread_id].value.store(length);
  }
  /**
   * @brief insert the key into local sketch
   * @param p local sketch
   * @param q_group parallel ringbuffer group for communication between local structures and global partitions
   * @param key the key to be inserted
   * @param thread_id id of the thread owning the local sketch
   * @param sketch_sub_section the partitioned global sketch
   */
  void InsertChild(MyChild_CM<Key> *p, myQueue q_group[][thread_num],
                   const Key &key,
                   uint64_t thread_id,
                   struct GlobalSketchSubSection *sketch_sub_section)
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
      // update local sketch
      sketch[hashPos][pos[hashPos]] += 1;
      if (__builtin_expect(sketch[hashPos][pos[hashPos]] >= COUNTERMAX, 0))
      {
        sketch[hashPos][pos[hashPos]] = 0;
        // if the counter is larger than the predefined threshold
        // insert the counter into remote stash for batching
        InsertStash(q_group, key, thread_id, sketch_sub_section, hashPos, pos[hashPos]);
      }
    }
  }

  inline void InsertStash(myQueue q_group[][thread_num],
                          const Key &key,
                          uint64_t thread_id,
                          struct GlobalSketchSubSection *sketch_sub_section, uint16_t hashPos, uint16_t pos)
  {
    bool update_flag = false;
    uint64_t minVal = 0x7fffffff;
    uint32_t minPos = 0;
    // ensure the update events from the same counter will be mapped into the same bucket
    uint64_t idx = pos % FILTER_BUCKET_LENGTH;
    for (uint32_t j = 0; j < COUNTER_PER_BUCKET_FILTER; j++)
    {
      if (this->stash[thread_id].buckets[hashPos][idx].ID[j] ==
          key)
      {
        this->stash[thread_id].buckets[hashPos][idx].count[j]++;
        if (this->stash[thread_id].buckets[hashPos][idx].count[j] == 9)
        {
          uint64_t push_sec_id = pos / sub_sketch_length;
          uint16_t sec_inner_index = pos % sub_sketch_length;
          // if the batched counter is large enough, push it into the global sketch via the ringbuffer
          while (__builtin_expect(!q_group[push_sec_id][thread_id].enqueue_fast(
                                      CM_Entry<Key>(key, hashPos, sec_inner_index,
                                                    9 * COUNTERMAX)),
                                  0))
          {
            ProcessQueue(q_group, sketch_sub_section, thread_id);
          }
          this->stash[thread_id].buckets[hashPos][idx].count[j] = 0;
          this->stash[thread_id].buckets[hashPos][idx].ID[j] = 0;
          this->stash[thread_id].buckets[hashPos][idx].pos[j] = 0;
        }
        update_flag = true;
      }
      if (this->stash[thread_id].buckets[hashPos][idx].count[j] <
          minVal)
      {
        minPos = j;
        minVal = this->stash[thread_id].buckets[hashPos][idx].count[j];
      }
    }
    if (update_flag)
      return;
    uint64_t insert_value = 1;
    uint64_t insert_pos = pos;
    uint64_t insert_key = key;
    // use vote strategy to simulate the TTL of each bucket
    // if the vote is large enough, the TTL of the stored element is smaller than zero
    // the stored element with TTL smaller than zero needs to be replaced

    // if the incoming element locates in the same counter as the stored element, the stored element needs to be replaced
    // this guarantees sequential consistency of updates within a single counter.
    if ((this->stash[thread_id].buckets[hashPos][idx].vote + 1) >= minVal * 8 || insert_pos == this->stash[thread_id].buckets[hashPos][idx].pos[minPos])
    {
      insert_value = minVal;
      insert_key = this->stash[thread_id].buckets[hashPos][idx].ID[minPos];
      insert_pos = this->stash[thread_id].buckets[hashPos][idx].pos[minPos];
      this->stash[thread_id].buckets[hashPos][idx].vote = 0;
      this->stash[thread_id].buckets[hashPos][idx].ID[minPos] =
          key;
      this->stash[thread_id].buckets[hashPos][idx].count[minPos] =
          1;
      this->stash[thread_id].buckets[hashPos][idx].pos[minPos] =
          pos;
    }
    else
    {
      this->stash[thread_id].buckets[hashPos][idx].vote++;
    }
    if (insert_value == 0)
      return;
    uint64_t push_sec_id = insert_pos / sub_sketch_length;
    uint16_t sec_inner_index = insert_pos % sub_sketch_length;
    // push the evicted element into the global sketc
    while (__builtin_expect(!q_group[push_sec_id][thread_id].enqueue_fast(CM_Entry<Key>(
                                insert_key, hashPos, sec_inner_index, insert_value * COUNTERMAX)),
                            0))
    {
      ProcessQueue(q_group, sketch_sub_section, thread_id);
    }
  }
  /**
   * @brief merge the updates into the global sketch partition
   * @param p local sketch
   * @param q_group parallel ringbuffer group for communication between local structures and global partitions
   * @param key the key to be inserted
   * @param thread_id id of the thread owning the partition
   * @param sketch_sub_section the partitioned global sketch
   */
  inline void ProcessQueue(myQueue q_group[][thread_num],
                           struct GlobalSketchSubSection *sketch_sub_section,
                           uint64_t thread_id = 0)
  {
    CM_Entry<Key> temp;
    uint64_t cnt = 0;
    for (uint64_t i = 0; i < thread_num; i++)
    {
      while (q_group[thread_id][i].try_dequeue(temp))
      {
        // merge the update into the global partition
        const uint16_t &hash_pos = temp.hashPos;
        const uint16_t &pos = temp.pos;
        uint64_t counter_pos = hash_pos * sub_sketch_length + pos;
        sketch_sub_section[thread_id].counter[counter_pos].fetch_add(
            temp.value);
        uint64_t counter_val =
            sketch_sub_section[thread_id].counter[counter_pos].load();

        // if the counter value is large enough, try to push the key into HH keeper
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

          // get the estimated frequency of the key
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

          // update HH keeper
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
#ifdef ONLINEQUERY
    // update the snapshot for the query thread
    if (cnt <= 10)
      return;
    UpdateSnapshot(thread_id);
#endif
  }
#ifdef ONLINEQUERY
  inline void UpdateSnapshot(uint64_t thread_id)
  {
    // a lightweight RCU based strategy for single writer to update the heavy hitter candidates snapshot

    // retrieve the snapshot buffer not being written according to the reader counter
    // occupy the free snapshot buffer by setting the writer flag
    uint32_t expected[2] = {0x00000100, 0x01010000};
    uint32_t swap[2] = {0x01000100, 0x00010000};
    uint32_t local_round = this->thread_snapshot_round[thread_id].value;
    uint32_t exchange = this->query_flag[thread_id].value.compare_exchange_strong(expected[local_round], swap[local_round]);
    local_round = local_round ^ exchange;
    uint64_t idx = 0;
    for (uint64_t i = 0; i < BUCKET_LENGTH; i++)
    {
      for (uint64_t j = 0; j < COUNTER_PER_BUCKET; j++)
      {
        this->threads_outcome[thread_id].outcome_cache.double_outcome_cache[local_round].array[idx++] = std::pair<Key, uint64_t>(
            this->global_buckets[thread_id].buckets[i].ID[j],
            this->global_buckets[thread_id].buckets[i].count[j]);
      }
    }
    // reverse the writer flag to exchange the snapshot for reader to access
    this->thread_snapshot_round[thread_id].value ^= 1;
    this->query_flag[thread_id].value.fetch_xor(0x01000000);
  }
#endif
};
#endif