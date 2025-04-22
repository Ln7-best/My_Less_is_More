#ifndef OURS_H
#define OURS_H
#include <algorithm>
#include <atomic>
#include <chrono>
#include "util.h"
#include <fstream>
#include "config.h"
#include "UnivMon.h"
#include <barrier>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <ratio>
#include <vector>
#include <unordered_map>
#include "readerwriterqueue.h"
#include <unordered_set>
#include "hash.h"
#define MEASURETIME
#define ONLINEQUERY
const uint64_t PROCESSGAP = 50000;
const uint64_t COUTERGAP = PROCESSGAP / 100;

struct alignas(64) GlobalSketchSubSection
{
  std::atomic<int64_t> counter[MAX_LEVEL][HASH_NUM * sub_sketch_length];
  char padding[64 - (sizeof(counter) % 64)];
};

template <typename Key>
struct HH_Bucket
{
  uint64_t vote;
  Key ID[COUNTER_PER_BUCKET];
  int64_t count[COUNTER_PER_BUCKET];
};

template <typename Key>
struct Stash_Bucket
{
  uint64_t vote;
  Key ID[COUNTER_PER_BUCKET_FILTER];
  int64_t count[COUNTER_PER_BUCKET_FILTER];
  uint16_t pos[COUNTER_PER_BUCKET_FILTER];
};

template <typename Key>
struct alignas(64) GlobalBucketsSubSection
{
  HH_Bucket<Key> buckets[MAX_LEVEL][BUCKET_LENGTH];
  char padding[64 - sizeof(buckets) % 64];
};

template <typename Key>
struct alignas(64) ChildBucketsSubSection
{
  Stash_Bucket<Key> buckets[MAX_LEVEL][HASH_NUM][FILTER_BUCKET_LENGTH];
  char padding[64 - sizeof(buckets) % 64];
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
  DoubleCacheArray<T> outcome_cache[MAX_LEVEL];
  uint64_t process_snapshot[NUM_OUTCOME];
  char padding0[64 - ((sizeof(process_snapshot) + sizeof(outcome_cache)) % 64)];
  T *outcome[NUM_OUTCOME][MAX_LEVEL];
  char padding[64 - (sizeof(outcome)) % 64];

  QueryOutcome()
  {
    for (int i = 0; i < NUM_OUTCOME; i++)
    {
      for (int j = 0; j < MAX_LEVEL; j++)
      {
        outcome[i][j] = static_cast<T *>(
            aligned_alloc(64, sizeof(T) * ARRAY_SIZE));
        if (!outcome[i][j])
        {
          std::cerr << "Memory allocation failed for outcome[" << i << "]\n";
          std::exit(EXIT_FAILURE);
        }
        volatile char *touch = reinterpret_cast<volatile char *>(outcome[i][j]);
        for (size_t j = 0; j < sizeof(T) * ARRAY_SIZE; j += 4096)
        {
          touch[j] = 0;
        }
      }
    }
  }

  ~QueryOutcome()
  {
    for (int i = 0; i < NUM_OUTCOME; i++)
    {
      for(int j=0; j < MAX_LEVEL; j++)
        std::free(outcome[i][j]); 
    }
  }
};
#endif

template <typename Key, uint32_t thread_num>
class Ours
{
public:
  typedef std::unordered_map<Key, int32_t> HashMap;
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
  const uint64_t filter_len[6] = {64, 32, 16, 8, 4, 2};
  GlobalSketchSubSection global_sketch[thread_num];
  double running_time[thread_num];
  uint32_t dataset_size[thread_num];
  std::atomic<uint32_t> partition_num;
  Paddedatomic process_counter[thread_num];
  Paddedint operation_cnt[thread_num];
  GlobalBucketsSubSection<Key> global_buckets[thread_num];
  ChildBucketsSubSection<Key> stash[thread_num];
  typedef ReaderWriterQueue<Univ_Entry<Key>> myQueue;
  struct Paddedint global_cnt[thread_num];
  struct Paddedint round[thread_num];
  myQueue que[thread_num][thread_num];
  static void HHCompare(HashMap test, HashMap real, int32_t threshold, std::ofstream *outputFile = nullptr)
  {
    double estHH = 0, HH = 0, both = 0;
    double CR = 0, PR = 0, AAE = 0, ARE = 0;
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

  /**
   * The thread of the aggregator
   */
  HashMap GetHHCandidates()
  {
    HashMap ret;

    for (uint64_t l = 0; l < MAX_LEVEL; l++)
    {
      for (uint32_t i = 0; i < thread_num; i++)
      {
        for (uint32_t j = 0; j < BUCKET_LENGTH; j++)
        {
          for (uint32_t k = 0; k < COUNTER_PER_BUCKET; k++)
          {
            if (global_buckets[i].buckets[l][j].count[k] == 0)
              continue;
            if (ret.find(global_buckets[i].buckets[l][j].ID[k]) == ret.end())
              ret[global_buckets[i].buckets[l][j].ID[k]] = global_buckets[i].buckets[l][j].count[k];
          }
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
    std::thread measure_thd;
    partition_num.store(0);
    for (uint64_t l = 0; l < MAX_LEVEL; l++)
    {
      for (uint64_t i = 0; i < thread_num; i++)
      {
        for (uint64_t j = 0; j < HASH_NUM * sub_sketch_length; j++)
        {
          global_sketch[i].counter[l][j] = 0;
        }
      }
    }
    for (uint64_t l = 0; l < MAX_LEVEL; l++)
    {
      for (size_t t = 0; t < thread_num; ++t)
      {
        for (size_t h = 0; h < HASH_NUM; ++h)
        {
          for (size_t b = 0; b < FILTER_BUCKET_LENGTH; ++b)
          {
            stash[t].buckets[l][h][b].vote = 0;
            for (uint64_t i = 0; i < COUNTER_PER_BUCKET_FILTER; i++)
            {
              stash[t].buckets[l][h][b].ID[i] = 0;
              stash[t].buckets[l][h][b].count[i] = 0;
              stash[t].buckets[l][h][b].pos[i] = 0;
            }
          }
        }
      }
    }
    for (uint64_t l = 0; l < MAX_LEVEL; l++)
    {
      for (uint64_t i = 0; i < thread_num; i++)
      {
        for (uint64_t j = 0; j < BUCKET_LENGTH; j++)
        {
          global_buckets[i].buckets[l][j].vote = 0;
          for (uint64_t k = 0; k < COUNTER_PER_BUCKET; k++)
          {
            global_buckets[i].buckets[l][j].ID[k] = 0;
            global_buckets[i].buckets[l][j].count[k] = 0;
          }
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
      global_cnt[i].value = 0;
      round[i].value = 0;
      process_counter[i].value.store(0);
      thd[i] = std::thread(&Ours::ChildThread, this, &(thd[i]), i, start, size,
                           &finish);
    }
    while (partition_num < thread_num)
    {
    }
#ifdef ONLINEQUERY
    Query(finish);
#endif
    for (uint32_t i = 0; i < thread_num; ++i)
    {
      thd[i].join();
    }
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
                << " running time: " << running_time[i]
                << std::endl;
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
    for (uint64_t i = 0; i < thread_num; i++)
    {
      double max = 0;
      double min = 1e9;
      double tot = 0;
      uint64_t expansion_cnt = 0;
      for (uint64_t j = 0; j < thread_num; j++)
      {

        for (auto time : que[j][i].make_block_time)
        {
          max = std::max(time, max);
          min = std::min(time, min);
          tot += time;
          expansion_cnt += que[j][i].expansion_cnt;
        }
      }

      std::cout << " thread: " << i << " expansion count: " << expansion_cnt
                << " max: " << std::fixed << std::setprecision(10) << max
                << " min: " << std::fixed << std::setprecision(10) << min
                << " tot: " << std::fixed << std::setprecision(10) << tot
                << std::endl;
    }

    for (uint64_t i = 0; i < thread_num; i++)
    {
      double max = 0;
      double min = 1e9;
      double tot = 0;
      uint64_t expansion_cnt = 0;
      for (uint64_t j = 0; j < thread_num; j++)
      {
        for (auto time : que[i][j].make_block_time)
        {
          max = std::max(time, max);
          min = std::min(time, min);
          tot += time;
          expansion_cnt += que[i][j].expansion_cnt;
        }
      }
      std::cout << " thread: " << i << " expansion count: " << expansion_cnt
                << " max: " << std::fixed << std::setprecision(10) << max
                << " min: " << std::fixed << std::setprecision(10) << min
                << " tot: " << std::fixed << std::setprecision(10) << tot
                << std::endl;
    }
    // HashMap ret = query_all();
    // HHCompare(ret, (*mp), size / sizeof(Key) * ALPHA);
  }

  /**
   * The thread of each worker
   */
  void ChildThread(std::thread *thisThd, uint32_t thread_id, void *start,
                   uint64_t size, std::atomic<int32_t> *finish)
  {
#ifdef __linux__
    if (!setaffinity(thisThd, thread_id))
      return;
#endif
    MyChild_Univ<Key> *sketch = initialize_child();
    uint64_t local_min = 0;
    std::vector<Key> dataset;
    Partition((Key *)start, size / sizeof(Key), thread_id,
              dataset);
    partition_num.fetch_add(1);
    while (partition_num < thread_num)
    {
    }
#ifdef MEASURETIME
    auto start_time = std::chrono::high_resolution_clock::now();
#endif
    insert(dataset, sketch, que, thread_id);

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
    // a lightweight RCU strategy for single reader to get the Heavy Hitter candidates from the HH keeper
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
        for (uint64_t k = 0; k < MAX_LEVEL; k++)
        {
          for (uint64_t j = 0; j < BUCKET_LENGTH * COUNTER_PER_BUCKET; j++)
            threads_outcome[i].outcome[round][k][j] = this->threads_outcome[i].outcome_cache[k].double_outcome_cache[idx].array[j];
        }
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
    uint64_t sum = 0;
    while (finish < thread_num)
    {
      issue_cnt.value++;
      std::this_thread::sleep_for(std::chrono::microseconds(1));
      IssueQuery();
    }
    std::cout << "issue count:" << issue_cnt.value << std::endl;
  }
#endif

  void insert(const std::vector<Key> &dataset, MyChild_Univ<Key> *sketch,
              myQueue queue_group[][thread_num], uint32_t thread_id)
  {
    uint32_t length = dataset.size();
    for (uint32_t i = 0; i < length; ++i)
    {
      insert_child(sketch, queue_group, dataset[i], thread_id, global_sketch);
      if (i % COUTERGAP == 0)
        process_counter[thread_id].value.store(i + 1, std::memory_order_relaxed);
      if (i % PROCESSGAP == 0)
        process_queue(queue_group, global_sketch,
                      thread_id);
#ifdef ONLINEQUERY
      if (i % 500 == 0)
        real_value_for_query[thread_id].value.store(i + 1, std::memory_order_relaxed);
#endif
    }
    process_counter[thread_id].value.store(length);
  }
  MyChild_Univ<Key> *initialize_child()
  {
    return new MyChild_Univ<Key>();
  }

  void insert_child(MyChild_Univ<Key> *p, myQueue q_group[][thread_num],
                    const Key &key,
                    uint64_t thread_id,
                    struct GlobalSketchSubSection *sketch_sub_section)
  {
    auto sketch = p->sketch;

    uint32_t polar = hash(key, 199);
    uint32_t max_level = MIN(MAX_LEVEL - 1, __builtin_clz(polar));

    for (uint32_t level = 0; level <= max_level; ++level)
    {

      uint32_t pos[HASH_NUM];
      int32_t incre[HASH_NUM];

      for (uint32_t hashPos = 0; hashPos < HASH_NUM; ++hashPos)
      {
        uint32_t hashNum = hash(key, level * HASH_NUM + hashPos);
        pos[hashPos] = (hashNum >> 1) % LENGTH;
        incre[hashPos] = increment[hashNum & 1];
      }

      for (uint32_t hashPos = 0; hashPos < HASH_NUM; ++hashPos)
      {
        sketch[level][hashPos][pos[hashPos]] += incre[hashPos];

        if (sketch[level][hashPos][pos[hashPos]] * incre[hashPos] >= PROMASK)
        {

          sketch[level][hashPos][pos[hashPos]] = 0;
          InsertStash(q_group, key, thread_id, sketch_sub_section, hashPos, pos[hashPos], level, incre);
        }
      }
    }
  }

  inline void InsertStash(myQueue q_group[][thread_num],
                          const Key &key,
                          uint64_t thread_id,
                          struct GlobalSketchSubSection *sketch_sub_section, uint16_t hashPos, uint16_t pos, uint16_t level, int32_t *incre)
  {
    bool update_flag = false;
    int64_t minVal = 0xffffffff;
    uint32_t minPos = 0;
    uint64_t idx = pos % filter_len[level];
    for (uint32_t j = 0; j < COUNTER_PER_BUCKET_FILTER; j++)
    {
      if (this->stash[thread_id].buckets[level][hashPos][idx].ID[j] ==
          key)
      {
        this->stash[thread_id].buckets[level][hashPos][idx].count[j] +=
            incre[hashPos];
        if (this->stash[thread_id].buckets[level][hashPos][idx].count[j] *
                incre[hashPos] ==
            9)
        {
          uint64_t push_sec_id = pos / sub_sketch_length;
          uint16_t sec_inner_index = pos % sub_sketch_length;
          while (__builtin_expect(
              !q_group[push_sec_id][thread_id].enqueue_fast(
                  Univ_Entry<Key>(key, level, hashPos, sec_inner_index,
                                  9 * PROMASK * incre[hashPos])),
              0))
          {
            process_queue(q_group, sketch_sub_section, thread_id);
          }
          this->stash[thread_id].buckets[level][hashPos][idx].count[j] = 0;
          this->stash[thread_id].buckets[level][hashPos][idx].ID[j] = 0;
          this->stash[thread_id].buckets[level][hashPos][idx].pos[j] = 0;
        }
        update_flag = true;
      }
      minPos = j;
      minVal = abs(this->stash[thread_id].buckets[level][hashPos][idx].count[j]);
    }
    if (update_flag)
      return;
    int64_t insert_value = incre[hashPos];
    uint64_t insert_pos = pos;
    uint64_t insert_key = key;
    if ((this->stash[thread_id].buckets[level][hashPos][idx].vote + 1) >=
            minVal * 8 ||
        insert_pos == this->stash[thread_id].buckets[level][hashPos][idx].pos[minPos])
    {
      insert_value = this->stash[thread_id]
                         .buckets[level][hashPos][idx]
                         .count[minPos];
      insert_key =
          this->stash[thread_id].buckets[level][hashPos][idx].ID[minPos];
      insert_pos =
          this->stash[thread_id].buckets[level][hashPos][idx].pos[minPos];
      this->stash[thread_id].buckets[level][hashPos][idx].vote = 0;
      this->stash[thread_id].buckets[level][hashPos][idx].ID[minPos] =
          key;
      this->stash[thread_id].buckets[level][hashPos][idx].count[minPos] =
          incre[hashPos];
      this->stash[thread_id].buckets[level][hashPos][idx].pos[minPos] = pos;
    }
    else
    {
      this->stash[thread_id].buckets[level][hashPos][idx].vote++;
    }
    if (insert_value == 0)
      return;
    uint64_t push_sec_id = insert_pos / sub_sketch_length;
    uint16_t sec_inner_index = insert_pos % sub_sketch_length;

    while (__builtin_expect(
        !q_group[push_sec_id][thread_id].enqueue_fast(Univ_Entry<Key>(
            insert_key, level, hashPos, sec_inner_index, insert_value * PROMASK)),
        0))
    {
      process_queue(q_group, sketch_sub_section, thread_id);
    }
  }

  inline void
  process_queue(myQueue q_group[][thread_num],
                struct GlobalSketchSubSection *sketch_sub_section,
                uint64_t thread_id = 0)
  {
    uint64_t process_cnt = 0;
    Univ_Entry<Key> temp;
    uint64_t que_cnt = 0;
    uint64_t candidate_cnt = 0;
    uint64_t cnt = 0;
    bool level_update_flag[MAX_LEVEL]={false};
    for (uint64_t i = 0; i < thread_num; i++)
    {
      while (q_group[thread_id][i].try_dequeue(temp))
      {
        process_cnt++;
        const uint16_t &level = temp.level;
        const uint16_t &hash_pos = temp.hashPos;
        const uint16_t &pos = temp.pos;
        uint64_t counter_pos = hash_pos * sub_sketch_length + pos;
        sketch_sub_section[thread_id].counter[level][counter_pos].fetch_add(temp.value);
        int64_t counter_val = sketch_sub_section[thread_id].counter[level][counter_pos].load();
        if ((double)abs(counter_val) >
            (double)this->global_cnt[thread_id].value * (ALPHA / 2))
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
          candidate_cnt++;
          int32_t count[HASH_NUM] = {0};
          for (uint32_t tempHash = 0; tempHash < HASH_NUM; ++tempHash)
          {
            uint32_t hashNum = hash(temp.key, tempHash);
            uint32_t tempPos = (hashNum >> 1) % LENGTH;
            int32_t incre = increment[hashNum & 1];
            uint64_t push_sec_id = tempPos / sub_sketch_length;
            uint16_t sec_inner_index = tempPos % sub_sketch_length;
            count[tempHash] = sketch_sub_section[push_sec_id]
                                  .counter[level][tempHash * sub_sketch_length + sec_inner_index]
                                  .load() *
                              incre;
          }
          int64_t minimum = MEDIAN3(count);
          if (minimum <= 0)
          {
            continue;
          }
          que_cnt++;
          int64_t minVal = 0xffffffff;
          uint32_t minPos = 0;
          uint64_t idx = hash(temp.key, 101 + level) % BUCKET_LENGTH;
          bool update_flag = 0;
          level_update_flag[level] = true;
          for (uint32_t j = 0; j < COUNTER_PER_BUCKET; j++)
          {
            if (this->global_buckets[thread_id].buckets[level][idx].ID[j] ==
                temp.key)
            {
              update_flag = 1;
              this->global_buckets[thread_id].buckets[level][idx].count[j] = minimum;
              break;
            }
            if (this->global_buckets[thread_id].buckets[level][idx].count[j] <
                minVal)
            {
              minPos = j;
              minVal = this->global_buckets[thread_id].buckets[level][idx].count[j];
            }
          }
          if (!update_flag)
          {
            if ((this->global_buckets[thread_id].buckets[level][idx].vote + minimum) >=
                minVal * 8)
            {
              this->global_buckets[thread_id].buckets[level][idx].vote = 0;
              this->global_buckets[thread_id].buckets[level][idx].ID[minPos] =
                  temp.key;
              this->global_buckets[thread_id].buckets[level][idx].count[minPos] =
                  minimum;
            }
            else
            {
              this->global_buckets[thread_id].buckets[level][idx].vote += minimum;
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
    UpdateSnapshot(thread_id, level_update_flag);
#endif
  }

#ifdef ONLINEQUERY
  inline void UpdateSnapshot(uint64_t thread_id, bool* level_update_flag)
  {
    // a lightweight RCU strategy for single writer to update the heavy hitter candidates snapshot

    // retrieve the snapshot buffer not being written according to the reader counter
    // occupy the free snapshot buffer by setting the writer flag
    uint32_t expected[2] = {0x00000100, 0x01010000};
    uint32_t swap[2] = {0x01000100, 0x00010000};
    uint32_t local_round = this->thread_snapshot_round[thread_id].value;
    uint32_t exchange = this->query_flag[thread_id].value.compare_exchange_strong(expected[local_round], swap[local_round]);
    local_round = local_round ^ exchange;
    uint64_t idx = 0;
    for(uint64_t k = 0; k<MAX_LEVEL; k++)
    {
      if(!level_update_flag[k])
        continue;
      uint64_t idx = 0;
      for (uint64_t i = 0; i < BUCKET_LENGTH; i++) {
        for (uint64_t j = 0; j < COUNTER_PER_BUCKET; j++) {
          this->threads_outcome[thread_id].outcome_cache[k].double_outcome_cache[local_round].array[idx++] = std::pair<Key, uint64_t>(
            this->global_buckets[thread_id].buckets[k][i].ID[j],
            this->global_buckets[thread_id].buckets[k][i].count[j]);
        }
      }
    }
    // reverse the writer flag to exchange the snapshot for reader to access
    this->thread_snapshot_round[thread_id].value ^= 1;
    this->query_flag[thread_id].value.fetch_xor(0x01000000);
  }
#endif

  void ProcessQuery(void *start, uint64_t size)
  {
    HashMap real;
    std::unordered_set<uint64_t> key_set;
    // for (auto it = keyset.begin(); it != keyset.end(); it++)
    //   real[*it] = 0;
    std::vector<Key> dataset[thread_num];
    for (uint32_t i = 0; i < thread_num; i++)
    {
      Partition((Key *)start, size / sizeof(Key), i, dataset[i]);
    }

#ifdef QLENGTHONLY
    std::ofstream outputFile(
        "/home/ln7/OctoSketch/eval/OctoSketchQLENSKIPHHLargeSet" +
        std::to_string(thread_num) + std::to_string(PROMASK) + ".txt");
    for (uint32_t i = 0; i < queue_lengths.size(); i++)
    {
      uint64_t tot = 0;
      for (uint32_t j = 0; j < thread_num; j++)
      {
        tot += thd_packets[j][i];
      }
      if (tot == 0)
        continue;
      outputFile << tot << " " << queue_lengths[i] << std::endl;
    }
#else
#ifdef AGGRESSIVE_PUSH
std::ofstream outputFile("OctoSketchACCLargeSetSuperACC" +
  std::to_string(thread_num) +
  std::to_string(PROMASK) + ".txt");
#else
    std::ofstream outputFile0("test.txt");
    // std::ofstream outputFile1("OctoSketchACCLargeSet0.002Snapshot" +
    //                           std::to_string(thread_num) +
    //                           std::to_string(PROMASK) + ".txt");
    // std::ofstream outputFile2("OctoSketchACCLargeSet0.02Snapshot" +
    //                             std::to_string(thread_num) +
    //                             std::to_string(PROMASK) + ".txt");
    // std::ofstream outputFile3("OctoSketchACCLargeSet0.0016" +
    //                           std::to_string(thread_num) +
    //                           std::to_string(PROMASK) + ".txt");     
    // std::ofstream outputFile4("OctoSketchACCLargeSet0.0032" +
    //                             std::to_string(thread_num) +
    //                             std::to_string(PROMASK) + ".txt");     
    // std::ofstream outputFile1("OctoSketchACCLargeSet0.02" +
    //                             std::to_string(thread_num) +
    //                             std::to_string(PROMASK) + ".txt");              
#endif
    for (uint64_t i = 0; i < issue_cnt.value - 1; i++)
    {
      uint64_t tot = 0;
      for (uint64_t j = 0; j < thread_num; j++)
      {
        tot += threads_outcome[j].process_snapshot[i];
        if (i == 0)
        {
          for (uint64_t k = 0; k < threads_outcome[j].process_snapshot[i];
               k++)
          {
            if (key_set.find(dataset[j][k]) != key_set.end())
            {
              real[dataset[j][k]] += 1;
            }
            else
            {
              key_set.insert(dataset[j][k]);
              real[dataset[j][k]] = 1;
            }
          }
        }
        else
        {
          for (uint64_t k = threads_outcome[j].process_snapshot[i - 1];
               k < threads_outcome[j].process_snapshot[i]; k++)
            if (key_set.find(dataset[j][k]) != key_set.end())
            {
              real[dataset[j][k]] += 1;
            }
            else
            {
              key_set.insert(dataset[j][k]);
              real[dataset[j][k]] = 1;
            }
        }
      }
      // std::cout << "packets num: " << tot << std::endl;
      if(i<=1e9 || i % 100 == 0)
      {
      if (tot == 0)
        continue;
       
      HashMap ret;
      for (uint64_t j = 0; j < thread_num; j++)
      {
        for(uint64_t t = 0 ;t < MAX_LEVEL; t++)
        {
          for (uint64_t k = 0; k < BUCKET_LENGTH * COUNTER_PER_BUCKET; k++)
          {
            if (threads_outcome[j].outcome[i][t][k].second == 0)
              continue;
            if(ret.find(threads_outcome[j].outcome[i][t][k].first) == ret.end())  
              ret[threads_outcome[j].outcome[i][t][k].first] = threads_outcome[j].outcome[i][t][k].second;
          }
        }   
      }
      outputFile0 << "Tot:" << std::to_string(tot) << std::endl;
      HHCompare(ret, real, tot * ALPHA, &outputFile0);
      // outputFile1 << "Tot:" << std::to_string(tot) << std::endl;
      // HHCompare(ret, real, tot * 0.002, &outputFile1);
      // outputFile2 << "Tot:" << std::to_string(tot) << std::endl;
      // HHCompare(ret, real, tot * 0.02, &outputFile2);
      // outputFile3 << "Tot:" << std::to_string(tot) << std::endl;
      // HHCompare(ret, real, tot * 0.0016, &outputFile3);
      // outputFile4 << "Tot:" << std::to_string(tot) << std::endl;
      // HHCompare(ret, real, tot * 0.0032, &outputFile4);
      }

      // if(tot>200050000)
      // {
      //   return;
      // }
    }
#endif
    outputFile0.close();
    // outputFile1.close();
    // outputFile2.close();
    // outputFile3.close();
    // outputFile4.close();
  }
};
#endif