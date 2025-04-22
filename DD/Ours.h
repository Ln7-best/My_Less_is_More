#ifndef OURS_H
#define OURS_H
#include <algorithm>
#include <atomic>
#include <chrono>
#include "readerwriterqueue.h"
#include "config.h"
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <ratio>
#include <vector>
#include "hash.h"
#include "util.h"
#include <unordered_map>
#define MEASURETIME
#define ONLINEQUERY
const uint64_t PROCESSGAP = 17000;
struct alignas(64) GlobalSketchSubSection
{
  uint64_t counter[sub_sketch_length];
  char padding[64 - (sizeof(counter) % 64)];
};

template <typename Key>
struct Bucket_
{
  uint64_t vote;
  Key ID[COUNTER_PER_BUCKET_FILTER];
  uint16_t count[COUNTER_PER_BUCKET_FILTER];
  uint16_t pos[COUNTER_PER_BUCKET_FILTER];
};

template <typename Key>
struct alignas(64) ChildBucketsSubSection
{
  Bucket_<Key> buckets[FILTER_BUCKET_LENGTH];
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

template <typename Key, uint32_t thread_num>
class Ours
{
public:
  typedef std::unordered_map<Key, uint64_t> HashMap;
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
  QueryOutcome<uint64_t> threads_outcome[thread_num];
#endif
  struct GlobalSketchSubSection global_sketch[thread_num];
  std::atomic<uint32_t> partition_num;
  ChildBucketsSubSection<Key> stash[thread_num];
  double running_time[thread_num];
  uint32_t dataset_size[thread_num];
  Paddedint operation_cnt[thread_num];
  typedef ReaderWriterQueue<DD_Entry<Key>> myQueue;
  myQueue que[thread_num][thread_num];
  MyChild_DD<Key> *initialize_child() { return new MyChild_DD<Key>(); }
  static void Partition(Key *start, uint64_t size, uint32_t id, std::vector<Key> &vec)
  {
    uint64_t interval = size / thread_num;
    for (uint64_t i = interval * id; i < interval * (id + 1); i++)
    {
      vec.push_back(start[i]);
    }
  }
  static void QuantCompare(HashMap test, HashMap &real)
  {

    for (int i = 0; i < 8; i++)
    {
      std::cout << i << ": " << ((double)test[i] / (double)real[i]) - 1 << std::endl;
    }
  }

  /**
   * The thread of the aggregator
   */
  HashMap GetQuantiles()
  {
    HashMap ret;
    std::vector<double> quantiles = {0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 0.999};
    double sum = 0;
    for (uint32_t i = 0; i < LENGTH; ++i)
    {
      uint64_t secid = i / sub_sketch_length;
      uint64_t inner_index = i % sub_sketch_length;
      sum += global_sketch[secid].counter[inner_index];
    }

    uint32_t index = 0;
    double preSum = 0;
    for (uint32_t i = 0; i < 8; ++i)
    {
      while (preSum <= quantiles[i] * sum)
      {
        uint64_t secid = index / sub_sketch_length;
        uint64_t inner_index = index % sub_sketch_length;
        preSum += global_sketch[secid].counter[inner_index];
        index += 1;
      }
      std::cout << index << std::endl;
      uint64_t secid = (index - 1) / sub_sketch_length;
      uint64_t inner_index = (index - 1) % sub_sketch_length;
      std::cout << global_sketch[secid].counter[inner_index] << std::endl;
      ret[i] = 2 * pow(GAMMA, index) / (1 + GAMMA);
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
    for (uint64_t i = 0; i < thread_num; i++)
    {
      for (uint64_t j = 0; j < sub_sketch_length; j++)
      {
        global_sketch[i].counter[j] = 0;
      }
    }
    for (size_t t = 0; t < thread_num; ++t)
    {
      for (size_t b = 0; b < FILTER_BUCKET_LENGTH; ++b)
      {
        stash[t].buckets[b].vote = 0;
        for (uint64_t i = 0; i < COUNTER_PER_BUCKET_FILTER; i++)
        {
          stash[t].buckets[b].ID[i] = 0;
          stash[t].buckets[b].count[i] = 0;
          stash[t].buckets[b].pos[i] = 0;
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
  }

  void ChildThread(std::thread *thisThd, uint32_t thread_id, void *start,
                   uint64_t size, std::atomic<int32_t> *finish)
  {
#ifdef __linux__
    if (!setaffinity(thisThd, thread_id))
      return;
#endif
    MyChild_DD<Key> *sketch = initialize_child();
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
        // get the snapshot of the buffer
        // any further queries for quantiles can be completed by reading this buffer
        for (uint64_t j = 0; j < sub_sketch_length; j++)
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

  void insert(const std::vector<Key> &dataset, MyChild_DD<Key> *sketch,
              myQueue queue_group[][thread_num], uint32_t thread_id)
  {
    uint32_t length = dataset.size();
    for (uint32_t i = 0; i < length; ++i)
    {
      insert_child(sketch, queue_group, dataset[i], thread_id, global_sketch);
      if (i % PROCESSGAP == 0)
        process_queue(queue_group, global_sketch,
                      thread_id);
#ifdef ONLINEQUERY
      if ((i + 1) % 50 == 0)
        real_value_for_query[thread_id].value.store(i + 1, std::memory_order_relaxed);
#endif
    }
  }

  void insert_child(MyChild_DD<Key> *p, myQueue q_group[][thread_num],
                    const Key &key,
                    uint64_t thread_id,
                    struct GlobalSketchSubSection *sketch_sub_section)
  {

    auto sketch = p->sketch;
    auto log_gamma = p->log_gamma;
    uint32_t pos = log(key) / log_gamma;
    if (pos >= LENGTH)
    {
      pos = LENGTH - 1;
    }
    sketch[pos] += 1;
    if (__builtin_expect(sketch[pos] >= PROMASK, 0))
    {
      sketch[pos] = 0;
      InsertStash(q_group, key, thread_id, sketch_sub_section, pos);
    }
  }

  inline void InsertStash(myQueue q_group[][thread_num],
                          const Key &key,
                          uint64_t thread_id,
                          struct GlobalSketchSubSection *sketch_sub_section, uint16_t pos)
  {
    bool update_flag = false;
    int64_t minVal = 0xffffffff;
    uint32_t minPos = 0;
    uint64_t idx = pos % FILTER_BUCKET_LENGTH;
    for (uint32_t j = 0; j < COUNTER_PER_BUCKET_FILTER; j++)
    {
      if (this->stash[thread_id].buckets[idx].ID[j] ==
          key)
      {
        this->stash[thread_id].buckets[idx].count[j]++;
        if (this->stash[thread_id].buckets[idx].count[j] == 9)
        {
          uint64_t push_sec_id = pos / sub_sketch_length;
          uint16_t sec_inner_index = pos % sub_sketch_length;
          while (__builtin_expect(
              !q_group[push_sec_id][thread_id].enqueue_fast(
                  DD_Entry<Key>(sec_inner_index, 9 * PROMASK)),
              0))
          {
            process_queue(q_group, sketch_sub_section, thread_id);
          }
          this->stash[thread_id].buckets[idx].count[j] = 0;
          this->stash[thread_id].buckets[idx].ID[j] = 0;
          this->stash[thread_id].buckets[idx].pos[j] = 0;
        }
        update_flag = true;
      }
      if (this->stash[thread_id].buckets[idx].count[j] < minVal)
      {
        minPos = j;
        minVal = this->stash[thread_id].buckets[idx].count[j];
      }
    }
    if (update_flag)
      return;
    int64_t insert_value = 1;
    uint64_t insert_pos = pos;
    if ((this->stash[thread_id].buckets[idx].vote + 1) >=
            minVal * 8 ||
        insert_pos == this->stash[thread_id].buckets[idx].pos[minPos])
    {
      insert_value = this->stash[thread_id]
                         .buckets[idx]
                         .count[minPos];
      insert_pos =
          this->stash[thread_id].buckets[idx].pos[minPos];
      this->stash[thread_id].buckets[idx].vote = 0;
      this->stash[thread_id].buckets[idx].ID[minPos] =
          key;
      this->stash[thread_id].buckets[idx].count[minPos] = 1;
      this->stash[thread_id].buckets[idx].pos[minPos] = pos;
    }
    else
    {
      this->stash[thread_id].buckets[idx].vote++;
    }
    if (insert_value == 0)
      return;
    uint64_t push_sec_id = insert_pos / sub_sketch_length;
    uint16_t sec_inner_index = insert_pos % sub_sketch_length;

    while (__builtin_expect(
        !q_group[push_sec_id][thread_id].enqueue_fast(DD_Entry<Key>(sec_inner_index, insert_value * PROMASK)),
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
    DD_Entry<Key> temp;
    uint64_t que_cnt = 0;
    uint64_t candidate_cnt = 0;
    uint64_t cnt = 0;
    for (uint64_t i = 0; i < thread_num; i++)
    {
      while (q_group[thread_id][i].try_dequeue(temp))
      {
        process_cnt++;
        const uint16_t &pos = temp.pos;
        sketch_sub_section[thread_id].counter[pos] += (temp.value);
      }
    }
#ifdef ONLINEQUERY
    if (cnt <= 10)
      return;
    UpdateSnapshot(thread_id, sketch_sub_section);
#endif
  }
#ifdef ONLINEQUERY
  inline void UpdateSnapshot(uint64_t thread_id,  struct GlobalSketchSubSection *sketch_sub_section)
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
    for (uint64_t j = 0; j < sub_sketch_length; j++)
    {
      this->threads_outcome[thread_id].outcome_cache.double_outcome_cache[local_round].array[j] = sketch_sub_section[thread_id].counter[j];
    }
    // reverse the writer flag to exchange the snapshot for reader to access
    this->thread_snapshot_round[thread_id].value ^= 1;
    this->query_flag[thread_id].value.fetch_xor(0x01000000);
  }
#endif
};
#endif