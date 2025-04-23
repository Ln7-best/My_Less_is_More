#ifndef OURS_H
#define OURS_H
#include <algorithm>
#include <atomic>
#include <chrono>
#include "config.h"
#include "hash.h"
#include "util.h"
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <ratio>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include "readerwriterqueue.h"
// #define ONLINEQUERY
const uint64_t PROCESSGAP = 25000;

struct alignas(64) GlobalSketchSubSection
{

  SS_bucket counter[sub_sketch_length * HASH_NUM];
  char padding[64 - (sizeof(counter) % 64)];
};
#ifdef ONLINEQUERY
template <typename T>
struct alignas(64) CacheArray
{
  T array[HASH_NUM * sub_sketch_length];
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

    size_t sub_sketch_length_round = (sub_sketch_length + 63) / 64 * 64;
    for (int i = 0; i < NUM_OUTCOME; i++)
    {
      outcome[i] = static_cast<SS_bucket *>(aligned_alloc(64, sizeof(SS_bucket) * HASH_NUM * sub_sketch_length_round));
      if (!outcome[i])
      {
        std::cerr << "Memory allocation failed for outcome[" << i << "]\n";
        std::exit(EXIT_FAILURE);
      }
      volatile char *touch = reinterpret_cast<volatile char *>(outcome[i]);
      for (size_t j = 0; j < sizeof(uint64_t) * HASH_NUM * sub_sketch_length_round; j += 4096)
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
  typedef std::unordered_map<uint64_t, int64_t> HashMap;
  void update(void *start, uint64_t size)
  {
    size = size;
    std::thread parent;
    parent = std::thread(&Ours::ParentThread, this, &parent, start, size);
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
  QueryOutcome<SS_bucket> threads_outcome[thread_num];
#endif
  struct  GlobalSketchSubSection global_sketch[thread_num];
  std::atomic<uint32_t> partition_num;
  uint32_t dataset_size[thread_num];
  Paddedint operation_cnt[thread_num];
  typedef ReaderWriterQueue<SS_Entry<Key>> myQueue;
  myQueue que[thread_num][thread_num];
  MyChild_SS<Key> *initialize_child() { return new MyChild_SS<Key>(); }

  static void Partition(Key *start, uint64_t size, uint32_t id, std::vector<Key> &vec)
  {
    uint64_t interval = size / thread_num;
    for (uint64_t i = interval * id; i < interval * (id + 1); i++)
    {
      vec.push_back(start[i]);
    }
  }
  HashMap GetSSCandidates()
  {
    HashMap ret;
    uint32_t card[HASH_NUM][LENGTH];
    uint64_t cnt = 0;
    std::unordered_set<uint32_t> srcs;
    for (uint64_t hashpos = 0; hashpos < HASH_NUM; hashpos++)
    {
      for (uint64_t pos = 0; pos < LENGTH; pos++)
      {
        uint64_t push_sec_id = pos / sub_sketch_length;
        uint64_t sec_inner_index = pos % sub_sketch_length;
        double sum = 0;
        uint32_t counter_pos = sec_inner_index + hashpos * sub_sketch_length;
        for (uint32_t i = 0; i < HLL_LEN; i++)
        {
          sum += (1.0 / ((unsigned long long)1 << (unsigned long long)global_sketch[push_sec_id].counter[counter_pos].hll_rank[i]));
          if (global_sketch[push_sec_id].counter[counter_pos].hll_rank[i] > 20)
          {
            cnt++;
          }
        }
        card[hashpos][pos] = 0.7213 * HLL_LEN / (HLL_LEN + 1.079) * HLL_LEN / sum * HLL_LEN;
        srcs.insert(global_sketch[push_sec_id].counter[counter_pos].ss_candidate);
      }
    }
    for (auto src : srcs)
    {
      uint32_t minimum = 0xffffffff;
      for (uint64_t hashPos = 0; hashPos < HASH_NUM; hashPos++)
      {
        uint32_t pos = hash(src, hashPos) % LENGTH;
        minimum = std::min(minimum, card[hashPos][pos]);
      }
      ret[src] = minimum;
    }

    return ret;
  }

  void ParentThread(std::thread *thisThd, void *start, uint64_t size)
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
      for (uint64_t j = 0; j < HASH_NUM * sub_sketch_length; j++)
      {
        global_sketch[i].counter[j].level = 0;
        global_sketch[i].counter[j].ss_candidate = 0;
        for (uint64_t k = 0; k < HLL_LEN; k++)
        {
          global_sketch[i].counter[j].hll_rank[k] = 0;
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
    uint64_t tot_keys = 0;
    for (uint32_t i = 0; i < thread_num; ++i)
    {
      thd[i].join();
      tot_keys += dataset_size[i];
    }
    std::cout << "Insert " << tot_keys << " keys." << std::endl;
  }

  void ChildThread(std::thread *thisThd, uint32_t thread_id, void *start,
                   uint64_t size, std::atomic<int32_t> *finish)
  {
#ifdef __linux__
    if (!setaffinity(thisThd, thread_id))
      return;
#endif
    MyChild_SS<Key> *sketch = initialize_child();
    uint64_t local_min = 0;
    std::vector<Key> dataset;
    Partition((Key *)start, size / sizeof(Key), thread_id,
              dataset);
    partition_num.fetch_add(1);
    while (partition_num < thread_num)
    {
    }
    insert(dataset, sketch, que, thread_id);
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
        std::memcpy(
            threads_outcome[i].outcome[round],
            this->threads_outcome[i].outcome_cache.double_outcome_cache[idx].array,
            sizeof(SS_bucket) * HASH_NUM * sub_sketch_length);
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
  void insert(const std::vector<Key> &dataset, MyChild_SS<Key> *sketch,
              myQueue queue_group[][thread_num], uint32_t thread_id)
  {
    uint32_t length = dataset.size();
    for (uint32_t i = 0; i < length; ++i)
    {
      insert_child(sketch, queue_group, dataset[i], thread_id, global_sketch);
      if (i % PROCESSGAP == 0)
        process_queue(queue_group, global_sketch, thread_id);
#ifdef ONLINEQUERY
      if ((i + 1) % 50 == 0)
        real_value_for_query[thread_id].value.store(i + 1, std::memory_order_relaxed);
#endif
    }
  }

  void insert_child(MyChild_SS<Key> *p, myQueue q_group[][thread_num],
                    const Key &packet,
                    uint64_t thread_id,
                    struct GlobalSketchSubSection *sketch_sub_section)
  {

    auto sketch = ((MyChild_SS<Key> *)p)->sketch;
    uint32_t pos[HASH_NUM];
    uint32_t src = packet >> 32;
    uint32_t temp = hash<uint64_t>(packet, 101);
    uint8_t rank = MIN(31, __builtin_clz(temp) + 1);
    uint32_t hll_pos = temp % HLL_LEN;
    for (uint32_t hashPos = 0; hashPos < HASH_NUM; ++hashPos)
    {
      pos[hashPos] = hash(src, hashPos) % LENGTH;
    }
    for (uint32_t hashPos = 0; hashPos < HASH_NUM; ++hashPos)
    {
      bool modify = false;
      if (sketch[hashPos][pos[hashPos]].level < rank)
      {
        sketch[hashPos][pos[hashPos]].level = rank;
        sketch[hashPos][pos[hashPos]].ss_candidate = src;
        modify = true;
      }
      if (sketch[hashPos][pos[hashPos]].hll_rank[hll_pos] < rank)
      {
        sketch[hashPos][pos[hashPos]].hll_rank[hll_pos] = rank;
        modify = true;
      }
      // due to the unique properties of HLL, buckets with large distinct counts are pushed less frequently
      // therefore, batching updates with small ranks is sufficient
      // there's no need to use a stash
      if (modify && rank > 4)
      {
        uint64_t push_sec_id = pos[hashPos] / sub_sketch_length;
        uint16_t sec_inner_index = pos[hashPos] % sub_sketch_length;
        while (__builtin_expect(!q_group[push_sec_id][thread_id].enqueue_fast(
                                    SS_Entry<Key>(src, hashPos, sec_inner_index, temp)),
                                0))
        {
          process_queue(q_group, sketch_sub_section, thread_id);
        }
      }
    }
  }

  inline void
  process_queue(myQueue q_group[][thread_num],
                struct GlobalSketchSubSection *sketch_sub_section,
                uint64_t thread_id = 0)
  {
    SS_Entry<Key> temp;
    uint64_t cnt = 0;
    for (uint64_t i = 0; i < thread_num; i++)
    {
      while (q_group[thread_id][i].try_dequeue(temp))
      {
        cnt++;
        const uint32_t &hashpos = temp.hashPos;
        const uint32_t &pos = temp.pos;
        uint8_t rank = MIN(31, __builtin_clz(temp.value) + 1);
        uint32_t hll_pos = temp.value % HLL_LEN;
        uint64_t counter_pos = hashpos * sub_sketch_length + pos;
        if (sketch_sub_section[thread_id].counter[counter_pos].level < rank)
        {
          sketch_sub_section[thread_id].counter[counter_pos].level = rank;
          sketch_sub_section[thread_id].counter[counter_pos].ss_candidate = temp.key;
        }
        if (sketch_sub_section[thread_id].counter[counter_pos].hll_rank[hll_pos] < rank)
        {
          sketch_sub_section[thread_id].counter[counter_pos].hll_rank[hll_pos] = rank;
        }
      }
    }
#ifdef ONLINEQUERY
    // update the snapshot for the query thread
    if (cnt <= 10)
      return;
    UpdateSnapshot(thread_id, sketch_sub_section);
#endif
  }
#ifdef ONLINEQUERY
  inline void UpdateSnapshot(uint64_t thread_id, struct  GlobalSketchSubSection *sketch_sub_section)
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
    std::memcpy(
        this->threads_outcome[thread_id].outcome_cache.double_outcome_cache[local_round].array,
        sketch_sub_section[thread_id].counter,
        sizeof(SS_bucket) * HASH_NUM * sub_sketch_length);
    // reverse the writer flag to exchange the snapshot for reader to access
    this->thread_snapshot_round[thread_id].value ^= 1;
    this->query_flag[thread_id].value.fetch_xor(0x01000000);
  }
#endif
};

#endif