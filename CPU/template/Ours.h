#ifndef OURS_H
#define OURS_H
#include "heap.h"
#include <algorithm>
#include <atomic>
#include <chrono>
// mark
//  #define MEASUREACC
//  #include "config.h"
//  #include "hash.h"
#include "config.h"
#include "sketch.h"
#include <barrier>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <ratio>
#include <vector>

#define MEASURETIME
#include "Abstract.h"
const uint64_t PROCESSGAP = 100000;
const uint64_t COUTERGAP = PROCESSGAP / 100;
std::atomic<uint32_t> partition_num;
template <typename T> void time_process(std::vector<T> time_vector) {
  T max_time = 0;
  T min_time = 1e9;
  T tot_time = 0;
  for (auto time : time_vector) {
    max_time = std::max(max_time, time);
    min_time = std::min(min_time, time);
    tot_time += time;
  }
  std::cout << "max: " << max_time << " min: " << min_time
            << " avg: " << (double)tot_time / time_vector.size()
            << " tot: " << tot_time << "size: " << time_vector.size()
            << std::endl;
}

struct alignas(64) global_sketch_sub_section {
  std::atomic<uint64_t> counter[HASH_NUM * sub_sketch_length];
  char padding[64 - (sizeof(counter) % 64)];
};

template <typename Key> struct Bucket {
  uint64_t vote;
  Key ID[COUNTER_PER_BUCKET];
  uint64_t count[COUNTER_PER_BUCKET];
};


template <typename Key> struct Bucket_ {
  uint64_t vote;
  Key ID[COUNTER_PER_BUCKET_FILTER];
  uint64_t count[COUNTER_PER_BUCKET_FILTER];
  uint16_t pos[COUNTER_PER_BUCKET_FILTER];
};



template <typename Key> struct alignas(64) global_buckets_sub_section {
  Bucket<Key> buckets[BUCKET_LENGTH];
  char padding[64 - ((sizeof(Bucket<Key>) * BUCKET_LENGTH) % 64)];
};

template<typename Key> struct alignas(64) child_buckets_sub_section{
  Bucket_<Key> buckets[HASH_NUM][FILTER_BUCKET_LENGTH];
  char padding[64 - sizeof(buckets) % 64];
};
// global_buckets_sub_section<Key> global_buckets[thread_num];
// child_buckets_sub_section<Key> child_filters[thread_num];
template <typename Key, typename Entry, uint32_t thread_num>
class Ours : public Abstract {
private:
  struct global_sketch_sub_section global_sketch[thread_num];
  double running_time[thread_num];
  uint32_t dataset_size[thread_num];
public:
  Paddedatomic process_counter[thread_num];
  Paddedint operation_cnt[thread_num];
  global_buckets_sub_section<Key> global_buckets[thread_num];
  child_buckets_sub_section<Key> child_filters[thread_num];
  Heap<Key, int32_t> *heap;
  std::unordered_map<Key, Entry> thread_map[3 * thread_num];
  typedef ReaderWriterQueue<Entry> myQueue;
  struct Paddedint global_cnt[thread_num];
  struct Paddedint round[thread_num];
  myQueue que[thread_num][thread_num];
  virtual Sketch<Key> *initialize_child() = 0;
  virtual void
  insert_child(Sketch<Key> *sketch, myQueue q_group[][thread_num],
               const Key &packet,
               uint64_t thread_id,
               struct global_sketch_sub_section *sketch_sub_section) = 0;
  virtual void
  process_queue(myQueue q_group[][thread_num],
                struct global_sketch_sub_section *sketch_sub_section,
                uint64_t thread_id) = 0;
  void update(void *start, uint64_t size, HashMap mp,
              double *throughput = nullptr) {
    size = size;
    std::thread parent;
    parent = std::thread(&Ours::ParentThread, this, &parent, start, size, &mp,
                         throughput);
    parent.join();
  }


  HashMap query_all(){
    HashMap ret;
    for(uint32_t i=0;i<thread_num;i++)
    {
      for(uint32_t j=0;j<BUCKET_LENGTH;j++)
      {
        for(uint32_t k=0;k<COUNTER_PER_BUCKET;k++)
        {
          if(global_buckets[i].buckets[j].count[k]==0)
            continue;
          ret[global_buckets[i].buckets[j].ID[k]]=global_buckets[i].buckets[j].count[k];
        }
      }
    }
    return ret;
  }

  void ParentThread(std::thread *thisThd, void *start, uint64_t size,
                    HashMap *mp, double *throughput) {
#ifdef __linux__
    if (!setaffinity(thisThd, thread_num))
      return;
#endif
    std::atomic<int32_t> finish(0);
    std::thread thd[thread_num];
    partition_num.store(0);
    for (uint64_t i = 0; i < thread_num; i++) {
      for (uint64_t j = 0; j < HASH_NUM * sub_sketch_length; j++) {
        global_sketch[i].counter[j] = 0;
      }
    }
    for (size_t t = 0; t < thread_num; ++t) {
      for (size_t h = 0; h < HASH_NUM; ++h) {
          for (size_t b = 0; b < FILTER_BUCKET_LENGTH; ++b) {
              child_filters[t].buckets[h][b].vote = 0;
              for(uint64_t i=0;i<COUNTER_PER_BUCKET_FILTER;i++)
              {
                child_filters[t].buckets[h][b].ID[i] = 0;
                child_filters[t].buckets[h][b].count[i] = 0;
                child_filters[t].buckets[h][b].pos[i] = 0;
              }
          }
      }
  }
    for (uint64_t i = 0; i < thread_num; i++) {
      for (uint64_t j = 0; j < BUCKET_LENGTH; j++) {
        global_buckets[i].buckets[j].vote = 0;
        for (uint64_t k = 0; k < COUNTER_PER_BUCKET; k++) {
          global_buckets[i].buckets[j].ID[k] = 0;
          global_buckets[i].buckets[j].count[k] = 0;
        }
      }
    }
  
    for (uint32_t i = 0; i < thread_num; ++i) {
      operation_cnt[i].value = 0;
      global_cnt[i].value = 0;
      round[i].value = 0;
      thd[i] = std::thread(&Ours::ChildThread, this, &(thd[i]), i, start, size,
                           &finish, global_sketch + i);
    }
    while (partition_num < thread_num) {
    }
    for (uint32_t i = 0; i < thread_num; ++i) {
      thd[i].join();
    }
    // collect(heap, finish);
    while (finish.load(std::memory_order_seq_cst) < thread_num) {
    }
    double avg_time = 0;
    double max_time = 0;
    double avg_numa1 = 0;
    uint64_t tot_size = 0;
    double avg_numa2 = 0;
    for (uint32_t i = 0; i < thread_num; i++) {
      avg_time += running_time[i];
      max_time = std::max(max_time, running_time[i]);
      tot_size += dataset_size[i];
      std::cout << "thread id: " << i << " dataset size: " << dataset_size[i]
                << " running time: " << running_time[i]<<std::endl;
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
                   struct global_sketch_sub_section *sketch_sub_section) {
#ifdef __linux__
    if (!setaffinity(thisThd, thread_id))
      return;
#endif
    Sketch<Key> *sketch = initialize_child();
    std::vector<Key> dataset;
    // uint64_t round = 0;
    Partition<Key, thread_num>((Key *)start, size / sizeof(Key), thread_id,
                               dataset);
    partition_num.fetch_add(1);
    while (partition_num < thread_num) {
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
    (*finish).fetch_add(1,std::memory_order_seq_cst);
    while ((*finish).load(std::memory_order_seq_cst) < thread_num) {
      process_queue(que, global_sketch, thread_id);

    }
    delete sketch;
  }

  void insert(const std::vector<Key> &dataset, Sketch<Key> *sketch,
                      myQueue queue_group[][thread_num], uint32_t thread_id,
                      struct global_sketch_sub_section *sketch_sub_section) {
    uint32_t length = dataset.size();
    for (uint32_t i = 0; i < length; ++i) {
      insert_child(sketch, queue_group, dataset[i], thread_id, global_sketch);
      if (i % COUTERGAP == 0)
        process_counter[thread_id].value.store(i + 1,std::memory_order_relaxed);
      if (i % PROCESSGAP == 0)
        process_queue(queue_group, global_sketch, thread_id);
    }
    process_counter[thread_id].value.store(length);
  }


};
#endif