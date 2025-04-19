#ifndef CM_OURS_H
#define CM_OURS_H
// #define RAW
#define ATOMIC
// #define ATOMIC_OPT
#include "../template/Ours.h"
#include "CM.h"
#include "config.h"
#include <chrono>
#include <cstdint>
#include <numa.h>
#include <sys/types.h>
#include <unordered_map>
/**
 * OctoSketch for the Count-Min sketch
 */
// implementation for direct write
// template <typename Key>
// struct Hash_Entry{
//   Key key;
//   uint64_t counter;

// };
template <typename Key, uint32_t thread_num>
class CM_Ours : public Ours<Key, CM_Entry<Key>, thread_num>
{
public:
  typedef ReaderWriterQueue<CM_Entry<Key>> myQueue;
  typedef ReaderWriterQueue<Heap_entry<Key>> heapQueue;
  Sketch<Key> *initialize_parent() { return new MyCM<Key>(); }

  Sketch<Key> *initialize_child() { return new MyChild_CM<Key>(); }

  void insert_child(Sketch<Key> *p, myQueue q_group[][thread_num],
                    const Key &packet,
                    //  Heap<Key, int32_t>& heap,
                    //  uint64_t test,
                    uint64_t thread_id,
                    struct global_sketch_sub_section *sketch_sub_section)
  {
    double hit_rate;
    uint64_t local_min;
    auto sketch = ((MyChild_CM<Key> *)p)->sketch;
    uint32_t pos[HASH_NUM];
    for (uint32_t hashPos = 0; hashPos < HASH_NUM; ++hashPos)
    {
      pos[hashPos] = hash(packet, hashPos) % LENGTH;
    }

    for (uint32_t hashPos = 0; hashPos < HASH_NUM; ++hashPos)
    {
      sketch[hashPos][pos[hashPos]] += 1;
      if (__builtin_expect(sketch[hashPos][pos[hashPos]] >= PROMASK,0))
      // if (sketch[hashPos][pos[hashPos]] >= PROMASK)
      {
        // this->full_cnt[thread_id].value += 1;
        // uint64_t push_sec_id = pos[hashPos] / sub_sketch_length;
        // uint16_t sec_inner_index = pos[hashPos] % sub_sketch_length;
        // while( __builtin_expect(!q_group[push_sec_id][thread_id].enqueue_fast(CM_Entry<Key>(
        //   packet, hashPos, sec_inner_index, sketch[hashPos][pos[hashPos]])), 0))
        // {
        //   // this->full_cnt[thread_id].value += 1;
        //   process_queue(q_group, sketch_sub_section, this->heap, &(this->heap_que[thread_id]), thread_id);
        // }
        // sketch[hashPos][pos[hashPos]] = 0;
       // --------------------------------------------------------------------
        bool update_flag = false;
        uint64_t minVal = 0x7fffffff;
        uint32_t minPos = 0;
        // uint64_t temp_val = 1;
        sketch[hashPos][pos[hashPos]] = 0;
        uint64_t idx = hash(pos[hashPos], 101) % FILTER_BUCKET_LENGTH;
        for (uint32_t j = 0; j < COUNTER_PER_BUCKET_FILTER; j++)
        {
          if (this->child_filters[thread_id].buckets[hashPos][idx].ID[j] ==
              packet)
          {
            this->child_filters[thread_id].buckets[hashPos][idx].count[j]++;
            if (this->child_filters[thread_id].buckets[hashPos][idx].count[j] == 9)
            // if (__builtin_expect(this->child_filters[thread_id].buckets[hashPos][idx].count[j] == 9,0))
            {
              // this->full_cnt[thread_id].value += 1;
              uint64_t push_sec_id = pos[hashPos] / sub_sketch_length;
              uint16_t sec_inner_index = pos[hashPos] % sub_sketch_length;
              while (__builtin_expect(!q_group[push_sec_id][thread_id].enqueue_fast(
                  CM_Entry<Key>(packet, hashPos, sec_inner_index,
                                9*PROMASK)),0))
              // while (!q_group[push_sec_id][thread_id].enqueue_fast(
              //   CM_Entry<Key>(packet, hashPos, sec_inner_index,
              //                 this->child_filters[thread_id]
              //                     .buckets[hashPos][idx]
              //                     .count[j])))
              {
                // this->full_cnt[thread_id].value += 1;
                // measure_time([&](){
                process_queue(q_group, sketch_sub_section, this->heap,
                              &(this->heap_que[thread_id]), thread_id);
                // },this->merge_time[thread_id]);

              }
              // q_group[push_sec_id][thread_id].enqueue_fast(
              //   CM_Entry<Key>(packet, hashPos, sec_inner_index,
              //                 this->child_filters[thread_id]
              //                     .buckets[hashPos][idx]
              //                     .count[j]));
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
        if(update_flag)
          continue;
        uint64_t insert_value = 1;
        uint64_t insert_pos = pos[hashPos];
        uint64_t insert_key = packet;
        // if (__builtin_expect((this->child_filters[thread_id].buckets[hashPos][idx].vote + 1) >=
        //     minVal * 8,0))
        if ((this->child_filters[thread_id].buckets[hashPos][idx].vote + 1) >= minVal * 8 || insert_pos == this->child_filters[thread_id].buckets[hashPos][idx].pos[minPos])
        {
          // this->full_cnt[thread_id].value += 1;
          insert_value = minVal;
          insert_key = this->child_filters[thread_id].buckets[hashPos][idx].ID[minPos];
          insert_pos = this->child_filters[thread_id].buckets[hashPos][idx].pos[minPos];
          this->child_filters[thread_id].buckets[hashPos][idx].vote = 0;
          this->child_filters[thread_id].buckets[hashPos][idx].ID[minPos] =
              packet;
          this->child_filters[thread_id].buckets[hashPos][idx].count[minPos] =
              1;
          this->child_filters[thread_id].buckets[hashPos][idx].pos[minPos] =
              pos[hashPos];
        }
        else
        {
          // this->full_cnt[thread_id].value += 1;
          this->child_filters[thread_id].buckets[hashPos][idx].vote++;
        }
        if (insert_value == 0)
          continue;
        uint64_t push_sec_id = insert_pos / sub_sketch_length;
        uint16_t sec_inner_index = insert_pos % sub_sketch_length;
        // if(push_sec_id >= thread_num)
        // {
        //   std::cout << pos[hashPos]<<std::endl;
        //   std::cout << sub_sketch_length;
        //   std::cout <<"access error"<<std::endl;
        // }
        
        while (__builtin_expect(!q_group[push_sec_id][thread_id].enqueue_fast(CM_Entry<Key>(
            insert_key, hashPos, sec_inner_index, insert_value*PROMASK)),0))
        // while (!q_group[push_sec_id][thread_id].enqueue_fast(CM_Entry<Key>(
        //   insert_key, hashPos, sec_inner_index, insert_value)))
        {
          // this->full_cnt[thread_id].value += 1;
          // this->full_cnt[thread_id].value += 1;
          // measure_time([&](){
          process_queue(q_group, sketch_sub_section, this->heap,
                        &(this->heap_que[thread_id]), thread_id);
          // },this->merge_time[thread_id]);

        }
        // // q_group[push_sec_id][thread_id].enqueue_fast(CM_Entry<Key>(
        //   insert_key, hashPos, sec_inner_index, insert_value));
      }
    }
    // process_queue(q_group, sketch_sub_section, thread_id);
  }

  // inline void direct_write(struct global_sketch_sub_section
  // *sketch_sub_section, uint16_t&  )

  inline void process_queue(myQueue q_group[][thread_num],
                     struct global_sketch_sub_section *sketch_sub_section,
                     Heap<Key, int32_t> *heap, heapQueue *heap_que,
                     //  uint64_t& local_min ,uint64_t& round,
                     uint64_t thread_id = 0)
  {
    uint64_t process_cnt = 0;
    CM_Entry<Key> temp;
    uint64_t que_cnt = 0;
    uint64_t candidate_cnt = 0;
    uint64_t cnt = 0;
    // uint64_t global_cnt = 0;
    // std::cout<<"yes"<<std::endl;
    // uint64_t local_min = 0;
    for (uint64_t i = 0; i < thread_num; i++)
    {
      while (q_group[thread_id][i].try_dequeue(temp))
      {
        process_cnt++;
        // direct write to global sketch

        const uint16_t &hash_pos = temp.hashPos;
        const uint16_t &pos = temp.pos;
        uint64_t counter_pos = hash_pos * sub_sketch_length + pos;
#ifdef RAW
        sketch_sub_section[thread_id].counter[counter_pos] += temp.value;
        if (sketch_sub_section[thread_id].counter[counter_pos] >
            this->local_min[thread_id].value)
#elif defined(ATOMIC_OPT)
        sketch_sub_section[thread_id].counter[counter_pos] += temp.value;
        sketch_sub_section[thread_id].stale_time[counter_pos]++;
        if (sketch_sub_section[thread_id].stale_time[counter_pos] >= 10)
        {
          sketch_sub_section[thread_id].stale_time[counter_pos] = 0;
          sketch_sub_section[thread_id].shadow_counter[counter_pos].store(
              sketch_sub_section[thread_id].counter[counter_pos]);
        }
#else
        // this->full_cnt[thread_id].value++;
        sketch_sub_section[thread_id].counter[counter_pos].fetch_add(
            temp.value);
        uint64_t counter_val =
            sketch_sub_section[thread_id].counter[counter_pos].load();
        // uint64_t counter_val += temp.value;
        // sketch_sub_section[thread_id].counter[counter_pos].store(counter_val);
        if ((double)counter_val >
            this->global_cnt[thread_id].value * (ALPHA / 2))
#endif
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
          // else
          candidate_cnt++;
#ifdef RAW
          int32_t minimum = sketch_sub_section[thread_id].counter[counter_pos];
#elif defined(ATOMIC_OPT)
          int32_t minimum = sketch_sub_section[thread_id].counter[counter_pos];
#else
          int32_t minimum = counter_val;
#endif
          int32_t origin_min = minimum;
          // uint64_t candidate = 0;
          // uint64_t min_hash_pos = hash_pos;
          for (uint32_t tempHash = 0; tempHash < HASH_NUM; ++tempHash)
          {
            uint32_t tempPos = hash(temp.key, tempHash) % LENGTH;
            uint64_t push_sec_id = tempPos / sub_sketch_length;
            uint16_t sec_inner_index = tempPos % sub_sketch_length;
            // candidate = sketch_sub_section[push_sec_id].counter[tempHash *
            // sub_sketch_length + sec_inner_index]; if(minimum > candidate)
            // {
            //   minimum = candidate;
            //   min_hash_pos = tempHash;
            // }
#ifdef RAW
            minimum = MIN(
                minimum,
                sketch_sub_section[push_sec_id]
                    .counter[tempHash * sub_sketch_length + sec_inner_index]);
#elif defined(ATOMIC_OPT)
            minimum =
                MIN(minimum, sketch_sub_section[push_sec_id]
                                 .shadow_counter[tempHash * sub_sketch_length +
                                                 sec_inner_index]
                                 .load());
#else
            minimum =
                MIN(minimum,
                    sketch_sub_section[push_sec_id]
                        .counter[tempHash * sub_sketch_length + sec_inner_index]
                        .load());
#endif
          }
          // }
          // this->my_minimum[thread_id].value = minimum;
          // if(min_hash_pos == hash_pos)
          //   heap_que->enqueue(Heap_entry<Key>(temp.key, minimum));

          // if (__builtin_expect(minimum == origin_min, 1))
          // {
          que_cnt++;
          int32_t minVal = 0x7fffffff;
          uint32_t minPos = 0;
          // uint64_t idx =
          //     (sub_sketch_length * thread_id + pos) % BUCKET_LENGTH;
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
          // }
          cnt++;
        }
      }
    }
    this->round[thread_id].value++;
    // this->process_cnt[thread_id].push_back(process_cnt);
    // if(candidate_cnt!=0)
    // this->process_cnt[thread_id].push_back((double)que_cnt/candidate_cnt);
  }
  void merge(Sketch<Key> *p, CM_Entry<Key> temp)
  {
    ((MyCM<Key> *)p)->Merge(temp);
  }
};

#endif