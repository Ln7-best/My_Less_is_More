#ifndef CONFIG_H
#define CONFIG_H

#include "DD.h"
#define THREAD_NUM 64
#define HASH_NUM 3
#define LENGTH (1 << 16)
#define COUNTER_PER_BUCKET 4
#define COUNTER_PER_BUCKET_FILTER 1
#define HEAP_SIZE 1024
#define BUCKET_LENGTH (HEAP_SIZE * 3 / COUNTER_PER_BUCKET/THREAD_NUM)
#define FILTER_BUCKET_LENGTH 64
constexpr size_t ceil_div(size_t dividend, size_t divisor) {
    return static_cast<size_t>(std::ceil(static_cast<double>(dividend) / divisor));
}
const size_t sub_sketch_length = ceil_div(LENGTH, THREAD_NUM);
#define NUM_OUTCOME 100000
#define ARRAY_SIZE 1152
template<typename Key>
using MyChild_DD = Child_DD<Key, LENGTH>;

#define PROMASK 0x7f

#endif
