#ifndef CONFIG_H
#define CONFIG_H

#include "CM.h"
#include <cstddef>
#include <cmath>
#define THREAD_NUM 64
// Configuration for the sketch
#define LENGTH (1<<16)
#define HASH_NUM 3
#define COUNTER_PER_BUCKET 4
#define COUNTER_PER_BUCKET_FILTER 1
#define FILTER_BUCKET_LENGTH 64
#define BUCKET_LENGTH (1024 * 3 / COUNTER_PER_BUCKET/THREAD_NUM)
constexpr size_t ceil_div(size_t dividend, size_t divisor) {
    return static_cast<size_t>(std::ceil(static_cast<double>(dividend) / divisor));
}
const size_t sub_sketch_length = ceil_div(LENGTH, THREAD_NUM);
template<typename Key>
using MyChild_CM = Child_CM<Key, HASH_NUM, LENGTH>;
#define ALPHA 0.0002 // Threshold for finding heavy hitters
#define COUNTERMAX 0x7f
#define NUM_OUTCOME 100000
#define ARRAY_SIZE 128 // Must be a multiple of 64
#endif
