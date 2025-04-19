#ifndef CONFIG_H
#define CONFIG_H

#include "CM.h"
#include <cstddef>

#define THREAD_NUM 64
// Configuration for the sketch
#define LENGTH (1<<16)
#define HASH_NUM 3
#define LEN_BITS 16
#define HEAP_SIZE 1024
#define COUNTER_PER_BUCKET 4
#define COUNTER_PER_BUCKET_FILTER 1
#define FILTER_BUCKET_LENGTH 64
#define BUCKET_LENGTH (HEAP_SIZE * 3 / COUNTER_PER_BUCKET/THREAD_NUM)
// #define BUCKET_LENGTH 14
constexpr size_t ceil_div(size_t dividend, size_t divisor) {
    return static_cast<size_t>(std::ceil(static_cast<double>(dividend) / divisor));
}
const size_t sub_sketch_length = ceil_div(LENGTH, THREAD_NUM);
// const size_t sub_sketch_length=1640;
// const size_t real_sub_sketch_length = 1640;
// const size_t sub_sketch_length = 2048;
// const size_t LENGTH = sub_sketch_length * THREAD_NUM;
#define INTERVAL 20000
#define PRINTLENGTH 
template<typename Key>
using IdealCM = CM<Key, HASH_NUM, LENGTH*(THREAD_NUM/4+1), HEAP_SIZE>;

template<typename Key>
using MyCM = CM<Key, HASH_NUM, LENGTH, HEAP_SIZE>;

template<typename Key>
using MyChild_CM = Child_CM<Key, HASH_NUM, LENGTH>;
#define UPPER_LENGTH 0xa0
#define TARGET_LENGTH 0x80
#define LOWER_LENGTH 0x60
#define THP_TIME 20

#define ALPHA 0.0002 // Threshold for finding heavy hitters

#define PROMASK 0x7f

//#define QUEUELENGTH

#endif
