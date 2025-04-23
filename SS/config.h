#ifndef CONFIG_H
#define CONFIG_H
#include <cmath>
#include <cstddef>
#include "SS.h"
#define THREAD_NUM 64
#define HASH_NUM 3
#define LENGTH (1 << 12)
constexpr size_t ceil_div(size_t dividend, size_t divisor) {
    return static_cast<size_t>(std::ceil(static_cast<double>(dividend) / divisor));
}
const size_t sub_sketch_length = ceil_div(LENGTH, THREAD_NUM);



template<typename Key>
using MyChild_SS = Child_SS<Key, HASH_NUM, LENGTH>;

#define ALPHA 0.002

#define PROMASK 0x7f

#define NUM_OUTCOME 40000

//#define QUEUELENGTH

#endif
