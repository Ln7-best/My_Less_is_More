#ifndef CM_H
#define CM_H
#include <cstdint>
#include <ratio>
#include <chrono>
#include <cstring>
// #define QTIMESTAMP
#define HLL_LEN 63
typedef uint8_t Value;
template <typename Key> struct SS_Entry {
  uint32_t key;
  uint16_t hashPos;
  uint16_t pos;
  uint32_t value;
  SS_Entry(uint32_t _key = 0, uint64_t _hashPos = 0, uint64_t _pos = 0,
           uint64_t _value = 0)
      : key(_key), hashPos(_hashPos), pos(_pos), value(_value){};
};
struct SS_bucket{
  uint8_t hll_rank[HLL_LEN];
  uint32_t ss_candidate;
  uint8_t level; 
};  
template <typename Key, uint32_t _HASH_NUM, uint32_t _LENGTH>
class Child_SS {
public:
  SS_bucket sketch[_HASH_NUM][_LENGTH];
  
  Child_SS() { memset(sketch, 0, sizeof(SS_bucket) * _HASH_NUM * _LENGTH);}

  ~Child_SS() {}

};

#endif