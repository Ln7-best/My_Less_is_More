#ifndef CM_H
#define CM_H
#include <cstdint>
#define SKIPHH
#include "sketch.h"
#include <chrono>
// #define QTIMESTAMP

typedef uint8_t Value;

template <typename Key> struct CM_Entry {
  Key key;
  uint16_t hashPos;
  uint16_t pos;
  uint16_t value;
  // char padding[3];

  CM_Entry(Key _key = 0, uint64_t _hashPos = 0, uint64_t _pos = 0,
           uint64_t _value = 0)
      : key(_key), hashPos(_hashPos), pos(_pos), value(_value){};
};

template <typename Key, uint32_t _HASH_NUM, uint32_t _LENGTH>
class Child_CM : public Sketch<Key> {
public:
  Value sketch[_HASH_NUM][_LENGTH];
  Child_CM() { 
    memset(sketch, 0, sizeof(Value) * _HASH_NUM * _LENGTH); 
  }

  ~Child_CM() {}

  void insert_one(const Key &packet) {}

  HashMap query_all() {
    HashMap ret;
    return ret;
  }
};
#endif