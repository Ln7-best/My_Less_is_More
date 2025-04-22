#ifndef UNIVMON_H
#define UNIVMON_H

#include <cstdint>
#include <cstring>
typedef int8_t Value;

constexpr static int32_t increment[2] = {1, -1};

template<typename Key>
struct Univ_Entry{
    Key key;
    uint16_t level;
    uint16_t hashPos;
    uint16_t pos;
    int16_t value;

    Univ_Entry(Key _key = 0, uint16_t _level = 0, uint16_t _hashPos = 0,
            uint16_t _pos = 0, int16_t _value = 0):
            key(_key), level(_level), hashPos(_hashPos), pos(_pos), value(_value){};
};

template<typename Key, uint32_t MAX_LEVEL, uint32_t HASH_NUM, uint32_t LENGTH>
class Child_UnivMon {
public:

    Value sketch[MAX_LEVEL][HASH_NUM][LENGTH];

    Child_UnivMon(){
        memset(sketch, 0, sizeof(Value) * MAX_LEVEL * HASH_NUM * LENGTH);
    }

    ~Child_UnivMon(){}

};
#endif