#ifndef COUNT_H
#define COUNT_H
#include <cstdint>
#include <cstring>
#include <chrono>
typedef int8_t Value;

constexpr static int32_t increment[2] = {1, -1};

template<typename Key>
struct Count_Entry{
    Key key;
    uint16_t hashPos;
    uint16_t pos;
    int16_t value;

    Count_Entry(Key _key = 0, uint16_t _hashPos = 0, uint16_t _pos = 0, int16_t _value = 0):
            key(_key), hashPos(_hashPos), pos(_pos), value(_value){};
};

template<typename Key, uint32_t HASH_NUM, uint32_t LENGTH>
class Child_Count{
public:

    Value sketch[HASH_NUM][LENGTH];

    Child_Count(){
        memset(sketch, 0, sizeof(Value) * HASH_NUM * LENGTH);
    }

    ~Child_Count(){}
};

#endif