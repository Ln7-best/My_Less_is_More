#ifndef DD_H
#define DD_H
#include <cmath>
#include <cstdint>
#include <cstring>
#include <chrono>
#define GAMMA 1.001
// #define GAMMA_ADD 0.020202020202020204

typedef uint8_t Value;

template<typename Key>
struct DD_Entry{
    uint16_t pos;
    uint16_t value;

    DD_Entry(uint16_t _pos = 0, uint16_t _value = 0):
            pos(_pos), value(_value){};
};

template<typename Key, uint32_t LENGTH>
class Child_DD {
public:

    Value sketch[LENGTH];
    double log_gamma;
    // CubicallyInterpolatedMapping mapping;
    // double multiplier;
    Child_DD(){
        log_gamma = log(GAMMA);
        memset(sketch, 0, sizeof(Value) * LENGTH);
    }

    ~Child_DD(){}
};
#endif