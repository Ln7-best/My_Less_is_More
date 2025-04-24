#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <vector>

#include "loader.h"
#include "Ours.h"
#include <unordered_map>

template <typename Key>
class Benchmark
{
public:
    Benchmark(const char *PATH)
    {
        result = Load(PATH);
        uint64_t length = result.length / sizeof(Key);
        Key *dataset = (Key *)result.start;
    }
    ~Benchmark()
    {
        UnLoad(result);
    }

    void Bench()
    {
        auto alg = new Ours<Key, THREAD_NUM>;
        alg->Update(result.start, result.length);
        delete alg;
    }

private:
    LoadResult result;
};

#endif
