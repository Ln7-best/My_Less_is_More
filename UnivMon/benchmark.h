#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <vector>

#include "loader.h"
#include "Ours.h"

// #define THREAD_NUM 16

template <typename Key>
class Benchmark
{
public:
    typedef std::unordered_map<Key, int32_t> HashMap;

    Benchmark(const char *PATH, uint64_t benchtimes)
    {
        result = Load(PATH);
        this->benchtimes = benchtimes;
        uint64_t length = result.length / sizeof(Key);
        Key *dataset = (Key *)result.start;

        // for(uint64_t i = 0;i < length;++i){
        //     mp[dataset[i]] += 1;
        // }
    }

    ~Benchmark()
    {
        UnLoad(result);
    }

    void Bench()
    {
        double tot_throughput = 0;
        for (uint64_t i = 0; i < benchtimes; i++)
        {
            auto alg = new Ours<Key, THREAD_NUM>;
            double throughput;
            alg->update(result.start, result.length, mp, &throughput);
            delete alg;
            tot_throughput += throughput;
        }
        std::cout << "average throughput:" << tot_throughput / benchtimes
                  << std::endl;
    }

private:
    uint64_t benchtimes;
    LoadResult result;
    HashMap mp;
};

#endif
