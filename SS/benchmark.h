#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <vector>

#include "loader.h"
#include <unordered_map>
// #include "Ideal.h"
// #include "Merge.h"
#include "Ours.h"

// #define THREAD_NUM 16

template <typename Key>
class Benchmark
{
public:
    typedef std::unordered_map<uint64_t, int64_t> HashMap;

    Benchmark(const char *PATH, uint64_t benchtimes)
    {
        result = Load(PATH);
        this->benchtimes = benchtimes;
        uint64_t length = result.length / sizeof(Key);
        Key *dataset = (Key *)result.start;
        std::unordered_set<uint64_t> st;
        // for(uint64_t i = 0;i < length;++i){
        //     // if(st.find(dataset[i]) == st.end()){
        //     //     st.insert(dataset[i]);
        //         uint64_t src = dataset[i] >> 32;
        //     //     if(mp.find(src) == mp.end())
        //     //         mp[src] = 1;
        //     //     else
        //             mp[src] += 1;
        //     // }
        //     // std::cout<<i<<std::endl;
        //     if(i%100000000==0)
        //         std::cout<<i<<std::endl;
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
