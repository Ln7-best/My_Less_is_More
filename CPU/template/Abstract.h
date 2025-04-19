#ifndef ABSTRACT_H
#define ABSTRACT_H
// #include <gperftools/profiler.h>

#include "sketch.h"
#include "readerwriterqueue.h"
#include <condition_variable>
#include <thread>
#include <chrono>
#include <unordered_map>
#include <fstream>
#include <unordered_set>
#include <random>
typedef std::unordered_map<uint64_t, int32_t> HashMap;
unsigned int Random_Generate()
{
    unsigned int x = rand();
    unsigned int h = rand();

    return x ^ ((h & 1) << 31);
}
class Abstract
{

public:
    std::atomic<uint32_t> processed_packets_num;
    virtual void update(void *start, uint64_t size, HashMap mp, double* throughput) = 0;
    virtual ~Abstract() {};

    template <typename Key, uint32_t thread_num>
    static void Partition(Key *start, uint64_t size, uint32_t id, std::vector<Key> &vec)
    {
        uint64_t interval = size / thread_num;
        for (uint64_t i = interval * id; i < interval * (id + 1); i++)
        {
            vec.push_back(start[i]);
        }
        // for (uint32_t i = 0; i < size; ++i)
        // {
        //     if (hash(start[i], 101) % thread_num == id)
        //     {
        //         vec.push_back(start[i]);
        //     }
        // }
    }
    template <typename Key>
    static void Initkeyset(Key *start, uint32_t size, std::unordered_set<Key> &keyset)
    {
        std::cout << "init key set" << std::endl;
        for (uint32_t i = 0; i < size; ++i)
        {
            keyset.insert(start[i]);
        }
        std::cout << "init end" << std::endl;
    }
    static void HHCompare(HashMap test, HashMap real, int32_t threshold, std::ofstream *outputFile = nullptr)
    {
        double estHH = 0, HH = 0, both = 0;
        double CR = 0, PR = 0, AAE = 0, ARE = 0;
        std::cout << "real size" << real.size() << std::endl;
        std::cout << "test size" << test.size() << std::endl;
        uint64_t cnt =0;
        for (auto it = test.begin(); it != test.end(); ++it)
        {            
            if (it->second > threshold)
            {
                estHH += 1;
                if (real[it->first] > threshold)
                {
                    both += 1;
                    AAE += abs(real[it->first] - it->second);
                    ARE += abs(real[it->first] - it->second) / (double)real[it->first];
                }
            }
            else
            if ((int64_t)threshold - (int64_t)it->second < 10000)
            {
                cnt++;
            }
        }
        uint64_t test_hit = 0;
        for (auto it = real.begin(); it != real.end(); ++it)
        {
            if (it->second > threshold)
            {
                HH += 1;
                if(test.find(it->first)!=test.end())
                {
                    std::cout << (int64_t)threshold - (int64_t) test[it->first]<<std::endl; 
                    test_hit++;
                }
                
            }
        }
        if (!outputFile)
            std::cout<<"real HH:"<<HH<<std::endl
                    << "correct HH: "<<both<<std::endl
                    <<"test hit: "<<test_hit<<std::endl
                    <<"num: "<<cnt<<std::endl
                     << "CR: " << both / HH << std::endl
                      << "PR: " << both / estHH << std::endl
                      << "AAE: " << AAE / both << std::endl
                      << "ARE: " << ARE / both << std::endl
                      << std::endl;
        else
            *outputFile << "CR: " << std::to_string(both / HH) << std::endl
                        << "PR: " << std::to_string(both / estHH) << std::endl
                        << "AAE: " << std::to_string(AAE / both) << std::endl
                        << "ARE: " << std::to_string(ARE / both) << std::endl;
    }
};

#endif