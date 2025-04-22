#ifndef BENCHMARK_H
#define BENCHMARK_H


#include <vector>

#include "loader.h"

// #include "Ideal.h"
// #include "Merge.h"
#include "Ours.h"
#include <unordered_map>

// #define THREAD_NUM 16

template<typename Key>
class Benchmark{
public:
    typedef std::unordered_map<Key, uint64_t> HashMap;

    Benchmark(const char *PATH, uint64_t benchtimes){
        result = Load(PATH);
        this->benchtimes = benchtimes;
        uint64_t length = result.length  / sizeof(Key);
        Key* dataset = (Key*)result.start;

        // mp =calculateQuantiles(dataset,length);
    }
    HashMap calculateQuantiles(const Key* data, uint64_t length) {
        std::vector<double> quantiles = {0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875,0.999};
        HashMap quantileMap;
        if (length == 0) {
            return quantileMap; // 空数据直接返回
        }
    
        // 复制数据并排序（避免修改原数据）
        std::vector<Key> sortedData(data, data + length);
        std::sort(sortedData.begin(), sortedData.end());
    
        // quantileMap[0] = sortedData[0];
        // for (int k = 1; k < 1000; ++k) {
        //     // 计算位置（0-based）
        //     double pos = (length - 1) * k / 1000.0;
        //     uint64_t index = static_cast<uint64_t>(pos); // 直接取整
        //     quantileMap[k] = sortedData[index];
        // }
        for (int k = 0; k < 8; ++k) {
            // 计算位置（0-based）
            double pos = (length - 1)* quantiles[k];
            uint64_t index = static_cast<uint64_t>(pos); // 直接取整
            quantileMap[k] = sortedData[index];
        }
        return quantileMap;
    }
    ~Benchmark(){
    	UnLoad(result);
    }

    void Bench(){
        double tot_throughput = 0;
        for (uint64_t i = 0; i < benchtimes; i++) {
            double throughput;
            auto alg = new Ours<Key, THREAD_NUM>;
            alg->update(result.start, result.length, mp, &throughput);
            delete alg;
    }
}
    
private:
    uint64_t benchtimes;
    LoadResult result;
    HashMap mp;
};

#endif
