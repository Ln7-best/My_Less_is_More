#include <cstdint>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <map>
#include <unordered_set>
#include "dirty_zipfian_int_distribution.h" // 包含你的 Zipf 分布类
#define MASK 0xffffffff00000000
// #define KV_GENERATE
std::map<uint32_t, uint64_t> my_map;
std::vector<uint64_t> output_data;
std::unordered_set<uint64_t> my_set;
int main(int argc, char *argv[])
{
    std::string outputFileName = argv[1];


    std::vector<uint64_t> data;
    uint64_t value;
    std::default_random_engine generator;
    dirtyzipf::dirty_zipfian_int_distribution<int> distribute(1, 1e7, 0.99);
    std::ofstream outFile(outputFileName, std::ios::binary);
    if (!outFile)
    {
        std::cerr << "Fail to open output file: " << outputFileName << std::endl;
        return 1;
    }

    const uint64_t queryCount = 3.2e9;
    for (uint64_t i = 0; i < queryCount; i++)
    {
        int zipfdata = distribute(generator);
        outFile.write(reinterpret_cast<const char *>(&data[zipfdata]), sizeof(uint64_t));
    }
    outFile.close();
    return 0;
}
