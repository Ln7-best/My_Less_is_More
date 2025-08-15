#ifndef BENCHMARK_H
#define BENCHMARK_H
#include <cstdint>
#include <unordered_map>
#include <vector>
#include "loader.h"
#include "Ours.h"

template <typename Key>
class Benchmark
{
public:
  /**
   * Load a binary file (dataset) and get the statistics
   */
  Benchmark(const char *PATH)
  {
    result = Load(PATH);
    uint64_t length = result.length / sizeof(Key);
    Key *dataset = (Key *)result.start;
    // for (uint64_t i = 0; i < length; ++i) {
    //   mp[dataset[i]] += 1;
    // }
  }

  ~Benchmark() { UnLoad(result); }

  void Bench()
  {
    auto alg = new Ours<Key, THREAD_NUM>;
    alg->Update(result.start, result.length ,mp);
    delete alg;
  }

private:
  LoadResult result;
  std::unordered_map<Key,int32_t> mp;
};

#endif
