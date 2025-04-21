#ifndef BENCHMARK_H
#define BENCHMARK_H
#define MEASUREAVGTHROUGHPUT
// #define BENCHTIMES 2
#include <cstdint>
#include <vector>

#include "loader.h"

#include "Ours.h"

template <typename Key>
class Benchmark
{
public:
  typedef std::unordered_map<Key, int32_t> HashMap;

  /**
   * Load a binary file (dataset) and get the statistics
   */
  Benchmark(const char *PATH)
  {
    result = Load(PATH);
    uint64_t length = result.length / sizeof(Key);
    Key *dataset = (Key *)result.start;
  }

  ~Benchmark() { UnLoad(result); }

  void Bench()
  {
    auto alg = new Ours<Key, THREAD_NUM>;
    alg->Update(result.start, result.length);
    delete alg;
  }

private:
  LoadResult result;
  HashMap mp;
};

#endif
