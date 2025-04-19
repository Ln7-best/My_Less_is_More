#ifndef BENCHMARK_H
#define BENCHMARK_H
#define MEASUREAVGTHROUGHPUT
// #define BENCHTIMES 2
#include <cstdint>
#include <vector>

#include "loader.h"

#include "Ideal.h"
#include "Merge.h"
#include "Ours.h"

template <typename Key> class Benchmark {
public:
  typedef std::unordered_map<Key, int32_t> HashMap;

  /**
   * Load a binary file (dataset) and get the statistics
   */
  Benchmark(const char *PATH, uint64_t benchtimes) {
    result = Load(PATH);
    this->benchtimes = benchtimes;
    uint64_t length = result.length / sizeof(Key);
    Key *dataset = (Key *)result.start;

    // for (uint64_t i = 0; i < length; ++i) {
    //   mp[dataset[i]] += 1;
    // }
  }

  ~Benchmark() { UnLoad(result); }

  void Bench() {
    double tot_throughput = 0;
    for (uint64_t i = 0; i < benchtimes; i++) {

      std::vector<Abstract *> absVec = {
          // new CM_Ideal<Key>(),
          // new CM_Merge<Key, THREAD_NUM>(),
          new CM_Ours<Key, THREAD_NUM>(),
      };
      for (auto alg : absVec) {
        double throughput;
        std::cout<<*(Key*)result.start;
        alg->update(result.start, result.length, mp, &throughput);
        delete alg;
        tot_throughput += throughput;
      }
    }
#ifdef MEASUREAVGTHROUGHPUT
    std::cout << "average throughput:" << tot_throughput / benchtimes
              << std::endl;
#endif
  }

private:
  uint64_t benchtimes;
  LoadResult result;
  HashMap mp;
};

#endif
