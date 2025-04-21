#include "../common/util.h"
#include "benchmark.h"
#include <cstdlib>
#include <pthread.h>
#include <string>
int main(int argc, char *argv[]) {
  std::ios::sync_with_stdio(false);
  if (argc < 3) {
    std::cerr << "File name needed" << std::endl;
    return -1;
  }
  Benchmark<uint64_t> bench(argv[1],std::stoll(argv[2]) );
  bench.Bench();
  return 0;
}
