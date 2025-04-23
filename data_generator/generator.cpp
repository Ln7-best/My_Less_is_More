#include <cstdint>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <map>
#include <unordered_set>
#define MASK 0xffffffff00000000
// #define KV_GENERATE
std::map<uint32_t, uint64_t> my_map;
std::vector<uint64_t> output_data;
std::unordered_set<uint64_t> my_set;
struct GenRange {
    uint64_t min;
    uint64_t max;
  };
  
class Generator {
    public:
     Generator(GenRange gen_range)
         : gen_range_(gen_range), rd_(), gen_(rd_()) {}
     virtual uint64_t Next() = 0;
     virtual ~Generator() {}
   
    protected:
     inline double RandomDouble(double min = 0.0, double max = 1.0) {
       static std::uniform_real_distribution<double> uniform(min, max);
       return uniform(gen_);
     }
   
    protected:
     GenRange gen_range_;
   
     std::random_device rd_;
     std::mt19937 gen_;
   };
class ZipfGenerator : public Generator {
    public:
     ZipfGenerator(GenRange gen_range, double theta)
         : Generator(gen_range),
           theta_(theta),
           n_(gen_range.max - gen_range.min),
           dist_(nullptr) {
       Init();
     }
     ~ZipfGenerator() override {
       if (dist_) delete dist_;
     }
   
     uint64_t Next() override { return (*dist_)(gen_) + gen_range_.min; }
   
    private:
     inline double ZipfPMF(uint64_t k) { return 1.0 / pow(k, theta_); }
   
     double ComputeNormalizationConstant() {
       double normalization_constant = 0.0;
       for (uint64_t i = 1; i <= n_; ++i) {
         normalization_constant += ZipfPMF(i);
       }
       return normalization_constant;
     }
   
     void Init() {
       const auto normalization_constant = ComputeNormalizationConstant();
   
       probs_.resize(n_);
       for (uint64_t i = 1; i <= n_; ++i) {
         probs_[i - 1] = ZipfPMF(i) / normalization_constant;
       }
   
       dist_ = new std::discrete_distribution<uint64_t>(probs_.begin(), probs_.end());
     }
   
    private:
     /* zipf param */
     double theta_;
   
     /* param for std::discrete_distribution */
     uint64_t n_;
     std::vector<double> probs_;
   
     /* gen */
     std::discrete_distribution<uint64_t>* dist_;
   };
 
int main(int argc, char *argv[])
{
    std::string outputFileName = argv[1];
    uint64_t value;
    std::default_random_engine generator;
    ZipfGenerator zipf_generator({1, 10000000}, 0.99);
    std::ofstream outFile(outputFileName, std::ios::binary);
    if (!outFile)
    {
        std::cerr << "Fail to open output file: " << outputFileName << std::endl;
        return 1;
    }
    const uint64_t queryCount = 3.2e9;
    for (uint64_t i = 0; i < queryCount; i++)
    {
        int zipfdata = zipf_generator.Next();
        outFile.write(reinterpret_cast<const char *>(&zipfdata), sizeof(uint64_t));
    }
    outFile.close();
    return 0;
}
