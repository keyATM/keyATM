#include "sampler.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

namespace sampler{

  double slice_uniform(const double lower, const double upper)
  {
    return lower + (upper - lower) * R::unif_rand();
  }


  std::vector<int> shuffled_indexes(const int m)
  {
    // Returns a vector of shuffled indexes for sampling
    std::vector<int> v(m);
    std::iota(v.begin(), v.end(), 0);
    std::random_shuffle(v.begin(), v.end(), sampler::rand_wrapper);
    return v;
  }


  int rcat(Eigen::VectorXd &prob, const int size)
  { 
    double u = R::unif_rand();
    double temp = 0.0;
    int index = 0;

    for (int ii = 0; ii < size; ii++) {
      temp += prob(ii);
      if (u < temp) {
        index = ii;
        break;
      }
    }
    return index;
  }


  int rcat_without_normalize(Eigen::VectorXd &prob, const double total, const int size)
  { 
    // Draw from a categorial distribution
    // This function does not requiare a normalized probability vector.
    double u = R::unif_rand() * total;
    double temp = 0.0;
    int index = 0;
    for (int ii = 0; ii < size; ii++) {
      temp += prob(ii);
      if (u < temp) {
        index = ii;
        break;
      }
    }
    return index;
  }


  int rcat_eqsize(const int size)
  {
    double u = R::unif_rand(); 
    double temp = 0.0;
    int index = 0;
    double prob = 1.0 / size;

    for (int ii = 0; ii < size; ii++) {
      temp += prob;
      if (u < temp) {
        index = ii;
        break;
      }
    }
    return index;
  }


  int rcat_eqprob(const double prob, const int size)
  {
    double u = R::unif_rand(); 
    double temp = 0.0;
    int index = 0;

    for (int ii = 0; ii < size; ii++) {
      temp += prob;
      if (u < temp) {
        index = ii;
        break;
      }
    }
    return index; 
  }
}
