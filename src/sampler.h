#ifndef __sampler__INCLUDED__
#define __sampler__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Eigen;

namespace sampler
{
  // Defines sampler used in keyATM
 
  inline int rand_wrapper(const int n) { return floor(R::unif_rand() * n); }

  double slice_uniform(const double lower, const double upper);

  std::vector<int> shuffled_indexes(const int m);

  int rcat(VectorXd &prob, const int size);
  int rcat_without_normalize(VectorXd &prob, const double total, const int size);

  int rcat_eqsize(const int size);
  int rcat_eqprob(const double prob, const int size);
}

#endif
