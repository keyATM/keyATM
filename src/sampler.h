#ifndef __sampler__INCLUDED__
#define __sampler__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Eigen;

namespace sampler
{
  // Defines sampler used in keyATM
 
  inline int rand_wrapper(const int n) { return floor(R::unif_rand() * n); }

  double slice_uniform(double& lower, double& upper);

  std::vector<int> shuffled_indexes(int m);

  int rcat(VectorXd &prob, int &size);
  int rcat_without_normalize(VectorXd &prob, double &total, int &size);

}

#endif
