#ifndef __LDA_weight__INCLUDED__
#define __LDA_weight__INCLUDED__
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#include "LDA_base.h"
#include "keyATM_base.h"
#include "sampler.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class LDAweight : public LDAbase, public keyATMbase {
public:
  //
  // Parameters
  //
  int estimate_alpha;
  int store_alpha;

  // Slice Sampling
  double start, end, previous_p, new_p, newlikelihood, slice_;
  std::vector<int> topic_ids;
  VectorXd keep_current_param;
  double store_loglik;
  double newalphallk;

  // in alpha_loglik
  MatrixXd ndk_a;

  //
  // Functions
  //

  // Constructor
  LDAweight(List model_)
      : keyATMmeta(model_), LDAbase(model_), keyATMbase(model_) {};

  // Iteration
  virtual void iteration_single(int it) override final;
  virtual double loglik_total() override final;
};

#endif
