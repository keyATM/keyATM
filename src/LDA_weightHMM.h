#ifndef __keyATM_weightHMM__INCLUDED__
#define __keyATM_weightHMM__INCLUDED__
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#include "LDA_base.h"
#include "keyATM_HMM.h"
#include "sampler.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class LDAhmm : public LDAbase, public keyATMhmm {
public:
  // Constructor
  LDAhmm(List model_)
      : keyATMmeta(model_), LDAbase(model_), keyATMhmm(model_) {};

  // Functions
  void iteration_single(int it) override final;
  double loglik_total() override final;
};

#endif
