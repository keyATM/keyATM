#ifndef __keyATM_weightCov__INCLUDED__
#define __keyATM_weightCov__INCLUDED__
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "LDA_base.h"
#include "keyATM_cov.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class LDAcov : public LDAbase, public keyATMcov
{
  public:
    // Constructor
    LDAcov(List model_) :
      keyATMmeta(model_),
      LDAbase(model_),
      keyATMcov(model_) {};

    // Functions
    void iteration_single(int it) override final;
    double loglik_total() override final;
};

#endif

