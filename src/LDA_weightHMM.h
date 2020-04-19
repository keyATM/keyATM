#ifndef __keyATM_weightHMM__INCLUDED__
#define __keyATM_weightHMM__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "LDA_base.h"
#include "keyATM_HMM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class LDAhmm : public LDAbase, public keyATMhmm
{
  public:
    // Constructor
    LDAhmm(List model_, const int iter_) :
      keyATMmeta(model_, iter_),
      LDAbase(model_, iter_),
      keyATMhmm(model_, iter_) {};

    // Functions
    void iteration_single(int it) final;
    double loglik_total() final;

};

#endif


