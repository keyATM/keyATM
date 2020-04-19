#ifndef __LDA_weight__INCLUDED__
#define __LDA_weight__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM_base.h"
#include "LDA_base.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class LDAweight : public LDAbase, public keyATMbase
{
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
    LDAweight(List model_, const int iter_) :
      keyATMmeta(model_, iter_),
      LDAbase(model_, iter_),
      keyATMbase(model_, iter_) {};

    // Iteration
    void iteration_single(int it) final;
    double loglik_total() final;
};

#endif
