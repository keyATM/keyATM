#ifndef __keyATM_cov__INCLUDED__
#define __keyATM_cov__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM_meta.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMcov : virtual public keyATMmeta
{
  public:  
    //
    // Parameters
    //
    MatrixXd Alpha;
    int num_cov;
    MatrixXd Lambda;
    MatrixXd C;

    int mh_use;
    double mu;
    double sigma;

    // During the sampling
      std::vector<int> topic_ids;
      std::vector<int> cov_ids;

      // Slice sampling
      double val_min;
      double val_max;
    
    //
    // Functions
    //

    // Constructor
    keyATMcov(List model_, const int iter_) :
      keyATMmeta(model_, iter_) {};

    // Read data
    void read_data_specific() final;

    // Initialization
    void initialize_specific() final;

    // Iteration
    virtual void iteration_single(int it);
    void sample_parameters(int it);
    void sample_lambda();
    void sample_lambda_mh();
    void sample_lambda_slice();
    double alpha_loglik();
    virtual double loglik_total();
    double loglik_total_label();

    double likelihood_lambda(int k, int t);
    void proposal_lambda(int k);
};


#endif


