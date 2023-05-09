#ifndef __keyATM_cov__INCLUDED__
#define __keyATM_cov__INCLUDED__
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

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
    keyATMcov(List model_) :
      keyATMmeta(model_) {};

    // Read data
    virtual void read_data_specific() override final;

    // Initialization
    virtual void initialize_specific() override final;

    // Resume
    virtual void resume_initialize_specific() override final;

    // Iteration
    virtual void iteration_single(int it) override;
    virtual void sample_parameters(int it) override final;
    void sample_lambda();
    void sample_lambda_mh();
    void sample_lambda_slice();
    double alpha_loglik();
    virtual double loglik_total() override;

    double likelihood_lambda(int k, int t);
    void proposal_lambda(int k);
};


#endif
