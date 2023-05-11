#ifndef __keyATM_covPG__INCLUDED__
#define __keyATM_covPG__INCLUDED__
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "utils.h"
#include "keyATM_meta.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMcovPG : virtual public keyATMmeta
{
  public:
    //
    // Parameters
    //
    int num_cov;
    List PG_params;
    MatrixXd theta;

    // During the sampling
      std::vector<int> topic_ids;
      std::vector<int> cov_ids;

    // Constructor
    keyATMcovPG(List model_) :
      keyATMmeta(model_) {};

    // Read data
    void read_data_specific() override final;

    // Initialization
    void initialize_specific() override final;

    // Iteration
    virtual void iteration_single(int it) override final;
    virtual void sample_parameters(int it) override final;
    void sample_PG(int it);
    int sample_z_PG(int z, int s, int w, int doc_id);
    virtual double loglik_total() override final;
};


#endif
