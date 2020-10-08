#ifndef __keyATM_covPG__INCLUDED__
#define __keyATM_covPG__INCLUDED__

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
    MatrixXd Alpha;   // remove later
    int num_cov;
    MatrixXd Lambda;  // remvoe later
    MatrixXd C;   // remove later
    List PG_params;
    MatrixXd theta;

    int mh_use;  // remove later
    double mu;  // remove later
    double sigma;  // remove later

    // During the sampling
      std::vector<int> topic_ids;
      std::vector<int> cov_ids;

      // Slice sampling
      double val_min;
      double val_max;
    
    // Constructor
    keyATMcovPG(List model_, const int iter_) :
      keyATMmeta(model_, iter_) {};

    // Read data
    void read_data_specific() final;

    // Initialization
    void initialize_specific() final;

    // Iteration
    virtual void iteration_single(int it);
    void sample_parameters(int it);
    void sample_PG();
    int sample_z_PG(VectorXd &alpha, int z, int s, int w, int doc_id); 
    void sample_lambda();
    void sample_lambda_mh();
    void sample_lambda_slice();
    double alpha_loglik();
    virtual double loglik_total();

    double likelihood_lambda(int k, int t);
    void proposal_lambda(int k);
};


#endif


