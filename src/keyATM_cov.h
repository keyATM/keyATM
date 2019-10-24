#ifndef __keyATM_cov__INCLUDED__
#define __keyATM_cov__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMcov : virtual public keyATMbase
{
  public:  
    //
    // Parameters
    //
    MatrixXd Alpha;
    int num_cov;
    MatrixXd Lambda;
    MatrixXd C;

    double mu;
    double sigma;

    // During the sampling
      std::vector<int> topic_ids;
      std::vector<int> cov_ids;

      double Lambda_current;
      double llk_current;
      double llk_proposal;
      double diffllk;
      double r, u;

      // Slice sampling
      double start, end, previous_p, new_p, newlikelihood, slice_, current_lambda;
      double store_loglik;
      double newlambdallk;
    
    //
    // Functions
    //

    // Constructor
    keyATMcov(List model_, const int iter_, const int output_per_) :
      keyATMbase(model_, iter_, output_per_) {};

    // Read data
    void read_data_specific();

    // Initialization
    void initialize_specific();

    // Iteration
    void iteration_single(int &it);
    void sample_parameters(int &it);
    void sample_lambda();
    double alpha_loglik();
    double loglik_total();

    void sample_lambda_slice();
    double likelihood_lambda();
    void proposal_lambda(int& k);

};


#endif


