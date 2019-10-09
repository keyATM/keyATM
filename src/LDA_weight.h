#ifndef __LDA_weight__INCLUDED__
#define __LDA_weight__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class LDAweight : public keyATMbase
{
  public:  
    //
    // Parameters
    //

    // Slice Sampling
    double start, end, previous_p, new_p, newlikelihood, slice_;
    std::vector<int> topic_ids;
    VectorXd keep_current_param;
    double store_loglik;
    double newalphallk;

    double loglik;
    double fixed_part;

      // in alpha_loglik
      MatrixXd ndk_a;

    // Latent variables for LDA weight
    MatrixXd n_kv;
      // n_dk is already defined in parent class
    VectorXd n_k;
    VectorXd n_k_noWeight;
    
    //
    // Functions
    //

    // Constructor
    LDAweight(List model_, const int iter_, const int output_per_);

    // Read data
    void read_data_specific();

    // Initialization
    void initialize_specific();

    // Iteration
    void iteration_single(int &it);
    void sample_parameters(int &it);
    void sample_alpha();
    int sample_z(VectorXd &alpha, int &z, int &x,
                 int &w, int &doc_id);
    double alpha_loglik();
    double loglik_total();
};

#endif
