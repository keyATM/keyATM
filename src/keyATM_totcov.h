#ifndef __keyATM_totcov__INCLUDED__
#define __keyATM_totcov__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMtotcov : public keyATMbase
{

  public:  
    // Settings
      // TOT
      int logsumexp_approx;
      int use_mom;

    //
    // Parameters
    //
      // COV
      VectorXd timestamps;  // time stamps (document share the same time)
      MatrixXd beta_params;  // parameter for time Beta, K \times 2

      VectorXd beta_tg;  // apply tgamma 
      VectorXd beta_lg;  // apply lgamma 
      VectorXd beta_tg_base;  // apply tgamma 
      VectorXd beta_lg_base;  // apply lgamma 

      // TOT
      MatrixXd Alpha;
      int num_cov;
      MatrixXd Lambda;
      MatrixXd C;

    // MH sampling
    double mu;
    double sigma;
    double mh_sigma;

    // Sampling info
    std::vector<int> mh_info {0,0};

    // During the sampling
      // Sample Beta parameters
      vector<vector<double>> store_t;  
        // store time stamp to estiamte Beta parameters
      VectorXd timestamps_k;  // time stamps for topic k

      // mh_single
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

        double loglik;
        double fixed_part;

        // in alpha_loglik
        MatrixXd ndk_a;

      // In sampling z
      double beta_a;
      double beta_b;
      double check_frac;
      double timestamp_d;
      int use_log;

      // In sampling betaparam
      double beta_mean;
      double beta_var;
      double current_param;
      double temp_beta_loglik;
      double ts_g1;  // parameters for gamma
      double ts_g2;  // parameters for gamma
    
    //
    // Functions
    //

    // Constructor
    keyATMtotcov(List model_, const int iter_, const int output_per_);

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

    void sample_lambda_mh_single();
    void sample_lambda_mh();
    void sample_lambda_slice();
    double likelihood_lambda();
    void proposal_lambda(int& k);

    int sample_z(VectorXd &alpha, int &z, int &x, int &w, int &doc_id);
    int sample_z_log(VectorXd &alpha, int &z, int &x, int &w, int &doc_id);
    void sample_betaparam();
    double beta_loglik(const int &k, const int &i);
    void verbose_special(int &r_index);  // store sampled beta param
};

# endif
