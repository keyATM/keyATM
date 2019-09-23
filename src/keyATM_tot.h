#ifndef __keyATM_tot__INCLUDED__
#define __keyATM_tot__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMtot : public keyATMbase
{
  public:
    // Settings  
    int logsumexp_approx;

    // Parameters
    VectorXd timestamps;  // time stamps (document share the same time)
    VectorXi time_topics;  // consider time trends only for topics in this vector
    int num_topics_time;  // number of topics to sample beta
    MatrixXd beta_params;  // parameter for time Beta, K \times 2

    VectorXd beta_tg;  // apply tgamma 
    VectorXd beta_lg;  // apply lgamma 
    VectorXd beta_tg_base;  // apply tgamma 
    VectorXd beta_lg_base;  // apply lgamma 
    
    // Constructor
    keyATMtot(List model_, const int iter_, const int output_per_);

    // During sampling
      // Sample Beta parameters
      vector<vector<double>> store_t;  
        // store time stamp to estiamte Beta parameters
      VectorXd timestamps_k;  // time stamps for topic k


      // Sample alpha

        double start, end, previous_p, new_p, newlikelihood, slice_;
        std::vector<int> topic_ids;
        VectorXd keep_current_param;
        double store_loglik;
        double newalphallk;

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
  
    // Read data and Initialize
    // Read data
    void read_data_specific();

    // Initialization
    void initialize_specific();

    // Iteration
    void iteration_single(int &it);
    int sample_z(VectorXd &alpha, int &z, int &x, int &w, int &doc_id);
    int sample_z_log(VectorXd &alpha, int &z, int &x, int &w, int &doc_id);
    void sample_parameters(int &it);
    void sample_betaparam();
    void sample_alpha();
    double alpha_loglik();
    double beta_loglik(const int &k, const int &i);
    double loglik_total();
    void verbose_special(int &r_index);  // store sampled beta param
};

#endif


