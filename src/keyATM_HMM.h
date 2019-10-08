#ifndef __keyATM_HMM__INCLUDED__
#define __keyATM_HMM__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMhmm : public keyATMbase
{
  public:
    // Data
    VectorXi time_index;
		int time_num;
    VectorXi time_doc_start;
    VectorXi time_doc_end;

    // Parameters
    int num_states;
    int index_states;  // num_states - 1
    MatrixXd Psk;    // (num_doc, num_states)
    VectorXi S_est;  // stores state index, (num_doc)
    VectorXi S_count;  // stores the count of each state
    MatrixXd P_est;  // (num_states, num_states)
    MatrixXd alphas;  // (num_states, num_topics)
    double loglik;
    double fixed_part;

    // Constructor
    keyATMhmm(List model_, const int iter_, const int output_per_);

    // During sampling
      // sample_forward()
      VectorXd logfy;  // (num_states)
      VectorXd st_1l;
      VectorXd st_k;
      VectorXd logst_k;
      double logsum;
      int added;

      int state_id;
      VectorXd state_prob_vec;
      double pii;

      // Sample alpha
      VectorXi states_start;
      VectorXi states_end;

        // Slice Sampling
        double start, end, previous_p, new_p, newlikelihood, slice_;
        std::vector<int> topic_ids;
        VectorXd keep_current_param;
        double store_loglik;
        double newalphallk;
        MatrixXd ndk_a;
  
    // 
    // Functions
    //
  
    // Read data and Initialize
    void read_data_specific();
    void initialize_specific();
  
    // Iteration
    void iteration_single(int &it);
    void sample_parameters(int &it);

    void sample_alpha();
    void sample_alpha_state(int &state, int &state_start, int &state_end);
    double alpha_loglik(int &state_start, int &state_end);

    void sample_forward();  // calculate Psk
    void sample_backward();  // sample S_est
    void sample_P();  // sample P_est
    void store_S_est();

    double polyapdfln(int &doc_id, VectorXd &alpha);
    double loglik_total();
    void verbose_special(int &r_index);
};

#endif

