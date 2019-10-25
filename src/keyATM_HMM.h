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

class keyATMhmm : virtual public keyATMbase
{
  public:
    // Data
    VectorXi time_index;
    int num_time;  // number of time segment
    VectorXi time_doc_start;
    VectorXi time_doc_end;

    // Parameters
    // In this implementation, we should consider
    // Y_n in Chib (1998) as a 'block' of documents that
    // share the same time index.
    // 'Time-Blocked' Change Point
    int num_states;
    int index_states;  // num_states - 1
    int store_transition_matrix;
    MatrixXd Psk;    // (num_time, num_states)
    VectorXi S_est;  // stores state index, (num_time)
    VectorXi S_count;  // stores the count of each state
      // Sec 2.3 in Chib (1998)
      // "the elements of p_ij of P may be simulated from P|S_n
      // without regard to the sampling model for the data"
      // so we just store the count of each state
    MatrixXd P_est;  // (num_states, num_states)
    MatrixXd alphas;  // (num_states, num_topics)
    double loglik;
    double fixed_part;

    // Constructor
    keyATMhmm(List model_, const int iter_, const int output_per_) :
      keyATMbase(model_, iter_, output_per_) {};

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
    // Utilities
    int get_state_index(const int &doc_id);
  
    // Read data and Initialize
    void read_data_specific() final;
    void initialize_specific() final;
  
    // Iteration
    virtual void iteration_single(int &it);
    void sample_parameters(int &it);

    void sample_alpha();
    void sample_alpha_state(int &state, int &state_start, int &state_end);
    double alpha_loglik(int &state_start, int &state_end);

    void sample_forward();  // calculate Psk
    void sample_backward();  // sample S_est
    void sample_P();  // sample P_est
    void store_S_est();
    void store_P_est();

    double polyapdfln(int &t, VectorXd &alpha);
    virtual double loglik_total();
    void verbose_special(int &r_index);
};

#endif

