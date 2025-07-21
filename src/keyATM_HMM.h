#ifndef __keyATM_HMM__INCLUDED__
#define __keyATM_HMM__INCLUDED__
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#include "keyATM_meta.h"
#include "sampler.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMhmm : virtual public keyATMmeta {
public:
  // Data
  VectorXi time_index;
  int num_time; // number of time segment
  VectorXi time_doc_start;
  VectorXi time_doc_end;

  // Parameters
  // In this implementation, we should consider
  // Y_n in Chib (1998) as a 'block' of documents that
  // share the same time index.
  // 'Time-Blocked' Change Point
  int num_states;
  int index_states; // num_states - 1
  int store_transition_matrix;
  MatrixXd Prk;     // (num_time, num_states)
  VectorXi R_est;   // stores state index, (num_time)
  VectorXi R_count; // stores the count of each state
                    // Sec 2.3 in Chib (1998)
                    // "the elements of p_ij of P may be simulated from P|R_n
                    // without regard to the sampling model for the data"
                    // so we just store the count of each state
  MatrixXd P_est;  // (num_states, num_states)
  MatrixXd alphas; // (num_states, num_topics)

  // Constructor
  keyATMhmm(List model_) : keyATMmeta(model_) {};

  // During sampling
  // sample_forward()
  VectorXd logfy; // (num_states)
  VectorXd rt_1l;
  VectorXd rt_k;
  VectorXd logrt_k;

  VectorXd state_prob_vec;

  // Sample alpha
  VectorXi states_start;
  VectorXi states_end;

  // Slice Sampling
  std::vector<int> topic_ids;
  VectorXd keep_current_param;
  MatrixXd ndk_a;

  //
  // Functions
  //
  // Utilities
  int get_state_index(const int doc_id);

  // Read data and Initialize
  virtual void read_data_specific() override final;
  virtual void initialize_specific() override final;

  // Resume
  virtual void resume_initialize_specific() override final;

  // Iteration
  virtual void iteration_single(int it) override;
  virtual void sample_parameters(int it) override final;

  void sample_alpha();
  void sample_alpha_state(int state, int state_start, int state_end);
  double alpha_loglik(int k, int state_start, int state_end);

  void sample_forward();  // calculate Prk
  void sample_backward(); // sample R_est
  void sample_P();        // sample P_est
  void store_R_est();     // store state
  void store_P_est();     // store the transition matrix
  void keep_P_est();      // keep the latest transition matrix

  double polyapdfln(int t, VectorXd &alpha);
  virtual double loglik_total() override;
  virtual void verbose_special(int r_index) override;
};

#endif
