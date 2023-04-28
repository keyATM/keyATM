#ifndef __keyATM_multi_base__INCLUDED__
#define __keyATM_multi_base__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM_meta.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMmultibase : virtual public keyATMmeta
{
public:  
  
  // Data
  int num_corpus;
  IntegerVector corpus_id;
  std::vector< int > total_words_corpus;
  ListOf<IntegerVector> global_id;
  IntegerVector num_doc_all;
  
  //Corpus-specific weights
  MatrixXd vocab_weights_corpus;
  std::vector< int > num_vocab_corpus;
  
  //
  // Sufficient stats
  //
  std::vector<MatrixXd> n_s0_kv_all;
  MatrixXd n_s0_k_all;
  
  
  //
  // Parameters
  //
  int estimate_alpha;
  int store_alpha;
  //std::vector<MatrixXd> Alpha_all;
  MatrixXd Alpha;
  
  std::vector<int> topic_ids;
  VectorXd keep_current_param;
  int g_doc_id;
  double store_loglik;
  double newalphallk;
  
  // in alpha_loglik
  MatrixXd ndk_a;
  
  //
  // Functions
  //
  
  // Constructor
  keyATMmultibase(List model_, const int iter_) :
    keyATMmeta(model_, iter_) {};
  
  // Read data
  virtual void read_data_specific();
  
  // Initialization
  virtual void initialize_specific();
  
  // Iteration
  virtual void iteration_single(int it);
  virtual void sample_parameters(int it);
  int sample_z(Ref<VectorXd> alpha, int z, int s,
                                int w, int doc_id);
  int sample_s(int z, int s, int w, int doc_id);
  void sample_alpha(int corpus_);
  virtual double loglik_total();
  double alpha_loglik(int k, int corpus_);
  double loglik_total_label();
  void store_pi_iter(int r_index);
  List return_model();
};


#endif

