#ifndef __keyATM_multi_base__INCLUDED__
#define __keyATM_multi_base__INCLUDED__
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

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
    //
    // Data
    //
    List C;
    int num_corpus;
    IntegerVector corpus_id;
    std::vector< int > total_words_corpus;

    // Corpus-specific weights
    MatrixXd vocab_weights_corpus;
    std::vector< int > num_vocab_corpus;

    // Sufficient Stats
    std::vector<MatrixXd> n_s0_c0_kv_all;

    MatrixXd n_s0_c0_k_all;
    MatrixXd n_c0_k_all;
    MatrixXd n_c1_k_all;

    MatrixXd n_s0_c1_kv;
    VectorXd n_s0_c1_k;

    VectorXd n_c0_k;
    VectorXd n_c1_k;

    SparseMatrix<double,RowMajor> n_s1_kv_multi;
    
    std::vector<double> doc_each_len_weighted_multi;

    IntegerVector doc_c;

    //
    // Parameters
    //
    int estimate_alpha;
    int store_alpha;
    VectorXd alpha;
    double beta_c;
    double Vbeta_c;
    double Vbeta_s;
    MatrixXd prior_omega;

    std::vector<int> topic_ids;
    VectorXd keep_current_param;
    double store_loglik;
    double newalphallk;
    

    // in alpha_loglik
    MatrixXd ndk_a;

    //
    // Functions
    //

    // Constructor
    keyATMmultibase(List model_) :
      keyATMmeta(model_) {};

    // Initialize
    virtual void read_data_specific() override final;
    virtual void initialize_specific() override final;

    // Resume
    virtual void resume_initialize_specific() override final;

    void create_sufficient_stats();

    // Iteration
    virtual void iteration_single(int it) override;
    virtual void sample_parameters(int it) override final;
    virtual int sample_z(VectorXd &alpha, int z, int s, int w, int c, int doc_id);
    virtual int sample_s(int z, int s, int w, int c, int doc_id);
    virtual int sample_c(int z, int s, int w, int c, int doc_id);
    void sample_alpha();
    double alpha_loglik(int k);
    virtual double loglik_total() override;
};


#endif
