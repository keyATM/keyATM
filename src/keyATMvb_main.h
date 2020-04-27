#ifndef __keyATMvb_main__INCLUDED__
#define __keyATMvb_main__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include <algorithm>
#include "sampler.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMvb
{
  // Variational Bayes for keyATM models
  //  VB for z and s to avoid label switching issues

  public:
    //
    // Data
    //
    List model;
    List W, Z, S;
    std::string model_name;
    StringVector vocab;
    List keywords_list;
    List priors_list;
    List options_list;
    List vb_options;
    List Perplexity;
    NumericVector Perplexity_value;
    NumericVector Perplexity_iter;
    int use_weight;
    List stored_values;

    //
    // Parameters
    //
    int num_doc;
    int num_topics;
    int num_vocab;
    int regular_k, keyword_k;
    int total_words;
    double total_words_weighted;
    double beta, Vbeta, beta_s;
    MatrixXd prior_gamma;
    std::vector<int> doc_each_len;
    std::vector<double> doc_each_len_weighted;

    MatrixXd alphas;
    VectorXd alpha_d;  // document level sum of alpha

    std::vector< std::unordered_set<int> > keywords;
    std::unordered_set<int> keywords_all;
    std::vector<int> keywords_num;

    VectorXd vocab_weights;

    std::vector<std::vector<vector<double>>> qz;
    std::vector<std::vector<vector<double>>> qs;

    MatrixXd n_s0_kv;
    MatrixXd n_s1_kv;
    MatrixXd n_dk;
    VectorXd n_s0_k;
    VectorXd n_s1_k;


    // During the iteration
    std::vector<int> doc_indexes;
    VectorXd z_prob_vec;
    VectorXd s_prob_vec;
    VectorXd s0_temp;
    VectorXd s1_temp;



    List doc_w;
    List doc_z;
    List doc_s;

    // For calculating perplexity
    std::vector<int> ppl_doc_indexes;
    int num_doc_perp;
    double ppl_words;

    //
    // Functions
    //
    keyATMvb(List model_);
    ~keyATMvb();
    void fit();
    List return_model();
      virtual void get_QOI();

    // Read
    void read_data();
    void read_data_common();
      void read_data_common_alpha();
        void read_data_common_alpha_base();
        void read_data_common_alpha_cov();
        void read_data_common_alpha_hmm();
      virtual void read_data_words();
        void read_data_documents();
        void read_data_keywords();
    virtual void read_data_specific();
    
    // Initialize
    void initialize();
    void initialize_common();
      virtual void initialize_common_MCMCcount();
      virtual void initialize_common_q();
      virtual void initialize_common_qz(int doc_id, int w,
                                        int z, int s, vector<double> &qzdk);
      virtual void initialize_common_qs(int doc_id, int w, 
                                        int z, int s, vector<double> &qsds);
      void initialize_common_expectation();
      void initialize_weightedlen();
    void initialize_specific();
    
    // Iteration
    void iteration();
      void iteration_single();
      virtual void update_q();
      virtual void update_decrese_count(int doc_id, int w_position, int v);
      virtual void update_increase_count(int doc_id, int w_position, int v);

    virtual double calc_perplexity(int iter);
      void store_perplexity(int iter, double ppl);
};

#endif
