#ifndef __keyATM__INCLUDED__
#define __keyATM__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMbase
{
  public:
    //
    // Parameters
    //
    int iter;
    int output_per;
    double eta_1;
    double eta_2;
    double eta_1_regular;
    double eta_2_regular;
    int use_weight;
    int store_theta;
    int thinning;

    double slice_A; // parameter for slice sampling 

    // Data
    List model;
    List W, Z, X;
    StringVector vocab;
    NumericVector nv_alpha;
    MatrixXd prior_gamma;
    double beta, beta_s;
    int regular_k, keyword_k;
    List keywords_list;
    List model_fit;
    std::vector<int> doc_each_len;
    
    int num_vocab, num_doc, total_words;

    List options_list;
    List Z_tables;
    List priors_list;
    List model_settings;
    List stored_values;
    NumericMatrix Z_table;

    // alpha
    int num_topics;
    VectorXd alpha;

    // Keywords
    std::vector< std::unordered_set<int> > keywords;
    std::vector<int> keywords_num;

    // Latent Variables
    MatrixXd n_x0_kv;
    SparseMatrix<double,RowMajor> n_x1_kv;
    typedef Eigen::Triplet<double> Triplet;
    MatrixXd n_dk;
    VectorXd n_x0_k;
    VectorXd n_x0_k_noWeight;
    VectorXd n_x1_k;
    VectorXd n_x1_k_noWeight;
    VectorXd vocab_weights;

    // Use during the iteration
      // Declaration
      std::vector<int> doc_indexes;
      int doc_id_;
      std::vector<int> token_indexes;
      IntegerVector doc_x, doc_z, doc_w;
      int w_position;
      int x_, z_, w_;
      int doc_length;
      // int new_z, new_x;  // defined in sample_z and sample_x
  
      // sample_z
      VectorXd z_prob_vec;
      int new_z;
      double numerator, denominator;
      double sum;

      // sample_x
      int new_x;
      double x0_prob;
      double x1_prob;
      int k;

      // sample alpha
      double alpha_sum_val;
      double min_v;
      double max_v;
      int max_shrink_time;

      // gammaln_sum
      double gammaln_val;


    //
    // Functions
    //
    keyATMbase(List model_, const int iter_, const int output_per_);
    ~keyATMbase();
    void fit();

    // Reading and Initialization
    void read_data();
    virtual void read_data_common();
    virtual void read_data_specific() = 0;

    void initialize();
    virtual void initialize_common();
    virtual void initialize_specific() = 0;

    // Sampling
    void iteration();
    virtual void iteration_single(int &it) = 0;
    virtual void sample_parameters(int &it) = 0;

    virtual int sample_z(VectorXd &alpha, int &z, int &x,
                         int &w, int &doc_id);
    int sample_x(VectorXd &alpha, int &z, int &x,
                   int &w, int &doc_id);

    void sampling_store(int &r_index);
    void store_theta_iter(int &r_index);

    virtual void verbose_special(int &r_index);

    virtual double loglik_total() = 0;

    // Utilities
    double vmax, vmin;

    double gammapdfln(const double &x, const double &a, const double &b);
    double betapdf(const double &x, const double &a, const double &b);
    double betapdfln(const double &x, const double &a, const double &b);
    NumericVector alpha_reformat(VectorXd& alpha, int& num_topics);

    double gammaln_frac(const double &value, const int &count);

    // Inline functions

    double expand(double &p, const double &A)
    {
      return (-(1.0/A) * log((1.0/p) - 1.0));
    };
    double shrink(double &x, const double &A)
    {
      return (1.0 / (1.0 + exp(-A*x)));
    };
    

    double logsumexp(double x, double y, bool flg)
    {
      if (flg) return y; // init mode
      if (x == y) return x + 0.69314718055; // log(2)
      double vmin = std::min (x, y);
      double vmax = std::max (x, y);
      if (vmax > vmin + 50) {
        return vmax;
      } else {
        return vmax + std::log (std::exp (vmin - vmax) + 1.0);
      }
    };

    double logsumexp_Eigen(VectorXd &vec, const int size){
      vmax = vec.maxCoeff();
      sum = 0.0;

      for(int i = 0; i < size; i++){
        sum += exp(vec(i) - vmax);
      }
      
      return vmax + log(sum);
    }

    double mylgamma(const double &x){
      // gammaln_val = 0.0;
      // gammaln_val = lgamma(x);
      
      // Good approximation when x > 1
      //    x > 1: max abs err: 2.272e-03
      //    x > 0.5: 0.012
      //    x > 0.6: 0.008
      // Abramowitz and Stegun p.257
      
      if(x < 0.6)
        return (lgamma(x));
      else
        return ((x-0.5)*log(x) - x + 0.91893853320467 + 1/(12*x));
    };

    // These approximations are not used
    double mypow(const double &a, const double &b){
      // Reference: https://github.com/ekmett/approximate/blob/master/cbits/fast.c
      // Probably not good to use if b>1.0
      
      if(b > 1.0)
        return(pow(a,b));

      union { double d; long long x; } u = { a };
      u.x = (long long)(b * (u.x - 4606921278410026770LL) + 4606921278410026770LL);
      return u.d;
    };

    double myexp(const double &a){
      // Seems to be not very good
      union { double d; long long x; } u, v;
      u.x = (long long)(3248660424278399LL * a + 0x3fdf127e83d16f12LL);
      v.x = (long long)(0x3fdf127e83d16f12LL - 3248660424278399LL * a);
      return u.d / v.d;
    };

    double mylog(const double &a){
      // Looks fine even with large a
      union { double d; long long x; } u = { a };
      return (u.x - 4606921278410026770) * 1.539095918623324e-16;  
    };

    List return_model();
  
};

#endif
