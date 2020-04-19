#ifndef __keyATM_label__INCLUDED__
#define __keyATM_label__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM_meta.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMlabel : virtual public keyATMmeta
{
  public:  
    //
    // Parameters
    //
    int estimate_alpha;
    int store_alpha;

    double start, end, previous_p, new_p, newlikelihood, slice_;
    std::vector<int> topic_ids;
    VectorXd keep_current_param;
    double store_loglik;
    double newalphallk;

    double loglik;
    double fixed_part;
          // 
    MatrixXd Alpha;
    MatrixXd label_dk;
    VectorXd Alpha_sum_vec;


    // use during initialization
    int doc_label;
    // std::vector<int> label_vec;
    IntegerVector label_vec;

    // in alpha_loglik
    MatrixXd ndk_a;
    //
    // Functions
    //

    // Constructor
    keyATMlabel(List model_, const int iter_) :
      keyATMmeta(model_, iter_) {};

    // Read data
    void read_data_specific();

    // Initialization
    void initialize_specific();

    // Iteration
    virtual void iteration_single(int it);
    void sample_parameters(int it);
    void sample_alpha();
    double alpha_loglik_label(int k);
    virtual double loglik_total();
    double loglik_total_label();
};


#endif

