#ifndef __keyATM_basic__INCLUDED__
#define __keyATM_basic__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMbasic : virtual public keyATMbase
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

      // in alpha_loglik
      MatrixXd ndk_a;
    
    //
    // Functions
    //

    // Constructor
    keyATMbasic(List model_, const int iter_, const int output_per_) :
      keyATMbase(model_, iter_, output_per_) {};

    // Read data
    void read_data_specific();

    // Initialization
    void initialize_specific();

    // Iteration
    void iteration_single(int &it);
    void sample_parameters(int &it);
    void sample_alpha();
    double alpha_loglik();
    double loglik_total();
};


#endif

