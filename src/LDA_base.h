#ifndef __LDAbase__INCLUDED__
#define __LDAbase__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM_meta.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class LDAbase : virtual public keyATMmeta
{
  // Base function for the Weighted LDA models
  // This inherits keyATMmeta class.
  
  public:
    // Constructor
    LDAbase(List model_, const int iter_) :
      keyATMmeta(model_, iter_) {};
  
    // Variables
    MatrixXd n_kv;
    VectorXd n_k;
    VectorXd n_k_noWeight;

    // Functions
    // In LDA, we do not need to read and initialize X
    virtual void read_data_common();
    virtual void initialize_common();
    virtual void iteration_single(int it) = 0;
    void parameters_store(int r_index); 
    virtual int sample_z(VectorXd &alpha, int z, int s,
                         int w, int doc_id) final;
    virtual double loglik_total() = 0;
};

#endif

