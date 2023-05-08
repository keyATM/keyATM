#ifndef __LDAbase__INCLUDED__
#define __LDAbase__INCLUDED__
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

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
    LDAbase(List model_) :
      keyATMmeta(model_) {};

    // Variables
    MatrixXd n_kv;
    VectorXd n_k;
    VectorXd n_k_noWeight;

    // Functions
    // In LDA, we do not need to read and initialize X
    virtual void read_data_common() override final;
    virtual void initialize_common() override final;
    virtual void parameters_store(int r_index) override final;
    virtual int sample_z(VectorXd &alpha, int z, int s,
                         int w, int doc_id) override final;
};

#endif

