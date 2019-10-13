#ifndef __keyATM_weightHMM__INCLUDED__
#define __keyATM_weightHMM__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM_HMM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class LDAhmm : public keyATMhmm
{
  public:
    // Constructor
    LDAhmm(List model_, const int iter_, const int output_per_);

    // Parameters in LDA HMM
    MatrixXd n_kv;
    VectorXd n_k;
    VectorXd n_k_noWeight;

    // Functions
		void read_data_specific();
    void initialize_specific();
    int sample_z(VectorXd &alpha, int &z, int &x,
                 int &w, int &doc_id);
    void iteration_single(int &it);
    double loglik_total();

};

#endif


