#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

/* LOGGING NAMESPACE */
namespace iter_log {
  NumericMatrix iter_alpha;
  NumericMatrix iter_n_x0_k;
  NumericMatrix iter_n_x1_k;
  NumericMatrix iter_p_k;

  template <class EigenVec>
  void store_values(NumericMatrix &rmat, int &iter_num, EigenVec &vec);

  void initialize(int &iter_num, int &num_topics){
  static int set = 1;

  if(set){
    NumericMatrix temp(iter_num, num_topics);
    iter_alpha = clone(temp);
    iter_n_x0_k = clone(temp);
    iter_n_x1_k = clone(temp);
    iter_p_k = clone(temp);

    set = 0; // initialize only once
  }
}

void store_alpha(int &iter_num, VectorXd &vec){
  store_values(iter_alpha, iter_num, vec);
}

void store_n_x0_k(int &iter_num, VectorXi &vec){
  store_values(iter_n_x0_k, iter_num, vec);
}

void store_n_x1_k(int &iter_num, VectorXi &vec){
  store_values(iter_n_x1_k, iter_num, vec);
}

void store_p_k(int &iter_num, VectorXd &vec){
  store_values(iter_p_k, iter_num, vec);
}

template <class EigenVec>
  void store_values(NumericMatrix &rmat, int &iter_num, EigenVec &vec){
    int num = vec.size();

    for(int i=0; i<num; i++){
      rmat(iter_num, i) = vec(i);
    }
  }
}
