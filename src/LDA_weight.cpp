#include "LDA_weight.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;


void LDAweight::iteration_single(int it)
{ // Single iteration
  int doc_id_;
  int doc_length;
  int w_, z_, s_;
  int new_z;
  int w_position;

  s_ = -1;  // we do not use x_ in LDA weight
  doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle

  for (int ii = 0; ii < num_doc; ii++) {
    doc_id_ = doc_indexes[ii];
    doc_z = Z[doc_id_], doc_w = W[doc_id_];
    doc_length = doc_each_len[doc_id_];
    
    token_indexes = sampler::shuffled_indexes(doc_length); //shuffle
    
    // Iterate each word in the document
    for (int jj = 0; jj < doc_length; jj++) {
      w_position = token_indexes[jj];
      z_ = doc_z[w_position], w_ = doc_w[w_position];
    
      new_z = sample_z(alpha, z_, s_, w_, doc_id_);
      doc_z[w_position] = new_z;
    }
    
    Z[doc_id_] = doc_z;
  }
  sample_parameters(it);

}


double LDAweight::loglik_total()
{
  double loglik = 0.0;
  for (int k = 0; k < num_topics; k++) {
    for (int v = 0; v < num_vocab; v++) { // word
      loglik += mylgamma(beta + n_kv(k, v)) - mylgamma(beta);
    }
    // word normalization
    loglik += mylgamma( beta * (double)num_vocab ) - mylgamma(beta * (double)num_vocab + n_k(k) );

  }

  // z
  for (int d = 0; d < num_doc; d++) {
    loglik += mylgamma( alpha.sum() ) - mylgamma( doc_each_len_weighted[d] + alpha.sum() );
    for (int k = 0; k < num_topics; k++) {
      loglik += mylgamma( n_dk(d,k) + alpha(k) ) - mylgamma( alpha(k) );
    }
  }
  return loglik;
}

