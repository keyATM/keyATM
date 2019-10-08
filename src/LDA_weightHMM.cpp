#include "LDA_weightHMM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */

LDAhmm::LDAhmm(List model_, const int iter_, const int output_per_) :
  keyATMhmm(model_, iter_, output_per_) // pass to parent!
{
  // Nothing to add
}


void LDAhmm::initialize_specific()
{
  // Initialize Psk
  Psk = MatrixXd::Zero(num_time, num_states);

  // Initialize S_est
  // Use multinomial distribution (with flat probability)
  // to decide the number of each state
  // and push it into S_est.
  VectorXi S_est_num = VectorXi::Zero(num_states);
  VectorXd S_est_temp = VectorXd::Zero(num_states);
  double cumulative = 1.0 / num_states;
  double u;
  int index;
  for(int i=0; i<num_states; i++){
    S_est_temp(i) = cumulative * (i+1);
  }

  for(int j=0; j<num_time; j++){
    u = R::runif(0, 1);
    for(int i=0; i<num_states; i++){
      if(u < S_est_temp(i)){
        index = i;
        break;
      }
    }
    S_est_num(index) += 1;
  }


  S_est = VectorXi::Zero(num_time);
  S_count = S_est_num;
  int count;
  index = 0;
  for(int i=0; i<num_states; i++){
    count = S_est_num(i);
    for(int j=0; j<count; j++){
      S_est(index) = i;
      index += 1;
    }
  }

  // Initializae P_est
  P_est = MatrixXd::Zero(num_states, num_states);
  double prob;
  for(int i=0; i<=(index_states-1); i++){
    prob = R::rbeta(1.0, 1.0);
    P_est(i, i) = prob;
    P_est(i, i+1) = 1-prob;
  }
  P_est(index_states, index_states) = 1;


  // Initialize alphas;
  alphas = MatrixXd::Constant(num_states, num_topics, 50.0/num_topics);

  // Initialize variables we use in the sampling
  logfy = VectorXd::Zero(num_states);
  st_k = VectorXd::Zero(num_states);
  logst_k = VectorXd::Zero(num_states);
  state_prob_vec = VectorXd::Zero(num_states);
  
  states_start = VectorXi::Zero(num_states);
  states_end = VectorXi::Zero(num_states);


  // Initialization for LDA weights 

  n_kv = MatrixXd::Zero(num_topics, num_vocab);
  n_dk = MatrixXd::Zero(num_doc, num_topics);
  n_k = VectorXd::Zero(num_topics);
  n_k_noWeight = VectorXd::Zero(num_topics);

  int z, w;
  int doc_len;
  IntegerVector doc_z, doc_w;


  // Construct data matrices
  for(int doc_id = 0; doc_id < num_doc; doc_id++){
    doc_z = Z[doc_id], doc_w = W[doc_id];
    doc_len = doc_each_len[doc_id];

    for(int w_position = 0; w_position < doc_len; w_position++){
      z = doc_z[w_position], w = doc_w[w_position];

      n_kv(z, w) += vocab_weights(w);
      n_k(z) += vocab_weights(w);
      n_k_noWeight(z) += 1.0;
      n_dk(doc_id, z) += 1.0;
    }
  }
  

  // Use during the iteration
  z_prob_vec = VectorXd::Zero(num_topics);


}


void LDAhmm::iteration_single(int &it)
{
  x_ = -1;  // we do not use x_ in LDA HMM
  doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle

  for (int ii = 0; ii < num_doc; ii++){
    doc_id_ = doc_indexes[ii];
    doc_z = Z[doc_id_], doc_w = W[doc_id_];
    doc_length = doc_each_len[doc_id_];

    alpha = alphas.row(S_est(doc_id_)).transpose(); // select alpha for this document
    
    token_indexes = sampler::shuffled_indexes(doc_length); //shuffle
    
    // Iterate each word in the document
    for (int jj = 0; jj < doc_length; jj++){
      w_position = token_indexes[jj];
      z_ = doc_z[w_position], w_ = doc_w[w_position];
    
      new_z = sample_z(alpha, z_, x_, w_, doc_id_);
      doc_z[w_position] = new_z;
    }
    
    Z[doc_id_] = doc_z;
  }

  sample_parameters(it);
}


int LDAhmm::sample_z(VectorXd &alpha, int &z, int &x,
                         int &w, int &doc_id)
{
  // remove data
  n_kv(z, w) -= vocab_weights(w);
  n_k(z) -= vocab_weights(w);
  n_k_noWeight(z) -= 1.0;
  n_dk(doc_id, z) -= 1;

  new_z = -1; // debug


  for (int k = 0; k < num_topics; ++k){

    numerator = (beta + n_kv(k, w)) *
      (n_dk(doc_id, k) + alpha(k));

    denominator = ((double)num_vocab * beta + n_k(k)) ;

    z_prob_vec(k) = numerator / denominator;
  }

  sum = z_prob_vec.sum(); // normalize
  new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample


  // add back data counts
  n_kv(new_z, w) += vocab_weights(w);
  n_k(new_z) += vocab_weights(w);
  n_k_noWeight(new_z) += 1.0;
  n_dk(doc_id, new_z) += 1;

  return new_z;
}

double LDAhmm::loglik_total()
{
  loglik = 0.0;
  for (int k = 0; k < num_topics; k++){
    for (int v = 0; v < num_vocab; v++){ // word
      loglik += mylgamma(beta + n_kv(k, v) / vocab_weights(v) ) - mylgamma(beta);
    }

    // word normalization
    loglik += mylgamma( beta * (double)num_vocab ) - mylgamma(beta * (double)num_vocab + n_k_noWeight(k) );
  }


  for (int d = 0; d < num_doc; d++){
    // z
    alpha = alphas.row(S_est(doc_id_)).transpose(); // Doc alpha, column vector  
    
    loglik += mylgamma( alpha.sum() ) - mylgamma( doc_each_len[d] + alpha.sum() );
    for (int k = 0; k < num_topics; k++){
      loglik += mylgamma( n_dk(d,k) + alpha(k) ) - mylgamma( alpha(k) );
    }

    // HMM part
    state_id = S_est(d);
    loglik += log( P_est(state_id, state_id) );
  }


  return loglik;
}




