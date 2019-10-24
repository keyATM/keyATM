#include "LDA_weightHMM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */


void LDAhmm::read_data_specific()
{
  model_settings = model["model_settings"];
  num_states = model_settings["num_states"];
  index_states = num_states - 1;

  IntegerVector time_index_R = model_settings["time_index"];
  time_index = Rcpp::as<Eigen::VectorXi>(time_index_R);
  num_time = time_index.maxCoeff();
  time_index = time_index.array() - 1;  // adjust index
  time_doc_start = VectorXi::Zero(num_time);
  time_doc_end = VectorXi::Zero(num_time);

  // Prepare start and end of time
  int index_prev = -1;
  int index;
  int store_index = 0;
  for(int d = 0; d < num_doc; d++){
     index = time_index[d]; 
    if(index != index_prev){
      time_doc_start[store_index] = d;
      index_prev = index;
      store_index += 1;
    }
  }

  for(int s = 0; s < num_time - 1; s++){
    time_doc_end(s) = time_doc_start(s+1) - 1;  
  }
  time_doc_end(num_time - 1) = num_doc - 1;

  store_transition_matrix = options_list["store_transition_matrix"];
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
  for(int i = 0; i < num_states; i++){
    S_est_temp(i) = cumulative * (i+1);
  }

  for(int j = 0; j < num_time; j++){
    u = R::runif(0, 1);
    for(int i = 0; i < num_states; i++){
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
  for(int i = 0; i < num_states; i++){
    count = S_est_num(i);
    for(int j = 0; j < count; j++){
      S_est(index) = i;
      index += 1;
    }
  }

  // Initializae P_est
  P_est = MatrixXd::Zero(num_states, num_states);
  double prob;
  for(int i = 0; i <= (index_states - 1); i++){
    prob = R::rbeta(1.0, 1.0);
    P_est(i, i) = prob;
    P_est(i, i + 1) = 1 - prob;
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

}


void LDAhmm::iteration_single(int &it)
{
  x_ = -1;  // we do not use x_ in LDA HMM
  doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle

  for (int ii = 0; ii < num_doc; ii++){
    doc_id_ = doc_indexes[ii];
    doc_z = Z[doc_id_], doc_w = W[doc_id_];
    doc_length = doc_each_len[doc_id_];

    alpha = alphas.row(get_state_index(doc_id_)).transpose(); // select alpha for this document
    
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
    alpha = alphas.row(get_state_index(doc_id_)).transpose(); // Doc alpha, column vector  
    
    loglik += mylgamma( alpha.sum() ) - mylgamma( doc_each_len[d] + alpha.sum() );
    for (int k = 0; k < num_topics; k++){
      loglik += mylgamma( n_dk(d,k) + alpha(k) ) - mylgamma( alpha(k) );
    }

  }

  // HMM part
  for(int t = 0; t < num_time; t++){
    state_id = S_est(t);
    loglik += log( P_est(state_id, state_id) );
  }

  return loglik;
}




