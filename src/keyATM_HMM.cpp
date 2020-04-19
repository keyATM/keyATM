#include "keyATM_HMM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

void keyATMhmm::read_data_specific()
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
  for (int d = 0; d < num_doc; ++d) {
     index = time_index[d]; 
    if (index != index_prev) {
      time_doc_start[store_index] = d;
      index_prev = index;
      store_index += 1;
    }
  }

  for (int s = 0; s < num_time-1; ++s) {
    time_doc_end(s) = time_doc_start(s + 1) - 1;  
  }
  time_doc_end(num_time-1) = num_doc-1;

  store_transition_matrix = options_list["store_transition_matrix"];
}


void keyATMhmm::initialize_specific()
{
  // Initialize Prk
  Prk = MatrixXd::Zero(num_time, num_states);

  // Initialize R_est
  // Use multinomial distribution (with flat probability)
  // to decide the number of each state
  // and push it into R_est.
  VectorXi R_est_num = VectorXi::Constant(num_states, 1);
  VectorXd R_est_temp = VectorXd::Zero(num_states);
  double cumulative = 1.0 / num_states;
  double u;
  int index;
  for (int i = 0; i < num_states; i++) {
    R_est_temp(i) = cumulative * (i + 1);
  }

  for (int j = 0; j < num_time-num_states; j++) {
    u = R::runif(0, 1);
    for (int i = 0; i < num_states; i++) {
      if (u < R_est_temp(i)) {
        index = i;
        break;
      }
    }
    R_est_num(index) += 1;
  }


  R_est = VectorXi::Zero(num_time);
  R_count = R_est_num;
  int count;
  index = 0;
  for (int i = 0; i < num_states; i++) {
    count = R_est_num(i);
    for (int j = 0; j < count; j++) {
      R_est(index) = i;
      index += 1;
    }
  }

  // Initializae P_est
  P_est = MatrixXd::Zero(num_states, num_states);
  double prob;
  for (int i = 0; i <= (index_states-1); i++) {
    prob = R::rbeta(1.0, 1.0);
    P_est(i, i) = prob;
    P_est(i, i + 1) = 1-prob;
  }
  P_est(index_states, index_states) = 1;

  // cout << R_est_num.transpose() << endl;  //debug
  // cout << R_est.transpose() << endl;  //debug
  // cout << P_est << endl;  //debug

  // Initialize alphas;
  alphas = MatrixXd::Constant(num_states, num_topics, 50.0/num_topics);

  // Initialize variables we use in the sampling
  logfy = VectorXd::Zero(num_states);
  rt_k = VectorXd::Zero(num_states);
  logrt_k = VectorXd::Zero(num_states);
  state_prob_vec = VectorXd::Zero(num_states);
  
  states_start = VectorXi::Zero(num_states);
  states_end = VectorXi::Zero(num_states);
}


int keyATMhmm::get_state_index(const int doc_id)
{
  // Which time segment the document belongs to
  int t;
  for (t = 0; t < num_time; ++t) {
    if (time_doc_start(t) <= doc_id && doc_id <= time_doc_end(t)) {
      break;  
    }
  }
  return(R_est(t));
}


void keyATMhmm::iteration_single(int it)
{ // Single iteration
  int doc_id_;
  int doc_length;
  int w_, z_, s_;
  int new_z, new_s;
  int w_position;

  doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle

  for (int ii = 0; ii < num_doc; ++ii){
    doc_id_ = doc_indexes[ii];
    doc_s = S[doc_id_], doc_z = Z[doc_id_], doc_w = W[doc_id_];
    doc_length = doc_each_len[doc_id_];

    alpha = alphas.row(get_state_index(doc_id_)).transpose(); // select alpha for this document
    
    token_indexes = sampler::shuffled_indexes(doc_length); //shuffle
    
    // Iterate each word in the document
    for (int jj = 0; jj < doc_length; ++jj){
      w_position = token_indexes[jj];
      s_ = doc_s[w_position], z_ = doc_z[w_position], w_ = doc_w[w_position];
    
      new_z = sample_z(alpha, z_, s_, w_, doc_id_);
      doc_z[w_position] = new_z;
    
      if (keywords[new_z].find(w_) == keywords[new_z].end())	
        continue;
  
      z_ = doc_z[w_position]; // use updated z
      new_s = sample_s(alpha, z_, s_, w_, doc_id_);
      doc_s[w_position] = new_s;
    }
    
    Z[doc_id_] = doc_z;
    S[doc_id_] = doc_s;
  }

  sample_parameters(it);
}


void keyATMhmm::verbose_special(int r_index)
{
  // If there is anything special to show, write here.
}


void keyATMhmm::sample_parameters(int it)
{
  // alpha
  sample_alpha();
  
  // HMM
  sample_forward();  // calculate Prk
  sample_backward();  // sample R_est
  sample_P();  // sample P_est  

  // Store alpha and S
  int r_index = it + 1;
  if (r_index % thinning == 0 || r_index == 1 || r_index == iter) {
    Rcpp::NumericMatrix alphas_R = Rcpp::wrap(alphas);
    List alpha_iter = stored_values["alpha_iter"];
    alpha_iter.push_back(alphas_R);
    stored_values["alpha_iter"] = alpha_iter;  

    // Store S
    store_R_est();

    // Store transition matrix
    if (store_transition_matrix)
      store_P_est();
  }

}


void keyATMhmm::sample_alpha()
{

  // Retrieve start and end indexes of states in documents
  int index_start, index_end;
  for (int r = 0; r < num_states; ++r) {
    if (r == 0) {
      // First state
      // Which time segment correspond to s = 0
      index_start = 0;
      index_end = R_count(r) - 1;

      // Index of documents that belong to s = 0
      states_start(r) = time_doc_start(index_start);
      states_end(r) = time_doc_end(index_end);
      continue;
    }  
    
    index_start = index_end + 1;
    index_end = index_start + R_count(r) - 1;
    states_start(r) = time_doc_start(index_start);
    states_end(r) = time_doc_end(index_end);
  }

  for (int r = 0; r < num_states; ++r) {
    sample_alpha_state(r, states_start(r),
                          states_end(r));  
  }
  
}


void keyATMhmm::sample_alpha_state(int state, int state_start, int state_end)
{

  double start, end, previous_p, new_p, newlikelihood, slice_;
  double store_loglik;
  double newalphallk;

  keep_current_param = alpha;
  topic_ids = sampler::shuffled_indexes(num_topics);
  newalphallk = 0.0;
  int k;

  alpha = alphas.row(state).transpose();  // select alpha to update

  
  for (int i = 0; i < num_topics; ++i) {
    k = topic_ids[i];
    store_loglik = alpha_loglik(k, state_start, state_end);
    start = min_v ; // shrinked with shrinkp()
    end = max_v;  // shrinked with shrinkp()

    previous_p = alpha(k) / (1.0 + alpha(k)); // shrinkp
    slice_ = store_loglik - 2.0 * log(1.0 - previous_p) 
            + log(unif_rand()); // <-- using R random uniform
    
    for (int shrink_time = 0; shrink_time < max_shrink_time; ++shrink_time){
      new_p = sampler::slice_uniform(start, end); // <-- using R function above
      alpha(k) = new_p / (1.0 - new_p); // expandp
    
      newalphallk = alpha_loglik(k, state_start, state_end);
      newlikelihood = newalphallk - 2.0 * log(1.0 - new_p);
    
      if (slice_ < newlikelihood){
        break;
      } else if (previous_p < new_p){
        end = new_p;
      } else if (new_p < previous_p){
        start = new_p;
      } else {
        Rcpp::stop("Something goes wrong in sample_lambda_slice(). Adjust `A_slice`.");
        alpha(k) = keep_current_param(k);
        break;
      }
    }
  }

  // Set new alpha
  alphas.row(state) = alpha.transpose();

}


double keyATMhmm::alpha_loglik(int k, int state_start, int state_end)
{
  double loglik = 0.0;
  double fixed_part = 0.0;

  ndk_a = n_dk.rowwise() + alpha.transpose(); // Use Eigen Broadcasting
  double alpha_sum_val = alpha.sum();


  fixed_part += mylgamma(alpha_sum_val); // first term numerator
  fixed_part -= mylgamma(alpha(k)); // first term denominator
  // Add prior
  if (k < keyword_k) {
    loglik += gammapdfln(alpha(k), eta_1, eta_2);
  } else {
    loglik += gammapdfln(alpha(k), eta_1_regular, eta_2_regular);
  }

  for (int d = state_start; d <= state_end; ++d) {
    loglik += fixed_part;

    // second term numerator
    loglik += mylgamma(ndk_a(d,k));

    // second term denominator
    loglik -= mylgamma(doc_each_len_weighted[d] + alpha_sum_val);
  }
  return loglik;
}


void keyATMhmm::sample_forward()
{ // Calculate Prk (num_doc, num_states)
  double logsum;
  int added;
  double loglik;

  // Prk = MatrixXd::Zero(num_time, num_states);

  for (int t = 0; t < num_time; ++t) {
    if (t == 0) {
      // First time segment should be the first state
      Prk(0, 0) = 1.0;
      continue;
    }  

    // Prepare f in Eq.(6) of Chib (1998)
    for (int r = 0; r < num_states; ++r) {
      // f(y_t | ...) in the numerator
      alpha = alphas.row(r).transpose();
      logfy(r) = polyapdfln(t, alpha);
    }  

    // Prepare Pst
    rt_1l = Prk.row(t-1);  // previous time block
    rt_k = (rt_1l.transpose() * P_est); 
        // p(s_{t} = k), summation is done as matrix calculation
        // Note that P has a lot of 0 elements
        // This is a first term of the numerator in Eq.(6)

    // Format numerator and calculate denominator at the same time
    logsum = 0.0;
    added = 0;
    for (int r = 0; r < num_states; ++r) {
      if (rt_k(r) != 0.0) {
        loglik = log(rt_k(r)) + logfy(r);
        logrt_k(r) = loglik;
        logsum = logsumexp(logsum, loglik, (added == 0));
        added += 1;
      } else {
        logrt_k(r) = 0.0;  // place holder
      }
    }

    for (int r = 0; r < num_states; ++r) {
      if (rt_k(r) != 0.0) {
        Prk(t, r) = exp(logrt_k(r) - logsum);  
      } else {
        Prk(t, r) = 0.0;
      }  
    }

  }

}


double keyATMhmm::polyapdfln(int t, VectorXd &alpha)
{ // Polya distribution: log-likelihood
  double loglik = 0.0;

  int doc_start, doc_end;
  doc_start = time_doc_start(t);  // starting doc index of time segment t
  doc_end = time_doc_end(t);

  for (int d = doc_start; d <= doc_end; ++d) {
    loglik += mylgamma( alpha.sum() ) - mylgamma( doc_each_len_weighted[d] + alpha.sum() );
    for (int k = 0; k < num_topics; ++k) {
      loglik += mylgamma( n_dk(d,k) + alpha(k) ) - mylgamma( alpha(k) );
    }  
  }

  return loglik;
}


void keyATMhmm::sample_backward()
{
  int state_id;

  // sample R_est
  // num_time - 2, because time segment index is (num_time - 1)
  // and we want to start from (time_index - 1)

  R_count = VectorXi::Zero(num_states); // reset counter

  // Last document
  R_est(num_time - 1) = index_states;
  R_count(index_states) += 1;  // last document

  for (int t = (num_time - 2); 0 <= t; --t) {
    state_id = R_est(t + 1);

    state_prob_vec.array() = Prk.row(t).transpose().array() * P_est.col(state_id).array(); 
    state_prob_vec.array() = state_prob_vec.array() / state_prob_vec.sum();

    state_id = sampler::rcat(state_prob_vec, num_states); // new state id
    R_est(t) = state_id;
    R_count(state_id) += 1;
  }

}


void keyATMhmm::sample_P()
{
  double pii;

  // sample P_est
  // iterate until index_state - 2
  for (int r = 0; r <= (num_states - 2); ++r) {
    pii = R::rbeta(R_count(r), 2);  
      // First value is 1 + R_count(s) - 1. 
      // R_count(s) - 1: the number of transitions from state
      // s to state s in the sequence of state
      // ----------------------------------------------------
      // "-1" because the first count in R_count is
      // the transition from s-1 to s
      // prior is Beta(1,1)

    P_est(r, r) = pii;
    P_est(r, r + 1) = 1.0 - pii;
  }
}


void keyATMhmm::store_R_est()
{
  // Store state
  Rcpp::NumericVector state_R = Rcpp::wrap(R_est);
  List R_iter = stored_values["R_iter"];
  R_iter.push_back(state_R);
  stored_values["R_iter"] = R_iter;
}


void keyATMhmm::store_P_est()
{
  // Store state
  Rcpp::NumericMatrix mat_R = Rcpp::wrap(P_est);
  List P_iter = stored_values["P_iter"];
  P_iter.push_back(mat_R);
  stored_values["P_iter"] = P_iter;
}


double keyATMhmm::loglik_total()
{
  double loglik = 0.0;
  int state_id;

  for (int k = 0; k < num_topics; ++k) {
    for (int v = 0; v < num_vocab; ++v) { // word
      loglik += mylgamma(beta + n_s0_kv(k, v)) - mylgamma(beta);
    }



    // word normalization
    loglik += mylgamma( beta * (double)num_vocab ) - mylgamma(beta * (double)num_vocab + n_s0_k(k) );

    if (k < keyword_k) {
      // For keyword topics
      
      // n_s1_kv
      for (SparseMatrix<double,RowMajor>::InnerIterator it(n_s1_kv, k); it; ++it) {
        loglik += mylgamma(beta_s + it.value() / vocab_weights(it.index()) ) - mylgamma(beta_s);
      }
      loglik += mylgamma( beta_s * (double)keywords_num[k] ) - mylgamma(beta_s * (double)keywords_num[k] + n_s1_k(k) );
      
      // Normalization
      loglik += mylgamma( prior_gamma(k, 0) + prior_gamma(k, 1)) - mylgamma( prior_gamma(k, 0)) - mylgamma( prior_gamma(k, 1));

      // s
      loglik += mylgamma( n_s0_k(k) + prior_gamma(k, 1) ) 
                - mylgamma(n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1))
                + mylgamma(n_s1_k(k) + prior_gamma(k, 0) );  
    }
  }


  for (int d = 0; d < num_doc; ++d) {
    // z
    alpha = alphas.row(get_state_index(d)).transpose(); // Doc alpha, column vector  
    
    loglik += mylgamma( alpha.sum() ) - mylgamma( doc_each_len_weighted[d] + alpha.sum() );
    for (int k = 0; k < num_topics; ++k) {
      loglik += mylgamma( n_dk(d,k) + alpha(k) ) - mylgamma( alpha(k) );
    }

  }

  // HMM part
  for (int t = 0; t < num_time; ++t) {
    state_id = R_est(t);
    loglik += log( P_est(state_id, state_id) );
  }

  return loglik;
}


double keyATMhmm::loglik_total_label()
{
  double loglik = 0.0;
  int state_id;

  for (int k = 0; k < num_topics; ++k) {
    for (int v = 0; v < num_vocab; ++v) { // word
      loglik += mylgamma(beta_s0kv(k, v) + n_s0_kv(k, v) ) - mylgamma(beta_s0kv(k, v));
    }

    // word normalization
    loglik += mylgamma( Vbeta_k(k) ) - mylgamma(Vbeta_k(k) + n_s0_k(k) );

    if (k < keyword_k) {
      // For keyword topics
      
      // n_s1_kv
      for (SparseMatrix<double,RowMajor>::InnerIterator it(n_s1_kv, k); it; ++it) {
        loglik += mylgamma(beta_s1kv.coeffRef(k, it.index()) + it.value()) - mylgamma(beta_s1kv.coeffRef(k, it.index()));
      }
      loglik += mylgamma( Lbeta_sk(k) ) - mylgamma(Lbeta_sk(k) + n_s1_k(k) );
      
      // Normalization
      loglik += mylgamma( prior_gamma(k, 0) + prior_gamma(k, 1)) - mylgamma( prior_gamma(k, 0)) - mylgamma( prior_gamma(k, 1));

      // s
      loglik += mylgamma( n_s0_k(k) + prior_gamma(k, 1) ) 
                - mylgamma(n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1))
                + mylgamma(n_s1_k(k) + prior_gamma(k, 0) );  
    }
  }


  for (int d = 0; d < num_doc; ++d) {
    // z
    alpha = alphas.row(get_state_index(d)).transpose(); // Doc alpha, column vector  
    
    loglik += mylgamma( alpha.sum() ) - mylgamma( doc_each_len_weighted[d] + alpha.sum() );
    for (int k = 0; k < num_topics; ++k) {
      loglik += mylgamma( n_dk(d,k) + alpha(k) ) - mylgamma( alpha(k) );
    }

  }

  // HMM part
  for (int t = 0; t < num_time; ++t) {
    state_id = R_est(t);
    loglik += log( P_est(state_id, state_id) );
  }

  return loglik;
}

