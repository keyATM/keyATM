#include "keyATM_HMM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */

keyATMhmm::keyATMhmm(List model_, const int iter_, const int output_per_) :
  keyATMbase(model_, iter_, output_per_) // pass to parent!
{

}


void keyATMhmm::read_data_specific()
{
  model_settings = model["model_settings"];
  num_states = model_settings["num_states"];
  index_states = num_states - 1;

  prior_gamma = MatrixXd::Zero(num_topics, 2);
  NumericMatrix RMatrix = priors_list["gamma"];
  prior_gamma = Rcpp::as<Eigen::MatrixXd>(RMatrix);
  beta_s = priors_list["beta_s"];

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
  for (int d = 0; d < num_doc; d++) {
     index = time_index[d]; 
    if (index != index_prev) {
      time_doc_start[store_index] = d;
      index_prev = index;
      store_index += 1;
    }
  }

  for (int s = 0; s < num_time-1; s++) {
    time_doc_end(s) = time_doc_start(s + 1) - 1;  
  }
  time_doc_end(num_time-1) = num_doc-1;

  store_transition_matrix = options_list["store_transition_matrix"];
}


void keyATMhmm::initialize_specific()
{
  // Initialize Psk
  Psk = MatrixXd::Zero(num_time, num_states);

  // Initialize S_est
  // Use multinomial distribution (with flat probability)
  // to decide the number of each state
  // and push it into S_est.
  VectorXi S_est_num = VectorXi::Constant(num_states, 1);
  VectorXd S_est_temp = VectorXd::Zero(num_states);
  double cumulative = 1.0 / num_states;
  double u;
  int index;
  for (int i = 0; i < num_states; i++) {
    S_est_temp(i) = cumulative * (i + 1);
  }

  for (int j = 0; j < num_time-num_states; j++) {
    u = R::runif(0, 1);
    for (int i = 0; i < num_states; i++) {
      if (u < S_est_temp(i)) {
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
  for (int i = 0; i < num_states; i++) {
    count = S_est_num(i);
    for (int j = 0; j < count; j++) {
      S_est(index) = i;
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

  // cout << S_est_num.transpose() << endl;  //debug
  // cout << S_est.transpose() << endl;  //debug
  // cout << P_est << endl;  //debug

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


int keyATMhmm::get_state_index(const int &doc_id)
{
  // Which time segment the document belongs to
  int t;
  for (t = 0; t < num_time; t++) {
    if (time_doc_start(t) <= doc_id && doc_id <= time_doc_end(t)) {
      break;  
    }
  }
  return(S_est(t));
}


void keyATMhmm::iteration_single(int &it)
{
  doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle

  for (int ii = 0; ii < num_doc; ii++){
    doc_id_ = doc_indexes[ii];
    doc_x = X[doc_id_], doc_z = Z[doc_id_], doc_w = W[doc_id_];
    doc_length = doc_each_len[doc_id_];

    alpha = alphas.row(get_state_index(doc_id_)).transpose(); // select alpha for this document
    
    token_indexes = sampler::shuffled_indexes(doc_length); //shuffle
    
    // Iterate each word in the document
    for (int jj = 0; jj < doc_length; jj++){
      w_position = token_indexes[jj];
      x_ = doc_x[w_position], z_ = doc_z[w_position], w_ = doc_w[w_position];
    
      new_z = sample_z(alpha, z_, x_, w_, doc_id_);
      doc_z[w_position] = new_z;
    
  
      z_ = doc_z[w_position]; // use updated z
      new_x = sample_x(alpha, z_, x_, w_, doc_id_);
      doc_x[w_position] = new_x;
    }
    
    Z[doc_id_] = doc_z;
    X[doc_id_] = doc_x;
  }

  sample_parameters(it);
}

void keyATMhmm::verbose_special(int &r_index)
{
  // If there is anything special to show, write here.
}

void keyATMhmm::sample_parameters(int &it)
{
  // alpha
  sample_alpha();
  
  // HMM
  sample_forward();  // calculate Psk
  sample_backward();  // sample S_est
  sample_P();  // sample P_est  

  // Store alpha and S
  int r_index = it + 1;
  if (r_index % thinning == 0 || r_index == 1 || r_index == iter) {
    Rcpp::NumericMatrix alphas_R = Rcpp::wrap(alphas);
    List alpha_iter = stored_values["alpha_iter"];
    alpha_iter.push_back(alphas_R);
    stored_values["alpha_iter"] = alpha_iter;  

    // Store S
    store_S_est();

    // Store transition matrix
    if (store_transition_matrix)
      store_P_est();
  }

}


void keyATMhmm::sample_alpha()
{

  // Retrieve start and end indexes of states in documents
  int index_start, index_end;
  for (int s = 0; s < num_states; s++) {
    if (s == 0) {
      // First state
      // Which time segment correspond to s = 0
      index_start = 0;
      index_end = S_count(s) - 1;

      // Index of documents that belong to s = 0
      states_start(s) = time_doc_start(index_start);
      states_end(s) = time_doc_end(index_end);
      continue;
    }  
    
    index_start = index_end + 1;
    index_end = index_start + S_count(s) - 1;
    states_start(s) = time_doc_start(index_start);
    states_end(s) = time_doc_end(index_end);
  }

  // // Debug
  // cout << S_est.transpose() << endl;
  // cout << S_count.transpose() << endl;
  // cout << states_start.transpose() << endl;
  // cout << states_end.transpose() << endl;

  for (int s = 0; s < num_states; s++) {
    sample_alpha_state(s, states_start(s),
                          states_end(s));  
  }
  
}


void keyATMhmm::sample_alpha_state(int &state, int &state_start, int &state_end)
{

  // start, end, previous_p, new_p, newlikelihood, slice_;
  keep_current_param = alpha;
  topic_ids = sampler::shuffled_indexes(num_topics);
  newalphallk = 0.0;
  int k;

  alpha = alphas.row(state).transpose();  // select alpha to update
  store_loglik = alpha_loglik(state_start, state_end);

  
  for (int i = 0; i < num_topics; i++) {
    k = topic_ids[i];
    start = min_v / (1.0 + min_v); // shrinkp
    end = 1.0;
    // end = shrinkp(max_v);
    previous_p = alpha(k) / (1.0 + alpha(k)); // shrinkp
    slice_ = store_loglik - 2.0 * log(1.0 - previous_p) 
            + log(unif_rand()); // <-- using R random uniform
    
    for (int shrink_time = 0; shrink_time < max_shrink_time; shrink_time++){
      new_p = sampler::slice_uniform(start, end); // <-- using R function above
      alpha(k) = new_p / (1.0 - new_p); // expandp
    
      newalphallk = alpha_loglik(state_start, state_end);
      newlikelihood = newalphallk - 2.0 * log(1.0 - new_p);
    
      if (slice_ < newlikelihood){
        store_loglik = newalphallk;
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


double keyATMhmm::alpha_loglik(int &state_start, int &state_end)
{

  loglik = 0.0;
  fixed_part = 0.0;
  ndk_a = n_dk.rowwise() + alpha.transpose(); // Use Eigen Broadcasting
  alpha_sum_val = alpha.sum();


  fixed_part += mylgamma(alpha_sum_val); // first term numerator
  for (int k = 0; k < num_topics; k++) {
    fixed_part -= mylgamma(alpha(k)); // first term denominator
    // Add prior
    if (k < keyword_k) {
      loglik += gammapdfln(alpha(k), eta_1, eta_2);
    } else {
      loglik += gammapdfln(alpha(k), eta_1_regular, eta_2_regular);
    }

  }
  for (int d = state_start; d <= state_end; d++) {
    loglik += fixed_part;
    // second term numerator
    for (int k = 0; k < num_topics; k++) {
      loglik += mylgamma(ndk_a(d,k));
    }
    // second term denominator
    loglik -= mylgamma(doc_each_len[d] + alpha_sum_val);

  }
  return loglik;

}


void keyATMhmm::sample_forward()
{ // Calculate Psk (num_doc, num_states)

  // Psk = MatrixXd::Zero(num_doc, num_states);

  for (int t = 0; t < num_time; t++) {
    if (t == 0) {
      // First time segment should be the first state
      Psk(0, 0) = 1.0;
      continue;
    }  

    // Prepare f in Eq.(6) of Chib (1998)
    for (int s = 0; s < num_states; s++) {
      alpha = alphas.row(s).transpose();
      logfy(s) = polyapdfln(t, alpha);
    }  

    // Prepare Pst
    st_1l = Psk.row(t-1);  // previous observation
    st_k = (st_1l.transpose() * P_est);

    // Format numerator and calculate denominator at the same time
    logsum = 0.0;
    added = 0;
    for (int s = 0; s < num_states; s++) {
      if (st_k(s) != 0.0) {
        loglik = log(st_k(s)) + logfy(s);
        logst_k(s) = loglik;
        logsum = logsumexp(logsum, loglik, (added == 0));
        added += 1;
      } else {
        logst_k(s) = 0.0;  // place holder
      }
    }

    for (int s = 0; s < num_states; s++) {
      if (st_k(s) != 0.0) {
        Psk(t, s) = exp( logst_k(s) - logsum );  
      } else {
        Psk(t, s) = 0.0;  
      }  
    }

  }

}


double keyATMhmm::polyapdfln(int &t, VectorXd &alpha)
{ // Polya distribution log-likelihood
  loglik = 0.0;

  int doc_start, doc_end;
  doc_start = time_doc_start(t);  // starting doc index of time segment t
  doc_end = time_doc_end(t);

  for (int d = doc_start; d <= doc_end; d++) {
    loglik += mylgamma( alpha.sum() ) - mylgamma( doc_each_len[d] + alpha.sum() );
    for (int k = 0; k < num_topics; k++){
      loglik += mylgamma( n_dk(d,k) + alpha(k) ) - mylgamma( alpha(k) );
    }  
  }

  return loglik;
}


void keyATMhmm::sample_backward()
{
  // sample S_est
  // num_time - 2, because time segment index is (num_time - 1)
  // and we want to start from (time_index - 1)

  S_count = VectorXi::Zero(num_states); // reset counter

  // Last document
  S_est(num_time-1) = index_states;
  S_count(index_states) += 1;  // last document

  for (int t=(num_time-2); 0<= t; --t) {
    state_id = S_est(t + 1);

    state_prob_vec.array() = Psk.row(t).transpose().array() * P_est.col(state_id).array(); 
    state_prob_vec.array() = state_prob_vec.array() / state_prob_vec.sum();

    state_id = sampler::rcat(state_prob_vec, num_states); // new state id
    S_est(t) = state_id;
    S_count(state_id) += 1;
  }

}


void keyATMhmm::sample_P()
{
  // sample P_est
  // iterate until index_state - 2
  for (int s = 0; s <= (num_states - 2); ++s) {
    pii = R::rbeta(1 + S_count(s), 1);  

    P_est(s, s) = pii;
    P_est(s, s + 1) = 1.0 - pii;
  }
}


void keyATMhmm::store_S_est()
{
  // Store state
  Rcpp::NumericVector state_R = Rcpp::wrap(S_est);
  List S_iter = stored_values["S_iter"];
  S_iter.push_back(state_R);
  stored_values["S_iter"] = S_iter;
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
  loglik = 0.0;
  for (int k = 0; k < num_topics; k++){
    for (int v = 0; v < num_vocab; v++){ // word
      loglik += mylgamma(beta + n_x0_kv(k, v) / vocab_weights(v) ) - mylgamma(beta);
      // loglik += mylgamma(beta_s + n_x1_kv.coeffRef(k, v) / vocab_weights(v) ) - mylgamma(beta_s);
    }



    // word normalization
    loglik += mylgamma( beta * (double)num_vocab ) - mylgamma(beta * (double)num_vocab + n_x0_k_noWeight(k) );

    if (k < keyword_k) {
      // For keyword topics
      
      // n_x1_kv
      for (SparseMatrix<double,RowMajor>::InnerIterator it(n_x1_kv, k); it; ++it){
        loglik += mylgamma(beta_s + it.value() / vocab_weights(it.index()) ) - mylgamma(beta_s);
      }
      loglik += mylgamma( beta_s * (double)keywords_num[k] ) - mylgamma(beta_s * (double)keywords_num[k] + n_x1_k_noWeight(k) );
      
      // Normalization
      loglik += mylgamma( prior_gamma(k, 0) + prior_gamma(k, 1)) - mylgamma( prior_gamma(k, 0)) - mylgamma( prior_gamma(k, 1));

      // x
      loglik += mylgamma( n_x0_k_noWeight(k) + prior_gamma(k, 1) ) 
                -  mylgamma(n_x1_k_noWeight(k) + prior_gamma(k, 0) + n_x0_k_noWeight(k) + prior_gamma(k, 1))
                + mylgamma( n_x1_k_noWeight(k) + prior_gamma(k, 0) );  
    }
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
  for (int t = 0; t < num_time; t++) {
    state_id = S_est(t);
    loglik += log( P_est(state_id, state_id) );
  }

  return loglik;
}

