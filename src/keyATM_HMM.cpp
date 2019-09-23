#include "keyATM_HMM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */

keyATMhmm::keyATMhmm(List model_, const int iter_, const int output_per_) :
  keyATMbase(model_, iter_, output_per_) // pass to parent!
{
  // Constructor
  read_data();
  initialize();
  iteration();
}


void keyATMhmm::read_data_specific()
{
  num_states = model["num_states"];
  index_states = num_states - 1;
}


void keyATMhmm::initialize_specific()
{
  // Initialize Psk
  Psk = MatrixXd::Zero(num_doc, num_states);

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

  for(int j=0; j<num_doc; j++){
    u = R::runif(0, 1);
    for(int i=0; i<num_states; i++){
      if(u < S_est_temp(i)){
        index = i;
        break;
      }
    }
    S_est_num(index) += 1;
  }


  S_est = VectorXi::Zero(num_doc);
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


}


void keyATMhmm::iteration_single(int &it)
{
  doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle

  for (int ii = 0; ii < num_doc; ii++){
    doc_id_ = doc_indexes[ii];
    doc_x = X[doc_id_], doc_z = Z[doc_id_], doc_w = W[doc_id_];
    doc_length = doc_each_len[doc_id_];

    alpha = alphas.row(S_est(doc_id_)).transpose(); // select alpha for this document
    
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

  sample_parameters(int &it);
}

void keyATMhmm::verbose_special(int &r_index){
  // If there is anything special to show, write here.
  if(floor(iter/2) < r_index-1)
    cout << "  Sampled S: " << S_count.transpose() << endl;

  // Store S
  store_S_est();
}

void keyATMhmm::sample_parameters(int &it)
{
  if(floor(iter/2) < it)
    return;

  // alpha
  sample_alpha();
  
  // HMM
  sample_forward();  // calculate Psk
  sample_backward();  // sample S_est
  sample_P();  // sample P_est
}


void keyATMhmm::sample_alpha()
{

  // Retrieve start and end indexes of states in documents
  for(int s=0; s<num_states; s++){
    if(s == 0){
      // First state
      states_start(s) = 0;
      states_end(s) = S_count(s) - 1;
      continue;
    }  
    
    states_start(s) = states_end(s-1) + 1;
    states_end(s) = states_start(s) + S_count(s) - 1;
  }

  
  for(int s=0; s<num_states; s++){
    sample_alpha_state(s, states_start(s),
                          states_end(s));  
  }
  
  // Store alphas
  Rcpp::NumericMatrix alphas_R = Rcpp::wrap(alphas);
  List alpha_iter = model["alpha_iter"];
  alpha_iter.push_back(alphas_R);
  model["alpha_iter"] = alpha_iter;
  
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

  
  for(int i = 0; i < num_topics; i++){
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
  alphas.row(state) = alpha.transpose();           // Use this line!!

}


double keyATMhmm::alpha_loglik(int &state_start, int &state_end)
{

  loglik = 0.0;
  fixed_part = 0.0;
  ndk_a = n_dk.rowwise() + alpha.transpose(); // Use Eigen Broadcasting
  alpha_sum_val = alpha.sum();


  fixed_part += mylgamma(alpha_sum_val); // first term numerator
  for(int k = 0; k < num_topics; k++){
    fixed_part -= mylgamma(alpha(k)); // first term denominator
    // Add prior
    if(k < k_seeded){
      loglik += gammapdfln(alpha(k), eta_1, eta_2);
    }else{
      loglik += gammapdfln(alpha(k), eta_1_regular, eta_2_regular);
    }

  }
  for(int d = state_start; d <= state_end; d++){
    loglik += fixed_part;
    // second term numerator
    for(int k = 0; k < num_topics; k++){
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

  for(int d=0; d<num_doc; d++){
    if(d == 0){
      // First document should be the first state
      Psk(0, 0) = 1.0;
      continue;
    }  

    // Prepare f in Eq.(6) of Chib (1998)
    for(int s=0; s<num_states; s++){
      alpha = alphas.row(s).transpose();
      logfy(s) = polyapdfln(d, alpha);
    }  

    // Prepare Pst
    st_1l = Psk.row(d-1);  // previous observation
    st_k = (st_1l.transpose() * P_est);

    // Format numerator and calculate denominator at the same time
    logsum = 0.0;
    added = 0;
    for(int s=0; s<num_states; s++){
      if(st_k(s) != 0.0){
        loglik = log(st_k(s)) + logfy(s);
        logst_k(s) = loglik;
        logsum = logsumexp(logsum, loglik, (added == 0));
        added += 1;
      }else{
        logst_k(s) = 0.0;  // place holder
      }
    }

    for(int s=0; s<num_states; s++){
      if(st_k(s) != 0.0){
        Psk(d, s) = exp( logst_k(s) - logsum );  
      }else{
        Psk(d, s) = 0.0;  
      }  
    }

  }

}

double keyATMhmm::polyapdfln(int &d, VectorXd &alpha)
{ // Polya distribution log-likelihood
  loglik = 0.0;

  loglik += mylgamma( alpha.sum() ) - mylgamma( doc_each_len[d] + alpha.sum() );
  for (int k = 0; k < num_topics; k++){
    loglik += mylgamma( n_dk(d,k) + alpha(k) ) - mylgamma( alpha(k) );
  }

  return loglik;
}


void keyATMhmm::sample_backward()
{
  // sample S_est
  // num_doc - 2, because doc_index is (num_doc - 1)
  // and we want to start from (doc_index - 1)

  S_count = VectorXi::Zero(num_states); // reset counter

  // Last document
  S_est(num_doc-1) = index_states;
  S_count(index_states) += 1;  // last document

  for(int d=(num_doc-2); 0<= d; --d){
    state_id = S_est(d+1);

    state_prob_vec.array() = Psk.row(d).transpose().array() * P_est.col(state_id).array(); 
    state_prob_vec.array() = state_prob_vec.array() / state_prob_vec.sum();

    state_id = sampler::rcat(state_prob_vec, num_states); // new state id
    S_est(d) = state_id;
    S_count(state_id) += 1;
  }

}


void keyATMhmm::sample_P()
{
  // sample P_est
  // iterate until index_state - 2
  for(int s=0; s<=(num_states-2); ++s){
    pii = R::rbeta(1+S_count(s), 1);  

    P_est(s, s) = pii;
    P_est(s, s+1) = 1.0 - pii;
  }
}


void keyATMhmm::store_S_est()
{
  // Store state
  Rcpp::NumericVector state_R = Rcpp::wrap(S_est);
  List S_iter = model["S_iter"];
  S_iter.push_back(state_R);
  model["S_iter"] = S_iter;
}


double keyATMhmm::loglik_total()
{
  loglik = 0.0;
  for (int k = 0; k < num_topics; k++){
    for (int v = 0; v < num_vocab; v++){ // word
      loglik += mylgamma(beta + n_x0_kv(k, v) / vocab_weights(v) ) - mylgamma(beta);
      // loglik += mylgamma(beta_s + n_x1_kv.coeffRef(k, v) / vocab_weights(v) ) - mylgamma(beta_s);
    }

    // n_x1_kv
    for (SparseMatrix<double,RowMajor>::InnerIterator it(n_x1_kv, k); it; ++it){
      loglik += mylgamma(beta_s + it.value() / vocab_weights(it.index()) ) - mylgamma(beta_s);
    }

    // word normalization
    loglik += mylgamma( beta * (double)num_vocab ) - mylgamma(beta * (double)num_vocab + n_x0_k_noWeight(k) );
    loglik += mylgamma( beta_s * (double)num_vocab ) - mylgamma(beta_s * (double)num_vocab + n_x1_k_noWeight(k) );
    // x
    loglik += mylgamma( n_x0_k_noWeight(k) + x_prior(k, 1) ) - mylgamma(n_x1_k_noWeight(k) + x_prior(k, 0) + n_x0_k_noWeight(k) + x_prior(k, 1))
      + mylgamma( n_x1_k_noWeight(k) + x_prior(k, 0) ) ;
    
    // x normalization
    loglik += mylgamma(x_prior(k, 0) + x_prior(k, 1)) - mylgamma(x_prior(k, 0)) - mylgamma(x_prior(k, 1));
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

