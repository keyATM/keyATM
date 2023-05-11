#include "keyATM_multi_base.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;


void keyATMmultibase::read_data_specific()
{
  num_corpus = model_settings["num_corpus"]; //Nr. of corpora
  corpus_id = model_settings["corpus_id"]; // Vector of ints of corpus id for each document
  global_id = model_settings["global_id"]; //List num_corpus of integer vectors, each with global id of document in corpus
  num_doc_all = model_settings["docs_per_corpus"]; //Vector of ints, nr. of documents per corpus
  
  Alpha = Rcpp::as<Eigen::MatrixXd>(priors_list["alpha"]);
  
  estimate_alpha = options_list["estimate_alpha"];
  if (estimate_alpha == 0) {
    store_alpha = 0;
  } else {
    store_alpha = 1;
  }
}


void keyATMmultibase::initialize_specific()
{
  int s_, z_, w_, corpus_, doc_len_;
  n_s0_kv_all.resize(num_corpus);
  n_s0_k_all = MatrixXd::Zero(num_topics, num_corpus);
  vocab_weights_corpus = MatrixXd::Constant(num_corpus, num_vocab, 1.0);
  total_words_corpus.resize(num_corpus);
  
  
  // storage for sufficient statistics and their margins
  n_s1_kv.resize(num_topics, num_vocab);
  n_dk = MatrixXd::Zero(num_doc, num_topics);
  n_s1_k = VectorXd::Zero(num_topics);
  
  
  for(int corpus = 0; corpus < num_corpus; ++corpus){
    //Corpus specific sufficient statistics init
    n_s0_kv_all[corpus] = MatrixXd::Zero(num_topics, num_vocab);
    
    //Corpus specific weights
    //vocab_weights_corpus[corpus] = VectorXd::Constant(num_vocab, 1.0);
    total_words_corpus[corpus] = 0.0;
  }
  
  for (int doc_id = 0; doc_id < num_doc; doc_id++) {
    corpus_ = corpus_id[doc_id];
    doc_w = W[doc_id];
    doc_len_ = doc_each_len[doc_id];
    for (int w_position = 0; w_position < doc_len_; w_position++) {
      w_ = doc_w[w_position];
      vocab_weights_corpus(corpus_, w_)++;
      (total_words_corpus[corpus_])++;
    }
  }
  
  for (int corpus = 0; corpus < num_corpus; ++corpus){
    if (weights_type == "inv-freq" || weights_type == "inv-freq-normalized") {
      // Inverse frequency
      vocab_weights_corpus.row(corpus) = (double)total_words_corpus[corpus] / (vocab_weights_corpus.row(corpus)).array();
    } else if (weights_type == "information-theory" ||
      weights_type == "information-theory-normalized")
    {
      // Information theory
      vocab_weights_corpus.row(corpus) = (vocab_weights_corpus.row(corpus)).array() / (double)total_words_corpus[corpus];
      vocab_weights_corpus.row(corpus) = (vocab_weights_corpus.row(corpus)).array().log();
      vocab_weights_corpus.row(corpus) = -(vocab_weights_corpus.row(corpus)).array() / log(2);
    }
    
    // Do you want to use weights?
    if (use_weights == 0) {
      Rcpp::Rcout << "Not using weights!! Check `options$use_weights`." << std::endl;
      vocab_weights_corpus.row(corpus) = RowVectorXd::Constant(num_vocab, 1.0);
    }
  }
  
  // Normalize weights
  if (weights_type == "inv-freq-normalized" ||
      weights_type == "information-theory-normalized") {
    NumericVector total_weights(num_corpus);
    int doc_len;
    for (int doc_id = 0; doc_id < num_doc; ++doc_id) {
      corpus_ = corpus_id[doc_id];
      doc_w = W[doc_id];
      doc_len = doc_each_len[doc_id];
      for (int w_position = 0; w_position < doc_len; ++w_position) {
        w_ = doc_w[w_position];
        total_weights[corpus_] += vocab_weights_corpus(corpus_, w_);
      }
    }
    for(int corpus = 0; corpus< num_corpus; ++corpus){
     vocab_weights_corpus.row(corpus) = (vocab_weights_corpus.row(corpus)).array() * (double)total_words_corpus[corpus_] / total_weights[corpus_];
    }
  }
  
  
  //
  // Construct (corpus specific) data matrices
  //
  vector<Triplet> trip_s1;  // for a sparse matrix
  total_words_weighted = 0.0;
  double temp;
  
  for (int doc_id = 0; doc_id < num_doc; ++doc_id) {
    corpus_ = corpus_id[doc_id];
    doc_s = S[doc_id], doc_z = Z[doc_id], doc_w = W[doc_id];
    doc_len_ = doc_each_len[doc_id];
    
    for (int w_position = 0; w_position < doc_len_; ++w_position) {
      s_ = doc_s[w_position], z_ = doc_z[w_position], w_ = doc_w[w_position];
      if (s_ == 0){
        n_s0_kv_all[corpus_](z_, w_) += vocab_weights_corpus(corpus_, w_);
        n_s0_k_all(z_, corpus_) += vocab_weights_corpus(corpus_, w_);
      } else {
        trip_s1.push_back(Triplet(z_, w_, vocab_weights_corpus(corpus_, w_)));
        n_s1_k(z_) += vocab_weights_corpus(corpus_, w_);
      }
      n_dk(doc_id, z_) += vocab_weights_corpus(corpus_, w_);
    }
    
    temp = n_dk.row(doc_id).sum();
    doc_each_len_weighted.push_back(temp);
    total_words_weighted += temp;
  }
  n_s1_kv.setFromTriplets(trip_s1.begin(), trip_s1.end());
  
}


void keyATMmultibase::iteration_single(int it)
{ // Single iteration
  int doc_id_;
  int doc_length;
  int corpus_;
  int w_, z_, s_;
  int new_z, new_s;
  int w_position;
  
  doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle
  
  for (int ii = 0; ii < num_doc; ++ii) {
    doc_id_ = doc_indexes[ii];
    corpus_ = corpus_id[doc_id_];
    doc_s = S[doc_id_], doc_z = Z[doc_id_], doc_w = W[doc_id_];
    doc_length = doc_each_len[doc_id_];
    
    token_indexes = sampler::shuffled_indexes(doc_length); //shuffle
    
    // Iterate each word in the document
    for (int jj = 0; jj < doc_length; ++jj) {
      w_position = token_indexes[jj];
      s_ = doc_s[w_position], z_ = doc_z[w_position], w_ = doc_w[w_position];
      
      new_z = sample_z(Alpha.col(corpus_), z_, s_, w_, doc_id_);
      doc_z[w_position] = new_z;
      
      if (keywords[new_z].find(w_) == keywords[new_z].end())	
        continue;
      
      z_ = doc_z[w_position]; // use updated z
      new_s = sample_s(z_, s_, w_, doc_id_);
      doc_s[w_position] = new_s;
    }
    
    Z[doc_id_] = doc_z;
    S[doc_id_] = doc_s;
  }
  sample_parameters(it);
  
}

int keyATMmultibase::sample_z(Ref<VectorXd> alpha, int z, int s,
                              int w, int doc_id)
{
  int new_z;
  double numerator, denominator;
  double sum;
  
  // Get corpus id and choose containers
  int corpus_ = corpus_id[doc_id];
  MatrixXd& n_s0_kv_a = n_s0_kv_all[corpus_];
  Ref<VectorXd> n_s0_k_a = n_s0_k_all.col(corpus_);
  
  // remove data
  if (s == 0) {
    n_s0_kv_a(z, w) -= vocab_weights_corpus(corpus_, w);
    n_s0_k_a(z) -= vocab_weights_corpus(corpus_, w);
  } else if (s == 1) {
    n_s1_kv.coeffRef(z, w) -= vocab_weights_corpus(corpus_, w);
    n_s1_k(z) -= vocab_weights_corpus(corpus_, w);
  } else {
    Rcerr << "Error at sample_z, remove" << std::endl;
  }
  
  n_dk(doc_id, z) -= vocab_weights_corpus(corpus_, w);
  n_dk_noWeight(doc_id, z) -= 1.0;
  
  new_z = -1; // debug
  if (s == 0) {
    for (int k = 0; k < num_topics; ++k) {
      
      numerator = (beta + n_s0_kv_a(k, w)) *
        (n_s0_k_a(k) + prior_gamma(k, 1)) *
        (n_dk(doc_id, k) + alpha(k));
      
      denominator = (Vbeta + n_s0_k_a(k)) *
        (n_s1_k(k) + prior_gamma(k, 0) + n_s0_k_a(k) + prior_gamma(k, 1));
      
      z_prob_vec(k) = numerator / denominator;
    }
    
    sum = z_prob_vec.sum(); // normalize
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample
    
  } else {
    for (int k = 0; k < num_topics; ++k) {
      if (keywords[k].find(w) == keywords[k].end()) {
        z_prob_vec(k) = 0.0;
        continue;
      } else {
        numerator = (beta_s + n_s1_kv.coeffRef(k, w)) *
          (n_s1_k(k) + prior_gamma(k, 0)) *
          (n_dk(doc_id, k) + alpha(k));
        denominator = (Lbeta_sk(k) + n_s1_k(k) ) *
          (n_s1_k(k) + prior_gamma(k, 0) + n_s0_k_a(k) + prior_gamma(k, 1));
        
        z_prob_vec(k) = numerator / denominator;
      }
    }
    
    sum = z_prob_vec.sum();
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample
    
  }
  
  // add back data counts
  if (s == 0) {
    n_s0_kv_a(new_z, w) += vocab_weights_corpus(corpus_, w);
    n_s0_k_a(new_z) += vocab_weights_corpus(corpus_, w);
  } else if (s == 1) {
    n_s1_kv.coeffRef(new_z, w) += vocab_weights_corpus(corpus_, w);
    n_s1_k(new_z) += vocab_weights_corpus(corpus_, w);
  } else {
    Rcerr << "Error at sample_z, add" << std::endl;
  }
  n_dk(doc_id, new_z) += vocab_weights_corpus(corpus_, w);
  n_dk_noWeight(doc_id, new_z) += 1.0;
  
  return new_z;
}



int keyATMmultibase::sample_s(int z, int s, int w, int doc_id)
{
  int new_s;
  double numerator, denominator;
  double s0_prob;
  double s1_prob;
  double sum;
  
  // Get corpus id and choose containers
  int corpus_ = corpus_id[doc_id];
  MatrixXd& n_s0_kv_a = n_s0_kv_all[corpus_];
  Ref<VectorXd> n_s0_k_a = n_s0_k_all.col(corpus_);
  
  // remove data
  if (s == 0) {
    n_s0_kv_a(z, w) -= vocab_weights_corpus(corpus_, w);
    n_s0_k_a(z) -= vocab_weights_corpus(corpus_, w);
  } else {
    n_s1_kv.coeffRef(z, w) -= vocab_weights_corpus(corpus_, w);
    n_s1_k(z) -= vocab_weights_corpus(corpus_, w);
  }
  
  // newprob_s1()
  
  numerator = (beta_s + n_s1_kv.coeffRef(z, w)) *
    ( n_s1_k(z) + prior_gamma(z, 0) );
  denominator = (Lbeta_sk(z) + n_s1_k(z) );
  s1_prob = numerator / denominator;
  
  // newprob_s0()
  numerator = (beta + n_s0_kv_a(z, w)) *
    (n_s0_k_a(z) + prior_gamma(z, 1));
  
  denominator = (Vbeta + n_s0_k_a(z) );
  s0_prob = numerator / denominator;
  
  // Normalize
  sum = s0_prob + s1_prob;
  
  s1_prob = s1_prob / sum;
  new_s = R::runif(0,1) <= s1_prob;  //new_s = Bern(s0_prob, s1_prob);
  
  // add back data counts
  if (new_s == 0) {
    n_s0_kv_a(z, w) += vocab_weights_corpus(corpus_, w);
    n_s0_k_a(z) += vocab_weights_corpus(corpus_, w);
  } else {
    n_s1_kv.coeffRef(z, w) += vocab_weights_corpus(corpus_, w);
    n_s1_k(z) += vocab_weights_corpus(corpus_, w);
  }
  
  return new_s;
}





void keyATMmultibase::sample_parameters(int it)
{
  if (estimate_alpha){
    for(int corpus = 0; corpus < num_corpus; ++corpus){
      sample_alpha(corpus);
    }
  }
  
  // Store alpha
  if (store_alpha) {
    int r_index = it + 1;
    if (r_index % thinning == 0 || r_index == 1 || r_index == iter) {
      SEXP tmp = Rcpp::wrap(Alpha);
      NumericMatrix alpha_rmat(tmp);
      List alpha_iter = stored_values["alpha_iter"];
      alpha_iter.push_back(alpha_rmat);
      stored_values["alpha_iter"] = alpha_iter;  
    }
  }
}


void keyATMmultibase::sample_alpha(int corpus_)
{
  
  double start, end, previous_p, new_p, newlikelihood, slice_;
  Ref<VectorXd> alpha_c = Alpha.col(corpus_);
  keep_current_param = alpha_c.transpose();
  topic_ids = sampler::shuffled_indexes(num_topics);
  newalphallk = 0.0;
  int k;
  
  for (int i = 0; i < num_topics; ++i) {
    k = topic_ids[i];
    store_loglik = alpha_loglik(k, corpus_);
    start = min_v ; // shrinked with shrinkp()
    end = max_v;  // shrinked with shrinkp()
    
    previous_p = alpha_c(k) / (1.0 + alpha_c(k)); // shrinkp
    slice_ = store_loglik - 2.0 * log(1.0 - previous_p) 
      + log(unif_rand()); // <-- using R random uniform
    
    for (int shrink_time = 0; shrink_time < max_shrink_time; ++shrink_time) {
      new_p = sampler::slice_uniform(start, end); // <-- using R function above
      alpha_c(k) = new_p / (1.0 - new_p); // expandp
      
      newalphallk = alpha_loglik(k, corpus_);
      newlikelihood = newalphallk - 2.0 * log(1.0 - new_p);
      
      if (slice_ < newlikelihood) {
        break;
      } else if (previous_p < new_p) {
        end = new_p;
      } else if (new_p < previous_p) {
        start = new_p;
      } else {
        Rcerr << "Something went wrong in sample_lambda_slice()." << std::endl;
        alpha_c(k) = keep_current_param(k);
        break;
      }
    }
  }
}


double keyATMmultibase::alpha_loglik(int k, int corpus_)
{
  double loglik = 0.0;
  double fixed_part = 0.0;
  int num_doc_a = num_doc_all[corpus_];
  const VectorXd& alpha_c = Alpha.col(corpus_);
  
  
  ndk_a = n_dk.rowwise() + alpha_c.transpose(); // Use Eigen Broadcasting
  double alpha_sum_val = alpha_c.sum();
  
  
  fixed_part += mylgamma(alpha_sum_val); // first term numerator
  fixed_part -= mylgamma(alpha_c(k)); // first term denominator
  // Add prior
  if (k < keyword_k) {
    loglik += gammapdfln(alpha_c(k), eta_1, eta_2);
  } else {
    loglik += gammapdfln(alpha_c(k), eta_1_regular, eta_2_regular);
  }
  
  for (int d = 0; d < num_doc_a; ++d) {
    g_doc_id = global_id[corpus_][d];
    loglik += fixed_part;
    
    // second term numerator
    loglik += mylgamma(ndk_a(g_doc_id, k));
    
    // second term denominator
    loglik -= mylgamma(doc_each_len_weighted[g_doc_id] + alpha_sum_val);
  }
  
  return loglik;
}


double keyATMmultibase::loglik_total()
{
  int doc_id_ = 0;
  double loglik = 0.0;
  double fixed_part = 0.0;
  for(int c = 0; c < num_corpus; ++c){
    Ref<VectorXd> alpha_c = Alpha.col(c);
    for (int k = 0; k < num_topics; ++k) {
      for (int v = 0; v < num_vocab; ++v) { // word
        loglik += mylgamma(beta + n_s0_kv_all[c](k, v) ) - mylgamma(beta);
      }
      
      // word normalization
      loglik += mylgamma( beta * (double)num_vocab ) - mylgamma(beta * (double)num_vocab + n_s0_k_all(k, c) );
      
      if (k < keyword_k) {
        // For keyword topics
        
        // n_s1_kv
        for (SparseMatrix<double,RowMajor>::InnerIterator it(n_s1_kv, k); it; ++it) {
          loglik += mylgamma(beta_s + it.value()) - mylgamma(beta_s);
        }
        loglik += mylgamma( beta_s * (double)keywords_num[k] ) - mylgamma(beta_s * (double)keywords_num[k] + n_s1_k(k) );
        
        // Normalization
        loglik += mylgamma( prior_gamma(k, 0) + prior_gamma(k, 1)) - mylgamma( prior_gamma(k, 0)) - mylgamma( prior_gamma(k, 1));
        
        // s
        loglik += mylgamma( n_s0_k_all(k, c) + prior_gamma(k, 1) ) 
          - mylgamma(n_s1_k(k) + prior_gamma(k, 0) + n_s0_k_all(k, c) + prior_gamma(k, 1))
          + mylgamma(n_s1_k(k) + prior_gamma(k, 0));  
      }
    }
    
    // z
    fixed_part = alpha_c.sum();
    for (int d = 0; d < num_doc_all[c]; ++d) {
      doc_id_ = global_id[c][d];
      loglik += mylgamma( fixed_part ) - mylgamma( doc_each_len_weighted[doc_id_] + fixed_part );
      
      for (int k = 0; k < num_topics; ++k) {
        loglik += mylgamma( n_dk(doc_id_,k) + alpha_c(k) ) - mylgamma( alpha_c(k) );
      }
    }
  }
  
  return loglik;
}

double keyATMmultibase::loglik_total_label() // NEEDS TO CHANGE IF LABELS ARE IMPLEMENTED
{
  int doc_id_ = 0;
  double loglik = 0.0;
  double fixed_part = 0.0;
  
  
  for(int c = 0; c < num_corpus; ++c){
    const VectorXd& alpha_c= Alpha.col(c);
    for (int k = 0; k < num_topics; ++k) {
      for (int v = 0; v < num_vocab; ++v) { // word
        loglik += mylgamma(beta_s0kv(k, v) + n_s0_kv_all[c](k, v) ) - mylgamma(beta_s0kv(k, v));
      }
      
      // word normalization
      loglik += mylgamma( Vbeta_k(k) ) - mylgamma(Vbeta_k(k) + n_s0_k_all(k, c) );
      
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
        loglik += mylgamma( n_s0_k_all(k, c) + prior_gamma(k, 1) ) 
          - mylgamma(n_s1_k(k) + prior_gamma(k, 0) + n_s0_k_all(k, c) + prior_gamma(k, 1))
          + mylgamma(n_s1_k(k) + prior_gamma(k, 0));  
      }
    }
    
    // z
    fixed_part = alpha_c.sum();
    for (int d = 0; d < num_doc; ++d) {
      doc_id_ = global_id[c][d];
      loglik += mylgamma( fixed_part ) - mylgamma( doc_each_len_weighted[doc_id_] + fixed_part );
      
      for (int k = 0; k < num_topics; ++k) {
        loglik += mylgamma( n_dk(doc_id_,k) + alpha_c(k) ) - mylgamma( alpha_c(k) );
      }
    }
  }
  
  return loglik;
}


void keyATMmultibase::store_pi_iter(int r_index)
{
  List pi_vectors = stored_values["pi_vectors"];
  // calculate
  VectorXd numer = n_s1_k.array() + prior_gamma.col(0).array();
  VectorXd n_s0_k_row = n_s0_k_all.rowwise().sum();
  VectorXd denom = n_s0_k_row.array() + prior_gamma.col(1).array() + numer.array();
  VectorXd pi = numer.array() / denom.array();
  
  // store
  NumericVector pi_vector = Rcpp::wrap(pi);
  pi_vectors.push_back(pi_vector);
  stored_values["pi_vectors"] = pi_vectors;
}


List keyATMmultibase::return_model()
{
  // Return output to R
  
  if (use_labels) {
    // Return prior to use in R
    NumericMatrix R_betas0 = Rcpp::wrap(beta_s0kv);
    SEXP R_betas1 = Rcpp::wrap(beta_s1kv);
    
    priors_list.push_back(R_betas0, "beta_s0");
    priors_list.push_back(R_betas1, "beta_s1");
    model["priors"] = priors_list;
  }
  
  // Vocabulary weights
  SEXP tmp = wrap(vocab_weights_corpus.colwise().mean());
  NumericVector vocab_weights_R(tmp);
  stored_values["vocab_weights"] = vocab_weights_R;
  model["stored_values"] = stored_values;
  
  model["stored_values"] = stored_values;
  
  return model;
}
