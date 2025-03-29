#include "keyATM_multibase.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

void keyATMmultibase::read_data_specific()
{
  num_corpus = model_settings["num_corpus"]; // Nr. of corpora
  corpus_id = model_settings["corpus_id"]; // Vector of ints of corpus id for each document

  C = model["C"];

  NumericMatrix RMatrix = priors_list["omega"];
  prior_omega = Rcpp::as<Eigen::MatrixXd>(RMatrix);

  beta_c = priors_list["beta_c"];
}

void keyATMmultibase::create_sufficient_stats()
{
  Vbeta_c = double(num_vocab) * beta_c;
  Vbeta_s = double(num_vocab) * beta_s;

  int s_, z_, w_, c_, corpus_, doc_len_;
  n_dk = MatrixXd::Zero(num_doc, num_topics);

  n_s0_c0_kv_all.resize(num_corpus);
  
  n_s0_c0_k_all = MatrixXd::Zero(num_topics, num_corpus);
  n_s0_c1_k = VectorXd::Zero(num_topics);
  n_s0_c1_kv = MatrixXd::Zero(num_topics, num_vocab);
  n_c0_k = VectorXd::Zero(num_topics);
  n_c1_k = VectorXd::Zero(num_topics);
  n_c0_k_all = MatrixXd::Zero(num_topics, num_corpus);
  n_c1_k_all = MatrixXd::Zero(num_topics, num_corpus);
  
  n_s0_kv = MatrixXd::Zero(num_topics, num_vocab);
  n_s1_kv_multi.resize(num_topics, num_vocab);
  n_s0_k = VectorXd::Zero(num_topics);
  n_s1_k = VectorXd::Zero(num_topics);

  vocab_weights_corpus = MatrixXd::Constant(num_corpus, num_vocab, 1.0);
  total_words_corpus.resize(num_corpus);

  for (int corpus = 0; corpus < num_corpus; ++corpus) {
    n_s0_c0_kv_all[corpus] = MatrixXd::Zero(num_topics, num_vocab);
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

  for (int corpus = 0; corpus < num_corpus; ++corpus) {
    if (weights_type == "inv-freq" || weights_type == "inv-freq-normalized") {
      vocab_weights_corpus.row(corpus) = (double)total_words_corpus[corpus] / (vocab_weights_corpus.row(corpus)).array();
    } else if (weights_type == "information-theory" || weights_type == "information-theory-normalized") {
      vocab_weights_corpus.row(corpus) = (vocab_weights_corpus.row(corpus)).array() / (double)total_words_corpus[corpus];
      vocab_weights_corpus.row(corpus) = (vocab_weights_corpus.row(corpus)).array().log();
      vocab_weights_corpus.row(corpus) = -(vocab_weights_corpus.row(corpus)).array() / log(2);
    }

    if (use_weights == 0) {
      Rcpp::Rcout << "Not using weights!! Check `options$use_weights`." << std::endl;
      vocab_weights_corpus.row(corpus) = RowVectorXd::Constant(num_vocab, 1.0);
    }
  }

  if (weights_type == "inv-freq-normalized" || weights_type == "information-theory-normalized") {
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
      for (int corpus = 0; corpus < num_corpus; ++corpus) {
        vocab_weights_corpus.row(corpus) = (vocab_weights_corpus.row(corpus)).array() * (double)total_words_corpus[corpus_] / total_weights[corpus_];
      }
    }
  }

  //
  // Construct (corpus specific) data matrices
  //
  vector<Triplet> trip_s1_multi;  // for a sparse matrix
  total_words_weighted = 0.0;
  double temp;
  
  for (int doc_id = 0; doc_id < num_doc; ++doc_id) {
    corpus_ = corpus_id[doc_id];
    doc_s = S[doc_id], doc_z = Z[doc_id], doc_w = W[doc_id], doc_c = C[doc_id];
    doc_len_ = doc_each_len[doc_id];
    
    for (int w_position = 0; w_position < doc_len_; ++w_position) {
      s_ = doc_s[w_position], z_ = doc_z[w_position], w_ = doc_w[w_position], c_ = doc_c[w_position];
      if (s_ == 0){
        if (c_ == 0) {
          n_s0_c0_kv_all[corpus_](z_, w_) += vocab_weights_corpus(corpus_, w_);
          n_c0_k(z_) += vocab_weights_corpus(corpus_, w_);
          n_s0_c0_k_all(z_, corpus_) += vocab_weights_corpus(corpus_, w_);
          n_c0_k_all(z_, corpus_) += vocab_weights_corpus(corpus_, w_);
        } else {
          n_c1_k(z_) += vocab_weights_corpus(corpus_, w_);
          n_s0_c1_k(z_) += vocab_weights_corpus(corpus_, w_);
          n_s0_c1_kv(z_, w_) += vocab_weights_corpus(corpus_, w_);
          n_c1_k_all(z_, corpus_) += vocab_weights_corpus(corpus_, w_);
        }
        n_s0_kv(z_, w_) += vocab_weights_corpus(corpus_, w_);
        n_s0_k(z_) += vocab_weights_corpus(corpus_, w_);
      } else {
        trip_s1_multi.push_back(Triplet(z_, w_, vocab_weights_corpus(corpus_, w_)));
        n_s1_k(z_) += vocab_weights_corpus(corpus_, w_);
        if (c_ == 0) {
          n_c0_k(z_) += vocab_weights_corpus(corpus_, w_);
          n_c0_k_all(z_, corpus_) += vocab_weights_corpus(corpus_, w_);
        } else {
          n_c1_k(z_) += vocab_weights_corpus(corpus_, w_);
          n_c1_k_all(z_, corpus_) += vocab_weights_corpus(corpus_, w_);
        }
      }
      n_dk(doc_id, z_) += vocab_weights_corpus(corpus_, w_);
    }
    
    temp = n_dk.row(doc_id).sum();
    doc_each_len_weighted_multi.push_back(temp);
    total_words_weighted += temp;
  }
  n_s1_kv_multi.setFromTriplets(trip_s1_multi.begin(), trip_s1_multi.end());
}

void keyATMmultibase::initialize_specific()
{
  alpha = Rcpp::as<Eigen::VectorXd>(priors_list["alpha"]);
  
  estimate_alpha = options_list["estimate_alpha"];
  if (estimate_alpha == 0) {
    store_alpha = 0;
  } else {
    store_alpha = 1;
  }

  create_sufficient_stats();
}


void keyATMmultibase::resume_initialize_specific()
{
  estimate_alpha = options_list["estimate_alpha"];
  if (estimate_alpha == 0) {
    nv_alpha = priors_list["alpha"];
    alpha = Rcpp::as<Eigen::VectorXd>(nv_alpha);
    store_alpha = 0;
  } else {
    List alpha_iter = stored_values["alpha_iter"];
    NumericVector alpha_rvec = alpha_iter[alpha_iter.size() - 1];  // last estimated alpha
    alpha = Rcpp::as<Eigen::VectorXd>(alpha_rvec);
    store_alpha = 1;
  }

  create_sufficient_stats();
}


void keyATMmultibase::iteration_single(int it)
{ // Single iteration
  int doc_id_;
  int doc_length;
  int w_, z_, s_, c_;
  int new_z, new_s, new_c;
  int w_position;

  doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle

  for (int ii = 0; ii < num_doc; ++ii) {
    doc_id_ = doc_indexes[ii];
    doc_s = S[doc_id_], doc_z = Z[doc_id_], doc_w = W[doc_id_], doc_c = C[doc_id_];
    doc_length = doc_each_len[doc_id_];

    token_indexes = sampler::shuffled_indexes(doc_length); //shuffle

    // Iterate each word in the document
    for (int jj = 0; jj < doc_length; ++jj) {
      w_position = token_indexes[jj];
      s_ = doc_s[w_position], z_ = doc_z[w_position], w_ = doc_w[w_position], c_ = doc_c[w_position];

      new_z = keyATMmultibase::sample_z(alpha, z_, s_, w_, c_, doc_id_);
      doc_z[w_position] = new_z;

      if (keywords[new_z].find(w_) == keywords[new_z].end())
        continue;

      z_ = doc_z[w_position]; // use updated z
      new_s = keyATMmultibase::sample_s(z_, s_, w_, c_, doc_id_);
      doc_s[w_position] = new_s;

      s_ = doc_s[w_position]; // use updated s
      new_c = keyATMmultibase::sample_c(z_, s_, w_, c_, doc_id_);
      doc_c[w_position] = new_c;
    }

    Z[doc_id_] = doc_z;
    S[doc_id_] = doc_s;
    C[doc_id_] = doc_c;
  }

  sample_parameters(it);
}

void keyATMmultibase::sample_parameters(int it)
{
  if (estimate_alpha)
    sample_alpha();

  if (store_alpha) {
    int r_index = it + 1;
    if (r_index % thinning == 0 || r_index == 1 || r_index == iter) {
      NumericVector alpha_rvec = alpha_reformat(alpha, num_topics);
      List alpha_iter = stored_values["alpha_iter"];
      alpha_iter.push_back(alpha_rvec);
      stored_values["alpha_iter"] = alpha_iter;
    }
  }
}

int keyATMmultibase::sample_z(VectorXd &alpha, int z, int s, int w, int c, int doc_id)
{
  int new_z;
  double numerator, denominator;
  double sum;

  int corpus_ = corpus_id[doc_id];
  Ref<VectorXd> n_s0_c0_k_a = n_s0_c0_k_all.col(corpus_);

  // remove data
  if (s == 0) {
    n_s0_c0_k_a(z) -= vocab_weights_corpus(corpus_, w);
    n_s0_k(z) -= vocab_weights_corpus(corpus_, w);
    n_s0_kv(z, w) -= vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_s0_c0_kv_all[corpus_](z, w) -= vocab_weights_corpus(corpus_, w);
      n_c0_k(z) -= vocab_weights_corpus(corpus_, w);
      n_c0_k_all(z, corpus_) -= vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) -= vocab_weights_corpus(corpus_, w);
      n_s0_c1_k(z) -= vocab_weights_corpus(corpus_, w);
      n_s0_c1_kv(z, w) -= vocab_weights_corpus(corpus_, w);
      n_c1_k_all(z, corpus_) -= vocab_weights_corpus(corpus_, w);
    }
  } else if (s == 1) {
    n_s1_kv_multi.coeffRef(z, w) -= vocab_weights_corpus(corpus_, w);
    n_s1_k(z) -= vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_c0_k(z) -= vocab_weights_corpus(corpus_, w);
      n_c0_k_all(z, corpus_) -= vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) -= vocab_weights_corpus(corpus_, w);
      n_c1_k_all(z, corpus_) -= vocab_weights_corpus(corpus_, w);
    }
  } else {
    Rcerr << "Error at sample_z, remove" << std::endl;
  }

  n_dk(doc_id, z) -= vocab_weights_corpus(corpus_, w);

  new_z = -1; // initialize
  if (s == 0 && c == 0) {
    for (int k = 0; k < num_topics; ++k) {

      numerator = (beta + n_s0_c0_kv_all[corpus_](k, w)) *
        (n_s0_k(k) + prior_gamma(k, 1)) *
        (n_c0_k(k) + prior_omega(k, 1)) *
        (n_dk(doc_id, k) + alpha(k));

      denominator = (Vbeta + n_s0_c0_k_a(k)) *
        (n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1)) *
        (n_c1_k(k) + prior_omega(k, 0) + n_c0_k(k) + prior_omega(k, 1));

      z_prob_vec(k) = numerator / denominator;
    }

    sum = z_prob_vec.sum(); // normalize
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample

  } else if (s == 0 && c == 1) {
    for (int k = 0; k < num_topics; ++k) {

      numerator = (beta_c + n_s0_c1_kv(k, w)) *
        (n_s0_k(k) + prior_gamma(k, 1)) *
        (n_c1_k(k) + prior_omega(k, 0)) *
        (n_dk(doc_id, k) + alpha(k));

      denominator = (Vbeta_c + n_s0_c1_k(k)) *
        (n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1)) *
        (n_c1_k(k) + prior_omega(k, 0) + n_c0_k(k) + prior_omega(k, 1));

      z_prob_vec(k) = numerator / denominator;
    }

    sum = z_prob_vec.sum(); // normalize
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample
    
  } else if (s == 1 && c == 0) {
    for (int k = 0; k < num_topics; ++k) {
      if (keywords[k].find(w) == keywords[k].end()) {
        z_prob_vec(k) = 0.0;
        continue;
      } else {
        numerator = (beta_s + n_s1_kv_multi.coeffRef(k, w)) *
          (n_s1_k(k) + prior_gamma(k, 0)) *
          (n_c0_k(k) + prior_omega(k, 1)) *
          (n_dk(doc_id, k) + alpha(k));

        denominator = (Lbeta_sk(k) + n_s1_k(k)) *
          (n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1)) *
          (n_c1_k(k) + prior_omega(k, 0) + n_c0_k(k) + prior_omega(k, 1));

        z_prob_vec(k) = numerator / denominator;
      }
    }

    sum = z_prob_vec.sum(); // normalize
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample

  } else {
    for (int k = 0; k < num_topics; ++k) {
      if (keywords[k].find(w) == keywords[k].end()) {
        z_prob_vec(k) = 0.0;
        continue;
      } else {
        numerator = (beta_s + n_s1_kv_multi.coeffRef(k, w)) *
          (n_s1_k(k) + prior_gamma(k, 0)) *
          (n_c1_k(k) + prior_omega(k, 0)) *
          (n_dk(doc_id, k) + alpha(k));
        denominator = (Lbeta_sk(k) + n_s1_k(k) ) *
          (n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1)) *
          (n_c1_k(k) + prior_omega(k, 0) + n_c0_k(k) + prior_omega(k, 1));

        z_prob_vec(k) = numerator / denominator;
      }
    }

    sum = z_prob_vec.sum();
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample

  }

  // add back data counts
  if (s == 0) {
    n_s0_c0_k_a(new_z) += vocab_weights_corpus(corpus_, w);
    n_s0_k(new_z) += vocab_weights_corpus(corpus_, w);
    n_s0_kv(new_z, w) += vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_s0_c0_kv_all[corpus_](new_z, w) += vocab_weights_corpus(corpus_, w);
      n_c0_k(new_z) += vocab_weights_corpus(corpus_, w);
      n_c0_k_all(new_z, corpus_) += vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(new_z) += vocab_weights_corpus(corpus_, w);
      n_s0_c1_k(new_z) += vocab_weights_corpus(corpus_, w);
      n_s0_c1_kv(new_z, w) += vocab_weights_corpus(corpus_, w);
      n_c1_k_all(new_z, corpus_) += vocab_weights_corpus(corpus_, w);
    }
  } else if (s == 1) {
    n_s1_kv_multi.coeffRef(new_z, w) += vocab_weights_corpus(corpus_, w);
    n_s1_k(new_z) += vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_c0_k(new_z) += vocab_weights_corpus(corpus_, w);
      n_c0_k_all(new_z, corpus_) += vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(new_z) += vocab_weights_corpus(corpus_, w);
      n_c1_k_all(new_z, corpus_) += vocab_weights_corpus(corpus_, w);
    }
  } else {
    Rcerr << "Error at sample_z, add" << std::endl;
  }
  n_dk(doc_id, new_z) += vocab_weights_corpus(corpus_, w);

  return new_z;
}

int keyATMmultibase::sample_s(int z, int s, int w, int c, int doc_id)
{
  int new_s;
  double numerator, denominator;
  double s0_prob;
  double s1_prob;
  double sum;

  // Get corpus id and choose containers
  int corpus_ = corpus_id[doc_id];
  Ref<VectorXd> n_s0_c0_k_a = n_s0_c0_k_all.col(corpus_);

  // remove data
  if (s == 0) {
    n_s0_c0_k_a(z) -= vocab_weights_corpus(corpus_, w);
    n_s0_k(z) -= vocab_weights_corpus(corpus_, w);
    n_s0_kv(z, w) -= vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_s0_c0_kv_all[corpus_](z, w) -= vocab_weights_corpus(corpus_, w);
      n_c0_k(z) -= vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) -= vocab_weights_corpus(corpus_, w);
      n_s0_c1_k(z) -= vocab_weights_corpus(corpus_, w);
      n_s0_c1_kv(z, w) -= vocab_weights_corpus(corpus_, w);
    }
  } else if (s == 1) {
    n_s1_kv_multi.coeffRef(z, w) -= vocab_weights_corpus(corpus_, w);
    n_s1_k(z) -= vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_c0_k(z) -= vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) -= vocab_weights_corpus(corpus_, w);
    }
  } else {
    Rcerr << "Error at sample_z, remove" << std::endl;
  }

  // newprob_s1()
  numerator = (beta_s + n_s1_kv_multi.coeffRef(z, w)) *
    ( n_s1_k(z) + prior_gamma(z, 0) );
  denominator = (Lbeta_sk(z) + n_s1_k(z) );
  s1_prob = numerator / denominator;

  // newprob_s0()
  if (c == 0) {
    numerator = (beta + n_s0_c0_kv_all[corpus_](z, w)) *
        (n_s0_k(z) + prior_gamma(z, 1));

    denominator = (Vbeta + n_s0_c0_k_a(z) );
  } else {
    numerator = (beta_c + n_s0_c1_kv(z, w)) *
        (n_s0_k(z) + prior_gamma(z, 1));

    denominator = (Vbeta_c + n_s0_c1_k(z) );
  }
  s0_prob = numerator / denominator;
  
  // Normalize
  sum = s0_prob + s1_prob;
  s1_prob = s1_prob / sum;

  new_s = R::runif(0,1) <= s1_prob;  //new_s = Bern(s0_prob, s1_prob);

  // add back data counts
  if (new_s == 0) {
    n_s0_c0_k_a(z) += vocab_weights_corpus(corpus_, w);
    n_s0_k(z) += vocab_weights_corpus(corpus_, w);
    n_s0_kv(z, w) += vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_s0_c0_kv_all[corpus_](z, w) += vocab_weights_corpus(corpus_, w);
      n_c0_k(z) += vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) += vocab_weights_corpus(corpus_, w);
      n_s0_c1_k(z) += vocab_weights_corpus(corpus_, w);
      n_s0_c1_kv(z, w) += vocab_weights_corpus(corpus_, w);
    }
  } else if (new_s == 1) {
    n_s1_kv_multi.coeffRef(z, w) += vocab_weights_corpus(corpus_, w);
    n_s1_k(z) += vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_c0_k(z) += vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) += vocab_weights_corpus(corpus_, w);
    }
  }

  return new_s;
}

int keyATMmultibase::sample_c(int z, int s, int w, int c, int doc_id)
{
  int new_c;
  double numerator, denominator;
  double c0_prob;
  double c1_prob;
  double sum;

  // Get corpus id and choose containers
  int corpus_ = corpus_id[doc_id];
  Ref<VectorXd> n_s0_c0_k_a = n_s0_c0_k_all.col(corpus_);

  // remove data
  if (s == 0) {
    n_s0_c0_k_a(z) -= vocab_weights_corpus(corpus_, w);
    n_s0_k(z) -= vocab_weights_corpus(corpus_, w);
    n_s0_kv(z, w) -= vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_s0_c0_kv_all[corpus_](z, w) -= vocab_weights_corpus(corpus_, w);
      n_c0_k(z) -= vocab_weights_corpus(corpus_, w);
      n_c0_k_all(z, corpus_) -= vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) -= vocab_weights_corpus(corpus_, w);
      n_s0_c1_k(z) -= vocab_weights_corpus(corpus_, w);
      n_s0_c1_kv(z, w) -= vocab_weights_corpus(corpus_, w);
      n_c1_k_all(z, corpus_) -= vocab_weights_corpus(corpus_, w);
    }
  } else if (s == 1) {
    n_s1_kv_multi.coeffRef(z, w) -= vocab_weights_corpus(corpus_, w);
    n_s1_k(z) -= vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_c0_k(z) -= vocab_weights_corpus(corpus_, w);
      n_c0_k_all(z, corpus_) -= vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) -= vocab_weights_corpus(corpus_, w);
      n_c1_k_all(z, corpus_) -= vocab_weights_corpus(corpus_, w);
    }
  } else {
    Rcerr << "Error at sample_c, remove" << std::endl;
  }

  // newprob_c1()
  if (s == 0) {
    numerator = (beta_c + n_s0_c1_kv(z, w)) *
      ( n_c1_k(z) + prior_omega(z, 0) );
    denominator = (Vbeta_c + n_s0_c1_k(z) );
  } else {
    numerator = (beta_s + n_s1_kv_multi.coeffRef(z, w)) *
      ( n_c1_k(z) + prior_omega(z, 0) );
    denominator = (Lbeta_sk(z) + n_s1_k(z) );
  }
  c1_prob = numerator / denominator;

  // newprob_c0()
  if (s == 0) {
    numerator = (beta + n_s0_c0_kv_all[corpus_](z, w)) *
        (n_c0_k(z) + prior_omega(z, 1));

    denominator = (Vbeta + n_s0_c0_k_a(z) );
  } else {
    numerator = (beta_c + n_s0_c1_kv(z, w)) *
        (n_c0_k(z) + prior_omega(z, 1));

    denominator = (Vbeta_c + n_s0_c1_k(z) );
  }
  c0_prob = numerator / denominator;
  
  // Normalize
  sum = c0_prob + c1_prob;
  c1_prob = c1_prob / sum;

  new_c = R::runif(0,1) <= c1_prob;

  // add back data counts
  if (s == 0) {
    n_s0_c0_k_a(z) += vocab_weights_corpus(corpus_, w);
    n_s0_k(z) += vocab_weights_corpus(corpus_, w);
    n_s0_kv(z, w) += vocab_weights_corpus(corpus_, w);
    if (new_c == 0) {
      n_s0_c0_kv_all[corpus_](z, w) += vocab_weights_corpus(corpus_, w);
      n_c0_k(z) += vocab_weights_corpus(corpus_, w);
      n_c0_k_all(z, corpus_) += vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) += vocab_weights_corpus(corpus_, w);
      n_s0_c1_k(z) += vocab_weights_corpus(corpus_, w);
      n_s0_c1_kv(z, w) += vocab_weights_corpus(corpus_, w);
      n_c1_k_all(z, corpus_) += vocab_weights_corpus(corpus_, w);
    }
  } else if (s == 1) {
    n_s1_kv_multi.coeffRef(z, w) += vocab_weights_corpus(corpus_, w);
    n_s1_k(z) += vocab_weights_corpus(corpus_, w);
    if (new_c == 0) {
      n_c0_k(z) += vocab_weights_corpus(corpus_, w);
      n_c0_k_all(z, corpus_) += vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) += vocab_weights_corpus(corpus_, w);
      n_c1_k_all(z, corpus_) += vocab_weights_corpus(corpus_, w);
    }
  }

  return new_c;
}

void keyATMmultibase::sample_alpha()
{
  double start, end, previous_p, new_p, newlikelihood, slice_;
  keep_current_param = alpha;
  topic_ids = sampler::shuffled_indexes(num_topics);
  newalphallk = 0.0;
  int k;

  for (int i = 0; i < num_topics; ++i) {
    k = topic_ids[i];
    store_loglik = alpha_loglik(k);
    start = min_v ; // shrinked with shrinkp()
    end = max_v;  // shrinked with shrinkp()

    previous_p = alpha(k) / (1.0 + alpha(k)); // shrinkp
    slice_ = store_loglik - 2.0 * log(1.0 - previous_p)
            + log(unif_rand()); // <-- using R random uniform

    for (int shrink_time = 0; shrink_time < max_shrink_time; ++shrink_time) {
      new_p = sampler::slice_uniform(start, end); // <-- using R function above
      alpha(k) = new_p / (1.0 - new_p); // expandp

      newalphallk = alpha_loglik(k);
      newlikelihood = newalphallk - 2.0 * log(1.0 - new_p);

      if (slice_ < newlikelihood) {
        break;
      } else if (previous_p < new_p) {
        end = new_p;
      } else if (new_p < previous_p) {
        start = new_p;
      } else {
        Rcpp::stop("Something goes wrong in sample_alpha().");
        alpha(k) = keep_current_param(k);
        break;
      }
    }
  }
}


double keyATMmultibase::alpha_loglik(int k)
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

  for (int d = 0; d < num_doc; ++d) {
    loglik += fixed_part;

    // second term numerator
    loglik += mylgamma(ndk_a(d,k));

    // second term denominator
    loglik -= mylgamma(doc_each_len_weighted_multi[d] + alpha_sum_val);
  }

  return loglik;
}


double keyATMmultibase::loglik_total()
{
  double loglik = 0.0;
  double fixed_part = 0.0;

  for (int k = 0; k < num_topics; ++k) {
    for (int v = 0; v < num_vocab; ++v) { // word
      loglik += mylgamma(beta_c + n_s0_c1_kv(k, v)) - mylgamma(beta_c);
      for (int j = 0; j < num_corpus; ++j) {
        loglik += mylgamma(beta + n_s0_c0_kv_all[j](k, v)) - mylgamma(beta);
      }
    }

    loglik += mylgamma(beta_c * (double)num_vocab) - mylgamma(beta_c * (double)num_vocab + n_s0_c1_k(k));
    for (int j = 0; j < num_corpus; ++j) {
      loglik += mylgamma(beta * (double)num_vocab + n_s0_c0_k_all(k, j)) - mylgamma(beta * (double)num_vocab + n_s0_c0_k_all(k, j));
    }

    if (k < keyword_k) {
      for (SparseMatrix<double,RowMajor>::InnerIterator it(n_s1_kv_multi, k); it; ++it) {
        loglik += mylgamma(beta_s + it.value()) - mylgamma(beta_s);
      }
      loglik += mylgamma(beta_s * (double)keywords_num[k]) - mylgamma(beta_s * (double)keywords_num[k] + n_s1_k(k));
    
      loglik += mylgamma(prior_gamma(k, 0) + prior_gamma(k, 1)) - mylgamma(prior_gamma(k, 0)) - mylgamma(prior_gamma(k, 1));

      loglik += mylgamma(n_s0_k(k) + prior_gamma(k, 1))
                - mylgamma(n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1))
                + mylgamma(n_s1_k(k) + prior_gamma(k, 0));
    }

    for (int j = 0; j < num_corpus; ++j) {
      loglik += mylgamma(prior_omega(k, 0) + prior_omega(k, 1)) - mylgamma(prior_omega(k, 0)) - mylgamma(prior_omega(k, 1));
      
      loglik += mylgamma(n_c0_k_all(k, j) + prior_omega(k, 1))
                - mylgamma(n_c1_k_all(k, j) + prior_omega(k, 0) + n_c0_k_all(k, j) + prior_omega(k, 1))
                + mylgamma(n_c1_k_all(k, j) + prior_omega(k, 0));
    }

    if (k < keyword_k) {
      loglik += gammapdfln(alpha(k), eta_1, eta_2);
    } else {
      loglik += gammapdfln(alpha(k), eta_1_regular, eta_2_regular);
    }
  }

  fixed_part = alpha.sum();
  for (int d = 0; d < num_doc; ++d) {
    loglik += mylgamma(fixed_part) - mylgamma(doc_each_len_weighted_multi[d] + fixed_part);

    for (int k = 0; k < num_topics; ++k) {
      loglik += mylgamma(n_dk(d,k) + alpha(k)) - mylgamma(alpha(k));
    }
  }

  return loglik;
}