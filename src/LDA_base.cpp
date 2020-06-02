#include "LDA_base.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;


void LDAbase::read_data_common()
{
  // Read data
  W = model["W"]; Z = model["Z"];
  vocab = model["vocab"];
  regular_k = model["no_keyword_topics"];
  model_fit = model["model_fit"];

  keyword_k = 0;
  num_topics = regular_k;

  // document-related constants
  num_vocab = vocab.size();
  num_doc = W.size();
  // alpha -> specific function  
  
  // Options
  options_list = model["options"];
  use_weights = options_list["use_weights"];
  slice_A = options_list["slice_shape"];
  store_theta = options_list["store_theta"];
  thinning = options_list["thinning"];
  llk_per = options_list["llk_per"];
  verbose = options_list["verbose"];
  weights_type = Rcpp::as<std::string>(options_list["weights_type"]);

  // Priors
  priors_list = model["priors"];
  beta = priors_list["beta"];

  // Stored values
  stored_values = model["stored_values"];

  // Slice Sampling
  // this is used except cov models
  model_settings = model["model_settings"];
  min_v = model_settings["slice_min"];
  min_v = shrinkp(min_v);

  max_v = model_settings["slice_max"];
  max_v = shrinkp(max_v);
}


void LDAbase::initialize_common()
{
  // Parameters
  eta_1 = 1.0;
  eta_2 = 1.0;
  eta_1_regular = 2.0;
  eta_2_regular = 1.0;

  // Slice sampling initialization
  max_shrink_time = 200;

  // No labels in LDA
  use_labels = 0;
 
  //
  // Vocabulary weights
  //
  vocab_weights = VectorXd::Constant(num_vocab, 1.0);

  int z, w;
  int doc_len;
  IntegerVector doc_z, doc_w;

  
  // Construct vocab weights
  for (int doc_id = 0; doc_id < num_doc; doc_id++) {
    doc_w = W[doc_id];
    doc_len = doc_w.size();
    doc_each_len.push_back(doc_len);
  
    for (int w_position = 0; w_position < doc_len; w_position++) {
      w = doc_w[w_position];
      vocab_weights(w) += 1.0;
    }
  }
  total_words = (int)vocab_weights.sum();

  if (weights_type == "inv-freq" || weights_type == "inv-freq-normalized") {
    // Inverse frequency
    weights_invfreq(); 
  } else if (weights_type == "information-theory" || weights_type == "information-theory-normalized") {
    // Information theory 
    weights_inftheory();
  }
    
  // Normalize weights
  if (weights_type == "inv-freq-normalized" || 
      weights_type == "information-theory-normalized") {
    weights_normalize_total(); 
  } 

  // Do you want to use weights?
  if (use_weights == 0) {
    Rcpp::Rcout << "Not using weights!! Check `options$use_weights`." << std::endl;
    vocab_weights = VectorXd::Constant(num_vocab, 1.0);
  }


  //
  // Construct data matrices
  //
  n_kv = MatrixXd::Zero(num_topics, num_vocab);
  n_dk = MatrixXd::Zero(num_doc, num_topics);
  n_dk_noWeight = MatrixXd::Zero(num_doc, num_topics);
  n_k = VectorXd::Zero(num_topics);

  total_words_weighted = 0.0;
  double temp;
  for(int doc_id = 0; doc_id < num_doc; doc_id++){
    doc_z = Z[doc_id], doc_w = W[doc_id];
    doc_len = doc_each_len[doc_id];

    for(int w_position = 0; w_position < doc_len; w_position++){
      z = doc_z[w_position], w = doc_w[w_position];

      n_kv(z, w) += vocab_weights(w);
      n_k(z) += vocab_weights(w);
      n_dk(doc_id, z) += vocab_weights(w);
      n_dk_noWeight(doc_id, z) += 1.0;
    }

    temp = n_dk.row(doc_id).sum();
    doc_each_len_weighted.push_back(temp);
    total_words_weighted += temp;
  }
  

  // Use during the iteration
  z_prob_vec = VectorXd::Zero(num_topics);

}


void LDAbase::parameters_store(int r_index)
{
  if (store_theta)
    store_theta_iter(r_index);
}


int LDAbase::sample_z(VectorXd &alpha, int z, int s, int w, int doc_id)
{
  int new_z;
  double numerator, denominator;
  double sum;

  // remove data
  n_kv(z, w) -= vocab_weights(w);
  n_k(z) -= vocab_weights(w);
  n_dk(doc_id, z) -= vocab_weights(w);
  n_dk_noWeight(doc_id, z) -= 1;

  new_z = -1; // debug


  for (int k = 0; k < num_topics; ++k) {

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
  n_dk(doc_id, new_z) += vocab_weights(w);
  n_dk_noWeight(doc_id, new_z) += 1;

  return new_z;
}

