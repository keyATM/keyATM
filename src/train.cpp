#include <RcppEigen.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;
using namespace Rcpp;

double logsumexp(double &x, double &y, bool flg){
  if (flg) return y; // init mode
  if (x == y) return x + 0.69314718055; // log(2)
  double vmin = std::min(x, y);
  double vmax = std::max(x, y);
  if (vmax > vmin + 50){
    return vmax;
  } else {
    return vmax + std::log(std::exp(vmin - vmax) + 1.0);
  }
}

double logsumexp_Eigen(VectorXd &vec){
  double sum = 0.0;
  for(int i = 0; i < vec.size(); i++){
    sum = logsumexp(sum, vec[i], (i == 0));
  }
  return sum;
}

double gammapdfln(double x, double a, double b){
  return a * log(b) - lgamma(a) + (a - 1.0) * log(x) - b * x;
}

int rcat(Eigen::VectorXd &prob){ // was called 'multi1'
  double u = R::runif(0, 1);
  double temp = 0.0;
  int index = 0;
  for (int ii = 0; ii < prob.size(); ii++){
    temp += prob(ii);
    if (u < temp){
      index = ii;
      break;
    }
  }
  return index;
}

double loglikelihood1(SparseMatrix<int, RowMajor>& n_x0_kv,
                      SparseMatrix<int, RowMajor>& n_x1_kv,
                      VectorXi& n_x0_k, VectorXi& n_x1_k,
                      MatrixXd& n_dk, VectorXd& alpha,
                      double beta, double beta_s,
                      double gamma_1, double gamma_2,
                      double lambda_1, double lambda_2,
                      int num_topics, int k_seeded, int num_vocab, int num_doc,
                      std::vector< std::unordered_map<int, double> > & phi_s) {
  double loglik = 0.0;
  for (int k = 0; k < num_topics; k++){
    for (int v = 0; v < num_vocab; v++){ // word
      loglik += lgamma(beta + (double)n_x0_kv.coeffRef(k, v) ) - lgamma(beta);
      loglik += lgamma(beta_s + (double)n_x1_kv.coeffRef(k, v) ) - lgamma(beta_s);
    }
    // word normalization
    loglik += lgamma( beta * (double)num_vocab ) - lgamma(beta * (double)num_vocab + (double)n_x0_kv.row(k).sum() );
    loglik += lgamma( beta_s * (double)num_vocab ) - lgamma(beta_s * (double)num_vocab + (double)n_x1_kv.row(k).sum() );
    // x
    loglik += lgamma( (double)n_x0_k(k) + gamma_2 ) - lgamma((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2)
      + lgamma( (double)n_x1_k(k) + gamma_1 ) ;

		// std::cout << (double)n_x0_k(k) << " / " << (double)n_x1_k(k) << std::endl; // debug

    // x normalization
    loglik += lgamma(gamma_1 + gamma_2) - lgamma(gamma_1) - lgamma(gamma_2);
  }
  // z
  for (int d = 0; d < num_doc; d++){
    loglik += lgamma( alpha.sum() ) - lgamma( n_dk.row(d).sum() + alpha.sum() );
    for (int k = 0; k < num_topics; k++){
      loglik += lgamma( n_dk.coeffRef(d,k) + alpha(k) ) - lgamma( alpha(k) );
    }
  }
  return loglik;
}

// takes the sufficient statistics and parameters of the model and
// returns a new value for z
int sample_z(SparseMatrix<int, RowMajor>& n_x0_kv,
             SparseMatrix<int, RowMajor>& n_x1_kv,
             VectorXi& n_x0_k, VectorXi& n_x1_k,
             MatrixXd& n_dk, VectorXd& alpha,
             std::vector<int>& seed_num,
             int x, int z, int w, int doc_id, int num_vocab,
             int num_topics, int k_seeded,
             double gamma_1, double gamma_2, double beta, double beta_s,
             std::vector< std::unordered_map<int, double> > & phi_s){
  // remove data
  if (x == 0){
    n_x0_kv.coeffRef(z, w) -= 1;
    n_x0_k(z) -= 1;
  } else {
    n_x1_kv.coeffRef(z, w) -= 1;
    n_x1_k(z) -= 1;
  }
  n_dk.coeffRef(doc_id, z) -= 1;

  VectorXd z_prob_vec = VectorXd::Zero(num_topics);
  int new_z = -1; // debug
  double numerator, denominator;
  if (x == 0){
    for (int k = 0; k < num_topics; ++k){
      numerator = log(beta + (double)n_x0_kv.coeffRef(k, w)) +
        log((double)n_x0_k(k) + gamma_2) +
        log((double)n_dk.coeffRef(doc_id, k) + alpha(k));
      denominator = log((double)num_vocab * beta + (double)n_x0_kv.row(k).sum()) +
        log((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);
      z_prob_vec(k) = numerator - denominator;
    }
    double sum = logsumexp_Eigen(z_prob_vec); // normalize
    for (int k = 0; k < num_topics; k++)
      z_prob_vec(k) = exp(z_prob_vec(k) - sum);
    new_z = rcat(z_prob_vec); // take a sample

  } else {
    std::vector<int> make_zero_later;
    for (int k = 0; k < num_topics; ++k){ 
      if (phi_s[k].find(w) == phi_s[k].end()){
        z_prob_vec(k) = 1.0;
        make_zero_later.push_back(k);
        continue;
      } else{ // w not one of the seeds
        numerator = log(beta_s + (double)n_x1_kv.coeffRef(k, w)) +
          log( ((double)n_x1_k(k) + gamma_1) ) +
          log( ((double)n_dk.coeffRef(doc_id, k) + alpha(k)) );
      }
      denominator = log((double)seed_num[k] * beta_s + (double)n_x1_kv.row(k).sum() ) +
        log((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);
      z_prob_vec(k) = numerator - denominator;
    }
    double sum = logsumexp_Eigen(z_prob_vec);

    for (int k = 0; k < num_topics; k++)
      z_prob_vec(k) = exp(z_prob_vec(k) - sum);

    for (int k = 0; k < make_zero_later.size(); k++){ // zero out elements
			int change_k = make_zero_later[k];
      z_prob_vec(change_k) = 0.0; // make it 0 explicitly
		}

    z_prob_vec = z_prob_vec / z_prob_vec.sum(); // and renormalize
    new_z = rcat(z_prob_vec); // take a sample

  }

  // add back data counts
  if (x == 0){
    n_x0_kv.coeffRef(new_z, w) += 1;
    n_x0_k(new_z) += 1;
  } else {
    n_x1_kv.coeffRef(new_z, w) += 1;
    n_x1_k(new_z) += 1;
  }
  n_dk.coeffRef(doc_id, new_z) += 1;

  return new_z;
}

// takes the sufficient statistics and parameters of the model and
// returns a new value for x
int sample_x(SparseMatrix<int, RowMajor>& n_x0_kv,
             SparseMatrix<int, RowMajor>& n_x1_kv,
             VectorXi& n_x0_k, VectorXi& n_x1_k,
             MatrixXd& n_dk, VectorXd& alpha,
             std::vector<int>& seed_num,
             int x, int z, int w, int doc_id, int num_vocab,
             int num_topics, int k_seeded,
             double gamma_1, double gamma_2, double beta, double beta_s,
             std::vector< std::unordered_map<int, double> > & phi_s){
  // remove data
  if (x == 0){
    n_x0_kv.coeffRef(z, w) -= 1;
    n_x0_k(z) -= 1;
  } else {
    n_x1_kv.coeffRef(z, w) -= 1;
    n_x1_k(z) -= 1;
  }
  n_dk.coeffRef(doc_id, z) -= 1; // not necessary to remove 

  // newprob_x1()
  double x1_logprob;
  int k = z;
  double numerator;
  double denominator;
  if ( phi_s[k].find(w) == phi_s[k].end() ){
       x1_logprob = -1.0;
  } else {
    numerator = log(beta_s + (double)n_x1_kv.coeffRef(k, w)) +
      log( ((double)n_x1_k(k) + gamma_1) );
    denominator = log((double)seed_num[k] * beta_s + (double)n_x1_kv.row(k).sum() ) +
      log((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);
    x1_logprob = numerator - denominator;
  }

  int new_x;
  if(x1_logprob == -1.0){
    // if probability of x_di = 1 case is 0, it should be x=0 (regular topic)
    new_x = 0;
  } else {
    // newprob_x0()
    int k = z;
    numerator = log(beta + (double)n_x0_kv.coeffRef(k, w)) +
      log((double)n_x0_k(k) + gamma_2);

    denominator = log((double)num_vocab * beta + (double)n_x0_kv.row(k).sum() ) +
      log((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);
    double x0_logprob = numerator - denominator;

    // Normalize
    double sum = 0.0;
    double prob_array[] = {x0_logprob, x1_logprob};
    for (int i = 0; i < 2; i++)
      sum = logsumexp(sum, prob_array[i], (i == 0));

    double x0_prob = exp( x0_logprob - sum);
    double x1_prob = exp( x1_logprob - sum);
    new_x = R::runif(0,1) <= x1_prob;  //new_x = Bern(x0_prob, x1_prob);
  }
  // add back data counts
  if (new_x == 0){
    n_x0_kv.coeffRef(z, w) += 1;
    n_x0_k(z) += 1;
  } else {
    n_x1_kv.coeffRef(z, w) += 1;
    n_x1_k(z) += 1;
  }
  n_dk.coeffRef(doc_id, z) += 1;

  return new_x;
}

// Wrapper around R's uniform random number generator for std shuffles
inline int rand_wrapper(const int n) { return floor(unif_rand() * n); }

// Turns K into random permutation of 0 to K-1 using R's RNG
std::vector<int> shuffled_indexes(int m) {
  std::vector<int> v(m);
  std::iota(v.begin(), v.end(), 0);
  std::random_shuffle(v.begin(), v.end(), rand_wrapper);
  return v;
}

// a version of randgen::uniform(lower, upper) that uses R's RNG
double slice_uniform(double lower, double upper){
  return lower + (upper - lower) * unif_rand();
}

// This used to differentiate between input_alpha and alpha
// for reasons I could not understand.  TODO: check this is still correct?
double alpha_loglik(VectorXd &alpha,  MatrixXd& n_dk,
                    int num_topics, int num_doc){
  double loglik = 0.0;
  double fixed_part = 0.0;
  VectorXd ndk_ak;

  fixed_part += lgamma(alpha.sum()); // first term numerator
  for(int k = 0; k < num_topics; k++){
    fixed_part -= lgamma(alpha(k)); // first term denominator
    // Add prior
    loglik += gammapdfln(alpha(k), 1.0, 2.0);
  }
  for(int d = 0; d < num_doc; d++){
    loglik += fixed_part;
    ndk_ak = n_dk.row(d) + alpha;
    // second term numerator
    for(int k = 0; k < num_topics; k++){
      loglik += lgamma(ndk_ak(k));
    }
    // second term denominator
    loglik -= lgamma(ndk_ak.sum());
  }
  return loglik;
}

// updates alpha in place, currently just hands back alpha (by reference)
VectorXd& slice_sample_alpha(VectorXd& alpha, MatrixXd& n_dk,
                             int num_topics, int num_doc,
                             double min_v = 1e-9, double max_v = 100.0,
                             int max_shrink_time = 3000){

  double start, end, previous_p, new_p, newlikelihood, slice_;
  VectorXd keep_current_param = alpha;
  std::vector<int> topic_ids = shuffled_indexes(num_topics);
  for(int i = 0; i < num_topics; i++){
    int k = topic_ids[i];
    start = min_v / (1.0 + min_v); // shrinkp
    end = 1.0;
    // end = shrinkp(max_v);
    previous_p = alpha(k) / (1.0 + alpha(k)); // shrinkp
    slice_ = alpha_loglik(alpha, n_dk, num_topics, num_doc)
              - 2.0 * log(1.0 - previous_p)
              + log(unif_rand()); // <-- using R random uniform

    for (int shrink_time = 0; shrink_time < max_shrink_time; shrink_time++){
      new_p = slice_uniform(start, end); // <-- using R function above
      alpha(k) = new_p / (1.0 - new_p); // expandp
      newlikelihood = alpha_loglik(alpha, n_dk, num_topics, num_doc)
                      - 2.0 * log(1.0 - new_p);
      if (slice_ < newlikelihood){
        break;
      } else if (previous_p < new_p){
        end = new_p;
      } else if (new_p < previous_p){
        start = new_p;
      } else {
        alpha(k) = keep_current_param(k);
      }
    }
  }
  return alpha;
}

//' Run the Gibbs sampler
//'
//' @param model A model, from \code{init} or a previous invocation of \code{train}
//' @param alpha_k A starting value for alpha (will be removed when alpha updates are back in)
//' @param iter Required number of iterations
//'
//' @export
// [[Rcpp::export]]
List topicdict_train(List model, double alpha_k, int iter = 0){

	// Data
  List W = model["W"], Z = model["Z"], X = model["X"];
  StringVector files = model["files"], vocab = model["vocab"];
  int k_free = model["extra_k"];
  List seeds = model["seeds"];
  int k_seeded = seeds.size();

  // alpha-related constants
  int num_topics = k_seeded + k_free;
	alpha_k /= (double)num_topics; // recommended alpha initialization in Griffiths and Steyvers (2004)
  VectorXd alpha = VectorXd::Constant(num_topics, alpha_k);


	// phi_s
  std::vector< std::unordered_map<int, double> > phi_s(num_topics);
  std::vector<int> seed_num(num_topics);
  for (int ii = 0; ii < k_seeded; ii++){
    IntegerVector wd_ids = seeds[ii];
    seed_num[ii] = wd_ids.size();
    std::unordered_map<int, double> phi_sk;
    for (int jj = 0; jj < wd_ids.size(); jj++)
      phi_sk[wd_ids(jj)] = 1.0 / wd_ids.size();
    phi_s[ii] = phi_sk;
  }
	for(int i=k_seeded; i<num_topics; i++){
		std::unordered_map<int, double> phi_sk{ {-1, -1.0} };
		seed_num[i] = 0;
		phi_s[i] = phi_sk;
	}


  // document-related constants
  int num_vocab = vocab.size(), num_doc = files.size();

  // distributional constants
  double gamma_1 = 1.0, gamma_2 = 1.0;
  double lambda_1 = 1.0, lambda_2 = 2.0;
  double beta = 0.01, beta_s = 0.1;

  // storage for sufficient statistics and their margins
  SparseMatrix<int, RowMajor> n_x0_kv = SparseMatrix<int, RowMajor> (num_topics, num_vocab);
  SparseMatrix<int, RowMajor> n_x1_kv = SparseMatrix<int, RowMajor> (num_topics, num_vocab);
  MatrixXd n_dk = MatrixXd::Zero(num_doc, num_topics);
  MatrixXd theta_dk = MatrixXd::Zero(num_doc, num_topics);
  VectorXi n_x0_k = VectorXi::Zero(num_topics);
  VectorXi n_x1_k = VectorXi::Zero(num_topics);

  for(int doc_id = 0; doc_id < num_doc; doc_id++){
    IntegerVector doc_x = X[doc_id], doc_z = Z[doc_id], doc_w = W[doc_id];
    for(int w_position = 0; w_position < doc_x.size(); w_position++){
      int x = doc_x[w_position], z = doc_z[w_position], w = doc_w[w_position];
      if (x == 0){
        n_x0_kv.coeffRef(z, w) += 1;
        n_x0_k(z) += 1;
      } else {
        n_x1_kv.coeffRef(z, w) += 1;
        n_x1_k(z) += 1;
      }
      n_dk.coeffRef(doc_id, z) += 1.0;
    }
  }
  int total_words = (int)n_dk.sum();


  // Randomized update sequence
  for (int it = 0; it < iter; it++){
    std::vector<int> doc_indexes = shuffled_indexes(num_doc); // shuffle
    for (int ii = 0; ii < num_doc; ii++){
      int doc_id = doc_indexes[ii];
      IntegerVector doc_x = X[doc_id], doc_z = Z[doc_id], doc_w = W[doc_id];

      std::vector<int> token_indexes = shuffled_indexes(doc_x.size()); //shuffle
      for (int jj = 0; jj < doc_x.size(); jj++){
        int w_position = token_indexes[jj];
        int x = doc_x[w_position], z = doc_z[w_position], w = doc_w[w_position];

        doc_z[w_position] = sample_z(n_x0_kv, n_x1_kv, n_x0_k, n_x1_k, n_dk,
                                     alpha, seed_num, x, z, w, doc_id,
                                     num_vocab, num_topics, k_seeded, gamma_1,
                                     gamma_2, beta, beta_s, phi_s);

				// Debug
				for(int s=0; s<num_topics; ++s){
					if((double)n_x0_k(s) < 0 | (double)n_x1_k(s) < 0)
						std::cout << (double)n_x0_k(s) << " / " << (double)n_x1_k(s) << std::endl;
				}

				z = doc_z[w_position]; // use updated z

        doc_x[w_position] = sample_x(n_x0_kv, n_x1_kv, n_x0_k, n_x1_k, n_dk,
                                    alpha, seed_num, x, z, w, doc_id,
                                    num_vocab, num_topics, k_seeded, gamma_1,
                                    gamma_2, beta, beta_s, phi_s);
      }
    }
    slice_sample_alpha(alpha, n_dk, num_topics, num_doc);

    double loglik = loglikelihood1(n_x0_kv, n_x1_kv, n_x0_k, n_x1_k, n_dk, alpha,
                                   beta, beta_s, gamma_1, gamma_2, lambda_1, lambda_2,
                                   num_topics, k_seeded, num_vocab, num_doc, phi_s);
    double perplexity = exp(-loglik / (double)total_words);
    Rcerr << "log likelihood: " << loglik <<
             " (perplexity: " << perplexity << ")" << std::endl;

    checkUserInterrupt();
  }
  return model;
}

// int bern(double &prob0, double &prob1){
//  double value = R::runif(0,1);
//  if(value <= prob1)
//    return 1;
//  return 0;
//}
//
// This function replaced by:   R::runif(0,1) <= prob1

/*
 double loglikelihood2(SparseMatrix<int, RowMajor>& n_x0_kv,
                       SparseMatrix<int, RowMajor>& n_x1_kv,
VectorXi& n_x0_k, VectorXi& n_x1_k,
MatrixXd& n_dk, VectorXd& alpha,
double beta, double beta_s,
double gamma_1, double gamma_2,
double lambda_1, double lambda_2,
int num_topics, int num_doc) {
double prod_k = 0.0;// calc_loglik_prod_k();
for (int k = 0; k < num_topics; k++){
// (a) Seed Topic Part
prod_k += lgamma((double)num_vocab * beta_s); // first term numerator
prod_k -= lgamma(beta_s) * (double)num_vocab; // first term denominator
prod_k -= lgamma( (double)num_vocab * beta_s + (double)n_x1_kv.row(k).sum()  ); // second term denominator
// Regular Topic Part
prod_k += lgamma((double)num_vocab * beta); //probably constant
prod_k -= lgamma(beta) * (double)num_vocab; // probably constant
prod_k -= lgamma( (double)num_vocab * beta + (double)n_x0_kv.row(k).sum()  ); // last part denominator
// (b) second part
prod_k += lgamma((double)n_x1_k(k) + gamma_1) + lgamma((double)n_x0_k(k) + gamma_2);
prod_k -=  lgamma( (double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2); // probably constant
// (c) prior for alpha
prod_k += gammapdfln(alpha(k), lambda_1, lambda_2);
}
int v = 0, prod_v = 0.0; // calc_loglik_prod_v();
for (size_t v_sizet = 0; v_sizet < num_vocab; v_sizet++){
v = static_cast<int>(v_sizet);
for (int k = 0; k < k_seeded; k++){ // <--------- NOTE: using k_seeded
// (a), seed topic part, second part numerator
if (phi_s[k].find(v) == phi_s[k].end()){ // not a seed word
prod_v += lgamma(beta_s);
} else {
prod_v += lgamma(beta_s + (double)n_x1_kv.coeffRef(k, v) );
}
// (a) regular part numerator
prod_v += lgamma(beta + (double)n_x0_kv.coeffRef(k, v) );
}
}
double prod_d = 0.0; // calc_loglik_prod_d();
for (int d = 0; d < num_doc; d++){ // (c)
prod_d += lgamma(alpha.sum());
for (int k = 0; k < k_seeded; k++){ // <--------- NOTE: using k_seeded
prod_d += lgamma(n_dk.coeffRef(d,k) + alpha(k)); // second numerator
prod_d -= lgamma(alpha(k)); // first denominator
}
prod_d -= lgamma(n_dk.row(d).sum() + alpha.sum()); // second denominator
}
// (b) first part
double others = (double)num_topics *
(lgamma(gamma_1 + gamma_2) - lgamma(gamma_1) - lgamma(gamma_2)); // constant?
double loglik = prod_k + prod_v + prod_d + others;
return loglik;
}
*/

