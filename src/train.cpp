#include <RcppEigen.h>
#include <chrono>
#include <iostream>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]


using namespace Eigen;
using namespace Rcpp;
using namespace std;

double time_z_rcat = 0.0;
double time_z_prepare_vec = 0.0;
double time_z_denom = 0.0;
std::chrono::high_resolution_clock::time_point  time_start, time_end;
std::chrono::high_resolution_clock::time_point  time_start_z, time_end_z;
std::chrono::high_resolution_clock::time_point  time_start_z1, time_end_z1;

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

NumericVector alpha_reformat(VectorXd& alpha, int& num_topics){
	NumericVector alpha_rvec(num_topics);

	for(int i=0; i<num_topics; ++i){
		alpha_rvec[i] = alpha(i);
	}

	return alpha_rvec;
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

int rcat_without_normalize(Eigen::VectorXd &prob, double &total){ // was called 'multi1'
  double u = R::runif(0, 1) * total;
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

double loglikelihood1(MatrixXi& n_x0_kv,
                      MatrixXi& n_x1_kv,
                      VectorXi& n_x0_k, VectorXi& n_x1_k,
                      MatrixXd& n_dk, VectorXd& alpha,
                      double beta, double beta_s,
                      double gamma_1, double gamma_2,
                      int num_topics, int k_seeded, int num_vocab, int num_doc,
                      std::vector< std::unordered_map<int, double> > & phi_s) {
  double loglik = 0.0;
  for (int k = 0; k < num_topics; k++){
    for (int v = 0; v < num_vocab; v++){ // word
      loglik += lgamma(beta + (double)n_x0_kv(k, v) ) - lgamma(beta);
      loglik += lgamma(beta_s + (double)n_x1_kv(k, v) ) - lgamma(beta_s);
    }
    // word normalization
    loglik += lgamma( beta * (double)num_vocab ) - lgamma(beta * (double)num_vocab + (double)n_x0_k(k) );
    loglik += lgamma( beta_s * (double)num_vocab ) - lgamma(beta_s * (double)num_vocab + (double)n_x1_k(k) );
    // x
    loglik += lgamma( (double)n_x0_k(k) + gamma_2 ) - lgamma((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2)
      + lgamma( (double)n_x1_k(k) + gamma_1 ) ;

		// Rcout << (double)n_x0_k(k) << " / " << (double)n_x1_k(k) << std::endl; // debug

    // x normalization
    loglik += lgamma(gamma_1 + gamma_2) - lgamma(gamma_1) - lgamma(gamma_2);
  }
  // z
  for (int d = 0; d < num_doc; d++){
    loglik += lgamma( alpha.sum() ) - lgamma( n_dk.row(d).sum() + alpha.sum() );
    for (int k = 0; k < num_topics; k++){
      loglik += lgamma( n_dk(d,k) + alpha(k) ) - lgamma( alpha(k) );
    }
  }
  return loglik;
}

// takes the sufficient statistics and parameters of the model and
// returns a new value for z
int sample_z(MatrixXi& n_x0_kv,
             MatrixXi& n_x1_kv,
             VectorXi& n_x0_k, VectorXi& n_x1_k,
             MatrixXd& n_dk, VectorXd& alpha,
             std::vector<int>& seed_num,
             int x, int z, int w, int doc_id, int num_vocab,
             int num_topics, int k_seeded,
             double gamma_1, double gamma_2, double beta, double beta_s,
             std::vector< std::unordered_map<int, double> > & phi_s){
  // remove data
  if (x == 0){
    n_x0_kv(z, w) -= 1;
    n_x0_k(z) -= 1;
  } else if (x==1) {
    n_x1_kv(z, w) -= 1;
    n_x1_k(z) -= 1;
  } else {
		Rcerr << "Error at sample_z, remove" << std::endl;
	}

  n_dk(doc_id, z) -= 1;

  VectorXd z_prob_vec = VectorXd::Zero(num_topics);
  int new_z = -1; // debug
  double numerator, denominator;
  if (x == 0){
		time_start_z = std::chrono::high_resolution_clock::now();
    for (int k = 0; k < num_topics; ++k){
      numerator = (beta + (double)n_x0_kv(k, w)) *
        ((double)n_x0_k(k) + gamma_2) *
        ((double)n_dk(doc_id, k) + alpha(k));

			time_start_z1 = std::chrono::high_resolution_clock::now();
      denominator = ((double)num_vocab * beta + (double)n_x0_k(k)) *
        ((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);
			time_end_z1 = std::chrono::high_resolution_clock::now();
			time_z_denom += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_z1-time_start_z1).count();

      z_prob_vec(k) = numerator / denominator;
    }
		time_end_z = std::chrono::high_resolution_clock::now();
		time_z_prepare_vec += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_z-time_start_z).count();

    double sum = z_prob_vec.sum(); // normalize
		time_start_z = std::chrono::high_resolution_clock::now();
    new_z = rcat_without_normalize(z_prob_vec, sum); // take a sample
		time_end_z = std::chrono::high_resolution_clock::now();
		time_z_rcat += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_z-time_start_z).count();

  } else {
		time_start_z = std::chrono::high_resolution_clock::now();
    for (int k = 0; k < num_topics; ++k){
      if (phi_s[k].find(w) == phi_s[k].end()){
        z_prob_vec(k) = 0.0;
				continue;
      } else{ // w not one of the seeds
        numerator = (beta_s + (double)n_x1_kv(k, w)) *
          ( ((double)n_x1_k(k) + gamma_1) ) *
          ( ((double)n_dk(doc_id, k) + alpha(k)) );
      }
			time_start_z1 = std::chrono::high_resolution_clock::now();
      denominator = ((double)seed_num[k] * beta_s + (double)n_x1_k(k) ) *
        ((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);
			time_end_z1 = std::chrono::high_resolution_clock::now();
			time_z_denom += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_z1-time_start_z1).count();

      z_prob_vec(k) = numerator / denominator;
    }
		time_end_z = std::chrono::high_resolution_clock::now();
		time_z_prepare_vec += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_z-time_start_z).count();


		double sum = z_prob_vec.sum();
		time_start_z = std::chrono::high_resolution_clock::now();
    new_z = rcat_without_normalize(z_prob_vec, sum); // take a sample
		time_end_z = std::chrono::high_resolution_clock::now();
		time_z_rcat += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_z-time_start_z).count();

  }

  // add back data counts
  if (x == 0){
    n_x0_kv(new_z, w) += 1;
    n_x0_k(new_z) += 1;
  } else if (x==1) {
    n_x1_kv(new_z, w) += 1;
    n_x1_k(new_z) += 1;
  } else {
		Rcerr << "Error at sample_z, add" << std::endl;
	}
  n_dk(doc_id, new_z) += 1;

  return new_z;
}


// takes the sufficient statistics and parameters of the model and
// returns a new value for x
int sample_x(MatrixXi& n_x0_kv,
             MatrixXi& n_x1_kv,
             VectorXi& n_x0_k, VectorXi& n_x1_k,
             MatrixXd& n_dk, VectorXd& alpha,
             std::vector<int>& seed_num,
             int x, int z, int w, int doc_id, int num_vocab,
             int num_topics, int k_seeded,
             double gamma_1, double gamma_2, double beta, double beta_s,
             std::vector< std::unordered_map<int, double> > & phi_s){
  // remove data
  if (x == 0){
    n_x0_kv(z, w) -= 1;
    n_x0_k(z) -= 1;
  } else {
    n_x1_kv(z, w) -= 1;
    n_x1_k(z) -= 1;
  }
  n_dk(doc_id, z) -= 1; // not necessary to remove

  // newprob_x1()
  double x1_prob;
  int k = z;
  double numerator;
  double denominator;
  if ( phi_s[k].find(w) == phi_s[k].end() ){
       x1_prob = -1.0;
  } else {
    numerator = (beta_s + (double)n_x1_kv(k, w)) *
      ( ((double)n_x1_k(k) + gamma_1) );
    denominator = ((double)seed_num[k] * beta_s + (double)n_x1_k(k) ) *
      ((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);
    x1_prob = numerator / denominator;
  }

  int new_x;
  if(x1_prob == -1.0){
    // if probability of x_di = 1 case is 0, it should be x=0 (regular topic)
    new_x = 0;
  } else {
    // newprob_x0()
    int k = z;
    numerator = (beta + (double)n_x0_kv(k, w)) *
      ((double)n_x0_k(k) + gamma_2);

    denominator = ((double)num_vocab * beta + (double)n_x0_k(k) ) *
      ((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);
    double x0_prob = numerator / denominator;

    // Normalize
    double sum = x0_prob + x1_prob;

    x1_prob = x1_prob / sum;
    new_x = R::runif(0,1) <= x1_prob;  //new_x = Bern(x0_prob, x1_prob);
  }
  // add back data counts
  if (new_x == 0){
    n_x0_kv(z, w) += 1;
    n_x0_k(z) += 1;
  } else {
    n_x1_kv(z, w) += 1;
    n_x1_k(z) += 1;
  }
  n_dk(doc_id, z) += 1;

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
double alpha_loglik(VectorXd &alpha, MatrixXd& n_dk,
                    int num_topics, int num_doc){
  double loglik = 0.0;
  double fixed_part = 0.0;
	double eta_1 = 0.5;
	double eta_2 = 5;
  VectorXd ndk_ak;


  fixed_part += lgamma(alpha.sum()); // first term numerator
  for(int k = 0; k < num_topics; k++){
    fixed_part -= lgamma(alpha(k)); // first term denominator
    // Add prior
    loglik += gammapdfln(alpha(k), eta_1, eta_2);

  }
  for(int d = 0; d < num_doc; d++){
    loglik += fixed_part;
    ndk_ak = n_dk.row(d) + alpha.transpose();
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
                             int max_shrink_time = 1000){

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
				Rcerr << "Something goes wrong in slice_sample_alpha()" << std::endl;
        alpha(k) = keep_current_param(k);
      }
    }
  }
  return alpha;
}

//' Run the Gibbs sampler
//'
//' @param model A model, from \code{init} or a previous invocation of \code{train}
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @export
// [[Rcpp::export]]
List topicdict_train(List model, int iter = 0, int output_per = 10){
	auto start = std::chrono::high_resolution_clock::now();

	// Data
  List W = model["W"], Z = model["Z"], X = model["X"];
  StringVector files = model["files"], vocab = model["vocab"];
  NumericVector nv_alpha = model["alpha"];
  double gamma_1 = model["gamma_1"], gamma_2 = model["gamma_2"];
  double beta = model["beta"], beta_s = model["beta_s"];
  int k_free = model["extra_k"];
  List seeds = model["seeds"];
  int k_seeded = seeds.size();
	List model_fit = model["model_fit"];

  // alpha
  int num_topics = k_seeded + k_free;
	VectorXd alpha = Rcpp::as<Eigen::VectorXd>(nv_alpha);

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

  // storage for sufficient statistics and their margins
  MatrixXi n_x0_kv = MatrixXi::Zero(num_topics, num_vocab);
  MatrixXi n_x1_kv = MatrixXi::Zero(num_topics, num_vocab);
  MatrixXd n_dk = MatrixXd::Zero(num_doc, num_topics);
  MatrixXd theta_dk = MatrixXd::Zero(num_doc, num_topics);
  VectorXi n_x0_k = VectorXi::Zero(num_topics);
  VectorXi n_x1_k = VectorXi::Zero(num_topics);


  for(int doc_id = 0; doc_id < num_doc; doc_id++){
    IntegerVector doc_x = X[doc_id], doc_z = Z[doc_id], doc_w = W[doc_id];
    for(int w_position = 0; w_position < doc_x.size(); w_position++){
      int x = doc_x[w_position], z = doc_z[w_position], w = doc_w[w_position];
      if (x == 0){
        n_x0_kv(z, w) += 1;
        n_x0_k(z) += 1;
      } else {
        n_x1_kv(z, w) += 1;
        n_x1_k(z) += 1;
      }
      n_dk(doc_id, z) += 1.0;
    }
  }
  int total_words = (int)n_dk.sum();

	double prepare_data = std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::high_resolution_clock::now() - start).count();
	double time_z, time_x, time_alpha, time_loglik, time_getvector, rpushback, shuffle_time;
			time_z = 0.0 ; time_x = 0.0; time_alpha=0.0; time_loglik=0.0; time_getvector=0.0; rpushback=0.0; shuffle_time=0.0;


  // Randomized update sequence
  for (int it = 0; it < iter; it++){
		time_start = std::chrono::high_resolution_clock::now();
    std::vector<int> doc_indexes = shuffled_indexes(num_doc); // shuffle
		time_end = std::chrono::high_resolution_clock::now();
		shuffle_time += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end-time_start).count();

    for (int ii = 0; ii < num_doc; ii++){
			time_start = std::chrono::high_resolution_clock::now();
      int doc_id = doc_indexes[ii];
      IntegerVector doc_x = X[doc_id], doc_z = Z[doc_id], doc_w = W[doc_id];
			time_end = std::chrono::high_resolution_clock::now();
			time_getvector += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end-time_start).count();

			time_start = std::chrono::high_resolution_clock::now();
      std::vector<int> token_indexes = shuffled_indexes(doc_x.size()); //shuffle
			time_end = std::chrono::high_resolution_clock::now();
			shuffle_time += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end-time_start).count();

      for (int jj = 0; jj < doc_x.size(); jj++){
				time_start = std::chrono::high_resolution_clock::now();
        int w_position = token_indexes[jj];
        int x = doc_x[w_position], z = doc_z[w_position], w = doc_w[w_position];
				time_getvector += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end-time_start).count();

				time_start = std::chrono::high_resolution_clock::now();
        doc_z[w_position] = sample_z(n_x0_kv, n_x1_kv, n_x0_k, n_x1_k, n_dk,
                                     alpha, seed_num, x, z, w, doc_id,
                                     num_vocab, num_topics, k_seeded, gamma_1,
                                     gamma_2, beta, beta_s, phi_s);
				time_end = std::chrono::high_resolution_clock::now();
				time_z += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end-time_start).count();

				z = doc_z[w_position]; // use updated z

				time_start = std::chrono::high_resolution_clock::now();
        doc_x[w_position] = sample_x(n_x0_kv, n_x1_kv, n_x0_k, n_x1_k, n_dk,
                                    alpha, seed_num, x, z, w, doc_id,
                                    num_vocab, num_topics, k_seeded, gamma_1,
                                    gamma_2, beta, beta_s, phi_s);
				time_end = std::chrono::high_resolution_clock::now();
				time_x += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end-time_start).count();

      }

			time_start = std::chrono::high_resolution_clock::now();
			X[doc_id] = doc_x; // is doc_x not a pointer/ref to X[doc_id]?
			Z[doc_id] = doc_z;
			time_end = std::chrono::high_resolution_clock::now();
			time_getvector += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end-time_start).count();
    }
		time_start = std::chrono::high_resolution_clock::now();
    slice_sample_alpha(alpha, n_dk, num_topics, num_doc);
    model["alpha"] = alpha;
		time_end = std::chrono::high_resolution_clock::now();
		time_alpha += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end-time_start).count();

		// Store Alpha
		time_start = std::chrono::high_resolution_clock::now();
		NumericVector alpha_rvec = alpha_reformat(alpha, num_topics);
		List alpha_iter = model["alpha_iter"];
		alpha_iter.push_back(alpha_rvec);
		model["alpha_iter"] = alpha_iter;
		time_end = std::chrono::high_resolution_clock::now();
		rpushback += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end-time_start).count();

		// Log-likelihood and Perplexity
		int r_index = it + 1;
		if(r_index % output_per == 0 || r_index == 1 || r_index == iter ){
			time_start = std::chrono::high_resolution_clock::now();
			double loglik = loglikelihood1(n_x0_kv, n_x1_kv, n_x0_k, n_x1_k, n_dk, alpha,
																		 beta, beta_s, gamma_1, gamma_2,
																		 num_topics, k_seeded, num_vocab, num_doc, phi_s);
			double perplexity = exp(-loglik / (double)total_words);
			time_end = std::chrono::high_resolution_clock::now();
			time_loglik += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end-time_start).count();

			time_start = std::chrono::high_resolution_clock::now();
			NumericVector model_fit_vec;
			model_fit_vec.push_back(r_index);
			model_fit_vec.push_back(loglik);
			model_fit_vec.push_back(perplexity);
			model_fit.push_back(model_fit_vec);
			time_end = std::chrono::high_resolution_clock::now();
			rpushback += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end-time_start).count();

			Rcerr << "[" << r_index << "] log likelihood: " << loglik <<
							 " (perplexity: " << perplexity << ")" << std::endl;

		}

    checkUserInterrupt();
  }

	double end_time = std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::high_resolution_clock::now() - start).count();

	model["model_fit"] = model_fit;
	double devide = 1.0/100000;

	cout << "Preparation inside C++: " << prepare_data * devide  << endl;
	cout << "Sampling Z: " << time_z * devide << endl;
	cout << "      Rcat: " << time_z_rcat * devide << endl;
	cout << "  prep_vec: " << time_z_prepare_vec * devide << endl;
	cout << "          denominator: " << time_z_denom * devide << endl;
	cout << "Sampling X: " << time_x * devide << endl;
	cout << "Sampling alpha: " << time_alpha* devide << endl;
	cout << "Calculation loglik: " << time_loglik* devide << endl;
	cout << "Robj push_back: " << rpushback* devide << endl;
	cout << "Shuffle: " << shuffle_time* devide << endl;
	cout << "Total time: " << end_time* devide << endl;

  return model;
}

