#include "idealpoint.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */


// 今後の改善の方針
// 今は直にalphaを決めてしまっているけれど、a0と掛け合わせる形にして
// (CSTMみたいに)、トピックごとの標準の登場確率$alpha_0$を、
// イデオロギーの位置$alpha$でずらす形にするのが良いのでは。
//
// alpha_d = alpha0 * exp(...ここに今のがくる...)


IDEALPOINT::IDEALPOINT(List model_, List author_info_, const int iter_, const int output_iter_)
{
	// Get Info from R
	model = model_;

	// Speaker Info
	author_info = author_info_;
	author_ids = Rcpp::as<vector<int>> (author_info["author_ids"]);
	num_authors = Rcpp::as<int> (author_info["num_authors"]);
	num_authordoc = Rcpp::as<vector<int>> (author_info["num_authordoc"]);

	sigma_lambda = Rcpp::as<double> (author_info["sigma_lambda"]);
	sigma_psi = Rcpp::as<double> (author_info["sigma_psi"]);

	W = model["W"];
	Z = model["Z"];
	X = model["X"];
  gamma_1 = model["gamma_1"];
	gamma_2 = model["gamma_2"];

  beta = model["beta"];
	beta_s = model["beta_s"];
  k_free = model["extra_k"];
  seeds = model["seeds"];
  k_seeded = seeds.size();
	List model_fit = model["model_fit"];

	iter = iter_;
	output_iter = output_iter_;
	num_topics = k_seeded + k_free;
	
	vocab = model["vocab"];
	num_vocab = vocab.size();
	StringVector files = model["files"];
	num_doc = files.size();
	
	// Initialize
	initialize();
	
	// Iteration
	iteration();
	
}

void IDEALPOINT::initialize()
{
	sigma_lambda = 1.0;
	sigma_psi = 1.0;

	// Psi
	Psi = VectorXd::Zero(num_authors);
	for(int j=0; j<num_authors; j++){
		Psi(j) = R::rnorm(0.0, sigma_psi);
	}
	
	
	// Lambda
	Lambda = MatrixXd::Zero(num_authors, num_topics);
	for(int j=0; j<num_authors; j++){
		// Initialize with R random
		for(int k=0; k<num_topics; k++){
			Lambda.coeffRef(j, k) = R::rnorm(Psi(j), sigma_lambda);
		}
	}

	// Vector that stores seed words (words in dictionary)
  for (int ii = 0; ii < k_seeded; ii++){
    IntegerVector wd_ids = seeds[ii];
    seed_num.push_back(wd_ids.size());
		
    std::unordered_set<int> keywords_set;
    for (int jj = 0; jj < wd_ids.size(); jj++){
			int wd_id = wd_ids(jj);
			keywords_set.insert(wd_id);
		}
		
    keywords.push_back(keywords_set);
  }

	for(int i=k_seeded; i<num_topics; i++){
		std::unordered_set<int> keywords_set{ -1 };
	
		seed_num.push_back(0);
		keywords.push_back(keywords_set);
	}
	
	
  // storage for sufficient statistics and their margins
  n_x0_kv = MatrixXi::Zero(num_topics, num_vocab);
  n_x1_kv = MatrixXi::Zero(num_topics, num_vocab);
  n_dk = MatrixXd::Zero(num_doc, num_topics);
  n_x0_k = VectorXi::Zero(num_topics);
  n_x1_k = VectorXi::Zero(num_topics);
	
	
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
  total_words = (int)n_dk.sum();

}


void IDEALPOINT::iteration()
{
	// Iteration
	VectorXd alpha_ = VectorXd::Zero(num_topics);
	MatrixXd Alpha;
	int author_id;
	
	for(int it=0; it<iter; it++){
		std::vector<int> doc_indexes = shuffled_indexes(num_doc);
		
		Alpha = Lambda.array().exp();
		
    for (int ii = 0; ii < num_doc; ii++){
      int doc_id_ = doc_indexes[ii];
			author_id = author_ids[doc_id_];
		
      IntegerVector doc_x = X[doc_id_], doc_z = Z[doc_id_], doc_w = W[doc_id_];
			
      std::vector<int> token_indexes = shuffled_indexes(doc_z.size()); //shuffle
			
			// Prepare Alpha for the doc
			alpha_ = Alpha.row(author_id).transpose(); // take out alpha
			
			// Iterate each word in the document
      for (int jj = 0; jj < doc_z.size(); jj++){
        int w_position = token_indexes[jj];
        int x_ = doc_x[w_position], z_ = doc_z[w_position], w_ = doc_w[w_position];
			
				int new_z = sample_z(alpha_, z_, x_, w_, doc_id_);
        doc_z[w_position] = new_z;
			
		
				z_ = doc_z[w_position]; // use updated z
				int new_x = sample_x(alpha_, z_, x_, w_, doc_id_);
				doc_x[w_position] = new_x;
      }
			
			Z[doc_id_] = doc_z;
			X[doc_id_] = doc_x;
    }
		lambda_sample();
		lambda_store();

		psi_sample();
	
		int r_index = it + 1;
		if(r_index % output_iter == 0 || r_index == 1 || r_index == iter ){
			loglik_store(r_index);
		}


		checkUserInterrupt();
	}
}


void IDEALPOINT::psi_sample()
{
	double new_mu;
	double new_sigma;
	double new_psi = 0.0;
	double mu = 0;

	for(int j=0; j<num_authors; j++){
		// int N = num_authordoc[j];
		int N = num_topics;
		double sigma_c2 = pow(sigma_lambda, 2);
		double sigma_p2 = pow(sigma_psi, 2);

		// new_mu
		new_mu = (sigma_c2 * mu + sigma_p2 * Lambda.row(j).sum()) / (sigma_c2 + N * sigma_p2);

		// new_sigma
		new_sigma = (sigma_c2 * sigma_p2) / (sigma_c2 + N * sigma_p2);
			// 1/(numerator/denominator) = denominator/numerator

		// Sample
		new_psi = R::rnorm(new_mu, new_sigma);

		// Store
		Psi(j) = new_psi;	
	}

	// Save
	NumericVector Psi_vec = Rcpp::wrap(Psi);
	List Psi_iter = author_info["Psi_iter"];
	Psi_iter.push_back(Psi_vec);
	author_info["Psi_iter"] = Psi_iter;

}

double IDEALPOINT::loglik_lambda(int &author_id)
{
	double loglik = 0.0;
	double mu = Psi(author_id);
	VectorXd alpha = Lambda.row(author_id).transpose().array().exp();

	for(int d=0; d<num_doc; d++){

		if (author_id != author_ids[d])
			// Calculate only author_id matches
			continue;

		// alpha = Alpha.row(author_id).transpose(); // Doc alpha, column vector

		loglik += lgamma(alpha.sum()); 
				// the first term numerator in the first square bracket
		loglik -= lgamma( n_dk.row(d).sum() + alpha.sum() ); 
				// the second term denoinator in the first square bracket

		for(int k=0; k<num_topics; k++){
			loglik -= lgamma(alpha(k));
				// the first term denominator in the first square bracket
			loglik += lgamma( n_dk(d, k) + alpha(k) );
				// the second term numerator in the firist square bracket
		}
	}

	// Prior
	double prior_fixedterm = -0.5 * log(2.0 * PI_V * std::pow(sigma_lambda, 2.0) );
	for(int k=0; k<num_topics; k++){
		loglik += prior_fixedterm;

		loglik -= ( std::pow( (Lambda(author_id,k) - mu) , 2.0) / (2.0 * std::pow(sigma_lambda, 2.0)) );
	}

	return loglik;
}

void IDEALPOINT::lambda_sample()
{

	// Slice sampling for Lambda

	double start, end, previous_p, new_p, newlikelihood, slice_, current_lambda;
  std::vector<int> topic_ids = shuffled_indexes(num_topics);
	std::vector<int> authors_shuffled = shuffled_indexes(num_authors);
	static double A = 0.7; // important, 0.5-1.5
	static int max_shrink_time = 1000;

	double start_min = -8.0;
	double end_max = 8.0;

	double newlambdallk = 0.0;

	for(int aa=0; aa<num_authors; aa++){
		int author_id = authors_shuffled[aa];	


		for(int kk=0; kk<num_topics; kk++){
			int k = topic_ids[kk];

			start = shrink(start_min, A); // shrink
			end = shrink(end_max, A); // shrink

			current_lambda = Lambda(author_id, k);
			previous_p = shrink(current_lambda, A);
			slice_ = loglik_lambda(author_id) - std::log(A * previous_p * (1.0 - previous_p)) 
							+ log(unif_rand()); // <-- using R random uniform


			for (int shrink_time = 0; shrink_time < max_shrink_time; shrink_time++){
				new_p = sampler::slice_uniform(start, end); // <-- using R function above
				Lambda(author_id, k) = expand(new_p, A); // expand

				newlambdallk = loglik_lambda(author_id);

				newlikelihood = newlambdallk - std::log(A * new_p * (1.0 - new_p));

				if (slice_ < newlikelihood){
					break;
				} else if (previous_p < new_p){
					end = new_p;
				} else if (new_p < previous_p){
					start = new_p;
				} else {
					Rcerr << "Something goes wrong in sample_lambda_slice()" << std::endl;
					Lambda(author_id, k) = current_lambda;
					break;
				}

			} // for loop for shrink time

		} // for loop for num_cov
	} // for loop for num_topics

	model["model_fit"] = model_fit;

}

void IDEALPOINT::loglik_store(int& r_index)
{
	double loglik = loglik_calc();
	double perplexity = exp(-loglik / (double)total_words);

	NumericVector model_fit_vec;
	model_fit_vec.push_back(r_index);
	model_fit_vec.push_back(loglik);
	model_fit_vec.push_back(perplexity);
	model_fit.push_back(model_fit_vec);

	Rcerr << "[" << r_index << "] log likelihood: " << loglik <<
					 " (perplexity: " << perplexity << ")" << std::endl;
}

double IDEALPOINT::loglik_calc()
{
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

    // x normalization
    loglik += lgamma(gamma_1 + gamma_2) - lgamma(gamma_1) - lgamma(gamma_2);
  }


  // z
	MatrixXd Alpha = Lambda.array().exp();
	VectorXd alpha = VectorXd::Zero(num_topics);
	int author_id;

  for (int d = 0; d < num_doc; d++){
		author_id = author_ids[d];
		alpha = Alpha.row(author_id).transpose(); // Doc alpha, column vector	

    loglik += lgamma( alpha.sum() ) - lgamma( n_dk.row(d).sum() + alpha.sum() );
    for (int k = 0; k < num_topics; k++){
      loglik += lgamma( n_dk(d,k) + alpha(k) ) - lgamma( alpha(k) );
    }
  }
  return loglik;
}

void IDEALPOINT::lambda_store()
{
	Rcpp::NumericMatrix Lambda_R = Rcpp::wrap(Lambda);
	List Lambda_iter = model["Lambda_iter"];
	Lambda_iter.push_back(Lambda_R);
	model["Lambda_iter"] = Lambda_iter;
}

int IDEALPOINT::sample_z(Eigen::VectorXd &alpha, int &z, int &x,
				  	     int &w, int &doc_id)
{
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
    for (int k = 0; k < num_topics; ++k){
      numerator = (beta + (double)n_x0_kv(k, w)) *
        ((double)n_x0_k(k) + gamma_2) *
        ((double)n_dk(doc_id, k) + alpha(k));

      denominator = ((double)num_vocab * beta + (double)n_x0_k(k)) *
        ((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);

      z_prob_vec(k) = numerator / denominator;
    }

    double sum = z_prob_vec.sum(); // normalize
    new_z = rcat_without_normalize(z_prob_vec, sum); // take a sample

  } else {
    for (int k = 0; k < num_topics; ++k){
      if (keywords[k].find(w) == keywords[k].end()){
        z_prob_vec(k) = 0.0;
				continue;
      } else{ // w not one of the seeds
        numerator = (beta_s + (double)n_x1_kv(k, w)) *
          ( ((double)n_x1_k(k) + gamma_1) ) *
          ( ((double)n_dk(doc_id, k) + alpha(k)) );
      }
      denominator = ((double)seed_num[k] * beta_s + (double)n_x1_k(k) ) *
        ((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);

      z_prob_vec(k) = numerator / denominator;
    }


		double sum = z_prob_vec.sum();
    new_z = rcat_without_normalize(z_prob_vec, sum); // take a sample

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

int IDEALPOINT::sample_x(Eigen::VectorXd &alpha, int &z, int &x,
				  	     int &w, int &doc_id)
{
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
  if ( keywords[k].find(w) == keywords[k].end() ){
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

int IDEALPOINT::rcat_without_normalize(Eigen::VectorXd &prob,
		double &total)
{ // was called 'multi1'
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

vector<int> IDEALPOINT::shuffled_indexes(int m) {
  vector<int> v(m);
  iota(v.begin(), v.end(), 0);
  random_shuffle(v.begin(), v.end(), sampler::rand_wrapper);
  return v;
}


double IDEALPOINT::expand(double& p, double& A){
	double res = -(1.0/A) * log((1.0/p) - 1.0);
	return res;
}

double IDEALPOINT::shrink(double& x, double& A){
	double res = 1.0 / (1.0 + exp(-A*x));
	return res;
}
