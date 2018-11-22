#ifndef __IDEALPOINT__INCLUDED__
#define __IDEALPOINT__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"

class IDEALPOINT
{
	public:
		// Variables From R
		Rcpp::List model;
		Rcpp::List author_info;
		Rcpp::List W, Z, X;
		Rcpp::StringVector vocab;

		double gamma_1, gamma_2;
	 	double beta;
		double beta_s;

		int k_free;
		Rcpp::List seeds;
		int k_seeded;
		Rcpp::List model_fit;

		std::vector<int> author_ids;
		std::vector<int> num_authordoc;
		int num_authors; // number of authors

		std::vector< std::unordered_set<int> > keywords;
		std::vector<int> seed_num;

		int num_doc;
		int num_vocab;
		int num_topics;
		int num_cov;
		int total_words;
		int iter;
		int output_iter;
		std::vector<int> mh_info{0,0};

		// Ideal point related parameters
		double sigma_lambda = 1.0;  // sigma_c
		double sigma_psi = 1.0;  // sigma_p

		// Variables for the Model
	  Eigen::VectorXd Psi;  // speakers' ideal point

		Eigen::MatrixXi n_x0_kv;	
		Eigen::MatrixXi n_x1_kv;	
		Eigen::MatrixXd n_dk;
	  Eigen::VectorXi n_x0_k;
	  Eigen::VectorXi n_x1_k;
		Eigen::MatrixXd Lambda;
		
		// Constructor
		IDEALPOINT(Rcpp::List model_, Rcpp::List author_info_,
				const int iter_, const int output_iter_);

		// Main functions
		void initialize();
		void iteration();
		int sample_z(Eigen::VectorXd &alpha, int &z, int &x,
				  	     int &w, int &doc_id);
		int sample_x(Eigen::VectorXd &alpha, int &z, int &x,
				  	     int &w, int &doc_id);
		void lambda_sample();

		// Store functions
		void lambda_store();
		void loglik_store(int& r_index);
		double loglik_calc();
		double loglik_lambda(int &author_id);
		void psi_sample();

		// Sub functions
		std::vector<int> shuffled_indexes(int m);
		int rcat_without_normalize(Eigen::VectorXd &prob, double &total);
		double expand(double& p, double& A);
		double shrink(double& x, double& A);
};

#endif
