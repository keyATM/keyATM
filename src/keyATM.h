#ifndef __keyATM__INCLUDED__
#define __keyATM__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMbase
{
	public:
		//
		// Parameters
		//
		int iter;
		int output_per;
		double eta_1 = 1;
		double eta_2 = 1;
		double eta_1_regular = 2;
		double eta_2_regular = 1;

		double slice_A = 1.2; // parameter for slice sampling 

		// Data
		List model;
		List W, Z, X;
		StringVector files, vocab;
		NumericVector nv_alpha;
		double gamma_1, gamma_2;
		double beta, beta_s;
		int k_free, k_seeded;
		List seeds;
		List model_fit;
		
		int num_vocab, num_doc, total_words;

		// alpha
		int num_topics;
		VectorXd alpha;

		// Keywords
		std::vector< std::unordered_set<int> > keywords;
		std::vector<int> seed_num;

		// Latent Variables
		MatrixXi n_x0_kv;
		MatrixXi n_x1_kv;
		MatrixXd n_dk;
		VectorXi n_x0_k;
		VectorXi n_x1_k;

		// Track time
		std::chrono::time_point<std::chrono::high_resolution_clock> start;
		double prepare_data;

		//
		// Functions
		//
		keyATMbase(List model_, const int iter_, const int output_per_);

		// Reading and Initialization
		void read_data();
		void read_data_common();
		virtual void read_data_specific() = 0;

		void initialize();
		void initialize_common();
		virtual void initialize_specific() = 0;

		// Sampling
		void iteration();
		virtual void iteration_single() = 0;
		virtual void sample_parameters() = 0;

		int sample_z(VectorXd &alpha, int &z, int &x,
									 int &w, int &doc_id);
		int sample_x(VectorXd &alpha, int &z, int &x,
									 int &w, int &doc_id);

		void sampling_store(int &r_index);
		virtual double loglik_total() = 0;

		// Utilities
		double logsumexp(double &x, double &y, bool flg);
		double logsumexp_Eigen(VectorXd &vec);
		double gammapdfln(const double x, const double a, const double b);
		NumericVector alpha_reformat(VectorXd& alpha, int& num_topics);

		double expand(double &p, const double &A);
		double shrink(double &x, const double &A);

		List return_model();
	
};

#endif
