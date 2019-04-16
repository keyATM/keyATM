#ifndef __keyATM_cov__INCLUDED__
#define __keyATM_cov__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMcov : public keyATMbase
{
	public:	
		//
		// Parameters
		//
		MatrixXd Alpha;
		int num_cov;
		MatrixXd Lambda;
		MatrixXd C;

		// Slice Sampling
		double min_v = 1e-9;
		double max_v = 100.0;
		int max_shrink_time = 1000;

		// MH sampling
		double mu = 0.0;
		double sigma = 50.0;
		double mh_sigma = 0.05;

		// Sampling info
		std::vector<int> mh_info{0,0};
		
		//
		// Functions
		//

		// Constructor
		keyATMcov(List model_, const int iter_, const int output_per_);

		// Read data
		void read_data_specific();

		// Initialization
		void initialize_specific();

		// Iteration
		void iteration_single();
		void sample_parameters();
		void sample_lambda();
		double alpha_loglik();
		double loglik_total();

		void sample_lambda_mh_single();
		void sample_lambda_mh();
		void sample_lambda_slice();
		double likelihood_lambda();
		void proposal_lambda(int& k);

};


#endif


