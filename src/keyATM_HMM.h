#ifndef __keyATM_HMM__INCLUDED__
#define __keyATM_HMM__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMhmm : public keyATMbase
{
	public:
		// Parameters
		int num_states;
		int index_states;  // num_states - 1
		MatrixXd Psk;
		VectorXi S_est;
		MatrixXd P_est;
	
		// Constructor
		keyATMhmm(List model_, const int iter_, const int output_per_);
	
		// 
		// Functions
		//
	
		// Read data and Initialize
		void read_data_specific();
		void initialize_specific();
	
		// Iteration
		void iteration_single();
		void sample_parameters();
		double loglik_total();
};

#endif

