#ifndef __keyATM_tot__INCLUDED__
#define __keyATM_tot__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMtot : public keyATMbase
{
	public:
		// Parameters
		VectorXd timestamps;  // time stamps (document share the same time)
		MatrixXd beta_params;  // parameter for time Beta, K \times 2

		// Constructor
		keyATMtot(List model_, const int iter_, const int output_per_);

		// During sampling
			// Sample Beta parameters
			vector<vector<double>> store_t;  
				// store time stamp to estiamte Beta parameters
			VectorXd timestamps_k;  // time stamps for topic k


			// Sample alpha

				// Slice Sampling
				double min_v = 1e-9;
				double max_v = 100.0;
				int max_shrink_time = 1000;

				double start, end, previous_p, new_p, newlikelihood, slice_;
				std::vector<int> topic_ids;
				VectorXd keep_current_param;
				double store_loglik;
				double newalphallk;

				double loglik;
				double fixed_part;

					// in alpha_loglik
					MatrixXd ndk_a;

				// In sampling z
				double beta_a;
				double beta_b;
				double timestamp_d;
	
		// 
		// Functions
		//
	
		// Read data and Initialize
		// Read data
		void read_data_specific();

		// Initialization
		void initialize_specific();

		// Iteration
		void iteration_single(int &it);
		int sample_z(VectorXd &alpha, int &z, int &x, int &w, int &doc_id);
		void sample_parameters();
		void sample_betaparam();
		void sample_alpha();
		double alpha_loglik();
		double loglik_total();
};

#endif


