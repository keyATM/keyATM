#include "keyATM_HMM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */

keyATMhmm::keyATMhmm(List model_, const int iter_, const int output_per_) :
  keyATMbase(model_, iter_, output_per_) // pass to parent!
{
	// Constructor
	read_data();
	initialize();
	iteration();
}


void keyATMhmm::read_data_specific()
{
	num_states = model["num_states"];
	index_states = num_states - 1;
}


void keyATMhmm::initialize_specific()
{
	// Initialize Psk
	Psk = MatrixXd::Zero(num_doc, num_states);

	// Initialize S_est
	VectorXi S_est_num = VectorXi::Zero(num_states);
	VectorXd S_est_temp = VectorXd::Zero(num_states);
	double cumulative = 1.0 / num_states;
	double u;
	int index;
	for(int i=0; i<num_states; i++){
		S_est_temp(i) = cumulative * (i+1);
	}

	for(int j=0; j<num_doc; j++){
		u = R::runif(0, 1);
		for(int i=0; i<num_states; i++){
			if(u < S_est_temp(i)){
				index = i;
				break;
			}
		}
		S_est_num(index) += 1;
	}

	S_est = VectorXi::Zero(num_doc);
	int count;
	index = 0;
	for(int i=0; i<num_states; i++){
		count = S_est_num(i);
		for(int j=0; j<count; j++){
			S_est(index) = i;
			index += 1;
		}
	}

	// Initializae P_est
	P_est = MatrixXd::Zero(num_states, num_states);
	double prob;
	for(int i=0; i<=(index_states-1); i++){
		prob = R::rbeta(1.0, 1.0);
		P_est(i, i) = prob;
		P_est(i, i+1) = 1-prob;
	}
	P_est(index_states, index_states) = 1;
	
}


void keyATMhmm::iteration_single()
{

}


void keyATMhmm::sample_parameters()
{

}


double keyATMhmm::loglik_total()
{
	return 0.0;
}

