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
}


void keyATMhmm::initialize_specific()
{

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

