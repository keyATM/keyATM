#include "keyATM_HMM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;


keyATMhmm::keyATMhmm(List model_, const int iter_, const int output_per_) :
	keyATMbase(model_, iter_, output_per_) // pass to parent!
{
	// Constructor
	read_data();
	// initialize();
	// iteration();
}


void keyATMhmm::read_data_specific()
{
	num_states = model["num_states"];
	cout << num_states << endl;
}


void keyATMhmm::initialize_specific()
{

}


void keyATMhmm::sample_parameters()
{

}


double keyATMhmm::loglik_total()
{
	return 0.0;
}
