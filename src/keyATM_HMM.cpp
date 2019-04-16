#include "keyATM_HMM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;


keyATMhmm::keyATMhmm(List model_, const int iter_, const int output_per_) :
	keyATMbase(model_, iter_, output_per_) // pass to parent!
{
	// Constructor
	// read_data();
	// initialize();
	// iteration();
}
