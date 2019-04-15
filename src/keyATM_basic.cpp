#include "keyATM_basic.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;


keyATMbasic::keyATMbasic(List model_, const int iter_, const int output_per_) :
	keyATMbase(model_, iter_, output_per_) // pass to parent!
{
	// Constructor
	read_data();
}


void keyATMbasic::read_data_specific()
{
	cout << "specific" << endl;
	cout << iter << endl;
}
