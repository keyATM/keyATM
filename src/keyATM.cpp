#include "keyATM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */

keyATMbase::keyATMbase(List model_, const int iter_, const int output_per_)
{
	model = model_;
	iter = iter_;
	output_per = output_per_;
}

void keyATMbase::read_data()
{
	read_data_common();
	read_data_specific();
}

void keyATMbase::read_data_common()
{
	// Read data
}
