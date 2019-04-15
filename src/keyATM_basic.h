#ifndef __keyATM_basic__INCLUDED__
#define __keyATM_basic__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class keyATMbasic : public keyATMbase
{
	public:	
		
		//
		// Functions
		//

		// Constructor
		keyATMbasic(List model_, const int iter_, const int output_per_);

		// Read data
		void read_data_specific();
};


#endif

