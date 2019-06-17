#include "util.h"

using namespace Rcpp;
using namespace std;

namespace util{

	double mylgamma(double x){
  	double z=1.0/(x*x);

    x=x+6.0;

    z=(((-0.000595238095238*z+0.000793650793651)
			*z-0.002777777777778)*z+0.083333333333333)/x;

    z=(x-0.5)*log(x)-x+0.918938533204673+z-log(x-1.0)-
			log(x-2.0)-log(x-3.0)-log(x-4.0)-log(x-5.0)-log(x-6.0);

    return z;	
	}

}
