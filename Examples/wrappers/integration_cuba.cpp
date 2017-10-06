// ===============================================================================================
// Example code: how to integrate a multidimensional function using wrappers to the CUBA libraries
// ===============================================================================================

#include "Func.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;


int main () {

  try {

    // dimension of the integrand function
    const int ndim = 2;

    // function to be integrated
    auto integrand = [] (const vector<double> x)
      {
	double sum = 0.;
	for (size_t i=0; i<x.size(); i++)
	  sum += x[i];
	return sum;
      };
    
    // limits of the integral
    vector<vector<double>> integration_limits(2);
    integration_limits[0] = {0., 2.};
    integration_limits[1] = {0., 2.};

    // wrapper to CUBA libraries
    cosmobl::cuba::CUBAwrapper CW(integrand, ndim);

    // write the outputs
    cout << "The integral computed with Vegas is: " << CW.IntegrateVegas(integration_limits) << endl;
    cout << "The integral computed with Suave is: " << CW.IntegrateSuave(integration_limits) << endl;
    cout << "The integral computed with Divonne is: " << CW.IntegrateDivonne(integration_limits) << endl;
    cout << "The integral computed with Cuhre is: " << CW.IntegrateCuhre(integration_limits) << endl;
  
  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }

  return 0;
}

 
