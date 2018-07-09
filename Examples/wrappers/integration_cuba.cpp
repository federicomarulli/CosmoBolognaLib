// ===============================================================================================
// Example code: how to integrate a multidimensional function using wrappers to the CUBA libraries
// ===============================================================================================

#include "Func.h"

int main () {

  try {

    // dimension of the integrand function
    const int ndim = 2;

    // function to be integrated
    auto integrand = [] (const std::vector<double> x)
      {
	double sum = 0.;
	for (size_t i=0; i<x.size(); i++)
	  sum += x[i];
	return sum;
      };
    
    // limits of the integral
    std::vector<std::vector<double>> integration_limits(2);
    integration_limits[0] = {0., 2.};
    integration_limits[1] = {0., 2.};

    // wrapper to CUBA libraries
    cbl::cuba::CUBAwrapper CW(integrand, ndim);

    // write the outputs
    std::cout << "The integral computed with Vegas is: " << CW.IntegrateVegas(integration_limits) << std::endl;
    std::cout << "The integral computed with Suave is: " << CW.IntegrateSuave(integration_limits) << std::endl;
    std::cout << "The integral computed with Divonne is: " << CW.IntegrateDivonne(integration_limits) << std::endl;
    std::cout << "The integral computed with Cuhre is: " << CW.IntegrateCuhre(integration_limits) << std::endl;
  
  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }

  return 0;
}

 
