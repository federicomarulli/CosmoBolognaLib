// ===============================================================================
// Example code: how to generate random numbers from a given distribution function
// ===============================================================================

#include "Distribution.h"
#include "RandomNumbers.h"
#include "Histogram.h"


int main () {

  try {
    
    // Gaussian function
    const cbl::distribution_func func = [] (double x, const std::shared_ptr<void> modelInput, std::vector<double> parameter)
      {
	const std::vector<double> data = *std::static_pointer_cast<std::vector<double>>(modelInput);
	
	return 1./(parameter[1]*data[0])*exp(-data[1]*pow((x-parameter[0])/parameter[1], 2));
      };

    // fixed parameters to construct the Gaussian function
    const std::vector<double> vect = {sqrt(2.*cbl::par::pi), 0.5};
    auto data = std::make_shared<std::vector<double>>(vect);

    // mean, sigma 
    const std::vector<double> parameter = {0., 0.1}; 

    // construct the objects used to draw random numbers from a given distribution function
    cbl::random::CustomDistributionRandomNumbers ran1(func, data, parameter, 666, -10., 10.); // example 1
    cbl::random::CustomDistributionRandomNumbers ran2(&cbl::gaussian<double>, data, parameter, 666, -10., 10.); // example 2
   
    // extract the random numbers
    for (int i=0; i<10; ++i) std::cout << ran1() << "   " << ran1() << std::endl;
      
  }
  
  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}
