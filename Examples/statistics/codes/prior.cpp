// ======================================
// Example code: how to construct a prior
// ======================================

#include "PriorDistribution.h"

int main () {

  try {
    
    double xmin = -10., xmax = 10.;
    double mean = -1., sigma = 0.1;
    
    cbl::statistics::PriorDistribution priorDistribution { cbl::glob::DistributionType::_Gaussian_, {mean, sigma}, xmin, xmax };
    
    for (int i=0; i<100; ++i) {
      double value = priorDistribution.sample();
      std::cout << i << "  " << value << "  " << priorDistribution(value) << std::endl;
    }

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

 
