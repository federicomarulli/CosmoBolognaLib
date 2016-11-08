// ======================================
// Example code: how to construct a prior
// ======================================

#include "Prior.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;

int main () {

  try {
    
    double xmin = -10., xmax = 10.;
    double mean = -1., sigma = 0.1;
    
    cosmobl::statistics::Prior prior { cosmobl::statistics::PriorType::_GaussianPrior_, {mean, sigma}, xmin, xmax };
    
    for (int i=0; i<100; ++i) {
      double value = prior.sample(i);
      cout << i << "  " << value << "  " << prior(value) << endl;
    }

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}

 
