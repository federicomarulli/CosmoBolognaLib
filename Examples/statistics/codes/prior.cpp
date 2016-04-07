// ======================================
// Example code: how to construct a prior
// ======================================

#include "Prior.h"

using namespace cosmobl;
using namespace statistics;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;

int main () {

  double xmin = -10., xmax = 10;
  double mean = -1., sigma = 0.1;

  Prior prior { PriorType::_GaussianPrior_, {mean, sigma}, {xmin, xmax} };

  for (int i=0; i<100; i++) {
    double value = prior.sample(i);
    cout << i << "  " << value << "  " << prior(value) << endl;
  }

  return 0;
}

 
