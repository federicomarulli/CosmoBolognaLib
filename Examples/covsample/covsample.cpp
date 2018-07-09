// ===============================================================
// Example code: how to generate a sample from a covariance matrix
// ===============================================================

#include "Func.h"

int main () {

  try {

    // define a correlation matrix
    std::vector<std::vector<double> > correlation = { {1, 0.6, 0.3}, {0.6, 1, 0.5}, {0.3, 0.5, 1} };

    // define the mean and standard deviation of samples
    int npt = 3;
    std::vector<double> mean = {1., 3., 4.};
    std::vector<double> stdev = {0.1, 0.3, 0.5};
    std::vector<std::vector<double> > covariance(npt, std::vector<double>(npt, 0));

    for (int i=0; i<npt; i++)
      for (int j=0; j<npt; j++)
	covariance[i][j] = correlation[i][j]*stdev[i]*stdev[j];

    // generate nExtractions correlated samples
    int nExtractions = 10000;
    std::vector<std::vector<double> > sample = cbl::generate_correlated_data(nExtractions, mean, covariance, 231);

    // measure the covariance matrix from the samples
    std::vector<std::vector<double> > measured_covariance;
    cbl::covariance_matrix(sample, measured_covariance);
  
    for (size_t i=0; i<covariance.size(); i++)
      for (size_t j=0; j<covariance[i].size(); j++)
	std::cout << i << " " << j << " covariance: " << covariance[i][j] << ", measured covariance: " << measured_covariance[i][j] << std::endl;
    std::cout << std::endl;  
  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}
