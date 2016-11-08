// ===============================================================
// Example code: how to generate a sample from a covariance matrix
// ===============================================================

#include "Func.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;

int main () {

  try {

    // define a correlation matrix
    vector<vector<double> > correlation = { {1, 0.6, 0.3}, {0.6, 1, 0.5}, {0.3, 0.5, 1} };

    // define the mean and standard deviation of samples
    int npt = 3;
    vector<double> mean = {1., 3., 4.};
    vector<double> stdev = {0.1, 0.3, 0.5};
    vector<vector<double> > covariance(npt, vector<double>(npt, 0));

    for (int i=0; i<npt; i++)
      for (int j=0; j<npt; j++)
	covariance[i][j] = correlation[i][j]*stdev[i]*stdev[j];

    // generate nExtractions correlated samples
    int nExtractions = 10000;
    vector<vector<double> > sample = cosmobl::generate_correlated_data(nExtractions, mean, covariance, 231);

    // measure the covariance matrix from the samples
    vector<vector<double> > measured_covariance;
    cosmobl::covariance_matrix(sample, measured_covariance);
  
    for (int i=0; i<covariance.size(); i++)
      for (int j=0; j<covariance[i].size(); j++)
	cout << i << " " << j << " covariance: " << covariance[i][j] << ", measured covariance: " << measured_covariance[i][j] << endl;
    cout << endl;  
  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}
