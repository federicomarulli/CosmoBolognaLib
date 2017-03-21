// ============================================================================
// Example code: how to minimise a function using wrappers to the GSL libraries
// ============================================================================

#include "Func.h" 

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;


double chi2 (const double xx)
{
  double mean = 0.;
  double sigma = 1.;
  return pow(xx-mean, 2)/(sigma*sigma);
}

double chi2_2 (const vector<double> xx)
{
  double mean1 = 0.;
  double sigma1 = 1.;
  double mean2 = 1.;
  double sigma2 = 1.;
  return pow(xx[0]-mean1, 2)/(sigma1*sigma1)+pow(xx[1]-mean2, 2)/(sigma2*sigma2);
}


int main () {

  try {

    double starting_point = 99.;
    double min = -100.;
    double max = 100.;
    double res = cosmobl::gsl::GSL_minimize_1D(&chi2, starting_point, min, max);
    cout << res << endl;

    vector<double> starting_points = {5, 5};
    vector<double> result = cosmobl::gsl::GSL_minimize_nD(&chi2_2, starting_points);
    cosmobl::print(result);
    
  }
  
  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}
