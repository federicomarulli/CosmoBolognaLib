// =============================================================================
// Example code: how to integrate a function using wrappers to the GSL libraries
// =============================================================================

#include "Func.h" 

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;


double Gaussian (const double xx)
{
  const double mean = 0.;
  const double sigma = 1.;
  return 1./(sigma*sqrt(2.*cosmobl::par::pi))*exp(-pow(xx-mean, 2)/(2.*sigma*sigma));
}


int main () {

  try {

    const double mean = 0.;
    const double sigma = 1.;
    const shared_ptr<void> pp;

    // lower and upper limits of the integration
    const double lower_limit = -1.;
    const double upper_limit = 1.;

    // relative measurement precision (it raises a warning message if not satisfied)
    const double prec = 1.e-5;

    // default values (refer to the GSL documentation for different choiches)
    const int limit_size = 1000;
    const int rule = 6;

    // integrating a function of type: function<double<double, shared_ptr<void>, vector<double>>
    cout << "The integral is: " << cosmobl::gsl::GSL_integrate_qag(&cosmobl::gaussian<double>, pp, {mean, sigma}, lower_limit, upper_limit, prec, limit_size, rule) << endl;

    // integrating a function of type: function<double<double>>
    cout << "The integral is: " << cosmobl::gsl::GSL_integrate_qag(&Gaussian, lower_limit, upper_limit, prec, limit_size, rule) << endl;

  }
  
  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}
