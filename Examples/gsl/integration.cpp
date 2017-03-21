// =============================================================================
// Example code: how to integrate a function using wrappers to the GSL libraries
// =============================================================================

#include "Func.h" 

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;


double fixed_parameters_gaussian (const double xx)
{
  double mean = 0.;
  double sigma = 1.;
  return 1./(sigma*sqrt(2.*cosmobl::par::pi))*exp(-pow(xx-mean, 2)/(2.*sigma*sigma));
}


int main () {

  try {

    double mean = 0.;
    double sigma = 1.;
    shared_ptr<void> pp;

    // lower and upper limits of integration
    double lower_limit=-1;
    double upper_limit=1;

    // relative measurement precision (it raise a warning message if not satisfacted)
    double prec = 1.e-5;

    // default values (refer to the GSL documentation for different choiches)
    int limit_size = 1000;
    int rule = 6;

    // integrating a function of type: function<double<double, shared_ptr<void>, vector<double>>
    double Int = cosmobl::gsl::GSL_integrate_qag(&cosmobl::gaussian<double>, pp, {mean, sigma}, lower_limit, upper_limit, prec, limit_size, rule);
    cout << "The integral is: " << Int << endl;

    // integrating a function of type: function<double<double>>
    Int = cosmobl::gsl::GSL_integrate_qag(&fixed_parameters_gaussian, lower_limit, upper_limit, prec, limit_size, rule);
    cout << "The integral is: " << Int << endl;

  }
  
  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}
