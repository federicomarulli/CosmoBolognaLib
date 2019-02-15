// =============================================================================
// Example code: how to integrate a function using wrappers to the GSL libraries
// =============================================================================

#include "Func.h" 


double Gaussian (const double xx)
{
  const double mean = 0.;
  const double sigma = 1.;
  return 1./(sigma*sqrt(2.*cbl::par::pi))*exp(-pow(xx-mean, 2)/(2.*sigma*sigma));
}


int main () {

  try {

    const double mean = 0.;
    const double sigma = 1.;
    const std::shared_ptr<void> pp;

    // lower and upper limits of the integration
    const double lower_limit = -1.;
    const double upper_limit = 1.;

    // relative measurement precision (it raises a warning message if not satisfied)
    const double prec = 1.e-5;

    // default values (refer to the GSL documentation for different choiches)
    const int limit_size = 1000;
    const int rule = 6;

    // integrating a function of type: function<double<double, std::shared_ptr<void>, std::vector<double>>
    std::cout << "The integral is: " << cbl::gsl::GSL_integrate_qag(&cbl::gaussian<double>, pp, {mean, sigma}, lower_limit, upper_limit, prec, limit_size, rule) << std::endl;

    // integrating a function of type: function<double<double>>
    std::cout << "The integral is: " << cbl::gsl::GSL_integrate_qag(&Gaussian, lower_limit, upper_limit, prec, limit_size, rule) << std::endl;

    // integrating using a lambda function
    std::function<double(double)> GaussianL = [&] (const double xx) { return 1./(sigma*sqrt(2.*cbl::par::pi))*exp(-pow(xx-mean, 2)/(2.*sigma*sigma)); };
    std::cout << "The integral is: " << cbl::gsl::GSL_integrate_qag(GaussianL, lower_limit, upper_limit, prec, limit_size, rule) << std::endl;
  
  }
  
  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}
