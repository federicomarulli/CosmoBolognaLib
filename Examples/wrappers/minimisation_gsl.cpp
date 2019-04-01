// ============================================================================
// Example code: how to minimise a function using wrappers to the GSL libraries
// ============================================================================

#include "Func.h" 

int main () {

  try {

    // define the (lambda) function to minimise
    
    auto chi2 = [] (const double xx) 
    {
      const double mean = 0.;
      const double sigma = 1.;
      return pow(xx-mean, 2)/(sigma*sigma);
    };

    
    // run the 1D minimisation GSL wrapper
    
    const double starting_point = 99.;
    const double min = -100.;
    const double max = 100.;
    const double res = cbl::wrapper::gsl::GSL_minimize_1D(chi2, starting_point, min, max);
    std::cout << "result of the minimisation: " << res << std::endl;

    
    // define a 2D (lambda) function to minimise

    auto chi2_2 = [] (std::vector<double> &xx)
   { 
      const double mean1 = 0.;
      const double sigma1 = 1.;
      const double mean2 = 1.;
      const double sigma2 = 1.;
      return pow(xx[0]-mean1, 2)/(sigma1*sigma1)+pow(xx[1]-mean2, 2)/(sigma2*sigma2);
    };

    
    // run the nD minimisation GSL wrapper
    
    std::vector<double> starting_points = {5, 5};
    std::vector<std::vector<double>> point_limits = {{-10, 10}, {-10, 10}};
    const int maxiter = 10000;
    std::vector<double> result = cbl::wrapper::gsl::GSL_minimize_nD(chi2_2, starting_points, point_limits, maxiter);
    cbl::Print(result);
    
  }
  
  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}
