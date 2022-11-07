// =================================================================================
// Example code: how to use the NormalRandomNumbers class to generate random numbers
// =================================================================================

#include "Distribution.h" 

int main () {

  try {
  
    // choose a seed
    int seed = 1231;

    // choose a range 
    double xmin = 0.;
    double xmax = 10.;

    // choose the parameters for the Normal distribution
    double mean = 5.;
    double sigma = 2.;

    // construct the object used to extract random numbers from a Normal distribution
    cbl::random::NormalRandomNumbers ran(mean, sigma, seed, xmin, xmax);

    // extract the random numbers
    std::vector<double> numbers(1000);
    for (auto && num : numbers) num = ran();

    // compute the distribution of the extracted random numbers
    std::vector<double> xx, fx, error;
    cbl::distribution(xx, fx, error, numbers, {}, 10);
  
    // show the results
    for (size_t i=0; i<xx.size(); ++i)
      std::cout << std::setprecision(1) << std::fixed << xx[i] << "  " << std::string(fx[i]*100./numbers.size(), '*') << std::endl;

  }
  
  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}
