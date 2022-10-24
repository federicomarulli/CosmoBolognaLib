// ============================================
// Example code: how to use the class Histogram
// ============================================

#include "Distribution.h" 
#include "Histogram.h" 

int main () {

  try {
  
    // choose a seed
    int seed = 1231;

    // choose a range 
    double xmin = 0.;
    double xmax = 10.;

    // choose the parameters for a Normal distribution
    double mean = 5.;
    double sigma = 2.;

    // construct the object used to extract random numbers from the Normal distribution
    cbl::random::NormalRandomNumbers ran(mean, sigma, seed, xmin, xmax);

    // extract the random numbers
    std::vector<double> numbers(1000);
    for (auto && num : numbers) num = ran();

    // construct the object used to extract random numbers from an Uniform distribution
    cbl::random::UniformRandomNumbers uniran(0.9, 1.1, seed);

    // set a dummy weight vector
    std::vector<double> weights(1000);
    for (auto && num : weights) num = uniran();

    // create an histogram and write the output
    cbl::glob::Histogram1D hist(numbers, weights, 10);
    hist.write("./", "histogram.dat", cbl::glob::HistogramType::_N_V_);

  }
  
  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}
