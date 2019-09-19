// ========================================================================================
// Example code: how to generate correlated random samples using the Cholesky decomposition
// ========================================================================================

#include "Func.h" 

int main () {

  try {
    
    // the mean values
    const std::vector<double> mean = {1., 2.};

    // the standard deviation values
    const std::vector<double> std = {0.1, 0.3};

    // the sample correlation
    const double correlation_factor = 0.5;

    // the covariance matrix
    std::vector<std::vector<double>> covariance = {{std[0]*std[0], std[0]*std[1]*correlation_factor},
                                                   {std[0]*std[1]*correlation_factor, std[1]*std[1]}};
    
    // the seed
    int seed = 666;

    // the number of samples to be extracted
    int nExtractions = 100000;

    // construct the object to generate correlated samples; dataset has size nExtractions x nData
    std::vector<std::vector<double>> dataset = cbl::generate_correlated_data(nExtractions, mean, covariance, seed);

    // traspose it, to simplify the output verification
    std::vector<std::vector<double>> datasetT = cbl::transpose(dataset);

    // write data statistics
    for (size_t i=0; i<datasetT.size(); i++) {
      std::cout << "Data" << i+1 << " mean = " << mean[i] << "  " <<  cbl::Average(datasetT[i]) << std::endl;
      std::cout << "     " << " std  = " << std[i] << "  " << cbl::Sigma(datasetT[i]) << std::endl;
    }

    // write the input covariance
    std::cout << std::endl << "Input covariance : " << std::endl;
    cbl::Print(covariance);

    // write the output
    std::cout << std::endl << "Measured covariance : " << std::endl;
    std::vector<std::vector<double>> measured_covariance;
    cbl::covariance_matrix (dataset, measured_covariance);
    cbl::Print(measured_covariance);

  }
  
  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}
