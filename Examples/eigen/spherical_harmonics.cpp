// ===================================================================
// Test the differences in performances using Eigen std::vectorization
// ===================================================================

#include "RandomNumbers.h"
#include "SphericalHarmonics_Coefficients.h"

int main () {
  

  try {
   
    // create random points in a cubic volume
    
    const double BoxSize = 500.;
    const long nran = 10000;
    const int seed = 666;
    
    cbl::random::UniformRandomNumbers random(0, BoxSize, seed);

    std::vector<double> x(nran), y(nran), z(nran);

    for (long i=0; i<nran; ++i) {
      x[i] = random();
      y[i] = random();
      z[i] = random();
      //v[i] = {x[i], y[i], z[i], 0.};
    }

    time_t timer_start, timer_end;

    time(&timer_start); 

    std::cout << "Distance..." << std::endl;

    double rr;
    for (long i=0; i<nran; ++i) {
      for (long j=0; j<nran; ++j) {

	const double xx = x[i]-x[j];
	const double yy = y[i]-y[j];
	const double zz = z[i]-z[j];
	rr = sqrt(xx*xx+yy*yy+zz*zz);
      }
    }

    time(&timer_end); 
    std::cout << "Time for distances: " << difftime( timer_end, timer_start) << " seconds" << " " << rr << std::endl;

    
    // Using standard method

    time_t timer_standard_start, timer_standard_end;

    time(&timer_standard_start); 

    std::cout << "Starting the standard method..." << std::endl;

    auto sph_standard = cbl::glob::SphericalHarmonicsArray::factory(cbl::glob::SpHarMethod::_STANDARD_, 20);
    (void)sph_standard;

    for (long i=0; i<nran; ++i) {
      for (long j=i+1; j<nran; ++j) {

	const double xx = x[i]-x[j];
	const double yy = y[i]-y[j];
	const double zz = z[i]-z[j];
	double rr = sqrt(xx*xx+yy*yy+zz*zz);

	sph_standard->compute(xx/rr, yy/rr, zz/rr, 1.);
	(void)rr;
	
      }
    }

    sph_standard->print();

    time(&timer_standard_end); 

    std::cout << "Standard approach: " << difftime( timer_standard_end, timer_standard_start) << " seconds" << std::endl;

    // Using GSL
   
    std::cout << "Starting the GSL method..." << std::endl;

    auto sph_gsl = cbl::glob::SphericalHarmonicsArray::factory(cbl::glob::SpHarMethod::_GSL_, 20);

    time(&timer_standard_start); 

    for (long i=0; i<nran; ++i) {
      for (long j=i+1; j<nran; ++j) {

	const double xx = x[i]-x[j];
	const double yy = y[i]-y[j];
	const double zz = z[i]-z[j];
	double rr = sqrt(xx*xx+yy*yy+zz*zz);
	//sph_gsl->compute(xx/rr, yy/rr, zz/rr, 1.);
	
	(void)rr;
	
      }
    }

    time(&timer_standard_end); 

    std::cout << "GSL approach: " << difftime( timer_standard_end, timer_standard_start) << " seconds" << std::endl;

    
    // using Eigen vectorization

    time_t timer_eigen_start,  timer_eigen_end;
    auto sph_eigen = cbl::glob::SphericalHarmonicsArray::factory(cbl::glob::SpHarMethod::_EIGEN_, 20);

    time(&timer_eigen_start); 

    std::cout << "Starting the Eigen method..." << std::endl;

    for (long i=0; i<nran; ++i) {
      for (long j=i+1; j<nran; ++j) {
	
	const double xx = x[i]-x[j];
	const double yy = y[i]-y[j];
	const double zz = z[i]-z[j];
	double rr = sqrt(xx*xx+yy*yy+zz*zz);
	(void)rr;
	sph_eigen->compute(xx/rr, yy/rr, zz/rr, 1.);
      }
    }
    sph_eigen->print();

    time(&timer_eigen_end); 

    std::cout << "Eigen approach: " << difftime( timer_eigen_end, timer_eigen_start) << " seconds" << std::endl;

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }

  return 0;
}

 
