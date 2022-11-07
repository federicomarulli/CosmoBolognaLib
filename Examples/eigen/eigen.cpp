// ===================================================================
// Test the differences in performances using Eigen std::vectorization
// ===================================================================

#include "RandomNumbers.h"

int main () {

  try {
   
    // create random points in a cubic volume
    
    const double BoxSize = 500.;
    const long nran = 200;
    const int seed = 666;
    
    cbl::random::UniformRandomNumbers random(0, BoxSize, seed);

    std::vector<double> x(nran), y(nran), z(nran);

    std::vector<cbl::Vector3D> v(nran);

    for (long i=0; i<nran; ++i) {
      x[i] = random();
      y[i] = random();
      z[i] = random();
      v[i] = {x[i], y[i], z[i]};
    }

    
    // test the standard way to compute the norm of a std::vector

    time_t timer_standard_start, timer_standard_end;

    time(&timer_standard_start); 

    std::cout << "Starting the standard method..." << std::endl;

    for (long i=0; i<nran; ++i) {
      for (long j=0; j<nran; ++j) {

	const double xx = x[i]-x[j];
	const double yy = y[i]-y[j];
	const double zz = z[i]-z[j];
	double rr = sqrt(xx*xx+yy*yy+zz*zz);
	(void)rr;
	
      }
    }

    time(&timer_standard_end); 

    std::cout << "Standard approach: " << difftime( timer_standard_end, timer_standard_start) << " seconds" << std::endl;

    
    // test the Eigen std::vectorization

    time_t timer_eigen_start,  timer_eigen_end;

    time(&timer_eigen_start); 

    std::cout << "Starting the Eigen method..." << std::endl;

    for (long i=0; i<nran; ++i) {
      for (long j=0; j<nran; ++j) {
	
	double rr = (v[i]-v[j]).norm();
	(void)rr;

      }
    }

    time(&timer_eigen_end); 

    std::cout << "Eigen approach: " << difftime( timer_eigen_end, timer_eigen_start) << " seconds" << std::endl;

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }

  return 0;
}

 
