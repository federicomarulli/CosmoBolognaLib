// ==============================================================
// Test the differences in performances using Eigen vectorization
// ==============================================================

#include "Func.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;

int main () {

  try {

    // create random points in a cubic volume
    
    const double BoxSize = 500.;
    
    const long nran = 40000;
    const int seed = 666;
    
    cosmobl::random::UniformRandomNumbers random(0, BoxSize, seed);

    vector<double> x(nran), y(nran), z(nran);

    vector<cosmobl::Vector3D> v(nran);

    for (long i=0; i<nran; i++) {
      const double _x = random();
      const double _y = random();
      const double _z = random();

      x[i] = _x;
      y[i] = _y;
      z[i] = _z;
      v[i] << _x, _y, _z;
    }

    
    // test the standard way to compute the norm of a vector

    time_t timer_standard_start, timer_standard_end;

    time(&timer_standard_start); 

    cout << "Starting the standard method..." << endl;

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

    cout << "Standard approach: " << difftime( timer_standard_end, timer_standard_start) << " seconds" << endl;

    
    // test the Eigen vectorization

    time_t timer_eigen_start,  timer_eigen_end;

    time(&timer_eigen_start); 

    cout << "Starting the Eigen method..." << endl;

    for (long i=0; i<nran; ++i) {
      for (long j=0; j<nran; ++j) {
	
	double rr = (v[j]-v[i]).norm();
	(void)rr;

      }
    }

    time(&timer_eigen_end); 

    cout << "Eigen approach: " << difftime( timer_eigen_end, timer_eigen_start) << " seconds" << endl;

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }

  return 0;
}

 
