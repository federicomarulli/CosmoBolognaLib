// ========================================================
// Example code: how to fit a function with a generic model
// ========================================================

#include "Sampler.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;


double func (vector<double> &x)
{
  return -(100.*pow(x[1]-x[0]*x[0], 2)+pow(1.-x[0], 2))/20.;
}

vector<vector<double>> starting_values(const int nwalkers, const vector<double> center, const double radius)
{
  int npar = center.size();
  vector<vector<double>> start(nwalkers, vector<double>(npar, 0));
  cosmobl::random::NormalRandomNumbers ran(0, radius, 43213234, -100, 100);
  for(int i=0; i<nwalkers;i++)
    for(int j=0;j<npar;j++)
      start[i][j] = center[j]+ran();
  return start;
}


int main () {
  
  try {

    int npar = 2;
    int nwalkers = 1000;
    int chain_size = 10000;

    double radius = 1.e-3;
    vector<double> center= {0,0};
    vector<vector<double>> start = starting_values(nwalkers, center, radius);
  
    string dir_output = "../output/";
    cosmobl::statistics::Sampler sampler(npar, &func);

    // sample the likelihood
    sampler.sample_stretch_move(chain_size, nwalkers, start, 4314234, 2);

    // write the chain output
    sampler.write_chain(dir_output, "chain.dat", 100, 10);

    // sample the likelihood in parallel
    sampler.sample_stretch_move_parallel(chain_size, nwalkers, start, 4314234, 2);
  
    // write the chain output
    sampler.write_chain(dir_output, "chain_parallel.dat", 100, 10);

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }

  return 0;
}

 
