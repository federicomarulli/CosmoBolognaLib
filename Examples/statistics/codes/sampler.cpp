// ========================================================
// Example code: how to fit a function with a generic model
// ========================================================

#include "Sampler.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;


double func (const vector<double> x)
{
  return -(100.*pow(x[1]-x[0]*x[0], 2)+pow(1.-x[0], 2))/20.;
}


double func2 (const vector<double> x)
{
  double ff1 = 0.5*cosmobl::gaussian(x[0], NULL, {-4, 1, 1})*cosmobl::gaussian(x[1], NULL, {30, 0.5, 1});
  double ff2 = 0.5*cosmobl::gaussian(x[0], NULL, {4, 1, 1})*cosmobl::gaussian(x[1], NULL, {30, 0.5, 1});
  return log(ff1+ff2+0.1*exp(func(x)));
}


int main () {

  int npar = 2;
  int nchains = 1000;
  int chain_size = 10000;
  vector<double> start = {0, 0};
  
  string dir_output = "../output/";
  cosmobl::statistics::Sampler sampler(npar, &func2);

  // sample the likelihood
  sampler.sample_stretch_move(nchains, chain_size, start, 1., 4314234, 2);

  // write the chain output
  sampler.write_chain(dir_output, "chain.dat", 100, nchains, 10);

  // sample the likelihood in parallel
  sampler.sample_stretch_move_parallel(nchains, chain_size, start, 1.e-4, 4314234, 2);
  
  // write the chain output
  sampler.write_chain(dir_output, "chain_parallel.dat", 100, nchains, 10);

  return 0;
}

 
