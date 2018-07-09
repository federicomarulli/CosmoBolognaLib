// ========================================================
// Example code: how to fit a function with a generic model
// ========================================================

#include "Sampler.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
std::string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;


double func (std::vector<double> &x)
{
  return -(100.*pow(x[1]-x[0]*x[0], 2)+pow(1.-x[0], 2))/20.;
}

std::vector<std::vector<double>> starting_values (const int nwalkers, const std::vector<double> center, const double radius)
{
  int npar = center.size();
  std::vector<std::vector<double>> start(nwalkers, std::vector<double>(npar, 0));
  cbl::random::NormalRandomNumbers ran(0, radius, 43213234, -100, 100);
  for(int i=0; i<nwalkers;i++)
    for(int j=0;j<npar;j++)
      start[i][j] = center[j]+ran();
  return start;
}


int main () {
  
  try {

    int npar = 2;
    int nwalkers = 10;
    int chain_size = 1000;

    double radius = 1.e-3;
    std::vector<double> center= {0,0};
    std::vector<std::vector<double>> start = starting_values(nwalkers, center, radius);
  
    std::string dir_output = cbl::par::DirLoc+"../output/";
    cbl::statistics::Sampler sampler(npar, &func);

    // sample the likelihood
    sampler.sample_stretch_move(chain_size, nwalkers, start, 4314234, 2);

    // write the chain output
    sampler.write_chain(dir_output, "chain.dat", 100, 10);

    // sample the likelihood in parallel
    sampler.sample_stretch_move_parallel(chain_size, nwalkers, start, 4314234, 2);
  
    // write the chain output
    sampler.write_chain(dir_output, "chain_parallel.dat", 100, 10);

    // sample the likelihood in parallel
    sampler.sample_stretch_move_parallel(chain_size, nwalkers, start, 434, 2);

    cbl::statistics::Sampler sampler2(npar, &func);

    // sample the likelihood
    sampler2.sample_stretch_move(chain_size, nwalkers, start, 4314234, 2);

    // write the chain output
    sampler2.write_chain(dir_output, "chain2.dat", 100, 10);

    // sample the likelihood in parallel
    sampler2.sample_stretch_move_parallel(chain_size, nwalkers, start, 4314234, 2);
  
    // write the chain output
    sampler2.write_chain(dir_output, "chain_parallel2.dat", 100, 10);


  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }

  return 0;
}

 
