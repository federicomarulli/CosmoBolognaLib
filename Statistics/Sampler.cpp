/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
 *  federico.marulli3@unibo.it                                      *
 *                                                                  *
 *  This program is free software; you can redistribute it and/or   *
 *  modify it under the terms of the GNU General Public License as  *
 *  published by the Free Software Foundation; either version 2 of  *
 *  the License, or (at your option) any later version.             *
 *                                                                  *
 *  This program is distributed in the hope that it will be useful, *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   *
 *  GNU General Public License for more details.                    *
 *                                                                  *
 *  You should have received a copy of the GNU General Public       *
 *  License along with this program; if not, write to the Free      *
 *  Software Foundation, Inc.,                                      *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.       *
 ********************************************************************/

/** @file Statistics/Sampler.cpp
 *
 *  @brief Methods of the class Sampler 
 *
 *  This file contains the implementation of the methods of the class
 *  Sampler, used for bayesian analyses
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Sampler.h"

using namespace std;

using namespace cbl;


// ============================================================================================


shared_ptr<random::DistributionRandomNumbers> cbl::statistics::Sampler::m_set_gz (const int seed, const double aa)
{
  if (aa<=1)
    ErrorCBL("Error in cbl::statistics::Sampler::m_set_gz() of Sampler.cpp: the stretch move parameter must be >1!");

  double zmin = 1./aa;
  double zmax = aa;
  vector<double> zz = linear_bin_vector(1000,zmin,zmax), gzz;
  for (auto &&zzz : zz)
    gzz.push_back(1./sqrt(zzz));

  return make_shared<random::DistributionRandomNumbers>(random::DistributionRandomNumbers(zz, gzz, "Spline", seed));
}


// ============================================================================================



void cbl::statistics::Sampler::m_initialize_chains(const vector< vector<double> > start)
{
  for (int i=0; i<m_nwalkers; i++) { 
    vector<double> pp = start[i];
    m_function_chain[0][i] = m_function(pp);
    m_chains[0][i] = pp;
  }
}


// ============================================================================================


void cbl::statistics::Sampler::set_chain(const int npar, const int npar_free, const int chain_size, const int nwalkers)
{
  m_npar = npar;
  m_npar_free = npar_free;

  m_nwalkers = nwalkers;
  m_chain_size = chain_size;

  m_chains.erase(m_chains.begin(), m_chains.end());
  m_chains.resize(m_chain_size, vector<vector<double>>(m_nwalkers, vector<double>(m_npar,0)));

  m_acceptance.erase(m_acceptance.begin(), m_acceptance.end());
  m_acceptance.resize(m_nwalkers, 0.);

  m_function_chain.erase(m_function_chain.begin(), m_function_chain.end());
  m_function_chain.resize(m_chain_size, vector<double>(m_nwalkers, 0.));
}


// ============================================================================================


void cbl::statistics::Sampler::set_function (const function<double(vector<double> &)> function)
{
  m_function = function;
  m_use_python = false;
}


// ============================================================================================

/*
void cbl::statistics::Sampler::set_function (PyObject *pp)
{
  PyCallback pc(pp);
  m_function = bind(&PyCallback::V2D, pc, placeholders::_1);
  m_use_python = true;
}
*/

// ============================================================================================


void cbl::statistics::Sampler::sample_stretch_move (const int chain_size, const int nwalkers, const vector<vector<double>> start, const int seed, const double aa, const string outputFile)
{
  set_chain(m_npar, m_npar_free, chain_size, nwalkers);

  function<void(int, int, const vector<double>, const double)>  write_func;
  ofstream fout;

  write_func = [&] (const int nc, const int nw, const vector<double> parameters, const double func) 	
  {(void)nc; (void)nw; (void)parameters; (void)func;};

  if (outputFile!=par::defaultString) {
    fout.open(outputFile.c_str());
    write_func = [&] (const int nc, const int nw, const vector<double> parameters, const double func)
      {
	fout << nc << " " << nw << " ";
	for (int p=0; p<m_npar; p++)
	  fout << parameters[p] << " ";
	fout << func << endl;
      };
  }

  // set seed for random_numbers
  random::UniformRandomNumbers_Int seeds(0, 23412432, seed);

  const int chain_seed = seeds();
  const int MH_seed = seeds();
  const int gz_seed = seeds();

  coutCBL << "Starting seed = " << seed << endl;
  coutCBL << "Seed for random walker choice = " << chain_seed << endl;
  coutCBL << "Seed for Metropolis-Hastings algorithm = " << MH_seed << endl;
  coutCBL << "Seed for random extraction from g(z) = " << gz_seed << endl;
 
  random::UniformRandomNumbers_Int chains(0, m_nwalkers-1, chain_seed);
  random::UniformRandomNumbers MH_random(0., 1., MH_seed);
  auto gz = m_set_gz(gz_seed, aa);

  // initialize chains
  m_initialize_chains(start);
  
  for (int n=1; n<m_chain_size; n++) {
    for (int i=0; i<m_nwalkers; i++)
    {
      int kk = i;
      while (kk==i)
	kk = chains();
        
      vector<double> parameters_i;
      vector<double> parameters_k;

      parameters_i = m_chains[n-1][i];
      parameters_k = m_chains[n-1][kk];

      vector<double> parameters(m_npar, 0);
      double proposed_function_chains = par::defaultDouble;

      double gen_z = gz->operator()();

      gen_z = gz->operator()();
      for (int p=0; p<m_npar; p++)
	parameters[p] = parameters_k[p] + gen_z*(parameters_i[p]-parameters_k[p]);
      proposed_function_chains = m_function(parameters);

      parameters_i = parameters;

      double ratio = min(1.,pow(gen_z, m_npar_free-1)*exp(proposed_function_chains-m_function_chain[n-1][i]));

      if (MH_random() <ratio) {
	m_function_chain[n][i] = proposed_function_chains;
	m_chains[n][i] = parameters_i;
	m_acceptance[i] += 1./m_chain_size;
      }
      else {
	m_function_chain[n][i] = m_function_chain[n-1][i];
	m_chains[n][i] = m_chains[n-1][i]; 
      }

      write_func(n, i, m_chains[n][i], m_function_chain[n][i]);
    }
    
    double progress = double((n+1)*m_nwalkers)/(m_nwalkers*m_chain_size)*100;
    coutCBL << setprecision(2) << setiosflags(ios::fixed) << setw(8) << progress << "% \r"; cout.flush();
  }
  
  cout << endl;
  coutCBL << "Done!" << endl;
}


// ============================================================================================


void cbl::statistics::Sampler::m_sample_stretch_move_parallel_cpp (const int chain_size, const int nwalkers, const vector<vector<double>> start, const int seed, const double aa)
{
  set_chain(m_npar, m_npar_free, chain_size, nwalkers);

  // set seeds and random_numbers
  random::UniformRandomNumbers_Int seeds(0, 23412432, seed);

  vector<double> chains_index = linear_bin_vector(m_nwalkers, 0., m_nwalkers-1.);
  vector<double> chains_weights(m_nwalkers,1);

  const int chain_seed1 = seeds();
  const int chain_seed2 = seeds();
  const int MH_seed = seeds();
  const int gz_seed = seeds();

  coutCBL << "Starting seed = " << seed << endl;
  coutCBL << "Seed for random walker choice, first half = " << chain_seed1 << endl;
  coutCBL << "Seed for random walker choice, second half = " << chain_seed2 << endl;
  coutCBL << "Seed for Metropolis-Hastings algorithm = " << MH_seed << endl;
  coutCBL << "Seed for random extraction from g(z) = " << gz_seed << endl;

  int half = m_nwalkers/2;
  random::UniformRandomNumbers_Int chains_first_half(0, half-1, chain_seed1);
  random::UniformRandomNumbers_Int chains_second_half(half, m_nwalkers-1, chain_seed2);

  vector<random::UniformRandomNumbers_Int> chains_sample = {chains_second_half, chains_first_half};

  random::UniformRandomNumbers MH_random(0., 1., MH_seed);

  auto gz = m_set_gz(gz_seed, aa);

  // initialise the chains
  m_initialize_chains(start);
  
  // stretch-move

  for (int n=1; n<m_chain_size; n++) {
    for (int ss=0; ss<2; ss++) {
#pragma omp parallel for schedule(dynamic)
      for (int ii=0; ii<half; ii++)
      {
	int nn = n-1+ss; 
	int i = half*ss+ii;
	int kk = int(chains_sample[ss]());;

	vector<double> parameters_i= m_chains[n-1][i];
	vector<double> parameters_k = m_chains[nn][kk];;

        vector<double> parameters(m_npar, 0);
        double proposed_function_chains = par::defaultDouble;

        double gen_z = gz->operator()();

	gen_z = gz->operator()();
	for (int p=0; p<m_npar; p++)
	  parameters[p] = parameters_k[p] + gen_z*(parameters_i[p]-parameters_k[p]);
	proposed_function_chains = m_function(parameters);

        parameters_i = parameters;

	double ratio = min(1.,pow(gen_z, m_npar_free-1)*exp(proposed_function_chains-m_function_chain[n-1][i]));

        if (MH_random() <ratio) {
          m_function_chain[n][i] = proposed_function_chains;
	  m_chains[n][i] = parameters_i;
	  m_acceptance[i] += 1./m_chain_size;
	}
	else {
	  m_function_chain[n][i] = m_function_chain[n-1][i];
	  m_chains[n][i] = m_chains[n-1][i]; 
	}

      }
    }
    
    double progress = double((n+1)*m_nwalkers)/(m_nwalkers*m_chain_size)*100;
    coutCBL << setprecision(2) << setiosflags(ios::fixed) << setw(8) << progress << "% \r"; cout.flush();
  }

  cout << endl;
  coutCBL << "Done!" << endl;
}


// ============================================================================================


void cbl::statistics::Sampler::m_sample_stretch_move_parallel_py (const int chain_size, const int nwalkers, const vector<vector<double>> start, const int seed, const double aa)
{
  (void)chain_size; (void)nwalkers; (void)start; (void)seed; (void)aa;
  cbl::ErrorCBL("Work in progress", glob::ExitCode::_workInProgress_);


  /*
  set_chain(m_npar, m_npar_free, chain_size, nwalkers);

  // set seeds and random_numbers
  random::UniformRandomNumbers seeds(0, 23412432, seed);

  vector<double> chains_index = linear_bin_vector(m_nwalkers, 0., m_nwalkers-1.);
  vector<double> chains_weights(m_nwalkers,1);

  int half = m_nwalkers/2;
  random::DiscreteRandomNumbers chains_first_half(chains_index, chains_weights, int(seeds()), 0, half-1);
  random::DiscreteRandomNumbers chains_second_half(chains_index, chains_weights, int(seeds()), half, m_nwalkers-1);

  vector<random::DiscreteRandomNumbers> chains_sample = {chains_second_half, chains_first_half};

  random::UniformRandomNumbers MH_random(0., 1., int(seeds()));

  auto gz = m_set_gz(int(seeds()), aa);

  // initialize chains
  m_initialize_chains(start);

  // stretch-move

  Py_BEGIN_ALLOW_THREADS

  for (int n=1; n<m_chain_size; n++) {
    for (int ss=0; ss<2; ss++) {
#pragma omp parallel for schedule(dynamic)
      for (int ii=0; ii<half; ii++)
      {
	PyGILState_STATE gstate = PyGILState_Ensure();

	int nn = n-1+ss; 
	int i = half*ss+ii;
	int kk = int(chains_sample[ss]());;

	vector<double> parameters_i= m_chains[n-1][i];
	vector<double> parameters_k = m_chains[nn][kk];;

        double gen_z = gz->operator()();

        for (int p=0; p<m_npar; p++)
          parameters_i[p] = parameters_k[p] + gen_z*(parameters_i[p]-parameters_k[p]);

        double proposed_function_chains = m_function(parameters_i);

        double ratio = min(1.,pow(gen_z, m_npar-1)*exp(proposed_function_chains-m_function_chain[n-1][i]));

        if (MH_random() <ratio) {
          m_function_chain[n][i] = proposed_function_chains;
          m_chains[n][i] = parameters_i;
          m_acceptance[i] += 1./m_chain_size;
        }
        else{
          m_function_chain[n][i] = m_function_chain[n-1][i];
          m_chains[n][i] = m_chains[n-1][i]; 
        }

        PyGILState_Release(gstate);

      }
    }

    double progress = double((n+1)*m_nwalkers)/(m_nwalkers*m_chain_size)*100;
    coutCBL << setprecision(2) << setiosflags(ios::fixed) << setw(8) << progress << "% \r"; cout.flush();
  }

  cout << endl;

  Py_END_ALLOW_THREADS
    */
}


// ============================================================================================


void cbl::statistics::Sampler::sample_stretch_move_parallel (const int chain_size, const int nwalkers, const vector<vector<double>> start, const int seed, const double aa)
{
  if (nwalkers%2 != 0)
    ErrorCBL("Error in cbl::statistics::Sampler::sample_stretch_move_parallel(): the number of walkers must be an even integer!");

  if (m_use_python)
    m_sample_stretch_move_parallel_py(chain_size, nwalkers, start, seed, aa);
  else
    m_sample_stretch_move_parallel_cpp(chain_size, nwalkers, start, seed, aa);
}


// ============================================================================================


void cbl::statistics::Sampler::get_chain_function_acceptance (vector<vector<double>> &chains, vector<double> &function, vector<double> &acceptance, const int start, const int thin)
{
  chains.resize(m_npar);
  function.erase(function.begin(), function.end());
  acceptance.erase(acceptance.begin(), acceptance.end());

  for (int i=start; i<m_chain_size; i+=thin) {
    for (int j=0; j<m_nwalkers; j++) {
      function.push_back(m_function_chain[i][j]);
      acceptance.push_back(m_acceptance[j]);
      for (int k=0; k<m_npar; k++)
	chains[k].push_back(m_chains[i][j][k]);
    }
  } 
}

  
// ============================================================================================


void cbl::statistics::Sampler::write_chain (const string dir_output, const string file, const int start, const int thin)
{
  string file_out = dir_output+file;
  ofstream fout(file_out.c_str());
  fout.precision(10);

  for (int i=start;i<m_chain_size; i+=thin) {
    for (int j=0; j<m_nwalkers; j++) {
      fout << i*m_nwalkers+j << "  ";
      for (int k=0; k<m_npar; k++) 
	fout <<  m_chains[i][j][k] << "  ";   
      fout << m_function_chain[i][j] << "  " << m_acceptance[j] << endl;
    }
  } 
  fout.clear(); fout.close();

}
