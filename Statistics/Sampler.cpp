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

using namespace cosmobl;


// ============================================================================================


shared_ptr<random::DistributionRandomNumbers> cosmobl::statistics::Sampler::m_set_gz(const int seed, const double aa)
{
  if(aa<=1)
    ErrorCBL("Error in Sampler, the stretch move parameter must be >1!");

  double zmin = 1./aa;
  double zmax = aa;
  vector<double> zz = linear_bin_vector(1000,zmin,zmax), gzz;
  for (auto &&zzz : zz)
    gzz.push_back(1./sqrt(zzz));

  return make_shared<random::DistributionRandomNumbers>(random::DistributionRandomNumbers(zz, gzz, "Spline", seed));
}


// ============================================================================================



void cosmobl::statistics::Sampler::m_initialize_chains(const int seed, const vector<double> start, const double radius)
{
  random::NormalRandomNumbers normal_random(0., 1., seed, -100, 100);
  for (int i=0; i<m_nchains; i++) 
    for(int j=0; j<m_npar; j++)
    m_chains[0][i][j] = start[j]+radius*normal_random();
  
  for (int i=0; i<m_nchains; i++) 
    m_function_chain[0][i] = m_function(m_chains[0][i]);
}


// ============================================================================================


void cosmobl::statistics::Sampler::set_chain(const int npar, const int nchains, const int chain_size)
{
  m_npar = npar;
  m_nchains = nchains;
  m_chain_size = chain_size;

  m_chains.resize(m_chain_size, vector<vector<double>>(m_nchains, vector<double>(m_npar,0)));

  m_acceptance.resize(m_nchains, 0.);

  m_function_chain.resize(m_chain_size, vector<double>(m_nchains, 0.));
}


// ============================================================================================


void cosmobl::statistics::Sampler::set_function (const function<double(vector<double>)> function)
{
  m_function = function;
}


// ============================================================================================


void cosmobl::statistics::Sampler::sample_stretch_move (const int nchains, const int chain_size, const vector<double> start, double radius, const int seed, const double aa)
{
  coutCBL << "Sampling..." << endl;

  set_chain(m_npar, nchains, chain_size);

  // set seed for random_numbers
  random::UniformRandomNumbers prior_seeds(0, 23412432, seed);
 
  vector<double> chains_index = linear_bin_vector(m_nchains, 0., m_nchains-1.);
  vector<double> chains_weights(m_nchains,1);
  random::DiscreteRandomNumbers chains(chains_index, chains_weights, int(prior_seeds()), 0, m_nchains-1);

  random::UniformRandomNumbers MH_random(0., 1., int(prior_seeds()));

  auto gz = m_set_gz(int(prior_seeds()), aa);

  // initialize chains
  m_initialize_chains(int(prior_seeds()), start, radius);
  
  for (int n=1; n<m_chain_size; n++) {
    for (int i=0; i<m_nchains; i++)
    {
      int kk = i;
      while (kk==i)
	kk = int(chains());
        
      vector<double> parameters_i;
      vector<double> parameters_k;

      double gen_z = gz->operator()();

      parameters_i = m_chains[n-1][i];
      parameters_k = m_chains[n-1][kk];

      for(int p=0;p<m_npar;p++)
	parameters_i[p] = parameters_k[p] + gen_z*(parameters_i[p]-parameters_k[p]);

      double proposed_function_chains = m_function(parameters_i);

      double ratio = min(1.,pow(gen_z, m_npar-1)*exp(proposed_function_chains-m_function_chain[n-1][i]));

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
    double progress = double((n+1)*m_nchains)/(m_nchains*m_chain_size)*100;
    coutCBL << setprecision(2) << setiosflags(ios::fixed) << setw(8) << progress << "% \r"; cout.flush();

  }
  cout << endl;
}


// ============================================================================================


void cosmobl::statistics::Sampler::sample_stretch_move_parallel (const int nchains, const int chain_size, const vector<double> start, double radius, const int seed, const double aa)
{
  coutCBL << "Sampling..." << endl;

  set_chain(m_npar, nchains, chain_size);

  // set seeds and random_numbers
  random::UniformRandomNumbers prior_seeds(0, 23412432, seed);

  vector<double> chains_index = linear_bin_vector(m_nchains, 0., m_nchains-1.);
  vector<double> chains_weights(m_nchains,1);

  int half = m_nchains/2;
  random::DiscreteRandomNumbers chains_first_half(chains_index, chains_weights, int(prior_seeds()), 0, half-1);
  random::DiscreteRandomNumbers chains_second_half(chains_index, chains_weights, int(prior_seeds()), half, m_nchains-1);

  vector<random::DiscreteRandomNumbers> chains_sample = {chains_second_half, chains_first_half};

  random::UniformRandomNumbers MH_random(0., 1., int(prior_seeds()));

  auto gz = m_set_gz(int(prior_seeds()), aa);

  // initialize chains
  m_initialize_chains(int(prior_seeds()), start, radius);

  // stretch-move

  for (int n=1; n<m_chain_size; n++) {
    for(int ss=0; ss<2; ss++){
#pragma omp parallel for schedule(dynamic)
      for (int ii=0; ii<half; ii++)
      {
	int nn = n-1+ss; 
	int i = half*ss+ii;
	int kk = int(chains_sample[ss]());;

	vector<double> parameters_i= m_chains[n-1][i];
	vector<double> parameters_k = m_chains[nn][kk];;

	double gen_z = gz->operator()();


	for(int p=0;p<m_npar;p++)
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

      }
    }

    double progress = double((n+1)*m_nchains)/(m_nchains*m_chain_size)*100;
    coutCBL << setprecision(2) << setiosflags(ios::fixed) << setw(8) << progress << "% \r"; cout.flush();
  }

  cout << endl;
}


// ============================================================================================


void cosmobl::statistics::Sampler::write_chain(const string dir_output, const string file, const int start, const int stop, const int thin)
{
  string file_out = dir_output+file;
  ofstream fout(file_out.c_str());
  fout.precision(10);

  for(int i=start;i<stop; i+=thin){
    for(int j=0;j<m_nchains;j++){
      fout << i*m_nchains+j << " ";
      for(int k=0;k<m_npar;k++){
	fout <<  m_chains[i][j][k] << "  ";
      }
      fout << m_function_chain[i][j] << " " << m_acceptance[j] << endl;
    }
  } 
  fout.clear(); fout.close();

}
