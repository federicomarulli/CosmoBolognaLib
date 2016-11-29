/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo  *
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

/**
 *  @file CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation1D.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation1D used to
 *  measure the monopole of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation1D used to measure the monopole of the
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation1D.h"
#include "Data1D_extra.h"

using namespace cosmobl;
using namespace catalogue;
using namespace chainmesh;
using namespace pairs;
using namespace twopt;


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::write_pairs (const shared_ptr<pairs::Pair> PP, const string dir, const string file) const 
{  
  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {}
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  
  // ----- standard info: scales at the bin centre + number of pairs -----

  if (PP->pairInfo()==_standard_)
    for (int i=0; i<PP->nbins(); i++) 
      fout << setiosflags(ios::fixed) << setprecision(5) << setw(8) << i << "   " << setw(8) << PP->scale(i) << "   " << setw(8) << PP->PP1D(i) << endl;

  
  // ----- standard + extra info -----
  
  else if (PP->pairInfo()==_extra_) 
    for (int i=0; i<PP->nbins(); i++)
      fout << setiosflags(ios::fixed) << setprecision(5) << setw(8) << i << "   " << setw(8) << PP->scale(i) << "   " << setw(8) << PP->PP1D(i) << "   " << setw(8) << PP->scale_mean(i) << "   " << setw(8) << PP->scale_sigma(i) << "   " << setw(8) << PP->z_mean(i) << "   " << setw(8) << PP->z_sigma(i) << endl;

  else
      ErrorCBL("Error in write_pairs() of TwoPointCorrelation1D.cpp: no such pairInfo!");
  
  fout.clear(); fout.close(); coutCBL << "I wrote the file " << file_out << endl;
  
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::read_pairs (shared_ptr<pairs::Pair> PP, const vector<string> dir, const string file) const
{
  if (dir.size()==0)
    ErrorCBL("Error in cosmobl::twopt::TwoPointCorrelation1D::read_pairs of TwoPointCorrelation1D.cpp! dir.size()=0!");

  int i;
  double rad, pairs;
  
  
  // ----- standard info: scales at the bin centre + number of pairs -----
  
  if (PP->pairInfo()==_standard_)
  
    for (size_t dd=0; dd<dir.size(); dd++) {
      
      string file_in = dir[dd]+file; 
      ifstream fin(file_in.c_str()); checkIO(fin, file_in);
      
      while (fin >> i >> rad >> pairs) 
	PP->add_data1D(i, {pairs});

      fin.clear(); fin.close(); coutCBL << "I read the file " << file_in << endl;
    }
  
    
  // ----- standard + extra info -----

  else if (PP->pairInfo()==_extra_)

    for (size_t dd=0; dd<dir.size(); dd++) {
      
      string file_in = dir[dd]+file; 
      ifstream fin(file_in.c_str()); checkIO(fin, file_in);
      
      double scale_mean, scale_sigma, z_mean, z_sigma;
      
      while (fin >> i >> rad >> pairs >> scale_mean >> scale_sigma >> z_mean >> z_sigma) 
	PP->add_data1D(i, {pairs, scale_mean, scale_sigma, z_mean, z_sigma});

      fin.clear(); fin.close(); coutCBL << "I read the file " << file_in << endl;
    }
    
  else
    ErrorCBL("Error in read_pairs() of TwoPointCorrelation1D.cpp: no such pairInfo!");
    
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::write_pairs (const vector<shared_ptr<pairs::Pair>> PP, const string dir, const string file) const 
{  
  size_t nRegions = m_data->region_list().size();
  
  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {}
  
  string file_out = dir+file;  
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

   
  // ----- standard info: scales at the bin centre + number of pairs -----

  if (PP[0]->pairInfo()==_standard_)
    
    for (size_t i=0; i<nRegions; i++) 
      for (size_t j=i; j<nRegions; j++) {
	int index = i*nRegions+j-(i-1)*i/2-i;
	for (int r1=0; r1<PP[index]->nbins(); r1++)
	  if (PP[index]->PP1D(r1)>0)
	    fout << setiosflags(ios::fixed) << setprecision(5) << setw(8) << i << "   " << setw(8) << j << "   " << setw(8) << r1 << "   " << setw(8) << PP[index]->scale(r1) << "   " << setw(8) << PP[index]->PP1D(r1) << endl;
      }

  
  // ----- standard + extra info -----

  else if (PP[0]->pairInfo()==_extra_)
    
    for (size_t i=0; i<nRegions; i++) 
      for (size_t j=i; j<nRegions; j++) {
	int index = i*nRegions+j-(i-1)*i/2-i;
	for (int r1=0; r1<PP[index]->nbins(); r1++)
	  if (PP[index]->PP1D(r1)>0)
	    fout << setiosflags(ios::fixed) << setprecision(5) << setw(8) << i << "   " << setw(8) << j << "   " << setw(8) << r1 << "   " << setw(8) << PP[index]->scale(r1) << "   " << setw(8) << PP[index]->PP1D(r1) << "   " << setw(8) << PP[index]->scale_mean(r1) << "   " << setw(8) << PP[index]->scale_sigma(r1) << "   " << setw(8) << PP[index]->z_mean(r1) << "   " << setw(8) << PP[index]->z_sigma(r1) << endl;
      }
      
  else
    ErrorCBL("Error in write_pairs() of TwoPointCorrelation1D.cpp: no such pairInfo!");
      
  fout.clear(); fout.close();
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::read_pairs (vector<shared_ptr<pairs::Pair>> PP, const vector<string> dir, const string file) const
{
  size_t nRegions = m_data->region_list().size();

  int i, j, bin, index;
  double rad, pairs;
    
   
  // ----- standard info: scales at the bin centre + number of pairs -----

  if (PP[0]->pairInfo()==_standard_)
    
    for (size_t dd=0; dd<dir.size(); dd++) {
      string ff = dir[dd]+file; 
      ifstream fin(ff.c_str()); checkIO(fin, ff);
      
      while (fin >> i >> j >> bin >> rad >> pairs) {
	index = i*nRegions-(i-1)*i/2+j-i;
	PP[index]->add_data1D(bin, {pairs});
      }
      
      fin.clear(); fin.close(); coutCBL << "I read the file " << ff << endl;
    }
  

  // ----- standard + extra info -----

  else if (PP[0]->pairInfo()==_extra_)

    for (size_t dd=0; dd<dir.size(); dd++) {
      string ff = dir[dd]+file; 
      ifstream fin(ff.c_str()); checkIO(fin, ff);
      
      double scale_mean, scale_sigma, z_mean, z_sigma;

      while (fin >> i >> j >> bin >> rad >> pairs >> scale_mean >> scale_sigma >> z_mean >> z_sigma) {
	index = i*nRegions-(i-1)*i/2+j-i;
	PP[index]->add_data1D(bin, {pairs, scale_mean, scale_sigma, z_mean, z_sigma});
      }
      
      fin.clear(); fin.close(); coutCBL << "I read the file " << ff << endl;
    }
  
  else
    ErrorCBL("Error in read_pairs() of TwoPointCorrelation1D.cpp: no such pairInfo!");
}


// ============================================================================


shared_ptr<data::Data> cosmobl::twopt::TwoPointCorrelation1D::data_with_extra_info (const shared_ptr<pairs::Pair> dd, const vector<double> rad, const vector<double> xi, const vector<double> error) const
{
  vector<vector<double>> extra(4);
  
  for (int i=0; i<dd->nbins(); ++i) {
    extra[0].push_back(dd->scale_mean(i));
    extra[1].push_back(dd->scale_sigma(i));
    extra[2].push_back(dd->z_mean(i));
    extra[3].push_back(dd->z_sigma(i));
  }
  
  return move(unique_ptr<data::Data1D_extra>(new data::Data1D_extra(rad, xi, error, extra)));
}


// ============================================================================


shared_ptr<data::Data> cosmobl::twopt::TwoPointCorrelation1D::NaturalEstimator (const shared_ptr<pairs::Pair> dd, const shared_ptr<pairs::Pair> rr, const int nData, const int nRandom)
{
  vector<double> rad(m_dd->nbins()), xi(m_dd->nbins(), -1.), error(m_dd->nbins(), 1000.);
  
  double norm = double(nRandom)*double(nRandom-1)/(double(nData)*double(nData-1));

  for (int i=0; i<dd->nbins(); i++) {

    rad[i] = dd->scale(i);
    
    if (dd->PP1D(i)>0) {

      if (rr->PP1D(i)<1.e-30) 
	ErrorCBL("Error in NaturalEstimator() of TwoPointCorrelation1D.cpp: there are no random objects in the bin "+conv(i, par::fINT)+"; please, either increase the total number of random objects or enlarge the bin size!");
      
      xi[i] = max(-1., norm*dd->PP1D(i)/rr->PP1D(i)-1.);
      error[i] = PoissonError(dd->PP1D(i), rr->PP1D(i), 0, nData, nRandom); 
    }
    
  }

  return (!m_compute_extra_info) ? move(unique_ptr<data::Data1D>(new data::Data1D(rad, xi, error))) : data_with_extra_info(dd, rad, xi, error);
}


// ============================================================================


shared_ptr<data::Data> cosmobl::twopt::TwoPointCorrelation1D::LandySzalayEstimator (const shared_ptr<pairs::Pair> dd, const shared_ptr<pairs::Pair> rr, const shared_ptr<pairs::Pair> dr, const int nData, const int nRandom)
{
  vector<double> rad(m_dd->nbins()), xi(m_dd->nbins(), -1.), error(m_dd->nbins(), 1000.);
  
  double norm1 = double(nRandom)*double(nRandom-1)/(double(nData)*double(nData-1));
  double norm2 = double(nRandom-1)/double(nData);
  
  for (int i=0; i<dd->nbins(); i++) {
  
    rad[i] = dd->scale(i);
    
    if (dd->PP1D(i)>0) {

      if (rr->PP1D(i)<1.e-30) 
	ErrorCBL("Error in LandySzalayEstimator() of TwoPointCorrelation1D.cpp: there are no random objects in the bin "+conv(i, par::fINT)+"; please, either increase the total number of random objects or enlarge the bin size!");

      xi[i] = max(-1., norm1*dd->PP1D(i)/rr->PP1D(i)-norm2*dr->PP1D(i)/rr->PP1D(i)+1.);
      error[i] = PoissonError(dd->PP1D(i), rr->PP1D(i), dr->PP1D(i), nData, nRandom);
      
    }
  }
  
  return (!m_compute_extra_info) ? move(unique_ptr<data::Data1D>(new data::Data1D(rad, xi, error))) : data_with_extra_info(dd, rad, xi, error);
}


// ============================================================================


vector<shared_ptr<data::Data>> cosmobl::twopt::TwoPointCorrelation1D::XiJackknife (const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr)
{
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<shared_ptr<data::Data>> data;

  for (size_t i=0; i<nRegions; i++) {

    coutCBL << "analysing region: " << i << " of " << nRegions << "\r"; cout.flush();
    
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin(), m_dd->sMax(), m_dd->nbins(), m_dd->shift(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin(), m_rr->sMax(), m_rr->nbins(), m_rr->shift(), m_rr->angularUnits(), m_rr->angularWeight());
    
    double nData_SS = m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    double nRandom_SS = m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    
    vector<int> w(nRegions, 1);
    w[i] = 0;
    
    for (size_t j=0; j<nRegions; j++) 
      for (size_t k=j; k<nRegions; k++) {
	int index = j*nRegions-(j-1)*j/2+k-j;
	double ww = w[j]*w[k];
	if (ww>0)
	  for (int bin=0; bin<dd_SS->nbins(); bin++) {
	    dd_SS->add_data1D(bin, dd[index]);
	    rr_SS->add_data1D(bin, rr[index]);
	  }
      }

    // finalise the computation of the extra information
    if (m_compute_extra_info) {
      dd_SS->finalise();
      rr_SS->finalise();
    }
    
    data.push_back(move(NaturalEstimator(dd_SS, rr_SS, nData_SS, nRandom_SS)));
  } 
  
  return data;
}


// ============================================================================


vector<shared_ptr<data::Data>> cosmobl::twopt::TwoPointCorrelation1D::XiJackknife (const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr, const vector<shared_ptr<pairs::Pair>> dr)
{
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<shared_ptr<data::Data>> data;

  
  for (size_t i=0; i<nRegions; i++) {
    
    coutCBL << "analysing region: " << i << " of " << nRegions << "\r"; cout.flush();
      
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin(), m_dd->sMax(), m_dd->nbins(), m_dd->shift(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin(), m_rr->sMax(), m_rr->nbins(), m_rr->shift(), m_rr->angularUnits(), m_rr->angularWeight());
    auto dr_SS = Pair::Create(m_dr->pairType(), m_dr->pairInfo(), m_dr->sMin(), m_dr->sMax(), m_dr->nbins(), m_dr->shift(), m_dr->angularUnits(), m_dr->angularWeight());

    double nData_SS = m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    double nRandom_SS = m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);

    vector<int> w(nRegions, 1);
    w[i] = 0;

    for (size_t j=0; j<nRegions; j++) 
      for (size_t k=j; k<nRegions; k++) {
	int index = j*nRegions+k-(j-1)*j/2-j;
	double ww = w[j]*w[k];
	if (ww>0) 
	  for (int bin=0; bin<dd_SS->nbins(); bin++) {
	    dd_SS->add_data1D(bin, dd[index]);
	    rr_SS->add_data1D(bin, rr[index]);
	    dr_SS->add_data1D(bin, dr[index]);
	  }
      }

    // finalise the computation of the extra information
    if (m_compute_extra_info) {
      dd_SS->finalise();
      rr_SS->finalise();
      dr_SS->finalise();
    }
    
    data.push_back(move(LandySzalayEstimator(dd_SS, rr_SS, dr_SS, nData_SS, nRandom_SS)));	
  }

  return data;
}


// ============================================================================


vector<shared_ptr<data::Data>> cosmobl::twopt::TwoPointCorrelation1D::XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr)
{
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<shared_ptr<data::Data>> data;
  vector<double> nData_reg, nRandom_reg;

  for (int i=0; i<nMocks; i++) {
    nData_reg.push_back(m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg.push_back(m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
  }

  uniform_int_distribution<int> uni(0, nRegions-1);
  default_random_engine rng;
  int val = 1; // see Norberg et al. 2009

    
  for (int i=0; i<nMocks; i++) {

    coutCBL << "analysing mock: " << i << " of " << nMocks << "\r"; cout.flush();
    
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin(), m_dd->sMax(), m_dd->nbins(), m_dd->shift(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin(), m_rr->sMax(), m_rr->nbins(), m_rr->shift(), m_rr->angularUnits(), m_rr->angularWeight());

    double nData_SS=0, nRandom_SS=0;

    vector<int> w(nRegions, 0);
    for (size_t n=0; n<val*nRegions; n++)
      w[uni(rng)] ++;

    for (size_t j=0; j<nRegions; j++) {
      nData_SS += w[j]*nData_reg[j];
      nRandom_SS += w[j]*nRandom_reg[j];

      for (size_t k=j; k<nRegions; k++) {
	int index = j*nRegions+k-(j-1)*j/2-j;
	double ww = (k==j) ? w[k] : w[j]*w[k];
	if (ww>0) 
	  for (int bin=0; bin<dd_SS->nbins(); bin++) {
	    dd_SS->add_data1D(bin, dd[index], ww);
	    rr_SS->add_data1D(bin, rr[index], ww);
	  }
      }
    }

    // finalise the computation of the extra information
    if (m_compute_extra_info) {
      dd_SS->finalise();
      rr_SS->finalise();
    }
      
    data.push_back(move(NaturalEstimator(dd_SS, rr_SS, nData_SS, nRandom_SS)));
  }
  
  return data;
}


// ============================================================================


vector<shared_ptr<data::Data>> cosmobl::twopt::TwoPointCorrelation1D::XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr, const vector<shared_ptr<pairs::Pair>> dr)
{
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<shared_ptr<data::Data>> data;
  vector<double> nData_reg, nRandom_reg;

  for (int i=0; i<nMocks; i++) {
    nData_reg.push_back(m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg.push_back(m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
  }

  uniform_int_distribution<int> uni(0, nRegions-1);
  default_random_engine rng;
  int val = 1; // see Norberg et al. 2009

  
  for (int i=0; i<nMocks; i++) {

    coutCBL << "analysing mock: " << i << " of " << nMocks << "\r"; cout.flush();
    
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin(), m_dd->sMax(), m_dd->nbins(), m_dd->shift(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin(), m_rr->sMax(), m_rr->nbins(), m_rr->shift(), m_rr->angularUnits(), m_rr->angularWeight());
    auto dr_SS = Pair::Create(m_dr->pairType(), m_dr->pairInfo(), m_dr->sMin(), m_dr->sMax(), m_dr->nbins(), m_dr->shift(), m_dr->angularUnits(), m_dr->angularWeight());

    double nData_SS=0, nRandom_SS=0;

    vector<int> w(nRegions, 0);
    for (size_t n=0; n<val*nRegions; n++)
      w[uni(rng)] ++;

    
    for (size_t j=0; j<nRegions; j++) {
      nData_SS += w[j]*nData_reg[j];
      nRandom_SS += w[j]*nRandom_reg[j];
      
      for (size_t k=j; k<nRegions; k++) {
	int index = j*nRegions+k-(j-1)*j/2-j;
	double ww = (k==j) ? w[k] : w[j]*w[k];
	if (ww>0) 
	  for (int bin=0; bin<dd_SS->nbins(); bin++) {
	    dd_SS->add_data1D(bin, dd[index], ww);
	    rr_SS->add_data1D(bin, rr[index], ww);
	    dr_SS->add_data1D(bin, dr[index], ww);
	  }
      }
    }

    // finalise the computation of the extra information
    if (m_compute_extra_info) {
      dd_SS->finalise();
      rr_SS->finalise();
      dr_SS->finalise();
    }
    
    data.push_back(move(LandySzalayEstimator(dd_SS, rr_SS, dr_SS, nData_SS, nRandom_SS)));
  }

  return data;
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::read_covariance_matrix (const string dir, const string file)
{
  m_dataset->set_covariance(dir+file);
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::write_covariance_matrix (const string dir, const string file) const
{
  m_dataset->write_covariance(dir, file, "r");
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::compute_covariance_matrix (const vector<shared_ptr<data::Data>> xi_collection, const bool doJK)
{
  vector<vector<double>> xi;

  for (size_t i=0; i<xi_collection.size(); i++)
    xi.push_back(xi_collection[i]->fx());

  vector<vector<double>> cov_mat;
  cosmobl::covariance_matrix(xi, cov_mat, doJK);
  
  m_dataset->set_covariance(cov_mat);
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::compute_covariance_matrix (const vector<string> file_xi, const bool doJK)
{
  vector<double> rad, mean;
  vector<vector<double>> cov_mat;

  cosmobl::covariance_matrix(file_xi, rad, mean, cov_mat, doJK); 
  m_dataset->set_covariance(cov_mat);
}

