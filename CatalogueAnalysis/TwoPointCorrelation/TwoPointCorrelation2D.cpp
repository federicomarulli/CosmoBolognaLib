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
 *  @file CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation2D.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation2D used to
 *  measure the monopole of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation2D used to measure the monopole of the
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "TwoPointCorrelation2D.h"
#include "Data2D_extra.h"

using namespace cosmobl;
using namespace catalogue;
using namespace chainmesh;
using namespace data;
using namespace pairs;
using namespace twopt;


// ===================================================================================================


void cosmobl::twopt::TwoPointCorrelation2D::write_pairs (const shared_ptr<pairs::Pair> PP, const string dir, const string file) const 
{  
  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {}
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  
  // ----- standard info: scales at the bin centre + number of pairs -----
  
  if (PP->pairInfo()==_standard_)
    for (int i=0; i<PP->nbins_D1(); i++)
      for (int j=0; j<PP->nbins_D2(); j++) 
	fout << setiosflags(ios::fixed) << setprecision(5) << setw(8) << i << "   " << setw(8) << j << "   " << setw(8) << PP->scale_D1(i) << "   " << setw(8) << PP->scale_D2(j) << "   " << setw(8) << PP->PP2D(i, j) << endl;


  // ----- standard + extra info -----
  
  else if (PP->pairInfo()==_extra_)
    for (int i=0; i<PP->nbins_D1(); i++)
      for (int j=0; j<PP->nbins_D2(); j++) 
	fout << setiosflags(ios::fixed) << setprecision(5) << setw(8) << i << "   " << setw(8) << j << "   " << setw(8) << PP->scale_D1(i) << "   " << setw(8) << PP->scale_D2(j) << "   " << setw(8) << PP->PP2D(i, j) << "   " << setw(8) << PP->scale_D1_mean(i, j) << "   " << setw(8) << PP->scale_D1_sigma(i, j) << "   " << setw(8) << PP->scale_D2_mean(i, j) << "   " << setw(8) << PP->scale_D2_mean(i, j) << "   " << setw(8) << PP->z_mean(i, j) << "   " << setw(8) << PP->z_sigma(i, j) << endl;

  else
      ErrorCBL("Error in write_pairs() of TwoPointCorrelation2D.cpp: no such pairInfo!");
  
  fout.clear(); fout.close(); coutCBL << "I wrote the file " << file_out << endl;
  
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation2D::read_pairs (shared_ptr<pairs::Pair> PP, const vector<string> dir, const string file) const
{
  if (dir.size()==0)
    ErrorCBL("Error in cosmobl::twopt::TwoPointCorrelation2D::read_pairs of TwoPointCorrelation2D.cpp! dir.size()=0!");

  int i, j;
  double rad1, rad2, pairs;
  

  // ----- standard info: scales at the bin centre + number of pairs -----
  
  if (PP->pairInfo()==_standard_)
  
    for (size_t dd=0; dd<dir.size(); dd++) {
        
      string file_in = dir[dd]+file; 
      coutCBL << "I'm reading the pair file: " << file_in << endl;
      
      ifstream fin(file_in.c_str()); checkIO(fin, file_in);

      while (fin >> i >> j >> rad1 >> rad2 >> pairs) 
	PP->add_data2D(i, j, {pairs});
      
      fin.clear(); fin.close(); coutCBL << "I read the file " << file_in << endl;
    }

  
  // ----- standard + extra info -----

  else if (PP->pairInfo()==_extra_)
  
    for (size_t dd=0; dd<dir.size(); dd++) {
        
      string file_in = dir[dd]+file; 
      coutCBL << "I'm reading the pair file: " << file_in << endl;
      
      ifstream fin(file_in.c_str()); checkIO(fin, file_in);

      double scale_D1_mean, scale_D1_sigma, scale_D2_mean, scale_D2_sigma, z_mean, z_sigma;
      
      while (fin >> i >> j >> rad1 >> rad2 >> pairs >> scale_D1_mean >> scale_D1_sigma >> scale_D2_mean >> scale_D2_sigma >> z_mean >> z_sigma) 
	PP->add_data2D(i, j, {pairs, scale_D1_mean, scale_D1_sigma, scale_D2_mean, scale_D2_sigma, z_mean, z_sigma});

      fin.clear(); fin.close(); coutCBL << "I read the file " << file_in << endl;
      
    }
  
  else
    ErrorCBL("Error in read_pairs() of TwoPointCorrelation2D.cpp: no such pairInfo!");
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation2D::write_pairs (const vector<shared_ptr<pairs::Pair>> PP, const string dir, const string file) const 
{  
  size_t nRegions = m_data->region_list().size();

  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {}
  
  string file_out = dir+file;  
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  
  // ----- standard info: scales at the bin centre + number of pairs -----
  
  if (PP[0]->pairInfo()==_standard_) 

    for (size_t i=0; i<nRegions; i++) 
      for (size_t j=i; j<nRegions; j++) {
	int index = i*nRegions-(i-1)*i/2+j-i;
	for (int r1=0; r1<PP[index]->nbins_D1(); r1++)
	  for (int r2=0; r2<PP[index]->nbins_D2(); r2++)
	    if (PP[index]->PP2D(r1, r2)>0)
	      fout << setiosflags(ios::fixed) << setprecision(5) << setw(8) << i << "   " << setw(8) << j << "   " << setw(8) << r1 << "   " << setw(8) << r2 << "   " << setw(8) << PP[index]->scale_D1(r1) << setw(8) << "   " << PP[index]->scale_D2(r2) << "   " << setw(8) << PP[index]->PP2D(r1, r2) << endl;
      }

  
  // ----- standard + extra info -----
  
  else if (PP[0]->pairInfo()==_extra_) 

    for (size_t i=0; i<nRegions; i++) 
      for (size_t j=i; j<nRegions; j++) {
	int index = i*nRegions-(i-1)*i/2+j-i;
	for (int r1=0; r1<PP[index]->nbins_D1(); r1++)
	  for (int r2=0; r2<PP[index]->nbins_D2(); r2++)
	    if (PP[index]->PP2D(r1, r2)>0)
	      fout << setiosflags(ios::fixed) << setprecision(5) << setw(8) << i << "   " << setw(8) << j << "   " << setw(8) << r1 << "   " << setw(8) << r2 << "   " << setw(8) << PP[index]->scale_D1(r1) << setw(8) << "   " << PP[index]->scale_D2(r2) << "   " << setw(8) << PP[index]->PP2D(r1, r2) << "   " << setw(8) << PP[index]->scale_D1_mean(r1, r2) << "   " << setw(8) << PP[index]->scale_D1_sigma(r1, r2) << "   " << setw(8) << PP[index]->scale_D2_mean(r1, r2) << "   " << setw(8) << PP[index]->scale_D2_sigma(r1, r2) << "   " << setw(8) << PP[index]->z_mean(r1, r2) << "   " << setw(8) << PP[index]->z_sigma(r1, r2) << endl;
      }

  else
    ErrorCBL("Error in write_pairs() of TwoPointCorrelation2D.cpp: no such pairInfo!");
  
  fout.clear(); fout.close();
  
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation2D::read_pairs (vector<shared_ptr<pairs::Pair>> PP, const vector<string> dir, const string file) const
{
  size_t nRegions = m_data->region_list().size();

  int i, j, r1, r2, index;
  double rad1, rad2, pairs;

    
  // ----- standard info: scales at the bin centre + number of pairs -----

  if (PP[0]->pairInfo()==_standard_)
	
    for (size_t dd=0; dd<dir.size(); dd++) {

      string ff = dir[dd]+file; 
      coutCBL << "I'm reading the pair file: " << ff << endl;
      ifstream fin(ff.c_str()); checkIO(fin, ff);
      
      while (fin >> i >> j >> r1 >> r2 >> rad1 >> rad2 >> pairs) {
	index = i*nRegions-(i-1)*i/2+j-i;
	PP[index]->add_data2D(r1, r2, {pairs});
      }
      
      fin.clear(); fin.close(); coutCBL << "I read the file " << ff << endl;
    }
  

  // ----- standard + extra info -----
  
  else if (PP[0]->pairInfo()==_extra_) 

    for (size_t dd=0; dd<dir.size(); dd++) {

      string ff = dir[dd]+file; 
      coutCBL << "I'm reading the pair file: " << ff << endl;
      ifstream fin(ff.c_str()); checkIO(fin, ff);
    
      double scale_D1_mean, scale_D1_sigma, scale_D2_mean, scale_D2_sigma, z_mean, z_sigma;
      
      while (fin >> i >> j >> r1 >> r2 >> rad1 >> rad2 >> pairs >> scale_D1_mean >> scale_D1_sigma >> scale_D2_mean >> scale_D2_sigma >> z_mean >> z_sigma) {
	index = i*nRegions-(i-1)*i/2+j-i;
	PP[index]->add_data2D(r1, r2, {pairs, scale_D1_mean, scale_D1_sigma, scale_D2_mean, scale_D2_sigma, z_mean, z_sigma});
      }
      
      fin.clear(); fin.close(); coutCBL << "I read the file " << ff << endl;
    }
  
  else
    ErrorCBL("Error in read_pairs() of TwoPointCorrelation2D.cpp: no such pairInfo!");

}


// ============================================================================


shared_ptr<data::Data> cosmobl::twopt::TwoPointCorrelation2D::data_with_extra_info (const shared_ptr<pairs::Pair> dd, const vector<double> scale_D1, const vector<double> scale_D2, const vector<vector<double>> xi, const vector<vector<double>> error) const
{
  vector<vector<double>> extra(6);
  
  for (int i=0; i<dd->nbins_D1(); ++i) 
    for (int j=0; j<dd->nbins_D2(); ++j) {
      extra[0].push_back(dd->scale_D1_mean(i, j));
      extra[1].push_back(dd->scale_D1_sigma(i, j));
      extra[2].push_back(dd->scale_D2_mean(i, j));
      extra[3].push_back(dd->scale_D2_sigma(i, j));
      extra[4].push_back(dd->z_mean(i, j));
      extra[5].push_back(dd->z_sigma(i, j));
    }

  return move(unique_ptr<data::Data2D_extra>(new data::Data2D_extra(scale_D1, scale_D2, xi, error, extra)));
}


// ============================================================================


shared_ptr<Data> cosmobl::twopt::TwoPointCorrelation2D::NaturalEstimator (const shared_ptr<pairs::Pair> dd, const shared_ptr<pairs::Pair> rr, const int nData, const int nRandom)
{
  vector<double> scale_D1, scale_D2;
  vector<vector<double>> xi, error;
  scale_D1.resize(m_dd->nbins_D1()); 
  scale_D2.resize(m_dd->nbins_D2()); 

  xi.resize(m_dd->nbins_D1(), vector<double>(m_dd->nbins_D2(), 0));
  error.resize(m_dd->nbins_D1(), vector<double>(m_dd->nbins_D2(), 0));

  double norm = double(nRandom)*double(nRandom-1)/(double(nData)*double(nData-1));

  for (int i=0; i<dd->nbins_D1(); i++) {
    scale_D1[i] = dd->scale_D1(i);
    for (int j=0; j<dd->nbins_D2(); j++) {
      scale_D2[j] = dd->scale_D2(j);

      xi[i][j] = -1.;
      error[i][j] = 1000.;

      if (dd->PP2D(i, j)>0) {

	if (rr->PP2D(i, j)<1.e-30) {
	  string MSG = "Error in cosmobl::twopt::TwoPointCorrelation2D::NaturalEstimator: there are no random objects in the bin ("+conv(i, par::fINT)+", "+conv(j, par::fINT)+"); please, either increase the total number of random objects or enlarge the bin size!";
	  ErrorCBL(MSG);
	}
	
	xi[i][j] = max(-1., norm*dd->PP2D(i, j)/rr->PP2D(i, j)-1.);

	error[i][j]= PoissonError(dd->PP2D(i, j), rr->PP2D(i, j), 0, nData, nRandom); 
      }
    }
  }

  return (!m_compute_extra_info) ? move(unique_ptr<Data2D>(new Data2D(scale_D1, scale_D2, xi, error))) : data_with_extra_info(dd, scale_D1, scale_D2, xi, error);
}


// ============================================================================


shared_ptr<Data> cosmobl::twopt::TwoPointCorrelation2D::LandySzalayEstimator (const shared_ptr<pairs::Pair> dd, const shared_ptr<pairs::Pair> rr, shared_ptr<pairs::Pair> dr, int nData, int nRandom)
{
  vector<double> scale_D1, scale_D2;
  vector<vector<double>> xi, error;
  scale_D1.resize(m_dd->nbins_D1()); 
  scale_D2.resize(m_dd->nbins_D2()); 

  xi.resize(m_dd->nbins_D1(), vector<double>(m_dd->nbins_D2(), 0));
  error.resize(m_dd->nbins_D1(), vector<double>(m_dd->nbins_D2(), 0));
  
  double norm1 = double(nRandom)*double(nRandom-1)/(double(nData)*double(nData-1));
  double norm2 = double(nRandom-1)/double(nData);

  for (int i=0; i<dd->nbins_D1(); i++) {
    scale_D1[i] = dd->scale_D1(i);
    for (int j=0; j<dd->nbins_D2(); j++) {
      scale_D2[j] = dd->scale_D2(j);

      xi[i][j] = -1.;
      error[i][j] = 1000.;

      if (dd->PP2D(i, j)>0) {

	if (rr->PP2D(i, j)<1.e-30) {
	  string MSG = "Error in cosmobl::twopt::TwoPointCorrelation2D::LandySzalayEstimator: there are no random objects in the bin ("+conv(i, par::fINT)+", "+conv(j, par::fINT)+"); please, either increase the total number of random objects or enlarge the bin size!";
	  ErrorCBL(MSG);
	}
	
	xi[i][j] = max(-1., norm1*dd->PP2D(i, j)/rr->PP2D(i, j)-norm2*dr->PP2D(i, j)/rr->PP2D(i, j)+1.);

	error[i][j]= PoissonError(dd->PP2D(i, j), rr->PP2D(i, j), dr->PP2D(i, j), nData, nRandom); 
      }
    }
  }
 
  return (!m_compute_extra_info) ? move(unique_ptr<Data2D>(new Data2D(scale_D1, scale_D2, xi, error))) : data_with_extra_info(dd, scale_D1, scale_D2, xi, error);
}


// ============================================================================


vector<shared_ptr<Data>> cosmobl::twopt::TwoPointCorrelation2D::XiJackknife (const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr)
{
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<shared_ptr<Data>> data;

  for (size_t i=0; i<nRegions; i++) {

    coutCBL << "analysing region: " << i << " of " << nRegions << "\r"; cout.flush();
    
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2(), m_rr->angularUnits(), m_rr->angularWeight());

    double nData_SS = m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    double nRandom_SS = m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);

    vector<int> w(nRegions, 1);
    w[i] = 0;

    for (size_t j=0; j<nRegions; j++) 
      for (size_t k=j; k<nRegions; k++) {
	int index = j*nRegions-(j-1)*j/2+k-j;
	double ww = w[j]*w[k];
	if (ww>0) 
	  for (int bin1=0; bin1<dd_SS->nbins_D1(); bin1++) 
	    for (int bin2=0; bin2<dd_SS->nbins_D2(); bin2++) {
	      dd_SS->add_data2D(bin1, bin2, dd[index]);
	      rr_SS->add_data2D(bin1, bin2, rr[index]);
	    }
      }
    
    data.push_back(move(NaturalEstimator(dd_SS, rr_SS, nData_SS, nRandom_SS)));
  }
  
  return data;
}


// ============================================================================


vector<shared_ptr<Data>> cosmobl::twopt::TwoPointCorrelation2D::XiJackknife (const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr, const vector<shared_ptr<pairs::Pair>> dr)
{
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<shared_ptr<Data>> data;

  for (size_t i=0; i<nRegions; i++) {

    coutCBL << "analysing region: " << i << " of " << nRegions << "\r"; cout.flush();
    
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2(), m_rr->angularUnits(), m_rr->angularWeight());
    auto dr_SS = Pair::Create(m_dr->pairType(), m_dr->pairInfo(), m_dr->sMin_D1(), m_dr->sMax_D1(), m_dr->nbins_D1(), m_dr->shift_D1(), m_dr->sMin_D2(), m_dr->sMax_D2(), m_dr->nbins_D2(), m_dr->shift_D2(), m_dr->angularUnits(), m_dr->angularWeight());

    double nData_SS = m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    double nRandom_SS = m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);

    vector<int> w(nRegions, 1);
    w[i] = 0;

    for (size_t j=0; j<nRegions; j++) 
      for (size_t k=j; k<nRegions; k++) {
	int index = j*nRegions-(j-1)*j/2+k-j;
	double ww = w[j]*w[k];
	if (ww>0) 
	  for (int bin1=0; bin1<dd_SS->nbins_D1(); bin1++) 
	    for (int bin2=0; bin2<dd_SS->nbins_D2(); bin2++) {
	      dd_SS->add_data2D(bin1, bin2, dd[index]);
	      rr_SS->add_data2D(bin1, bin2, rr[index]);
	      dr_SS->add_data2D(bin1, bin2, dr[index]);
	    }
      }
      
    data.push_back(move(LandySzalayEstimator(dd_SS, rr_SS, dr_SS, nData_SS, nRandom_SS)));  
  }
  
  return data;
}


// ============================================================================


vector<shared_ptr<Data>> cosmobl::twopt::TwoPointCorrelation2D::XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr)
{
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<shared_ptr<Data>> data;
  vector<double> nData_reg, nRandom_reg;

  for (int i=0; i<nMocks; i++) {
    nData_reg.push_back(m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg.push_back(m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
  }

  uniform_int_distribution<int> uni(0, nRegions-1);
  default_random_engine rng;
  int val = 3; // see Norberg et al. 2009

  for (int i=0; i<nMocks; i++) {

    coutCBL << "analysing mock: " << i << " of " << nMocks << "\r"; cout.flush();
    
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2(), m_rr->angularUnits(), m_rr->angularWeight());
    
    double nData_SS = 0., nRandom_SS = 0.;

    vector<int> w(nRegions, 0);
    for (size_t n=0; n<val*nRegions; n++)
      w[uni(rng)] +=1;

    for (size_t j=0; j<nRegions; j++) {
      nData_SS += w[j]*nData_reg[j];
      nRandom_SS += w[j]*nRandom_reg[j];

      for (size_t k=j; k<nRegions; k++) {
	int index = j*nRegions-(j-1)*j/2+k-j;
	double ww = (k==j) ? w[k] : w[j]*w[k];
	if (ww>0) 
	  for (int bin1=0; bin1<dd_SS->nbins_D1(); bin1++) 
	    for (int bin2=0; bin2<dd_SS->nbins_D2(); bin2++) {
	      dd_SS->add_data2D(bin1, bin2, dd[index], ww);
	      rr_SS->add_data2D(bin1, bin2, rr[index], ww);
	    }
      }
    }
    
    data.push_back(move(NaturalEstimator(dd_SS, rr_SS, nData_SS, nRandom_SS)));
  }
  
  return data;
}


// ============================================================================


vector<shared_ptr<Data>> cosmobl::twopt::TwoPointCorrelation2D::XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr, const vector<shared_ptr<pairs::Pair>> dr)
{
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<shared_ptr<Data>> data;
  vector<double> nData_reg, nRandom_reg;

  for (int i=0; i<nMocks; i++) {
    nData_reg.push_back(m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg.push_back(m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
  }

  uniform_int_distribution<int> uni(0, nRegions-1);
  default_random_engine rng;
  int val = 3; // see Norberg et al. 2009

  for (int i=0; i<nMocks; i++) {

    coutCBL << "analysing mock: " << i << " of " << nMocks << "\r"; cout.flush();
    
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2(), m_rr->angularUnits(), m_rr->angularWeight());
    auto dr_SS = Pair::Create(m_dr->pairType(), m_dr->pairInfo(), m_dr->sMin_D1(), m_dr->sMax_D1(), m_dr->nbins_D1(), m_dr->shift_D1(), m_dr->sMin_D2(), m_dr->sMax_D2(), m_dr->nbins_D2(), m_dr->shift_D2(), m_dr->angularUnits(), m_dr->angularWeight());

    double nData_SS = 0., nRandom_SS = 0.;

    vector<int> w(nRegions, 0);
    for (size_t n=0; n<val*nRegions; n++)
      w[uni(rng)] +=1;

    for (size_t j=0; j<nRegions; j++) {
      nData_SS += w[j]*nData_reg[j];
      nRandom_SS += w[j]*nRandom_reg[j];

      for (size_t k=j; k<nRegions; k++) {
	int index = j*nRegions-(j-1)*j/2+k-j;
	double ww = (k==j) ? w[k] : w[j]*w[k];
	if (ww>0) 
	  for (int bin1=0;bin1<dd_SS->nbins_D1();bin1++) 
	    for (int bin2=0;bin2<dd_SS->nbins_D2();bin2++) {
	      dd_SS->add_data2D(bin1, bin2, dd[index], ww);
	      rr_SS->add_data2D(bin1, bin2, rr[index], ww);
	      dr_SS->add_data2D(bin1, bin2, dr[index], ww);
	    }
      }
    }

    data.push_back(move(LandySzalayEstimator(dd_SS, rr_SS, dr_SS, nData_SS, nRandom_SS)));
  }
  
  return data;
}


