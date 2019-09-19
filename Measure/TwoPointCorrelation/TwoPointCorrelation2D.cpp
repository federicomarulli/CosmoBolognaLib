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
 *  @file
 *  CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation2D.cpp
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

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace chainmesh;
using namespace data;
using namespace pairs;
using namespace measure;
using namespace twopt;


// ===================================================================================================


void cbl::measure::twopt::TwoPointCorrelation2D::write_pairs (const std::shared_ptr<pairs::Pair> PP, const std::string dir, const std::string file) const 
{  
  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {}
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  
  // ----- standard info: scales at the bin centre + number of pairs -----
  
  if (PP->pairInfo()==PairInfo::_standard_)
    for (int i=0; i<PP->nbins_D1(); i++)
      for (int j=0; j<PP->nbins_D2(); j++) 
	fout << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << i
	     << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << j
	     << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->scale_D1(i)
	     << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->scale_D2(j)
	     << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->PP2D(i, j)
	     << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->PP2D_weighted(i, j) << endl;
  

  // ----- standard + extra info -----
  
  else if (PP->pairInfo()==PairInfo::_extra_)
    for (int i=0; i<PP->nbins_D1(); i++)
      for (int j=0; j<PP->nbins_D2(); j++) 
	fout << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << i
	     << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << j
	     << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->scale_D1(i)
	     << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->scale_D2(j)
	     << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->PP2D(i, j)
	     << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->PP2D_weighted(i, j)
	     << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->scale_D1_mean(i, j)
	     << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->scale_D1_sigma(i, j)
	     << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->scale_D2_mean(i, j)
	     << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->scale_D2_sigma(i, j)
	     << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->z_mean(i, j)
	     << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->z_sigma(i, j) << endl;

  else
    ErrorCBL("no such pairInfo!", "write_pairs", "TwoPointCorrelation2D.cpp");
  
  fout.clear(); fout.close(); coutCBL << "I wrote the file " << file_out << endl;
  
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation2D::read_pairs (std::shared_ptr<pairs::Pair> PP, const std::vector<std::string> dir, const std::string file) const
{
  if (dir.size()==0)
    ErrorCBL("dir.size()=0!", "read_pairs", "TwoPointCorrelation2D.cpp");

  string line;
  int i, j;
  double pairs, weighted_pairs;
  

  // ----- standard info: scales at the bin centre + number of pairs -----
  
  if (PP->pairInfo()==PairInfo::_standard_)
  
    for (size_t dd=0; dd<dir.size(); dd++) {
        
      string file_in = dir[dd]+file; 
      ifstream fin(file_in.c_str()); checkIO(fin, file_in);

      while (getline(fin, line)) {
	
	stringstream ss(line);
	vector<double> num;
	double val;
	while (ss >> val) num.emplace_back(val);

	if (num.size()!=6)
	  ErrorCBL("the number of lines in the input pair file: "+file_in+" must be 6!", "read_pairs", "TwoPointCorrelation2D.cpp");

	i = int(num[0]);
	j = int(num[1]);
	pairs = num[4];
	weighted_pairs = num[5];
	
	PP->add_data2D(i, j, {pairs, weighted_pairs});

      }
      
      fin.clear(); fin.close(); coutCBL << "I read the file " << file_in << endl;
    }

  
  // ----- standard + extra info -----

  else if (PP->pairInfo()==PairInfo::_extra_)
  
    for (size_t dd=0; dd<dir.size(); dd++) {
        
      string file_in = dir[dd]+file; 
      ifstream fin(file_in.c_str()); checkIO(fin, file_in);

      double scale_D1_mean, scale_D1_sigma, scale_D2_mean, scale_D2_sigma, z_mean, z_sigma;

      while (getline(fin, line)) {
	
	stringstream ss(line);
	vector<double> num;
	double val;
	while (ss >> val) num.emplace_back(val);

	if (num.size()!=12)
	  ErrorCBL("the number of lines in the input pair file: "+file_in+" must be 12!", "read_pairs", "TwoPointCorrelation2D.cpp");

	i = int(num[0]);
	j = int(num[1]);
	pairs = num[4];
	weighted_pairs = num[5];
	scale_D1_mean = num[6];
	scale_D1_sigma = num[7];
	scale_D2_mean = num[8];
	scale_D2_sigma = num[9];
	z_mean = num[10];
	z_sigma = num[11];
	
	PP->add_data2D(i, j, {pairs, weighted_pairs, scale_D1_mean, pow(scale_D1_sigma, 2)*weighted_pairs, scale_D2_mean, pow(scale_D2_sigma, 2)*weighted_pairs, z_mean, pow(z_sigma, 2)*weighted_pairs});

      }
      
      fin.clear(); fin.close(); coutCBL << "I read the file " << file_in << endl;
      
    }
  
  else
    ErrorCBL("no such pairInfo!", "read_pairs", "TwoPointCorrelation2D.cpp");
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation2D::write_pairs (const std::vector<std::shared_ptr<pairs::Pair>> PP, const std::string dir, const std::string file) const 
{  
  size_t nRegions = m_data->region_list().size();

  bool cross = (PP.size()==nRegions*nRegions) ? true : false;

  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {}
  
  string file_out = dir+file;  
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  
  // ----- standard info: scales at the bin centre + number of pairs -----
  
  if (PP[0]->pairInfo()==PairInfo::_standard_) 

    for (size_t i=0; i<nRegions; i++) 
      for (size_t j=(cross) ? 0 : i; j<nRegions; j++) {
	int index = (cross) ? i*nRegions+j : i*nRegions+j-(i-1)*i/2-i;
	for (int r1=0; r1<PP[index]->nbins_D1(); r1++)
	  for (int r2=0; r2<PP[index]->nbins_D2(); r2++)
	    if (PP[index]->PP2D(r1, r2)>0)
	      fout << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << i
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << j
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << r1
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << r2
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->scale_D1(r1)
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->scale_D2(r2)
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->PP2D(r1, r2)
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->PP2D_weighted(r1, r2) << endl;
      }

  
  // ----- standard + extra info -----
  
  else if (PP[0]->pairInfo()==PairInfo::_extra_) 

    for (size_t i=0; i<nRegions; i++) 
      for (size_t j=(cross) ? 0 : i; j<nRegions; j++) {
	int index = (cross) ? i*nRegions+j : i*nRegions+j-(i-1)*i/2-i;
	for (int r1=0; r1<PP[index]->nbins_D1(); r1++)
	  for (int r2=0; r2<PP[index]->nbins_D2(); r2++)
	    if (PP[index]->PP2D(r1, r2)>0)
	      fout << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << i
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << j
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << r1
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << r2
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->scale_D1(r1) 
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->scale_D2(r2)
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->PP2D(r1, r2)
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->PP2D_weighted(r1, r2)
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->scale_D1_mean(r1, r2)
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->scale_D1_sigma(r1, r2)
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->scale_D2_mean(r1, r2)
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->scale_D2_sigma(r1, r2)
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->z_mean(r1, r2)
		   << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->z_sigma(r1, r2) << endl;
      }

  else
    ErrorCBL("no such pairInfo!", "write_pairs", "TwoPointCorrelation2D.cpp");
  
  fout.clear(); fout.close();
  
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation2D::read_pairs (std::vector<std::shared_ptr<pairs::Pair>> PP, const std::vector<std::string> dir, const std::string file) const
{
  size_t nRegions = m_data->region_list().size();
  
  bool cross = (PP.size() == nRegions*nRegions) ? true : false;

  int i, j, r1, r2, index;
  double rad1, rad2, pairs, weighted_pairs;

    
  // ----- standard info: scales at the bin centre + number of pairs -----

  if (PP[0]->pairInfo()==PairInfo::_standard_)
	
    for (size_t dd=0; dd<dir.size(); dd++) {

      string ff = dir[dd]+file; 
      coutCBL << "I'm reading the pair file: " << ff << endl;
      ifstream fin(ff.c_str()); checkIO(fin, ff);
      
      while (fin >> i >> j >> r1 >> r2 >> rad1 >> rad2 >> pairs >> weighted_pairs) {
	index = (cross) ? i*nRegions+j : i*nRegions-(i-1)*i/2+j-i;
	PP[index]->add_data2D(r1, r2, {pairs, weighted_pairs});
      }
      
      fin.clear(); fin.close(); coutCBL << "I read the file " << ff << endl;
    }
  

  // ----- standard + extra info -----
  
  else if (PP[0]->pairInfo()==PairInfo::_extra_) 

    for (size_t dd=0; dd<dir.size(); dd++) {

      string ff = dir[dd]+file; 
      coutCBL << "I'm reading the pair file: " << ff << endl;
      ifstream fin(ff.c_str()); checkIO(fin, ff);
    
      double scale_D1_mean, scale_D1_sigma, scale_D2_mean, scale_D2_sigma, z_mean, z_sigma;
      
      while (fin >> i >> j >> r1 >> r2 >> rad1 >> rad2 >> pairs >> weighted_pairs >> scale_D1_mean >> scale_D1_sigma >> scale_D2_mean >> scale_D2_sigma >> z_mean >> z_sigma) {
	index = (cross) ? i*nRegions+j : i*nRegions-(i-1)*i/2+j-i;
	PP[index]->add_data2D(r1, r2, {pairs, weighted_pairs, scale_D1_mean, pow(scale_D1_sigma, 2)*weighted_pairs, scale_D2_mean, pow(scale_D2_sigma, 2)*weighted_pairs, z_mean, pow(z_sigma, 2)*weighted_pairs});
      }

      fin.clear(); fin.close(); coutCBL << "I read the file " << ff << endl;
    }
  
  else
    ErrorCBL("no such pairInfo!", "read_pairs", "TwoPointCorrelation2D.cpp");

}


// ============================================================================


std::shared_ptr<data::Data> cbl::measure::twopt::TwoPointCorrelation2D::data_with_extra_info (const std::shared_ptr<pairs::Pair> dd, const std::vector<double> scale_D1, const std::vector<double> scale_D2, const std::vector<std::vector<double>> xi, const std::vector<std::vector<double>> error) const
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


std::shared_ptr<Data> cbl::measure::twopt::TwoPointCorrelation2D::correlation_NaturalEstimator (const std::shared_ptr<pairs::Pair> dd, const std::shared_ptr<pairs::Pair> rr, const int nData, const double nData_weighted, const int nRandom, const double nRandom_weighted)
{
  vector<double> scale_D1, scale_D2;
  vector<vector<double>> xi, error;
  scale_D1.resize(m_dd->nbins_D1()); 
  scale_D2.resize(m_dd->nbins_D2()); 

  xi.resize(m_dd->nbins_D1(), vector<double>(m_dd->nbins_D2(), 0));
  error.resize(m_dd->nbins_D1(), vector<double>(m_dd->nbins_D2(), 0));

  // number of objects in the data catalogue
  int nD = (nData>0) ? nData : m_data->nObjects();

  // weighted number of objects in the data catalogue
  double nDw = (nData_weighted>0) ? nData_weighted : m_data->weightedN();

  // number of objects in the random catalogue
  int nR = (nRandom>0) ? nRandom : m_random->nObjects();

  // weighted number of objects in the random catalogue
  double nRw = (nRandom_weighted>0) ? nRandom_weighted : m_random->weightedN();

  // inverse of the total number of data-data pairs
  double nDDi = 1./(nDw*(nDw-1.)*0.5);

  // inverse of the total number of random-random pairs
  double nRRi = 1./(nRw*(nRw-1.)*0.5);
  
  for (int i=0; i<dd->nbins_D1(); i++) {
    scale_D1[i] = dd->scale_D1(i);
    for (int j=0; j<dd->nbins_D2(); j++) {
      scale_D2[j] = dd->scale_D2(j);

      xi[i][j] = -1.;
      error[i][j] = 1000.;

      if (dd->PP2D_weighted(i, j)>0) {

	if (rr->PP2D_weighted(i, j)<1.e-30) 
	  ErrorCBL("there are no random objects in the bin "+conv(i, par::fINT)+","+conv(j, par::fINT)+"; please, either increase the total number of random objects or enlarge the bin size! (dd="+conv(dd->PP2D_weighted(i, j), par::fDP3)+", rr="+conv(rr->PP2D_weighted(i, j), par::fDP3)+")!", "correlation_NaturalEstimator", "TwoPointCorrelation2D.cpp");

	// normalised number of data-data weighted pairs
	double DD_norm = dd->PP2D_weighted(i, j)*nDDi;
	
	// normalised number of random-random weighted pairs
	double RR_norm = rr->PP2D_weighted(i, j)*nRRi;

	// natural estimator
	xi[i][j] = max(-1., DD_norm/RR_norm-1.);
	
	// Poisson error
	error[i][j]= PoissonError(Estimator::_natural_, dd->PP2D(i, j), rr->PP2D(i, j), 0, nD, nR); 
      }
    }
  }

  return (!m_compute_extra_info) ? move(unique_ptr<Data2D>(new Data2D(scale_D1, scale_D2, xi, error))) : data_with_extra_info(dd, scale_D1, scale_D2, xi, error);
}


// ============================================================================


std::shared_ptr<Data> cbl::measure::twopt::TwoPointCorrelation2D::correlation_LandySzalayEstimator (const std::shared_ptr<pairs::Pair> dd, const std::shared_ptr<pairs::Pair> rr, const std::shared_ptr<pairs::Pair> dr, const int nData, const double nData_weighted, const int nRandom, const double nRandom_weighted)
{
  vector<double> scale_D1, scale_D2;
  vector<vector<double>> xi, error;
  scale_D1.resize(m_dd->nbins_D1()); 
  scale_D2.resize(m_dd->nbins_D2()); 

  xi.resize(m_dd->nbins_D1(), vector<double>(m_dd->nbins_D2(), 0));
  error.resize(m_dd->nbins_D1(), vector<double>(m_dd->nbins_D2(), 0));
  
  // number of objects in the data catalogue
  int nD = (nData>0) ? nData : m_data->nObjects();

  // weighted number of objects in the data catalogue
  double nDw = (nData_weighted>0) ? nData_weighted : m_data->weightedN();

  // number of objects in the random catalogue
  int nR = (nRandom>0) ? nRandom : m_random->nObjects();

  // weighted number of objects in the random catalogue
  double nRw = (nRandom_weighted>0) ? nRandom_weighted : m_random->weightedN();

  // inverse of the total number of data-data pairs
  double nDDi = 1./(nDw*(nDw-1.)*0.5);

  // inverse of the total number of random-random pairs
  double nRRi = 1./(nRw*m_random_dilution_fraction*(nRw*m_random_dilution_fraction-1.)*0.5); 

  // inverse of the total number of data-random pairs
  double nDRi = 1./(nDw*nRw);

  for (int i=0; i<dd->nbins_D1(); i++) {
    scale_D1[i] = dd->scale_D1(i);
    for (int j=0; j<dd->nbins_D2(); j++) {
      scale_D2[j] = dd->scale_D2(j);
      
      xi[i][j] = -1.;
      error[i][j] = 1000.;

      if (dd->PP2D_weighted(i, j)>0) {

	if (rr->PP2D_weighted(i, j)<1.e-30) 
	  ErrorCBL("there are no random objects in the bin "+conv(i, par::fINT)+","+conv(j, par::fINT)+"; please, either increase the total number of random objects or enlarge the bin size! (dd="+conv(dd->PP2D_weighted(i, j), par::fDP3)+", rr="+conv(rr->PP2D_weighted(i, j), par::fDP3)+")!", "correlation_LandySzalayEstimator", "TwoPointCorrelation2D.cpp");
		   
	// normalised number of data-data weighted pairs
	double DD_norm = dd->PP2D_weighted(i, j)*nDDi;
	
	// normalised number of random-random weighted pairs
	double RR_norm = rr->PP2D_weighted(i, j)*nRRi;
	
	// normalised number of data-random weighted pairs
	double DR_norm = dr->PP2D_weighted(i, j)*nDRi;

	// Landy & Szalay estimator
	xi[i][j] = max(-1., (DD_norm-2.*DR_norm)/RR_norm+1.);

	// Poisson error
	error[i][j]= PoissonError(Estimator::_LandySzalay_, dd->PP2D(i, j), rr->PP2D(i, j), dr->PP2D(i, j), nD, nR);
	
      }
    }
  }
  
  return (!m_compute_extra_info) ? move(unique_ptr<Data2D>(new Data2D(scale_D1, scale_D2, xi, error))) : data_with_extra_info(dd, scale_D1, scale_D2, xi, error);
}


// ============================================================================


std::vector<std::shared_ptr<Data>> cbl::measure::twopt::TwoPointCorrelation2D::XiJackknife (const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr)
{
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<shared_ptr<Data>> data;

  for (size_t i=0; i<nRegions; i++) {

    coutCBL << "analysing region: " << i << " of " << nRegions << "\r"; cout.flush();
    
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2(), m_rr->angularUnits(), m_rr->angularWeight());

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
    

    int nData_SS = m_data->nObjects_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    double nData_SS_weighted = m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    int nRandom_SS = m_random->nObjects_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    double nRandom_SS_weighted = m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    
    data.push_back(move(correlation_NaturalEstimator(dd_SS, rr_SS, nData_SS, nData_SS_weighted, nRandom_SS, nRandom_SS_weighted)));
  }
  
  return data;
}


// ============================================================================


std::vector<std::shared_ptr<Data>> cbl::measure::twopt::TwoPointCorrelation2D::XiJackknife (const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr, const std::vector<std::shared_ptr<pairs::Pair>> dr)
{
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<shared_ptr<Data>> data;

  for (size_t i=0; i<nRegions; i++) {

    coutCBL << "analysing region: " << i << " of " << nRegions << "\r"; cout.flush();
    
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2(), m_rr->angularUnits(), m_rr->angularWeight());
    auto dr_SS = Pair::Create(m_dr->pairType(), m_dr->pairInfo(), m_dr->sMin_D1(), m_dr->sMax_D1(), m_dr->nbins_D1(), m_dr->shift_D1(), m_dr->sMin_D2(), m_dr->sMax_D2(), m_dr->nbins_D2(), m_dr->shift_D2(), m_dr->angularUnits(), m_dr->angularWeight());
    
    vector<int> w(nRegions, 1);
    w[i] = 0;

    for (size_t j=0; j<nRegions; j++){

      if(w[j]>0){
        for (size_t k=j; k<nRegions; k++) { // auto pairs
          if(w[k]>0){
            int index = j*nRegions+k-(j-1)*j/2-j;
            for (int bin1=0; bin1<dd_SS->nbins_D1(); bin1++) 
              for (int bin2=0; bin2<dd_SS->nbins_D2(); bin2++) {
                dd_SS->add_data2D(bin1, bin2, dd[index]);
                rr_SS->add_data2D(bin1, bin2, rr[index]);
              }
          }
        }

        for (size_t k=0; k<nRegions; k++) {// cross pairs
          if (w[k]>0){
            int index = j*nRegions+k;
            for (int bin1=0; bin1<dd_SS->nbins_D1(); bin1++) 
              for (int bin2=0; bin2<dd_SS->nbins_D2(); bin2++) {
                dr_SS->add_data2D(bin1, bin2, dr[index]);
              }
          }
        }
      }
    }

    double nData_SS = m_data->nObjects_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    double nData_SS_weighted = m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    double nRandom_SS = m_random->nObjects_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    double nRandom_SS_weighted = m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);

    data.push_back(move(correlation_LandySzalayEstimator(dd_SS, rr_SS, dr_SS, nData_SS, nData_SS_weighted, nRandom_SS, nRandom_SS_weighted)));  
  }
  
  return data;
}


// ============================================================================


std::vector<std::shared_ptr<Data>> cbl::measure::twopt::TwoPointCorrelation2D::XiBootstrap (const int nMocks, const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr, const int seed)
{
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<shared_ptr<Data>> data;
  vector<double> nData_reg, nData_reg_weighted, nRandom_reg, nRandom_reg_weighted;

  for (size_t i=0; i<nRegions; i++) {
    nData_reg.push_back(m_data->nObjects_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nData_reg_weighted.push_back(m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg.push_back(m_random->nObjects_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg_weighted.push_back(m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
  }

  random::UniformRandomNumbers_Int ran(0., nRegions-1, seed);

  int val = 3; // see Norberg et al. 2009

  int nbins_D1 = m_dd->nbins_D1();
  int nbins_D2 = m_dd->nbins_D2();

  for (int i=0; i<nMocks; i++) {

    coutCBL << "analysing mock: " << i << " of " << nMocks << "\r"; cout.flush();

    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2(), m_rr->angularUnits(), m_rr->angularWeight());

    double nData_SS = 0., nData_SS_weighted = 0., nRandom_SS = 0., nRandom_SS_weighted = 0.;

    vector<int> w(nRegions, 0);
    for (size_t n=0; n<val*nRegions; n++)
      w[ran()] ++;

    for (size_t j=0; j<nRegions; j++) {

      if (w[j]>0) {

	nData_SS += w[j]*nData_reg[j];
	nData_SS_weighted += w[j]*nData_reg_weighted[j];
	nRandom_SS += w[j]*nRandom_reg[j];
	nRandom_SS_weighted += w[j]*nRandom_reg_weighted[j];

	for (size_t k=j; k<nRegions; k++) {
	  if(w[k]>0){
	    int index = j*nRegions-(j-1)*j/2+k-j;
	    double ww = w[j]*w[k];
	    for (int bin1=0; bin1<nbins_D1; bin1++) 
	      for (int bin2=0; bin2<nbins_D2; bin2++) {
		dd_SS->add_data2D(bin1, bin2, dd[index], ww);
		rr_SS->add_data2D(bin1, bin2, rr[index], ww);
	      }
	  }
	}
      }
    }

    data.push_back(move(correlation_NaturalEstimator(dd_SS, rr_SS, nData_SS, nData_SS_weighted, nRandom_SS, nRandom_SS_weighted)));
  }


  return data;
}


// ============================================================================


std::vector<std::shared_ptr<Data>> cbl::measure::twopt::TwoPointCorrelation2D::XiBootstrap (const int nMocks, const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr, const std::vector<std::shared_ptr<pairs::Pair>> dr, const int seed)
{
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<shared_ptr<Data>> data;
  vector<double> nData_reg, nData_reg_weighted, nRandom_reg, nRandom_reg_weighted;

  for (size_t i=0; i<nRegions; i++) {
    nData_reg.push_back(m_data->nObjects_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nData_reg_weighted.push_back(m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg.push_back(m_random->nObjects_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg_weighted.push_back(m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
  }

  random::UniformRandomNumbers_Int ran(0., nRegions-1, seed);

  int val = 3; // see Norberg et al. 2009

  int nbins_D1 = m_dd->nbins_D1();
  int nbins_D2 = m_dd->nbins_D2();

  for (int i=0; i<nMocks; i++) {

    coutCBL << "analysing mock: " << i << " of " << nMocks << "\r"; cout.flush();

    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2(), m_rr->angularUnits(), m_rr->angularWeight());
    auto dr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2(), m_rr->angularUnits(), m_rr->angularWeight());

    double nData_SS = 0., nData_SS_weighted = 0., nRandom_SS = 0., nRandom_SS_weighted = 0.;

    vector<int> w(nRegions, 0);
    for (size_t n=0; n<val*nRegions; n++)
      w[ran()] ++;

    for (size_t j=0; j<nRegions; j++) {

      if (w[j]>0) {

	nData_SS += w[j]*nData_reg[j];
	nData_SS_weighted += w[j]*nData_reg_weighted[j];
	nRandom_SS += w[j]*nRandom_reg[j];
	nRandom_SS_weighted += w[j]*nRandom_reg_weighted[j];

	for (size_t k=j; k<nRegions; k++) {
	  if(w[k]>0){
	    int index = j*nRegions-(j-1)*j/2+k-j;
	    double ww = w[j]*w[k];
	    for (int bin1=0; bin1<nbins_D1; bin1++) 
	      for (int bin2=0; bin2<nbins_D2; bin2++) {
		dd_SS->add_data2D(bin1, bin2, dd[index], ww);
		rr_SS->add_data2D(bin1, bin2, rr[index], ww);
	      }
	  }
	}
      }

      for (size_t k=0; k<nRegions; k++) {
	if (w[k]>0) {
	  int index = j*nRegions+k;
	  double ww = w[j]*w[k];
	  for (int bin1=0; bin1<nbins_D1; bin1++) 
	    for (int bin2=0; bin2<nbins_D2; bin2++){
	      dr_SS->add_data2D(bin1, bin2, dr[index], ww);
	    }
	}

      }
    }
    data.push_back(move(correlation_LandySzalayEstimator(dd_SS, rr_SS, dr_SS, nData_SS, nData_SS_weighted, nRandom_SS, nRandom_SS_weighted)));
  }

  return data;
}


