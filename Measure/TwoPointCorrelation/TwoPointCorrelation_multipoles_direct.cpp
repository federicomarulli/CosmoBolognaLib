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
 *  CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation_multipoles_direct.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation_multipoles_direct
 *  used to measure the multipoles of the two-point correlation
 *  function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation_multipoles_direct used to measure the
 *  multipoles of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "Data1D_extra.h"
#include "TwoPointCorrelation_multipoles_direct.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace chainmesh;
using namespace data;
using namespace pairs;
using namespace measure;
using namespace twopt;

// ============================================================================


std::shared_ptr<data::Data> cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::data_with_extra_info (const std::shared_ptr<pairs::Pair> dd, const std::vector<double> rad, const std::vector<double> xi, const std::vector<double> error) const
{
  vector<vector<double>> extra(4);
  
  for (int l=0; l<3; l++)
    for (int i=0; i<dd->nbins(); ++i) {
      extra[0].push_back(dd->scale_mean(i));
      extra[1].push_back(dd->scale_sigma(i));
      extra[2].push_back(dd->z_mean(i));
      extra[3].push_back(dd->z_sigma(i));
    }

  
  return move(unique_ptr<data::Data1D_extra>(new data::Data1D_extra(rad, xi, error, extra)));
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::write_pairs (const std::shared_ptr<pairs::Pair> PP, const std::string dir, const std::string file) const 
{  
  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {}
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  
  // ----- standard info: scales at the bin centre + number of pairs -----

  if (PP->pairInfo()==PairInfo::_standard_)
    for (int i=0; i<PP->nbins(); i++) 
      fout << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << i
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->scale(i)
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->PP1D(i)
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->PP1D_weighted(i)
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->PP1D(i+PP->nbins()+1)
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->PP1D_weighted(i+PP->nbins()+1)
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->PP1D(i+2*(PP->nbins()+1))
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->PP1D_weighted(i+2*(PP->nbins()+1)) << endl;
  else if (PP->pairInfo()==PairInfo::_extra_) 
    for (int i=0; i<PP->nbins(); i++)
      fout << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << i
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->scale(i)
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->PP1D(i)
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->PP1D_weighted(i)
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->PP1D(i+PP->nbins()+1)
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->PP1D_weighted(i+PP->nbins()+1)
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->PP1D(i+2*(PP->nbins()+1))
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->PP1D_weighted(i+2*(PP->nbins()+1))
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->scale_mean(i)
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->scale_sigma(i)
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->z_mean(i)
	<< "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP->z_sigma(i) << endl;
  else
      ErrorCBL("Error in write_pairs() of TwoPointCorrelation_multipoles_direct.cpp: no such pairInfo!");

  fout.clear(); fout.close(); coutCBL << "I wrote the file " << file_out << endl;
}
  
// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::read_pairs (std::shared_ptr<pairs::Pair> PP, const std::vector<std::string> dir, const std::string file) const
{
  if (dir.size()==0)
    ErrorCBL("Error in cbl::measure::twopt::TwoPointCorrelation1D::read_pairs of TwoPointCorrelation_multipoles_direct! dir.size()=0!");

  string line;
  int i;
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
	
	if (num.size()!=8)
	  ErrorCBL("Error in read_pairs() of TwoPointCorrelation1D.cpp: the number of lines in the input pair file: "+file_in+" must be 8!");
	
	i = int(num[0]);
	pairs = num[2];
	weighted_pairs = num[3];
	
	PP->add_data1D(i, {pairs, weighted_pairs});

	pairs = num[4];
	weighted_pairs = num[5];
	
	PP->add_data1D(PP->nbins()+1+i, {pairs, weighted_pairs});

	pairs = num[6];
	weighted_pairs = num[7];
	
	PP->add_data1D(2*(PP->nbins()+1)+i, {pairs, weighted_pairs});
      }
	  
      fin.clear(); fin.close(); coutCBL << "I read the file " << file_in << endl;
    }

  // ----- standard + extra info -----

  else if (PP->pairInfo()==PairInfo::_extra_)

    for (size_t dd=0; dd<dir.size(); dd++) {

      string file_in = dir[dd]+file; 
      ifstream fin(file_in.c_str()); checkIO(fin, file_in);

      double scale_mean, scale_sigma, z_mean, z_sigma;

      while (getline(fin, line)) {

	stringstream ss(line);
	vector<double> num;
	double val;
	while (ss >> val) num.emplace_back(val);

	if (num.size()!=12)
	  ErrorCBL("Error in read_pairs() of TwoPointCorrelation_multipoles_direct.cpp: the number of lines in the input pair file: "+file_in+" must be 12!");

	i = int(num[0]);

	scale_mean = num[8];
	scale_sigma = num[9];
	z_mean = num[10];
	z_sigma = num[11];

	pairs = num[2];
	weighted_pairs = num[3];
	
	PP->add_data1D(i, {pairs, weighted_pairs, scale_mean, pow(scale_sigma, 2)*weighted_pairs, z_mean, pow(z_sigma, 2)*weighted_pairs});

	pairs = num[4];
	weighted_pairs = num[5];
	
	PP->add_data1D(PP->nbins()+1+i, {pairs, weighted_pairs, scale_mean, pow(scale_sigma, 2)*weighted_pairs, z_mean, pow(z_sigma, 2)*weighted_pairs});

	pairs = num[6];
	weighted_pairs = num[7];
	
	PP->add_data1D(2*(PP->nbins()+1)+i, {pairs, weighted_pairs, scale_mean, pow(scale_sigma, 2)*weighted_pairs, z_mean, pow(z_sigma, 2)*weighted_pairs});
      }

      fin.clear(); fin.close(); coutCBL << "I read the file " << file_in << endl;
    }

  else
    ErrorCBL("Error in read_pairs() of TwoPointCorrelation_multipoles_direct.cpp: no such pairInfo!");
    
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::write_pairs (const std::vector<std::shared_ptr<pairs::Pair>> PP, const std::string dir, const std::string file) const 
{  
  size_t nRegions = m_data->region_list().size();

  bool cross = (PP.size() == nRegions*nRegions) ? true : false;
  
  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {}
  
  string file_out = dir+file;  
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  // ----- standard info: scales at the bin centre + number of pairs -----

  if (PP[0]->pairInfo()==PairInfo::_standard_)
    
    for (size_t i=0; i<nRegions; i++) 
      for (size_t j=(cross) ? 0 : i; j<nRegions; j++) {
	int index = (cross) ? i*nRegions+j : i*nRegions+j-(i-1)*i/2-i;
	for (int r1=0; r1<PP[index]->nbins(); r1++)
	  if (PP[index]->PP1D(r1)>0)
	    fout << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << i
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << j
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << r1
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->scale(r1)
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->PP1D(r1)
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->PP1D_weighted(r1)
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->PP1D(r1+PP[index]->nbins()+1)
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->PP1D_weighted(r1+PP[index]->nbins()+1)
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->PP1D(r1+2*(PP[index]->nbins()+1))
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->PP1D_weighted(r1+2*(PP[index]->nbins()+1)) << endl;
      }


  // ----- standard + extra info -----

  else if (PP[0]->pairInfo()==PairInfo::_extra_)

    for (size_t i=0; i<nRegions; i++) 
      for (size_t j=(cross) ? 0 : i; j<nRegions; j++) {
	int index = (cross) ? i*nRegions+j : i*nRegions+j-(i-1)*i/2-i;
	for (int r1=0; r1<PP[index]->nbins(); r1++)
	  if (PP[index]->PP1D(r1)>0)
	    fout << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << i
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << j
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << r1
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->scale(r1)
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->PP1D(r1)
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->PP1D_weighted(r1)
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->PP1D(r1+PP[index]->nbins()+1)
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->PP1D_weighted(r1+PP[index]->nbins()+1)
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->PP1D(r1+2*(PP[index]->nbins()+1))
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->PP1D_weighted(r1+2*(PP[index]->nbins()+1))
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->scale_mean(r1)
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->scale_sigma(r1)
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->z_mean(r1)
	      << "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << PP[index]->z_sigma(r1) << endl;
      }

  else
    ErrorCBL("Error in write_pairs() of TwoPointCorrelation_multipoles_direct.cpp: no such pairInfo!");
      
  fout.clear(); fout.close();
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::read_pairs (std::vector<std::shared_ptr<pairs::Pair>> PP, const std::vector<std::string> dir, const std::string file) const
{
  size_t nRegions = m_data->region_list().size();

  bool cross = (PP.size() == nRegions*nRegions) ? true : false;

  int i, j, bin, index;
  double rad, pairs0, weighted_pairs0,  pairs2, weighted_pairs2, pairs4, weighted_pairs4;
    
  // ----- standard info: scales at the bin centre + number of pairs -----

  if (PP[0]->pairInfo()==PairInfo::_standard_)
    
    for (size_t dd=0; dd<dir.size(); dd++) {
      string ff = dir[dd]+file; 
      ifstream fin(ff.c_str()); checkIO(fin, ff);
      
      while (fin >> i >> j >> bin >> rad >> pairs0 >> weighted_pairs0 >> pairs2 >> weighted_pairs2 >> pairs4 >> weighted_pairs4) {
	index = (cross) ? i*nRegions+j : i*nRegions+j-(i-1)*i/2-i;
	PP[index]->add_data1D(bin, {pairs0, weighted_pairs0});
	PP[index]->add_data1D(bin+PP[index]->nbins()+1, {pairs2, weighted_pairs2});
	PP[index]->add_data1D(bin+2*(PP[index]->nbins()+1), {pairs4, weighted_pairs4});
      }
      
      fin.clear(); fin.close(); coutCBL << "I read the file " << ff << endl;
    }
  

  // ----- standard + extra info -----

  else if (PP[0]->pairInfo()==PairInfo::_extra_)

    for (size_t dd=0; dd<dir.size(); dd++) {
      string ff = dir[dd]+file; 
      ifstream fin(ff.c_str()); checkIO(fin, ff);
      
      double scale_mean, scale_sigma, z_mean, z_sigma;

      while (fin >> i >> j >> bin >> rad >> pairs0 >> weighted_pairs0 >> pairs2 >> weighted_pairs2 >> pairs4 >> weighted_pairs4 >> scale_mean >> scale_sigma >> z_mean >> z_sigma) {
	index = (cross) ? i*nRegions+j : i*nRegions+j-(i-1)*i/2-i;
	PP[index]->add_data1D(bin, {pairs0, weighted_pairs0, scale_mean, pow(scale_sigma, 2)*weighted_pairs0, z_mean, pow(z_sigma, 2)*weighted_pairs0});
	PP[index]->add_data1D(bin+PP[index]->nbins()+1, {pairs2, weighted_pairs2, scale_mean, pow(scale_sigma, 2)*weighted_pairs0, z_mean, pow(z_sigma, 2)*weighted_pairs0});
	PP[index]->add_data1D(bin+2*(PP[index]->nbins()+1), {pairs4, weighted_pairs4, scale_mean, pow(scale_sigma, 2)*weighted_pairs0, z_mean, pow(z_sigma, 2)*weighted_pairs0});
      } 
      
      fin.clear(); fin.close(); coutCBL << "I read the file " << ff << endl;
    }
  
  else
    ErrorCBL("Error in read_pairs() of TwoPointCorrelation_multipoles_direct.cpp: no such pairInfo!");
}


// ============================================================================


std::shared_ptr<data::Data> cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::correlation_NaturalEstimator (const std::shared_ptr<pairs::Pair> dd, const std::shared_ptr<pairs::Pair> rr, const int nData, const double nData_weighted, const int nRandom, const double nRandom_weighted)
{
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
  
  vector<double> rad(3*m_dd->nbins()), xi(3*m_dd->nbins(), -1.), error(3*m_dd->nbins(), 1000.);

  for (int l=0; l<3; l++) {
    for (int i=0; i<dd->nbins(); i++) {
      int idx = i+l*(dd->nbins()+1);

      rad[i+l*dd->nbins()] = dd->scale(i+l*dd->nbins());

      if (abs(dd->PP1D_weighted(idx))>0) {

        if (abs(rr->PP1D_weighted(i))<1.e-30) 
          ErrorCBL("Error in correlation_NaturalEstimator() of TwoPointCorrelation_multipoles_direct.cpp: there are no random objects in the bin "+conv(i, par::fINT)+"; please, either increase the total number of random objects or enlarge the bin size! (dd="+conv(dd->PP1D_weighted(i), par::fDP3)+", rr="+conv(rr->PP1D_weighted(i), par::fDP3)+")");

        // normalised number of data-data weighted pairs
        double DD_norm = dd->PP1D_weighted(idx)*nDDi;

        // normalised number of random-random weighted pairs
        double RR_iso = rr->PP1D_weighted(i)*nRRi;

        // Landy & Szalay estimator
        xi[i+l*dd->nbins()] = (DD_norm/RR_iso)-1;

        // Poisson error
        error[i+l*dd->nbins()] = PoissonError(Estimator::_LandySzalay_, dd->PP1D(i), rr->PP1D(i), 0, nD, nR); ///CHECK!!!!

      }
    }
  }

  return (!m_compute_extra_info) ? move(unique_ptr<data::Data1D>(new data::Data1D(rad, xi, error))) : data_with_extra_info(dd, rad, xi, error);

}


// ============================================================================


std::shared_ptr<data::Data> cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::correlation_LandySzalayEstimator (const std::shared_ptr<pairs::Pair> dd, const std::shared_ptr<pairs::Pair> rr, const std::shared_ptr<pairs::Pair> dr, const int nData, const double nData_weighted, const int nRandom, const double nRandom_weighted)
{
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
  
  vector<double> rad(3*m_dd->nbins()), xi(3*m_dd->nbins(), -1.), error(3*m_dd->nbins(), 1000.);

  for (int l=0; l<3; l++) {
    for (int i=0; i<dd->nbins(); i++) {
      int idx = i+l*(dd->nbins()+1);

      rad[i+l*dd->nbins()] = dd->scale(i+l*dd->nbins());

//if (abs(dd->PP1D_weighted(i))>0) {

        if (abs(rr->PP1D_weighted(i))<1.e-30) 
          ErrorCBL("Error in correlation_LandySzalayEstimator() of TwoPointCorrelation_multipoles_direct.cpp: there are no random objects in the bin "+conv(i, par::fINT)+"; please, either increase the total number of random objects or enlarge the bin size! (dd="+conv(dd->PP1D_weighted(idx), par::fDP3)+", rr="+conv(rr->PP1D_weighted(idx), par::fDP3)+")");

        // normalised number of data-data weighted pairs
        double DD_norm = dd->PP1D_weighted(idx)*nDDi;

        // normalised number of random-random weighted pairs
        double RR_norm = rr->PP1D_weighted(idx)*nRRi;

        // normalised number of random-random weighted pairs
        double RR_iso = rr->PP1D_weighted(i)*nRRi;

        // normalised number of data-random weighted pairs
        double DR_norm = dr->PP1D_weighted(idx)*nDRi;

        // Landy & Szalay estimator
        xi[i+l*dd->nbins()] = (DD_norm+RR_norm-2.*DR_norm)/RR_iso;

        // Poisson error
        error[i+l*dd->nbins()] = PoissonError(Estimator::_LandySzalay_, dd->PP1D(i), rr->PP1D(i), dr->PP1D(i), nD, nR); ///CHECK!!!!

      //}

    }
  }

  return (!m_compute_extra_info) ? move(unique_ptr<data::Data1D>(new data::Data1D(rad, xi, error))) : data_with_extra_info(dd, rad, xi, error);

}


// ============================================================================


std::vector<std::shared_ptr<data::Data>> cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::XiJackknife (const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr)
{
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<shared_ptr<data::Data>> data;

  for (size_t i=0; i<nRegions; i++) {

    coutCBL << "analysing region: " << i << " of " << nRegions << "\r"; cout.flush();
    
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin(), m_dd->sMax(), m_dd->nbins(), m_dd->shift(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin(), m_rr->sMax(), m_rr->nbins(), m_rr->shift(), m_rr->angularUnits(), m_rr->angularWeight());
    
    vector<int> w(nRegions, 1);
    w[i] = 0;
    
    for (size_t j=0; j<nRegions; j++) 
      for (size_t k=j; k<nRegions; k++) {
	int index = j*nRegions-(j-1)*j/2+k-j;
	double ww = w[j]*w[k];
	if (ww>0)
	  for (size_t bin=0; bin<dd_SS->PP1D().size(); bin++) {
	    dd_SS->add_data1D(bin, dd[index]);
	    rr_SS->add_data1D(bin, rr[index]);
	  }
      }

    double nData_SS = m_data->nObjects_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    double nData_SS_weighted = m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    double nRandom_SS = m_random->nObjects_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    double nRandom_SS_weighted = m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    
    data.push_back(move(correlation_NaturalEstimator(dd_SS, rr_SS, nData_SS, nData_SS_weighted, nRandom_SS, nRandom_SS_weighted)));
  } 
  
  return data;
}


// ============================================================================


std::vector<std::shared_ptr<data::Data>> cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::XiJackknife (const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr, const std::vector<std::shared_ptr<pairs::Pair>> dr)
{
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<shared_ptr<data::Data>> data;

  for (size_t i=0; i<nRegions; i++) {
    
    coutCBL << "analysing region: " << i << " of " << nRegions << "\r"; cout.flush();
      
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin(), m_dd->sMax(), m_dd->nbins(), m_dd->shift(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin(), m_rr->sMax(), m_rr->nbins(), m_rr->shift(), m_rr->angularUnits(), m_rr->angularWeight());
    auto dr_SS = Pair::Create(m_dr->pairType(), m_dr->pairInfo(), m_dr->sMin(), m_dr->sMax(), m_dr->nbins(), m_dr->shift(), m_dr->angularUnits(), m_dr->angularWeight());

    vector<int> w(nRegions, 1);
    w[i] = 0;

    for (size_t j=0; j<nRegions; j++) {

      if (w[j]>0) {
        for (size_t k=j; k<nRegions; k++) { // auto pairs
          if (w[k]>0) {
            int index = j*nRegions+k-(j-1)*j/2-j;
            for (size_t bin=0; bin<dd_SS->PP1D().size(); bin++) {
              dd_SS->add_data1D(bin, dd[index]);
              rr_SS->add_data1D(bin, rr[index]);
            }
          }
        }

        for (size_t k=0; k<nRegions; k++) {// cross pairs
          if (w[k]>0) {
            int index = j*nRegions+k;
            for (size_t bin=0; bin<dd_SS->PP1D().size(); bin++) 
              dr_SS->add_data1D(bin, dr[index]);
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


std::vector<std::shared_ptr<data::Data>> cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::XiBootstrap (const int nMocks, const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr, const int seed)
{
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<shared_ptr<data::Data>> data;
  vector<double> nData_reg, nData_reg_weighted, nRandom_reg, nRandom_reg_weighted;

  for (size_t i=0; i<nRegions; i++) {
    nData_reg.push_back(m_data->nObjects_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nData_reg_weighted.push_back(m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg.push_back(m_random->nObjects_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg_weighted.push_back(m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
  }

  random::UniformRandomNumbers_Int ran(0., nRegions-1, seed);
  
  int val = 3; // see Norberg et al. 2009

  for (int i=0; i<nMocks; i++) {

    coutCBL << "analysing mock: " << i << " of " << nMocks << "\r"; cout.flush();
    
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin(), m_dd->sMax(), m_dd->nbins(), m_dd->shift(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin(), m_rr->sMax(), m_rr->nbins(), m_rr->shift(), m_rr->angularUnits(), m_rr->angularWeight());

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
          if (w[k]>0) {
            int index = j*nRegions+k-(j-1)*j/2-j;
            double ww = w[j]*w[k]; //(k==j) ? w[k] : w[j]*w[k];
            for (size_t bin=0; bin<dd_SS->PP1D().size(); bin++) {
              dd_SS->add_data1D(bin, dd[index], ww);
              rr_SS->add_data1D(bin, rr[index], ww);
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


std::vector<std::shared_ptr<data::Data>> cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::XiBootstrap (const int nMocks, const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr, const std::vector<std::shared_ptr<pairs::Pair>> dr, const int seed)
{
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<shared_ptr<data::Data>> data;
  vector<double> nData_reg, nData_reg_weighted, nRandom_reg, nRandom_reg_weighted;

  for (size_t i=0; i<nRegions; i++) {
    nData_reg.push_back(m_data->nObjects_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nData_reg_weighted.push_back(m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg.push_back(m_random->nObjects_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg_weighted.push_back(m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
  }

  random::UniformRandomNumbers_Int ran(0., nRegions-1, seed);
  
  int val = 3; // see Norberg et al. 2009

  for (int i=0; i<nMocks; i++) {

    coutCBL << "analysing mock: " << i << " of " << nMocks << "\r"; cout.flush();
    
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin(), m_dd->sMax(), m_dd->nbins(), m_dd->shift(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin(), m_rr->sMax(), m_rr->nbins(), m_rr->shift(), m_rr->angularUnits(), m_rr->angularWeight());
    auto dr_SS = Pair::Create(m_dr->pairType(), m_dr->pairInfo(), m_dr->sMin(), m_dr->sMax(), m_dr->nbins(), m_dr->shift(), m_dr->angularUnits(), m_dr->angularWeight());

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
          if (w[k]>0) {
            int index = j*nRegions+k-(j-1)*j/2-j;
            double ww = w[j]*w[k]; //(k==j) ? w[k] : w[j]*w[k];
            for (size_t bin=0; bin<dd_SS->PP1D().size(); bin++) {
              dd_SS->add_data1D(bin, dd[index], ww);
              rr_SS->add_data1D(bin, rr[index], ww);
            }
          }
        }

        for (size_t k=0; k<nRegions; k++) {
          if (w[k]>0) {
            int index = j*nRegions+k;
            double ww = w[j]*w[k]; //(k==j) ? w[k] : w[j]*w[k];
            for (size_t bin=0; bin<dd_SS->PP1D().size(); bin++) 
              dr_SS->add_data1D(bin, dr[index], ww);
          }
        }
      }
    }
    
    data.push_back(move(correlation_LandySzalayEstimator(dd_SS, rr_SS, dr_SS, nData_SS, nData_SS_weighted, nRandom_SS, nRandom_SS_weighted)));
  }

  return data;
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::set_parameters (const BinType binType, const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight, const bool compute_extra_info) 
{
  if (!compute_extra_info) 
    m_dd = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_comoving_multipoles_log_, PairInfo::_standard_, rMin, rMax, nbins, shift, angularUnits, angularWeight))
      : move(Pair::Create(PairType::_comoving_multipoles_lin_, PairInfo::_standard_, rMin, rMax, nbins, shift, angularUnits, angularWeight));
  else 
    m_dd = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_comoving_multipoles_log_, PairInfo::_extra_, rMin, rMax, nbins, shift, angularUnits, angularWeight))
      : move(Pair::Create(PairType::_comoving_multipoles_lin_, PairInfo::_extra_, rMin, rMax, nbins, shift, angularUnits, angularWeight));
  
  m_rr = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_comoving_multipoles_log_, PairInfo::_standard_, rMin, rMax, nbins, shift, angularUnits))
    : move(Pair::Create(PairType::_comoving_multipoles_lin_, PairInfo::_standard_, rMin, rMax, nbins, shift, angularUnits));
  
  m_dr = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_comoving_multipoles_log_, PairInfo::_standard_, rMin, rMax, nbins, shift, angularUnits))
    : move(Pair::Create(PairType::_comoving_multipoles_lin_, PairInfo::_standard_, rMin, rMax, nbins, shift, angularUnits));
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::set_parameters (const BinType binType, const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight, const bool compute_extra_info)
{
  if (!compute_extra_info) 
    m_dd = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_comoving_multipoles_log_, PairInfo::_standard_, rMin, rMax, binSize, shift, angularUnits, angularWeight))
      : move(Pair::Create(PairType::_comoving_multipoles_lin_, PairInfo::_standard_, rMin, rMax, binSize, shift, angularUnits, angularWeight));
  else 
    m_dd = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_comoving_multipoles_log_, PairInfo::_extra_, rMin, rMax, binSize, shift, angularUnits, angularWeight))
      : move(Pair::Create(PairType::_comoving_multipoles_lin_, PairInfo::_extra_, rMin, rMax, binSize, shift, angularUnits, angularWeight));
    
  m_rr = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_comoving_multipoles_log_, PairInfo::_standard_, rMin, rMax, binSize, shift, angularUnits))
    : move(Pair::Create(PairType::_comoving_multipoles_lin_, PairInfo::_standard_, rMin, rMax, binSize, shift, angularUnits));
  
  m_dr = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_comoving_multipoles_log_, PairInfo::_standard_, rMin, rMax, binSize, shift, angularUnits))
    : move(Pair::Create(PairType::_comoving_multipoles_lin_, PairInfo::_standard_, rMin, rMax, binSize, shift, angularUnits));
}


// ============================================================================================


std::vector<double> cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::xx () const
{
  vector<double> rad, xx;
  m_dataset->xx(xx);

  for (size_t i=0; i<xx.size()/3; i++)
    rad.push_back(xx[i]);

  return rad;
}


// ============================================================================================


std::vector<double> cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::xiMonopole () const
{
  vector<double> vv; 
  m_dataset->data(vv);

  size_t sz = vv.size();

  vector<double> xi0;
  for (size_t i=0; i<sz/3; i++)
    xi0.push_back(vv[i]);

  return xi0;
}


// ============================================================================================


std::vector<double> cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::errorMonopole () const
{
  vector<double> vv; 
  m_dataset->error(vv);

  size_t sz = vv.size();

  vector<double> error_xi0;

  for (size_t i=0; i<sz/3; i++)
    error_xi0.push_back(vv[i]);

  return error_xi0;
}


// ============================================================================================


std::vector<double> cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::xiQuadrupole () const 
{
  vector<double> vv; 
  m_dataset->data(vv);

  size_t sz = vv.size();

  vector<double> xi2;

  for (size_t i=sz/3; i<2*sz/3; i++)
    xi2.push_back(vv[i]);

  return xi2;
}


// ============================================================================================


std::vector<double> cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::errorQuadrupole () const 
{
  vector<double> vv; 
  m_dataset->error(vv);

  size_t sz = vv.size();

  vector<double> error_xi2;

  for (size_t i=sz/3; i<2*sz/3; i++)
    error_xi2.push_back(vv[i]);

  return error_xi2;
}


// ============================================================================================


std::vector<double> cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::xiHexadecapole () const
{
  vector<double> vv; 
  m_dataset->data(vv);

  size_t sz = vv.size();

  vector<double> xi4;

  for (size_t i=2*sz/3; i<sz; i++)
    xi4.push_back(vv[i]);

  return xi4;
}


// ============================================================================================


std::vector<double> cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::errorHexadecapole () const 
{
  vector<double> vv; 
  m_dataset->error(vv);

  size_t sz = vv.size();

  vector<double> error_xi4;
  
  for (size_t i=2*sz/3; i<sz; i++)
    error_xi4.push_back(vv[i]);

  return error_xi4;
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation_multipoles_direct::write (const std::string dir, const std::string file, const int rank) const 
{
  (void)rank;
  
  vector<double> rad; m_dataset->xx(rad);
  vector<double> xil = m_dataset->data();
  vector<double> error = m_dataset->error();

  checkDim(rad, m_dd->nbins()*3, "rad");
  
  const string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  const int precision = 5;
 
  string header = "[1] separation at the bin centre # [2] monopole # [3] error on the monopole # [4] quadrupole # [5] error on the quadrupole # [6] hexadecapole # [7] error on the hexadecapole";
  
  if (m_compute_extra_info) header += " # [8] mean separation # [9] standard deviation of the separation distribution # [10] mean redshift # [11] standard deviation of the redshift distribution";
  
  fout << "### " << header << " ###" <<endl;

  for (int i=0; i<m_dd->nbins(); i++) {
    fout << setiosflags(ios::fixed) << setprecision(precision) << setw(10) << right << rad[i]
	 << "  "  << setiosflags(ios::fixed) << setprecision(precision) << setw(10) << right << xil[i]
	 << "  "  << setiosflags(ios::fixed) << setprecision(precision) << setw(10) << right << error[i]
	 << "  "  << setiosflags(ios::fixed) << setprecision(precision) << setw(10) << right << xil[i+m_dd->nbins()]
	 << "  "  << setiosflags(ios::fixed) << setprecision(precision) << setw(10) << right << error[i+m_dd->nbins()]
	 << "  "  << setiosflags(ios::fixed) << setprecision(precision) << setw(10) << right << xil[i+2*m_dd->nbins()]
	 << "  "  << setiosflags(ios::fixed) << setprecision(precision) << setw(10) << right << error[i+2*m_dd->nbins()];
    if (m_compute_extra_info)
      for (size_t ex=0; ex<m_dataset->extra_info().size(); ++ex)
	fout << setiosflags(ios::fixed) << setprecision(precision) << setw(10) << right << m_dataset->extra_info(ex, i);
    fout << endl;
  }
   
  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl << endl;
}
