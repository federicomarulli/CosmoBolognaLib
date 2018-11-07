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
 *  CatalogueAnalysis/ThreePointCorrelation/ThreePointCorrelation_comoving_connected.cpp
 *
 *  @brief Methods of the class
 *  ThreePointCorrelation_comoving_connected used to measure the
 *  connected three-point correlation function in comoving coordinates
 *
 *  This file contains the implementation of the methods of the class
 *  ThreePointCorrelation_comoving_connected used to measure the
 *  connected three-point correlation function in comoving coordinates
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "Data1D.h"
#include "ThreePointCorrelation_comoving_connected.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace triplets;
using namespace measure;
using namespace threept;
using namespace glob;


// ============================================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_connected::set_parameters (const triplets::TripletType tripletType, const double side_s, const double side_u, const double perc_increase, const int nbins) 
{
  double r12 = side_s;
  double r12_binSize = r12*2.*perc_increase;
  double r13 = side_s*side_u;
  double r13_binSize = r13*2.*perc_increase;

  m_ddd = move(Triplet::Create(tripletType, r12, r12_binSize, r13, r13_binSize, nbins));
  m_rrr = move(Triplet::Create(tripletType, r12, r12_binSize, r13, r13_binSize, nbins));
  m_ddr = move(Triplet::Create(tripletType, r12, r12_binSize, r13, r13_binSize, nbins));
  m_drr = move(Triplet::Create(tripletType, r12, r12_binSize, r13, r13_binSize, nbins));
}


// ============================================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_connected::set_parameters (const triplets::TripletType tripletType, const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int nbins) 
{
  m_ddd = move(Triplet::Create(tripletType, r12, r12_binSize, r13, r13_binSize, nbins));
  m_rrr = move(Triplet::Create(tripletType, r12, r12_binSize, r13, r13_binSize, nbins));
  m_ddr = move(Triplet::Create(tripletType, r12, r12_binSize, r13, r13_binSize, nbins));
  m_drr = move(Triplet::Create(tripletType, r12, r12_binSize, r13, r13_binSize, nbins));
}


// ============================================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_connected::measure (const std::string dir_output_triplets, const std::vector<std::string> dir_input_triplets, const bool count_ddd, const bool count_rrr, const bool count_ddr, const bool count_drr, const bool tcount, const int seed) 
{
  (void)seed;
  
  // -------- count the data-data-data, random-random-random, data-data-random and data-random-random triplets, or read them from file -------- 
  
  count_allTriplets(dir_output_triplets, dir_input_triplets, count_ddd, count_rrr, count_ddr, count_drr, tcount);
  

  // ----------- compute the three-point correlation function ----------- 

  m_scale.resize(m_ddd->nbins()); m_zeta.resize(m_ddd->nbins()); m_error.resize(m_ddd->nbins());
 
  double nGal = m_data->weightedN();
  double nRan = m_random->weightedN();

  coutCBL << "# Galaxies: " << nGal << " ; # Random: " << nRan << endl;
  
  double norm1 = (double(nGal)*double(nGal-1)*double(nGal-2))/6.;
  double norm2 = (double(nGal)*double(nGal-1)*double(nRan))*0.5;
  double norm3 = (double(nGal)*double(nRan)*double(nRan-1))*0.5;
  double norm4 = (double(nRan)*double(nRan-1)*double(nRan-2))/6.;
  
  for (int i=0; i<m_ddd->nbins(); i++) {
    
    m_scale[i] = m_ddd->scale(i);
    
    if (m_ddd->TT1D(i)>0 && m_rrr->TT1D(i)>0) {  
      m_zeta[i] = ((m_ddd->TT1D(i)/norm1)/(m_rrr->TT1D(i)/norm4))-3.*((m_ddr->TT1D(i)/norm2)/(m_rrr->TT1D(i)/norm4))+3.*(((m_drr->TT1D(i)/norm3)/(m_rrr->TT1D(i)/norm4)))-1.;
      m_error[i] = -1.; // work in progress...
    }

  }

  m_dataset = move(unique_ptr<data::Data1D>(new data::Data1D(m_scale, m_zeta, m_error)));
}


// ============================================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_connected::measure (const std::vector<std::vector<double>> weight, const bool doJK, const std::string dir_output_triplets, const std::vector<std::string> dir_input_triplets, const bool count_ddd, const bool count_rrr, const bool count_ddr, const bool count_drr, const bool tcount, const int seed) 
{
  (void)seed;
  
  // -------- count the data-data-data, random-random-random, data-data-random and data-random-random triplets, or read them from file -------- 
  
  count_allTriplets_region (weight, dir_output_triplets, dir_input_triplets, count_ddd, count_rrr, count_ddr, count_drr, tcount);
  

  // ----------- compute the three-point correlation function ----------- 

  m_scale.resize(m_ddd->nbins()); m_zeta.resize(m_ddd->nbins()); m_error.resize(m_ddd->nbins());
 
  double nGal = m_data->weightedN();
  double nRan = m_random->weightedN();
  
  double norm1 = (double(nGal)*double(nGal-1)*double(nGal-2))/6.;
  double norm2 = (double(nGal)*double(nGal-1)*double(nRan))*0.5;
  double norm3 = (double(nGal)*double(nRan)*double(nRan-1))*0.5;
  double norm4 = (double(nRan)*double(nRan-1)*double(nRan-2))/6.;

  
  for (int i=0; i<m_ddd->nbins(); i++) {

    m_scale[i] = m_ddd->scale(i);
    
    if (m_ddd->TT1D(i)>0 && m_rrr->TT1D(i)>0) 
      m_zeta[i] = ((m_ddd->TT1D(i)/norm1)/(m_rrr->TT1D(i)/norm4))-3.*((m_ddr->TT1D(i)/norm2)/(m_rrr->TT1D(i)/norm4))+3.*(((m_drr->TT1D(i)/norm3)/(m_rrr->TT1D(i)/norm4)))-1.;

  }

  /// Compute resamplings and covariance matrix
  
  vector<vector<double>> resampling_threept(weight.size(), vector<double>(m_ddd->nbins(), 0));
  vector<long> region_list = m_data->region_list();

  vector<double> nData_reg_weighted, nRandom_reg_weighted;

  const int nRegions = m_data->nRegions(); 

  for (int i=0; i<nRegions; i++) {
    nData_reg_weighted.push_back(m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg_weighted.push_back(m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
  }


  for (size_t i=0; i<weight.size(); i++) {

    double nData = 0;
    double nRan = 0;

    for (size_t j=0; j<weight[i].size(); j++) {
      nData += weight[i][j]*nData_reg_weighted[j];
      nRan += weight[i][j]*nRandom_reg_weighted[j];
    }

    double norm1 = (double(nData)*double(nData-1)*double(nData-2))/6.;
    double norm2 = (double(nData)*double(nData-1)*double(nRan))*0.5;
    double norm3 = (double(nData)*double(nRan)*double(nRan-1))*0.5;
    double norm4 = (double(nRan)*double(nRan-1)*double(nRan-2))/6.;

    for (int j=0; j<m_ddd->nbins(); j++) 
      if (m_ddd_regions[i]->TT1D(j)>0 && m_rrr_regions[i]->TT1D(j)>0) 
	resampling_threept[i][j] = ((m_ddd_regions[i]->TT1D(j)/norm1)/(m_rrr_regions[i]->TT1D(j)/norm4))-3.*((m_ddr_regions[i]->TT1D(j)/norm2)/(m_rrr_regions[i]->TT1D(j)/norm4))+3.*(((m_drr_regions[i]->TT1D(j)/norm3)/(m_rrr_regions[i]->TT1D(j)/norm4)))-1.;
  }

  vector<vector<double>> cov_mat;
  cbl::covariance_matrix(resampling_threept, cov_mat, doJK);
  m_dataset = move(unique_ptr<data::Data1D>(new data::Data1D(m_scale, m_zeta, cov_mat)));

  m_dataset->error(m_error);
}
  
  

// ============================================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_connected::measure (const cbl::measure::ErrorType errorType, const std::string dir_output_triplets, const std::vector<std::string> dir_input_triplets, const int nResamplings, const bool count_ddd, const bool count_rrr, const bool count_ddr, const bool count_drr, const bool tcount, const int seed) 
{  

  switch(errorType) {
    
    case cbl::measure::ErrorType::_None_:
      {
	measure(dir_output_triplets, dir_input_triplets, count_ddd, count_rrr, count_ddr, count_drr, tcount);
	break;
      }
    
    case cbl::measure::ErrorType::_Jackknife_:
      {
	const int nRegions = m_data->nRegions();

	vector<vector<double>> weight(nRegions, vector<double>(nRegions, 1));
	for (int i=0; i<nRegions; i++)
	  weight[i][i] = 0;

	measure(weight, true, dir_output_triplets, dir_input_triplets, count_ddd, count_rrr, count_ddr, count_drr, tcount);
	break;
      }

    case cbl::measure::ErrorType::_Bootstrap_:
      {
	const int nRegions = m_data->nRegions();

	random::UniformRandomNumbers_Int ran(0., nRegions-1, seed);
	
	int val = 2; // see Norberg et al. 2009

	vector<vector<double>> weight(nResamplings, vector<double>(nRegions, 0));
	for (int i=0; i<nResamplings; i++)
	  for (int j=0; j<val*nRegions; j++)
	    weight[i][ran()] ++;

	measure(weight, false, dir_output_triplets, dir_input_triplets, count_ddd, count_rrr, count_ddr, count_drr, tcount);
	break;
      }

    default:
      ErrorCBL("Error in measure() of ThreePointCorrelation_comoving_connected, no such kind of error type!");
  }

}


// ============================================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_connected::write (const std::string dir, const std::string file) const
{      
  checkDim(m_scale, m_ddd->nbins(), "scale");

  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  fout << "# scale  zeta  error(work in progress)" << endl;
  
  for (size_t i=0; i<m_scale.size(); i++) 
    fout << setiosflags(ios::fixed) << setprecision(4) << setw(10) << right << m_scale[i]
	 << "   " << setiosflags(ios::fixed) << setprecision(4) << setw(10) << right << m_zeta[i]
	 << "   " << setiosflags(ios::fixed) << setprecision(4) << setw(10) << right << m_error[i] << endl;
    
  fout.close(); coutCBL << endl << "I wrote the file: " << file_out << endl << endl;
}  


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_connected::write_covariance (const std::string dir, const std::string file) const
{
  m_dataset->write_covariance(dir, file);
}


