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
 *  @file CatalogueAnalysis/ThreePointCorrelation/ThreePointCorrelation_comoving_connected.cpp
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


#include "ThreePointCorrelation_comoving_connected.h"

using namespace cosmobl;
using namespace catalogue;
using namespace triplets;
using namespace threept;


// ============================================================================================


void cosmobl::threept::ThreePointCorrelation_comoving_connected::set_parameters (const TripletType tripletType, const double side_s, const double side_u, const double perc_increase, const int nbins) 
{
  m_ddd = move(Triplet::Create(tripletType, side_s, side_u, perc_increase, nbins));
  m_rrr = move(Triplet::Create(tripletType, side_s, side_u, perc_increase, nbins));
  m_ddr = move(Triplet::Create(tripletType, side_s, side_u, perc_increase, nbins));
  m_drr = move(Triplet::Create(tripletType, side_s, side_u, perc_increase, nbins));
}


// ============================================================================================


void cosmobl::threept::ThreePointCorrelation_comoving_connected::measure (const string dir_output_triplets, const vector<string> dir_input_triplets, const int count_ddd, const int count_rrr, const int count_ddr, const int count_drr, const bool tcount) 
{  

  // -------- count the data-data-data, random-random-random, data-data-random and data-random-random triplets, or read them from file -------- 
  
  count_allTriplets(dir_output_triplets, dir_input_triplets, count_ddd, count_rrr, count_ddr, count_drr, tcount);
  

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
    
    if (m_ddd->TT1D(i)>0 && m_rrr->TT1D(i)>0) {  
      m_zeta[i] = ((m_ddd->TT1D(i)/norm1)/(m_rrr->TT1D(i)/norm4))-3.*((m_ddr->TT1D(i)/norm2)/(m_rrr->TT1D(i)/norm4))+3.*(((m_drr->TT1D(i)/norm3)/(m_rrr->TT1D(i)/norm4)))-1.;
      m_error[i] = -1.; // work in progress...
    }

  }
}


// ============================================================================================


void cosmobl::threept::ThreePointCorrelation_comoving_connected::write (const string dir, const string file) const
{      
  checkDim(m_scale, m_ddd->nbins(), "scale");

  string file_out = dir+file;
  ofstream fout (file_out.c_str()); checkIO(file_out, 0);

  fout << "# scale  zeta  error(work in progress)" << endl;
  
  for (size_t i=0; i<m_scale.size(); i++) 
    fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << m_scale[i] << "  " << setw(8) << m_zeta[i] << "  " << setw(8) << m_error[i] << endl;
    
  fout.close(); cout << endl << "I wrote the file: " << file_out << endl << endl;
}  
