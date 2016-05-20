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
 *  @file CatalogueAnalysis/ThreePointCorrelation/ThreePointCorrelation_comoving_reduced.cpp
 *
 *  @brief Methods of the class ThreePointCorrelation_comoving_reduced used to
 *  measure the monopole of the three-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  ThreePointCorrelation_comoving_reduced used to measure the monopole of the
 *  three-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "ThreePointCorrelation_comoving_reduced.h"
#include "TwoPointCorrelation1D_monopole.h"

using namespace cosmobl;
using namespace catalogue;
using namespace triplets;
using namespace threept;


// ============================================================================


void cosmobl::threept::ThreePointCorrelation_comoving_reduced::measure (const string dir_output_triplets, const string dir_output_2pt, const vector<string> dir_input_triplets, const int count_ddd, const int count_rrr, const int count_ddr, const int count_drr, const bool tcount) 
{   

  // ----------- compute the connected three-point correlation function -----------
  
  ThreePointCorrelation_comoving_connected::measure(dir_output_triplets, dir_input_triplets, count_ddd, count_rrr, count_ddr, count_drr, tcount);
  m_scale = ThreePointCorrelation_comoving_connected::m_scale;
  
  
  // ----------- compute the two-point correlation function ----------- 

  Catalogue data = *m_data;
  Catalogue random = *m_random;
  double rMin = m_ddd->side_s();
  double rMax = (1+m_ddd->side_u())*m_ddd->side_s()*(1.+2.*m_ddd->perc_increase());
  double binSize = 0.05;
  double shift = 0.5;
  twopt::TwoPointCorrelation1D_monopole TwoP {data, random, _logarithmic_, rMin, rMax, binSize, shift};
  TwoP.measure(twopt::ErrorType::_Poisson_,dir_output_triplets, {}, par::defaultString,0, 1, 1, 1, tcount);

  vector<double> log_r(TwoP.dd()->nbins()), log_xi(TwoP.dd()->nbins());
  for (int i=0; i<TwoP.dd()->nbins(); i++) {
    log_r[i] = log10(TwoP.xx()[i]);
    log_xi[i] = log10(1.+TwoP.xi1D()[i]);
  }
  
  vector<double> values_interp(m_ddd->nbins()+2, 0.), xi_real_lin(m_ddd->nbins()+2, 0.);
  
  values_interp[0] = log10(m_ddd->side_s());
  values_interp[1] = log10(m_ddd->side_u()*m_ddd->side_s());
  
  for (int i=0; i<m_ddd->nbins(); i++) {
    double tmp_value = (m_ddd->side_s()+((i+0.5)*m_ddd->binSize()));
    values_interp[i+2] = log10(tmp_value);
  }

  string file_2pt = dir_output_2pt+"2ptCorrelation_3pt.dat";
  ofstream fout (file_2pt.c_str()); checkIO(file_2pt, 0);
 
  for (size_t i=0; i<values_interp.size(); i++) {
    xi_real_lin[i] = pow(10., interpolated(values_interp[i], log_r, log_xi, "Linear"))-1.;
    //cout << pow(10.,values_interp[i]) << " --- " << xi_real_lin[i] << endl;
    fout << pow(10.,values_interp[i]) << "     " << xi_real_lin[i] << endl;
  }
  
  fout.clear(); fout.close(); cout <<"I wrote the file "<<file_2pt<<endl;

  
  // ----------- compute the reduced three-point correlation function -----------

  m_QQ.resize(m_ddd->nbins()); m_error.resize(m_ddd->nbins());
  
  for (int i=0; i<m_ddd->nbins(); i++) {
    m_QQ[i] = m_zeta[i]/((xi_real_lin[0]*xi_real_lin[1])+(xi_real_lin[0]*xi_real_lin[i+2])+(xi_real_lin[1]*xi_real_lin[i+2]));
    m_error[i] = -1.; // work in progress...
  }
}


// ============================================================================


void cosmobl::threept::ThreePointCorrelation_comoving_reduced::write (const string dir, const string file, bool connected) const
{      
  checkDim(m_scale, m_ddd->nbins(), "scale");
  
  string file_out = dir+file;
  ofstream fout (file_out.c_str()); checkIO(file_out, 0);

  if (!connected) {  
    fout << "# scale  Q  error(work in progress)" << endl;
    for (size_t i=0; i<m_scale.size(); i++) 
      fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << m_scale[i] << "  " << setw(8) << m_QQ[i] << "  " << setw(8) << m_error[i] << endl;
  }
  
  else {
    fout << "# scale  z  error(work in progrss)  Q  error(work in progress)" << endl;
    for (size_t i=0; i<m_scale.size(); i++) 
      fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << m_scale[i] << "  " << setw(8) << m_zeta[i] << "  " << setw(8) << ThreePointCorrelation_comoving_connected::m_error[i] << setw(8) << m_QQ[i] << "  " << setw(8) << m_error[i] << endl;
  }
  
  fout.close(); cout << endl << "I wrote the file: " << file_out << endl << endl;
}  
