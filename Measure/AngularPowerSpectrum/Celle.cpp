/*******************************************************************
 *  Copyright (C) 2015 by Federico Marulli                         *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file Measure/AngularPowerSpectrum/Celle.cpp
 *
 *  @brief Methods of the class PowerSpectrum_angular  
 *
 *  This file contains the implementation of all the methods of the
 *  classes PowerSpectrum_angular*, used to measure the angular power spectrum
 *
 *  @author Federico Marulli, Massimiliano Romanello
 *
 *  @author federico.marulli3@unibo.it massimilia.romanell2@unibo.it
 */

#include "PowerSpectrum_Angular.h"
#include "Data1D.h"
#include <math.h>
#include <boost/math/special_functions/bessel.hpp>
#include <fstream>

using namespace std;

using namespace cbl;
using namespace measure;
using namespace twopt;
using namespace boost::math;


// ============================================================================================


cbl::measure::angularpk::PowerSpectrum_angular::PowerSpectrum_angular (const catalogue::Catalogue data, const catalogue::Catalogue random, const double ell_min, const double ell_max, const int Nell, const BinType binType, const double thetaMin, const double thetaMax, const int nbins, const double shift, const CoordinateUnits angularUnits)
{
  m_ell_min = ell_min;   
  m_ell_max = ell_max;   
  m_Nell = Nell;
  m_binType = binType;
  m_angularUnits = angularUnits;
  set_ell_output(binType);
  cbl::measure::twopt::TwoPointCorrelation1D_angular TwoPointCorrelation1D_angular(data, random, binType, thetaMin, thetaMax, nbins, shift, angularUnits);
  m_TwoPointCorrelation1D_angular = std::make_shared<twopt::TwoPointCorrelation1D_angular>(twopt::TwoPointCorrelation1D_angular(std::move(TwoPointCorrelation1D_angular)));
}


// ============================================================================================


double cbl::measure::angularpk::PowerSpectrum_angular::m_dtheta_theta (const BinType binType, int i)
{     //Returns dtheta_i * theta_i
  double value=0;
  if (binType==cbl::BinType::_linear_)      //dtheta * theta
    value = (m_theta[1]-m_theta[0])*m_theta[i];   

  else if (binType==cbl::BinType::_logarithmic_) // d log theta * theta^2 = d theta / theta * theta^2 = dtheta * theta
    value = (log(m_theta[1])-log(m_theta[0]))*m_theta[i]*m_theta[i];
  
  return value;
}


// ============================================================================================


void cbl::measure::angularpk::PowerSpectrum_angular::set_ell_output (const BinType binType)
{
  int size_bin;
  if (binType==cbl::BinType::_linear_) {    //multipole scale is linear 
    size_bin = (m_ell_max-m_ell_min)/m_Nell;   //if linear binning
    for(int i=0; i<m_Nell; i++){
      m_ell.emplace_back(m_ell_min+size_bin*i);
    }
    
    /* int i = 0;
    m_ell.push_back(m_ell_min);
    do{
      i++;
      m_ell.push_back(m_ell_min+size_bin*i);
      }while(m_ell_min+size_bin*i<m_ell_max);*/
    }
  
  if (binType==cbl::BinType::_logarithmic_) {    //multipole scale is logarithmic
    const double ell_min_log = std::log10(m_ell_min);
    const double ell_max_log = std::log10(m_ell_max);
    const double log_increment = (ell_max_log - ell_min_log) / m_Nell;
    
    for (int i=0; i<m_Nell+1; ++i) 
      m_ell.emplace_back(std::pow(10, ell_min_log+i*log_increment));
  }
}


// ============================================================================================


void cbl::measure::angularpk::PowerSpectrum_angular::write (const std::string dir, const std::string file)
{

  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {}
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);
     
  fout << "#l \t Celle/2/pi \t Celle_error"<<endl;

  for (size_t j=0; j<m_Celle.size()/2; j++) {
    fout<<setiosflags(ios::fixed) << setprecision(5) << setw(10) << right<<m_ell[j]<<
      "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right <<m_Celle[j] <<
      "   " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right <<m_Celle[m_Celle.size()/2+j]<<endl;
  }
  
  fout.clear(); fout.close(); coutCBL << "I wrote the file " << file_out << endl;

}


// ============================================================================================

std::vector<double> cbl::measure::angularpk::PowerSpectrum_angular::m_error (const BinType binType, std::vector<double> w_err)
{
  std::vector<double> Celle_err(m_ell.size(), 0.0); 
  double dlogell = log(m_ell[1])-log(m_ell[0]);    //logarithmic bin amplitude
  double p0 = 1/dlogell;     
  double d, pref, Gp, fp;
  for (size_t j=0; j<Celle_err.size(); j++) {    
    for (size_t i=0; i<w_err.size();i++) {   
      d = m_dtheta_theta(binType, i);
      pref = p0*d/(m_theta[i]*m_theta[i]);
      Gp = m_theta[i]*m_ell[j]*sqrt(std::exp(dlogell))*cyl_bessel_j(1,m_theta[i]*m_ell[j]*sqrt(std::exp(dlogell)))-m_theta[i]*m_ell[j]/sqrt(std::exp(dlogell))*cyl_bessel_j(1,m_theta[i]*m_ell[j]/sqrt(std::exp(dlogell)));
      fp = m_w_err[i]*Gp;    
      Celle_err[j] += pref*pref*fp*fp;     //square for error propagation	  
    }
  }
  for (size_t i=0; i<Celle_err.size(); i++) {
    Celle_err[i] = sqrt(Celle_err[i]);
  }
  return Celle_err;
}


// ============================================================================================


void cbl::measure::angularpk::PowerSpectrum_angular::m_convert_angular_units (cbl::CoordinateUnits inputUnits)
{
  if (inputUnits==cbl::CoordinateUnits::_degrees_) {
    const double fact = 1./180.*M_PI;
    for(size_t i=0; i<m_theta.size(); i++){
      m_theta[i] *= fact;
    }
  }
  if (inputUnits==cbl::CoordinateUnits::_arcminutes_) {
    const double fact = 1./60./180.*M_PI;
    for(size_t i=0; i<m_theta.size(); i++){
      m_theta[i] *= fact;
    }
  }
  if (inputUnits==cbl::CoordinateUnits::_arcseconds_) {
    const double fact = 1./60./60./180.*M_PI;
    for (size_t i=0; i<m_theta.size(); i++)
      m_theta[i] *= fact;
  }
}


// ============================================================================================


void cbl::measure::angularpk::PowerSpectrum_angular::measure (const Estimator estimator, const string dir_correlation_input, const string file_correlation_input, const int n_lines_header, CoordinateUnits inputUnits, BinType input_binType, cbl::measure::ErrorType errorType, const string dir_correlation_output, const string file_correlation_output, const int nMocks)
{ 
  if (estimator==Estimator::_Fast_) {
    if (file_correlation_input=="") {  //if file is not specified the program measures the angular correlation function
      coutCBL<<"I'm measuring the Angular Two Point Correlation"<<endl;
      m_TwoPointCorrelation1D_angular->measure(errorType, dir_correlation_output, {}, par::defaultString, nMocks);
      m_TwoPointCorrelation1D_angular->write(dir_correlation_output, file_correlation_output);
      inputUnits = m_angularUnits;      //unità scritte dalle cbl possono essere diverse da quelle di input da file
      input_binType = m_binType;        //bin del file prodotto dalle cbl può essere diverso da quello di input da file

      coutCBL<<"I'm reading the file: "<<dir_correlation_output+file_correlation_output<<endl;
      const cbl::data::Data1D AngularCorrelation_data(dir_correlation_output+file_correlation_output, n_lines_header,{1},{2},{3});
      AngularCorrelation_data.Print(); cout << endl;
      m_theta = AngularCorrelation_data.xx();
      AngularCorrelation_data.get_data(m_w);
      AngularCorrelation_data.get_error(m_w_err);    
    }
    else {
      coutCBL<<"I'm reading the file: "<<dir_correlation_input+file_correlation_input<<endl;
      const cbl::data::Data1D AngularCorrelation_data(dir_correlation_input+file_correlation_input, n_lines_header,{1},{2},{3});
      AngularCorrelation_data.Print(); cout << endl;
      m_theta = AngularCorrelation_data.xx();
      AngularCorrelation_data.get_data(m_w);
      AngularCorrelation_data.get_error(m_w_err);
    }      
  }
 
  else if (estimator==Estimator::_SphericalArmonic_) {
    ErrorCBL("SphericalArmonic_Estimator not still implemented!","cbl::measure::PowerSpectrum_angular::measure", "Celle.cpp" );
  }

   else ErrorCBL("Invalid estimator!","cbl::measure::PowerSpectrum_angular::measure", "Celle.cpp" );
  

  m_convert_angular_units(inputUnits);    //if read from file, converts the theta units to radians
  std::vector<double> Celle(m_ell.size(), 0.0),Celle_err(m_ell.size(), 0.0); 
  double dlogell = log10(m_ell[1])-log10(m_ell[0]);     //logarithmic bin amplitude
  double p0 = 1/dlogell;     
  double d, pref, Gp, fp;
  for (size_t j=0; j<Celle.size(); j++){   
    for (size_t i=0; i<m_w.size();i++){   
      d = m_dtheta_theta(input_binType, i);
      pref = p0*d/(m_theta[i]*m_theta[i]);
      Gp = m_theta[i]*m_ell[j]*sqrt(std::exp(dlogell))*cyl_bessel_j(1,m_theta[i]*m_ell[j]*sqrt(std::exp(dlogell)))-m_theta[i]*m_ell[j]/sqrt(std::exp(dlogell))*cyl_bessel_j(1,m_theta[i]*m_ell[j]/sqrt(std::exp(dlogell)));
      fp = m_w[i] * Gp;
      Celle[j] += pref*fp;
    }
  }
  m_Celle.insert(std::end(m_Celle), std::begin(Celle), std::end(Celle));
  Celle_err=m_error(input_binType,m_w_err);
  m_Celle.insert(std::end(m_Celle), std::begin(Celle_err), std::end(Celle_err));
}


// ============================================================================================
