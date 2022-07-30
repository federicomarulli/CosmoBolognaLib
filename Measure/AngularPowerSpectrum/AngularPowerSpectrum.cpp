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
 *  @file Measure/AngularPowerSpectrum/AngularPowerSpectrum.cpp
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
#include "Catalogue.h"
#include "Data1D.h"
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include "LegendrePolynomials.h"

using namespace std;
using namespace cbl;   
using namespace measure;
using namespace twopt;
using namespace boost::math;


// ============================================================================================


void cbl::measure::angularpk::PowerSpectrum_angular::m_read_angular_mask(const string mask_file, std::string mask_type, double pixel_area, int n_lines_header){

  double ra_mask, theta_mask;  //RA and colatitute of the pixel center of the mask
  
  if(mask_type=="Healpix"){  //binary mask produced with pixelfunc.pix2ang in healpix. Colatitude and RA of the centers of the occupied pixel
    m_pixel_area=pixel_area;
    ifstream fin(mask_file.c_str()); checkIO(fin, mask_file);
    string line;     
    // skip the first nlines (in case of header lines)
    if (n_lines_header>0)
      for (int i=0; i<n_lines_header; ++i)
	getline(fin, line); 
    
    // read the file lines
    while (getline(fin, line)) {
      
      stringstream ss(line);
      ss >> theta_mask; 
      ss >> ra_mask; 
      m_theta_mask.emplace_back(theta_mask);
      m_RA_mask.emplace_back(ra_mask);
      
    }
    fin.close();
    
    m_survey_area=m_theta_mask.size()*pixel_area;     //measure total area of the survey
        
  }else{     //mask obtained from the fits file. Colatitude and RA of the pixel centers and edges
      
    double ra_mask_min, ra_mask_max;
    double theta_mask_min, theta_mask_max, pixel_area_file;
    std::vector<double> area_pixel; //vector of pixel areas
    ifstream fin(mask_file.c_str()); checkIO(fin, mask_file);
    
    string line;     
    // skip the first nlines (in case of header lines)
    if (n_lines_header>0)
      for (int i=0; i<n_lines_header; ++i)
	getline(fin, line);
    
    // read the file lines
    while (getline(fin, line)) {
      
      stringstream ss(line);
      ss >> theta_mask; 
      ss >> theta_mask_min; 
      ss >> theta_mask_max; 
      ss >> ra_mask; 
      ss >> ra_mask_min; 
      ss >> ra_mask_max; 
      ss >> pixel_area_file;
      
      m_theta_mask.emplace_back(theta_mask);
      m_theta_mask_min.emplace_back(theta_mask_min);
      m_theta_mask_max.emplace_back(theta_mask_max);
      m_RA_mask.emplace_back(ra_mask);
      m_RA_mask_min.emplace_back(ra_mask_min);
      m_RA_mask_max.emplace_back(ra_mask_max);
      area_pixel.emplace_back(pixel_area_file);
      
    }
    fin.close();
    m_survey_area=0;             //measure total area of the survey
    for(size_t i=0; i<m_theta_mask.size(); i++) {
      m_survey_area+=area_pixel[i];
    }

    double area_overlap=0;     //remove overlap between pixels
    //pixel_area=(m_RA_mask_max[0]-m_RA_mask_min[0])*(m_theta_mask_max[0]-m_theta_mask_min[0]);  //effective pixel area, not normalized with colatitude
    std::vector<double> theta_side, RA_side;
    for (size_t i=0; i<m_theta_mask.size(); ++i){  //in case of rectangular and slightly different pixels
      theta_side.emplace_back(abs(m_theta_mask_max[i]-m_theta_mask_min[i]));
      RA_side.emplace_back(abs(m_RA_mask_max[i]-m_RA_mask_min[i]));
    }
    for(size_t i=0; i<m_theta_mask.size(); i++) {
      for(size_t j=i; j<m_theta_mask.size(); j++){
	if(i!=j){
	  if(abs(m_theta_mask[i]-m_theta_mask[j])<cbl::Average(theta_side) && abs(m_RA_mask[i]-m_RA_mask[j])<cbl::Average(RA_side)){
	    area_overlap+=(cbl::Average(theta_side)-abs(m_theta_mask[i]-m_theta_mask[j]))*(cbl::Average(RA_side)-abs(m_RA_mask[i]-m_RA_mask[j]))*sin((m_theta_mask[i]+m_theta_mask[j])/2);
	  }
	}
      }
    }
    coutCBL<<"Area overlap= "<<area_overlap<<endl;
    m_survey_area-=area_overlap;    
  }
  coutCBL<<"Survey area= "<<m_survey_area<<endl;

}


// ============================================================================================


double cbl::measure::angularpk::PowerSpectrum_angular::angular_mask(double theta, double RA){

  if(m_theta_mask.size()!=0){
    bool square=true;
    if(m_theta_mask_min.size()==0 || m_theta_mask_max.size()==0 || m_RA_mask_min.size()==0 || m_RA_mask_max.size()==0)  square=false;
    double masked=0.;
    if(!square){
      double radius=sqrt(m_pixel_area/M_PI);
      for(size_t i=0; i<m_theta_mask.size(); ++i)
	if(pow(theta-m_theta_mask[i], 2)+pow(RA-m_RA_mask[i], 2)<radius*radius) {
	  masked=1.;
	  break;
	}
    }
    else{
      for(size_t i=0; i<m_theta_mask.size(); ++i)
	if(m_theta_mask_min[i]<theta && m_theta_mask_max[i]>theta)
	  if(m_RA_mask_min[i]<RA && m_RA_mask_max[i]>RA){
	    masked=1.;
	    break;
	  }
    }
    return masked;
  }
  else return 1;
}


//====================================================================================================


cbl::measure::angularpk::PowerSpectrum_angular::PowerSpectrum_angular (const catalogue::Catalogue data, const catalogue::Catalogue random, const double l_min, const double l_max, const int Nl, const BinType binType, const double thetaMin, const double thetaMax, const int nbins, const double shift, const CoordinateUnits angularUnits, const BinType correlation_binType)
{
  m_l_min = l_min;   
  m_l_max = l_max;
  
  if(int(Nl)!=Nl || Nl<2)  ErrorCBL("In spherical armonic estimator the Nl should be an >=2!","cbl::measure::PowerSpectrum_angular::PowerSpectrum_angular", "AngularPowerSpectrum.cpp" );
  else m_Nl=Nl;
  m_binType = correlation_binType;
  m_angularUnits = angularUnits;
  m_nbins=nbins;
  set_l_output(binType);
  cbl::measure::twopt::TwoPointCorrelation1D_angular TwoPointCorrelation1D_angular(data, random, correlation_binType, thetaMin, thetaMax, nbins, shift, angularUnits);
  m_TwoPointCorrelation1D_angular = std::make_shared<twopt::TwoPointCorrelation1D_angular>(twopt::TwoPointCorrelation1D_angular(std::move(TwoPointCorrelation1D_angular)));
}


// ============================================================================================


cbl::measure::angularpk::PowerSpectrum_angular::PowerSpectrum_angular (const catalogue::Catalogue data, const double l_min, const double l_max, const int bandwidth, const string mask_file, const string mask_type, double pixel_area, int n_lines_header)
{
  m_l_min = l_min;   
  m_l_max = l_max;
  
  if(int(bandwidth)!=bandwidth || int(bandwidth)<1)  ErrorCBL("In spherical armonic estimator the bandwidth (l_max-l_min)/Nl should be an integer>=2!","cbl::measure::PowerSpectrum_angular::PowerSpectrum_angular", "AngularPowerSpectrum.cpp" );
  else if(int(bandwidth)==1){
  }else m_Nl = (l_max-l_min)/bandwidth;
  set_catalogue(data);
  
  if(mask_file!="") m_read_angular_mask(mask_file, mask_type, pixel_area, n_lines_header);
  else {
    double pi_2=cbl::par::pi*0.5;
    m_survey_area=(cos(pi_2-m_catalogue->Max(cbl::catalogue::Var::_Dec_))-cos(pi_2-m_catalogue->Min(cbl::catalogue::Var::_Dec_)))*(m_catalogue->Max(cbl::catalogue::Var::_RA_)-m_catalogue->Min(cbl::catalogue::Var::_RA_));
  }
}


// ============================================================================================


void cbl::measure::angularpk::PowerSpectrum_angular::set_catalogue (cbl::catalogue::Catalogue catalogue)
{
  m_catalogue = std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(catalogue)));
  if(!m_catalogue->isSetVar(cbl::catalogue::Var::_RA_) || !m_catalogue->isSetVar(cbl::catalogue::Var::_Dec_) || !m_catalogue->isSetVar(cbl::catalogue::Var::_Dc_)) ErrorCBL("The catalogue must be in observed coordinates!","cbl::measure::PowerSpectrum_angular::set_catalogue", "AngularPowerSpectrum.cpp" );

}
  

//==================================================================================


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


void cbl::measure::angularpk::PowerSpectrum_angular::set_l_output (const BinType binType)
{
  if (m_l_min==0) m_l_min=1;
  double size_bin;
  if (binType==cbl::BinType::_linear_) {    //multipole scale is linear 
    size_bin = (m_l_max-m_l_min)/m_Nl;   //if linear binning
    for(int i=0; i<m_Nl; i++){
      m_l.emplace_back(m_l_min+size_bin*i);
    }
  }
  
  if (binType==cbl::BinType::_logarithmic_) {    //multipole scale is logarithmic
    const double l_min_log = std::log10(m_l_min);
    const double l_max_log = std::log10(m_l_max);
    const double log_increment = (l_max_log - l_min_log)/ m_Nl;

    for (int i=0; i<m_Nl+1; ++i) {
      m_l.emplace_back(std::pow(10, l_min_log+i*log_increment));
    }
  }
}


// ============================================================================================


void cbl::measure::angularpk::PowerSpectrum_angular::write (const std::string dir, const std::string file)
{

  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {}
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);
     
  fout << "#l \t AngularPowerSpectrum \t AngularPowerSpectrum_error"<<endl;

  for (size_t j=0; j<m_AngularPowerSpectrum.size()/2; j++) {
    fout<<setiosflags(ios::fixed) << setprecision(5) << setw(10) << right<<m_l[j]<<
      "   " << setiosflags(ios::fixed) << setprecision(9) << setw(15) << right <<m_AngularPowerSpectrum[j] <<
      "   " << setiosflags(ios::fixed) << setprecision(9) << setw(15) << right <<m_AngularPowerSpectrum[m_AngularPowerSpectrum.size()/2+j]<<endl;
  }
  
  fout.clear(); fout.close(); coutCBL << "I wrote the file " << file_out << endl;

}


// ============================================================================================


void cbl::measure::angularpk::PowerSpectrum_angular::write_average (const std::string dir, const std::string file)
{
  if(isSet(m_Nl)){
    string MK = "mkdir -p "+dir; if (system (MK.c_str())) {}
  
    string file_out = dir+file;
    ofstream fout(file_out.c_str()); checkIO(fout, file_out);
    
    fout << "#l \t AngularPowerSpectrum average \t AngularPowerSpectrum_error average"<<endl;
    
    for (size_t j=0; j<m_AngularPowerSpectrum_av.size()/2; j++) {
      fout<<setiosflags(ios::fixed) << setprecision(5) << setw(10) << right<<m_l_av[j]<<
	"   " << setiosflags(ios::fixed) << setprecision(9) << setw(15) << right <<m_AngularPowerSpectrum_av[j] <<
	"   " << setiosflags(ios::fixed) << setprecision(9) << setw(15) << right <<m_AngularPowerSpectrum_av[m_AngularPowerSpectrum_av.size()/2+j]<<endl;
    }
    
    fout.clear(); fout.close(); coutCBL << "I wrote the file " << file_out << endl;
  }else{
    ErrorCBL("m_Nl is not set! Specify an integer bandwidth>=2 and use measure with SphericalArmonic estimator!","cbl::measure::PowerSpectrum_angular::measure", "AngularPowerSpectrum.cpp" );
  }
}


// ===========================================================================================


void  cbl::measure::angularpk::PowerSpectrum_angular::write_mixing_matrix(const std::string dir, const std::string file, bool store_window) {

  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {}
 
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);
 
  fout << "#l \t R_ll' "<<endl;
  for(size_t i=0; i<m_l.size(); ++i){
    fout<<setiosflags(ios::fixed) << setprecision(5) << setw(10) << right<<m_l[i]<<"   ";
    for(size_t j=0; j<m_l.size(); ++j){
      fout<<setiosflags(ios::fixed) << setprecision(9) << setw(10) << right<<m_RR[i][j]<<"   ";
    }
    fout<<setiosflags(ios::fixed) << setprecision(9) << setw(10) << right<<endl;
  }
  fout.clear(); fout.close(); coutCBL << "I wrote the file " << file_out << endl;

  if(store_window){
    string file_out = dir+"window.dat";
    ofstream fout(file_out.c_str()); checkIO(fout, file_out);
    
    fout << "#l \t W_l"<<endl;
    for(size_t i=0; i<m_l.size(); ++i){
      fout<<setiosflags(ios::fixed) << setprecision(5) << setw(10) << right<<m_l[i]<<"   " << setiosflags(ios::fixed) << setprecision(9) << setw(10) << right <<m_Wl[i]<<endl;
    }
    fout.clear(); fout.close(); coutCBL << "I wrote the file " << file_out << endl;
  }
}


// ===========================================================================================


void cbl::measure::angularpk::PowerSpectrum_angular::read (const std::string dir, const std::string file, const std::vector<int> column_x, const std::vector<int> column_y, const std::vector<int> column_error, const int n_lines_header)
{
  coutCBL<<"I'm reading the file: "<<dir+file<<endl;
  const cbl::data::Data1D data(dir+file, n_lines_header,column_x, column_y, column_error);
  m_x = data.xx();
  data.get_data(m_y);
  data.get_error(m_y_err);

}


// ============================================================================================


void cbl::measure::angularpk::read_mixing_matrix (const std::string dir_matrix, const std::string file_matrix, std::vector<double> &ll, std::vector<std::vector<double>> &matrix)
{
  std::ifstream in;
  in.open(dir_matrix+file_matrix);
  std::string line;
  std::vector<double> row;
  if(in.is_open())
    {
      int skipped_lines=1;
      for (int i=0; i<skipped_lines; ++i) getline(in, line);
      while(std::getline(in, line))                     //get 1 row as a string
        {
	  std::istringstream iss(line);                 //put line into stringstream
	  double element;
	  int i=0;
	  while(iss >> element)                         //read double by double
	    {
	      if (i==0) {i++;                           //skip first column (multipole)
		ll.emplace_back(element);
		continue;
	      }
	      row.emplace_back(element);
	      i++;
	    }
	  matrix.push_back(row);
	  row.erase(row.begin(), row.end());
	}
    }
}


// ======================================================================================================


std::vector<double> cbl::measure::angularpk::PowerSpectrum_angular::m_error (const BinType binType, std::vector<double> w_err)
{
  std::vector<double> AngularPowerSpectrum_err(m_l.size(), 0.0); 
  double dlogl = log(m_l[1])-log(m_l[0]);    //logarithmic bin amplitude
  double p0 = 2*cbl::par::pi/dlogl;     
  double d, pref, Gp, fp;
  for (size_t j=0; j<AngularPowerSpectrum_err.size(); j++) {    
    for (size_t i=0; i<w_err.size();i++) {   
      d = m_dtheta_theta(binType, i);
      pref = p0*d/(m_theta[i]*m_theta[i]);
      Gp = m_theta[i]*m_l[j]*sqrt(std::exp(dlogl))*cyl_bessel_j(1,m_theta[i]*m_l[j]*sqrt(std::exp(dlogl)))-m_theta[i]*m_l[j]/sqrt(std::exp(dlogl))*cyl_bessel_j(1,m_theta[i]*m_l[j]/sqrt(std::exp(dlogl)));
      fp = m_w_err[i]*Gp;    
      AngularPowerSpectrum_err[j] += pref*pref*fp*fp;     //square for error propagation	  
    }
  }
  for (size_t i=0; i<AngularPowerSpectrum_err.size(); i++) {
    AngularPowerSpectrum_err[i] = sqrt(AngularPowerSpectrum_err[i])/m_l[i]/m_l[i];
  }
  return AngularPowerSpectrum_err;
}


// ============================================================================================


void cbl::measure::angularpk::PowerSpectrum_angular::m_convert_angular_units (cbl::CoordinateUnits inputUnits)
{
  if (inputUnits==cbl::CoordinateUnits::_radians_) {}
  
  if (inputUnits==cbl::CoordinateUnits::_degrees_) {
    const double fact = 1./180.*cbl::par::pi;
    for(size_t i=0; i<m_theta.size(); i++){
      m_theta[i] *= fact;
    }
  }
  if (inputUnits==cbl::CoordinateUnits::_arcminutes_) {
    const double fact = 1./60./180.*cbl::par::pi;
    for(size_t i=0; i<m_theta.size(); i++){
      m_theta[i] *= fact;
    }
  }
  if (inputUnits==cbl::CoordinateUnits::_arcseconds_) {
    const double fact = 1./60./60./180.*cbl::par::pi;
    for (size_t i=0; i<m_theta.size(); i++)
      m_theta[i] *= fact;
  }
}

// ============================================================================================


void cbl::measure::angularpk::PowerSpectrum_angular::measure (const AngularEstimator estimator, const string dir_correlation_input, const string file_correlation_input, const int n_lines_header, CoordinateUnits inputUnits, BinType input_binType, cbl::measure::ErrorType errorType, const string dir_correlation_output, const string file_correlation_output, const int nMocks)
{

  switch (estimator) {
  case(AngularEstimator::_Fast_):
     measureFast(dir_correlation_input, file_correlation_input, n_lines_header, inputUnits, input_binType, errorType, dir_correlation_output, file_correlation_output, nMocks);
     break;
   case(AngularEstimator::_SphericalArmonic_):
     measureSphericalArmonic();
     break;
   default:
     ErrorCBL("Invalid AngularEstimator!","cbl::measure::PowerSpectrum_angular::measure", "AngularPowerSpectrum.cpp" );
  }
  
}


//=====================================================================================================


void cbl::measure::angularpk::PowerSpectrum_angular::measureFast (const string dir_correlation_input, const string file_correlation_input, const int n_lines_header, CoordinateUnits inputUnits, BinType input_binType, cbl::measure::ErrorType errorType, const string dir_correlation_output, const string file_correlation_output, const int nMocks)
{
  m_AngularPowerSpectrum.erase(std::begin(m_AngularPowerSpectrum), std::end(m_AngularPowerSpectrum));

  
  if (file_correlation_input=="") {  //if file is not specified the program measures the angular correlation function
    coutCBL<<"I'm measuring the Angular Two Point Correlation"<<endl;
    m_TwoPointCorrelation1D_angular->measure(errorType, dir_correlation_output, {}, par::defaultString, nMocks);
    m_TwoPointCorrelation1D_angular->write(dir_correlation_output, file_correlation_output); 
    inputUnits = m_angularUnits;     
    input_binType = m_binType;        
    read(dir_correlation_output, file_correlation_output, {1}, {2}, {3}, n_lines_header);
    m_theta=m_x;
    m_w=m_y;
    m_w_err=m_y_err;
    
  }
  else {
    
    read(dir_correlation_input, file_correlation_input, {1}, {2}, {3}, n_lines_header);
    m_theta=m_x;
    m_w=m_y;
    m_w_err=m_y_err;
    
  }
  
  m_convert_angular_units(inputUnits);    //if read from file, converts the theta units to radians
  std::vector<double> AngularPowerSpectrum(m_l.size(), 0.0),AngularPowerSpectrum_err(m_l.size(), 0.0); 
  double dlogl = log(m_l[1])-log(m_l[0]);     //logarithmic bin amplitude
  double p0 = 2*cbl::par::pi/dlogl;     
  double d, pref, Gp, fp;
  for (size_t j=0; j<AngularPowerSpectrum.size(); j++){   
    for (size_t i=0; i<m_w.size();i++){   
      d = m_dtheta_theta(input_binType, i);
      pref = p0*d/(m_theta[i]*m_theta[i]);
      Gp = m_theta[i]*m_l[j]*sqrt(std::exp(dlogl))*cyl_bessel_j(1,m_theta[i]*m_l[j]*sqrt(std::exp(dlogl)))-m_theta[i]*m_l[j]/sqrt(std::exp(dlogl))*cyl_bessel_j(1,m_theta[i]*m_l[j]/sqrt(std::exp(dlogl)));
      fp = m_w[i] * Gp;
      AngularPowerSpectrum[j] += pref*fp/m_l[j]/m_l[j];
    }
  }
  m_AngularPowerSpectrum.insert(std::end(m_AngularPowerSpectrum), std::begin(AngularPowerSpectrum), std::end(AngularPowerSpectrum));
  AngularPowerSpectrum_err=m_error(input_binType,m_w_err);
  m_AngularPowerSpectrum.insert(std::end(m_AngularPowerSpectrum), std::begin(AngularPowerSpectrum_err), std::end(AngularPowerSpectrum_err));
  
}


//==================================================================================


void cbl::measure::angularpk::PowerSpectrum_angular::measureSphericalArmonic ()
{
  double pi_2=cbl::par::pi*0.5;
  double _pi_4=1./cbl::par::pi*0.25;
  m_l.erase(std::begin(m_l), std::end(m_l));
  m_AngularPowerSpectrum.erase(std::begin(m_AngularPowerSpectrum), std::end(m_AngularPowerSpectrum));
  
  double survey_area=m_survey_area;
  double sigma0=m_catalogue->nObjects()/survey_area;
  double fsky=survey_area*_pi_4;    //fraction of the sky covered by the survey
  coutCBL<<"fsky= "<<fsky<<endl;
  std::vector<double>  W_ll, AngularPowerSpectrum_err;
  complex<double> Y_lm;
  cout<<endl;
  coutCBL<<" l \t C_l "<<endl;
  for(int l=m_l_min; l<m_l_max; l++){ //works only with l<=879
    double A_lm, C_lm=0, C_l=0, W_l=0;
    for(int m=-l; m<l+1; m++){
      Y_lm=0;
      for(size_t i=0; i<m_catalogue->nObjects(); i++){
	Y_lm+=conj(cbl::spherical_harmonics(cbl::CoordinateType::_observed_,l, m, m_catalogue->dec(i), m_catalogue->ra(i), 1.));
      }
      if(fsky>0.9){    //total survey
	A_lm=abs(Y_lm);
	C_lm=A_lm*A_lm-m_catalogue->nObjects()*_pi_4;
	C_l+=C_lm;
      }
      else{
	const int ndim = 2;
	auto integrand_real = [this, &l, &m] (std::vector<double> x) { return angular_mask(x[0], x[1])*pow(-1,m)*boost::math::spherical_harmonic(l, m, x[0], x[1]).real()*sin(x[0]); };
	auto integrand_im = [this, &l, &m] (std::vector<double> x) { return angular_mask(x[0], x[1])*pow(-1,m)*boost::math::spherical_harmonic(l, m, x[0], x[1]).imag()*sin(x[0]); };
	auto integrand_J_lm = [this, &l, &m] (std::vector<double> x) { return angular_mask(x[0], x[1])*pow(abs(boost::math::spherical_harmonic(l, m, x[0], x[1])),2)*sin(x[0]); };
	
	std::vector<std::vector<double>> integration_limits(2);
	integration_limits[0] = {pi_2-m_catalogue->Max(cbl::catalogue::Var::_Dec_), pi_2-m_catalogue->Min(cbl::catalogue::Var::_Dec_)};
	integration_limits[1] = {m_catalogue->Min(cbl::catalogue::Var::_RA_), m_catalogue->Max(cbl::catalogue::Var::_RA_)};
	
	// wrapper to CUBA libraries
	cbl::wrapper::cuba::CUBAwrapper CW_real(integrand_real, ndim);
	double I_lm_re=CW_real.IntegrateCuhre(integration_limits);
	cbl::wrapper::cuba::CUBAwrapper CW_im(integrand_im, ndim);
	double I_lm_im=CW_im.IntegrateCuhre(integration_limits);
	complex<double> I_lm (I_lm_re, -1.*I_lm_im);//complex conjugate
	cbl::wrapper::cuba::CUBAwrapper CW_J_lm(integrand_J_lm, ndim);
	double J_lm=CW_J_lm.IntegrateCuhre(integration_limits);
	C_lm=abs(Y_lm-m_catalogue->nObjects()/survey_area*I_lm);
	C_lm=C_lm*C_lm/J_lm-m_catalogue->nObjects()/survey_area;
	C_l+=C_lm;
	W_l+=pow(abs(I_lm),2);
      }
    }
    C_l=C_l/(2*l+1);
    C_l/=pow(sigma0,2);
    W_l=W_l/(2*l+1);
    m_l.emplace_back(l);
    m_Wl.emplace_back(W_l);
    m_AngularPowerSpectrum.emplace_back(C_l);
    AngularPowerSpectrum_err.emplace_back(sqrt(2/fsky/(2*l+1))*(C_l+1/sigma0));
    coutCBL<<l<<"   "<<C_l<<endl;

  }
  
  m_AngularPowerSpectrum.insert(std::end(m_AngularPowerSpectrum), std::begin(AngularPowerSpectrum_err), std::end(AngularPowerSpectrum_err));
  if(isSet(m_Nl)) m_averagePowerSpectrum();
  // if (fsky<0.9) compute_mixing_matrix();

    
}


// ============================================================================================


void cbl::measure::angularpk::PowerSpectrum_angular::m_averagePowerSpectrum()
{
  double size=(m_l_max-m_l_min)/m_Nl;
  std::vector<double>  AngularPowerSpectrum_error_av(m_Nl, 0);

  for(int i=0; i<m_Nl; i++){
    m_l_av.emplace_back(0);
    m_AngularPowerSpectrum_av.emplace_back(0);
  }
  
  for(size_t j=0; j<m_Nl; j++){
    for(int i=0; i<size; i++){
      if(i+j*size>=m_l.size()) break;
      m_l_av[j]+=m_l[j*size+i];
      m_AngularPowerSpectrum_av[j]+=m_AngularPowerSpectrum[j*size+i];
      AngularPowerSpectrum_error_av[j]+=m_AngularPowerSpectrum[m_AngularPowerSpectrum.size()/2+j*size+i];
    }
    m_l_av[j]/=size;
    m_AngularPowerSpectrum_av[j]/=size;
    AngularPowerSpectrum_error_av[j]/=size;
  }
  m_AngularPowerSpectrum_av.insert(std::end(m_AngularPowerSpectrum_av), std::begin(AngularPowerSpectrum_error_av), std::end(AngularPowerSpectrum_error_av));
}


// ============================================================================================


void cbl::measure::angularpk::PowerSpectrum_angular::compute_mixing_matrix(string dir_window_input, string file_window_input)
{
  double _pi_4=1./cbl::par::pi*0.25;
  if(dir_window_input!="" && file_window_input!="") {
    read(dir_window_input, file_window_input, {1}, {2}, {2}, 1);
    m_l=m_x;
    m_Wl=m_y;
  }
  
  std::vector<double> row;
  double R;
  for(size_t k=0; k<m_l.size(); ++k) { // l cycle
    for(size_t j=0; j<m_l.size(); ++j) { // l' cycle
      R=0;
      for(size_t i=0; i<m_l.size(); ++i)  { // l'' cycle
	R+=(2*m_l[i]+1)*m_Wl[i]*cbl::wigner3j(m_l[k], m_l[j], m_l[i], 0, 0, 0)*cbl::wigner3j(m_l[k], m_l[j], m_l[i], 0, 0, 0); 
      }
      R*=(2*m_l[j]+1)*_pi_4;
      row.emplace_back(R);
    }
    m_RR.emplace_back(row);
    row.erase(row.begin(), row.end());
  }
}
