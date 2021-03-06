/********************************************************************
 *  Copyright (C) 2020 by Giorgio Lesci and Federico Marulli        *
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
 *  @file StackedDensityProfile/StackedDensityProfile.cpp
 *
 *  @brief Methods of the class StackedDensityProfile 
 *
 *  This file contains the implementation of the methods of the class
 *  StackedDensityProfile, used to compute the stacked 
 *  density profile of galaxy clusters from lensing data
 *
 *  @author Giorgio Lesci (and Fabio Bellagamba, Federico Marulli)
 *
 *  @author giorgio.lesci2@unibo.it (and fabiobg83@gmail.com, federico.marulli3@unibo.it)
 */

#include "StackedDensityProfile.h"
#include "CCfits/CCfits"
#include "Data1D_extra.h"
#include<time.h>

using namespace cbl;
using namespace measure::stackprofile;


// ============================================================================

void cbl::measure::stackprofile::StackedDensityProfile::m_set_multiplicative_calib ()
{
  if (m_m_calib.size() > 0) {

    if (m_m_calib.size() != m_m_calib_zEdges.size())
      cbl::ErrorCBL("multiplicative_calibration_stats and multiplicative_calibration_zEdges must have the same size!","m_set_multiplicative_calib","StackedDensityProfile.cpp");
    
    for (size_t i=0; i<m_m_calib.size(); i++) {
      if (m_m_calib[i].size() != 2)
	cbl::ErrorCBL("All the vectors in multiplicative_calibration_stats must have size = 2","m_set_multiplicative_calib","StackedDensityProfile.cpp");
      if (m_m_calib_zEdges[i].size() != 2)
	cbl::ErrorCBL("All the vectors in multiplicative_calibration_zEdges must have size = 2","m_set_multiplicative_calib","StackedDensityProfile.cpp");
    }

    coutCBL<<"I'm setting the m parameter for each galaxy by extracting values from Gaussian distributions..."<<std::endl;
    std::vector<double> m_values((int)(m_galData->nObjects()));
    for (size_t j=0; j<m_galData->nObjects(); j++) {
      for (size_t k=0; k<m_m_calib_zEdges.size(); k++) {
	if (m_galData->redshift(j) >= m_m_calib_zEdges[k][0] && m_galData->redshift(j) <= m_m_calib_zEdges[k][1]) {
	  
	  srand(time(0));
	  cbl::random::NormalRandomNumbers random(m_m_calib[k][0], m_m_calib[k][1], rand());
	  m_values[j] = random();
	    
	}
      }
    }
    
    m_galData->set_var(cbl::catalogue::Var::_LensingCalib_, m_values);
    
  } else {
    if (m_galData->isSetVar(cbl::catalogue::Var::_LensingCalib_) != true)
      cbl::ErrorCBL("You must set the galaxy catalogue variable LensingCalib","m_set_multiplicative_calib","StackedDensityProfile.cpp");
  }
}
  

// ============================================================================

bool cbl::measure::stackprofile::StackedDensityProfile::m_colourSelection_Oguri (const int i_gal)
{
  bool is_good_gri = ((m_galData->magnitudeG(i_gal)-m_galData->magnitudeR(i_gal) < 0.3 || m_galData->magnitudeR(i_gal)-m_galData->magnitudeI(i_gal) > 1.3 || m_galData->magnitudeR(i_gal)-m_galData->magnitudeI(i_gal) > m_galData->magnitudeG(i_gal)-m_galData->magnitudeR(i_gal)) && m_galData->redshift(i_gal)>=0.2);
  return is_good_gri;
}

// ============================================================================

bool cbl::measure::stackprofile::StackedDensityProfile::m_zphotSelection_ODDS (const std::vector<int> x)
{
  // x[0] is the cluster index, while x[1] is the galaxy index
  bool is_good_photz = (m_galData->redshiftMin(x[1]) > m_cluData->redshift(x[0])+m_zphot_sel_pars[0] && m_galData->redshift(x[1])>=m_zphot_sel_pars[1] && m_galData->redshift(x[1])<=m_zphot_sel_pars[2] && m_galData->odds(x[1]) > m_zphot_sel_pars[3]);
  return is_good_photz;
}

// ============================================================================

bool cbl::measure::stackprofile::StackedDensityProfile::m_zphotSelection_noODDS (const std::vector<int> x)
{
  // x[0] is the cluster index, while x[1] is the galaxy index
  bool is_good_photz = (m_galData->redshiftMin(x[1]) > m_cluData->redshift(x[0])+m_zphot_sel_pars[0] && m_galData->redshift(x[1])>=m_zphot_sel_pars[1] && m_galData->redshift(x[1])<=m_zphot_sel_pars[2]);
  return is_good_photz;
}

// ============================================================================

void cbl::measure::stackprofile::StackedDensityProfile::m_set_logicSelection (const std::string sel)
{
  std::function<bool(const std::vector<bool>)> returned_fc;

  std::function<bool(const std::vector<bool>)> logicSelection_and = [] (const std::vector<bool> x){
								      // x[0]->colour selection, x[1]->photz selection
								      bool thebool = (x[0] && x[1]);
								      return thebool;
								    };

  std::function<bool(const std::vector<bool>)> logicSelection_or = [] (const std::vector<bool> x){
								     // x[0]->colour selection, x[1]->photz selection
								     bool thebool = (x[0] || x[1]);
								     return thebool;
								   };
  if (sel=="and"){
    returned_fc = logicSelection_and;
  }else if (sel=="or"){
    returned_fc = logicSelection_or;
  }else{
    cbl::ErrorCBL("The chosen logic operator between the selection criteria is not valid! Choose \'and\' or \'or\'.","StackedDensityProfile","StackedDensityProfile.cpp");
  }
  m_logicSelection = returned_fc;
}

// ============================================================================

void cbl::measure::stackprofile::StackedDensityProfile::m_set_selections ()
{    
  if (m_colour_sel=="Oguri"){
    m_colourSel = &cbl::measure::stackprofile::StackedDensityProfile::m_colourSelection_Oguri;
  }else{
    cbl::ErrorCBL("The chosen colour selection is not valid!","StackedDensityProfile","StackedDensityProfile.cpp");
  }
  if (m_zphot_sel=="Odds"){
    if (m_galData->isSetVar(cbl::catalogue::Var::_ODDS_) != true)
      cbl::ErrorCBL("You must set the ODDS galaxy catalogue variable!","m_set_selections","StackedDensityProfile.cpp");
    else
      m_photzSel = &cbl::measure::stackprofile::StackedDensityProfile::m_zphotSelection_ODDS;
  }else if (m_zphot_sel=="noOdds"){
    m_photzSel = &cbl::measure::stackprofile::StackedDensityProfile::m_zphotSelection_noODDS;
  }else{
    cbl::ErrorCBL("The chosen photo-z selection is not valid!","StackedDensityProfile","StackedDensityProfile.cpp");
  }    
  m_set_logicSelection(m_logic_sel);
}


// ============================================================================

void cbl::measure::stackprofile::StackedDensityProfile::m_check_catalogue_variables ()
{
  if (m_galData->isSetVar(cbl::catalogue::Var::_RA_) != true || m_galData->isSetVar(cbl::catalogue::Var::_Dec_) != true || m_galData->isSetVar(cbl::catalogue::Var::_Redshift_) != true || m_galData->isSetVar(cbl::catalogue::Var::_Shear1_) != true || m_galData->isSetVar(cbl::catalogue::Var::_Shear2_) != true || m_galData->isSetVar(cbl::catalogue::Var::_RedshiftMin_) != true || m_galData->isSetVar(cbl::catalogue::Var::_RedshiftMax_) != true || m_galData->isSetVar(cbl::catalogue::Var::_LensingWeight_) != true)
    cbl::ErrorCBL("You must set the following galaxy catalogue variables: RA, Dec, redshift, redshift_min, redshift_max, shear1, shear2, lensing_weight","m_check_catalogue_variables","StackedDensityProfile.cpp");
  if (m_cluData->isSetVar(cbl::catalogue::Var::_RA_) !=true || m_cluData->isSetVar(cbl::catalogue::Var::_Dec_) !=true || m_cluData->isSetVar(cbl::catalogue::Var::_Redshift_) !=true || m_cluData->isSetVar(cbl::catalogue::Var::_SN_) !=true || m_cluData->isSetVar(cbl::catalogue::Var::_MassProxy_) !=true)
    cbl::ErrorCBL("You must set the following cluster catalogue variables: RA, Dec, redshift, S/N, Mass Proxy","m_check_catalogue_variables","StackedDensityProfile.cpp");
}


// ============================================================================

void cbl::measure::stackprofile::StackedDensityProfile::m_resize (const double rad_min, const double rad_max, const int nRad, const bool log_rad)
{
  if (m_z_binEdges.size()-1 != m_proxy_binEdges.size())
      cbl::ErrorCBL("The number of proxy binnings must be equal to the number of redshift bins!","StackedDensityProfile","StackedDensityProfile.cpp");
  
  m_rad_arr.resize(nRad+1);
  if (log_rad == false){
    for (int i=0; i < nRad+1; i++)
      m_rad_arr[i] = rad_min + (rad_max-rad_min)/double(nRad)*i;
  }else{
    for (int i = 0; i<nRad+1; i++)
      m_rad_arr[i] = exp(log(rad_min) + (log(rad_max)-log(rad_min))/double(nRad)*i);
  }
  
  m_ngal_arr.resize(m_rad_arr.size()-1,0);
  m_deltasigma_t.resize(m_rad_arr.size()-1,0.);
  m_deltasigma_x.resize(m_rad_arr.size()-1,0.);
  m_deltasigma_err.resize(m_rad_arr.size()-1,0.);
  m_wetasquareSum.resize(m_rad_arr.size()-1,0.);
  m_rad_eff_arr.resize(m_rad_arr.size()-1,0.);
  m_rad_sigma_arr.resize(m_rad_arr.size()-1,0.);
  m_wSum.resize(m_rad_arr.size()-1,0.);
  m_wetaSum.resize(m_rad_arr.size()-1,0.);
  m_deltasigmaSum.resize(m_rad_arr.size()-1,0.);
  
  m_stacked_ngal_arr.resize(m_z_binEdges.size()-1);
  m_stacked_deltasigma_t.resize(m_z_binEdges.size()-1);
  m_stacked_deltasigma_x.resize(m_z_binEdges.size()-1);
  m_stacked_deltasigma_t_err.resize(m_z_binEdges.size()-1);
  m_stacked_deltasigma_x_err.resize(m_z_binEdges.size()-1);
  m_stacked_deltasigma_err.resize(m_z_binEdges.size()-1);
  m_stacked_wetasquareSum.resize(m_z_binEdges.size()-1);
  m_stacked_wSum.resize(m_z_binEdges.size()-1);
  m_stacked_wetaSum.resize(m_z_binEdges.size()-1);
  m_stacked_deltasigmaSum.resize(m_z_binEdges.size()-1);
  m_stacked_deltasigma_wei.resize(m_z_binEdges.size()-1);
  m_proxy_eff.resize(m_z_binEdges.size()-1);
  m_proxy_sigma.resize(m_z_binEdges.size()-1);
  m_stacked_rad_eff_arr.resize(m_z_binEdges.size()-1);
  m_stacked_rad_sigma_arr.resize(m_z_binEdges.size()-1);
  m_single_deltasigma_t.resize(m_z_binEdges.size()-1);
  m_single_deltasigma_x.resize(m_z_binEdges.size()-1);
  m_single_deltasigma_err.resize(m_z_binEdges.size()-1);
  m_single_ngal_arr.resize(m_z_binEdges.size()-1);
  m_nClu_inBin.resize(m_z_binEdges.size()-1);
  m_z_eff.resize(m_z_binEdges.size()-1);
  m_z_sigma.resize(m_z_binEdges.size()-1);
  for (size_t i=0; i<m_z_binEdges.size()-1; i++){
    m_stacked_ngal_arr[i].resize(m_proxy_binEdges[i].size()-1);
    m_stacked_deltasigma_t[i].resize(m_proxy_binEdges[i].size()-1);
    m_stacked_deltasigma_x[i].resize(m_proxy_binEdges[i].size()-1);
    m_stacked_deltasigma_t_err[i].resize(m_proxy_binEdges[i].size()-1);
    m_stacked_deltasigma_x_err[i].resize(m_proxy_binEdges[i].size()-1);
    m_stacked_deltasigma_err[i].resize(m_proxy_binEdges[i].size()-1);
    m_stacked_wetasquareSum[i].resize(m_proxy_binEdges[i].size()-1);
    m_stacked_wSum[i].resize(m_proxy_binEdges[i].size()-1);
    m_stacked_wetaSum[i].resize(m_proxy_binEdges[i].size()-1);
    m_stacked_deltasigmaSum[i].resize(m_proxy_binEdges[i].size()-1);
    m_stacked_deltasigma_wei[i].resize(m_proxy_binEdges[i].size()-1);
    m_proxy_eff[i].resize(m_proxy_binEdges[i].size()-1);
    m_proxy_sigma[i].resize(m_proxy_binEdges[i].size()-1);
    m_stacked_rad_eff_arr[i].resize(m_proxy_binEdges[i].size()-1);
    m_stacked_rad_sigma_arr[i].resize(m_proxy_binEdges[i].size()-1);
    m_single_deltasigma_t[i].resize(m_proxy_binEdges[i].size()-1);
    m_single_deltasigma_x[i].resize(m_proxy_binEdges[i].size()-1);
    m_single_deltasigma_err[i].resize(m_proxy_binEdges[i].size()-1);
    m_single_ngal_arr[i].resize(m_proxy_binEdges[i].size()-1);
    m_nClu_inBin[i].resize(m_proxy_binEdges[i].size()-1,0);
    m_z_eff[i].resize(m_proxy_binEdges[i].size()-1);
    m_z_sigma[i].resize(m_proxy_binEdges[i].size()-1);
    for (size_t j=0; j<m_proxy_binEdges[i].size()-1; j++){
      m_stacked_ngal_arr[i][j].resize(m_rad_arr.size()-1);
      m_stacked_deltasigma_t[i][j].resize(m_rad_arr.size()-1);
      m_stacked_deltasigma_x[i][j].resize(m_rad_arr.size()-1);
      m_stacked_deltasigma_t_err[i][j].resize(m_rad_arr.size()-1);
      m_stacked_deltasigma_x_err[i][j].resize(m_rad_arr.size()-1);
      m_stacked_deltasigma_err[i][j].resize(m_rad_arr.size()-1);
      m_stacked_wetasquareSum[i][j].resize(m_rad_arr.size()-1);
      m_stacked_wSum[i][j].resize(m_rad_arr.size()-1);
      m_stacked_wetaSum[i][j].resize(m_rad_arr.size()-1);
      m_stacked_deltasigmaSum[i][j].resize(m_rad_arr.size()-1);
      m_stacked_deltasigma_wei[i][j].resize(m_rad_arr.size()-1);
      m_stacked_rad_eff_arr[i][j].resize(m_rad_arr.size()-1);
      m_single_deltasigma_t[i][j].resize(m_rad_arr.size()-1);
      m_single_deltasigma_x[i][j].resize(m_rad_arr.size()-1);
      m_single_deltasigma_err[i][j].resize(m_rad_arr.size()-1);
      m_single_ngal_arr[i][j].resize(m_rad_arr.size()-1);
      m_stacked_rad_sigma_arr[i][j].resize(m_rad_arr.size()-1);
    }}

  m_deltasigma_cov_matr.resize(m_z_binEdges.size()-1);
  for (size_t i=0; i<m_z_binEdges.size()-1; i++){
    m_deltasigma_cov_matr[i].resize(m_proxy_binEdges[i].size()-1);
    for (size_t j=0; j<m_proxy_binEdges[i].size()-1; j++){
      m_deltasigma_cov_matr[i][j].resize(m_rad_arr.size()-1);
      for (size_t k=0; k<m_rad_arr.size()-1; k++){
	m_deltasigma_cov_matr[i][j][k].resize(m_rad_arr.size()-1);
      }}}
}

// ============================================================================

void cbl::measure::stackprofile::StackedDensityProfile::m_linked_list ()
{
  std::vector<double> ra_limits = {1000,-1000}, dec_limits = {1000,-1000};
  for (size_t i=0; i<m_galData->nObjects(); i++){
    if (m_galData->ra(i) < ra_limits[0]) ra_limits[0] = m_galData->ra(i);
    if (m_galData->ra(i) > ra_limits[1]) ra_limits[1] = m_galData->ra(i);
    if (m_galData->dec(i) < dec_limits[0]) dec_limits[0] = m_galData->dec(i);
    if (m_galData->dec(i) > dec_limits[1]) dec_limits[1] = m_galData->dec(i);
  }
  m_ra_start = ra_limits[0]; m_dec_start = dec_limits[0];
  m_nPix.resize(2);
  m_nPix[0] = (ra_limits[1]-m_ra_start)/m_pix_size +1;
  m_nPix[1] = (dec_limits[1]-m_dec_start)/m_pix_size +1;

  m_last.resize(m_nPix[0], std::vector<int>(m_nPix[1]));
  m_first.resize(m_nPix[0], std::vector<int>(m_nPix[1]));
  m_next.resize(m_galData->nObjects(),-1);

  for (size_t i=0; i<m_next.size(); i++)
    m_next[i] = -1;
  for (size_t i=0; i<m_first.size(); i++){
    for (size_t j=0; j<m_first[i].size(); j++){
      m_first[i][j] = -1;
    }
  }
  int pix[2];
  for (size_t i=0; i<m_next.size(); i++){
    // Pixel position
    pix[0] = std::floor((m_galData->ra(i)-m_ra_start)/m_pix_size);
    pix[1] = std::floor((m_galData->dec(i)-m_dec_start)/m_pix_size);
    if (m_first[pix[0]][pix[1]]==-1){ // If m_first==-1 then this is the first particle of the list
      m_first[pix[0]][pix[1]] = i;
      m_last[pix[0]][pix[1]] = i;
    }
    else{ // otherwise I search for the end of the list and put it there.                 
	m_next[m_last[pix[0]][pix[1]]] = i;
	m_last[pix[0]][pix[1]] = i;
    }
  }
}

// ============================================================================

void cbl::measure::stackprofile::StackedDensityProfile::m_add_galaxy (const int i_gal, const int clu_index, const double coscludec, const double clu_dist)
{ 
  std::vector<double> ang_dist(2);
  ang_dist[0] = (m_cluData->ra(clu_index)-m_galData->ra(i_gal))*coscludec;
  ang_dist[1] = m_galData->dec(i_gal)-m_cluData->dec(clu_index);

  double alpha1 = m_cluData->ra(clu_index)/180.*cbl::par::pi;
  double delta1 = m_cluData->dec(clu_index)/180.*cbl::par::pi;
  double alpha2 = m_galData->ra(i_gal)/180.*cbl::par::pi;
  double delta2 = m_galData->dec(i_gal)/180.*cbl::par::pi;

  // Consider a spherical triangle with sides a, b, c. We want to find c, the angular distance between two sources.
  // The law of haversines is: cos(c) = cos(a)cos(b) + sin(a)sin(b)cos(C), where C is the angle of the corner opposite c.
  // Suppose also that a and b intersect in the north pole of a unit sphere, so that a = 90°-dec1, b = 90°-dec2, and C = ra1-ra2. Thus the formula becomes:
  // cos(c) = sin(dec1)sin(dec2) + cos(dec1)cos(dec2)cos(ra1-ra2)
  // NOTE: this formula is valid only if c is small!
  double dist = acos(sin(delta1)*sin(delta2)+cos(delta1)*cos(delta2)*cos(alpha1-alpha2))*clu_dist/1.e+6;
  double dist2 = dist*dist;

  if (dist2 < m_rad_arr[m_rad_arr.size()-1]*m_rad_arr[m_rad_arr.size()-1] && dist2 >= m_rad_arr[0]*m_rad_arr[0] && m_galData->redshift(i_gal) > m_cluData->redshift(clu_index)+0.05){
    std::vector<int> idxs = {clu_index, i_gal};
    std::vector<bool> single_selections = { (this->*m_colourSel)(i_gal), (this->*m_photzSel)(idxs) };
    bool selection = m_logicSelection(single_selections);      
    if (selection){
      // locate correct annulus
      std::vector<double>::iterator it = upper_bound(m_rad_arr.begin(), m_rad_arr.end(), dist);
      int p = it-m_rad_arr.begin()-1;
      m_ngal_arr[p]++;

      double ds = 2.99792e+09*m_cosm->D_A(0.,m_galData->redshift(i_gal))/cbl::par::cc*100; // in pc/h
      double dds = 2.99792e+09*m_cosm->D_A(m_cluData->redshift(clu_index),m_galData->redshift(i_gal))/cbl::par::cc*100; // in pc/h 
      double eta = clu_dist*dds/ds; // (sigma_crit / (c*c/4/pi/G))^-1
      std::vector<double> ell = {m_galData->shear1(i_gal),m_galData->shear2(i_gal)};

      double phi = {atan2(ang_dist[1],ang_dist[0])}; // azimuthal angle between cluster and lensed source
      double cos2phi = cos(2*phi); double sin2phi = sin(2*phi);
      std::vector<double> new_ell = {-ell[0]*cos2phi-ell[1]*sin2phi, ell[0]*sin2phi-ell[1]*cos2phi}; // tang. shear = - (shear1*cos(2phi) + shear2*sin(2phi))
	  
      m_wSum[p] += m_galData->lensingWeight(i_gal);
      m_wetaSum[p] += m_galData->lensingWeight(i_gal)*eta;
      m_wetasquareSum[p] += m_galData->lensingWeight(i_gal)*eta*eta;	  
      m_deltasigmaSum[p] += m_galData->lensingWeight(i_gal)*eta*eta*m_galData->lensingCalib(i_gal); // uncertainty on m_wetasquareSum[p]	  
      m_deltasigma_t[p] += new_ell[0]*m_galData->lensingWeight(i_gal)*eta;
      m_deltasigma_x[p] += new_ell[1]*m_galData->lensingWeight(i_gal)*eta;
      m_rad_eff_arr[p] += m_galData->lensingWeight(i_gal)*eta*eta*std::pow(dist,-m_rad_alpha);
      m_rad_sigma_arr[p] += m_galData->lensingWeight(i_gal)*eta*eta*std::pow(dist,2.);
    }
  }
}


// ============================================================================

void cbl::measure::stackprofile::StackedDensityProfile::m_profile (const int clu_index)
{
  for (size_t i=0; i<m_rad_arr.size()-1; i++){
    m_ngal_arr[i]=0; m_deltasigma_t[i]=0; m_deltasigma_x[i]=0; m_deltasigma_err[i]=0; m_wetasquareSum[i]=0; m_rad_eff_arr[i]=0; m_rad_sigma_arr[i]=0; m_wSum[i]=0; m_wetaSum[i]=0; m_deltasigmaSum[i]=0;
  }
  
  const double clu_dist = 2.99792e+09*m_cosm->D_A(0.,m_cluData->redshift(clu_index))/cbl::par::cc*100; // in pc/h
  const double coscludec = cos(m_cluData->dec(clu_index)*cbl::par::pi/180.); // by multiplying by pi/180, we convert degrees to radians
  std::vector<double> clu_pix = {std::floor((m_cluData->ra(clu_index)-m_ra_start)/m_pix_size), std::floor((m_cluData->dec(clu_index)-m_dec_start)/m_pix_size)};
  double max_rad = m_rad_arr[m_rad_arr.size()-1]/clu_dist*1.e+6*180/cbl::par::pi;

  // limits of relevant pixels
  std::vector<int> pix_ra_lim(2); std::vector<int> pix_dec_lim(2);
  pix_ra_lim[0] = std::min(std::max(0,int(clu_pix[0] - max_rad/coscludec/m_pix_size)),int(m_nPix[0]-1));
  pix_ra_lim[1] = std::max(std::min(int(m_nPix[0]-1),int(clu_pix[0] + max_rad/coscludec/m_pix_size + 1)),0);
  pix_dec_lim[0] = std::min(std::max(0,int(clu_pix[1] - max_rad/m_pix_size)),int(m_nPix[1]-1));
  pix_dec_lim[1] = std::max(std::min(int(m_nPix[1]-1),int(clu_pix[1] + max_rad/m_pix_size + 1)),0);

  std::vector<int> pix(2);
  for (pix[0] = pix_ra_lim[0]; pix[0] <= pix_ra_lim[1]; pix[0]++){
    for (pix[1] = pix_dec_lim[0]; pix[1] <= pix_dec_lim[1]; pix[1]++){
      int igal = m_first[pix[0]][pix[1]];
      while ( igal!=-1 ){
	m_add_galaxy(igal,clu_index,coscludec,clu_dist);
	igal = m_next[igal];
      }
    }
  }
  
  // Normalize
  for (size_t i=0; i<m_rad_arr.size()-1; i++){
    double deltasigmaMean = m_deltasigmaSum[i]/m_wetasquareSum[i];
    m_deltasigma_t[i] /= (1.+deltasigmaMean)*(m_wetasquareSum[i]/m_sigma_fac); // final units: h*M_sun/(pc^2)
    m_deltasigma_x[i] /= (1.+deltasigmaMean)*(m_wetasquareSum[i]/m_sigma_fac);
    m_deltasigma_err[i] = 1./sqrt(m_wetasquareSum[i])*m_sigma_fac; // statistical uncertainty from weighted mean
    m_rad_eff_arr[i] = std::pow(m_rad_eff_arr[i]/m_wetasquareSum[i],-1./m_rad_alpha);
    m_rad_sigma_arr[i] = std::pow(m_rad_sigma_arr[i]/m_wetasquareSum[i] - std::pow(m_rad_eff_arr[i],2.), 0.5);
  }
}

// ============================================================================

void cbl::measure::stackprofile::StackedDensityProfile::m_stacker ()
{

  m_set_multiplicative_calib();
  
  std::cout<<std::endl; coutCBL<<"Performing the stacking..."<<std::endl;
  cbl::WarningMsgCBL("If you do not set \'radians\' as the coordinate units in the galaxy and cluster catalogues, the code is much slower! If not already done, set \'radians\'.","m_stacker","StackedDensityProfile.cpp");
  
  // This is the inverse of Fac_eta_ToInv_Sigma_Cr in Mauro Sereno's code
  m_sigma_fac = 1.6628e+12; // c^2/(4piG) in Msun/pc
  
  // linked_list pix_size could depend on cluster min z and rad_max
  m_linked_list();

  // --------------------------------------------------------------------------------------------
  // ---------------------------------- Perform the stacking ------------------------------------
  // --------------------------------------------------------------------------------------------
  std::vector<std::vector<double>> wetasquareSum_j(m_z_binEdges.size()-1,std::vector<double>(m_proxy_binEdges[0].size()-1)); // Array used for the stacked proxy and redshift
  std::vector<std::vector<double>> wetasquareSum_tot(m_z_binEdges.size()-1,std::vector<double>(m_proxy_binEdges[0].size()-1)); // Array used for the stacked proxy and redshift

  double nClu = (double)(m_cluData->nObjects());
  for (int i=0; i<nClu; i++){
    if (m_cluData->redshift(i) >= m_z_binEdges[0] && m_cluData->redshift(i) < m_z_binEdges[m_z_binEdges.size()-1] && m_cluData->sn(i) > m_SN_min){
      std::vector<double>::iterator it_p = std::upper_bound(m_z_binEdges.begin(), m_z_binEdges.end(), m_cluData->redshift(i)); // find redshift and mass proxy bin where to put cluster
      int p = it_p-m_z_binEdges.begin()-1;
	  
      if (m_cluData->mass_proxy(i) >= m_proxy_binEdges[p][0] && m_cluData->mass_proxy(i) < m_proxy_binEdges[p][m_proxy_binEdges[p].size()-1]){
	std::vector<double>::iterator it_q = std::upper_bound(m_proxy_binEdges[p].begin(), m_proxy_binEdges[p].end(), m_cluData->mass_proxy(i));
	int q = it_q-m_proxy_binEdges[p].begin()-1;
	  
	if (p>=0 && p < (int)(m_z_binEdges.size()) && q>=0 && q < (int)(m_proxy_binEdges[p].size())){
	  m_nClu_inBin[p][q] ++;

	  m_profile(i);
	  for (size_t j=0; j<m_rad_arr.size()-1; j++){
	    if (m_ngal_arr[j] > 0){
	      double err = m_deltasigma_err[j]; double wei = 1./err/err;
	      m_stacked_deltasigma_t[p][q][j] += m_deltasigma_t[j]*wei; m_stacked_deltasigma_x[p][q][j] += m_deltasigma_x[j]*wei;
	      m_stacked_deltasigma_wei[p][q][j] += wei; m_stacked_ngal_arr[p][q][j] += m_ngal_arr[j];
	      m_single_deltasigma_t[p][q][j].emplace_back(m_deltasigma_t[j]); m_single_deltasigma_x[p][q][j].emplace_back(m_deltasigma_x[j]); // For bootstrap/resampling in general
	      m_single_deltasigma_err[p][q][j].emplace_back(m_deltasigma_err[j]); m_single_ngal_arr[p][q][j].emplace_back(m_ngal_arr[j]); // For bootstrap/resampling in general
	      wetasquareSum_j[p][q] += m_wetasquareSum[j];
	      m_stacked_wetasquareSum[p][q][j] += m_wetasquareSum[j];
	      m_stacked_rad_eff_arr[p][q][j] += pow(m_rad_eff_arr[j],m_obs_gamma)*m_wetasquareSum[j];
	      if (m_ngal_arr[j] > 1) m_stacked_rad_sigma_arr[p][q][j] += pow(m_rad_sigma_arr[j],m_obs_gamma*2.)*m_wetasquareSum[j];
	    }
	  }
	  m_proxy_eff[p][q] += pow(m_cluData->mass_proxy(i),m_obs_gamma)*wetasquareSum_j[p][q]; // wetasquaresum = \Sum (weight * (sigma_crit / (c*c/4/pi/G))^-2)
	  m_proxy_sigma[p][q] += pow(m_cluData->mass_proxy(i),2*m_obs_gamma)*wetasquareSum_j[p][q]; 
	  m_z_eff[p][q] += pow(m_cluData->redshift(i),m_obs_gamma)*wetasquareSum_j[p][q];
	  m_z_sigma[p][q] += pow(m_cluData->redshift(i),2*m_obs_gamma)*wetasquareSum_j[p][q];
	  wetasquareSum_tot[p][q] += wetasquareSum_j[p][q];
	  wetasquareSum_j[p][q] = 0;
	}
      }
    }
    const int progress = (int)((i+1)/nClu*100);
    coutCBL << progress << "% \r"; std::cout.flush();
  }
  for (size_t p=0; p<m_z_binEdges.size()-1; p++){
    for (size_t q=0; q<m_proxy_binEdges[p].size()-1; q++){
      m_proxy_eff[p][q] = pow(m_proxy_eff[p][q]/wetasquareSum_tot[p][q],1./m_obs_gamma);
      m_proxy_sigma[p][q] = pow((m_proxy_sigma[p][q]/wetasquareSum_tot[p][q] - pow(m_proxy_eff[p][q],2*m_obs_gamma))/m_nClu_inBin[p][q],1./2./m_obs_gamma);
      m_z_eff[p][q] = pow(m_z_eff[p][q]/wetasquareSum_tot[p][q],1./m_obs_gamma);            
      m_z_sigma[p][q] = pow((m_z_sigma[p][q]/wetasquareSum_tot[p][q] - pow(m_z_eff[p][q],2*m_obs_gamma))/m_nClu_inBin[p][q],1./2./m_obs_gamma);
      for (size_t j=0; j<m_rad_arr.size()-1; j++){
	m_stacked_deltasigma_t[p][q][j] /= m_stacked_deltasigma_wei[p][q][j];
	m_stacked_deltasigma_x[p][q][j] /= m_stacked_deltasigma_wei[p][q][j];
	m_stacked_deltasigma_t_err[p][q][j] = sqrt(1./m_stacked_deltasigma_wei[p][q][j]);
	m_stacked_deltasigma_x_err[p][q][j] = sqrt(1./m_stacked_deltasigma_wei[p][q][j]);
	// Make m_rad_arr
	// for radius we consider total scatter of galaxy data in bin (from weighted average of scatters), and not sample-average error
	m_stacked_rad_eff_arr[p][q][j] = pow(m_stacked_rad_eff_arr[p][q][j]/m_stacked_wetasquareSum[p][q][j],1./m_obs_gamma);;
	m_stacked_rad_sigma_arr[p][q][j] = pow(m_stacked_rad_sigma_arr[p][q][j]/m_stacked_wetasquareSum[p][q][j],1./2./m_obs_gamma);
      }}}
}


// ============================================================================

std::shared_ptr<data::Data> cbl::measure::stackprofile::StackedDensityProfile::m_make_bootstrap (const std::vector<int> z_proxy_bin)
{
  m_stacker(); // Perform the stacking
  
  std::cout<<std::endl; coutCBL<<"Computing the bootstrap covariance matrix..."<<std::endl<<std::endl;  
  std::default_random_engine generator;
  for (size_t p=0; p<m_z_binEdges.size()-1; p++){
    for (size_t q=0; q<m_proxy_binEdges[p].size()-1; q++){
      std::uniform_int_distribution<int> distribution(0,m_single_ngal_arr[p][q][0].size()-1);
      std::vector< std::vector<double> > deltasigma_t_boot(m_n_resampling,std::vector<double>(m_rad_arr.size()-1,0.));
      std::vector< std::vector<double> > deltasigma_x_boot(m_n_resampling,std::vector<double>(m_rad_arr.size()-1,0.));
      std::vector< std::vector<double> > deltasigma_wei_boot(m_n_resampling,std::vector<double>(m_rad_arr.size()-1,0.));
      std::vector<double> deltasigma_t_boot_mean(m_rad_arr.size()-1,0.);

      for (int i=0; i<m_n_resampling; i++){
	for (size_t j=0; j<m_single_ngal_arr[p][q][0].size(); j++){
	  int clu_num = distribution(generator); // random extraction of a cluster
	  for (size_t k=0; k<m_rad_arr.size()-1; k++){
	    if (m_single_ngal_arr[p][q][k][clu_num] > 0){
	      deltasigma_t_boot[i][k] += m_single_deltasigma_t[p][q][k][clu_num]/m_single_deltasigma_err[p][q][k][clu_num]/m_single_deltasigma_err[p][q][k][clu_num];
	      deltasigma_x_boot[i][k] += m_single_deltasigma_x[p][q][k][clu_num]/m_single_deltasigma_err[p][q][k][clu_num]/m_single_deltasigma_err[p][q][k][clu_num];
	      deltasigma_wei_boot[i][k] += 1./m_single_deltasigma_err[p][q][k][clu_num]/m_single_deltasigma_err[p][q][k][clu_num];
	    }}}
	for (size_t k=0; k<m_rad_arr.size()-1; k++){
	  deltasigma_t_boot[i][k] /= deltasigma_wei_boot[i][k];
	  deltasigma_x_boot[i][k] /= deltasigma_wei_boot[i][k];
	} 
      }
      for (size_t k=0; k<m_rad_arr.size()-1; k++){
	for (int i=0; i<m_n_resampling; i++){
	  deltasigma_t_boot_mean[k] += deltasigma_t_boot[i][k];
	}
	deltasigma_t_boot_mean[k] /= double(m_n_resampling);
      }
      for (size_t k=0; k<m_rad_arr.size()-1; k++){
	for (size_t l=0; l<m_rad_arr.size()-1; l++){
	  for (int i=0; i<m_n_resampling; i++){
	    m_deltasigma_cov_matr[p][q][k][l] += (deltasigma_t_boot[i][k]-deltasigma_t_boot_mean[k])*(deltasigma_t_boot[i][l]-deltasigma_t_boot_mean[l]);
	  }
	  m_deltasigma_cov_matr[p][q][k][l] /= double(m_n_resampling);
	}}
    }}

  std::vector<double> bins = m_stacked_rad_eff_arr[z_proxy_bin[0]][z_proxy_bin[1]];
  std::vector<double> profile = m_stacked_deltasigma_t[z_proxy_bin[0]][z_proxy_bin[1]];
  std::vector<double> error; // Set the diagonal of the covariance matrix as the error
  std::vector<std::vector<double>> extra_info (4);
  for (size_t i=0; i<m_rad_arr.size()-1; i++){
    error.emplace_back(sqrt(m_deltasigma_cov_matr[z_proxy_bin[0]][z_proxy_bin[1]][i][i]));
    extra_info[0].emplace_back(m_z_eff[z_proxy_bin[0]][z_proxy_bin[1]]);
    extra_info[1].emplace_back(m_z_sigma[z_proxy_bin[0]][z_proxy_bin[1]]);
    extra_info[2].emplace_back(m_proxy_eff[z_proxy_bin[0]][z_proxy_bin[1]]);
    extra_info[3].emplace_back(m_proxy_sigma[z_proxy_bin[0]][z_proxy_bin[1]]);
  }
  auto dataset = std::make_shared<data::Data1D_extra> (data::Data1D_extra(bins, profile, error, extra_info));
  
  return dataset;
}

// ============================================================================

bool cbl::measure::stackprofile::StackedDensityProfile::m_check_file (const std::string checked_file, const std::vector<int> z_proxy_bin)
{
  bool out = false;
  
  m_inputs_to_str.clear();
  
  m_inputs_to_str.emplace_back(cbl::conv(m_cosm->OmegaM(),cbl::par::fDP5));
  m_inputs_to_str.emplace_back(cbl::conv(m_cosm->Omega_baryon(),cbl::par::fDP5));
  m_inputs_to_str.emplace_back(cbl::conv(m_cosm->Omega_k(),cbl::par::fDP5));
  m_inputs_to_str.emplace_back(cbl::conv(m_cosm->Omega_neutrinos(),cbl::par::fDP5));
  m_inputs_to_str.emplace_back(cbl::conv(pow(10,9)*m_cosm->scalar_amp(),cbl::par::fDP3));
  m_inputs_to_str.emplace_back(cbl::conv(m_cosm->tau(),cbl::par::fDP5));
  m_inputs_to_str.emplace_back(cbl::conv(m_cosm->n_spec(),cbl::par::fDP5));
  m_inputs_to_str.emplace_back(cbl::conv(m_cosm->hh(),cbl::par::fDP5));

  m_inputs_to_str.emplace_back(cbl::conv(m_galData->nObjects(),cbl::par::fINT));
  m_inputs_to_str.emplace_back(cbl::conv(m_cluData->nObjects(),cbl::par::fINT));

  for (size_t i=0; i<m_m_calib.size(); i++) {
    m_inputs_to_str.emplace_back(cbl::conv(m_m_calib[i][0],cbl::par::fDP5));
    m_inputs_to_str.emplace_back(cbl::conv(m_m_calib[i][1],cbl::par::fDP5));
    m_inputs_to_str.emplace_back(cbl::conv(m_m_calib_zEdges[i][0],cbl::par::fDP5));
    m_inputs_to_str.emplace_back(cbl::conv(m_m_calib_zEdges[i][1],cbl::par::fDP5));
  }
  
  m_inputs_to_str.emplace_back(m_colour_sel);
  m_inputs_to_str.emplace_back(m_zphot_sel);
  for (size_t j=0; j<m_zphot_sel_pars.size(); j++){
    m_inputs_to_str.emplace_back(cbl::conv(m_zphot_sel_pars[j],cbl::par::fDP5));
  }
  m_inputs_to_str.emplace_back(m_logic_sel);
  for (size_t j=0; j<m_z_binEdges.size(); j++){
    m_inputs_to_str.emplace_back(cbl::conv(m_z_binEdges[j],cbl::par::fDP5));
  }
  for (size_t i=0; i<m_proxy_binEdges.size(); i++){
    for (size_t j=0; j<m_proxy_binEdges[i].size(); j++){
      m_inputs_to_str.emplace_back(cbl::conv(m_proxy_binEdges[i][j],cbl::par::fDP5));
    }
  }
  for (size_t j=0; j<m_rad_arr.size(); j++){
    m_inputs_to_str.emplace_back(cbl::conv(m_rad_arr[j],cbl::par::fDP5));
  }
  m_inputs_to_str.emplace_back(cbl::conv(m_SN_min,cbl::par::fDP5));
  m_inputs_to_str.emplace_back(cbl::conv(m_pix_size,cbl::par::fDP5));
  m_inputs_to_str.emplace_back(cbl::conv(m_n_resampling,cbl::par::fINT));
  m_inputs_to_str.emplace_back(cbl::conv(m_rad_alpha,cbl::par::fDP5));
  m_inputs_to_str.emplace_back(cbl::conv(m_obs_gamma,cbl::par::fDP5));

  std::string check_string = "# ";
  for (size_t i=0; i<m_inputs_to_str.size(); i++){
    check_string += m_inputs_to_str[i];
    check_string += "  ";
  }
  std::ifstream file_input(checked_file.c_str(),std::ios::in);
  std::vector<double> rad, DSigma, error;
  std::vector<std::vector<double>> extra_info (4);
  if (file_input.std::ifstream::is_open()){
    std::string line;
    std::getline(file_input,line);
    if (line == check_string){
      out = true;
      while (std::getline(file_input,line)){
	if (line.at(0) != '#'){
	  std::stringstream ss(line);            
	  std::vector<double> num; double NUM;
	  while (ss>>NUM) num.emplace_back(NUM);
	  if ((int)(num[0])==z_proxy_bin[0] && (int)(num[1])==z_proxy_bin[1]){
	    rad.emplace_back(num[6]);
	    DSigma.emplace_back(num[8]);  
	    error.emplace_back(num[14]);
	    extra_info[0].emplace_back(num[2]);
	    extra_info[1].emplace_back(num[3]);
	    extra_info[2].emplace_back(num[4]);
	    extra_info[3].emplace_back(num[5]);
	  }
	}
      }
      if (rad.size()==0)
	cbl::ErrorCBL("The chosen redshift and proxy bins are not valid!","StackedDensityProfile","StackedDensityProfile.cpp");
    }
  }
  file_input.clear(); file_input.close();

  m_dataset = std::make_shared<data::Data1D_extra> (data::Data1D_extra(rad, DSigma, error, extra_info));
  
  return out;
}

// ============================================================================

void cbl::measure::stackprofile::StackedDensityProfile::m_write(const std::string output_dir, const std::string output_file)
{
  const string mkdir = "mkdir -p "+output_dir; if (system(mkdir.c_str())) {}

  std::ofstream out_stack (output_dir+"/"+output_file+".dat",std::ios::out);
  
  // Write the inputs in the header
  out_stack << "# ";
  for (size_t i=0; i<m_inputs_to_str.size(); i++)
    out_stack << m_inputs_to_str[i] << "  ";
  out_stack << std::endl;
  
  // Write the rest of the file
  out_stack << "# Redshift array:" << std::endl << "# ";
  for (size_t i=0; i<m_z_binEdges.size(); i++)
    out_stack << m_z_binEdges[i] << "  ";
  out_stack << std::endl;
  out_stack << "# Proxy arrays:";
  for (size_t i=0; i<m_proxy_binEdges.size(); i++){
    out_stack << std::endl << "# ";
    for (size_t j=0; j<m_proxy_binEdges[i].size(); j++){
      out_stack << m_proxy_binEdges[i][j] << "  ";
      if (i==m_proxy_binEdges.size()-1 && j==m_proxy_binEdges[i].size()-1)
	out_stack << std::endl;
    }
  }  
  out_stack << "# 1) redshift bin" << std::endl;
  out_stack << "# 2) proxy bin" << std::endl;
  out_stack << "# 3) mean redshift" << std::endl;
  out_stack << "# 4) redshift sigma" << std::endl;
  out_stack << "# 5) mean proxy" << std::endl;
  out_stack << "# 6) proxy sigma" << std::endl;
  out_stack << "# 7) mean radius" << std::endl;
  out_stack << "# 8) radius sigma" << std::endl;
  out_stack << "# 9) deltasigma_+ [h*MSun/pc^2]" << std::endl;
  out_stack << "# 10) deltasigma_+ error [h*MSun/pc^2]" << std::endl;
  out_stack << "# 11) deltasigma_x [h*MSun/pc^2]" << std::endl;
  out_stack << "# 12) deltasigma_x error [h*MSun/pc^2]" << std::endl;
  out_stack << "# 13) number of clusters in the bins of z and proxy" << std::endl;
  out_stack << "# 14) total weight" << std::endl;
  out_stack << "# 15) bootstrap error on deltasigma_+" << std::endl;
  out_stack << "# 16) rad min" << std::endl;
  out_stack << "# 17) rad max" << std::endl;
  out_stack << "# 18) number of background galaxies in the bins of z, proxy and radius" << std::endl;

  for (size_t i=0; i<m_z_binEdges.size()-1; i++){
    for (size_t j=0; j<m_proxy_binEdges[i].size()-1; j++){
      for (size_t k=0; k<m_rad_arr.size()-1; k++){
	out_stack << i << "  " << j << "  " << m_z_eff[i][j] << "  " << m_z_sigma[i][j] << "  " << m_proxy_eff[i][j] << "  " << m_proxy_sigma[i][j] << "  " << m_stacked_rad_eff_arr[i][j][k] << "  " << m_stacked_rad_sigma_arr[i][j][k] << "  " << m_stacked_deltasigma_t[i][j][k] << "  " << m_stacked_deltasigma_t_err[i][j][k] << "  " << m_stacked_deltasigma_x[i][j][k] << "  " << m_stacked_deltasigma_x_err[i][j][k] << "  " << m_nClu_inBin[i][j] << "  " << m_stacked_wetasquareSum[i][j][k] << "  " << sqrt(m_deltasigma_cov_matr[i][j][k][k]) << "  " << m_rad_arr[k] << "  " << m_rad_arr[k+1] << "  " << m_stacked_ngal_arr[i][j][k] << std::endl;
      }
    }
  }
  out_stack.clear(); out_stack.close(); coutCBL << "I wrote the file: " << output_dir+"/"+output_file+".dat" << std::endl << std::endl;
}

// ============================================================================

void cbl::measure::stackprofile::StackedDensityProfile::measure (const std::vector<int> z_proxy_bin, const std::string output_dir, const std::string output_file, const ErrorType errorType)
{
  bool file_exists = m_check_file(output_dir+"/"+output_file+".dat", z_proxy_bin);

  if (file_exists){
    coutCBL << "I've just read the already existing stacking file!" << std::endl; 
  }else{
    switch (errorType){
    case (ErrorType::_Bootstrap_):
      m_dataset = m_make_bootstrap(z_proxy_bin);
      break;  
    default:
      ErrorCBL("The input ErrorType is not allowed!", "measure", "StackedDensityProfile.cpp");
    }
    m_write(output_dir, output_file);
  }
}

// ============================================================================

void cbl::measure::stackprofile::StackedDensityProfile::write (const std::string dir, const std::string file)
{
  std::string mkdir = "mkdir -p "+dir;
  if (system(mkdir.c_str())) {}
  
  std::string header = "bin_center  profile  error  mean_redshift  mean_redshift_error  mean_mass_proxy  mean_mass_proxy_error";
  m_dataset->write(dir, file, header, 4, 8);
}

// ============================================================================

cbl::measure::stackprofile::StackedDensityProfile::StackedDensityProfile (cosmology::Cosmology cosm, const catalogue::Catalogue gal_cat, const catalogue::Catalogue clu_cat, const std::string colour_sel, const std::string zphot_sel, std::vector<double> zphot_sel_pars, const std::string logic_sel, std::vector<double> z_binEdges, std::vector<std::vector<double>> proxy_binEdges, const double rad_min, const double rad_max, const int nRad, const bool log_rad, const double SN_min, const double pix_size, const std::vector<std::vector<double>> multiplicative_calibration_stats, const std::vector<std::vector<double>> multiplicative_calibration_zEdges, const int n_resampling, const double rad_alpha, const double obs_gamma)
{
  // WARNING: if you change the arguments of this constructor, you must edit m_check_file !!!

  m_cosm = std::make_shared<cosmology::Cosmology>(cosmology::Cosmology(std::move(cosm)));
  m_cosm->set_unit(true); // Force cosmological units
  m_galData = std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(gal_cat)));
  m_cluData = std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(clu_cat)));
  m_rad_alpha = rad_alpha;
  m_obs_gamma = obs_gamma;
  m_SN_min = SN_min;
  m_pix_size = pix_size;
  m_m_calib = multiplicative_calibration_stats;
  m_m_calib_zEdges = multiplicative_calibration_zEdges;
  m_n_resampling = n_resampling;
  m_z_binEdges = z_binEdges;
  m_proxy_binEdges = proxy_binEdges;
  m_colour_sel = colour_sel;
  m_zphot_sel = zphot_sel;
  m_zphot_sel_pars = zphot_sel_pars;
  m_logic_sel = logic_sel;

  m_check_catalogue_variables();
  m_resize(rad_min,rad_max,nRad,log_rad);
  m_set_selections();
}
