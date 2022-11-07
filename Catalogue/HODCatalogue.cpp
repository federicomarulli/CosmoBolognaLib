/********************************************************************
 *  Copyright (C) 2020 by Federico Zangrandi, Luca Stabellini and   *
 *  Sofia Contarini, Federico Marulli                               *
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
 *  @file Catalogue/Catalogue.cpp
 *
 *  @brief Methods of the class Catalogue 
 *
 *  This file contains the implementation of a construtor of class
 *  Catalogue, which populate dark matter halos following different
 *  HOD algorithms
 *
 *  @author Federico Zangrandi, Luca Stabellini, Sofia Contarini,
 *  Federico Marulli
 *
 *  @author federico.zangrandi2@studio.unibo.it,
 *  luca.stabellini@studio.unibo.it, sofia.contarini3@unibo.it
 */

#include "Catalogue.h"
#include <chrono>

using namespace std;
using namespace cbl;
using namespace glob;


// ============================================================================


double cbl::catalogue::Average_c_Zehavi_2011 (const double x, const double logMmin, const double logsigma_c)
{
  return 0.5*(1.+erf((log10(x)-logMmin)/logsigma_c));
}


// ============================================================================


double cbl::catalogue::Average_s_Zehavi_2011 (const double x, const double M0, const double M1, const double alpha)
{
  return 0.5*pow((x-M0)/M1, alpha);
}


// ============================================================================


double cbl::catalogue::Average_c_Zehavi_2005 (const double M_h, const double M_min)
{
  return (M_h>M_min) ? 1. : 0.;
}


// ============================================================================


double cbl::catalogue::Average_s_Zehavi_2005 (const double x, const double M1, const double alpha)
{
  return pow(x/M1, alpha);
}


// ============================================================================


cbl::catalogue::Catalogue::Catalogue (const Catalogue halo_catalogue, const cosmology::Cosmology &cosm, const HODType HOD_Type, const double threshold, const bool substructures, std::vector<double> parameter)
{
  // MOSTER 10 model:

  // Compute sigma_c scatter for lognormal distribution Moster2010 for central galaxies, eq. 4.3.7
  auto fsigma_c = [](double x)
		  {
		    // All the parameters below depend on Halo Mass
		    const double sigma_infinity = 0.1592;	        // 0.0569;
		    const double sigma_1 = 0.0460;			// 0.1204;
		    const double chsi = 4.2503;		     	// 6.3020;
		    const double M_2 = pow(10, 11.8045);       	// pow(10.,11.9652);
		    //sigma_c SHMR:
		    return sigma_infinity + sigma_1 * (1. - (2. / cbl::par::pi) * atan(chsi * log10(x / M_2)));
		  };

  // Phi_s normalization of modified Schechter for satellite galaxies, Moster2010
  auto fPhi_s = [](double x)
		{
		  const double lamda = 0.8032;		       //0.8285;
		  const double phi_0 = pow(10., -10.8924);         //pow(10.,-11.1622);

		  return phi_0 * pow(x, lamda);
		};

  // compute alpha index for the modified Schechter
  auto falpha_s = [](double x)
		  {
		    const double alpha_infinity = -1.3676;       //-1.3740;
		    const double alpha_1 = -0.0524;		   //-0.0309;
		    const double zeta = 9.5727;		   //4.3629;
		    const double M_3 = pow(10, 12.3646);         //pow(10,12.5730);

		    return alpha_infinity + alpha_1 * (1. - (2. / cbl::par::pi) * atan(zeta * log10(x / M_3)));
		  };

  // compute minumum stellar mass, in h unit i.e. the minimum stellar
  // mass of the galaxy sample eq.5.1.6
  auto mag_mstar = [&cosm](double x)
		   {
		     const double MoverL = 1.;
		     return MoverL*pow(10., (4.77-x)/2.5);
		   };

  auto startTimer = chrono::steady_clock::now();

  // threshold of stellar masses
  double starMass_threshold = 0;

  // ZEHAVI PARAMETERS INITIALIZATION:
  // Zehavi11
  double logMmin = 0.;
  double logsigma_c = 0.;
  double M0 = 0.;
  double M1 = 0.;
  double alpha = 0.;

  // Zehavi05
  double M_min = 0.;

  // MOSTER10 PARAMETERS INITIALIZATION, CENTRAL:
  //Moster central
  double k_c = 0.;
  double M1_c = 0.;
  double beta_c = 0.;
  double gamma_c = 0.;
  double sigma_c = 0.;

  // MOSTER10 PARAMETERS INITIALIZATION, SATELLITE:
  double k_s = 0.;
  double M1_s = 0.;
  double beta_s = 0.;
  double gamma_s = 0.;
  double Phi_s = 0.;
  double alpha_s = 0.;

  // ----- ZEHAVI 05-11 -----

  // Paramters settings:
  /*	
	threshold -> input parameter in r-band Magnitude M_{r}.
	p.[0] = Log M_min: minimum halo mass for hosting central galaxy.
	p.[1] = Sigma_{log_M}: scatter applied to the smoothed step-function (used in mean occupation number of central galaxies).
	p.[2] and [3] = Log(M_{0}) and log(M'_{1}): relative to satellite galaxies.
	p.[3] = alpha: slope of the power-law used to obtain the mean occupation number of satellite galaxies. 
  */
  if (parameter.size() < 1)
    parameter = {-99., -99., -99., -99., -99., -99., -99., -99., -99., -99., -99., -99., -99., -99., -99., -99.};
	
  // Setting of first 4 parameters based on input threshold:
  // Zheavi models have minimum stellar mass set by threshold
  if (HOD_Type==HODType::_Zehavi11_) {	
    if (threshold == -22)
      starMass_threshold = mag_mstar(-22.);
    else if (threshold == -21.5)
      starMass_threshold = mag_mstar(-21.5);
    else if (threshold == -21.)
      starMass_threshold = mag_mstar(-21.);
    else if (threshold == -20.5)
      starMass_threshold = mag_mstar(-20.5);
    else if (threshold == -20.)
      starMass_threshold = mag_mstar(-20.);
    else if (threshold == -19.5)
      starMass_threshold = mag_mstar(-19.5);
    else if (threshold == -19.)
      starMass_threshold = mag_mstar(-19.);
    else if (threshold == -18.5)
      starMass_threshold = mag_mstar(-18.5);
    else if (threshold == -18.)
      starMass_threshold = mag_mstar(-18.);
    else
      ErrorCBL("Not default parameters available for this threshold", "Catalogue", "HODCatalogue.cpp");

    for (int ii = 0; ii < 5; ii++) {
      if (parameter[ii] == -99) {
	if (ii == 0) {
	  if (threshold == -22)
	    parameter[0] = 14.06;
	  else if (threshold == -21.5)
	    parameter[0] = 13.38;
	  else if (threshold == -21.)
	    parameter[0] = 12.78;
	  else if (threshold == -20.5)
	    parameter[0] = 12.14;
	  else if (threshold == -20.)
	    parameter[0] = 11.83;
	  else if (threshold == -19.5)
	    parameter[0] = 11.57;
	  else if (threshold == -19.)
	    parameter[0] = 11.45;
	  else if (threshold == -18.5)
	    parameter[0] = 11.33;
	  else if (threshold == -18.)
	    parameter[0] = 11.18;
	  else
	    ErrorCBL("Not default parameters available for this threshold", "Catalogue", "HODCatalogue.cpp");
	}
	else if (ii == 1) {
	  if (threshold == -22)
	    parameter[1] = 0.71;
	  else if (threshold == -21.5)
	    parameter[1] = 0.69;
	  else if (threshold == -21.)
	    parameter[1] = 0.68;
	  else if (threshold == -20.5)
	    parameter[1] = 0.17;
	  else if (threshold == -20.)
	    parameter[1] = 0.25;
	  else if (threshold == -19.5)
	    parameter[1] = 0.17;
	  else if (threshold == -19.)
	    parameter[1] = 0.19;
	  else if (threshold == -18.5)
	    parameter[1] = 0.26;
	  else if (threshold == -18.)
	    parameter[1] = 0.19;
	  else
	    ErrorCBL("Not default parameters available for this threshold", "Catalogue", "HODCatalogue.cpp");
	}
	else if (ii == 2) {
	  if (threshold == -22)
	    parameter[2] = pow(10., 13.72);
	  else if (threshold == -21.5)
	    parameter[2] = pow(10., 13.35);
	  else if (threshold == -21.)
	    parameter[2] = pow(10., 12.71);
	  else if (threshold == -20.5)
	    parameter[2] = pow(10, 11.62);
	  else if (threshold == -20.)
	    parameter[2] = pow(10., 12.35);
	  else if (threshold == -19.5)
	    parameter[2] = pow(10., 12.33);
	  else if (threshold == -19.)
	    parameter[2] = pow(10., 9.77);
	  else if (threshold == -18.5)
	    parameter[2] = pow(10., 8.99);
	  else if (threshold == -18.)
	    parameter[2] = pow(10., 9.81);
	  else
	    ErrorCBL("Not default parameters available for this threshold", "Catalogue", "HODCatalogue.cpp");
	}
	else if (ii == 3) {
	  if (threshold == -22)
	    parameter[3] = pow(10, 14.80);
	  else if (threshold == -21.5)
	    parameter[3] = pow(10., 14.20);
	  else if (threshold == -21.)
	    parameter[3] = pow(10., 13.76);
	  else if (threshold == -20.5)
	    parameter[3] = pow(10., 13.43);
	  else if (threshold == -20.)
	    parameter[3] = pow(10., 12.98);
	  else if (threshold == -19.5)
	    parameter[3] = pow(10., 12.75);
	  else if (threshold == -19.)
	    parameter[3] = pow(10., 12.63);
	  else if (threshold == -18.5)
	    parameter[3] = pow(10., 12.50);
	  else if (threshold == -18.)
	    parameter[3] = pow(10., 12.42);
	  else
	    ErrorCBL("Not default parameters available for this threshold", "Catalogue", "HODCatalogue.cpp");
	}
	else {
	  if (threshold == -22)
	    parameter[4] = 1.35;
	  else if (threshold == -21.5)
	    parameter[4] = 1.09;
	  else if (threshold == -21.)
	    parameter[4] = 1.15;
	  else if (threshold == -20.5)
	    parameter[4] = 1.15;
	  else if (threshold == -20.)
	    parameter[4] = 1.00;
	  else if (threshold == -19.5)
	    parameter[4] = 0.99;
	  else if (threshold == -19.)
	    parameter[4] = 1.02;
	  else if (threshold == -18.5)
	    parameter[4] = 1.02;
	  else if (threshold == -18.)
	    parameter[4] = 1.04;
	  else
	    ErrorCBL("Not default parameters available for this threshold", "Catalogue", "HODCatalogue.cpp");
	}
      }
    }// END of setting first 4 parameters based on input threshold.

    // Average naumbers Zehavi11 but stellar mass with Moster10

    // SETTING PARAMTERS [5] to [15]:
    /*
      p.[5] = (m_c / M)_0: normalization of SHMR for central galaxies (eq. 4.3.11).
      p.[6] = logM_{1c}: log10 of characteristic mass at which the SHMR changes slope for central galaxies.
      p.[7] and p[8] = beta and gamma coefficents for central galaxies (eq. 4.3.11).
      p.[9] = sigma: scatter for the log-normal distribution of central galaxies.
      p.[10] = same as [5] but for satellite galaxies
      p.[11] = same as [6] but for satellite galaxies
      p.[12] and [13]= same as [7] ans [8] but for satellite galaxies
      p.[14] = Phi_s: parametrization by Moster10, eq.4.3.9
      p.[15] = alpha_s: parametrization by Moster10, eq.4.3.10
    */

    if (parameter[5] == -99)
      parameter[5] = 0.0267; // norma SHMR_c
    if (parameter[6] == -99)
      parameter[6] = pow(10, 11.9347); //break SHMR_c
    if (parameter[7] == -99)
      parameter[7] = 1.0059; //beta_c SHMR_c
    if (parameter[8] == -99)
      parameter[8] = 0.5611; //gamma_c SHMR_c
    if (parameter[9] == -99)
      {
	// sigma_c, default value setted from halo mass
	//cout << par::col_green << "Using default sigma_c" << par::col_default << endl;
	//parameter[ii]=sigma_c(;
      }
    if (parameter[10] == -99)
      parameter[10] = 0.0186; // norma SHMR_s
    if (parameter[11] == -99)
      parameter[11] = pow(10, 12.1988); // break SHMR_s
    if (parameter[12] == -99)
      parameter[12] = 0.7817; // beta_s SHMR_s
    if (parameter[13] == -99)
      parameter[13] = 0.7334; // gamma_s SHMR_s
    if (parameter[14] == -99)
      { // PHI_s, default value setted from halo mass

	//cout << par::col_green << "Using default Phi_s" << par::col_default << endl;
	//parameter[ii]=Phi_s;
      }
    if (parameter[15] == -99)
      { 
	// alpha_s, default value setted from halo mass
	//
	//cout << par::col_green << "Using default alpha_s" << par::col_default << endl;
	// parameter[ii]=alpha_s;
      }

    //HOD parameter Zehavi11
    logMmin = parameter[0];
    logsigma_c = parameter[1];
    M0 = parameter[2];
    M1 = parameter[3];
    alpha = parameter[4];

    //Moster central stellar masses
    k_c = parameter[5];
    M1_c = parameter[6];
    beta_c = parameter[7];
    gamma_c = parameter[8];
    sigma_c = parameter[9];

    //Moster satallite
    k_s = parameter[10];
    M1_s = parameter[11];
    beta_s = parameter[12];
    gamma_s = parameter[13];
    Phi_s = parameter[14];
    alpha_s = parameter[15];

    cout << " Chosen parameters for HOD: " << endl;
    cout << " LogMmin = " << logMmin << endl;
    cout << " Logsigma_c = " << logsigma_c << endl;
    cout << " M0 = " << M0 << " M_sun/h" << endl;
    cout << " M1 = " << M1 << " M_sun/h" << endl;
    cout << " alpha = " << alpha << endl;
    cout << " Minimum stellar mass = " << starMass_threshold << " M_sun/h" << endl
	 << endl;

    cout << " Chosen parameters for stellar masses: " << endl;
    cout << " k_c = " << k_c << endl;
    cout << " M1_c = " << M1_c << " M_sun" << endl;
    cout << " beta_c = " << beta_c << endl;
    cout << " gamma_c = " << gamma_c << endl;
    cout << " sigma_c  = " << sigma_c << endl;

    cout << " k_s = " << k_s << endl;
    cout << " M1_s = " << M1_s << " M_sun" << endl;
    cout << " beta_s = " << beta_s << endl;
    cout << " gamma_s = " << gamma_s << endl;
    cout << " Phi_s = " << Phi_s << endl;
    cout << " alpha_s = " << alpha_s << endl
	 << endl;
  }

  else if (HOD_Type == HODType::_Zehavi05_) {
    if (threshold == -22)
      starMass_threshold = mag_mstar(-22.);
    else if (threshold == -21.5)
      starMass_threshold = mag_mstar(-21.5);
    else if (threshold == -21.)
      starMass_threshold = mag_mstar(-21.);
    else if (threshold == -20.5)
      starMass_threshold = mag_mstar(-20.5);
    else if (threshold == -20.)
      starMass_threshold = mag_mstar(-20.);
    else if (threshold == -19.5)
      starMass_threshold = mag_mstar(-19.5);
    else if (threshold == -19.)
      starMass_threshold = mag_mstar(-19.);
    else if (threshold == -18.5)
      starMass_threshold = mag_mstar(-18.5);
    else if (threshold == -18.)
      starMass_threshold = mag_mstar(-18.);
    else
      ErrorCBL("Not default parameters available for this threshold", "Catalogue", "HODCatalogue.cpp");

    for (int ii = 0; ii < 3; ii++) {
      if (parameter[ii] == -99) {
	if (ii == 0) {
	  if (threshold == -22)
	    parameter[0] = pow(10, 13.91);
	  else if (threshold == -21.5)
	    parameter[0] = pow(10, 13.27);
	  else if (threshold == -21.)
	    parameter[0] = pow(10, 12.72);
	  else if (threshold == -20.5)
	    parameter[0] = pow(10, 12.30);
	  else if (threshold == -20.)
	    parameter[0] = pow(10, 12.01);
	  else if (threshold == -19.5)
	    parameter[0] = pow(10, 11.76);
	  else if (threshold == -19.)
	    parameter[0] = pow(10, 11.59);
	  else if (threshold == -18.5)
	    parameter[0] = pow(10, 11.44);
	  else if (threshold == -18.)
	    parameter[0] = pow(10, 11.27);
	  else
	    ErrorCBL("Not default parameters available for this threshold", "Catalogue", "HODCatalogue.cpp");
	}
	else if (ii == 1) {
	  if (threshold == -22)
	    parameter[1] = pow(10, 14.92);
	  else if (threshold == -21.5)
	    parameter[1] = pow(10., 14.60);
	  else if (threshold == -21.)
	    parameter[1] = pow(10., 14.09);
	  else if (threshold == -20.5)
	    parameter[1] = pow(10., 13.67);
	  else if (threshold == -20.)
	    parameter[1] = pow(10., 13.42);
	  else if (threshold == -19.5)
	    parameter[1] = pow(10., 13.15);
	  else if (threshold == -19.)
	    parameter[1] = pow(10., 12.94);
	  else if (threshold == -18.5)
	    parameter[1] = pow(10., 12.77);
	  else if (threshold == -18.)
	    parameter[1] = pow(10., 12.57);
	  else
	    ErrorCBL("Not default parameter available for this threshold", "Catalogue", "HODCatalogue.cpp");
	}
	else {
	  if (threshold == -22)
	    parameter[2] = 1.43;
	  else if (threshold == -21.5)
	    parameter[2] = 1.94;
	  else if (threshold == -21.)
	    parameter[2] = 1.39;
	  else if (threshold == -20.5)
	    parameter[2] = 1.21;
	  else if (threshold == -20.)
	    parameter[2] = 1.16;
	  else if (threshold == -19.5)
	    parameter[2] = 1.13;
	  else if (threshold == -19.)
	    parameter[2] = 1.08;
	  else if (threshold == -18.5)
	    parameter[2] = 1.01;
	  else if (threshold == -18.)
	    parameter[2] = 0.92;
	  else
	    ErrorCBL("not default parameter available for this threshold", "Catalogue", "HODCatalogue.cpp");
	}
      }
    }
    //for (int ii=0; ii<10; ii++) {
    if (parameter[3] == -99)
      parameter[3] = 0.0267; //norm SHMR_c
    if (parameter[4] == -99)
      parameter[4] = pow(10, 11.9347); //break SHMR_c
    if (parameter[5] == -99)
      parameter[5] = 1.0059; //beta_c SHMR_c
    if (parameter[6] == -99)
      parameter[6] = 0.5611; //gamma_c SHMR_c
    if (parameter[7] == -99)
      {
	//cout << par::col_green << "If sigma_c=-99 it will be setted using the halo mass..." << par::col_default << endl;
	//parameter[ii]=sigma_c(;
      }
    if (parameter[8] == -99)
      parameter[8] = 0.0186; // norma SHMR_s
    if (parameter[9] == -99)
      parameter[9] = pow(10, 12.1988); // break SHMR_s
    if (parameter[10] == -99)
      parameter[10] = 0.7817; // beta_s SHMR_s
    if (parameter[11] == -99)
      parameter[11] = 0.7334; // gamma_s SHMR_s
    if (parameter[12] == -99)
      {
	// PHI_s, default value setted from halo mass
	//cout << par::col_green << "Using default PHI_s" << par::col_default << endl;
	//parameter[ii]=Phi_s;
      }
    if (parameter[13] == -99)
      { // alpha_s
	//cout << par::col_green << "If alpha_s=-99 it will be setted using the halo mass..." << par::col_default << endl;
	// parameter[ii]=alpha_s;
      }

    M_min = parameter[0];
    M1 = parameter[1];
    alpha = parameter[2];
    //Moster stellar mass centra
    k_c = parameter[3];
    M1_c = parameter[4];
    beta_c = parameter[5];
    gamma_c = parameter[6];
    sigma_c = parameter[7];
    //Moster satallite
    k_s = parameter[8];
    M1_s = parameter[9];
    beta_s = parameter[10];
    gamma_s = parameter[11];
    Phi_s = parameter[12];
    alpha_s = parameter[13];

    coutCBL << " Chosen parameters for  HOD: " << endl;
    coutCBL << " M_min = " << M_min << " M_sun/h" << endl;
    coutCBL << " M1 = " << M1 << " M_sun/h" << endl;
    coutCBL << " alpha = " << alpha << endl;
    coutCBL << " Minimum stellar mass = " << starMass_threshold << "M_sun/h" << endl
	 << endl;

    coutCBL << " Chosen parameters for stellar masses: " << endl;
    coutCBL << " k_c = " << k_c << endl;
    coutCBL << " M1_c = " << M1_c << " M_sun" << endl;
    coutCBL << " beta_c = " << beta_c << endl;
    coutCBL << " gamma_c = " << gamma_c << endl;
    coutCBL << " sigma_c  = " << sigma_c << endl;

    coutCBL << " k_s = " << k_s << endl;
    coutCBL << " M1_s = " << M1_s << " M_sun" << endl;
    coutCBL << " beta_s = " << beta_s << endl;
    coutCBL << " gamma_s = " << gamma_s << endl;
    coutCBL << " Phi_s = " << Phi_s << endl;
    coutCBL << " alpha_s = " << alpha_s << endl << endl;
  }

  
  // ----- Moster et al. 2010 -----

  else if (HOD_Type == HODType::_Moster10_)
    {
      // cout<<"Warning: threshold in Moster is minimum stellar mass!"<<endl;
      if (threshold < 0.)
	ErrorCBL("Warning: threshold in Moster is minimum stellar mass!", "Catalogue", "HODCatalogue.cpp");

      //for (int ii=0; ii<10; ii++) {
      //	if(parameter[ii]==-99){
      if (parameter[0] == -99)
	parameter[0] = 0.0297;   // 0.0267; //norm SHMR_c    //
      if (parameter[1] == -99)
	parameter[1] = pow(10., 11.9008);   // pow(10,11.9347);//break SHMR_c  //
      if (parameter[2] == -99)
	parameter[2] = 1.0757; //← w scatter// 1.0059; //beta_c SHMR_c  //
      if (parameter[3] == -99)
	parameter[3] = 0.6310; //← w scatter// 0.5611; //gamma_c SHMR_c   //
      if (parameter[4] == -99)
	{
	  //cout << par::col_green << "If sigma_c=-99 it will be setted using the halo mass..." << par::col_default << endl;
	  //parameter[ii]=sigma_c;
	}
      if (parameter[5] == -99)
	parameter[5] = 0.0198;   // 0.0186; // norma SHMR_s  //
      if (parameter[6] == -99)
	parameter[6] = pow(10., 12.0640);    //pow(10,12.1988); // break SHMR_s //
      if (parameter[7] == -99)
	parameter[7] = 0.8097;    //0.7817; // beta_s SHMR_s //
      if (parameter[8] == -99)
	parameter[8] = 0.6910; //← w scatter//0.7334; // gamma_s SHMR_s //
      if (parameter[9] == -99)
	{ // PHI_s

	  //cout << par::col_green << "If Phi_s= -99 it will be setted using the halo mass..." << par::col_default << endl;
	  //parameter[ii]=Phi_s;
	}
      if (parameter[10] == -99)
	{ // alpha_s
	  //cout << par::col_green << "If alpha_s=-99 it will be setted using the halo mass..." << par::col_default << endl;
	  // parameter[ii]=alpha_s;
	}

      starMass_threshold = threshold;

      k_c = parameter[0];
      M1_c = parameter[1];
      beta_c = parameter[2];
      gamma_c = parameter[3];
      sigma_c = parameter[4];
      //Moster satallite
      k_s = parameter[5];
      M1_s = parameter[6];
      beta_s = parameter[7];
      gamma_s = parameter[8];
      Phi_s = parameter[9];
      alpha_s = parameter[10];

      coutCBL << par::col_green << "Used parameters(-99->set default value from halo mass):" << par::col_default << endl;
      coutCBL << " k_c = " << k_c << endl;
      coutCBL << " M1_c = " << M1_c << " M_sun" << endl;
      coutCBL << " beta_c = " << beta_c << endl;
      coutCBL << " gamma_c = " << gamma_c << endl;
      coutCBL << " sigma_c  = " << sigma_c << endl;

      coutCBL << " k_s = " << k_s << endl;
      coutCBL << " M1_s = " << M1_s << " M_sun" << endl;
      coutCBL << " beta_s = " << beta_s << endl;
      coutCBL << " gamma_s = " << gamma_s << endl;
      coutCBL << " Phi_s = " << Phi_s << endl;
      coutCBL << " alpha_s = " << alpha_s << endl << endl;
    }

  else
    ErrorCBL("HOD_type not available...", "Catalogue", "HODCatalogue.cpp");


  // -----------------------------------------------------------
  // --------------- populate the halo catalogue ---------------
  // -----------------------------------------------------------
  
  coutCBL << par::col_green << "Creating HOD catalogue" << par::col_default << endl;
  // Vector which contain parameters to extract  the stellar mass for central galaxies
  vector<double> par_stellar_c = {k_c, M1_c, beta_c, gamma_c, sigma_c};

  // Conditional mass function Moster2010 for central galaxies, eq:4.3.6
  const cbl::distribution_func CMF_c = [&cosm, &fsigma_c, &starMass_threshold] (double x, const shared_ptr<void> modelInput, vector<double> par_stellar_c)
				       {
					 const vector<double> fix_par = *static_pointer_cast<std::vector<double>>(modelInput);
					 if (par_stellar_c[4] == -99)
					   par_stellar_c[4] = fsigma_c(fix_par[0]);
					 double m_c = 2. * fix_par[0] * par_stellar_c[0] * pow(pow(fix_par[0] / par_stellar_c[1], -par_stellar_c[2]) + pow(fix_par[0] / par_stellar_c[1], par_stellar_c[3]), -1);

					 return (1. / (par_stellar_c[4] * sqrt(2. * cbl::par::pi) * x * log(10))) * exp(-pow(log10(x / m_c) / (sqrt(2) * par_stellar_c[4]), 2.));
				       };

  // Moster et al 2010 distribution probability of stellar mass for satellite galaxies, IMF: Kroupa

  vector<double> par_stellar_s = {k_s, M1_s, beta_s, gamma_s, Phi_s, alpha_s};

  //Conditional mass function Moster2010 for satellite galaxies
  const cbl::distribution_func CMF_s = [&cosm, &fPhi_s, &falpha_s](double x, const shared_ptr<void> modelInput, vector<double> par_stellar_s)
				       {
					 const vector<double> fix_par = *static_pointer_cast<std::vector<double>>(modelInput);
					 if (par_stellar_s[4] == -99)
					   par_stellar_s[4] = fPhi_s(fix_par[0]);
					 if (par_stellar_s[5] == -99)
					   par_stellar_s[5] = falpha_s(fix_par[0]);
					 double m_s = 2. * fix_par[0] * par_stellar_s[0] * pow(pow(fix_par[0] / par_stellar_s[1], -par_stellar_s[2]) + pow(fix_par[0] / par_stellar_s[1], par_stellar_s[3]), -1.);

					 return (par_stellar_s[4] / m_s) * pow(x / m_s, par_stellar_s[5]) * exp(-pow((x / m_s), 2.));
				       };

  
  // -------------------------------------------------
  // --------------- average functions ---------------
  // -------------------------------------------------

  //Average numbers of galaxies from CMF for central and satellite galaxies Moster et. al 2010

  auto Average_c_Moster_2010 = [&cosm, &fsigma_c] (double x, double star_threshold, double k_c, double M1_c, double beta_c, double gamma_c, double sigma_c)
			       {
				 if (sigma_c == -99)
				   sigma_c = fsigma_c(x);
				 double m_c = 2. * x * k_c / (pow(x / M1_c, -beta_c) + pow(x / M1_c, gamma_c));

				 return 0.5 * (1 - erf(log10(star_threshold / m_c) / (sqrt(2) * sigma_c))); //threshold
			       };

  auto Average_s_Moster_2010 = [&fPhi_s, &falpha_s](double x, double star_threshold, double k_s, double M1_s, double beta_s, double gamma_s, double Phi_s, double alpha_s)
			       {
				 if (Phi_s == -99)
				   Phi_s = fPhi_s(x);
				 if (alpha_s == -99)
				   alpha_s = falpha_s(x);
				 double m_s = 2. * x * k_s * pow(pow(x / M1_s, -beta_s) + pow(x / M1_s, gamma_s), -1.);
				 double lower = pow(star_threshold / m_s, 2.);

				 return (Phi_s / 2.) * gsl_sf_gamma_inc(0.5 * alpha_s + 0.5, lower);
			       };

  
  //-------------------------------------------------------
  // --------------- sub-halo mass function ---------------
  //-------------------------------------------------------

  const double z = 0.;
  vector<double> parameter_NFW{z};

  
  // Giocoli et al. 2010/2011

  const double A = 9.33 * 1e-4;
  const double zres = sqrt(1 + z);
  const double alpha_SubHalo = -0.9;
  const double betha_SubHalo = -12.2715;

  vector<double> par_SubHalo{A, zres, alpha_SubHalo, betha_SubHalo};

  const cbl::distribution_func SubMass = [&] (double lx, const shared_ptr<void> modelInput, const vector<double> par_SubHalo)
					 {
					   const std::vector<double> fix_par_h = *static_pointer_cast<std::vector<double>>(modelInput);
					   double x = pow(10., lx);
					   return pow(x, par_SubHalo[2] + 1) * exp(par_SubHalo[3] * pow(x, 3));
					 };

  
  //-------------------------------------------------
  // --------------- galaxy positions ---------------
  //-------------------------------------------------

  const cbl::distribution_func NFW_profile = [&cosm] (double x, const shared_ptr<void> modelInput, const vector<double> parameter_NFW)
					     {
					       const vector<double> fix_par_h = *static_pointer_cast<std::vector<double>>(modelInput);
					       double r_vir = cosm.r_vir(fix_par_h[0], parameter_NFW[0]);
					       double c_vir = cosm.concentration_NFW_Duffy(fix_par_h[0], parameter_NFW[0]);
					       double r_s = r_vir / c_vir;
					       double A = pow(r_s + r_vir, 2) / pow(r_s, 2);

					       return A * x / (pow(1 + x * c_vir, 2));
					     };


  // -------------------------------------------
  // --------------- infall mass ---------------
  // -------------------------------------------
	
  auto InfallMass = [&cosm] (double x, const shared_ptr<void> modelInput, const vector<double> par_infall)
		    {
		      const vector<double> fix_par_h = *static_pointer_cast<std::vector<double>>(modelInput);
		      const double z = 0;
		      double M_infall = par_infall[0] * pow(pow(x / cosm.r_vir(fix_par_h[0], z), 2. / 3.), -1); //0.65!!

		      return M_infall;
		    };

  
  //-------------------------------------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------

  
  const double M_stellarmax = pow(10., 12.);

  // default values (refer to the GSL documentation for different choiches)
  const int limit_size = 1000;
  const int rule = 6;
  const double prec = 1.e-5;
  // set h of the input cosmology
  const double h = cosm.hh();

  const double x_NFW_min = 1.e-5; //r/r_vir
  const double x_NFW_max = 1.;	//r/r_vir
  //const double v_sat_min = -1e10; //km/s
  //const double v_sat_max = 1e10; //km/s

  const int N_halo = halo_catalogue.nObjects();

  vector<bool> occupied(N_halo, false);
  std::default_random_engine generator; //for poisson number
  std::vector<std::vector<std::shared_ptr<Object>>> support_tot(N_halo);

  // start the  parallelization
#pragma omp parallel num_threads(omp_get_max_threads())

  {
#pragma omp for schedule(dynamic)

    for (int i = 0; i < N_halo; ++i) {

      if (i > N_halo)
	ErrorCBL("number of iterations exeed number of halos, possible infinite loop", "Catalogue", "HODCatalogue.cpp");

#pragma omp critical // <--means: only a thread EXECUTE ONLY the NEXT LINE:
      coutCBL << "..." << int(double(i) / double(N_halo) * 100) << "% completed\r"; cout.flush();

      // Parallelization continue here:
			
      // vector used to support the parallelization

      std::vector<std::shared_ptr<Object>> support;
      const int seedtheta = 333 + i;
      const int seedphi = 8782 + i;
      const int seed = 34562 + i;

      // properties of halos
      double halo_mass = halo_catalogue.mass(i);
      const double halo_catalogue_x = halo_catalogue.xx(i);
      const double halo_catalogue_y = halo_catalogue.yy(i);
      const double halo_catalogue_z = halo_catalogue.zz(i);

      // pointer for functions
      const vector<double> vect{halo_mass / h}; // Msun unit
      const vector<double> vect_h{halo_mass};	  // Unit of h if the input catalogue is in unit of h

      auto fix_par = make_shared<vector<double>>(vect);	  // Used for Moster
      auto fix_par_h = make_shared<vector<double>>(vect_h); // Used for everything else

      // mass halo in solar unit
      const double halo_mass_wo_h = halo_mass / h;

      
      // ---------------------------------------------------------------------------------------------------
      // --------------- computing the average number of galaxy using function defined above ---------------
      // ---------------------------------------------------------------------------------------------------
      
      double Navg_c = 0.;
      if (HOD_Type == HODType::_Zehavi11_)
	Navg_c = Average_c_Zehavi_2011(halo_mass, logMmin, logsigma_c);
      if (HOD_Type == HODType::_Zehavi05_)
	Navg_c = Average_c_Zehavi_2005(halo_mass, M_min);
      if (HOD_Type == HODType::_Moster10_)
	Navg_c = Average_c_Moster_2010(halo_mass_wo_h, starMass_threshold, k_c, M1_c, beta_c, gamma_c, sigma_c);

      random::CustomDistributionRandomNumbers StellarMass_c(CMF_c, fix_par, par_stellar_c, seed, starMass_threshold, M_stellarmax);

      // ------------------------------------------------
      // --------------- central galaxies --------------- 
      // ------------------------------------------------
      
      // poiter for the central galaxies
      auto galaxy_c = std::make_shared<cbl::catalogue::Galaxy>();
      auto number_c = nearbyint(Navg_c);

      for (size_t j = 0; j < number_c; ++j) {
	occupied[i] = true;
	galaxy_c->set_xx(halo_catalogue_x);
	galaxy_c->set_yy(halo_catalogue_y);
	galaxy_c->set_zz(halo_catalogue_z);
	galaxy_c->set_ID(i);
	// galaxy_c -> set_IDHost(halo_catalogue_ID);
	galaxy_c->set_mstar(StellarMass_c());
	galaxy_c->set_massinfall(halo_mass); //used to write a file
	galaxy_c->set_mass(halo_mass);
	galaxy_c->set_galaxyTag(0);

	support_tot[i].emplace_back(galaxy_c);
      }

      
      // --------------------------------------------------
      // --------------- satellite galaxies ---------------
      // --------------------------------------------------
			
      if (occupied[i]) {
	double Navg_s = 0.;
	if (HOD_Type == HODType::_Zehavi11_)
	  Navg_s = Average_s_Zehavi_2011(halo_mass, M0, M1, alpha);
	else if (HOD_Type == HODType::_Zehavi05_)
	  Navg_s = Average_s_Zehavi_2005(halo_mass, M1, alpha);
	else if (HOD_Type == HODType::_Moster10_)
	  Navg_s = Average_s_Moster_2010(halo_mass_wo_h, starMass_threshold, k_s, M1_s, beta_s, gamma_s, Phi_s, alpha_s);

	random::CustomDistributionRandomNumbers StellarMass_s(CMF_s, fix_par, par_stellar_s, seed, starMass_threshold, M_stellarmax); //fix_par_s per massa sotto aloni
	random::CustomDistributionRandomNumbers Radius(NFW_profile, fix_par_h, parameter_NFW, seed, x_NFW_min, x_NFW_max);
	random::UniformRandomNumbers Theta(0., cbl::par::pi, seedtheta);
	random::UniformRandomNumbers Phi(0., 2. * cbl::par::pi, seedphi);

	std::poisson_distribution<int> Poisson_s(Navg_s);
	int number_s = Poisson_s(generator);

	// Extracting the stellar mass of satellite galaxies
	vector<double> stellarmass_s(number_s);
	for (int j = 0; j < number_s; ++j)
	  stellarmass_s[j] = StellarMass_s();

	//integration for the mass fraction
	if (number_s > 0) {

	  std::function<double(double)> shmf = [&par_SubHalo, &vect_h](const double lx)
					       {
						 double x = pow(10., lx);

						 return (x * vect_h[0]) * par_SubHalo[0] * par_SubHalo[1] * pow(x * vect_h[0], par_SubHalo[2]) * exp(par_SubHalo[3] * pow(x, 3)) * log(10);
					       };

	  double MinLogExtraction = log10(starMass_threshold / halo_mass);
	  //double f_min=log10(1.e10/halo_mass);
	  double f = cbl::wrapper::gsl::GSL_integrate_qag(shmf, -5, 0., prec, limit_size, rule);
	  double MaxLogExtraction = log10(f);

	  random::CustomDistributionRandomNumbers SubHaloMass(SubMass, fix_par_h, par_SubHalo, seed, MinLogExtraction, MaxLogExtraction);

	  const double a = 1.e-4;
	  double delta = 0.2 + a * number_s;
	  int retryCounter = 0;

	  vector<double> subhalomass = {};

	  
	  // if substructure is true, compute the right subhalo mass fraction

	  if (substructures) {
	    vector<double> temp_subtot;
	    double temp_sum;
	    int temp_count;
	    double comp_f;
	    int counter_s = 0;

	    while (true) {
	      retryCounter++;
	      temp_subtot = {};
	      temp_sum = 0;
	      temp_count = 0;

	      while (temp_sum < 1. || temp_count < number_s) {
		temp_subtot.emplace_back(pow(10., SubHaloMass()));
		temp_sum += temp_subtot[temp_subtot.size() - 1];
		temp_count++;
	      }

	      counter_s = number_s;
	      comp_f = 0;
	      for (int jj = 0; jj < number_s; jj++)
		comp_f += temp_subtot[jj];
	      if (comp_f < f + f * delta / 2.) {
		while (!(comp_f > f - f * delta / 2. && comp_f < f + f * delta / 2.)) {
		  counter_s++;
		  comp_f += temp_subtot[counter_s - 1];
		  if (comp_f > f + f * delta / 2.)
		    break;
		}

		if (comp_f > f - f * delta / 2. && comp_f < f + f * delta / 2.)
		  break;
	      }
	    }

	    for (int kk = 0; kk < counter_s; kk++)
	      subhalomass.emplace_back(temp_subtot[kk] * halo_mass);

	    // If you are in hurry we can save your time! But the substructures lose their meaning.
	  }
	  else
	    for (int j = 0; j < number_s; ++j)
	      subhalomass.emplace_back(pow(10, SubHaloMass()) * halo_mass);

	  double Sum = 0;
	  for (size_t j = 0; j < subhalomass.size(); j++)
	    Sum += subhalomass[j];

	  if (subhalomass.size() > 0) {
	    std::sort(subhalomass.begin(), subhalomass.end(), greater<double>());
	    std::sort(stellarmass_s.begin(), stellarmass_s.end(), greater<double>());

	    for (int j = 0; j < number_s; ++j) {
	      auto galaxy_s = std::make_shared<cbl::catalogue::Galaxy>();

	      vector<double> par_infall{subhalomass[j]};
	      const double radius = cosm.r_vir(halo_mass, z) * Radius();
	      const double theta = Theta();
	      const double phi = Phi();
	      double M_infall = InfallMass(radius, fix_par_h, par_infall);

	      galaxy_s->set_xx(halo_catalogue_x + radius * sin(theta) * cos(phi));
	      galaxy_s->set_yy(halo_catalogue_y + radius * sin(theta) * sin(phi));
	      galaxy_s->set_zz(halo_catalogue_z + radius * cos(theta));
	      galaxy_s->set_ID(i);
	      galaxy_s->set_mstar(stellarmass_s[j] * h);
	      galaxy_s->set_massinfall(M_infall);
	      galaxy_s->set_mass(subhalomass[j]);
	      galaxy_s->set_galaxyTag(1);

	      support_tot[i].emplace_back(galaxy_s);
	    }
	  }
	}
      }
    }
  }

  // End of parallelisation all the information are put toghether

  for (int i=0; i<N_halo; ++i) 
    for (size_t j = 0; j < support_tot[i].size(); j++)
      m_object.emplace_back(support_tot[i][j]);

  if (HOD_Type!=HODType::_Zehavi11_) {
    (void)(logMmin);
    (void)(logsigma_c);
    (void)(M0);
    (void)(M1);
    (void)(alpha);
  }

  if (HOD_Type!=HODType::_Zehavi05_) {
    (void)(M_min);
    (void)(M1);
    (void)(alpha);
  }

  cout.flush();
  coutCBL << "\n ...Done!" << endl;

  auto endTimer = chrono::steady_clock::now();

  float seconds = chrono::duration_cast<chrono::seconds>(endTimer - startTimer).count();


  coutCBL << "Time spent to compute: " << seconds << " seconds " << endl;
  //coutCBL << "Time spent to compute: " << seconds / 60. << " minutes " << endl;
  //coutCBL << "Time spent to compute: " << seconds / 3600. << " hours " << endl;
  coutCBL << "Writing time log file: timeLog.txt" << endl;
  cout << endl;

  // Create and open the timelog file
  ofstream timeLogFile("./output/timelog.txt");

  // Write to the file the timings
  timeLogFile << "Time spent to compute: " << seconds << " seconds " << endl;

  // Close the timeLog file
  timeLogFile.close();
  
}
