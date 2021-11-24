/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Cosmology/Lib/Bias.cpp
 *
 *  @brief Methods of the class Cosmology used to model the bias
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model the bias
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#include "Cosmology.h"

using namespace std;

using namespace cbl;


// =====================================================================================


vector<double> cbl::cosmology::Cosmology::bias_halo (const vector<double> Mass, const vector<double> Sigma, const double redshift, const std::string model_bias, const bool store_output, const std::string output_root, const std::string interpType, const double Delta, const double kk, const int norm, const double k_min, const double k_max, const double prec, const std::string method_SS, const std::string input_file, const bool is_parameter_file) 
{
  double D_N = DN(redshift);
  vector<double> bias(Mass.size());
  for (size_t i=0; i<Mass.size(); i++)
      bias[i] = m_bias_halo_generator(Sigma[i], redshift, D_N, model_bias, Delta); 

  if (m_fNL!=0)
      for (size_t i=0; i<Mass.size(); i++)
          bias[i] += bias_correction(kk, Mass[i], method_SS, store_output, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file)*Sigma[i]*pow(bias[i]-1, 2); // check!!!

  return bias; 
}


// =====================================================================================


double cbl::cosmology::Cosmology::m_bias_halo_generator (const double Sigma, const double redshift, const double D_N, const std::string author, const double Delta) const
{
  const double deltacz = deltac(redshift);
  const double sigmaz = Sigma*D_N;
  
  double bias = -1000.;

  if (author=="ST99") {
    double aa = 0.707;
    double pp = 0.3;
    double ni = pow(deltacz/sigmaz, 2); 
    bias = 1.+(aa*ni-1.)/deltacz+(2.*pp/deltacz)/(1.+pow(aa*ni,pp));
  }

  else if (author=="SMT01") {
    double aa = 0.707;
    double bb = 0.5;
    double cc = 0.6;
    double ni = deltacz/sigmaz; 
    bias = 1.+1./(sqrt(aa)*deltacz)*(sqrt(aa)*aa*pow(ni,2.)+sqrt(aa)*bb*pow(aa*pow(ni,2.),1.-cc)-pow(aa*pow(ni,2.),cc)/(pow(aa*pow(ni,2.),cc)+bb*(1.-cc)*(1.-cc*0.5)));
  }
  
  else if (author=="SMT01_WL04") {
    double aa = 0.707;
    double bb = 0.5;
    double cc = 0.6;
    double ni = deltacz/sigmaz; 
    double niI = sqrt(aa)*ni;
    bias = 1.+1./deltacz*(pow(niI,2.)+bb*pow(niI,2.*(1.-cc))-pow(niI,2.*cc)/sqrt(aa)/(pow(niI,2.*cc)+bb*(1.-cc)*(1.-cc*0.5)));
  }
  
  else if (author=="Tinker") { // Tinker et al. (2010)
    double yy = log10(Delta);
    double AA = 1.+0.24*yy*exp(-pow(4./yy,4));
    double aa = 0.44*yy-0.88;
    double BB = 0.183;
    double bb = 1.5;
    double CC = 0.019+0.107*yy+0.19*exp(-pow(4./yy,4));
    double ccc = 2.4;
    double ni = 1.686/sigmaz;
    bias = 1.-AA*pow(ni,aa)/(pow(ni,aa)+pow(1.686,aa))+BB*pow(ni,bb)+CC*pow(ni,ccc);
  }
  
  else
    ErrorCBL("author = " + author + "!", "m_bias_halo_generator", "Bias.cpp");
  
  return bias;
}

