/********************************************************************
 *  Copyright (C) 2023 by Federico Marulli                          *
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
 *  @file Cosmology/Lib/MassGrowth.cpp
 *
 *  @brief Methods of the class MassGrowth used to model the cosmic
 *  mass accretion history
 *
 *  This file contains the implementation of the methods of the class
 *  MassGrowth used to model the cosmic mass accretion history
 *
 *  @authors Federico Marulli, Carlo Giocoli 
 *
 *  @author federico.marulli3@unibo.it, carlo.giocoli@inaf.it
 */

#include "MassGrowth.h"

using namespace std;

using namespace cbl;


// =====================================================================================


double cbl::cosmology::MassGrowth::pw (const double ww, const double ff, const std::string author) const 
{
  if (author=="NS") {
    if (ff<0.5) coutCBL <<"Warning you are calling pw function for NS with f = "<<ff<<endl;
    return (1./ff-1.)*2.*ww*erfc(ww/sqrt(2.))+(2.-1./ff)*sqrt(2./par::pi)*exp(-ww*ww*0.5);
  }
  if (author=="GTS") {
    double alpha = 0.815*exp(-2.*ff*ff*ff)/pow(ff,0.707);
    double Denumerator = (exp(ww*ww*0.5)+alpha-1.);
    return alpha*ww*exp(ww*ww*0.5)/Denumerator/Denumerator;
  }
  else return -1.;
}


// =====================================================================================


double cbl::cosmology::MassGrowth::pz (const double m0, const double z0, const double frac, const double redshift, const std::string model_model, const std::string method_SS, const bool store_output, const std::string output_root) const
{
  double dcz0 = m_cosmology->deltac(z0)/m_cosmology->DN(z0);
  double dcz = m_cosmology->deltac(redshift)/m_cosmology->DN(redshift);

  Sigma SG(m_cosmology);
  
  double SS = SG.sigma2M({m0}, method_SS, redshift, store_output, output_root)[0]; 
  double mf = m0*frac;
  double SSf = SG.sigma2M({mf}, method_SS, redshift, store_output, output_root)[0];
  double ww = (dcz-dcz0)/sqrt(SSf-SS);
  if (model_model=="NS"){
    if(frac<0.5) coutCBL <<"Warning you are calling pw function for NS with frac = "<<frac<<endl;
    return (1./frac-1.)*2.*ww*erfc(ww/sqrt(2.))+(2.-1./frac)*sqrt(2./par::pi)*exp(-ww*ww*0.5);
  }
  if (model_model=="GTS"){
    double alpha = 0.815*exp(-2.*frac*frac*frac)/pow(frac,0.707);
    double Denumerator = (exp(ww*ww*0.5)+alpha-1.);
    return alpha*ww*exp(ww*ww*0.5)/Denumerator/Denumerator;
  }
  else return -1.;
}


// =====================================================================================

 
double cbl::cosmology::MassGrowth::cumPw (const double ww, const double ff, const std::string author) const
{
  if(author=="NS"){
    if(ff<0.5) coutCBL <<"Warning you are calling cumPw function for NS with f = "<<ff<<endl;
    double F0 = erf(ww/sqrt(2.))+ww*ww*erfc(ww/sqrt(2.))-sqrt(2./par::pi)*ww*exp(-ww*ww*0.5);
    double F1 = (1./ff-1.)-(1./ff-1.)*F0;
    double F2 = (2.-1./ff)*(1.-erf(ww/sqrt(2.)));
    return F1+F2;
  }
  if(author=="GTS"){
    double alpha = 0.815*exp(-2.*ff*ff*ff)/pow(ff,0.707);
    double Denumerator = (exp(0.5*ww*ww)+alpha-1.);
    return alpha/Denumerator;
  }
  else return -1.;
}


// =====================================================================================


void cbl::cosmology::MassGrowth::medianwf (const double ff, const std::string model_model, vector<double> &wf) const
{
  wf.resize(3);

  if (model_model=="NS") {
    int nn = 128;
    vector<double> ww = linear_bin_vector(nn, 0., 5.);
    vector<double> Pw(nn);
    for(int i=0; i<nn; i++) {
      double F0 = erf(ww[i]/sqrt(2.))+ww[i]*ww[i]*erfc(ww[i]/sqrt(2.))-sqrt(2./par::pi)*ww[i]*exp(-ww[i]*ww[i]*0.5);
      double F1 = (1./ff-1.)-(1./ff-1.)*F0;
      double F2 = (2.-1./ff)*(1.-erf(ww[i]/sqrt(2.)));
      Pw[i] = F1+F2;
    }
    string type = "Poly";
    wf[2] = interpolated(0.5, Pw, ww, type);
    wf[1] = interpolated(0.25, Pw, ww, type);
    wf[0] = interpolated(0.75, Pw, ww, type);
  }

  if (model_model=="GTS") {
    double alpha = 0.815*exp(-2*ff*ff*ff)/pow(ff, 0.707);
    wf[2] = sqrt(2.*log(4.*alpha-(alpha-1.)));
    wf[1] = sqrt(2.*log(alpha+1.));
    wf[0] = sqrt(2.*log(alpha/0.75-(alpha-1.)));
  }

}


// =====================================================================================


double cbl::cosmology::MassGrowth::wf (const double mm, const double redshift, const double ff, const double zf, const std::string method_SS, const bool store_output, const std::string output_root) const
{
  Sigma SG(m_cosmology);
  
  double deltacz = m_cosmology->deltac(redshift)/m_cosmology->DN(redshift);
  double deltaczf = m_cosmology->deltac(zf)/m_cosmology->DN(zf);
  double SS = SG.sigma2M({mm}, method_SS, redshift, store_output, output_root)[0]; 
  double mf = mm*ff;
  double SSf = SG.sigma2M({mf}, method_SS, redshift, store_output, output_root)[0]; 
  return (deltaczf-deltacz)/sqrt(SSf-SS);
}


// =====================================================================================


double cbl::cosmology::MassGrowth::Redshift (const double mm, const double redshift, const double ff, const std::string method_SS, const double wwf, const bool store_output, const std::string output_root) const
{
  Sigma SG(m_cosmology);
  
  int const nn = 128;
  vector<double> lzi = linear_bin_vector(nn, 0., 1.7);
  double dc0 = m_cosmology->deltac(redshift)/m_cosmology->DN(redshift);
  double SS = SG.sigma2M({mm}, method_SS, redshift, store_output, output_root)[0]; 
  double mf = mm*ff;
  double SSf = SG.sigma2M({mf}, method_SS, redshift, store_output, output_root)[0]; 
  double dd = wwf*sqrt(SSf-SS) + dc0;
  
  vector<double> dci(nn);

  for (int i=0; i<nn; i++) {
    double zi = -1 + pow(10.,lzi[i]);
    dci[i] = m_cosmology->deltac(zi)/m_cosmology->DN(zi);
  }
  
  return -1.+pow(10.,interpolated(dd, dci, lzi, "Poly"));
}


// =====================================================================================


void cbl::cosmology::MassGrowth::medianzf (const double ff, const double mass, const double z0, const std::string model_model, const std::string method_SS, std::vector<double> &zf, const bool store_output, const std::string output_root) const
{
  vector<double> wf;
  zf.resize(3);
  medianwf(ff, model_model, wf);

  zf[0] = Redshift(mass, z0, ff, method_SS, wf[0], store_output, output_root);
  zf[1] = Redshift(mass, z0, ff, method_SS, wf[1], store_output, output_root);
  zf[2] = Redshift(mass, z0, ff, method_SS, wf[2], store_output, output_root);
}

