/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Carlo Giocoli        *
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
 *  @brief Methods of the class Cosmology used to model the cosmic
 *  mass accretion history
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model the cosmic mass accretion history
 *
 *  @authors Federico Marulli, Carlo Giocoli 
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Cosmology.h"
using namespace cosmobl;


// =====================================================================================

/* ======== Carlo Giocoli ======== */


// -------- Halo Mass Growth History based on Giocoli et al. 2012, Nusser and Sheth, Lacey and Coles 1993 -------- 

// Differential rescaled and generalized formation redshift distribution

double cosmobl::Cosmology::pw (const double ww, const double ff, const string author) const 
{
  if (author=="NS") {
    if (ff<0.5) cout <<"Warning you are calling pw function for NS with f = "<<ff<<endl;
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


// Probability that a halo of a given mass m0 at redshift z0 make a mass fraction f at redshift z

double cosmobl::Cosmology::pz (const double m0, const double z0, const double frac, const double redshift, const string author_model, const string method_SS, const string output_root) const
{
  double zero = 0.;
  double dcz0 = deltac(z0)/DD(z0)*DD(zero);
  double dcz = deltac(redshift)/DD(redshift)*DD(zero);
  double SS = SSM_norm(m0,method_SS,redshift,output_root); 
  double mf = m0*frac;
  double SSf = SSM_norm(mf,method_SS,redshift,output_root);
  double ww = (dcz-dcz0)/sqrt(SSf-SS);
  if (author_model=="NS"){
    if(frac<0.5) cout <<"Warning you are calling pw function for NS with frac = "<<frac<<endl;
    return (1./frac-1.)*2.*ww*erfc(ww/sqrt(2.))+(2.-1./frac)*sqrt(2./par::pi)*exp(-ww*ww*0.5);
  }
  if (author_model=="GTS"){
    double alpha = 0.815*exp(-2.*frac*frac*frac)/pow(frac,0.707);
    double Denumerator = (exp(ww*ww*0.5)+alpha-1.);
    return alpha*ww*exp(ww*ww*0.5)/Denumerator/Denumerator;
  }
  else return -1.;
}


// =====================================================================================


// Cumulative rescaled and generalized formation redshift distribution
 
double cosmobl::Cosmology::cumPw (const double ww, const double ff, const string author) const
{
  if(author=="NS"){
    if(ff<0.5) cout <<"Warning you are calling cumPw function for NS with f = "<<ff<<endl;
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


void cosmobl::Cosmology::medianwf (const double ff, const string author_model, vector<double> &wf) const
{
  wf.resize(3);

  if (author_model=="NS") {
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

  if (author_model=="GTS") {
    double alpha = 0.815*exp(-2*ff*ff*ff)/pow(ff, 0.707);
    wf[2] = sqrt(2.*log(4.*alpha-(alpha-1.)));
    wf[1] = sqrt(2.*log(alpha+1.));
    wf[0] = sqrt(2.*log(alpha/0.75-(alpha-1.)));
  }

}


// =====================================================================================

// Conditional variable w = \delta_c(zf) - \delta_c(z)/\sqrt(s(fm)-s(m))
// we recall that \delta_c(z) = delta_c0(z)/D+(z)

double cosmobl::Cosmology::wf (const double mm, const double redshift, const double ff, const double zf, const string method_SS, const string output_root) const
{
  double zero = 0.;
  double deltacz = deltac(redshift)/DD(redshift)*DD(zero);
  double deltaczf = deltac(zf)/DD(zf)*DD(zero);
  double SS = SSM_norm(mm,method_SS,redshift,output_root); 
  double mf = mm*ff;
  double SSf = SSM_norm(mf,method_SS,redshift,output_root); 
  return (deltaczf-deltacz)/sqrt(SSf-SS);
}


// =====================================================================================

// with this routine you can estimate the redshift from w given the parent halo mass (at z=z_0), z_0 and its assembled fraction f

double cosmobl::Cosmology::Redshift (const double mm, const double redshift, const double ff, const string method_SS, const double wwf, const string output_root) const
{
  int const nn = 128;
  vector<double> lzi = linear_bin_vector(nn, 0., 1.7);
  double zero = 0.;
  double dc0 = deltac(redshift)/DD(redshift)*DD(zero);
  double SS = SSM_norm(mm, method_SS, redshift, output_root); 
  double mf = mm*ff;
  double SSf = SSM_norm(mf, method_SS, redshift, output_root); 
  double dd = wwf*sqrt(SSf-SS) + dc0;
  
  vector<double> dci(nn);

  for (int i=0; i<nn; i++) {
    double zi = -1 + pow(10.,lzi[i]);
    dci[i] = deltac(zi)/DD(zi)*DD(zero);
  }
  
  return -1.+pow(10.,interpolated(dd, dci, lzi, "Poly"));
}


// =====================================================================================


void cosmobl::Cosmology::medianzf (const double ff, const double mass, const double z0, const string author_model, const string method_SS, vector<double> &zf, const string output_root) const
{
  vector<double> wf;
  zf.resize(3);
  medianwf(ff, author_model, wf);

  zf[0] = Redshift(mass, z0, ff, method_SS, wf[0], output_root);
  zf[1] = Redshift(mass, z0, ff, method_SS, wf[1], output_root);
  zf[2] = Redshift(mass, z0, ff, method_SS, wf[2], output_root);
}


// ============================================================================


double cosmobl::Cosmology::concentration (const double Vmax, const double Rmax) const
{ 
  const int nn = 128;
  vector<double> xxi = linear_bin_vector(nn, 0.1, 50.);
  vector<double> yyi(nn);

  for (int i=0; i<nn; i++)
    // reset 200 for the spherical collapse model
    yyi[i] = 200./3.*pow(xxi[i],3.)/(log(1.+xxi[i])-xxi[i]/(1.+xxi[i]))-14.426*pow(Vmax/Rmax/m_H0,2.);
  
  double null = 0.;
  return interpolated(null, yyi, xxi, "Poly");
}

