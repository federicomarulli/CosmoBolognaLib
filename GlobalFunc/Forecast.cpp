/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file GlobalFunc/Forecast.cpp
 *
 *  @brief Functions to generate mock 2PCF/3PCF signal 
 *  and provide forecasts
 *
 *  This file contains the implementation of a set of functions to
 *  generate mock 2PCF/3PCF signal and provide forecasts 
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Data1D.h"
#include "GlobalFunc.h"

using namespace std;

using namespace cbl;


// ============================================================================


std::vector<double> cbl::fit_covariance_matrix_2PCF_monopole (const std::vector<double> mean, const std::vector<std::vector<double>> mock_xi0, const bool doJK, const cbl::cosmology::Cosmology cosmology, const double nObjects, const double Volume, const double bias, const double redshift, const double rMin, const double rMax, const int nbins, const cbl::BinType bin_type, const std::string method_Pk, const double sigma_NL, const bool NL)
{
  coutCBL << "Fitting the Covariance Matrix" << endl;
  const int nmock = mock_xi0.size();

  double term = (doJK) ? nmock-1 : 1.;
  vector<double> mock_mean(nbins, 0);

  if (mean.size()==0) {
    vector<vector<double>> vv = transpose(mock_xi0);
    for (int i=0; i<nbins; i++)
      mock_mean[i] = Average(vv[i]);
  }
  else
    for (int i=0; i<nbins; i++)
      mock_mean[i] = mean[i];

  vector<vector<double>> diff(nmock, vector<double>(nbins, 0));
  for (int i=0; i<nmock; i++)
    for (int j=0; j<nbins; j++)
      diff[i][j] = mock_xi0[i][j]-mock_mean[j];

  cbl::cosmology::Cosmology cosmo = cosmology;
  vector<double> kk = cbl::logarithmic_bin_vector(200, 1.e-4, 1.e2);
  vector<double> Pk0(kk.size()), Pk2(kk.size()), Pk4(kk.size());
  double linear_growth_rate = cosmo.linear_growth_rate(redshift, 1);
  double beta = linear_growth_rate/bias;
  
  double pk0_fact = bias*bias*(1+2.*beta/3+beta*beta/5);
  double pk2_fact = 4./3*beta+4./7*beta*beta;
  double pk4_fact = 8./35*beta+beta;

  if (sigma_NL==0 or NL==true) {
    for (size_t i=0; i<kk.size(); i++) {
      double pkDM = cosmo.Pk(kk[i], method_Pk, NL, redshift);
      Pk0[i] = pkDM*pk0_fact;
      Pk2[i] = pkDM*pk2_fact;
      Pk4[i] = pkDM*pk4_fact;
    }
  }
  else {
    vector<double> PkDM = cosmo.Pk_DM_DeWiggled(method_Pk, "EisensteinHu", kk, redshift, sigma_NL);
    for (size_t i=0; i<kk.size(); i++) {
      Pk0[i] = PkDM[i]*pk0_fact;
      Pk2[i] = PkDM[i]*pk2_fact;
      Pk4[i] = PkDM[i]*pk4_fact;
    }
  }

  vector<double> nObj = linear_bin_vector(20, 0.2*nObjects, nObjects);
  vector<double> Vol = linear_bin_vector(20, 0.2*Volume, Volume);

  auto func = [&] (vector<double> &parameters)
  {
    const double _nObjects = parameters[0];
    const double _Volume = parameters[1];
    if (_nObjects<0.2*nObjects || _nObjects>nObjects ) return 1.e30;
    if (_Volume<0.2*Volume || _Volume>Volume ) return 1.e30;

    vector<double> rr;
    vector<vector<double>> covariance, icovariance;
    cbl::Covariance_XiMultipoles (rr, covariance, nbins, rMin, rMax, _nObjects, _Volume, kk, {Pk0, Pk2, Pk4}, {0}, bin_type);
    for (int i=0; i<nbins; i++)
      for (int j=0; j<nbins; j++)
	covariance[i][j] /=term;

    invert_matrix(covariance, icovariance);
    double chi2 = 0.;
    for (int i=0; i<nmock; i++)
      chi2+=v_M_vt(diff[i], icovariance);
    double LL = nmock*nbins*nbins*log(2*par::pi)+nmock*log(determinant_matrix(covariance))+chi2;
    return LL;
  };
  
  double min_nObj = 0., min_Vol = 0.;
  double max_nObj = 50*nObjects, max_Vol = 50*Volume;
  double min_LL = 1.e30;
  for (size_t i=0; i<nObj.size(); i++)
    for (size_t j=0; j<Vol.size(); j++) {
      vector<double> pars = {nObj[i], Vol[j]};
      double LL = func(pars);
      if (LL<min_LL) {
	min_LL = LL;
	min_nObj = nObj[i];
	min_Vol = Vol[j];
      }
    }

  vector<double> start = {nObjects, Volume};
  vector<vector<double>> limits = { {min_nObj, max_nObj}, {min_Vol, max_Vol}};
  return cbl::wrapper::gsl::GSL_minimize_nD(func, start, limits);
}


// ============================================================================


std::shared_ptr<cbl::data::Data> cbl::generate_mock_2PCF_monopole (const cbl::cosmology::Cosmology cosmology, const double bias, const double nObjects, const double Volume, const double redshift, const double rMin, const double rMax, const int nbins, const cbl::BinType bin_type, const std::string method_Pk, const double sigma_NL, const bool NL)
{
  auto data = cbl::generate_mock_2PCF_multipoles(cosmology, bias, nObjects, Volume, redshift, rMin, rMax, nbins, bin_type, method_Pk, sigma_NL, NL);

  int ndata = data->ndata()/3;

  vector<double> rad(ndata, 0), xi0(ndata, 0);
  vector<vector<double>> covariance(ndata, vector<double>(ndata, 0));

  for (int i=0; i<ndata; i++) {
    rad[i] = data->xx(i);
    xi0[i] = data->data(i);
    for (int j=0; j<ndata; j++)
      covariance[i][j] = data->covariance(i, j);
  }

  return make_shared<cbl::data::Data1D> (cbl::data::Data1D(rad, xi0, covariance));
}

// ============================================================================


std::shared_ptr<cbl::data::Data> cbl::generate_mock_2PCF_multipoles (const cbl::cosmology::Cosmology cosmology, const double bias, const double nObjects, const double Volume, const double redshift, const double rMin, const double rMax, const int nbins, const cbl::BinType bin_type, const std::string method_Pk, const double sigma_NL, const bool NL)
{
  cbl::cosmology::Cosmology cosmo = cosmology;
  vector<double> kk = cbl::logarithmic_bin_vector(200, 1.e-4, 1.e2);
  vector<double> Pk0(kk.size()), Pk2(kk.size()), Pk4(kk.size());
  double linear_growth_rate = cosmo.linear_growth_rate(redshift, 1);
  double beta = linear_growth_rate/bias;
  
  double pk0_fact = bias*bias*(1+2.*beta/3+beta*beta/5);
  double pk2_fact = 4./3*beta+4./7*beta*beta;
  double pk4_fact = 8./35*beta+beta;

  if (sigma_NL==0 or NL==true) {
    for (size_t i=0; i<kk.size(); i++) {
      double pkDM = cosmo.Pk(kk[i], method_Pk, NL, redshift);
      Pk0[i] = pkDM*pk0_fact;
      Pk2[i] = pkDM*pk2_fact;
      Pk4[i] = pkDM*pk4_fact;
    }
  }
  else {
    vector<double> PkDM = cosmo.Pk_DM_DeWiggled(method_Pk, "EisensteinHu", kk, redshift, sigma_NL);
    for (size_t i=0; i<kk.size(); i++) {
      Pk0[i] = PkDM[i]*pk0_fact;
      Pk2[i] = PkDM[i]*pk2_fact;
      Pk4[i] = PkDM[i]*pk4_fact;
    }
  }


   vector<double> rr;
   vector<vector<double>> xil(3), covariance;

   cbl::Covariance_XiMultipoles (rr, covariance, nbins, rMin, rMax, nObjects, Volume, kk, {Pk0, Pk2, Pk4}, {0, 2, 4}, bin_type);

  xil[0] = cbl::wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk0, 0);
  xil[1] = cbl::wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk2, 2);
  xil[2] = cbl::wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk4, 4);

  vector<double> rad, xi;
  for (int i=0; i<3; i++) {
    for (int j=0; j<nbins; j++) {
      rad.push_back(rr[j]);
      xi.push_back(xil[i][j]);
    }
  }

  return make_shared<cbl::data::Data1D>(cbl::data::Data1D(rad, xi, covariance));
}

