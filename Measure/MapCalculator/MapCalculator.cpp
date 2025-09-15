/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli                          *
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
 *  @file MapCalculator.cpp
 *
 *  @brief Methods of the class MapCalculator 
 *
 *  This file contains the implementation of the methods of the class
 *  MapCalculator, used to produce sky maps with HEALPix
 *
 *  @author Matteo Santini
 *
 *  @author matteo.santini7@studio.unibo.it
 */
 
#include "MapCalculator.h"
#include "RandomNumbers.h"

std::vector<int> cbl::measure::MapCalculator::HealpixMapping () const
{
  Healpix_Base healpixBase(m_ns, m_scheme, SET_NSIDE);  // RING ordering scheme as default
  std::vector<int> ipx1;
  auto size = m_data->nObjects();
  //need to convert from RA,DEC into polar spherical coordinates
  for (size_t i=0; i<size; ++i) {
    //the angles must be between[0, pi[ and [0, 2pi[
    pointing ang((cbl::par::pi/2)-m_data->dec(i), m_data->ra(i));
    ipx1.push_back(healpixBase.ang2pix(ang));
  }
  return ipx1;
}

std::vector<std::vector<std::vector<double>>> cbl::measure::MapCalculator::ADF (std::vector<double> const& mskin) const
{

  std::vector<bool> cnzero(mskin.size()); // For the mask application
  
  std::transform(mskin.begin(), mskin.end(), cnzero.begin(), [](double x) 
  { 
    return x>0.0; 
  });
  // Now cnzero is a boolean array (0 and 1), values depending on the matching of the condition in the lambda function 

  const unsigned int nzobs = m_zobs.size();
  const unsigned int nsgz = m_sgzar.size();
  const unsigned int size = m_data->nObjects();

  auto ipx1 = HealpixMapping();

  std::vector<std::vector<std::vector<double>>> map_ADF (nzobs, std::vector<std::vector<double>>(nsgz, std::vector<double>(12 * m_ns * m_ns, 0.0)));
  std::vector<std::vector<std::vector<double>>> nocounts (nzobs, std::vector<std::vector<double>>(nsgz, std::vector<double>(12 * m_ns * m_ns, 0.0)));
  std::vector<std::vector<double>> nang (nzobs, std::vector<double>(nsgz, 0.0));
    
  for (size_t iz=0; iz<nzobs; ++iz) {
    for (size_t isgz=0; isgz<nsgz; ++isgz) {
      std::vector<int> ipxA;
      std::vector<double> W1a;
      std::vector<int> true_inda;
      std::vector<bool> csel1a(size, 0);
      for (size_t i=0; i<ipx1.size(); ++i) {                
        csel1a[i] = ((std::abs(m_data->redshift(i) - m_zobs[iz]) < (m_nsg * m_sgzar[isgz])) && (mskin[ipx1[i]] != 0));
        if (csel1a[i]==true) {
          true_inda.push_back(i);
        }
      }
      // m_data->var(cbl::catalogue::Var::_Redshift_) causes memory leak ....
      std::vector<double> red;
      for (size_t ind=0; ind<size; ++ind)
        red.push_back(m_data->redshift(ind));
      auto zin_sel = cbl::select<double>(red, true_inda);
      for (size_t i=0; i<true_inda.size(); ++i) 
        W1a.push_back(std::exp(-0.5*(std::pow((zin_sel[i]-m_zobs[iz]), 2)/std::pow(m_sgzar[isgz], 2))));
      auto ipx_sel = select<int>(ipx1, true_inda);
      for (size_t i=0; i<true_inda.size(); ++i)  
        ipxA.push_back(ipx_sel[i]);
      for (size_t ig=0; ig<true_inda.size(); ++ig) 
        map_ADF[iz][isgz][ipxA[ig]] += W1a[ig];
      for (size_t i=0; i<mskin.size(); ++i) { 
        if (cnzero[i]==true) {
          map_ADF[iz][isgz][i] = map_ADF[iz][isgz][i]/mskin[i];
          nocounts[iz][isgz][i] = map_ADF[iz][isgz][i];
          nang[iz][isgz] += map_ADF[iz][isgz][i];
        }
      }

      auto iter_n = std::next(cnzero.begin(), 12*m_ns*m_ns); //any difference?
      nang[iz][isgz] = nang[iz][isgz]/std::accumulate(cnzero.begin(), iter_n, 0.0);
      double meanADF = 0;
      for (size_t i=0; i<mskin.size(); ++i) {
	if (cnzero[i]==true) { 
	  map_ADF[iz][isgz][i] = map_ADF[iz][isgz][i]/nang[iz][isgz];
	  meanADF += map_ADF[iz][isgz][i];
	}
      }
      meanADF = meanADF/std::accumulate(cnzero.begin(), iter_n, 0.0);
      for (size_t i=0; i<mskin.size(); ++i) {
	if (cnzero[i]==true) 
	  map_ADF[iz][isgz][i] = map_ADF[iz][isgz][i] - meanADF;
      }
    }
  }
  return map_ADF;
}

std::tuple<std::vector<std::vector<std::vector<double>>> ,std::vector<std::vector<std::vector<double>>>> cbl::measure::MapCalculator::ADF_ARF (std::vector<double> const& mskin) const
{

  std::vector<bool> cnzero(mskin.size()); // For the mask application
    
  std::transform(mskin.begin(), mskin.end(), cnzero.begin(), [](double x) 
  { 
    return x>0.0; 
  });
  // Now cnzero is a boolean array (0 and 1), values depending on the matching of the condition in the lambda function 

  const unsigned int nzobs = m_zobs.size();
  const unsigned int nsgz = m_sgzar.size();
  const unsigned int size = m_data->nObjects();

  auto ipx1 = HealpixMapping();
    
  std::vector<std::vector<std::vector<double>>> map_ADF (nzobs, std::vector<std::vector<double>>(nsgz, std::vector<double>(12 * m_ns * m_ns, 0.0)));
  std::vector<std::vector<std::vector<double>>> nocounts (nzobs, std::vector<std::vector<double>>(nsgz, std::vector<double>(12 * m_ns * m_ns, 0.0)));
  std::vector<std::vector<std::vector<double>>> map_ARF (nzobs, std::vector<std::vector<double>>(nsgz, std::vector<double>(12 * m_ns * m_ns, 0.0)));
  std::vector<std::vector<double>> nang (nzobs, std::vector<double>(nsgz, 0.0));

  for (size_t iz=0; iz<nzobs; ++iz) {
    for (size_t isgz=0; isgz<nsgz; ++isgz) {
      std::vector<int> ipxA;
      std::vector<double> W1a;
      std::vector<int> true_inda;
      std::vector<bool> csel1a(size, 0);
      std::vector<double> zeff;
      std::vector<double> zav1a_vec;
      std::vector<bool> cnzero2(12*m_ns*m_ns, false);
      for (size_t i=0; i<ipx1.size(); ++i) {                
        csel1a[i] = ((std::abs(m_data->redshift(i)-m_zobs[iz])<(m_nsg*m_sgzar[isgz])) && (mskin[ipx1[i]]!=0));
        if (csel1a[i]==true) 
          true_inda.push_back(i);
      }
      // m_data->var(cbl::catalogue::Var::_Redshift_) causes memory leak ....
      std::vector<double> red;
      for (size_t ind=0; ind<size; ++ind)
        red.push_back(m_data->redshift(ind));
      auto zin_sel = cbl::select<double>(red, true_inda);
      for (size_t i=0; i<true_inda.size(); ++i)  
        W1a.push_back(std::exp(-0.5*(std::pow((zin_sel[i]-m_zobs[iz]), 2)/std::pow(m_sgzar[isgz], 2))));
      auto ipx_sel = select<int>(ipx1, true_inda);
      for (size_t i=0; i<true_inda.size(); ++i) {  
        ipxA.push_back(ipx_sel[i]);
        zav1a_vec.push_back((W1a[i]*zin_sel[i]));
        zeff.push_back(zin_sel[i]);
      }
            
      double zav1a = std::accumulate(zav1a_vec.begin(), zav1a_vec.end(), 0.0)/std::accumulate(W1a.begin(), W1a.end(), 0.0);
      for (size_t ig=0; ig<true_inda.size(); ++ig) {
        map_ADF[iz][isgz][ipxA[ig]] += W1a[ig];
        map_ARF[iz][isgz][ipxA[ig]] += W1a[ig]*(zeff[ig]-zav1a);
      }
      for (size_t i=0; i<mskin.size(); ++i) 
        cnzero2[i] = (map_ADF[iz][isgz][i]!=0); // assigning bools
      for (size_t i=0; i<mskin.size(); ++i) { 
        if (cnzero[i]==true) {
          map_ADF[iz][isgz][i] = map_ADF[iz][isgz][i]/mskin[i];
          nocounts[iz][isgz][i] = map_ADF[iz][isgz][i];
          nang[iz][isgz] += map_ADF[iz][isgz][i];
        }
        if (cnzero2[i]==true) 
          map_ARF[iz][isgz][i] = map_ARF[iz][isgz][i]/mskin[i];
      }

      auto iter_n = std::next(cnzero.begin(), 12*m_ns*m_ns);
      nang[iz][isgz] = nang[iz][isgz]/std::accumulate(cnzero.begin(), iter_n, 0.0);
      double meanADF = 0;
      double meanARF = 0;
      for (size_t i=0; i<mskin.size(); ++i) {
        if (cnzero[i]==true) { 
          map_ADF[iz][isgz][i] = map_ADF[iz][isgz][i]/nang[iz][isgz];
          meanADF += map_ADF[iz][isgz][i];
        }
        if (cnzero2[i]==true) { 
          map_ARF[iz][isgz][i] = map_ARF[iz][isgz][i]/nang[iz][isgz];
          meanARF += map_ARF[iz][isgz][i];
        }
      }
      meanADF = meanADF/std::accumulate(cnzero.begin(), iter_n, 0.0);
      auto iter_n2 = std::next(cnzero2.begin(), 12*m_ns*m_ns); //any difference?
      meanARF = meanARF/std::accumulate(cnzero2.begin(), iter_n2, 0.0);
      for (size_t i=0; i<mskin.size(); ++i) {
        if (cnzero[i]==true) 
          map_ADF[iz][isgz][i] = map_ADF[iz][isgz][i]-meanADF;
        if (cnzero2[i]==true) 
          map_ARF[iz][isgz][i]=map_ARF[iz][isgz][i]-meanARF;
      }
    }
  }
  return std::make_tuple(map_ADF, map_ARF);
}

std::vector<std::vector<std::vector<double>>> cbl::measure::MapCalculator::ADFr (std::vector<double> const& mskin, std::vector<double> const& zH) const 
{
  std::vector<bool> cnzero(mskin.size()); // For the mask application
    
  std::transform(mskin.begin(), mskin.end(), cnzero.begin(), [](double x) 
  { 
    return x>0.0; 
  });
  // Now cnzero is a boolean array (0 and 1), values depending on the matching of the condition in the lambda function 

  const unsigned int nzobs = m_zobs.size();
  const unsigned int nsgz = m_sgzar.size();
  const unsigned int size = m_data->nObjects();

  auto ipx1 = HealpixMapping();

  std::vector<std::vector<std::vector<double>>> map_ADFr (nzobs, std::vector<std::vector<double>>(nsgz, std::vector<double>(12 * m_ns * m_ns, 0.0)));
  std::vector<std::vector<double>> nangr (nzobs, std::vector<double>(nsgz, 0.0));

  for (size_t iz=0; iz<nzobs; ++iz) {
    for (size_t isgz=0; isgz<nsgz; ++isgz) {
      std::vector<int> ipxA;
      std::vector<double> W1b;
      std::vector<int> true_indb;
      std::vector<bool> csel1b(size, 0);
      std::vector<double> zeff;
      std::vector<bool> cnzero2(12*m_ns*m_ns, false);
      for (size_t i=0; i<ipx1.size(); ++i) {                
        csel1b[i] = (std::abs(zH[i]-m_zobs[iz])<m_nsg*m_sgzar[isgz]) && (mskin[ipx1[i]]!=0);
        if (csel1b[i]==true) 
          true_indb.push_back(i);
      }
      std::vector<double> red;
      for (size_t ind=0; ind<size; ++ind)
        red.push_back(m_data->redshift(ind));
      auto zH_sel = cbl::select<double>(red, true_indb);
      for (size_t i=0; i<true_indb.size(); ++i) 
        W1b.push_back(std::exp(-0.5*(std::pow((zH_sel[i]-m_zobs[iz]), 2)/std::pow(m_sgzar[isgz], 2))));
      auto ipx_sel = select<int>(ipx1, true_indb);
      for (size_t i=0; i<true_indb.size(); ++i)   
        ipxA.push_back(ipx_sel[i]);
      for (size_t ig=0; ig<true_indb.size(); ++ig) 
        map_ADFr[iz][isgz][ipxA[ig]] += W1b[ig];
      for (size_t i=0; i<mskin.size(); ++i) 
        cnzero2[i] = (map_ADFr[iz][isgz][i]!=0); // assigning bools
      for (size_t i=0; i<mskin.size(); ++i) { 
        if (cnzero[i]==true) {
          map_ADFr[iz][isgz][i] = map_ADFr[iz][isgz][i]/mskin[i];
          nangr[iz][isgz] += map_ADFr[iz][isgz][i];
        }
      }

      auto iter_n = std::next(cnzero.begin(), 12*m_ns*m_ns); 
      nangr[iz][isgz] = nangr[iz][isgz]/std::accumulate(cnzero.begin(), iter_n, 0.0);
      double meanADFr = 0;
      for (size_t i=0; i<mskin.size(); ++i) {
        if (cnzero[i]==true) { 
          map_ADFr[iz][isgz][i] = map_ADFr[iz][isgz][i]/nangr[iz][isgz];
          meanADFr += map_ADFr[iz][isgz][i];
        }
      }
      meanADFr=meanADFr/std::accumulate(cnzero.begin(), iter_n, 0.0);
      for (size_t i=0; i<mskin.size(); ++i) {
        if (cnzero[i]==true) 
          map_ADFr[iz][isgz][i] = map_ADFr[iz][isgz][i] - meanADFr;
      }
    }
  }
  return map_ADFr;
}

std::tuple<std::vector<std::vector<std::vector<double>>> ,std::vector<std::vector<std::vector<double>>>> cbl::measure::MapCalculator::ADFr_ARFr (std::vector<double> const& mskin, std::vector<double> const& zH) const
{
  std::vector<bool> cnzero(mskin.size()); // For the mask application
    
  std::transform(mskin.begin(), mskin.end(), cnzero.begin(), [](double x) 
  { 
    return x>0.0; 
  });
  // Now cnzero is a boolean array (0 and 1), values depending on the matching of the condition in the lambda function 

  const unsigned int nzobs = m_zobs.size();
  const unsigned int nsgz = m_sgzar.size();
  const unsigned int size = m_data->nObjects();

  auto ipx1 = HealpixMapping();

  std::vector<std::vector<std::vector<double>>> map_ADFr (nzobs, std::vector<std::vector<double>>(nsgz, std::vector<double>(12 * m_ns * m_ns, 0.0)));
  std::vector<std::vector<std::vector<double>>> map_ARFr (nzobs, std::vector<std::vector<double>>(nsgz, std::vector<double>(12 * m_ns * m_ns, 0.0)));
  std::vector<std::vector<double>> nangr (nzobs, std::vector<double>(nsgz, 0.0));

  for (size_t iz=0; iz<nzobs; ++iz) {
    for (size_t isgz=0; isgz<nsgz; ++isgz) {
      std::vector<int> ipxA;
      std::vector<double> W1b;
      std::vector<int> true_indb;
      std::vector<bool> csel1b(size, 0);
      std::vector<double> zeff;
      std::vector<double> zav1b_vec;
      std::vector<bool> cnzero2(12*m_ns*m_ns, false);
      for (size_t i=0; i<ipx1.size(); ++i) {                
        csel1b[i] = (std::abs(zH[i]-m_zobs[iz])<m_nsg*m_sgzar[isgz]) && (mskin[ipx1[i]]!=0);
        if (csel1b[i]==true) 
          true_indb.push_back(i);
      }
      std::vector<double> red;
      for (size_t ind=0; ind<size; ++ind)
        red.push_back(m_data->redshift(ind));
      auto zH_sel = cbl::select<double>(red, true_indb);
      for (size_t i=0; i<true_indb.size(); ++i) 
        W1b.push_back(std::exp(-0.5*(std::pow((zH_sel[i]-m_zobs[iz]), 2)/std::pow(m_sgzar[isgz], 2))));
      auto ipx_sel = select<int>(ipx1, true_indb);
      for (size_t i=0; i<true_indb.size(); ++i) {  
        zav1b_vec.push_back(W1b[i]*zH_sel[i]);
        ipxA.push_back(ipx_sel[i]);
        zeff.push_back(zH_sel[i]);
      }

      double zav1b = std::accumulate(zav1b_vec.begin(), zav1b_vec.end(), 0.0)/std::accumulate(W1b.begin(), W1b.end(), 0.0);
      for (size_t ig=0; ig<true_indb.size(); ++ig) {
        map_ADFr[iz][isgz][ipxA[ig]] += W1b[ig];
        map_ARFr[iz][isgz][ipxA[ig]] += W1b[ig]*(zeff[ig]-zav1b); 
      }
      for (size_t i=0; i<mskin.size(); ++i) 
        cnzero2[i] = (map_ADFr[iz][isgz][i]!=0); // assigning bools
      for (size_t i=0; i<mskin.size(); ++i) { 
        if (cnzero[i]==true) {
          map_ADFr[iz][isgz][i] = map_ADFr[iz][isgz][i]/mskin[i];
          nangr[iz][isgz] += map_ADFr[iz][isgz][i];
        }
        if (cnzero2[i]==true) 
          map_ARFr[iz][isgz][i] = map_ARFr[iz][isgz][i]/mskin[i];
      }

      auto iter_n = std::next(cnzero.begin(), 12*m_ns*m_ns); //any difference?
      nangr[iz][isgz] = nangr[iz][isgz]/std::accumulate(cnzero.begin(), iter_n, 0.0);
      double meanADFr = 0;
      double meanARFr = 0;
      for (size_t i=0; i<mskin.size(); ++i) {
        if (cnzero[i]==true) { 
          map_ADFr[iz][isgz][i] = map_ADFr[iz][isgz][i]/nangr[iz][isgz];
          meanADFr += map_ADFr[iz][isgz][i];
        }
        if (cnzero2[i]==true) { 
          map_ARFr[iz][isgz][i] = map_ARFr[iz][isgz][i]/nangr[iz][isgz]; // non nang2??
          meanARFr += map_ARFr[iz][isgz][i];
        }
      }
      meanADFr = meanADFr/std::accumulate(cnzero.begin(), iter_n, 0.0);
      auto iter_n2 = std::next(cnzero2.begin(), 12*m_ns*m_ns); //any difference?
      meanARFr = meanARFr/std::accumulate(cnzero2.begin(), iter_n2, 0.0);
      for (size_t i=0; i<mskin.size(); ++i) {
        if (cnzero[i]==true) 
          map_ADFr[iz][isgz][i] = map_ADFr[iz][isgz][i]-meanADFr;
        if (cnzero2[i]==true) 
          map_ARFr[iz][isgz][i] = map_ARFr[iz][isgz][i]-meanARFr;
      }
    }
  }
  return std::make_tuple(map_ADFr, map_ARFr);
}

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> cbl::measure::MapCalculator::SN (std::vector<double> const& mskin, const int nran) const
{   
  std::vector<bool> cnzero(mskin.size()); // For the mask application

  std::transform(mskin.begin(), mskin.end(), cnzero.begin(), [](double x) 
  { 
    return x>0.0; 
  });
  // Now cnzero is a boolean array (0 and 1), values depending on the matching of the condition in the lambda function 

  const unsigned int nzobs = m_zobs.size();
  const unsigned int nsgz = m_sgzar.size();
  const unsigned int size = m_data->nObjects();
  cbl::random::UniformRandomNumbers rnd(0., 1., m_seed);

  auto ipx1 = HealpixMapping();

  auto map_SN = std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> (
											 nzobs, std::vector<std::vector<std::vector<std::vector<double>>>>(
																			   nsgz, std::vector<std::vector<std::vector<double>>>(
																									       nran, std::vector<std::vector<double>>(
																														      2, std::vector<double>(12 * m_ns * m_ns, 0.0)
																														      )
																									       )
																			   )
											 );

  for (size_t iz=0; iz<nzobs; ++iz) {
    for (size_t isgz=0; isgz<nsgz; ++isgz) {
      std::vector<int> ipxA;
      std::vector<double> W1a;
      std::vector<int> true_inda;
      std::vector<bool> csel1a(size, 0);
      std::vector<double> zeff;
      std::vector<double> zav1a_vec;
      std::vector<bool> cnzero2(12*m_ns*m_ns, false);
      for (size_t i=0; i<ipx1.size(); ++i) {                
        csel1a[i] = ((std::abs(m_data->redshift(i)-m_zobs[iz])<(m_nsg*m_sgzar[isgz])) && (mskin[ipx1[i]]!=0));
        if (csel1a[i]==true) 
          true_inda.push_back(i);
      }
      // m_data->var(cbl::catalogue::Var::_Redshift_) causes memory leak ....
      std::vector<double> red;
      for (size_t ind=0; ind<size; ++ind)
        red.push_back(m_data->redshift(ind));
      auto zin_sel = cbl::select<double>(red, true_inda);
      for (size_t i=0; i<true_inda.size(); ++i)   
        W1a.push_back(std::exp(-0.5*(std::pow((zin_sel[i]-m_zobs[iz]), 2)/std::pow(m_sgzar[isgz], 2))));
      auto ipx_sel = select<int>(ipx1, true_inda);
      for (size_t i=0; i<true_inda.size(); ++i) {  
        ipxA.push_back(ipx_sel[i]);
        zav1a_vec.push_back((W1a[i]*zin_sel[i]));
        zeff.push_back(zin_sel[i]);
      }

      double zav1a = std::accumulate(zav1a_vec.begin(), zav1a_vec.end(), 0.0)/std::accumulate(W1a.begin(), W1a.end(), 0.0);
      double cumsum = 0;
      std::vector<int> true_cnzero;
      for (size_t i=0; i<cnzero.size(); ++i) {
        if (cnzero[i]==true)
          true_cnzero.push_back(i);
      }
      std::vector<double> cum = cbl::select<double>(mskin, true_cnzero);
      for (size_t i=0; i<cum.size(); ++i) {
        if (i>0) 
          cumsum = cum[i-1];
        cum[i]= cum[i]+cumsum;
      }
      std::vector<double> Cprobar; 
      for (size_t i=0; i<cum.size(); ++i) 
        Cprobar.push_back(cum[i]/std::accumulate(mskin.begin(), mskin.end(), 0.0));
      std::vector<double> jeff = cbl::arange(0, static_cast<unsigned int>(std::accumulate(cnzero.begin(), cnzero.end(), 0)), 1);
      assert(Cprobar.size()==jeff.size());
      std::vector<int> jindex;
      auto all_jindex = cbl::arange(12*m_ns*m_ns);
      for (size_t i=0; i<cnzero.size(); ++i) {
        if (cnzero[i]==true)
          jindex.push_back(all_jindex[i]);
      }
      std::vector<std::vector<double>> probg (nran, std::vector<double>(true_inda.size(), 0.0));
            
      // Filling with randoms
      for (int r=0; r<nran; ++r) {
        for (size_t c=0; c<true_inda.size(); ++c) 
          probg[r][c] = rnd();
      }
      std::vector<std::vector<int>> jgar0 (nran, std::vector<int>(true_inda.size(), 0.0));
      for (int r=0; r<nran; ++r) {
        for (size_t c=0; c<true_inda.size(); ++c) {
          //sets vanishing outrange values
          if (probg[r][c]>*std::max_element(Cprobar.begin(), Cprobar.end()) || probg[r][c]<*std::min_element(Cprobar.begin(), Cprobar.end())) 
            jgar0[r][c] = 0;
          else jgar0[r][c] = static_cast<int>(cbl::interpolated(probg[r][c], Cprobar, jeff, "Linear"));
        }
      }
      std::vector<std::vector<double>> jgar (nran, std::vector<double>(true_inda.size(), 0.0));
      for (int r=0; r<nran; ++r) {
        for (size_t c=0; c<true_inda.size(); ++c) 
          jgar[r][c] = jindex[jgar0[r][c]];
      }
      for (int jran=0; jran<nran; ++jran) {
        for (size_t ig=0; ig<true_inda.size(); ++ig) {
          int jg = jgar[jran][ig];
          map_SN[iz][isgz][jran][0][jg] += W1a[ig];
          map_SN[iz][isgz][jran][1][jg] += W1a[ig]*(zeff[ig]-zav1a);
        }
        for (size_t i=0; i<mskin.size(); ++i) 
          cnzero2[i] = map_SN[iz][isgz][jran][0][i]!=0; // assigning bools
        std::vector<double> nangT_vec;
        std::vector<double> nangT2_vec;
        for (size_t i=0; i<mskin.size(); ++i) { 
          if (cnzero[i]==true) {
            map_SN[iz][isgz][jran][0][i] = map_SN[iz][isgz][jran][0][i]/mskin[i];
            nangT_vec.push_back(map_SN[iz][isgz][jran][0][i]);
          }
        }
        for (size_t i=0; i<mskin.size(); ++i) { 
          if (cnzero2[i]==true) 
            map_SN[iz][isgz][jran][1][i] = map_SN[iz][isgz][jran][1][i]/mskin[i];
        }
        double nangT = std::accumulate(nangT_vec.begin(), nangT_vec.end(), 0.0)/nangT_vec.size();
        double meanSN0 = 0;
        double meanSN1 = 0;
        for (size_t i=0; i<mskin.size(); ++i) { 
          if (cnzero[i]==true) {
            map_SN[iz][isgz][jran][0][i] = map_SN[iz][isgz][jran][0][i]/nangT;
            meanSN0 += map_SN[iz][isgz][jran][0][i];
          }
          if (cnzero2[i]==true) {
            map_SN[iz][isgz][jran][1][i] = map_SN[iz][isgz][jran][1][i]/nangT;
            meanSN1 += map_SN[iz][isgz][jran][1][i];
          }
        }
        meanSN0 = meanSN0/std::accumulate(cnzero.begin(), cnzero.end(), 0.0);
        meanSN1 = meanSN1/std::accumulate(cnzero2.begin(), cnzero2.end(), 0.0);
        for (size_t i=0; i<mskin.size(); ++i) {
          if (cnzero2[i]==true) 
            map_SN[iz][isgz][jran][1][i] = map_SN[iz][isgz][jran][1][i]-meanSN1;
          if (cnzero[i]==true) 
            map_SN[iz][isgz][jran][0][i] = map_SN[iz][isgz][jran][0][i]-meanSN0;
        }
      }
    }
  }
  return map_SN;
}

void cbl::measure::MapCalculator::print_map (const std::vector<std::vector<std::vector<double>>>& map, const std::string& filename) const {

  std::ofstream outputFile(filename);

  if (!outputFile.is_open()) {
    ErrorCBL("An error occured in opening the file "+filename, "print_map", "MapCalculator.cpp");
    return;
  }
  // Printing map_ADF as 3D array
  for (size_t ra=0; ra<map.size(); ++ra) {
    for (size_t dec=0; dec<map[ra].size(); ++dec) {
      for (size_t z=0; z<map[ra][dec].size(); ++z) {
        outputFile << map[ra][dec][z] << " ";
      }
      outputFile << std::endl;
    }
    outputFile << std::endl;
  }
  // Closing file
  outputFile.close();
  coutCBL << "Data printed on " << filename << std::endl;
}
