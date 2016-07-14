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
 *  @file CatalogueAnalysis/TwoPointCorrelation/Pair.cpp
 *
 *  @brief Methods of the class Pair  
 *
 *  This file contains the implementation of all the methods of the
 *  classes Pair*, used to handle pairs of objects of any kind
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Pair1D.h"
#include "Pair2D.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;


// ============================================================================================


shared_ptr<Pair> cosmobl::pairs::Pair::Create (const PairType type, const double Min, const double Max, const int nbins, const double shift, const CoordUnits angularUnits, function<double(double)> angularWeight)
{
  if (type==_angular_lin_) return move(unique_ptr<Pair1D_angular_lin>{new Pair1D_angular_lin(Min, Max, nbins, shift, angularUnits, angularWeight)});
  else if (type==_angular_log_) return move(unique_ptr<Pair1D_angular_log>{new Pair1D_angular_log(Min, Max, nbins, shift, angularUnits, angularWeight)});
  else if (type==_comoving_lin_) return move(unique_ptr<Pair1D_comoving_lin>{new Pair1D_comoving_lin(Min, Max, nbins, shift, angularUnits, angularWeight)});
  else if (type==_comoving_log_) return move(unique_ptr<Pair1D_comoving_log>{new Pair1D_comoving_log(Min, Max, nbins, shift, angularUnits, angularWeight)});

  else ErrorMsg("Error in cosmobl::pairs::Create of Pairs.cpp: no such type of object!");
  
  return NULL;
}


// ============================================================================================


shared_ptr<Pair> cosmobl::pairs::Pair::Create (const PairType type, const double Min, const double Max, const double binSize, const double shift, const CoordUnits angularUnits, function<double(double)> angularWeight)
{
  if (type==_angular_lin_) return move(unique_ptr<Pair1D_angular_lin>{new Pair1D_angular_lin(Min, Max, binSize, shift, angularUnits, angularWeight)});
  else if (type==_angular_log_) return move(unique_ptr<Pair1D_angular_log>{new Pair1D_angular_log(Min, Max, binSize, shift, angularUnits, angularWeight)});
  else if (type==_comoving_lin_) return move(unique_ptr<Pair1D_comoving_lin>{new Pair1D_comoving_lin(Min, Max, binSize, shift, angularUnits, angularWeight)});
  else if (type==_comoving_log_) return move(unique_ptr<Pair1D_comoving_log>{new Pair1D_comoving_log(Min, Max, binSize, shift, angularUnits, angularWeight)});

  else ErrorMsg("Error in cosmobl::pairs::Create of Pairs.cpp: no such type of object!");
  
  return NULL;
}


// ============================================================================================


shared_ptr<Pair> cosmobl::pairs::Pair::Create (const PairType type, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const double Min_D2, const double Max_D2, const int nbins_D2, const double shift_D2, const CoordUnits angularUnits, function<double(double)> angularWeight)
{
  if (type==_comovingCartesian_linlin_) return move(unique_ptr<Pair2D_comovingCartesian_linlin>{new Pair2D_comovingCartesian_linlin(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});
  
  else if (type==_comovingCartesian_linlog_) return move(unique_ptr<Pair2D_comovingCartesian_linlog>{new Pair2D_comovingCartesian_linlog(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});

  else if (type==_comovingCartesian_loglin_) return move(unique_ptr<Pair2D_comovingCartesian_loglin>{new Pair2D_comovingCartesian_loglin(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});

  else if (type==_comovingCartesian_loglog_) return move(unique_ptr<Pair2D_comovingCartesian_loglog>{new Pair2D_comovingCartesian_loglog(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});

  else if (type==_comovingPolar_linlin_) return move(unique_ptr<Pair2D_comovingPolar_linlin>{new Pair2D_comovingPolar_linlin(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});
  
  else if (type==_comovingPolar_linlog_) return move(unique_ptr<Pair2D_comovingPolar_linlog>{new Pair2D_comovingPolar_linlog(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});

  else if (type==_comovingPolar_loglin_) return move(unique_ptr<Pair2D_comovingPolar_loglin>{new Pair2D_comovingPolar_loglin(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});

  else if (type==_comovingPolar_loglog_) return move(unique_ptr<Pair2D_comovingPolar_loglog>{new Pair2D_comovingPolar_loglog(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});
  
  else ErrorMsg("Error in cosmobl::pairs::Create of Pairs.cpp: no such type of object!");
 
  return NULL;
}


// ============================================================================================


shared_ptr<Pair> cosmobl::pairs::Pair::Create (const PairType type, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const double Min_D2, const double Max_D2, const double binSize_D2, const double shift_D2, const CoordUnits angularUnits, function<double(double)> angularWeight)
{
  if (type==_comovingCartesian_linlin_) return move(unique_ptr<Pair2D_comovingCartesian_linlin>{new Pair2D_comovingCartesian_linlin(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});
  
  else if (type==_comovingCartesian_linlog_) return move(unique_ptr<Pair2D_comovingCartesian_linlog>{new Pair2D_comovingCartesian_linlog(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

  else if (type==_comovingCartesian_loglin_) return move(unique_ptr<Pair2D_comovingCartesian_loglin>{new Pair2D_comovingCartesian_loglin(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

  else if (type==_comovingCartesian_loglog_) return move(unique_ptr<Pair2D_comovingCartesian_loglog>{new Pair2D_comovingCartesian_loglog(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

  else if (type==_comovingPolar_linlin_) return move(unique_ptr<Pair2D_comovingPolar_linlin>{new Pair2D_comovingPolar_linlin(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});
  
  else if (type==_comovingPolar_linlog_) return move(unique_ptr<Pair2D_comovingPolar_linlog>{new Pair2D_comovingPolar_linlog(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

  else if (type==_comovingPolar_loglin_) return move(unique_ptr<Pair2D_comovingPolar_loglin>{new Pair2D_comovingPolar_loglin(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

  else if (type==_comovingPolar_loglog_) return move(unique_ptr<Pair2D_comovingPolar_loglog>{new Pair2D_comovingPolar_loglog(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});
  
  else ErrorMsg("Error in cosmobl::pairs::Create of Pairs.cpp: no such type of object!");
 
  return NULL;
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair1D_angular_lin::m_set_parameters_nbins ()
{
  double binSize = ((m_thetaMax-m_thetaMin)/m_nbins);
  m_binSize_inv = 1./binSize;
  
  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = binSize*(i+m_shift)+m_thetaMin;
}


// ============================================================================================


void cosmobl::pairs::Pair1D_angular_lin::m_set_parameters_binSize ()
{
  m_nbins = nint((m_thetaMax-m_thetaMin)*m_binSize_inv);
  m_thetaMax = m_nbins/m_binSize_inv+m_thetaMin;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = (i+m_shift)/m_binSize_inv+m_thetaMin;
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair1D_angular_log::m_set_parameters_nbins ()
{
  if (m_thetaMin<1.e-30) ErrorMsg("Error in cosmobl::pairs::Pair1D_angular_log::m_set_parameters_nbins of Pair.cpp: Min must be >0!");
  
  double binSize = ((log10(m_thetaMax)-log10(m_thetaMin))/m_nbins);
  m_binSize_inv = 1./binSize;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = pow(10.,(i+m_shift)*binSize+log10(m_thetaMin));
}


// ============================================================================================


void cosmobl::pairs::Pair1D_angular_log::m_set_parameters_binSize ()
{
  if (m_thetaMin<1.e-30) ErrorMsg("Error in cosmobl::pairs::Pair1D_angular_log::m_set_parameters_binSize of Pair.cpp: Min must be >0!");

  m_nbins = nint((log10(m_thetaMax)-log10(m_thetaMin))*m_binSize_inv);
  m_thetaMax = pow(10.,m_nbins/m_binSize_inv+log10(m_thetaMin));

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = pow(10.,(i+m_shift)/m_binSize_inv+log10(m_thetaMin));
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair1D_comoving_lin::m_set_parameters_nbins ()
{
  double binSize = ((m_rMax-m_rMin)/m_nbins);
  m_binSize_inv = 1./binSize;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = (i+m_shift)*binSize+m_rMin;
}


// ============================================================================================


void cosmobl::pairs::Pair1D_comoving_lin::m_set_parameters_binSize ()
{
  m_nbins = nint((m_rMax-m_rMin)*m_binSize_inv);
  m_rMax = m_nbins/m_binSize_inv+m_rMin;
  
  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = (i+m_shift)/m_binSize_inv+m_rMin;
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair1D_comoving_log::m_set_parameters_nbins ()
{
  if (m_rMin<1.e-30) ErrorMsg("Error in cosmobl::pairs::Pair1D_comoving_log::m_set_parameters_nbins of Pair.cpp: Min must be >0!");
  
  double binSize = ((log10(m_rMax)-log10(m_rMin))/m_nbins);
  m_binSize_inv = 1./binSize;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = pow(10.,(i+m_shift)*binSize+log10(m_rMin));
}


// ============================================================================================


void cosmobl::pairs::Pair1D_comoving_log::m_set_parameters_binSize ()
{
  if (m_rMin<1.e-30) ErrorMsg("Error in cosmobl::pairs::Pair1D_comoving_log::m_set_parameters_binSize of Pair.cpp: Min must be >0!");
  
  m_nbins = nint((log10(m_rMax)-log10(m_rMin))*m_binSize_inv);
  m_rMax = pow(10.,m_nbins/m_binSize_inv+log10(m_rMin));

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = pow(10.,(i+m_shift)/m_binSize_inv+log10(m_rMin));
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_linlin::m_set_parameters_nbins ()
{
  double binSize_D1 = ((m_rpMax-m_rpMin)/m_nbins_D1);
  m_binSize_inv_D1 = 1./binSize_D1;

  double binSize_D2 = ((m_piMax-m_piMin)/m_nbins_D2);
  m_binSize_inv_D2 = 1./binSize_D2;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = (i+m_shift_D1)*binSize_D1+m_rpMin;
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = (i+m_shift_D2)*binSize_D2+m_piMin;
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_linlin::m_set_parameters_binSize ()
{
  m_nbins_D1 = nint((m_rpMax-m_rpMin)*m_binSize_inv_D1);
  m_rpMax = m_nbins_D1/m_binSize_inv_D1+m_rpMin;

  m_nbins_D2 = nint((m_piMax-m_piMin)*m_binSize_inv_D2);
  m_piMax = m_nbins_D2/m_binSize_inv_D2+m_piMin;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = (i+m_shift_D1)/m_binSize_inv_D1+m_rpMin;
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = (i+m_shift_D2)/m_binSize_inv_D2+m_piMin;
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_linlog::m_set_parameters_nbins ()
{
  if (m_piMin<1.e-30)
    ErrorMsg("Error in cosmobl::pairs::Pair2D_comovingCartesian_linlog::m_set_parameters_nbins of Pair.cpp: m_piMin must be >0!");
  
  double binSize_D1 = ((m_rpMax-m_rpMin)/m_nbins_D1);
  m_binSize_inv_D1 = 1./binSize_D1;

  double binSize_D2 = ((log10(m_piMax)-log10(m_piMin))/m_nbins_D2);
  m_binSize_inv_D2 = 1./binSize_D2;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = (i+m_shift_D1)*binSize_D1+m_rpMin;
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = pow(10.,(i+m_shift_D2)*binSize_D2+log10(m_piMin));
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_linlog::m_set_parameters_binSize ()
{
  if (m_piMin<1.e-30)
    ErrorMsg("Error in cosmobl::pairs::Pair2D_comovingCartesian_linlog::m_set_parameters_binSize of Pair.cpp: m_piMin must be >0!");
  
  m_nbins_D1 = nint((m_rpMax-m_rpMin)*m_binSize_inv_D1);
  m_rpMax = m_nbins_D1/m_binSize_inv_D1+m_rpMin;

  m_nbins_D2 = nint((log10(m_piMax)-log10(m_piMin))*m_binSize_inv_D2);
  m_piMax = pow(10.,m_nbins_D2/m_binSize_inv_D2+log10(m_piMin));

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = (i+m_shift_D1)/m_binSize_inv_D1+m_rpMin;
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = pow(10.,(i+m_shift_D2)/m_binSize_inv_D2+log10(m_piMin));
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_loglin::m_set_parameters_nbins ()
{
  if (m_rpMin<1.e-30)
    ErrorMsg("Error in cosmobl::pairs::Pair2D_comovingCartesian_loglin::m_set_parameters_nbins of Pair.cpp: m_rpMin must be >0!");
  
  double binSize_D1 = ((log10(m_rpMax)-log10(m_rpMin))/m_nbins_D1);
  m_binSize_inv_D1 = 1./binSize_D1;
  
  double binSize_D2 = ((m_piMax-m_piMin)/m_nbins_D2);
  m_binSize_inv_D2 = 1./binSize_D2;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = pow(10.,(i+m_shift_D1)*binSize_D1+log10(m_rpMin)); 
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = (i+m_shift_D2)*binSize_D2+m_piMin;
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_loglin::m_set_parameters_binSize ()
{
  if (m_rpMin<1.e-30)
    ErrorMsg("Error in cosmobl::pairs::Pair2D_comovingCartesian_linlin::m_set_parameters_binSize of Pair.cpp: m_rpMin must be >0!");
  
  m_nbins_D1 = nint((log10(m_rpMax)-log10(m_rpMin))*m_binSize_inv_D1);
  m_rpMax = pow(10.,m_nbins_D1/m_binSize_inv_D1+log10(m_rpMin));
  
  m_nbins_D2 = nint((m_piMax-m_piMin)*m_binSize_inv_D2);
  m_piMax = m_nbins_D2/m_binSize_inv_D2+m_piMin;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = pow(10.,(i+m_shift_D1)/m_binSize_inv_D1+log10(m_rpMin));
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = (i+m_shift_D2)/m_binSize_inv_D2+m_piMin;
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_loglog::m_set_parameters_nbins ()
{
  if (m_rpMin<1.e-30 || m_piMin<1.e-30)
    ErrorMsg("Error in cosmobl::pairs::Pair2D_comovingCartesian_loglog::m_set_parameters_nbins of Pair.cpp: m_rpMin and m_piMin must be >0!");
  
  double binSize_D1 = ((log10(m_rpMax)-log10(m_rpMin))/m_nbins_D1);
  m_binSize_inv_D1 = 1./binSize_D1;
  
  double binSize_D2 = ((log10(m_piMax)-log10(m_piMin))/m_nbins_D2);
  m_binSize_inv_D2 = 1./binSize_D2;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = pow(10.,(i+m_shift_D1)*binSize_D1+log10(m_rpMin)); 
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = pow(10.,(i+m_shift_D2)*binSize_D2+log10(m_piMin)); 
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_loglog::m_set_parameters_binSize ()
{
  if (m_rpMin<1.e-30 || m_piMin<1.e-30)
    ErrorMsg("Error in cosmobl::pairs::Pair2D_comovingCartesian_loglog::m_set_parameters_binSize of Pair.cpp: m_rpMin and m_piMin must be >0!");
  
  m_nbins_D1 = nint((log10(m_rpMax)-log10(m_rpMin))*m_binSize_inv_D1);
  m_rpMax = pow(10.,m_nbins_D1/m_binSize_inv_D1+log10(m_rpMin));
  
  m_nbins_D2 = nint((log10(m_piMax)-log10(m_piMin))*m_binSize_inv_D2);
  m_piMax = pow(10.,m_nbins_D2/m_binSize_inv_D2+log10(m_piMin));

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = pow(10.,(i+m_shift_D1)/m_binSize_inv_D1+log10(m_rpMin));
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = pow(10.,(i+m_shift_D2)/m_binSize_inv_D2+log10(m_piMin));
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_linlin::m_set_parameters_nbins ()
{ 
  double binSize_D1 = ((m_rMax-m_rMin)/m_nbins_D1);
  m_binSize_inv_D1 = 1./binSize_D1;

  double binSize_D2 = ((m_muMax-m_muMin)/m_nbins_D2);
  m_binSize_inv_D2 = 1./binSize_D2;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = (i+m_shift_D1)*binSize_D1+m_rMin;
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = (i+m_shift_D2)*binSize_D2+m_muMin;
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_linlin::m_set_parameters_binSize ()
{
  m_nbins_D1 = nint((m_rMax-m_rMin)*m_binSize_inv_D1);
  m_rMax = m_nbins_D1/m_binSize_inv_D1+m_rMin;

  m_nbins_D2 = nint((m_muMax-m_muMin)*m_binSize_inv_D2);
  m_muMax = m_nbins_D2/m_binSize_inv_D2+m_muMin;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = (i+m_shift_D1)/m_binSize_inv_D1+m_rMin;
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = (i+m_shift_D2)/m_binSize_inv_D2+m_muMin;
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_linlog::m_set_parameters_nbins ()
{
  if (m_muMin<1.e-30)
    ErrorMsg("Error in cosmobl::pairs::Pair2D_comovingPolar_linlog::m_set_parameters_nbins of Pair.cpp: m_muMin must be >0!");
  
  double binSize_D1 = ((m_rMax-m_rMin)/m_nbins_D1);
  m_binSize_inv_D1 = 1./binSize_D1;

  double binSize_D2 = ((log10(m_muMax)-log10(m_muMin))/m_nbins_D2);
  m_binSize_inv_D2 = 1./binSize_D2;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = (i+m_shift_D1)*binSize_D1+m_rMin;
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = pow(10.,(i+m_shift_D2)*binSize_D2+log10(m_muMin));
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_linlog::m_set_parameters_binSize ()
{
  if (m_muMin<1.e-30)
    ErrorMsg("Error in cosmobl::pairs::Pair2D_comovingPolar_linlog::m_set_parameters_linlog of Pair.cpp: m_muMin must be >0!");
  
  m_nbins_D1 = nint((m_rMax-m_rMin)*m_binSize_inv_D1);
  m_rMax = m_nbins_D1/m_binSize_inv_D1+m_rMin;

  m_nbins_D2 = nint((log10(m_muMax)-log10(m_muMin))*m_binSize_inv_D2);
  m_muMax = pow(10.,(m_nbins_D2-m_shift_D2)/m_binSize_inv_D2+log10(m_muMin));

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = (i+m_shift_D1)/m_binSize_inv_D1+m_rMin;
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = pow(10.,(i+m_shift_D2)/m_binSize_inv_D2+log10(m_muMin));
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_loglin::m_set_parameters_nbins ()
{
  if (m_rMin<1.e-30)
    ErrorMsg("Error in cosmobl::pairs::Pair2D_comovingPolar_loglin::m_set_parameters_nbins of Pair.cpp: m_rMin must be >0!");
  
  double binSize_D1 = ((log10(m_rMax)-log10(m_rMin))/m_nbins_D1);
  m_binSize_inv_D1 = 1./binSize_D1;
  
  double binSize_D2 = ((m_muMax-m_muMin)/m_nbins_D2);
  m_binSize_inv_D2 = 1./binSize_D2;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = pow(10.,(i+m_shift_D1)*binSize_D1+log10(m_rMin)); 
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = (i+m_shift_D2)*binSize_D2+m_muMin;
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_loglin::m_set_parameters_binSize ()
{
  if (m_rMin<1.e-30)
    ErrorMsg("Error in cosmobl::pairs::Pair2D_comovingPolar_loglin::m_set_parameters_binSize of Pair.cpp: m_rMin must be >0!");
  
  m_nbins_D1 = nint((log10(m_rMax)-log10(m_rMin))*m_binSize_inv_D1);
  m_rMax = pow(10.,m_nbins_D1/m_binSize_inv_D1+log10(m_rMin));
  
  m_nbins_D2 = nint((m_muMax-m_muMin)*m_binSize_inv_D2);
  m_muMax = m_nbins_D2/m_binSize_inv_D2+m_muMin;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = pow(10.,(i+m_shift_D1)/m_binSize_inv_D1+log10(m_rMin));
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = (i+m_shift_D2)/m_binSize_inv_D2+m_muMin;
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_loglog::m_set_parameters_nbins ()
{
  if (m_rMin<1.e-30 || m_muMin<1.e-30)
    ErrorMsg("Error in cosmobl::pairs::Pair2D_comovingPolar_loglog::m_set_parameters_nbins of Pair.cpp: m_rMin and m_muMin must be >0!");
  
  double binSize_D1 = ((log10(m_rMax)-log10(m_rMin))/m_nbins_D1);
  m_binSize_inv_D1 = 1./binSize_D1;
  
  double binSize_D2 = ((log10(m_muMax)-log10(m_muMin))/m_nbins_D2);
  m_binSize_inv_D2 = 1./binSize_D2;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = pow(10.,(i+m_shift_D1)*binSize_D1+log10(m_rMin)); 
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = pow(10.,(i+m_shift_D2)*binSize_D2+log10(m_muMin)); 
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_loglog::m_set_parameters_binSize ()
{
  if (m_rMin<1.e-30 || m_muMin<1.e-30)
    ErrorMsg("Error in cosmobl::pairs::Pair2D_comovingPolar_loglog::m_set_parameters_binSize of Pair.cpp: m_rMin and m_muMin must be >0!");
  
  m_nbins_D1 = nint((log10(m_rMax)-log10(m_rMin))*m_binSize_inv_D1);
  m_rMax = pow(10.,m_nbins_D1/m_binSize_inv_D1+log10(m_rMin));
  
  m_nbins_D2 = nint((log10(m_muMax)-log10(m_muMin))*m_binSize_inv_D2);
  m_muMax = pow(10.,m_nbins_D2/m_binSize_inv_D2+log10(m_muMin));

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = pow(10.,(i+m_shift_D1)/m_binSize_inv_D1+log10(m_rMin));
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = pow(10.,(i+m_shift_D2)/m_binSize_inv_D2+log10(m_muMin));
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair1D_angular_lin::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  double tt = (m_angularUnits==_radians_) ? angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()) 
	     : converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), _radians_, m_angularUnits);
  
  if (m_thetaMin < tt && tt < m_thetaMax) {

    int kk = max(0, min(int((tt-m_thetaMin)*m_binSize_inv), m_nbins));

    m_PP1D[kk] += obj1->weight()*obj2->weight();
    
  } 
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair1D_angular_log::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  double tt = (m_angularUnits==_radians_) ? angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()) 
	     : converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), _radians_, m_angularUnits);
  
  if (m_thetaMin < tt && tt < m_thetaMax) {

    int kk = max(0, min(int((log10(tt)-log10(m_thetaMin))*m_binSize_inv), m_nbins));

    m_PP1D[kk] += obj1->weight()*obj2->weight();
    
  }
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair1D_comoving_lin::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  double tt = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 
  
  if (m_rMin < tt && tt < m_rMax) {

    int kk = max(0, min(int((tt-m_rMin)*m_binSize_inv), m_nbins));
    
    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), _radians_, m_angularUnits));
    
    m_PP1D[kk] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair1D_comoving_log::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  double tt = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 
  
  if (m_rMin < tt && tt < m_rMax) {

    int kk = max(0, min(int((log10(tt)-log10(m_rMin))*m_binSize_inv), m_nbins));

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), _radians_, m_angularUnits));
  
    m_PP1D[kk] += obj1->weight()*obj2->weight()*angWeight;

  }
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_linlin::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  double rp = perpendicular_distance(obj1->ra(), obj2->ra(), obj1->dec(), obj2->dec(), obj1->dc(), obj2->dc());
  double pi = fabs(obj1->dc()-obj2->dc());
  
  if (m_rpMin<rp && rp<m_rpMax && m_piMin<pi && pi<m_piMax) {

    int ir = max(0, min(int((rp-m_rpMin)*m_binSize_inv_D1), m_nbins_D1));
    int jr = max(0, min(int((pi-m_piMin)*m_binSize_inv_D2), m_nbins_D2));

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), _radians_, m_angularUnits));
    
    m_PP2D[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_linlog::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  double rp = perpendicular_distance(obj1->ra(), obj2->ra(), obj1->dec(), obj2->dec(), obj1->dc(), obj2->dc());
  double pi = fabs(obj1->dc()-obj2->dc());
  
  if (m_rpMin<rp && rp<m_rpMax && m_piMin<pi && pi<m_piMax) {
    
    int ir = max(0, min(int((rp-m_rpMin)*m_binSize_inv_D1), m_nbins_D1));
    int jr = max(0, min(int((log10(pi)-log10(m_piMin))*m_binSize_inv_D2), m_nbins_D2));

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), _radians_, m_angularUnits));
    
    m_PP2D[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_loglin::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  double rp = perpendicular_distance(obj1->ra(), obj2->ra(), obj1->dec(), obj2->dec(), obj1->dc(), obj2->dc());
  double pi = fabs(obj1->dc()-obj2->dc());

  if (m_rpMin<rp && rp<m_rpMax && m_piMin<pi && pi<m_piMax) {

    int ir = max(0, min(int((log10(rp)-log10(m_rpMin))*m_binSize_inv_D1), m_nbins_D1)); 
    int jr = max(0, min(int((pi-m_piMin)*m_binSize_inv_D2), m_nbins_D2));

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), _radians_, m_angularUnits));
 
    m_PP2D[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_loglog::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  double rp = perpendicular_distance(obj1->ra(), obj2->ra(), obj1->dec(), obj2->dec(), obj1->dc(), obj2->dc());
  double pi = fabs(obj1->dc()-obj2->dc());
  
  if (m_rpMin<rp && rp<m_rpMax && m_piMin<pi && pi<m_piMax) {

    int ir = max(0, min(int((log10(rp)-log10(m_rpMin))*m_binSize_inv_D1), m_nbins_D1)); 
    int jr = max(0, min(int((log10(pi)-log10(m_piMin))*m_binSize_inv_D2), m_nbins_D2)); 

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), _radians_, m_angularUnits));
 
    m_PP2D[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_linlin::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  double rr = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz());
  double mu = fabs(obj1->dc()-obj2->dc())/rr;
    
  if (m_rMin<rr && rr<m_rMax && m_muMin<mu && mu<m_muMax) {
    
    int ir = max(0, min(int((rr-m_rMin)*m_binSize_inv_D1), m_nbins_D1));
    int jr = max(0, min(int((mu-m_muMin)*m_binSize_inv_D2), m_nbins_D2));

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), _radians_, m_angularUnits));
 
    m_PP2D[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_linlog::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  double rr = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz());
  double mu = fabs(obj1->dc()-obj2->dc())/rr;
  
  if (m_rMin<rr && rr<m_rMax && m_muMin<mu && mu<m_muMax) {
    
    int ir = max(0, min(int((rr-m_rMin)*m_binSize_inv_D1), m_nbins_D1));
    int jr = max(0, min(int((log10(mu)-log10(m_muMin))*m_binSize_inv_D2), m_nbins_D2)); 

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), _radians_, m_angularUnits));
    
    m_PP2D[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
  }
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_loglin::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  double rr = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz());
  double mu = fabs(obj1->dc()-obj2->dc())/rr;
  
  if (m_rMin<rr && rr<m_rMax && m_muMin<mu && mu<m_muMax) {

    int ir = max(0, min(int((log10(rr)-log10(m_rMin))*m_binSize_inv_D1), m_nbins_D1)); 
    int jr = max(0, min(int((mu-m_muMin)*m_binSize_inv_D2), m_nbins_D2));

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), _radians_, m_angularUnits));
    
    m_PP2D[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_loglog::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  double rr = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz());
  double mu = fabs(obj1->dc()-obj2->dc())/rr;
  
  if (m_rMin<rr && rr<m_rMax && m_muMin<mu && mu<m_muMax) {

    int ir = max(0, min(int((log10(rr)-log10(m_rMin))*m_binSize_inv_D1), m_nbins_D1)); 
    int jr = max(0, min(int((log10(mu)-log10(m_muMin))*m_binSize_inv_D2), m_nbins_D2)); 

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), _radians_, m_angularUnits));
    
    m_PP2D[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================
// ============================================================================================


void cosmobl::pairs::Pair1D::Sum (const shared_ptr<Pair> pp, const double ww)
{
  if (m_nbins != pp->nbins()) 
    ErrorMsg("Error in cosmobl::pairs::Pair1D::Sum of Pair.cpp: dimension problems!");
  
  for (int i=0; i<m_nbins; i++) 
    m_PP1D[i] += ww*pp->PP1D(i);
}


// ============================================================================
// ============================================================================================


void cosmobl::pairs::Pair2D::Sum (const shared_ptr<Pair> pp, const double ww)
{
  if (m_nbins_D1 != pp->nbins_D1() || m_nbins_D2 != pp->nbins_D2()) 
    ErrorMsg("Error in cosmobl::pairs::Pair2D::Sum of Pair.cpp: dimension problems!");
    
  for (int i=0; i<m_nbins_D1; i++)
    for (int j=0; j<m_nbins_D2; j++)
      m_PP2D[i][j] += ww*pp->PP2D(i,j);
}

