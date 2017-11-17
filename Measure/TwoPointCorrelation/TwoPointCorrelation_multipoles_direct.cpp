/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file
 *  CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation_multipoles_direct.cpp
 *
 *  @brief Methods of the class
 *  TwoPointCorrelation_multipoles_direct used to measure the
 *  multipoles of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation_multipoles_direct used to measure the
 *  multipoles of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation_multipoles_direct.h"

using namespace cosmobl;
using namespace catalogue;
using namespace chainmesh;
using namespace data;
using namespace pairs;
using namespace measure;
using namespace twopt;


// ============================================================================================


void cosmobl::measure::twopt::TwoPointCorrelation_multipoles_direct::set_parameters (const binType binType, const double rMin, const double rMax, const int nbins, const double shift, const CoordUnits angularUnits, function<double(double)> angularWeight, const bool compute_extra_info) 
{
  if (!compute_extra_info) 
    m_dd = (binType==_logarithmic_) ? move(Pair::Create(_comoving_multipoles_log_, _standard_, rMin, rMax, nbins, shift, angularUnits, angularWeight))
      : move(Pair::Create(_comoving_multipoles_lin_, _standard_, rMin, rMax, nbins, shift, angularUnits, angularWeight));
  else 
    m_dd = (binType==_logarithmic_) ? move(Pair::Create(_comoving_multipoles_log_, _extra_, rMin, rMax, nbins, shift, angularUnits, angularWeight))
      : move(Pair::Create(_comoving_multipoles_lin_, _extra_, rMin, rMax, nbins, shift, angularUnits, angularWeight));
    
  m_rr = (binType==_logarithmic_) ? move(Pair::Create(_comoving_multipoles_log_, _standard_, rMin, rMax, nbins, shift, angularUnits))
    : move(Pair::Create(_comoving_multipoles_lin_, _standard_, rMin, rMax, nbins, shift, angularUnits));
  
  m_dr = (binType==_logarithmic_) ? move(Pair::Create(_comoving_multipoles_log_, _standard_, rMin, rMax, nbins, shift, angularUnits))
    : move(Pair::Create(_comoving_multipoles_lin_, _standard_, rMin, rMax, nbins, shift, angularUnits));
}


// ============================================================================================


void cosmobl::measure::twopt::TwoPointCorrelation_multipoles_direct::set_parameters (const binType binType, const double rMin, const double rMax, const double binSize, const double shift, const CoordUnits angularUnits, function<double(double)> angularWeight, const bool compute_extra_info)
{
  if (!compute_extra_info) 
    m_dd = (binType==_logarithmic_) ? move(Pair::Create(_comoving_multipoles_log_, _standard_, rMin, rMax, binSize, shift, angularUnits, angularWeight))
      : move(Pair::Create(_comoving_multipoles_lin_, _standard_, rMin, rMax, binSize, shift, angularUnits, angularWeight));
  else 
    m_dd = (binType==_logarithmic_) ? move(Pair::Create(_comoving_multipoles_log_, _extra_, rMin, rMax, binSize, shift, angularUnits, angularWeight))
      : move(Pair::Create(_comoving_multipoles_lin_, _extra_, rMin, rMax, binSize, shift, angularUnits, angularWeight));
    
  m_rr = (binType==_logarithmic_) ? move(Pair::Create(_comoving_multipoles_log_, _standard_, rMin, rMax, binSize, shift, angularUnits))
    : move(Pair::Create(_comoving_multipoles_lin_, _standard_, rMin, rMax, binSize, shift, angularUnits));
  
  m_dr = (binType==_logarithmic_) ? move(Pair::Create(_comoving_multipoles_log_, _standard_, rMin, rMax, binSize, shift, angularUnits))
    : move(Pair::Create(_comoving_multipoles_lin_, _standard_, rMin, rMax, binSize, shift, angularUnits));
}


// ============================================================================================


vector<double> cosmobl::measure::twopt::TwoPointCorrelation_multipoles_direct::xx () const
{
  vector<double> rad, xx;
  m_dataset->xx(xx);

  for (size_t i=0; i<xx.size()/3; i++)
    rad.push_back(xx[i]);

  return rad;
}


// ============================================================================================


vector<double> cosmobl::measure::twopt::TwoPointCorrelation_multipoles_direct::xiMonopole () const
{
  vector<double> vv; 
  m_dataset->data(vv);

  size_t sz = vv.size();

  vector<double> xi0;
  for (size_t i=0; i<sz/3; i++)
    xi0.push_back(vv[i]);

  return xi0;
}


// ============================================================================================


vector<double> cosmobl::measure::twopt::TwoPointCorrelation_multipoles_direct::errorMonopole () const
{
  vector<double> vv; 
  m_dataset->error(vv);

  size_t sz = vv.size();

  vector<double> error_xi0;

  for (size_t i=0; i<sz/3; i++)
    error_xi0.push_back(vv[i]);

  return error_xi0;
}


// ============================================================================================


vector<double> cosmobl::measure::twopt::TwoPointCorrelation_multipoles_direct::xiQuadrupole () const 
{
  vector<double> vv; 
  m_dataset->data(vv);

  size_t sz = vv.size();

  vector<double> xi2;

  for (size_t i=sz/3; i<2*sz/3; i++)
    xi2.push_back(vv[i]);

  return xi2;
}


// ============================================================================================


vector<double> cosmobl::measure::twopt::TwoPointCorrelation_multipoles_direct::errorQuadrupole () const 
{
  vector<double> vv; 
  m_dataset->error(vv);

  size_t sz = vv.size();

  vector<double> error_xi2;

  for (size_t i=sz/3; i<2*sz/3; i++)
    error_xi2.push_back(vv[i]);

  return error_xi2;
}


// ============================================================================================


vector<double> cosmobl::measure::twopt::TwoPointCorrelation_multipoles_direct::xiHexadecapole () const
{
  vector<double> vv; 
  m_dataset->data(vv);

  size_t sz = vv.size();

  vector<double> xi4;

  for (size_t i=2*sz/3; i<sz; i++)
    xi4.push_back(vv[i]);

  return xi4;
}


// ============================================================================================


vector<double> cosmobl::measure::twopt::TwoPointCorrelation_multipoles_direct::errorHexadecapole () const 
{
  vector<double> vv; 
  m_dataset->error(vv);

  size_t sz = vv.size();

  vector<double> error_xi4;
  
  for (size_t i=2*sz/3; i<sz; i++)
    error_xi4.push_back(vv[i]);

  return error_xi4;
}


// ============================================================================================


void cosmobl::measure::twopt::TwoPointCorrelation_multipoles_direct::write (const string dir, const string file, const int rank) const 
{
  (void)rank;
  
  vector<double> rad; m_dataset->xx(rad);
  vector<double> xil = m_dataset->data();
  vector<double> error = m_dataset->error();

  checkDim(rad, m_dd->nbins()*3, "rad");
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  string header = "[1] separation at the bin centre # [2] monopole # [3] error on the monopole # [4] quadrupole # [5] error on the quadrupole # [6] hexadecapole # [7] error on the hexadecapole";
  
  if (m_compute_extra_info) header += " # [8] mean separation # [9] standard deviation of the separation distribution # [10] mean redshift # [11] standard deviation of the redshift distribution";
  
  fout << "### " << header << " ###" <<endl;

  for (int i=0; i<m_dd->nbins(); i++) {
    fout << setiosflags(ios::fixed) << setprecision(5) << setw(8) << rad[i] << "  " << setw(8) << xil[i] << "  " << setw(8) << error[i] << "  " << setw(8) << xil[i+m_dd->nbins()] << "  " << setw(8) << error[i+m_dd->nbins()] << "  " << setw(8) << xil[i+2*m_dd->nbins()] << "  " << setw(8) << error[i+2*m_dd->nbins()];
    if (m_compute_extra_info)
      for (size_t ex=0; ex<m_dataset->extra_info().size(); ++ex)
	fout << "   " << setw(8) << m_dataset->extra_info(ex, i);
    fout << endl;
  }
   
  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl << endl;
}
