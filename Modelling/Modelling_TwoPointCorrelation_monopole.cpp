/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Modelling/Modelling_TwoPointCorrelation_monopole.cpp
 *
 *  @brief Methods of the class
 *  Modelling_TwoPointCorrelation_monopole, used to model the monopole
 *  of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation_monopole
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "Modelling_TwoPointCorrelation_monopole.h"

using namespace cosmobl;


// ============================================================================================


void cosmobl::modelling::Modelling_TwoPointCorrelation_monopole::set_fiducial_xiDM ()
{
  coutCBL << "Setting up the fiducial two-point correlation function model" << endl;
  
  m_twop_parameters.fiducial_xiDM.erase(m_twop_parameters.fiducial_xiDM.begin(), m_twop_parameters.fiducial_xiDM.end());
  
  if (m_twop_parameters.sigmaNL==0) {
    for (size_t i=0; i<m_twop_parameters.fiducial_radDM.size(); i++) 
      m_twop_parameters.fiducial_xiDM.push_back(m_twop_parameters.cosmology->xi_DM(m_twop_parameters.fiducial_radDM[i], m_twop_parameters.method_Pk, m_twop_parameters.redshift, m_twop_parameters.output_root, m_twop_parameters.NL, m_twop_parameters.norm, m_twop_parameters.k_min, m_twop_parameters.k_max, m_twop_parameters.aa, m_twop_parameters.GSL, m_twop_parameters.prec, m_twop_parameters.file_par));
  }
  
  else {
    vector<double> kk = logarithmic_bin_vector(500, m_twop_parameters.k_min+1.e-4, m_twop_parameters.k_max), Pk;

    for (size_t i=0; i<kk.size(); i++)
      Pk.push_back(m_twop_parameters.cosmology->Pk_DeWiggle(kk[i], m_twop_parameters.redshift, m_twop_parameters.sigmaNL, m_twop_parameters.output_root, m_twop_parameters.norm, m_twop_parameters.k_min, m_twop_parameters.k_max, m_twop_parameters.aa, m_twop_parameters.prec));

    m_twop_parameters.fiducial_xiDM = Xi0(m_twop_parameters.fiducial_radDM, kk, Pk);
  }
  
  m_twop_parameters.func_xi = make_shared<classfunc::func_grid_GSL>(classfunc::func_grid_GSL(m_twop_parameters.fiducial_radDM, m_twop_parameters.fiducial_xiDM, "Spline"));

  m_twop_parameters.sigma8_z = m_twop_parameters.cosmology->sigma8_Pk(m_twop_parameters.method_Pk, m_twop_parameters.redshift, m_twop_parameters.output_root);
  
  m_twop_parameters.linear_growth_rate_z = m_twop_parameters.cosmology->linear_growth_rate(m_twop_parameters.redshift);

}
