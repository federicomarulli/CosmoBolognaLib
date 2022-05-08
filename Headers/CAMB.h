/*******************************************************************
 *  Copyright (C) 2022 by Federico Marulli and Giorgio Lesci       *
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
 *  @file Headers/CAMB.h
 *
 *  @brief class CAMB that wrap CAMB routines
 *
 *  This file contains the wrappers of CAMB routines
 *
 *  @author Giorgio Lesci
 *
 *  @author giorgio.lesci2@unibo.it
 */

#ifndef __CAMB__
#define __CAMB__

#include "Kernel.h"

namespace cbl {

  namespace wrapper {

    /**
     *  @brief The namespace of the <B> CAMB wrappers </B>
     *  
     *  The \e camb namespace contains all the wrapper functions of the
     *  CAMB routines
     */
    namespace camb {

      extern"C"
      {

	/**
	 *  @brief wrapper of the function GetCAMBparams
	 *  in CAMB/CAMBinterface.f90, creating an instance
	 *  of the CAMBparams type contained in
	 *  External/CAMB/fortran/model.f90
	 *
	 *  @return pointer to CAMBparams instance
	 */
	void* GetCAMBparams ();
	
	/**
	 *  wrapper of the function GetCAMBdata
	 *  in CAMB/CAMBinterface.f90, creating an instance
	 *  of the CAMBdata type contained in
	 *  External/CAMB/fortran/results.f90
	 *
	 *  @return pointer to CAMBdata instance
	 */
	void* GetCAMBdata ();
	
	/**
	 *  @brief wrapper of the subroutine SetCAMBparams
	 *  in CAMB/CAMBinterface.f90, setting CAMBparams
	 *  parameters
	 *
	 *  @param params pointer to CAMBparams instance
	 *
	 *  @param ombh2 \f$\Omega_{\rm b}h^2\f$
	 *
	 *  @param omch2 \f$\Omega_{\rm c}h^2\f$
	 *
	 *  @param omnuh2 \f$\Omega_{\nu}h^2\f$
	 *
	 *  @param massless_nu contribution of massless neutrinos
	 *  to \f$N_{\rm eff}\f$
	 *
	 *  @param massive_nu number of massive neutrinos
	 *
	 *  @param neutrino_hierarchy 1 \f$\rightarrow\f$ normal
	 *  hierarchy; 2 \f$\rightarrow\f$ inverted hierarchy;
	 *  3 \f$\rightarrow\f$ degenerate hierarchy
	 *
	 *  @param omk \f$\Omega_{\rm k}\f$
	 *
	 *  @param H0 \f$H_0\f$
	 *
	 *  @param dark_energy_model 0 \f$\rightarrow\f$ fluid;
	 *  1 \f$\rightarrow\f$ PPF; 2 \f$\rightarrow\f$ axion
	 *  effective fluid; 3 \f$\rightarrow\f$ early 
	 *  quintessence
	 *
	 *  @param w \f$w_0\f$
	 *
	 *  @param wa \f$w_{\rm a}\f$
	 *
	 *  @param tau \f$\tau\f$
	 *
	 *  @param cs2_lam constant comoving sound speed 
	 *  of the dark energy (1=quintessence), for fluid model
	 *
	 *  @param T_cmb CMB temperature
	 *
	 *  @param helium_fraction helium fraction
	 */
	void SetCAMBparams (void* params, const double ombh2, const double omch2, const double omnuh2, const double massless_nu, const int massive_nu, const int neutrino_hierarchy, const double omk, const double H0, const int dark_energy_model, const double w, const double wa, const double tau, const double cs2_lam, const double T_cmb, const double helium_fraction);

	/**
	 *  @brief wrapper of the subroutine SetCAMBPk
	 *  in CAMB/CAMBinterface.f90, setting the
	 *  power spectrum parameters
	 *
	 *  @param params pointer to CAMBparams instance
	 *
	 *  @param redshift the redshift
	 *
	 *  @param ns \f$n_s\f$
	 *
	 *  @param As \f$A_s\f$
	 *
	 *  @param pivot_scalar the pivot scalar
	 *
	 *  @param accurate_massive_nu whether you care 
	 *  about accuracy of the neutrino transfers themselves
	 *
	 *  @param kmax maximum \f$k\f$
	 *
	 *  @param nonlinear if true, \f$P(k)\f$ is nonlinear
	 */
	void SetCAMBPk (void* params, const double redshift, const double ns, const double As, const double pivot_scalar, const bool accurate_massive_nu, const double kmax, const bool nonlinear);

	/**
	 *  @brief wrapper of the subroutine GetCAMBresults
	 *  in CAMB/CAMBinterface.f90, calling the
	 *  subroutine CAMB_GetResults contained in
	 *  External/CAMB/fortran/camb.f90
	 *
	 *  @param params pointer to CAMBparams instance
	 *
	 *  @param data pointer to CAMBdata instance
	 */
	void GetCAMBresults (void* params, void* data);

	/**
	 *  @brief wrapper of the subroutine GetCAMBPk
	 *  in CAMB/CAMBinterface.f90, calling the
	 *  subroutine Transfer_GetMatterPowerD
	 *  contained in
	 *  External/CAMB/fortran/results.f90
	 *
	 *  @param data pointer to CAMBdata instance
	 *
	 *  @param Pk output power spectrum
	 * 
	 *  @param mink minimum \f$k\f$
	 * 
	 *  @param dlnkh log space between each \f$k\f$
	 * 
	 *  @param npoints number of points where \f$P(k)\f$
	 *  is computed
	 */
	void GetCAMBPk (void* data, double *Pk, const double mink, const double dlnkh, const int npoints);

	/**
	 *  @brief wrapper of the subroutine ReleaseCAMBparams
	 *  in CAMB/CAMBinterface.f90, deallocating 
	 *  CAMBparams instance
	 *
	 *  @param params pointer to CAMBparams instance
	 *
	 *  @return empty pointer
	 */
	void* ReleaseCAMBparams (void* params);

	/**
	 *  @brief wrapper of the subroutine ReleaseCAMBdata
	 *  in CAMB/CAMBinterface.f90, deallocating 
	 *  CAMBdata instance
	 *
	 *  @param data pointer to CAMBdata instance
	 *
	 *  @return empty pointer
	 */
	void* ReleaseCAMBdata (void* data);
	
      }
      
      /**
       *  @brief Get the matter power spectrum from CAMB
       *
       *  @param nonlinear if true, \f$P(k)\f$ is nonlinear
       *
       *  @param redshift the redshift
       * 
       *  @param kmin minimum \f$k\f$
       * 
       *  @param kmax maximum \f$k\f$
       * 
       *  @param npoints number of points where \f$P(k)\f$
       *  is computed
       *
       *  @param ombh2 \f$\Omega_{\rm b}h^2\f$
       *
       *  @param omch2 \f$\Omega_{\rm c}h^2\f$
       *
       *  @param omnuh2 \f$\Omega_{\nu}h^2\f$
       *
       *  @param massless_nu contribution of massless neutrinos
       *  to \f$N_{\rm eff}\f$
       *
       *  @param massive_nu number of massive neutrinos
       *
       *  @param omk \f$\Omega_{\rm k}\f$
       *
       *  @param H0 \f$H_0\f$
       *
       *  @param ns \f$n_s\f$
       *
       *  @param As \f$A_s\f$
       *
       *  @param pivot_scalar the pivot scalar
       *
       *  @param w \f$w_0\f$
       *
       *  @param wa \f$w_{\rm a}\f$
       *
       *  @param tau \f$\tau\f$
       *
       *  @param accurate_massive_nu whether you care 
       *  about accuracy of the neutrino transfers themselves
       *
       *  @param neutrino_hierarchy 1 \f$\rightarrow\f$ normal
       *  hierarchy; 2 \f$\rightarrow\f$ inverted hierarchy;
       *  3 \f$\rightarrow\f$ degenerate hierarchy
       *
       *  @param dark_energy_model 0 \f$\rightarrow\f$ fluid;
       *  1 \f$\rightarrow\f$ PPF; 2 \f$\rightarrow\f$ axion
       *  effective fluid; 3 \f$\rightarrow\f$ early 
       *  quintessence
       *
       *  @param cs2_lam constant comoving sound speed 
       *  of the dark energy (1=quintessence), for fluid model
       *
       *  @param T_cmb CMB temperature
       *
       *  @param helium_fraction helium fraction
       *
       *  @return CAMB power spectrum
       */
      std::vector<double> Pk_CAMB (const bool nonlinear, const double redshift, const double kmin, const double kmax, const int npoints, const double ombh2, const double omch2, const double omnuh2, const double massless_nu, const int massive_nu, const double omk, const double H0, const double ns, const double As, const double pivot_scalar=0.05, const double w=-1., const double wa=0., const double tau=2.1e-9, const bool accurate_massive_nu=false, const int neutrino_hierarchy=3, const int dark_energy_model=1, const double cs2_lam=1., const double T_cmb=2.7255, const double helium_fraction=0.24);
      
    }
  }
}
 
#endif
