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
 *  @file Headers/Modelling_TwoPointCorrelation_multipoles.h
 *
 *  @brief The class Modelling_TwoPointCorrelation_multipoles
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation_multipoles, used to model the
 *  multipoles two-point correlation function
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGTWOPOINTMULT__
#define __MODELLINGTWOPOINTMULT__


#include "Modelling_TwoPointCorrelation1D_monopole.h"
#include "ModelFunction_TwoPointCorrelation_multipoles.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    namespace twopt {

      /**
       *  @class Modelling_TwoPointCorrelation_multipoles
       *  Modelling_TwoPointCorrelation_multipoles.h
       *  "Headers/Modelling_TwoPointCorrelation_multipoles.h"
       *
       *  @brief The class Modelling_TwoPointCorrelation_multipoles
       *
       *  This file defines the interface of the base class
       *  Modelling_TwoPointCorrelation_multipoles, used to model the
       *  multipoles of the two-point correlation function
       *
       */
      class Modelling_TwoPointCorrelation_multipoles : public Modelling_TwoPointCorrelation1D_monopole {

      protected:

	/// number of multipoles
	int m_nmultipoles;

	/// number of multipoles used for the fit
	int m_nmultipoles_fit;

	/// vector containing the ordering of the data vector
	std::vector<int> m_multipoles_order;

	/// check if the model has been set
	bool m_ModelIsSet;
	
      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class
	 *  ModellingTwoPointCorrelation_multipoles
	 */
	Modelling_TwoPointCorrelation_multipoles () = default;
      
	/**
	 *  @brief constructor
	 *  
	 *  @param twop the two-point correlation function to model
	 *
	 *  @return object of type
	 *  Modelling_TwoPointCorrelation_multipoles
	 */
	Modelling_TwoPointCorrelation_multipoles (const std::shared_ptr<cbl::measure::twopt::TwoPointCorrelation> twop);

	/**
	 *  @brief constructor
	 *  
	 *  @param twop_dataset the dataset containing the two-point
	 *  correlation function to model
	 *
	 *  @param nmultipoles the number of multipoles
	 *
	 *  @return object of type
	 *  Modelling_TwoPointCorrelation_multipoles
	 */
	Modelling_TwoPointCorrelation_multipoles (const std::shared_ptr<data::Data> twop_dataset, const int nmultipoles);
      
	/**
	 *  @brief default destructor
	 *  @return none
	 */
	~Modelling_TwoPointCorrelation_multipoles () = default;
	
	///@}

	
	/**
	 *  @name Member functions used to set the model parameters
	 */
	///@{

	/**
	 *  @brief set the scale range used for the fit
	 *
	 *  @param xmin the minimum x value
	 *  @param xmax the maximum x value
	 *  @param nmultipoles the number of multipoles
	 *
	 *  @return none
	 */
	void set_fit_range (const double xmin, const double xmax, const int nmultipoles=-1);

	/**
	 *  @brief set the scale range used for the fit
	 *
	 *  @param fit_range vector containing the fitting range for
	 *  the multipoles
	 *
	 *  @return none
	 */
	void set_fit_range (std::vector<std::vector<double>> fit_range);

	/**
	 *  @brief set the fiducial model for the dark matter power
	 *  spectrum
	 *
	 *  @return none
	 */
	void set_fiducial_PkDM ();
	
	/**
	 *  @brief set the fiducial model for the dark matter
	 *  two-point correlation function and associated quantities
	 *
	 *  @return none
	 */
	void set_fiducial_xiDM ();

	/**
	 *  @brief set the model to fit the full shape of the
	 *  multipole moments of the two-point correlation function
	 *
	 *  the multipoles of the two-point correlation function will
	 *  be computed as follows:
	 *
	 *  \f[ \xi_l(s) = \frac{i^l}{2\pi^2} \int \mathrm{d} k P_l(k)
	 *  j_l(ks) \f]
	 *
	 *  where \f$j_l(ks)\f$ are the bessel functions and
	 *  \f$P_l(k)\f$ is computed by cbl::modelling::twopt::Pk_l
	 *
	 *  the model has 3+N parameters: 
	 *    - \f$\alpha_{\perp}\f$
	 *    - \f$\alpha_{\parallel}\f$
	 *    - \f$f(z)\sigma_8(z)\f$
	 *    - \f$b(z)\sigma_8(z)\f$
	 *    - \f$\Sigma_s\f$
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters
	 *
	 *  @param alpha_perpendicular_prior prior for the parameter 
	 *  \f$\alpha_{\perp}\f$
	 *
	 *  @param alpha_parallel_prior prior for the parameter 
	 *  \f$\alpha_{\parallel}\f$
	 *
	 *  @param fsigma8_prior prior for the parameter
	 *  \f$f(z)\sigma_8(z)\f$
	 *
	 *  @param bsigma8_prior prior for the parameter
	 *  \f$b(z)\sigma_8(z)\f$
	 *
	 *  @param SigmaS_prior prior for the parameters
	 *  \f$\Sigma_S\f$
	 *
	 *  @param compute_PkDM true \f$\rightarrow\f$ compute the
	 *  fiducial model of the dark matter power spectrum
	 *
	 *  @return none
	 */
	void set_model_fullShape (const statistics::PriorDistribution alpha_perpendicular_prior={}, const statistics::PriorDistribution alpha_parallel_prior={}, statistics::PriorDistribution fsigma8_prior={}, statistics::PriorDistribution bsigma8_prior={}, const statistics::PriorDistribution SigmaS_prior={}, const bool compute_PkDM=true);

	/**
	 *  @brief set the model to fit the full shape of the
	 *  multipole moments of the two-point correlation function
	 *
	 *  the multipoles of the two-point correlation function will
	 *  be computed as follows:
	 *
	 *  \f[ \xi_l(s) = \frac{i^l}{2\pi^2} \int \mathrm{d} k P_l(k)
	 *  j_l(ks); \f]
	 *
	 *  where \f$j_l(ks)\f$ are the bessel functions and
	 *  \f$P_l(k)\f$ is computed by cbl::modelling::twopt::Pk_l
	 *
	 *  the model has 2 parameters: 
	 *    - \f$\sigma_8(z)\f$
	 *    - \f$b(z)\f$
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters
	 *
	 *  @param sigma8_prior prior for the parameter
	 *  \f$\sigma_8(z)\f$
	 *
	 *  @param bias_prior prior for the parameter bias
	 *  \f$b(z)\f$
	 *
	 *  @return none
	 */
	void set_model_fullShape_sigma8_bias (const statistics::PriorDistribution sigma8_prior={}, const statistics::PriorDistribution bias_prior={});
	
	/**
	 *  @brief set the model to fit the monopole and quadrupole of
	 *  the two-point correlation function, used for anisotropic
	 *  BAO measurements
	 *
	 *  the monopole and quadrupole of the two-point correlation
	 *  function are computed as follows (Ross et al. 2017):
	 *
	 *  \f[ \xi_0(s) = B_0\xi_0(s, \alpha_{\perp},
	 *  \alpha_{\parallel})+A_0^0+\frac{A_0^1}{s}+\frac{A_0^2}{s^2}
	 *  \f]
	 *
	 *  \f[ \xi_2(s) = \frac{5}{2}\left[B_2\xi_{\mu2}(s,
	 *  \alpha_{\perp}, \alpha_{\parallel})-B_0\xi_0(s,
	 *  \alpha_{\perp}, \alpha_{\parallel})\right]
	 *  +A_2^0+\frac{A_2^1}{s}+\frac{A_2^2}{s^2} \f]
	 *
	 *  where \f$\xi_0(s, \alpha_{\perp}, \alpha_{\parallel})\f$
	 *  is the monopole computed at the fiducial cosmology,
	 *  \f$\xi_{\mu2}(s, \alpha_{\perp}, \alpha_{\parallel}) =
	 *  3\int_0^1\mathrm{d}\mu\mu^2\xi(s, \mu, \alpha_{\perp},
	 *  \alpha_{\parallel})\f$.
	 *
	 *  The function takes as inputs ten parameters
	 *    - \f$\alpha_{\perp}\f$
	 *    - \f$\alpha_{\parallel}\f$
	 *    - \f$B_0\f$
	 *    - \f$B_2\f$
	 *    - \f$A^0_0\f$
	 *    - \f$A^2_0\f$
	 *    - \f$A^0_1\f$
	 *    - \f$A^2_1\f$
	 *    - \f$A^0_2\f$
	 *    - \f$A^2_2\f$
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters
	 *
	 *  @param alpha_perpendicular_prior prior for the parameter 
	 *  \f$\alpha_{\perp}\f$
	 *
	 *  @param alpha_parallel_prior prior for the parameter
	 *  \f$\alpha_{\parallel}\f$
	 *  
	 *  @param B0_prior prior for the parameter \f$B_0\f$
	 *
	 *  @param B2_prior prior for the parameter \f$B_2\f$
	 *
	 *  @param A00_prior prior for the parameter \f$A^0_0\f$
	 *
	 *  @param A20_prior prior for the parameter \f$A^2_0\f$
	 *
	 *  @param A01_prior prior for the parameter \f$A^0_1\f$
	 *
	 *  @param A21_prior prior for the parameter \f$A^2_1\f$
	 *
	 *  @param A02_prior prior for the parameter \f$A^0_2\f$
	 *
	 *  @param A22_prior prior for the parameter \f$A^2_2\f$
	 *
	 *  @param compute_XiDM true \f$\rightarrow\f$ compute the
	 *  fiducial model of the dark matter two-point correlation
	 *  function
	 *
	 *  @return none
	 *
	 *  @warning the current implementation works only for
	 *  monopole and quadrupole
	 */
	void set_model_BAO (const statistics::PriorDistribution alpha_perpendicular_prior={}, const statistics::PriorDistribution alpha_parallel_prior={}, const statistics::PriorDistribution B0_prior={}, const statistics::PriorDistribution B2_prior={}, const statistics::PriorDistribution A00_prior={}, const statistics::PriorDistribution A20_prior={}, const statistics::PriorDistribution A01_prior={}, const statistics::PriorDistribution A21_prior={}, const statistics::PriorDistribution A02_prior={}, const statistics::PriorDistribution A22_prior={}, const bool compute_XiDM=true);
	
      };
    }
  }
}

#endif
