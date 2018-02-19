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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation_multipoles.h
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


namespace cosmobl {

  namespace modelling {

    namespace twopt {

      /**
       *  @class Modelling_TwoPointCorrelation_multipoles
       *  Modelling_TwoPointCorrelation_multipoles.h
       *  "Headers/Lib/Modelling_TwoPointCorrelation_multipoles.h"
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
	vector<int> m_multipoles_order;

	/// check if the model has been set
	bool m_ModelIsSet;
	
	/**
	 *  @name free parameters to model the two-point correlation function
	 */
	///@{
	
	/// \f$ \alpha_{\perp} \f$
	shared_ptr<statistics::Parameter> m_alpha_perpendicular;
	
	/// \f$ \alpha_{\parallel} \f$
	shared_ptr<statistics::Parameter> m_alpha_parallel;

	/// \f$b\sigma_8\f$
	shared_ptr<statistics::Parameter> m_bsigma8;

	/// \f$f\sigma_8\f$
	shared_ptr<statistics::Parameter> m_fsigma8;

	/// \f$f\Sigma_S\f$ the streaming scale
	shared_ptr<statistics::Parameter> m_SigmaS;

	/// \f$B_0\f$ nuisance parameter for anisotropic analysis 
	shared_ptr<statistics::Parameter> m_B0;

	/// \f$B_2\f$ nuisance parameter for anisotropic analysis 
	shared_ptr<statistics::Parameter> m_B2;
	
	/// \f$A_0^0\f$ nuisance parameter for anisotropic analysis 
	shared_ptr<statistics::Parameter> m_A00;
	
	/// \f$A_1^0\f$ nuisance parameter for anisotropic analysis 
	shared_ptr<statistics::Parameter> m_A01;
	
	/// \f$A_2^0\f$ nuisance parameter for anisotropic analysis 
	shared_ptr<statistics::Parameter> m_A02;
	
	/// \f$A_0^2\f$ nuisance parameter for anisotropic analysis 
	shared_ptr<statistics::Parameter> m_A20;
	
	/// \f$A_1^2\f$ nuisance parameter for anisotropic analysis 
	shared_ptr<statistics::Parameter> m_A21;

	/// \f$A_2^2\f$ nuisance parameter for anisotropic analysis 
	shared_ptr<statistics::Parameter> m_A22;

	///@}

	
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
	Modelling_TwoPointCorrelation_multipoles (const shared_ptr<cosmobl::measure::twopt::TwoPointCorrelation> twop);

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
	Modelling_TwoPointCorrelation_multipoles (const shared_ptr<data::Data> twop_dataset, const int nmultipoles);
      
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
	void set_fit_range (vector<vector<double>> fit_range);

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
	 *  @brief return the protected member 
	 *  alpha_perpendicular
	 *
	 *  @return pointer to parameter alpha_perpendicular
	 */
	shared_ptr<statistics::Parameter> alpha_perpendicular () { return m_alpha_perpendicular; }

	/**
	 *  @brief return the protected member 
	 *  m_alpha_parallel
	 *
	 *  @return pointer to parameter alpha_parallel
	 */
	shared_ptr<statistics::Parameter> alpha_parallel () { return m_alpha_parallel; }

	/**
	 *  @brief return the protected member 
	 *  m_bsigma8
	 *
	 *  @return pointer to parameter bsigma8
	 */
	shared_ptr<statistics::Parameter> bsigma8 () { return m_bsigma8; }

	/**
	 *  @brief return the protected member 
	 *  m_fsigma8
	 *
	 *  @return pointer to parameter fsigma8
	 */
	shared_ptr<statistics::Parameter> fsigma8 () { return m_fsigma8; }

	/**
	 *  @brief return the protected member 
	 *  m_SigmaS
	 *
	 *  @return pointer to parameter SigmaS
	 */
	shared_ptr<statistics::Parameter> SigmaS () { return m_SigmaS; }

	/**
	 *  @brief return the protected member 
	 *  m_B0
	 *
	 *  @return pointer to parameter B0
	 */
	shared_ptr<statistics::Parameter> B0 () { return m_B0; }

	/**
	 *  @brief return the protected member 
	 *  m_B0
	 *
	 *  @return pointer to parameter B0
	 */
	shared_ptr<statistics::Parameter> B2 () { return m_B2; }
	
	/**
	 *  @brief return the protected member 
	 *  m_A00
	 *
	 *  @return pointer to parameter A00
	 */
	shared_ptr<statistics::Parameter> A00 () { return m_A00; }
	
	/**
	 *  @brief return the protected member 
	 *  m_A00
	 *
	 *  @return pointer to parameter A01
	 */
	shared_ptr<statistics::Parameter> A01 () { return m_A01; }
	
	/**
	 *  @brief return the protected member 
	 *  m_A02
	 *
	 *  @return pointer to parameter A02
	 */
	shared_ptr<statistics::Parameter> A02 () { return m_A02; }
	
	/**
	 *  @brief return the protected member 
	 *  m_A20
	 *
	 *  @return pointer to parameter A20
	 */
	shared_ptr<statistics::Parameter> A20 () { return m_A20; }
	
	/**
	 *  @brief return the protected member 
	 *  m_A21
	 *
	 *  @return pointer to parameter A21
	 */
	shared_ptr<statistics::Parameter> A21 () { return m_A21; }

	/**
	 *  @brief return the protected member 
	 *  m_A22
	 *
	 *  @return pointer to parameter A22
	 */
	shared_ptr<statistics::Parameter> A22 () { return m_A22; }

	/**
	 *  @brief set the parameters used to model the full shape
	 *  of the multipoles of the two-point correlation function
	 *
	 *  the model is the following:
	 *
	 *  The functions computes the multipoles of the two-point
	 *  correlation function
	 *
	 *  \f[ \xi_l(s) = \frac{i^l}{2\pi^2} \int \mathrm{d} k P_l(k)
	 *  j_l(ks); \f]
	 *
	 *  where \f$j_l(ks)\f$ are the bessel functions.
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
	void set_model_fullShape (const statistics::Prior alpha_perpendicular_prior={}, const statistics::Prior alpha_parallel_prior={}, statistics::Prior fsigma8_prior={}, statistics::Prior bsigma8_prior={}, const statistics::Prior SigmaS_prior={}, const bool compute_PkDM=true);


	/**
	 * @brief set the parameters used to model the two-point 
	 * correlation function, intended for anisotropic BAO
	 * measurements (Ross et al. 2017). In only works with
	 * monopole and quadrupole.
	 *
	 * The functions computes the multipoles of the 
	 * two-point correlation function:
	 *
	 * \f[ \xi_0(s) = B_0\xi_0(s, \alpha_{\perp},
	 * \alpha_{\parallel})+A_0^0+\frac{A_0^1}{s}+\frac{A_0^2}{s^2}
	 * \f]
	 *
	 * \f[ \xi_2(s) = \frac{5}{2}\left[B_2\xi_{\mu2}(s,
	 * \alpha_{\perp}, \alpha_{\parallel})-B_0\xi_0(s,
	 * \alpha_{\perp}, \alpha_{\parallel})\right]
	 * +A_2^0+\frac{A_2^1}{s}+\frac{A_2^2}{s^2} \f]
	 *
	 * where \f$\xi_0(s, \alpha_{\perp}, \alpha_{\parallel})\f$ is
	 * the monopole computed at the fiducial cosmology,
	 * \f$\xi_{\mu2}(s, \alpha_{\perp}, \alpha_{\parallel}) =
	 * 3\int_0^1\mathrm{d}\mu\mu^2\xi(s, \mu, \alpha_{\perp},
	 * \alpha_{\parallel})\f$.
	 *
	 * The function takes as inputs ten parameters
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
	 */
	void set_model_BAO (const statistics::Prior alpha_perpendicular_prior={}, const statistics::Prior alpha_parallel_prior={}, const statistics::Prior B0_prior={}, const statistics::Prior B2_prior={}, const statistics::Prior A00_prior={}, const statistics::Prior A20_prior={}, const statistics::Prior A01_prior={}, const statistics::Prior A21_prior={}, const statistics::Prior A02_prior={}, const statistics::Prior A22_prior={}, const bool compute_XiDM=true);


	/**
	 *  @brief set the parameters used to model the full shape
	 *  of the multipoles of the two-point correlation function
	 *
	 *  the model is the following:
	 *
	 *  The functions computes the multipoles of the two-point
	 *  correlation function
	 *
	 *  \f[ \xi_l(s) = \frac{i^l}{2\pi^2} \int \mathrm{d} k P_l(k)
	 *  j_l(ks); \f]
	 *
	 *  where \f$j_l(ks)\f$ are the bessel functions.
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
	void set_model_sigma8_bias (const statistics::Prior sigma8_prior={}, const statistics::Prior bias_prior={});
	
        /**
         *  @brief compute the best-fit model from MCMC chains
         *
         *  @param xx vector of points at which the model is computed
	 *  @param start the starting position for each chain
	 *  @param thin the position step
         *  @param median_model the median model
         *  @param low_model the model at 16th percentile
         *  @param up_model the model at 84th percentile
         *  @return none
         */
        void compute_model_from_chains (const vector<double> xx, const int start, const int thin, vector<double> &median_model, vector<double> &low_model, vector<double> &up_model) override;	
	
	/**
	 *  @brief compute and write the model
	 *
	 *  @param dir the output directory
	 *
	 *  @param file the output file
	 *
	 *  @param xx vector of points at which the model is computed
	 *
	 *  @param parameter vector containing the input parameters
	 *  used to compute the model; if this vector is not provided,
	 *  the model will be computed using the best-fit parameters
	 *
	 *  @param start the starting position for each chain
	 *
	 *  @param thin the position step
	 *
	 *  @return none
	 */
	void write_model (const string dir, const string file, const vector<double> xx={}, const vector<double> parameter={}, const int start=0, const int thin=1);
	
      };
    }
  }
}

#endif
