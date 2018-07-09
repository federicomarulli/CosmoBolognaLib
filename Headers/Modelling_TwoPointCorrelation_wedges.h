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
 *  @file Headers/Modelling_TwoPointCorrelation_wedges.h
 *
 *  @brief The class Modelling_TwoPointCorrelatoin_wedges
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation_wedges, used to model the wedges of
 *  the two-point correlation function
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGTWOPOINTWED__
#define __MODELLINGTWOPOINTWED__


#include "Modelling_TwoPointCorrelation1D_monopole.h"


// ===================================================================================================


namespace cbl {
  
  namespace modelling {

    namespace twopt {
    
      /**
       *  @class Modelling_TwoPointCorrelation_wedges
       *  Modelling_TwoPointCorrelation_wedges.h
       *  "Headers/Modelling_TwoPointCorrelation_wedges.h"
       *
       *  @brief The class Modelling_TwoPointCorrelation_wedges
       *
       *  This file defines the interface of the base class
       *  Modelling_TwoPointCorrelation_wedges, used for modelling the
       *  wedges of the two-point correlation function
       *
       */
      class Modelling_TwoPointCorrelation_wedges : public Modelling_TwoPointCorrelation1D_monopole {

      protected:

	/// number of multipoles
	int m_nwedges;

	/// number of wedges used for the fit
	int m_nwedges_fit;

	/// vector containing the ordering of the data vector
	std::vector<int> m_wedges_order;

	/// the wedges aperture
	double m_deltamu;

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class Modelling_TwoPointCorrelation_wedges
	 */
	Modelling_TwoPointCorrelation_wedges () = default;

	/**
	 *  @brief constructor
	 *  
	 *  @param twop the two-point correlation function to model
	 *
	 *  @return object of type
	 *  Modelling_TwoPointCorrelation_wedges
	 */
	Modelling_TwoPointCorrelation_wedges (const std::shared_ptr<cbl::measure::twopt::TwoPointCorrelation> twop);

	/**
	 *  @brief constructor
	 *  
	 *  @param twop_dataset the dataset containing the two-point
	 *  correlation function to model
	 *  
	 *  @param nwedges the number of wedges
	 *
	 *  @return object of type
	 *  Modelling_TwoPointCorrelation_wedges
	 */
	Modelling_TwoPointCorrelation_wedges (const std::shared_ptr<data::Data> twop_dataset, const int nwedges);
      
	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_TwoPointCorrelation_wedges () = default;
	
	///@}

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
	 *  @brief set the scale range used for the fit
	 *
	 *  @param xmin the minimum x value
	 *  @param xmax the maximum x value
	 *  @param nwedges the number of wedges
	 *
	 *  @return none
	 */
	void set_fit_range (const double xmin, const double xmax, const int nwedges=-1);

	/**
	 *  @brief set the scale range used for the fit
	 *
	 *  @param fit_range vector containing the fitting range for
	 *  the wedges
	 *
	 *  @return none
	 */
	void set_fit_range (std::vector<std::vector<double>> fit_range);

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
	 *  the model has 3+N parameters: 
	 *    - \f$\alpha_{\perp}\f$
	 *    - \f$\alpha_{\parallel}\f$
	 *    - \f$f(z)\sigma_8(z)\f$
	 *    - \f$b(z)\sigma_8(z)\f$
	 *    - \f$\Sigma_s\f$*
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
	 *  @param compute_PkDM true \f$rightarrow\f$ compute the
	 *  fiducial model of the dark matter power spectrum
	 *
	 *  @return none
	 */
	void set_model_fullShape (const statistics::PriorDistribution alpha_perpendicular_prior={}, const statistics::PriorDistribution alpha_parallel_prior={}, statistics::PriorDistribution fsigma8_prior={}, statistics::PriorDistribution bsigma8_prior={}, const statistics::PriorDistribution SigmaS_prior={}, const bool compute_PkDM=true);

	/**
	 *  @brief set the model to fit the monopole and quadrupole of
	 *  the two-point correlation function, used for anisotropic
	 *  BAO measurements
	 *
	 *  the monopole and quadrupole of the two-point correlation
	 *  function are computed as follows (Ross et al. 2017):
	 *
	 *  \f[ \xi_{\perp}(s) = B_{perp}\xi_{\perp}(s,
	 *  \alpha_{\perp},
	 *  \alpha_{\parallel})+A_{perp}^0+\frac{A_{perp}^1}{s}+\frac{A_{perp}^2}{s^2};
	 *  \\ \xi_{\parallel}(s) = B_{parallel}\xi_{\parallel}(s,
	 *  \alpha_{\parallel},
	 *  \alpha_{\parallel})+A_{parallel}^0+\frac{A_{parallel}^1}{s}+\frac{A_{parallel}^2}{s^2};
	 *  \\ \f]
	 *
	 *  where \f$\xi_{\perp}\f$, \f$\xi_{\parallel}\f$ are the two
	 *  wedges of the two-point correlation function.
	 *
	 *  The function takes as inputs ten parameters
	 *    - \f$\alpha_{\perp}\f$
	 *    - \f$\alpha_{\parallel}\f$
	 *    - \f$B_0\f$
	 *    - \f$B_2\f$
	 *    - \f$A^0_0\f$
	 *    - \f$A^0_1\f$
	 *    - \f$A^0_2\f$
	 *    - \f$A^2_0\f$
	 *    - \f$A^2_1\f$
	 *    - \f$A^2_2\f$
	 *
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
	 *  @param Bperp_prior prior for the parameter
	 *  \f$B_perp\f$
	 *
	 *  @param Bpar_prior prior for the parameter
	 *  \f$B_par\f$
	 *
	 *  @param Aperp0_prior prior for the parameter
	 *  \f$A^0_0\f$
	 *
	 *  @param Apar0_prior prior for the parameter
	 *  \f$A^2_0\f$
	 *
	 *  @param Aperp1_prior prior for the parameter
	 *  \f$A^0_1\f$
	 *
	 *  @param Apar1_prior prior for the parameter
	 *  \f$A^2_1\f$
	 *
	 *  @param Aperp2_prior prior for the parameter
	 *  \f$A^0_2\f$
	 *
	 *  @param Apar2_prior prior for the parameter
	 *  \f$A^2_2\f$
	 *
	 *  @param compute_XiDM true \f$rightarrow\f$ compute the
	 *  fiducial model of the dark matter two-point correlation 
	 *  function
	 *
	 *  @return none
	 *
	 *  @warning the current implementation works only for
	 *  monopole and quadrupole
	 */
	void set_model_BAO (const statistics::PriorDistribution alpha_perpendicular_prior={}, const statistics::PriorDistribution alpha_parallel_prior={}, const statistics::PriorDistribution Bperp_prior={}, const statistics::PriorDistribution Bpar_prior={}, const statistics::PriorDistribution Aperp0_prior={}, const statistics::PriorDistribution Apar0_prior={}, const statistics::PriorDistribution Aperp1_prior={}, const statistics::PriorDistribution Apar1_prior={}, const statistics::PriorDistribution Aperp2_prior={}, const statistics::PriorDistribution Apar2_prior={}, const bool compute_XiDM=true);

      };
    }
  }
}

#endif
