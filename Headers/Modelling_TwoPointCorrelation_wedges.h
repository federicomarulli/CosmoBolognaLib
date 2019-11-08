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

	/// number of measured wedges in the input dataset
	int m_nWedges;

	/// the \f$\mu\f$ integral limits used to measure the wedges
	std::vector<std::vector<double>> m_mu_integral_limits;
	
	/// vector of booleans indicating the wedges to be modelled (m_use_wedge[i]=true -> the i-th wedge is modelled)
        std::vector<bool> m_use_wedge;

	/// vector containing the ordering of the data vector, which spcifies which wedge correponds to each data vector element
	std::vector<int> m_wedges_order;

	/// bolean to check if the model has been set
	bool m_ModelIsSet;
	
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
	 *  @param twop the wedges to model
	 *
	 *  @param nWedges the number of wedges
	 *	 
	 *  @param mu_integral_limits the \f$\mu\f$ integral limits
	 *  used to measure the wedges
	 *
	 *  @return object of type
	 *  Modelling_TwoPointCorrelation_wedges
	 */
	Modelling_TwoPointCorrelation_wedges (const std::shared_ptr<cbl::measure::twopt::TwoPointCorrelation> twop, const int nWedges=2, const std::vector<std::vector<double>> mu_integral_limits={{0., 0.5}, {0.5, 1.}});

	/**
	 *  @brief constructor
	 *  
	 *  @param twop_dataset the dataset containing the wedges to
	 *  model
	 *  
	 *  @param nWedges the number of wedges
	 *	 
	 *  @param mu_integral_limits the \f$\mu\f$ integral limits
	 *  used to measure the wedges
	 *
	 *  @return object of type
	 *  Modelling_TwoPointCorrelation_wedges
	 */
	Modelling_TwoPointCorrelation_wedges (const std::shared_ptr<data::Data> twop_dataset, const int nWedges, const std::vector<std::vector<double>> mu_integral_limits={{0., 0.5}, {0.5, 1.}});
      
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
	 *  @param nWedges the number of wedges
	 *
	 *  @return none
	 */
	void set_fit_range (const double xmin, const double xmax, const int nWedges=-1);

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
	 *  @brief set the de-wiggled model to fit the full shape of
	 *  the wedges of the two-point correlation function
	 *
	 *  the wedges of the two-point correlation function are
	 *  computed by cbl::modelling::twopt::xiWedges, where the
	 *  power spectrum is computed with the de-wiggled model by
	 *  cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fiducial_PkDM
	 *  with cbl::modelling::twopt::Pkmu_DeWiggled
	 *
	 *  the model has 7 parameters: 
	 *    - \f$\alpha_{\perp}\f$
	 *    - \f$\alpha_{\parallel}\f$
	 *    - \f$\Sigma_{NL,\perp}\f$ 
	 *    - \f$\Sigma_{NL,\parallel}\f$ 
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
	 *  @param SigmaNL_perpendicular_prior prior for the parameter
	 *  \f$\Sigma_{NL, \perp}\f$
	 *
	 *  @param SigmaNL_parallel_prior prior for the parameter 
	 *  \f$\Sigma_{NL, \parallel}\f$
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
	void set_model_fullShape_DeWiggled (const statistics::PriorDistribution alpha_perpendicular_prior={}, const statistics::PriorDistribution alpha_parallel_prior={}, const statistics::PriorDistribution SigmaNL_perpendicular_prior={}, const statistics::PriorDistribution SigmaNL_parallel_prior={}, statistics::PriorDistribution fsigma8_prior={}, statistics::PriorDistribution bsigma8_prior={}, const statistics::PriorDistribution SigmaS_prior={}, const bool compute_PkDM=true);

	/**
	 *  @brief set the Dispersion model to fit the wedges of the
	 *  two-point correlation function
	 *
	 *  the wedges of the two-point correlation function are
	 *  computed by cbl::modelling::twopt::xiWedges, where the
	 *  Dispersion power spectrum model \f$P(k, \mu)\f$ (see
	 *  e.g. Pezzotta et al. 2017,
	 *  https://arxiv.org/abs/1612.05645) is computed by
	 *  cbl:modelling::twopt::Pkmu_Dispersion
	 *
	 *  the model has 3 parameters:
	 *    - \f$f(z)\sigma_8(z)\f$
	 *    - \f$b(z)\sigma_8(z)\f$
	 *    - \f$\sigma_{12}\f$
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters
	 *
	 *  @author J.E. Garcia-Farieta
	 *  @author joegarciafa@unal.edu.co
	 *
	 *  @param fsigma8_prior prior for the parameter
	 *  \f$f(z)\sigma_8(z)\f$
	 *
	 *  @param bsigma8_prior prior for the parameter
	 *  \f$b(z)\sigma_8(z)\f$
	 *
	 *  @param sigma12_prior prior for the parameters
	 *  \f$\sigma_{12}\f$
	 *
	 *  @param DFoG true \f$\rightarrow\f$ Gaussian damping, false
	 *  \f$\rightarrow\f$ Lorentzian damping
	 *
	 *  @param compute_PkDM true \f$\rightarrow\f$ compute the
	 *  fiducial model of the dark matter power spectrum
	 *
	 *  @return none
	 */
	void set_model_Dispersion (statistics::PriorDistribution fsigma8_prior={}, statistics::PriorDistribution bsigma8_prior={}, const statistics::PriorDistribution sigma12_prior={}, const bool DFoG=true, const bool compute_PkDM=true); // [Jorge]

	/**
	 *  @brief set the Scoccimarro model to fit the wedges of the
	 *  two-point correlation function
	 *
	 *  the wedges of the two-point correlation function are
	 *  computed by cbl::modelling::twopt::xiWedges, where the
	 *  Scoccimarro power spectrum model \f$P(k, \mu)\f$ (see
	 *  e.g. Pezzotta et al. 2017,
	 *  https://arxiv.org/abs/1612.05645) is computed by
	 *  cbl:modelling::twopt::Pkmu_Scoccimarro_fitPezzotta
	 *
	 *  the model has 5 parameters:
	 *    - \f$f(z)\sigma_8(z)\f$
	 *    - \f$b(z)\sigma_8(z)\f$
	 *    - \f$\sigma_{12}\f$
	 *    - \f$k_\delta\f$
	 *    - \f$k_\theta\f$
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters
	 *
	 *  @author J.E. Garcia-Farieta
	 *  @author joegarciafa@unal.edu.co
	 *
	 *  @param fsigma8_prior prior for the parameter
	 *  \f$f(z)\sigma_8(z)\f$
	 *
	 *  @param bsigma8_prior prior for the parameter
	 *  \f$b(z)\sigma_8(z)\f$
	 *
	 *  @param sigma12_prior prior for the parameters
	 *  \f$\sigma_{12}\f$
	 *
	 *  @param kd_prior prior for the parameters \f$ k_\delta\f$
	 *
	 *  @param kt_prior prior for the parameters \f$k_t\f$
	 *
	 *  @param DFoG true \f$\rightarrow\f$ Gaussian damping, false
	 *  \f$\rightarrow\f$ Lorentzian damping
	 *
	 *  @param compute_PkDM true \f$\rightarrow\f$ compute the
	 *  fiducial model of the dark matter power spectrum
	 *
	 *  @return none
	 */
	void set_model_Scoccimarro_fitPezzotta (const statistics::PriorDistribution fsigma8_prior={}, const statistics::PriorDistribution bsigma8_prior={}, const statistics::PriorDistribution sigma12_prior={}, const statistics::PriorDistribution kd_prior={}, const statistics::PriorDistribution kt_prior={}, const bool DFoG=true, const bool compute_PkDM=true); 

	/**
	 *  @brief set the Scoccimarro model to fit the wedges of the
	 *  two-point correlation function
	 *
	 *  the wedges of the two-point correlation function are
	 *  computed by cbl::modelling::twopt::xiWedges, where the
	 *  Scoccimarro power spectrum model \f$P(k, \mu)\f$ (see
	 *  e.g. Bel et al. 2019, https://arxiv.org/abs/1809.09338) is
	 *  computed by cbl:modelling::twopt::Pkmu_Scoccimarro_fitBel
	 *
	 *  the model has 8 parameters:
	 *    - \f$f(z)\sigma_8(z)\f$
	 *    - \f$b(z)\sigma_8(z)\f$
	 *    - \f$\sigma_{12}\f$
	 *    - \f$k_\delta\f$
	 *    - \f$bb\f$
	 *    - \f$a1\f$
	 *    - \f$a2\f$
	 *    - \f$a3\f$
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters
	 *
	 *  @author J.E. Garcia-Farieta
	 *  @author joegarciafa@unal.edu.co
	 *
	 *  @param fsigma8_prior prior for the parameter
	 *  \f$f(z)\sigma_8(z)\f$
	 *
	 *  @param bsigma8_prior prior for the parameter
	 *  \f$b(z)\sigma_8(z)\f$
	 *
	 *  @param sigma12_prior prior for the parameters
	 *  \f$\sigma_{12}\f$
	 *
	 *  @param kd_prior prior for the parameters \f$ k_d\f$
	 *
	 *  @param bb_prior prior for the parameters \f$bb\f$
	 *
	 *  @param a1_prior prior for the parameters \f$a1\f$
	 *
	 *  @param a2_prior prior for the parameters \f$a2\f$
	 *
	 *  @param a3_prior prior for the parameters \f$a3\f$
	 *
	 *  @param DFoG true \f$\rightarrow\f$ Gaussian damping, false
	 *  \f$\rightarrow\f$ Lorentzian damping
	 *
	 *  @param compute_PkDM true \f$\rightarrow\f$ compute the
	 *  fiducial model of the dark matter power spectrum
	 *
	 *  @return none
	 */
	void set_model_Scoccimarro_fitBel (const statistics::PriorDistribution fsigma8_prior={}, const statistics::PriorDistribution bsigma8_prior={}, const statistics::PriorDistribution sigma12_prior={}, const statistics::PriorDistribution kd_prior={}, const statistics::PriorDistribution bb_prior={}, const statistics::PriorDistribution a1_prior={}, const statistics::PriorDistribution a2_prior={}, const statistics::PriorDistribution a3_prior={}, const bool DFoG=true, const bool compute_PkDM=true);

	/**
	 *  @brief set the Scoccimarro model to fit the wedges of the
	 *  two-point correlation function
	 *
	 *  the wedges of the two-point correlation function are
	 *  computed by cbl::modelling::twopt::xiWedges, where the
	 *  Scoccimarro power spectrum model \f$P(k, \mu)\f$ (see
	 *  e.g. Scoccimarro et al. 2004,
	 *  https://arxiv.org/abs/astro-ph/0407214) is computed by
	 *  cbl:modelling::twopt::Pkmu_Scoccimarro
	 *
	 *  the model has 3 parameters:
	 *    - \f$f(z)\sigma_8(z)\f$
	 *    - \f$b(z)\sigma_8(z)\f$
	 *    - \f$\sigma_{12}\f$
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters
	 *
	 *  @author J.E. Garcia-Farieta
	 *  @author joegarciafa@unal.edu.co
	 *
	 *  @param fsigma8_prior prior for the parameter
	 *  \f$f(z)\sigma_8(z)\f$
	 *
	 *  @param bsigma8_prior prior for the parameter
	 *  \f$b(z)\sigma_8(z)\f$
	 *
	 *  @param sigma12_prior prior for the parameters
	 *  \f$\sigma_{12}\f$
	 *
	 *  @param DFoG true \f$\rightarrow\f$ Gaussian damping, false
	 *  \f$\rightarrow\f$ Lorentzian damping
	 *
	 *  @param compute_PkDM true \f$\rightarrow\f$ compute the
	 *  fiducial model of the dark matter power spectrum
	 *
	 *  @return none
	 */
	void set_model_Scoccimarro (const statistics::PriorDistribution fsigma8_prior={}, const statistics::PriorDistribution bsigma8_prior={}, const statistics::PriorDistribution sigma12_prior={}, const bool DFoG=true, const bool compute_PkDM=true);

	/**
	 * @brief set the TNS (Taruya, Nishimichi and Saito) model to
	 * fit the wedges of the * two-point correlation function
	 *
	 *  the wedges of the two-point correlation function are
	 *  computed by cbl::modelling::twopt::xiWedges, where the TNS
	 *  power spectrum model \f$P(k, \mu)\f$ (see e.g. Taruya et
	 *  al. 2010, https://arxiv.org/abs/1006.0699) is computed by
	 *  cbl:modelling::twopt::Pkmu_TNS
	 *
	 *  the model has 3+N parameters:
	 *    - \f$f(z)\sigma_8(z)\f$
	 *    - \f$b(z)\sigma_8(z)\f$
	 *    - \f$\sigma_{12}\f$
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters
	 *
	 *  @author J.E. Garcia-Farieta
	 *  @author joegarciafa@unal.edu.co
	 *
	 *  @param fsigma8_prior prior for the parameter
	 *  \f$f(z)\sigma_8(z)\f$
	 *
	 *  @param bsigma8_prior prior for the parameter
	 *  \f$b(z)\sigma_8(z)\f$
	 *
	 *  @param sigma12_prior prior for the parameters
	 *  \f$\sigma_{12}\f$
	 *
	 *  @param DFoG true \f$\rightarrow\f$ Gaussian damping, false
	 *  \f$\rightarrow\f$ Lorentzian damping
	 *
	 *  @param compute_PkDM true \f$\rightarrow\f$ compute the
	 *  fiducial model of the dark matter power spectrum
	 *
	 *  @return none
	 */
	void set_model_TNS (const statistics::PriorDistribution fsigma8_prior={}, const statistics::PriorDistribution bsigma8_prior={}, const statistics::PriorDistribution sigma12_prior={}, const bool DFoG=true, const bool compute_PkDM=true);

	/**
	 *  @brief set the mode-coupling model to fit the full shape
	 *  of the wedges of the two-point correlation function
	 *
	 *  the wedges of the two-point correlation function are
	 *  computed by cbl::modelling::twopt::xiWedges, where the
	 *  power spectrum is computed with the mode-coupling model by
	 *  cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fiducial_PkDM
	 *  with cbl::modelling::twopt::Pkmu_ModeCoupling
	 *
	 *  the model has 6 parameters: 
	 *    - \f$\alpha_{\perp}\f$
	 *    - \f$\alpha_{\parallel}\f$
	 *    - \f$f(z)\sigma_8(z)\f$
	 *    - \f$b(z)\sigma_8(z)\f$
	 *    - \f$\sigma_v\f$
	 *    - \f$A_{MC}\f$
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
	 *  @param SigmaV_prior prior for the parameters
	 *  \f$\sigma_v\f$
	 *
	 *  @param AMC_prior prior for the parameters
	 *  \f$A_{MC}\f$
	 *
	 *  @param compute_PkDM true \f$\rightarrow\f$ compute the
	 *  fiducial model of the dark matter power spectrum
	 *
	 *  @return none
	 */
	void set_model_fullShape_ModeCoupling (const statistics::PriorDistribution alpha_perpendicular_prior={}, const statistics::PriorDistribution alpha_parallel_prior={}, const statistics::PriorDistribution fsigma8_prior={}, const statistics::PriorDistribution bsigma8_prior={}, const statistics::PriorDistribution SigmaV_prior={}, const statistics::PriorDistribution AMC_prior={}, const bool compute_PkDM=true);

	/**
	 *  @brief set the model to fit the wedges of the two-point
	 *  correlation function, used for anisotropic BAO
	 *  measurements
	 *
	 *  the wedges of the two-point correlation function are
	 *  computed by cbl::modelling::twopt::xiWedges_BAO
	 *
	 *  @param alpha_perpendicular_prior prior for the parameter
	 *  \f$\alpha_{\perp}\f$
	 *
	 *  @param alpha_parallel_prior prior for the parameter 
	 *  \f$\alpha_{\parallel}\f$
	 *
	 *  @param Bperp_prior prior for the parameter
	 *  \f$B_\perp\f$
	 *
	 *  @param Bpar_prior prior for the parameter
	 *  \f$B_\parallel\f$
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
	 *  @param compute_XiDM true \f$\rightarrow\f$ compute the
	 *  fiducial model of the dark matter two-point correlation 
	 *  function
	 *
	 *  @return none
	 */
	void set_model_BAO (const statistics::PriorDistribution alpha_perpendicular_prior={}, const statistics::PriorDistribution alpha_parallel_prior={}, const statistics::PriorDistribution Bperp_prior={}, const statistics::PriorDistribution Bpar_prior={}, const statistics::PriorDistribution Aperp0_prior={}, const statistics::PriorDistribution Apar0_prior={}, const statistics::PriorDistribution Aperp1_prior={}, const statistics::PriorDistribution Apar1_prior={}, const statistics::PriorDistribution Aperp2_prior={}, const statistics::PriorDistribution Apar2_prior={}, const bool compute_XiDM=true);

        /**
         *  @brief write the model at xx for given parameters
         *
         *  @param output_dir the output directory
	 *
         *  @param output_file the output file
	 *
         *  @param xx vector of points at which the model is computed
	 *
         *  @param parameters vector containing the input parameters
         *  used to compute the model; if this vector is not provided,
         *  the model will be computed using the best-fit parameters
         *
         *  @return none
         */
        void write_model (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> parameters);

        /**
         *  @brief write the model at xx 
         *  with best-fit parameters obtained from likelihood
         *  maximization
         *
         *  @param output_dir the output directory
         *  @param output_file the output file
         *  @param xx vector of points at which the model is computed,
         *
         *  @return none
         */
        void write_model_at_bestfit (const std::string output_dir, const std::string output_file, const std::vector<double> xx);

        /**
         *  @brief write the model at xx 
         *  computing 16th, 50th and 84th percentiles
         *  from the chains.
         *
         *  @param output_dir the output directory
         *  @param output_file the output file
         *  @param xx vector of points at which the model is computed,
         *  @param start the starting position for each chain
         *  @param thin the position step
         *
         *  @return none
         */
        void write_model_from_chains (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const int start=0, const int thin=1);

      };
    }
  }
}

#endif
