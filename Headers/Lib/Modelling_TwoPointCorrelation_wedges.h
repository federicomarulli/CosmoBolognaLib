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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation_wedges.h
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


namespace cosmobl {
  
  namespace modelling {

    namespace twopt {
    
      /**
       *  @class Modelling_TwoPointCorrelation_wedges
       *  Modelling_TwoPointCorrelation_wedges.h
       *  "Headers/Lib/Modelling_TwoPointCorrelation_wedges.h"
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
	vector<int> m_wedges_order;

	/// the wedges aperture
	double m_deltamu;

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

	/// \f$B_{\perp}\f$ nuisance parameter for anisotropic analysis 
	shared_ptr<statistics::Parameter> m_Bperp;

	/// \f$B_{\parallel}\f$ nuisance parameter for anisotropic analysis 
	shared_ptr<statistics::Parameter> m_Bpar;
	
	/// \f$A_0^{\perp}\f$ nuisance parameter for anisotropic analysis 
	shared_ptr<statistics::Parameter> m_Aperp0;
	
	/// \f$A_1^{\perp}\f$ nuisance parameter for anisotropic analysis 
	shared_ptr<statistics::Parameter> m_Aperp1;
	
	/// \f$A_2^{\perp}\f$ nuisance parameter for anisotropic analysis 
	shared_ptr<statistics::Parameter> m_Aperp2;
	
	/// \f$A_0^{\parallel}\f$ nuisance parameter for anisotropic analysis 
	shared_ptr<statistics::Parameter> m_Apar0;
	
	/// \f$A_1^{\parallel}\f$ nuisance parameter for anisotropic analysis 
	shared_ptr<statistics::Parameter> m_Apar1;

	/// \f$A_2^{\parallel}\f$ nuisance parameter for anisotropic analysis 
	shared_ptr<statistics::Parameter> m_Apar2;

	///@}
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
	Modelling_TwoPointCorrelation_wedges (const shared_ptr<cosmobl::measure::twopt::TwoPointCorrelation> twop);

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
	Modelling_TwoPointCorrelation_wedges (const shared_ptr<data::Data> twop_dataset, const int nwedges);
      
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
	 *  @param vector containing the fitting range for
	 *  the wedges
	 *
	 *  @return none
	 */
	void set_fit_range (vector<vector<double>> fit_range);

	/**
	 *  @brief return the protected member 
	 *  alpha_perpendicular
	 *
	 *  @return pointer to parameter alpha_perpendicular
	 */
	shared_ptr<statistics::Parameter> alpha_perpendicular () { return m_alpha_perpendicular;}

	/**
	 *  @brief return the protected member 
	 *  m_alpha_parallel
	 *
	 *  @return pointer to parameter alpha_parallel
	 */
	shared_ptr<statistics::Parameter> alpha_parallel () { return m_alpha_parallel;}

	/**
	 *  @brief return the protected member 
	 *  m_bsigma8
	 *
	 *  @return pointer to parameter bsigma8
	 */
	shared_ptr<statistics::Parameter> bsigma8 () { return m_bsigma8;}

	/**
	 *  @brief return the protected member 
	 *  m_fsigma8
	 *
	 *  @return pointer to parameter fsigma8
	 */
	shared_ptr<statistics::Parameter> fsigma8 () { return m_fsigma8;}

	/**
	 *  @brief return the protected member 
	 *  m_SigmaS
	 *
	 *  @return pointer to parameter SigmaS
	 */
	shared_ptr<statistics::Parameter> SigmaS () { return m_SigmaS;}

	/**
	 *  @brief return the protected member 
	 *  m_Bperp
	 *
	 *  @return pointer to parameter Bperp
	 */
	shared_ptr<statistics::Parameter> Bperp () { return m_Bperp;}

	/**
	 *  @brief return the protected member 
	 *  m_B0
	 *
	 *  @return pointer to parameter B0
	 */
	shared_ptr<statistics::Parameter> Bpar () { return m_Bpar;}
	
	/**
	 *  @brief return the protected member 
	 *  m_Aperp0
	 *
	 *  @return pointer to parameter Aperp0
	 */
	shared_ptr<statistics::Parameter> Aperp0 () { return m_Aperp0;}
	
	/**
	 *  @brief return the protected member 
	 *  m_Aperp0
	 *
	 *  @return pointer to parameter Aperp1
	 */
	shared_ptr<statistics::Parameter> Aperp1 () { return m_Aperp1;}
	
	/**
	 *  @brief return the protected member 
	 *  m_Aperp2
	 *
	 *  @return pointer to parameter Aperp2
	 */
	shared_ptr<statistics::Parameter> Aperp2 () { return m_Aperp2;}
	
	/**
	 *  @brief return the protected member 
	 *  m_Apar0
	 *
	 *  @return pointer to parameter Apar0
	 */
	shared_ptr<statistics::Parameter> Apar0 () { return m_Apar0;}
	
	/**
	 *  @brief return the protected member 
	 *  m_Apar1
	 *
	 *  @return pointer to parameter Apar1
	 */
	shared_ptr<statistics::Parameter> Apar1 () { return m_Apar1;}

	/**
	 *  @brief return the protected member 
	 *  m_Apar2
	 *
	 *  @return pointer to parameter Apar2
	 */
	shared_ptr<statistics::Parameter> Apar2 () { return m_Apar2;}

	/**
	 *  @brief set the parameters used to model the full shape
	 *  of the multipoles of the two-point correlation function
	 *
	 *  the model is the following:
	 *
	 * The functions computes the multipoles of the 
	 * two-point correlation function
	 *
	 * \f[
	 * \xi_l(s) = \frac{i^l}{2\pi^2} \int \mathrm{d} k P_l(k) j_l(ks);
	 * \f]
	 *
	 * where \f$j_l(ks)\f$ are the bessel functions.
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
	void set_model_fullShape (const statistics::Prior alpha_perpendicular_prior={}, const statistics::Prior alpha_parallel_prior={}, statistics::Prior fsigma8_prior={}, statistics::Prior bsigma8_prior={}, const statistics::Prior SigmaS_prior={}, const bool compute_PkDM=true);

	/**
	 * @brief set the parameters used to model the two-point 
	 * correlation function, intended for anisotropic BAO
	 * measurements (Ross et al. 2017). In only works with
	 * two two-point correlation function wedges.
	 *
	 * The functions computes the wedges of the 
	 * two-point correlation function:
	 *
	 * \f[
	 * \xi_{\perp}(s) =  B_{perp}\xi_{\perp}(s, \alpha_{\perp}, \alpha_{\parallel})+A_{perp}^0+\frac{A_{perp}^1}{s}+\frac{A_{perp}^2}{s^2}; \\
	 * \xi_{\parallel}(s) =  B_{parallel}\xi_{\parallel}(s, \alpha_{\parallel}, \alpha_{\parallel})+A_{parallel}^0+\frac{A_{parallel}^1}{s}+\frac{A_{parallel}^2}{s^2}; \\
	 * \f]
	 *
	 * where \f$\xi_{\perp}\f$, \f$\xi_{\parallel}\f$ are the two wedges of the two-point correlation function.
	 *
	 * The function takes as inputs ten parameters
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
	 */
	void set_model_BAO (const statistics::Prior alpha_perpendicular_prior={}, const statistics::Prior alpha_parallel_prior={}, const statistics::Prior Bperp_prior={}, const statistics::Prior Bpar_prior={}, const statistics::Prior Aperp0_prior={}, const statistics::Prior Apar0_prior={}, const statistics::Prior Aperp1_prior={}, const statistics::Prior Apar1_prior={}, const statistics::Prior Aperp2_prior={}, const statistics::Prior Apar2_prior={}, const bool compute_XiDM=true);

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

	void write_model (const string dir, const string file, const vector<double> xx, const vector<double> parameter, const int start, const int thin);

      };
    }
  }
}

#endif
