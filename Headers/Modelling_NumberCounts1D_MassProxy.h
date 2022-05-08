/********************************************************************
 *  Copyright (C) 2021 by Giorgio Lesci and Federico Marulli        *
 *  giorgio.lesci2@unibo.it, federico.marulli3@unibo.it             *
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
 *  @file Headers/Modelling_NumberCounts1D_MassProxy.h
 *
 *  @brief The class Modelling_NumberCounts1D_MassProxy
 *
 *  This file defines the interface of the class Modelling_NumberCounts1D_MassProxy, 
 *  used to  model number counts as a function of the mass proxy
 *
 *  @author Giorgio Lesci, Federico Marulli
 *
 *  @author giorgio.lesci2@unibo.it, federico.marulli3@unibo.it
 */

#ifndef __MODELLINGNCMP__
#define __MODELLINGNCMP__


#include "Modelling_NumberCounts1D.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    /**
     *  @brief The namespace of the <B> number counts
     *  modelling </B>
     *  
     *  The \e modelling::numbercounts namespace contains all the functions
     *  and classes to model number counts
     */
    namespace numbercounts {
    
      /**
       *  @class Modelling_NumberCounts1D_MassProxy
       *  Modelling_NumberCounts1D_MassProxy.h
       *  "Headers/Modelling_NumberCounts1D_MassProxy.h"
       *
       *  @brief The class Modelling_NumberCounts1D_MassProxy
       *
       *  This file defines the interface of the base class
       *  Modelling_NumberCounts1D_MassProxy, used for modelling 
       *  number counts measurements as a function of a mass proxy
       *
       */
      class Modelling_NumberCounts1D_MassProxy : public Modelling_NumberCounts1D
      {

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  _NumberCounts1D_MassProxy
	 */
	Modelling_NumberCounts1D_MassProxy () = default;
	
	/**
	 *  @brief constuctor
	 *  @param nc the number counts to model
	 */
	Modelling_NumberCounts1D_MassProxy (const std::shared_ptr<cbl::measure::numbercounts::NumberCounts> nc) 
	  : Modelling_NumberCounts1D (nc) {}
	
	/**
	 *  @brief constuctor
	 *  @param dataset the number counts dataset
	 *  @param hist_type the histogram type
	 *  @param fact the normalization factor
	 */
	Modelling_NumberCounts1D_MassProxy (const std::shared_ptr<cbl::data::Data> dataset, glob::HistogramType hist_type, double fact)
	  : Modelling_NumberCounts1D (dataset, hist_type, fact) {}
	
	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~Modelling_NumberCounts1D_MassProxy () = default;

	///@}


	/**
	 *  @name Member functions used to set the model parameters
	 *
	 */
	///@{

	/**
	 *  @brief Set the cosmological parameters used to model the 
	 *  number counts as a function of a mass proxy, here written
	 *  as \f$\lambda\f$, with the model expressed as follows:
	 *
	 *  \f$ \langle N(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle 
	 *  = w(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\,\,\Omega 
	 *  \int_{0}^{\infty} {\rm d} z_{\rm tr}\,\,
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{0}^{\infty}{\rm d}\lambda_{\rm tr}\,\,
	 *  P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\,
	 *  \int_{\Delta z_{\text{ob},j}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  \int_{\Delta\lambda_{\text{ob},i}}{\rm d} \lambda_{\rm ob} 
	 *  \,\,P(\lambda_{\rm ob}|\lambda_{\rm tr}), \f$
	 *
	 *  where \f$ w(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j}) \f$ is
	 *  the weight derived from the selection function (see e.g. Lesci et al. 2021),
	 *  and \f$\Omega\f$ is the survey effective area.
	 *
	 *  The distirbution \f$P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\f$ is expressed as:
	 *
	 *  \f$P(\lambda_{\rm tr}|M_{\rm tr},z_{\rm tr})= 
	 *  P(M_{\rm tr}|\lambda_{\rm tr},z_{\rm tr})\,
	 *  P(\lambda_{\rm tr}|z_{\rm tr})\,/\,P( M_{\rm tr}|z_{\rm tr}),\f$
	 *
	 *  where \f$P(M_{\rm tr}|\lambda_{\rm tr},z_{\rm tr})\f$ is
	 *  a log-normal whose mean is given by the mass-mass proxy
	 *  relation, i.e. 
	 *
	 *  \f$\log M = \alpha + \beta 
	 *  \log (\lambda/\lambda_{\rm piv}) + \gamma \log (f(z)),\f$
	 *
	 *  and whose standard deviation is given by the intrinsic scatter, \f$ \sigma_{\rm intr} \f$, expressed as:
	 *
	 *  \f$ \sigma_{\rm intr} = \sigma_0 + \sigma_{\lambda} 
	 *  \log (\lambda/\lambda_{\rm piv})^{e_{\lambda}} + \sigma_z \log (f(z))^{e_z}.\f$
	 *
	 *  The distribution \f$ P(\lambda_{\rm tr}|z_{\rm tr}) \f$ is derived
	 *  from mock catalogues, and it has the functional form of a power-law
	 *  with an exponential cut-off, i.e.:
	 *  \f$ P(\lambda_{\rm tr}|z_{\rm tr}) = a \, \lambda_{\rm tr}^{-b} \, e^{-c\lambda_{\rm tr}} \f$.
	 *
	 *  Furthermore, \f$P( M_{\rm tr}|z_{\rm tr}) = \int_{0}^{\infty}
	 *  {\rm d} \lambda_{\rm tr}\, 
	 *  P(M_{\rm tr}|\lambda_{\rm tr},z_{\rm tr})\,P(\lambda_{\rm tr}|z_{\rm tr}).\f$
	 *
	 *  Finally, \f$P(z_{\rm ob}|z_{\rm tr})\f$ and \f$P(\lambda_{\rm ob}|\lambda_{\rm tr})\f$
	 *  are Gaussian distributions.
	 *
	 *
	 *  @param scalrel_z_evo functional form of the redshift evolution
	 *  function in the scaling relation: "E_z" \f$\rightarrow\f$ 
	 *  \f$ f(z)=E(z)/E(z_{piv}) \f$, "direct" \f$\rightarrow\f$ 
	 *  \f$ f(z)=(1+z)/(1+z_{piv}) \f$
	 *
	 *  @param z_error_type if "absolute", set the absolute error
	 *  on redshift in the model; if "relative", set the relative error
	 *
	 *  @param proxy_error_type if "absolute", set the absolute error
	 *  on the mass proxy in the model; if "relative", set the relative error
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_param_prior vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @param alpha_prior prior on the normalization 
	 *  of the mass-mass proxy scaling relation, \f$\alpha\f$
	 *
	 *  @param beta_prior prior on the slope 
	 *  of the mass-mass proxy scaling relation, \f$\beta\f$
	 *
	 *  @param gamma_prior prior on the redshift evolution factor 
	 *  of the mass-mass proxy scaling relation, \f$\gamma\f$
	 *
	 *  @param scatter0_prior prior on the 
	 *  constant term of the intrinsic scatter, \f$ \sigma_0 \f$
	 *
	 *  @param scatterM_prior prior on the factor in the
	 *  proxy-dependent term of the intrinsic scatter, \f$ \sigma_{\lambda} \f$
	 *
	 *  @param scatterM_exponent_prior prior on the exponent in the
	 *  proxy-dependent term of the intrinsic scatter, \f$ e_{\lambda} \f$
	 *
	 *  @param scatterz_prior prior on the factor in the
	 *  redshift-dependent term of the intrinsic scatter, \f$ \sigma_z \f$
	 *
	 *  @param scatterz_exponent_prior prior on the exponent in the
	 *  redshift-dependent term of the intrinsic scatter, \f$ e_z \f$
	 *
	 *  @param z_bias_prior prior on the redshif bias,
	 *  governing the position of the mean of the Gaussian \f$ P(z_{\rm ob}|z_{\rm tr}) \f$.
	 *  For example, let's say that the redshift bias in the data is \f$\Delta z=0.02\f$:
	 *  the data redshifts would be corrected as \f$z_{\rm corr} = z - 0.02 (1+z)\f$, 
	 *  where \f$z\f$ is the not corrected redshift. If you correct the data in this way,
	 *  set this prior as constant and equal to zero. Otherwise you can choose to not
	 *  correct the data and set a non-zero prior on the redshift bias: in this case,
	 *  the mean of \f$ P(z_{\rm ob}|z_{\rm tr}) \f$ will be given by 
	 *  \f$z_{\rm corr} = z + \Delta z (1+z)\f$.
	 *
	 *  @param proxy_bias_prior prior on the mass proxy bias, analogous to z_bias_prior
	 *
	 *  @param z_error_prior prior on the redshift error, governing the distribution
	 *  \f$ P(z_{\rm ob}|z_{\rm tr}) \f$. Depending on the input parameter
	 *  z_error_type, this prior is related to the absolute error or to the relative error
	 *
	 *  @param proxy_error_prior prior on the mass proxy error, governing the distribution
	 *  \f$ P(\lambda_{\rm ob}|\lambda_{\rm tr}) \f$. Depending on the input parameter
	 *  proxy_error_type, this prior is related to the absolute error or to the relative error
	 *
	 *  @param Plambda_prior priors on the three parameters,
	 *  respectively \f$ a \f$, \f$ b \f$, \f$ c \f$, defining
	 *  the distribution \f$ P(\lambda_{\rm tr}|z_{\rm tr}) \f$
	 *  
	 */
	void set_model_NumberCounts_cosmology (const std::string scalrel_z_evo, const std::string z_error_type, const std::string proxy_error_type, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_param_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const statistics::PriorDistribution z_bias_prior, const statistics::PriorDistribution proxy_bias_prior, const statistics::PriorDistribution z_error_prior, const statistics::PriorDistribution proxy_error_prior, const std::vector<statistics::PriorDistribution> Plambda_prior);
	
	/**
	 *  @brief Set the cosmological parameters used to model the 
	 *  number counts as a function of a mass proxy, here written
	 *  as \f$\lambda\f$, with the model expressed as follows:
	 *
	 *  \f$ \langle N(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle 
	 *  = w(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\,\,\Omega 
	 *  \int_{0}^{\infty} {\rm d} z_{\rm tr}\,\,
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{0}^{\infty}{\rm d}\lambda_{\rm tr}\,\,
	 *  P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\,
	 *  \int_{\Delta z_{\text{ob},j}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  \int_{\Delta\lambda_{\text{ob},i}}{\rm d} \lambda_{\rm ob} 
	 *  \,\,P(\lambda_{\rm ob}|\lambda_{\rm tr}), \f$
	 *
	 *  where \f$ w(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j}) \f$ is
	 *  the weight derived from the selection function (see e.g. Lesci et al. 2021),
	 *  and \f$\Omega\f$ is the survey effective area.
	 *
	 *  The distirbution \f$P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\f$ is 
	 *  a log-normal whose mean is given by the mass-mass proxy
	 *  relation, i.e. 
	 *
	 *  \f$\log (\lambda/\lambda_{\rm piv}) = \alpha + \beta 
	 *  \log (M/M_{\rm piv}) + \gamma \log (f(z)),\f$
	 *
	 *  and whose standard deviation is given by the intrinsic scatter, \f$ \sigma_{\rm intr} \f$, expressed as:
	 *
	 *  \f$ \sigma_{\rm intr} = \sigma_0 + \sigma_{M} 
	 *  \log (M/M_{\rm piv})^{e_{M}} + \sigma_z \log (f(z))^{e_z}. \f$
	 *
	 *  Finally, \f$P(z_{\rm ob}|z_{\rm tr})\f$ and \f$P(\lambda_{\rm ob}|\lambda_{\rm tr})\f$
	 *  are Gaussian distributions.
	 *
	 *
	 *  @param scalrel_z_evo functional form of the redshift evolution
	 *  function in the scaling relation: "E_z" \f$\rightarrow\f$ 
	 *  \f$ f(z)=E(z)/E(z_{piv}) \f$, "direct" \f$\rightarrow\f$ 
	 *  \f$ f(z)=(1+z)/(1+z_{piv}) \f$
	 *
	 *  @param z_error_type if "absolute", set the absolute error
	 *  on redshift in the model; if "relative", set the relative error
	 *
	 *  @param proxy_error_type if "absolute", set the absolute error
	 *  on the mass proxy in the model; if "relative", set the relative error
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_param_prior vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @param alpha_prior prior on the normalization 
	 *  of the mass-mass proxy scaling relation, \f$\alpha\f$
	 *
	 *  @param beta_prior prior on the slope 
	 *  of the mass-mass proxy scaling relation, \f$\beta\f$
	 *
	 *  @param gamma_prior prior on the redshift evolution factor 
	 *  of the mass-mass proxy scaling relation, \f$\gamma\f$
	 *
	 *  @param scatter0_prior prior on the 
	 *  constant term of the intrinsic scatter, \f$ \sigma_0 \f$
	 *
	 *  @param scatterM_prior prior on the factor in the
	 *  proxy-dependent term of the intrinsic scatter, \f$ \sigma_{\lambda} \f$
	 *
	 *  @param scatterM_exponent_prior prior on the exponent in the
	 *  proxy-dependent term of the intrinsic scatter, \f$ e_{\lambda} \f$
	 *
	 *  @param scatterz_prior prior on the factor in the
	 *  redshift-dependent term of the intrinsic scatter, \f$ \sigma_z \f$
	 *
	 *  @param scatterz_exponent_prior prior on the exponent in the
	 *  redshift-dependent term of the intrinsic scatter, \f$ e_z \f$
	 *
	 *  @param z_bias_prior prior on the redshif bias,
	 *  governing the position of the mean of the Gaussian \f$ P(z_{\rm ob}|z_{\rm tr}) \f$.
	 *  For example, let's say that the redshift bias in the data is \f$\Delta z=0.02\f$:
	 *  the data redshifts would be corrected as \f$z_{\rm corr} = z - 0.02 (1+z)\f$, 
	 *  where \f$z\f$ is the not corrected redshift. If you correct the data in this way,
	 *  set this prior as constant and equal to zero. Otherwise you can choose to not
	 *  correct the data and set a non-zero prior on the redshift bias: in this case,
	 *  the mean of \f$ P(z_{\rm ob}|z_{\rm tr}) \f$ will be given by 
	 *  \f$z_{\rm corr} = z + \Delta z (1+z)\f$.
	 *
	 *  @param proxy_bias_prior prior on the mass proxy bias, analogous to z_bias_prior
	 *
	 *  @param z_error_prior prior on the redshift error, governing the distribution
	 *  \f$ P(z_{\rm ob}|z_{\rm tr}) \f$. Depending on the input parameter
	 *  z_error_type, this prior is related to the absolute error or to the relative error
	 *
	 *  @param proxy_error_prior prior on the mass proxy error, governing the distribution
	 *  \f$ P(\lambda_{\rm ob}|\lambda_{\rm tr}) \f$. Depending on the input parameter
	 *  proxy_error_type, this prior is related to the absolute error or to the relative error
	 *  
	 */
	void set_model_NumberCounts_cosmology (const std::string scalrel_z_evo, const std::string z_error_type, const std::string proxy_error_type, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_param_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const statistics::PriorDistribution z_bias_prior, const statistics::PriorDistribution proxy_bias_prior, const statistics::PriorDistribution z_error_prior, const statistics::PriorDistribution proxy_error_prior);

	///@}

      };
    }
  }
}

#endif
