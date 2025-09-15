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
 *  @file Headers/Modelling_PowerSpectrum_Angular.h
 *
 *  @brief The class Modelling_PowerSpectrum_angular
 *
 *  This file defines the interface of the class
 *  Modelling_PowerSpectrum_angular, used to model the angular power spectrum
 *
 *  @author Federico Marulli, Massimiliano Romanello
 *
 *  @author federico.marulli3@unibo.it, massimilia.romanell2@unibo.it
 */

#ifndef __MODELLINGPOWSPECTRUMANG__
#define __MODELLINGPOWSPECTRUMANG__

#include "Sigma.h"
#include "Modelling.h"
#include "PowerSpectrum_Angular.h"
#include "ModelFunction_PowerSpectrum_Angular.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    namespace angularpk {
      
      /**
       *  @class Modelling_PowerSpectrum_angular
       *  Modelling_PowerSpectrum_Angular.h
       *  "Headers/Modelling_PowerSpectrum_Angular.h"
       *
       *  @brief The class Modelling_PowerSpectrum_angular
       *
       *  This file defines the interface of the base class
       *  Modelling_PowerSpectrum_angular, used for modelling
       *  the angular power spectrum 
       */
      class Modelling_PowerSpectrum_angular : public Modelling {
      protected:
	
	/// the container of parameters for angular power spectrum model computation
	std::shared_ptr<modelling::angularpk::STR_data_model> m_data_model;
	
      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 */
	Modelling_PowerSpectrum_angular () = default;

	/**
	 *  @brief constructor
	 *  
	 *  @param Pow the angular power spectrum to model
	 */
	Modelling_PowerSpectrum_angular (const std::shared_ptr<cbl::measure::angularpk::PowerSpectrum_angular> Pow)
	{
	  m_data = Pow->dataset();
	};
	
	/**
	 *  @brief constructor
	 *  
	 *  @param dataset the dataset containing the angular power
	 *  spectrum to model
	 */
	Modelling_PowerSpectrum_angular (const std::shared_ptr<data::Data> dataset)
	{
	  m_data = dataset;
	};
	
	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~Modelling_PowerSpectrum_angular () = default;

	///@}
	

	/**
	 *  @name Member functions used to set the model parameters
	 */
	///@{

	/**
	 *  @brief set the model to fit the angular power spectrum
	 *
	 *  The model is the following:
	 *
	 *  \f[ C_l = \frac{2}{\pi}\int_{0}^{\infty} dk k^2 P_{mat}(k)\left(\psi_l(k)+\psi_l^r(k)\right)^2 \f]
	 *
	 *  where the kernel functions are:
	 *
	 *  \f[\psi_l(k) = \int_{0}^{\infty} dz \frac{1}{N}\frac{dN}{dz} b(z) f(z) D(z)j_l(kr(z)) \f]
	 *
	 *  \f[\psi_l^r(k) = \int_{0}^{\infty} dz
	 *  \frac{1}{N}\frac{dN}{dz} f(z)
	 *  D(z)\left[\frac{2l^2+2l-1}{(2l+3)(2l-1)}j_l(kr(z))-\frac{l(l-1)}{(2l+1)(2l-1)}j_{l-2}(kr(z))-\frac{(l+1)(l+2)}{(2l+1)(2l+3)}j_{l+2}(kr(z))\right]
	 *  \f] where the latter accounts for the contribution of
	 *  redshift-space distorsion (RSD), f(z) is the linear growth
	 *  rate, D(z) is the linear growth factor.  In the Limber
	 *  approximation this formula can be extended to l<2 and it
	 *  reduces to
	 *
	 *  \f[ C_l = \frac{b^2}{N^2}\int_{0}^{\infty}
	 *  \frac{dN}{dz}\frac{dN}{dz} P_{mat}
	 *  \left(\frac{l+1/2}{r(z)}\right)
	 *  \frac{H(z)}{c}\frac{1}{r^2(z)} dz\f]
	 *
	 *  the model has N cosmological parameters, plus 2 parameter: 
	 *    - \f$ b \f$
	 *    - \f$ shot noise \f$
	 *    
	 *  the parameters of the true dN/dz distribution ("gaussian",
	 *  "skewnormal" or "polynomial") in a given photometric
	 *  redshift bin are fixed by the user.  The angular power
	 *  spectrum is computed using the input cosmological
	 *  parameters
	 *
	 *  @param cosmo_param list of cosmological parameter for the model
	 *
	 *  @param cosmo_prior prior for the cosmological parameters
	 *
	 *  @param bias_prior prior for the parameter \f$ b \f$
	 *
	 *  @param shot_noise_prior prior for the parameter \f$ shot noise \f$
	 *  
	 */
	void set_model (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution bias_prior, const statistics::PriorDistribution shot_noise_prior);

	/**
	 *  @brief set the model to fit the cross correlation angular power spectrum
	 *
	 *  The model is the following:
	 *
	 *  \f[ C_l = \frac{2}{\pi}\int_{0}^{\infty} dk k^2 P_{mat}(k)\left(\psi_l(k)+\psi_l^r(k)\right)^2 \f]
	 *
	 *  where the kernel functions are:
	 *
	 *  \f[\psi_l(k) = \int_{0}^{\infty} dz \frac{1}{N}\frac{dN}{dz} b(z) f(z) D(z)j_l(kr(z)) \f]
	 *
	 *  \f[\psi_l^r(k) = \int_{0}^{\infty} dz \frac{1}{N}\frac{dN}{dz} f(z) D(z)\left[\frac{2l^2+2l-1}{(2l+3)(2l-1)}j_l(kr(z))-\frac{l(l-1)}{(2l+1)(2l-1)}j_{l-2}(kr(z))-\frac{(l+1)(l+2)}{(2l+1)(2l+3)}j_{l+2}(kr(z))\right] \f]
	 *  where the latter accounts for the contribution of redshift-space distorsion (RSD), f(z) is the linear growth rate, D(z) is the linear growth factor. 
	 *  In the Limber approximation this formula can be extended to l<2 and it reduces to
	 *
	 *  \f[ C_l = \frac{b^2}{N^2}\int_{0}^{\infty} \frac{dN}{dz}\frac{dN}{dz}
	 *  P_{mat} \left(\frac{l+1/2}{r(z)}\right) \frac{H(z)}{c}\frac{1}{r^2(z)} dz\f]
	 *
	 *
	 *  the model has N cosmological parameters plus 2 parameters: 
	 *    - \f$ b_i \f$
	 *    - \f$ b_j \f$
	 *    
	 *  the parameters of the true dN/dz distribution ("gaussian", "skewnormal" or "polynomial") in a given photometric redshift bin are fixed by the user.
	 *  the angular power spectrum is computed using the input cosmological parameters
	 *
	 *  @param cosmo_param list of cosmological parameter for the model
	 *
	 *  @param cosmo_prior prior for the cosmological parameters
	 *
	 *  @param bias_prior_bin1 prior for the parameter \f$ b_i \f$
	 *
	 *  @param bias_prior_bin2 prior for the parameter \f$ b_j \f$
	 *  
	 */
	void set_model_cross (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution bias_prior_bin1, const statistics::PriorDistribution bias_prior_bin2);

	/**
	 *  @brief set the model to fit the angular power spectrum
	 *
	 *  The model is the following:
	 *
	 *  \f[ C_l = \frac{2}{\pi}\int_{0}^{\infty} dk k^2 P_{mat}(k)\left(\psi_l(k)+\psi_l^r(k)\right)^2 \f]
	 *
	 *  where the kernel functions are:
	 *
	 *  \f[\psi_l(k) = \int_{0}^{\infty} dz \frac{1}{N}\frac{dN}{dz} b(z) f(z) D(z)j_l(kr(z)) \f]
	 *
	 *  \f[\psi_l^r(k) = \int_{0}^{\infty} dz \frac{1}{N}\frac{dN}{dz} f(z) D(z)\left[\frac{2l^2+2l-1}{(2l+3)(2l-1)}j_l(kr(z))-\frac{l(l-1)}{(2l+1)(2l-1)}j_{l-2}(kr(z))-\frac{(l+1)(l+2)}{(2l+1)(2l+3)}j_{l+2}(kr(z))\right] \f]
	 *  where the latter accounts for the contribution of redshift-space distorsion (RSD), f(z) is the linear growth rate, D(z) is the linear growth factor. 
	 *  In the Limber approximation this formula can be extended to l<2 and it reduces to
	 *
	 *  \f[ C_l = \frac{b^2}{N^2}\int_{0}^{\infty} \frac{dN}{dz}\frac{dN}{dz}
	 *  P_{mat} \left(\frac{l+1/2}{r(z)}\right) \frac{H(z)}{c}\frac{1}{r^2(z)} dz\f]
	 *
	 *  the model has N cosmological parameters plus 1 parameter: 
	 *    - \f$ b \f$
	 *    - \f$ shot noise \f$
	 *    
	 *  the dN/dz distribution in computed from the theoretical mass
	 *  function:
	 *
	 *  \f[ \frac{{\rm d}N }{{\rm d} z_{\rm tr}}(\Delta z_{\text{ob}})=\,\,\Omega 
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{\Delta z_{\text{ob}}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  , \f]
	 *  The angular power spectrum is computed
	 *  using the input cosmological parameters
	 *
	 *  @param cosmo_param list of cosmological parameter for the model
	 *
	 *  @param cosmo_prior prior for the cosmological parameters
	 *
	 *  @param bias_prior prior for the parameter \f$ b \f$
	 *
	 *  @param shot_noise_prior prior for the parameter \f$ shot noise \f$
	 *  
	 */
	void set_model_MF_theory (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution bias_prior, const statistics::PriorDistribution shot_noise_prior);
	
	/**
	 *  @brief set the model to fit the angular power spectrum
	 *
	 *  The model is the following:
	 *
	 *  \f[ C_l = \frac{2}{\pi}\int_{0}^{\infty} dk k^2 P_{mat}(k)\left(\psi_l(k)+\psi_l^r(k)\right)^2 \f]
	 *
	 *  where the kernel functions are:
	 *
	 *  \f[\psi_l(k) = \int_{0}^{\infty} dz \frac{1}{N}\frac{dN}{dz} b(z) f(z) D(z)j_l(kr(z)) \f]
	 *
	 *  \f[\psi_l^r(k) = \int_{0}^{\infty} dz \frac{1}{N}\frac{dN}{dz} f(z) D(z)\left[\frac{2l^2+2l-1}{(2l+3)(2l-1)}j_l(kr(z))-\frac{l(l-1)}{(2l+1)(2l-1)}j_{l-2}(kr(z))-\frac{(l+1)(l+2)}{(2l+1)(2l+3)}j_{l+2}(kr(z))\right] \f]
	 *  where the latter accounts for the contribution of redshift-space distorsion (RSD), f(z) is the linear growth rate, D(z) is the linear growth factor. 
	 *  In the Limber approximation this formula can be extended to l<2 and it reduces to
	 *
	 *  \f[ C_l = \frac{b^2}{N^2}\int_{0}^{\infty} \frac{dN}{dz}\frac{dN}{dz}
	 *  P_{mat} \left(\frac{l+1/2}{r(z)}\right) \frac{H(z)}{c}\frac{1}{r^2(z)} dz\f]
	 *
	 *  the model has N cosmological parameters
	 *    
	 *  the dN/dz distribution in computed from the theoretical mass
	 *  function:
	 *
	 *  \f[ \frac{{\rm d}N }{{\rm d} z_{\rm tr}}(\Delta z_{\text{ob}})=\,\,\Omega 
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{\Delta z_{\text{ob}}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  , \f]
	 *  The angular power spectrum is computed
	 *  using the input cosmological parameters
	 *
	 *  @param cosmo_param list of cosmological parameter for the model
	 *
	 *  @param cosmo_prior prior for the cosmological parameters
	 *
	 *  @param shot_noise_prior prior for the shot noise
	 *  
	 */
	void set_model_bias_MF_theory (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution shot_noise_prior);
	
/**
	 *  @brief set the model to fit the angular power spectrum
	 *
	 *  The model is the following:
	 *
	 *  \f[ C_l = \frac{2}{\pi}\int_{0}^{\infty} dk k^2 P_{mat}(k)\left(\psi_l(k)+\psi_l^r(k)\right)^2 \f]
	 *
	 *  where the kernel functions are:
	 *
	 *  \f[\psi_l(k) = \int_{0}^{\infty} dz \frac{1}{N}\frac{dN}{dz} b(z) f(z) D(z)j_l(kr(z)) \f]
	 *
	 *  \f[\psi_l^r(k) = \int_{0}^{\infty} dz \frac{1}{N}\frac{dN}{dz} f(z) D(z)\left[\frac{2l^2+2l-1}{(2l+3)(2l-1)}j_l(kr(z))-\frac{l(l-1)}{(2l+1)(2l-1)}j_{l-2}(kr(z))-\frac{(l+1)(l+2)}{(2l+1)(2l+3)}j_{l+2}(kr(z))\right] \f]
	 *  where the latter accounts for the contribution of redshift-space distorsion (RSD), f(z) is the linear growth rate, D(z) is the linear growth factor. 
	 *  In the Limber approximation this formula can be extended to l<2 and it reduces to
	 *
	 *  \f[ C_l = \frac{b^2}{N^2}\int_{0}^{\infty} \frac{dN}{dz}\frac{dN}{dz}
	 *  P_{mat} \left(\frac{l+1/2}{r(z)}\right) \frac{H(z)}{c}\frac{1}{r^2(z)} dz\f]
	 *
	 *  the model has N cosmological parameters, plus: 
	 *    - \f$\sigma}\f$ normalised shot noise 
	 *    - \f$\Sigma_{NL}\f$ damping at BAO
	 *    
	 *  the dN/dz distribution in computed from the theoretical mass
	 *  function:
	 *
	 *  \f[ \frac{{\rm d}N }{{\rm d} z_{\rm tr}}(\Delta z_{\text{ob}})=\,\,\Omega 
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{\Delta z_{\text{ob}}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  , \f]
	 *  The angular power spectrum is computed
	 *  using the input cosmological parameters
	 *
	 *  @param cosmo_param list of cosmological parameter for the model
	 *
	 *  @param cosmo_prior prior for the cosmological parameters
	 *
	 *  @param shot_noise_prior prior for the shot noise
	 *
	 *  @param sigmaNL_prior prior for the damping at BAO scale
	 *  
	 */
	void set_model_bias_MF_theory_BAO (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution shot_noise_prior, const statistics::PriorDistribution sigmaNL_prior);
	
	/**
	 *  @brief Set the parameters to model the monopole of the
	 *  two-point correlation function in redshift space, where
	 *  the effective bias can be computed in two ways, depending
	 *  on the cbl::modelling::twopt::Modelling_TwoPointCorrelation1D::set\_data\_model used.
	 *  In one case, the masses are given by a mass-observable scaling relation
	 *  with the following functional form:
	 *
	 *  \f[ \log \frac{M}{M_{\rm piv}} = \alpha+
	 *  \beta  \log \frac{\lambda}{\lambda_{\rm piv}}+
	 *  \gamma  \log f(z;z_{\rm piv}) + \sigma_{\rm intr}, \f]
	 *
	 *  where \f$\lambda\f$ is the mass proxy and
	 *  the intrinsic scatter, \f$\sigma_{\rm intr}\f$, should
	 *  be set as both negative and positive and it has the 
	 *  following functional form:
	 *
	 *  \f$ \sigma_{\rm intr} = \sigma_0 + \sigma_{M} 
	 *  \log (M/M_{\rm piv})^{e_{M}} + \sigma_z \log (f(z))^{e_z}. \f$
	 *
	 *  The bias is then computed for each object in correspondence of
	 *  such masses, and then averaged to obtain the effective bias.
	 *
	 *  In the second case, the effective bias is derived as follows:
	 *
	 *  \f$ b_{\rm eff} = \frac{1}{N}\sum_{i=1}^N b(z_{{\rm ob},i},\lambda_{{\rm ob},i}) = 
	 *  \frac{1}{N}\sum_{i=1}^N\int_0^\infty {\rm d}z\,\int_0^\infty{\rm d}\lambda\,\int_0^\infty{\rm d}M\,\, b(M,z) \, 
	 *  P(M|\lambda,z) \, P(z|z_{{\rm ob},i}) \, P(\lambda|\lambda_{{\rm ob},i}), \f$
	 *
	 *  where \f$N\f$ is the number of objects, \f$z_{{\rm ob},i}\f$ and \f$\lambda_{{\rm ob},i}\f$
	 *  are the observed redshift and mass proxy of the \f$i\f$th object, respectively,
	 *  \f$P(M|\lambda,z)\f$ is a log-normal whose mean is the mass proxy - mass relation and whose
	 *  rms is the intrinsic scatter of such relation (for details, see e.g. 
	 *  cbl::modelling::massobsrel::Modelling_MassObservableRelation), while \f$P(z|z_{{\rm ob},i})\f$
	 *  and \f$P(\lambda|\lambda_{{\rm ob},i})\f$ are Gaussian distributions whose standard deviations
	 *  are given by the errors on redshift and proxy, respectively.
	 *
	 *  The model is the following:
	 *
	 *  \f[ C_l = \frac{2}{\pi}\int_{0}^{\infty} dk k^2 P_{mat}(k)\left(\psi_l(k)+\psi_l^r(k)\right)^2 \f]
	 *
	 *  where the kernel functions are:
	 *
	 *  \f[\psi_l(k) = \int_{0}^{\infty} dz \frac{1}{N}\frac{dN}{dz} b(z) f(z) D(z)j_l(kr(z)) \f]
	 *
	 *  \f[\psi_l^r(k) = \int_{0}^{\infty} dz \frac{1}{N}\frac{dN}{dz} f(z) D(z)\left[\frac{2l^2+2l-1}{(2l+3)(2l-1)}j_l(kr(z))-\frac{l(l-1)}{(2l+1)(2l-1)}j_{l-2}(kr(z))-\frac{(l+1)(l+2)}{(2l+1)(2l+3)}j_{l+2}(kr(z))\right] \f]
	 *  where the latter accounts for the contribution of redshift-space distorsion (RSD), f(z) is the linear growth rate, D(z) is the linear growth factor. 
	 *  In the Limber approximation this formula can be extended to l<2 and it reduces to
	 *
	 *  \f[ C_l = \frac{b^2}{N^2}\int_{0}^{\infty} \frac{dN}{dz}\frac{dN}{dz}
	 *  P_{mat} \left(\frac{l+1/2}{r(z)}\right) \frac{H(z)}{c}\frac{1}{r^2(z)} dz\f]
	 *
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_prior vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @param alpha_prior prior on the scaling relation normalization
	 *
	 *  @param beta_prior prior on the scaling relation slope
	 *
	 *  @param gamma_prior prior on the redshift evolution factor of the scaling relation
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
	 *  @param shot_noise_prior prior on the shot noise
	 *
	 *  @param z_evo functional form of the redshift evolution
	 *  function in the scaling relation. See the documentation
	 *  of cbl::modelling::massobsrel::Modelling_MassObservableRelation
	 *  
	 */
	void set_model_scaling_relation_cosmology (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const statistics::PriorDistribution shot_noise_prior, const std::string z_evo);

	/**
	 *  @brief Set the parameters to model the angular power
	 *  spectrum in real or redshift space, where
	 *  the effective bias can be computed in two ways, depending
	 *  on the cbl::modelling::angularpk::Modelling_PowerSpectrum_angular::set\_data\_model used.
	 *  In one case, the masses are given by a mass-observable scaling relation
	 *  with the following functional form:
	 *
	 *  \f[ \log \frac{M}{M_{\rm piv}} = \alpha+
	 *  \beta  \log \frac{\lambda}{\lambda_{\rm piv}}+
	 *  \gamma  \log f(z;z_{\rm piv}) + \sigma_{\rm intr}, \f]
	 *
	 *  where \f$\lambda\f$ is the mass proxy and
	 *  the intrinsic scatter, \f$\sigma_{\rm intr}\f$, should
	 *  be set as both negative and positive and it has the 
	 *  following functional form:
	 *
	 *  \f$ \sigma_{\rm intr} = \sigma_0 + \sigma_{M} 
	 *  \log (M/M_{\rm piv})^{e_{M}} + \sigma_z \log (f(z))^{e_z}. \f$
	 *
	 *  The bias is then computed for each object in correspondence of
	 *  such masses, and then averaged to obtain the effective bias.
	 *
	 *  In the second case, the effective bias is derived as follows:
	 *
	 *  \f$ b_{\rm eff} = \frac{1}{N}\sum_{i=1}^N b(z_{{\rm ob},i},\lambda_{{\rm ob},i}) = 
	 *  \frac{1}{N}\sum_{i=1}^N\int_0^\infty {\rm d}z\,\int_0^\infty{\rm d}\lambda\,\int_0^\infty{\rm d}M\,\, b(M,z) \, 
	 *  P(M|\lambda,z) \, P(z|z_{{\rm ob},i}) \, P(\lambda|\lambda_{{\rm ob},i}), \f$
	 *
	 *  where \f$N\f$ is the number of objects, \f$z_{{\rm ob},i}\f$ and \f$\lambda_{{\rm ob},i}\f$
	 *  are the observed redshift and mass proxy of the \f$i\f$th object, respectively,
	 *  \f$P(M|\lambda,z)\f$ is a log-normal whose mean is the mass proxy - mass relation and whose
	 *  rms is the intrinsic scatter of such relation (for details, see e.g. 
	 *  cbl::modelling::massobsrel::Modelling_MassObservableRelation), while \f$P(z|z_{{\rm ob},i})\f$
	 *  and \f$P(\lambda|\lambda_{{\rm ob},i})\f$ are Gaussian distributions whose standard deviations
	 *  are given by the errors on redshift and proxy, respectively.
	 *
	 *  The model is the following:
	 *
	 *  \f[ C_l = \frac{2}{\pi}\int_{0}^{\infty} dk k^2 P_{mat}(k)\left(\psi_l(k)+\psi_l^r(k)\right)^2 \f]
	 *
	 *  where the kernel functions are:
	 *
	 *  \f[\psi_l(k) = \int_{0}^{\infty} dz \frac{1}{N}\frac{dN}{dz} b(z) f(z) D(z)j_l(kr(z)) \f]
	 *
	 *  \f[\psi_l^r(k) = \int_{0}^{\infty} dz \frac{1}{N}\frac{dN}{dz} f(z) D(z)\left[\frac{2l^2+2l-1}{(2l+3)(2l-1)}j_l(kr(z))-\frac{l(l-1)}{(2l+1)(2l-1)}j_{l-2}(kr(z))-\frac{(l+1)(l+2)}{(2l+1)(2l+3)}j_{l+2}(kr(z))\right] \f]
	 *  where the latter accounts for the contribution of redshift-space distorsion (RSD), f(z) is the linear growth rate, D(z) is the linear growth factor. 
	 *  In the Limber approximation this formula can be extended to l<2 and it reduces to
	 *
	 *  \f[ C_l = \frac{b^2}{N^2}\int_{0}^{\infty} \frac{dN}{dz}\frac{dN}{dz}
	 *  P_{mat} \left(\frac{l+1/2}{r(z)}\right) \frac{H(z)}{c}\frac{1}{r^2(z)} dz\f]
	 *  Finally, 
	 *  \f[ \frac{{\rm d}N }{{\rm d} z_{\rm tr}}(\Delta z_{\text{ob}, j})=\,\,\Omega 
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{0}^{\infty}{\rm d}\lambda_{\rm tr}\,\,
	 *  P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\,
	 *  \int_{\Delta z_{\text{ob},j}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  \int_{\Delta\lambda_{\text{ob},i}}{\rm d} \lambda_{\rm ob} 
	 *  \,\,P(\lambda_{\rm ob}|\lambda_{\rm tr})w(\lambda_{ob},\Delta z_{\text{ob},j}), \f]
	 *
	 *  as expressed in cbl::modelling::numbercounts::Modelling_NumberCounts1D_MassProxy::set_model_NumberCounts_cosmology
	 *
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_prior vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @param alpha_prior prior on the scaling relation normalization
	 *
	 *  @param beta_prior prior on the scaling relation slope
	 *
	 *  @param gamma_prior prior on the redshift evolution factor of the scaling relation
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
	 *  @param Plambda_prior priors on the three parameters,
	 *  respectively a, b, c, defining the distribution
	 *  \f$P(\lambda_{tr}|z_{tr})\f$
	 *
	 *  @param shot_noise_prior prior on the shot noise
	 *
	 *  @param z_evo functional form of the redshift evolution
	 *  function in the scaling relation. See the documentation
	 *  of cbl::modelling::massobsrel::Modelling_MassObservableRelation
	 *  
	 */
	void set_model_scaling_relation_cosmology_MF_theory (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const std::vector<statistics::PriorDistribution> Plambda_prior, const statistics::PriorDistribution shot_noise_prior, const std::string z_evo);
	
	/**
	 * @brief get the member \e m_data_model
	 * @return the container of parameters for two-point
	 * correlation function model computation
	 */
	std::shared_ptr<modelling::angularpk::STR_data_model> data_model () { return m_data_model; }
	
	/**
	 *  @brief Set the data used to construct models of
	 *  the angular power spectrum.
	 *
	 *  The model is based on the redshift distribution obtained directly from the data (e.g. galaxy catalogues). The true redshift distribution in a given redshift bin is obtained as: 

	 \f[\frac{dN}{dz_{true}} = \int_{0}^{\infty} dz_{phot} W(z_{phot}) \frac{dN}{dz_{phot}} P(z_{true}|z_{phot})  \f]
	 *
	 *  where \f$ W(z_{phot}) \f$ is a top-hat window function in the given photometric bin and \f$ P(z_{true}|z_{phot}) \f$ is assumed as gaussian with zero mean and std dev given by \f$ \sigma_{0,z}(1+z_{phot}) \f$. See Balaguera et al. (2018)
	 *
	 *  @param cosmology the cosmological model used to compute the angular power spectrum
	 *
	 *  @param z_min the minimum redshift of the photometric redshift slice
	 *
	 *  @param z_max the maximum redshift of the photometric redshift slice
	 *
	 *  @param method_Pk method used to compute the power spectrum
	 *  (i.e. the Boltzmann solver); valid choices for method_Pk
	 *  are: CAMB [http://camb.info/], CLASS
	 *  [http://class-code.net/], MPTbreeze-v1
	 *  [http://arxiv.org/abs/1207.1465], EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *
	 *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
	 *  non-linear power spectrum
	 *
	 *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
	 *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
	 *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
	 *
	 *  @param k_min minimum wave vector module up to which the
	 *  power spectrum is computed in order to estimate the power
	 *  spectrum normalisation; this parameter is used only if
	 *  either norm=1, or norm=-1 and sigma8 is set
	 *
	 *  @param k_max maximum wave vector module up to which the
	 *  power spectrum is computed to estimate the power spectrum
	 *  normalisation; this parameter is used only if norm=1
	 *
	 *  @param dN_par vector containing the parameters of the true distribution of objects dN/dz in a given photometric bin 
	 *
	 *  @param distribution_type the shape of the dN/dz distribution. It can be "polynomial", if you simply use a fit of photometric redshift distribution. In the case of convolution with P(z_true|z_phot) you should use "gaussian" or "skewnormal"
	 *
	 *  @param fsky the fraction of the sky covered by the survey
	 *
	 *  @param limber_approximation true->use the limber approximation (start from l=0), false->don't use the limber approximation (start from l=2)
	 *
	 *  @param ell_limber for l<ell_limber use exact formula,
	 *  for l>ell_limber use Limber approximation
	 *
	 *  @param ll the vector of the multipoles of the mixing matrix
	 *
	 *  @param mixing_matrix the mixing matrix
	 *
	 *  @param redshift_space true->redshift space, include
	 *  redshift space distorsion, false->real space
	 *
	 *  @param z_min_bin2 the minimum redshift of the second photometric redshift slice, for cross-correlation
	 *
	 *  @param z_max_bin2 the maximum redshift of the second redshift slice, for cross-correlation
	 *
	 *  @param dN_par_bin2 vector containing the parameters of the
	 *  true distribution of objects dN/dz in the second
	 *  photometric bin, for cross-correlation
	 *
	 *  @param interpType interpolation type
	 *
	 *  @param shot_noise_normalization the shot noise normalization. Shot noise il multiplyed for this value before summing with \f[C_\ell \f]
	 *
	 */
	void set_data_model (const std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const double z_min=par::defaultDouble, const double z_max=par::defaultDouble, const std::string method_Pk="CAMB", const bool NL=false, const int norm=-1, const double k_min=0.001, const double k_max=5., const std::vector<double> dN_par={}, const std::string distribution_type="skewnormal", const double fsky=1., const bool limber_approximation=true, const double ell_limber=20., std::vector<double> ll={}, std::vector<std::vector<double>> mixing_matrix={}, const bool redshift_space=false, const double z_min_bin2=par::defaultDouble, const double z_max_bin2=par::defaultDouble, const std::vector<double> dN_par_bin2={}, const std::string interpType="Spline", double shot_noise_normalization=1.);

	/**
	 *  @brief Set the data used to construct models of
	 *  the angular power spectrum, in which the redshift selection
	 *  function is computed through a theoretical mass function.
	 *
	 *  @param z_abs_err the redshift absolute error. Zero for spectroscopic survey, >0 for photometric survey.
	 *
	 *  @param cosmology the cosmological model used to compute the power spectrum
	 *
	 *  @param z_min the minimum redshift of the photometric redshift slice
	 *
	 *  @param z_max the maximum redshift of the photometric redshift slice
	 *
	 *  @param method_Pk method used to compute the power
	 *  spectrum; valid choices for method_Pk are: CAMB
	 *  [http://camb.info/], CLASS [http://class-code.net/],
	 *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465],
	 *  EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *
	 *  @param NL 0 &rarr; linear power spectrum; 1 &rarr;
	 *  non-linear power spectrum
	 *
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *
	 *  @param k_min minimum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *  
	 *  @param k_max maximum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *
	 *  @param area_degrees the area in deg. covered by the survey
	 *
	 *  @param limber_approximation true->use the limber approximation (ell will start from l=0), false->don't use the limber approximation (ell will start from l=2)
	 *
	 *  @param ell_limber for l<ell_limber use exact formula,
	 *  for l>ell_limber use Limber approximation
	 *
	 *  @param ll the vector of multipoles of the mixing matrix
	 *
	 *  @param mixing_matrix the mixing matrix
	 *
	 *  @param redshift_space true->redshift space, false->real space
	 *
	 *  @param interpType the interpolation type 
	 *
	 *  @param Delta \f$\Delta\f$, the overdensity
	 *
	 *  @param isDelta_critical \f$\rightarrow\f$ \f$\Delta\f$ is the
	 *  overdensity defined with respect to the critical density
	 *
	 *  @param model_MF author(s) who proposed the mass function;
	 *  valid authors are: PS (Press & Schechter), ST(Sheth & Tormen),
	 *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
	 *  ShenH (halo MF, Shen et al. 2006), ShenF
	 *  (filament MF, Shen et al. 2006), ShenS
	 *  (sheet MF, Shen et al. 2006), Pan (Pan 2007), Peacock
	 *  (Peacock at al. 2007), Reed (Reed et al. 2007), Tinker
	 *  (Tinker et al. 2008), Crocce (Crocce et al. 2010), Manera
	 *  (Manera et al. 2010), Courtin (Courtin et al. 2010),
	 *  Bhattacharya (Bhattacharya et al. 2011), Angulo_FOF
	 *  (FoF MF, Angulo et al. 2012), Angulo_Sub
	 *  (SUBFIND MF, Angulo et al. 2012), Watson_FOF
	 *  (FoF MF, Watson et al. 2012),
	 *  Watson_SOH (Spherical Overdensity halo MF, Watson et al. 2012)
	 *, Despali_Z0, Despali_AllZ, Despali_AllZAllCosmo,
	 *  Despali_HighM (Despali et al. 2016)
	 *
	 *  @param Mass_min minimum mass
	 *  
	 *  @param Mass_max maximum mass
	 *
	 *  @param Mass_step the number of mass step
	 *
	 *  @param shot_noise_normalization the shot noise normalization. Shot noise il multiplyed for this value before summing with \f[C_\ell \f]
	 *
	 */
	void set_data_model_MF_theory (const double z_abs_err, const std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const double z_min=par::defaultDouble, const double z_max=par::defaultDouble, const std::string method_Pk="CAMB", const bool NL=false, const int norm=-1, const double k_min=0.001, const double k_max=5., const double area_degrees=cbl::par::defaultDouble, const bool limber_approximation=true, const double ell_limber=20., std::vector<double> ll={}, std::vector<std::vector<double>> mixing_matrix={}, const bool redshift_space=false, const std::string interpType="Spline", const double Delta = 200., const bool isDelta_critical = true, const std::string model_MF = "Tinker", const double Mass_min = par::defaultDouble, const double Mass_max = par::defaultDouble, const int Mass_step = 100, double shot_noise_normalization=1.);
	
	/**
	 *  @brief Set the data used to construct models of
	 *  the angular power spectrum, in which the bias and the redshift selection
	 *  function are computed theorically computed
	 *
	 *  @param z_abs_err the redshift absolute error. Zero for spectroscopic survey, >0 for photometric survey.
	 *
	 *  @param cosmology the cosmological model used to compute the power spectrum
	 *
	 *  @param cluster_mass vector conteining the masses of the catalogue
	 *
	 *  @param cluster_redshift vector conteining the redshifts of the catalogue
	 *
	 *  @param z_min the minimum redshift of the photometric redshift slice
	 *
	 *  @param z_max the maximum redshift of the photometric redshift slice
	 *
	 *  @param method_Pk method used to compute the power
	 *  spectrum; valid choices for method_Pk are: CAMB
	 *  [http://camb.info/], CLASS [http://class-code.net/],
	 *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465],
	 *  EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *
	 *  @param NL 0 &rarr; linear power spectrum; 1 &rarr;
	 *  non-linear power spectrum
	 *
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *
	 *  @param k_min minimum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *  
	 *  @param k_max maximum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *
	 *  @param output_dir the output_dir directory
	 *  where the output of external codes are written
	 *
	 *  @param output_root output_root of the parameter file used
	 *  to compute the power spectrum and &sigma;(mass); it can be
	 *  any name
	 *
	 *  @param area_degrees the area in deg. covered by the survey
	 *
	 *  @param limber_approximation true->use the limber approximation (ell will start from l=0), false->don't use the limber approximation (ell will start from l=2)
	 *
	 *  @param ell_limber for l<ell_limber use exact formula,
	 *  for l>ell_limber use Limber approximation
	 *
	 *  @param ll the vector of multipoles of the mixing matrix
	 *
	 *  @param mixing_matrix the mixing matrix
	 *
	 *  @param redshift_space true->redshift space, false->real space
	 *
	 *  @param interpType the interpolation type 
	 *
	 *  @param Delta \f$\Delta\f$, the overdensity
	 *
	 *  @param isDelta_critical \f$\rightarrow\f$ \f$\Delta\f$ is the
	 *  overdensity defined with respect to the critical density
	 *
	 *  @param model_MF author(s) who proposed the mass function;
	 *  valid authors are: PS (Press & Schechter), ST(Sheth & Tormen),
	 *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
	 *  ShenH (halo MF, Shen et al. 2006), ShenF
	 *  (filament MF, Shen et al. 2006), ShenS
	 *  (sheet MF, Shen et al. 2006), Pan (Pan 2007), Peacock
	 *  (Peacock at al. 2007), Reed (Reed et al. 2007), Tinker
	 *  (Tinker et al. 2008), Crocce (Crocce et al. 2010), Manera
	 *  (Manera et al. 2010), Courtin (Courtin et al. 2010),
	 *  Bhattacharya (Bhattacharya et al. 2011), Angulo_FOF
	 *  (FoF MF, Angulo et al. 2012), Angulo_Sub
	 *  (SUBFIND MF, Angulo et al. 2012), Watson_FOF
	 *  (FoF MF, Watson et al. 2012),
	 *  Watson_SOH (Spherical Overdensity halo MF, Watson et al. 2012)
	 *, Despali_Z0, Despali_AllZ, Despali_AllZAllCosmo,
	 *  Despali_HighM (Despali et al. 2016)
	 *
	 *  @param model_bias author(s) who proposed the bias; valid
	 *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
	 *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
	 *  correction of Warren 2004), Tinker (Tinker et al. 2010)
	 *
	 *  @param Mass_step the number of mass step
	 *
	 *  @param bias_from_sample if true compute the eff bias as an average of the clusters bias from the catalogue
	 *
	 *  @param shot_noise_normalization the shot noise normalization. Shot noise il multiplyed for this value before summing with \f[C_\ell \f]
	 *
	 */
	void set_data_model_bias_MF_theory (const double z_abs_err, const std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const std::vector<double> cluster_mass, const std::vector<double> cluster_redshift, const double z_min=par::defaultDouble, const double z_max=par::defaultDouble, const std::string method_Pk="CAMB", const bool NL=false, const int norm=-1, const double k_min=0.001, const double k_max=100., const std::string output_dir=par::defaultString, const std::string output_root="test", const double area_degrees=cbl::par::defaultDouble, const bool limber_approximation=true, const double ell_limber=20., std::vector<double> ll={}, std::vector<std::vector<double>> mixing_matrix={}, const bool redshift_space=false, const std::string interpType="Spline", const double Delta = 200., const bool isDelta_critical = true, const std::string model_MF = "Tinker", const std::string model_bias = "Tinker", const int Mass_step = 100, const bool bias_from_sample=false, double shot_noise_normalization=1.);

	/**
	 *  @brief Set the data used to construct models of
	 *  the angular power spectrum, in which the effective bias
	 *  is computed through a mass-observable scaling relation.
	 *
	 *  In particular, the cluster masses are derived from the following
	 *  scaling relation:
	 *
	 *  \f$ \log (M/M_{\rm piv}) = \alpha + \beta\log(\lambda/\lambda_{\rm piv})+
	 *  \gamma \log f(z;z_{\rm piv}) \f$
	 *
	 *  See cbl::modelling::massobsrel::Modelling_MassObservableRelation for details.
	 *
	 *  The model is based on the redshift distribution obtained directly from the data (e.g. galaxy catalogues). The true redshift distribution in a given redshift bin is obtained as: 

	 \f[\frac{dN}{dz_{true}} = \int_{0}^{\infty} dz_{phot} W(z_{phot}) \frac{dN}{dz_{phot}} P(z_{true}|z_{phot})  \f]
	 *
	 *  where \f$ W(z_{phot}) \f$ is a top-hat window function in the given photometric bin and \f$ P(z_{true}|z_{phot}) \f$ is assumed as gaussian with zero mean and std dev given by \f$ \sigma_{0,z}(1+z_{phot}) \f$
	 *
	 *  @param z_abs_err the redshift absolute error. Zero for spectroscopic survey, >0 for photometric survey.
	 *
	 *  @param proxy_rel_err the proxy relative error
	 *
	 *  @param cosmology the cosmological model used to compute the power spectrum
	 *
	 *  @param redshift redshift
	 *
	 *  @param cluster_redshift vector containing the cluster
	 *  redshifts
	 *
	 *  @param cluster_mass_proxy vector containing the cluster
	 *  mass proxies
	 *
	 *  @param redshift_pivot redshift pivot in the scaling relation
	 *
	 *  @param proxy_pivot proxy pivot in the scaling relation
	 *
	 *  @param mass_pivot mass pivot in the scaling relation
	 *
	 *  @param log_base base of the logarithms in the scaling relation
	 *
	 *  @param method_Pk method used to compute the power
	 *  spectrum; valid choices for method_Pk are: CAMB
	 *  [http://camb.info/], CLASS [http://class-code.net/],
	 *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465],
	 *  EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *
	 *  @param NL 0 &rarr; linear power spectrum; 1 &rarr;
	 *  non-linear power spectrum
	 *
	 *  @param k_min minimum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *  
	 *  @param k_max maximum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *
	 *  @param dN_par the fitting parameter of the redshift
	 *  distribution
	 *
	 *  @param distribution_type the shape of the dN/dz distribution.
	 *  Can be "polynomial", for a fit of the photometric redshift
	 *  distribution, "gaussian" or "skewnormal", in the correct
	 *  case of convolution with P(z true|z_phot)
	 *
	 *  @param z_min the minimum redshift of the photometric redshift slice
	 *
	 *  @param z_max the maximum redshift of the photometric redshift slice
	 *
	 *  @param fsky the fraction of sky covered by the survey
	 *
	 *  @param limber_approximation true->use the limber approximation (ell will start from l=0), false->don't use the limber approximation (ell will start from l=2)
	 *
	 *  @param ell_limber for l<ell_limber use exact formula,
	 *  for l>ell_limber use Limber approximation
	 *
	 *  @param ll the vector of multipoles of the mixing matrix
	 *
	 *  @param mixing_matrix the mixing matrix
	 *
	 *  @param redshift_space true->redshift space, false->real space
	 *
	 *  @param interpType interpolation type
	 *
	 *  @param output_dir the output_dir directory
	 *  where the output of external codes are written
	 *
	 *  @param output_root output_root of the parameter file used
	 *  to compute the power spectrum and &sigma;(mass); it can be
	 *  any name
	 *  
	 *  @param prec accuracy of the GSL integration
	 *  
	 *  @param file_par name of the parameter file; if a parameter
	 *  file is provided (i.e. file_par!=NULL), it will be used,
	 *  ignoring the cosmological parameters of the object
	 *
	 *  @param Delta \f$\Delta\f$, the overdensity
	 *
	 *  @param isDelta_critical \f$\rightarrow\f$ \f$\Delta\f$ is the
	 *  overdensity defined with respect to the critical density
	 *
	 *  @param model_bias author(s) who proposed the bias; valid
	 *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
	 *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
	 *  correction of Warren 2004), Tinker (Tinker et al. 2010)
	 *
	 *  @param shot_noise_normalization the shot noise normalization. Shot noise il multiplyed for this value before summing with \f[C_\ell \f]
	 *
	 */
	void set_data_model_scaling_relation_cosmology (const double z_abs_err, const double proxy_rel_err, const std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const double redshift, const std::vector<double> cluster_redshift, const std::vector<double> cluster_mass_proxy, const double redshift_pivot, const double proxy_pivot, double mass_pivot, const double log_base, const std::string method_Pk="CAMB", const bool NL=true, const double k_min=1.e-4, const double k_max=100., const int norm=-1, std::vector<double> dN_par={}, const std::string distribution_type="skewnormal", const double z_min=par::defaultDouble, const double z_max=par::defaultDouble, const double fsky=1., const bool limber_approximation=true, const double ell_limber=20., std::vector<double> ll={}, std::vector<std::vector<double>> mixing_matrix={}, const bool redshift_space=false, const std::string interpType="Spline", const std::string output_dir=par::defaultString, const std::string output_root="test", const double prec=1.e-3, const std::string file_par=par::defaultString, const double Delta=200., const bool isDelta_critical=true, const std::string model_bias="Tinker", double shot_noise_normalization=1.);

	/**
	 *  @brief Set the data used to construct models of
	 *  the two-point correlation function, in which the effective bias
	 *  is computed through a mass-observable scaling relation.
	 *
	 *  In particular, the cluster masses are derived from the following
	 *  scaling relation:
	 *
	 *  \f$ \log (M/M_{\rm piv}) = \alpha + \beta\log(\lambda/\lambda_{\rm piv})+
	 *  \gamma \log f(z;z_{\rm piv}) \f$
	 *
	 *  See cbl::modelling::massobsrel::Modelling_MassObservableRelation for details.
	 *
	 *  @param z_abs_err the redshift absolute error
	 *
	 *  @param proxy_rel_err the proxy relative error
	 *
	 *  @param cosmology the cosmological model used to compute the power spectrum
	 *
	 *  @param redshift redshift
	 *
	 *  @param cluster_redshift vector containing the cluster
	 *  redshifts
	 *
	 *  @param cluster_mass_proxy vector containing the cluster
	 *  mass proxies
	 *
	 *  @param redshift_pivot redshift pivot in the scaling relation
	 *
	 *  @param proxy_pivot proxy pivot in the scaling relation
	 *
	 *  @param mass_pivot mass pivot in the scaling relation
	 *
	 *  @param log_base base of the logarithms in the scaling relation
	 *
	 *  @param method_Pk method used to compute the power
	 *  spectrum; valid choices for method_Pk are: CAMB
	 *  [http://camb.info/], CLASS [http://class-code.net/],
	 *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465],
	 *  EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *
	 *  @param NL 0 &rarr; linear power spectrum; 1 &rarr;
	 *  non-linear power spectrum
	 *    
	 *  @param k_min minimum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *  
	 *  @param k_max maximum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *
	 *  @param area_degrees the area in deg. covered by the survey
	 *
	 *  @param limber_approximation true->use the limber approximation (ell will start from l=0), false->don't use the limber approximation (ell will start from l=2)
	 *
	 *  @param ell_limber for l<ell_limber use exact formula,
	 *  for l>ell_limber use Limber approximation
	 *
	 *  @param z_min the minimum redshift of the photometric redshift slice
	 *
	 *  @param z_max the maximum redshift of the photometric redshift slice
	 *
	 *  @param proxy_min the minimum proxy of the photometric redshift slice
	 *
	 *  @param proxy_max the maximum proxy of the photometric redshift slice
	 *
	 *  @param ll the vector of the multipoles of the mixing matrix
	 *
	 *  @param mixing_matrix the mixing matrix
	 *
	 *  @param redshift_space true->redshift space, include
	 *  redshift space distorsion, false->real space
	 *
	 *  @param interpType interpolation type
	 *
	 *  @param output_dir the output_dir directory
	 *  where the output of external codes are written
	 *
	 *  @param output_root output_root of the parameter file used
	 *  to compute the power spectrum and &sigma;(mass); it can be
	 *  any name
	 *
	 *  @param prec accuracy of the GSL integration
	 *  
	 *  @param file_par name of the parameter file; if a parameter
	 *  file is provided (i.e. file_par!=NULL), it will be used,
	 *  ignoring the cosmological parameters of the object
	 *
	 *  @param Delta \f$\Delta\f$, the overdensity
	 *
	 *  @param isDelta_critical \f$\rightarrow\f$ \f$\Delta\f$ is the
	 *  overdensity defined with respect to the critical density
	 *
	 *  @param model_bias author(s) who proposed the bias; valid
	 *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
	 *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
	 *  correction of Warren 2004), Tinker (Tinker et al. 2010)
	 *
	 *  @param model_MF author(s) who proposed the mass function;
	 *  valid authors are: PS (Press & Schechter), ST(Sheth & Tormen),
	 *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
	 *  ShenH (halo MF, Shen et al. 2006), ShenF
	 *  (filament MF, Shen et al. 2006), ShenS
	 *  (sheet MF, Shen et al. 2006), Pan (Pan 2007), Peacock
	 *  (Peacock at al. 2007), Reed (Reed et al. 2007), Tinker
	 *  (Tinker et al. 2008), Crocce (Crocce et al. 2010), Manera
	 *  (Manera et al. 2010), Courtin (Courtin et al. 2010),
	 *  Bhattacharya (Bhattacharya et al. 2011), Angulo_FOF
	 *  (FoF MF, Angulo et al. 2012), Angulo_Sub
	 *  (SUBFIND MF, Angulo et al. 2012), Watson_FOF
	 *  (FoF MF, Watson et al. 2012),
	 *  Watson_SOH (Spherical Overdensity halo MF, Watson et al. 2012)
	 *, Despali_Z0, Despali_AllZ, Despali_AllZAllCosmo,
	 *  Despali_HighM (Despali et al. 2016)
	 *
	 *  @param SF_weights weights for completeness and purity of
	 *  number counts in lambda* redshift bins
	 *
	 *  @param shot_noise_normalization the shot noise normalization. Shot noise il multiplyed for this value before summing with \f[C_\ell \f]
	 *
	 */
	void set_data_model_scaling_relation_cosmology_MF_theory (const double z_abs_err, const double proxy_rel_err, const std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const double redshift, const std::vector<double> cluster_redshift, const std::vector<double> cluster_mass_proxy, const double redshift_pivot, const double proxy_pivot, double mass_pivot, const double log_base, const std::string method_Pk="CAMB", const bool NL=true, const double k_min=1.e-4, const double k_max=100., const int norm=-1,  const double area_degrees=cbl::par::defaultDouble, const bool limber_approximation=true, const double ell_limber=20., const double z_min=par::defaultDouble, const double z_max=par::defaultDouble, const double proxy_min=par::defaultDouble, const double proxy_max=par::defaultDouble, std::vector<double> ll={}, std::vector<std::vector<double>> mixing_matrix={}, const bool redshift_space=false, const std::string interpType="Spline", const std::string output_dir=par::defaultString, const std::string output_root="test", const double prec=1.e-3, const std::string file_par=par::defaultString, const double Delta=200., const bool isDelta_critical=true, const std::string model_bias="Tinker", const std::string model_MF = "Tinker", const std::vector<double> SF_weights={1.}, double shot_noise_normalization=1.);

	/**
	 *  @brief set the model to fit the covariance of the angular power spectrum
	 *
	 *  The model is the following:
	 *
	 *  \f[ \sigma_{C_\ell} = \frac{2}{(2\ell+1)f_{sky} \Delta \ell} \left(\beta C_\ell + \frac{1+\alpha}{\overline{n}} \right) \f]
	 *
	 *  where the angular power spectrum model is computed with set_model_bias_MF_theory.
	 *
	 *  The parameter of the model are:
	 *
	 *  - N cosmological parameters
	 *  - nuisance extra shot noise term (used only in set_model_bias_MF_theory, keep it fixed to zero) 
	 *  - \f$\beta\f$ halo bias inaccuracy
	 *  - \f$\alpha\f$ deviations from poissonian shot noise
	 *  - \f$b_{eff}\f$ the effective bias is a derived parameter (from cosmology)
	 *
	 *  @param cosmo_param list of cosmological parameter for the model
	 *
	 *  @param cosmo_prior prior for the cosmological parameters
	 *
	 *  @param shot_noise_prior prior for the shot noise
	 *  
	 *  @param beta_prior prior for the beta parameter
	 *
	 *  @param alpha_prior prior for the alpha parameter
 	 *
	 */
	void set_model_covariance_bias_MF_theory (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution shot_noise_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution alpha_prior);
	
        /**
         *  @brief Set the model to fit the angular power spectrum
         *
         *  The model is the following:
         *
         *  \f[ C_{ij}^{\rm gg}(\ell) =
              \int {\rm d}z
              \frac{W_{i}^{\rm g}(z)W_{j}^{\rm g}(z)}
              {H(z)r^2(z)}
              P^{\rm{photo}}_{\rm gg}
              \left[ \frac{\ell +1/2}{r(z)}, z \right] \f]
         *
         *  the model has n parameters:
         *    - N cosmological parameters
         *    
         *  the angular power spectrum is computed
         *  using the input cosmological parameters
         *
         *
         *  @param cosmo_param list of cosmological parameter for the model
         *
         *  @param cosmo_prior prior for the cosmological parameters
         */
   	void set_model_Cij_gg (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior);
   	
        /**
         *  @brief Set the model to fit the angular power spectrum
         *
         *  The model is the following:
         *
         *\f[ C_{ij}^{\rm GG}(\ell) = C_{ij}^{\rm gg}(\ell)
             + C_{ij}^{\rm g\mu}(\ell) + C_{ij}^{\mu\mu}(\ell) \f]
         *
         *  \f[ C_{ij}^{\rm gg}(\ell) = c
              \int {\rm d}z
              \frac{W_{i}^{\rm g}(z)W_{j}^{\rm g}(z)}
              {H(z)r^2(z)}
              P_{\rm gg}
              \left[ \frac{\ell +1/2}{r(z)}, z \right] \f]
         *
         *  \f[ C_{ij}^{\rm g\mu}(\ell) = c 
              \int {\rm d}z
              \frac{W_{i}^{\rm \mu}(z)W_{j}^{\rm g}(z) +
                    W_{i}^{\rm g}(z)W_{j}^{\rm \mu}(z)}
              {H(z)r^2(z)}
              P_{\rm g\delta}
              \left[\frac{\ell +1/2}{r(z)} , z \right] \f]
         *
         *  \f[  C_{ij}^{\rm \mu\mu}(\ell) = c
            \int {\rm d}z
            \frac{W_{i}^{\rm \mu}(z)W_{j}^{\rm \mu}(z)}
            {H(z)r^2(z)}
            P_{\delta \delta}
            \left[ \frac{\ell +1/2}{r(z)}, z \right] \f]
         *
         *  the model has n parameters:
         *    - N cosmological parameters
         *    
         *  the angular power spectrum is computed
         *  using the input cosmological parameters
         *
         *
         *  @param cosmo_param list of cosmological parameter for the model
         *
         *  @param cosmo_prior prior for the cosmological parameters
         *
         */
   	void set_model_Cij_GG (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior);
	
	/**
	 *  @brief set the parameters for the clustering window functions
	 *
	 *  @param nz_all vector that contains all the redshift distribution used to calculated the angular power spectrum
	 *
 	 *  @param z_bias vector that contains the redshift in order to interpolated the bias
 	 *
	 *  @param bias_all vector that contains the fiducial bias for each bin
	 *
	 *  @param mag_bias_all vector that contains magnification bias for each bin
	 *
	 */	
	void set_nz (const std::vector<std::shared_ptr<glob::FuncGrid>> nz_all, const std::vector<double> z_bias, const std::vector<double> bias_all, const std::vector<double> mag_bias_all);	
	
	///@}	
      };
    }
  }
}

#endif
