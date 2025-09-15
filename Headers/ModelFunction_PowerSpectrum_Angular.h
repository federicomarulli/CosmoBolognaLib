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
 *  @file Headers/ModelFunction_PowerSpectrum_Angular.h
 *
 *  @brief Global functions to model the angular power spectrum
 *
 *  This file contains all the prototypes of the functions used
 *  to model the angular power spectrum
 *  
 *  @author Federico Marulli, Massimiliano Romanello
 *
 *  @author federico.marulli3@unibo.it, massimilia.romanell2@unibo.it
 */

#ifndef __MODFUNCPOWSPECTRUMANG__
#define __MODFUNCPOWSPECTRUMANG__

#include "Modelling_MassObservableRelation.h"


// ============================================================================


namespace cbl {

  namespace modelling {

    /**
     *  @brief The namespace of the <B> angular power spectrum
     *  modelling </B>
     *  
     *  The \e modelling::angularpk namespace contains all the
     *  functions and classes to model the angular power spectrum
     */
    namespace angularpk {

      using ProxyFunction = std::function<double(double, double, std::shared_ptr<void>)>;
      
      struct STR_data_model {

	/// the cosmological model
	std::shared_ptr<cosmology::Cosmology> cosmology;

	/// cosmological parameters
	std::vector<std::string> Cpar;

	/// minimum redshift
	double z_min;

	/// maximum redshift
	double z_max;

	/// minimum redshift of the second redshift bin
	double z_min_bin2;

	/// maximum redshift of the second redshift bin
	double z_max_bin2;

	/// method used to compute the power spectrum (i.e. the Boltzmann solver)
	std::string method_Pk;

	/// linear or non-linear power spectrum
	bool NL;

	/// norm 0 \f$\rightarrow\f$ don't normalise the power spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum; -1 \f$\rightarrow\f$ normalise only if sigma8 is set
	int norm;

	/// minimum wave vector module
	double k_min;

	/// maximum wave vector module
	double k_max;

	/// the survey area (in square degrees)
	double area_degrees;

	/// fraction of the sky covered by the survey
	double fsky;

	/// use limber approximation
	bool limber_approximation;
	
	/// l upon which the Limber approximation is used
	double ell_limber;

	/// the mixing matrix
	std::vector<std::vector<double>> mixing_matrix;
	
	/// the mixing matrix used for the covariance matrix
	std::vector<std::vector<double>> mixing_matrix_cov;

	/// multipoles of the mixing matrix
	std::vector<double> ll;

	/// real or redshift space
	bool redshift_space;

	/// vector containing the parameters of dN/dz in a given photometric bin
	std::vector<double> dN_par;

	/// vector containing the parameters of dN/dz in the second photometric bin
	std::vector<double> dN_par_bin2;

	/// the distribution type of dN/dz
	std::string distribution_type;

	/// vector of cluster redshift
	std::vector<double> cluster_redshift;

	/// vector of cluster masses
	std::vector<double> cluster_mass;

	/// vector of cluster mass proxies
	std::vector<double> cluster_mass_proxy;

	/// redshift pivot in the scaling relation
	double redshift_pivot;

	/// proxy pivot in the scaling relation
	double proxy_pivot;

	/// mass pivot in the scaling relation
	double mass_pivot;

	/// base of the logarithms in the scaling relation
	double log_base;

	/// the output_dir directory where the output of external codes are written
	std::string output_dir;

	/// output_root of the parameter file 
	std::string output_root;

	/// prec accuracy of the GSL integration
	double prec;

	/// name of the parameter file
	std::string file_par;

	/// the overdensity
	double Delta;

	/// the input overdensity 
	double Delta_input;

	/// control if the overdensity is defined with respect to the critical density
	bool isDelta_critical;

	/// the model of the bias
	std::string model_bias;

	/// the model of mass function
	std::string model_MF;

	/// the minimum mass for the mass function
	double Mass_min;

	/// the minimum mass for the mass function
	double Mass_max;

	/// the number of mass steps
	int Mass_step;

	/// sigma8 at redshift z
	double sigma8_z;

	/// the scaling relation
	std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation;

	/// absolute photometric error
	double z_abs_err = -1;

	/// relatiove proxy error
	double proxy_rel_err = -1;

	/// interpolation type
	std::string interpType;

	/// the redshift evolution function in the scaling relation
	ProxyFunction fz;

	/// weights for completeness and purity of number counts in lambda* redshift bins 
	std::vector<double> SF_weights;

	/// minimum proxy in a given proxy bin
	double proxy_min;

	/// maximum proxy in a given proxy bin
	double proxy_max;
	
	/// the shot noise is rescaled with this factor
	double shot_noise_normalization; 
		  
	/// compute the bias from as an average of the clusters bias from catalogue
	bool bias_from_sample;
	
	/// pointer to a vector of FuncGrid class,
	std::vector<std::shared_ptr<glob::FuncGrid>> nz_all;

	/// the vector that contains the redshift to interpolated the bias
	std::vector<double> z_bias;
	
	/// the vector that contains the bias for each bin
	std::vector<double> bias_all;	
	
	/// the vector that contains the magnification bias for each bin
	std::vector<double> mag_bias_all;
	
	STR_data_model () = default;
      };


      /**
       *  @brief the angular power spectrum convolved with the mixing matrix 
       *
       *  the function computes:
       *
       *  \f[C_l^{mixed}=\sum_{l'} R_{ll'}C_{l'}\f]
       *
       *  @param l_mixing the vector of multipoles of the mixing matrix
       *
       *  @param mixing_matrix the mixing matrix
       *
       *  @param l the vector of multipoles of the power spectrum
       *
       *  @param Cl the angular power spectrum
       *
       *  @param fsky the fraction of sky covered by the survey
       *
       *  @param limber_approximation true->use the limber
       *  approximation (mixing starts from l=0), false->don't use the
       *  limber approximation (mixing starts from l=2)
       *
       *  @return the angular power spectrum convolved with the mixing matrix
       */
      std::vector<double> Cl_mixed (std::vector<double> l_mixing, std::vector<std::vector<double>> mixing_matrix, std::vector<double> l, std::vector<double> Cl, double fsky, bool limber_approximation);

      /**
       *  @brief compute the redshift selection function at given z
       *
       *  @param redshift redshift
       *
       *  @param parameter vector containing the fit parameter of the
       *  corresponding redshift bin
       *
       *  @param distribution_type the shape of the dN/dz distribution.
       *  Can be "polynomial", for a fit of the photometric redshift
       *  distribution, "gaussian" or "skewnormal", in the correct
       *  case of convolution with P(z true|z_phot)
       *
       *  @param redshift_min minimum redshift of photometric bin (for
       *  the normalization of polynomial distribution)
       *
       *  @param redshift_max maximum redshift of photometric bin (for
       *  the normalization of polynomial distribution)
       *
       *  @return the redshift selection function at given z
       */
      double redshift_selection_function (double redshift, std::vector<double> parameter, std::string distribution_type, double redshift_min, double redshift_max);
      
      /**
       *  @brief the integrand function for
       *  the calculus of angular power spectrum 
       *
       *  the function computes:
       *
       *  \f[C_l = \frac{1}{N^2}\int_{0}^{\infty}
       *  \frac{dN}{dz}\frac{dN}{dz} P_{mat}
       *  \left(\frac{l+1/2}{r(z)}\right)
       *  \frac{H(z)}{c}\frac{1}{r^2(z)} dz\f]
       *
       *  @param l the multipole at which angular power spectrum is
       *  computed
       *
       *  @param kk the k vector of the interpolation grid
       *
       *  @param pk0_interp the interpolated power spectrum at z=0
       *
       *  @param comoving_distance_interp the interpolated comoving distance
       *
       *  @param interp_DN the interpolated normalised amplitude of
       *  the growing mode
       *
       *  @param linear_growth_rate_interp the interpolated linear growth rate
       *
       *  @param redshift_counts_interp the interpolated dN/dz distribution
       *
       *  @param bias the bias
       *
       *  @param redshift_counts_interp_bin2 the interpolated dN/dz
       *  distribution of the second redshift bin
       *
       *  @param bias_bin2 the bias of the second redshift bin
       *
       *  @param cross true->compute cross correlation, false->compute
       *  autocorrelation
       *
       *  @param parameter the list of parameters, containing the
       *  cosmological parameters
       *
       *  @param input pointer to the structure that contains the
       *  power spectrum angular data model
       *
       *  @return the integral of the angular power spectrum in the
       *  Limber approximation
       */
      double integral_cl (double l, std::vector<double> kk, cbl::glob::FuncGrid pk0_interp, cbl::glob::FuncGrid comoving_distance_interp, cbl::glob::FuncGrid interp_DN, cbl::glob::FuncGrid linear_growth_rate_interp, cbl::glob::FuncGrid redshift_counts_interp, const double bias, cbl::glob::FuncGrid redshift_counts_interp_bin2, const double bias_bin2, bool cross, std::vector<double> parameter, std::shared_ptr<void> input);

      /**
       *  @brief the model for the angular power spectrum
       *
       *  the model is the following:
       *
       *  \f[ C_l = \frac{b^2}{N^2}\int_{0}^{\infty}
       *  \frac{dN}{dz}\frac{dN}{dz} P_{mat}
       *  \left(\frac{l+1/2}{r(z)}\right)
       *  \frac{H(z)}{c}\frac{1}{r^2(z)} dz\f]
       *
       *  the model has n+3 parameters:
       *    - \f$ n cosmological parameters, including \Omega_m, \Omega_b and \sigma_8 \f$
       *    - \f$ b \f$
       *    - \f$ offset \f$
       *    - \f$ slope \f$
       *    
       *  the angular power spectrum is computed using the input
       *  cosmological parameters
       *
       *  @param l the vector of multipoles at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  power spectrum angular data model
       *
       *  @param parameter 1D vector containing the linear bias and
       *  the offset and the slope of the normalized dN/dz
       *  distribution
       *  
       *  @return the angular power spectrum in the Limber
       *  approximation
       */
      std::vector<double> Cl (const std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       *  @brief the angular power spectrum, with mass/intrinsic richness scaling relation
       *
       *  The function computes the angular power spectrum 
       *
       *  with dN/dz measured directly from catalogue
       
       *  @param l the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector containing the linear bias
       *
       *  @return the angular power spectrum computed with mass/intrinsic richness scaling relation.
       */
      std::vector<double> Cl_scaling_relation_cosmology (const std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       *  @brief the angular power spectrum, with mass/intrinsic richness scaling relation, using a theoretical mass function to compute dN/dz
       *
       *  The function computes the angular power spectrum 
       *
       *  @param l the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector containing the linear bias
       *
       *  @return the angular power spectrum computed with mass/intrinsic richness scaling relation.
       */
      std::vector<double> Cl_scaling_relation_cosmology_MF_theory (const std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       *  @brief the angular power spectrum
       *
       *  The function computes the angular power spectrum, using a theoretical mass function to compute dN/dz
       *
       *  @param l the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the angular power spectrum
       *
       *  @param parameter 1D vector containing the cosmological parameters and the linear bias
       *
       *  @return the angular power spectrum
       */
      std::vector<double> Cl_MF_theory (const std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       *  @brief the angular power spectrum
       *
       *  The function computes the angular power spectrum, using a
       *  theoretical mass function to compute dN/dz and theoretical
       *  bias
       *
       *  @param l the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological parameters used to compute the angular power
       *  spectrum
       *
       *  @param parameter 1D vector containing the cosmological
       *  parameters
       *
       *  @return the angular power spectrum
       */
      std::vector<double> Cl_bias_MF_theory (const std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       *  @brief the angular power spectrum
       *
       *  The function computes the angular power spectrum, using a theoretical mass function to compute dN/dz and theoretical bias.
       *  It includes non-linear correction at BAO
       *
       *  @param l the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the angular power spectrum
       *
       *  @param parameter 1D vector containing the cosmological parameters 
       *
       *  @return the angular power spectrum
       */
      std::vector<double> Cl_bias_MF_theory_BAO (const std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

       /**
       * @brief compute the redshift distribution as function of mass
       * proxy
       *
       * @param redshift redshift
       *
       * @param alpha alpha
       *
       * @param beta beta
       *
       * @param gamma gamma
       *
       * @param scatter0 scatter0
       *
       * @param scatterM scatterM
       *
       * @param scatterM_exp scatterM_exp
       *
       * @param scatterz scatterz
       *
       * @param scatterz_exp scatterz_exp
       *
       * @param z_abs_err abs. error on z
       *
       * @param proxy_rel_err rel. error on proxy
       *
       * @param Plambda_a Plambda_a
       *
       * @param Plambda_b Plambda_b
       *
       * @param Plambda_c Plambda_c
       *
       * @param fz redshift evolution function in the scaling relation
       *
       * @param redshift_min minimum redshift
       *
       * @param redshift_max maximum redshift
       *
       * @param proxy_min minimum mass proxy
       *
       * @param proxy_max maximum mass proxy
       *
       * @param cosmology the cosmology 
       *
       * @param Area the area in degrees
       *
       * @param model_MF author(s) who proposed the mass function;
       * valid authors are: PS (Press & Schechter), ST (Sheth &
       * Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       * al. 2006), Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH
       * (halo MF by Shen et al. 2006), ShenF (filament MF by Shen et
       * al. 2006), ShenS (sheet MF by Shen et al. 2006), Tinker
       * (Tinker et al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF
       * (FoF MF by Angulo et al. 2012), Angulo_Sub (SUBFIND MF by
       * Angulo et al. 2012), Watson_FOF (FoF MF by Watson et
       * al. 2012), Watson_SOH (Spherical Overdensity halo MF by
       * Watson et al. 2012), Manera (Manera et al. 2010),
       * Bhattacharya (Bhattacharya et al. 2011), Courtin (Courtin et
       * al. 2010), Peacock (by Peacock at al. 2007)
       *
       * @param Delta \f$\Delta\f$, the overdensity
       *
       * @param isDelta_critical \f$\rightarrow\f$ \f$\Delta\f$ is the
       * overdensity defined with respect to the critical density
       *
       * @param interp_sigmaM interpolating function of
       * \f$\sigma(M)\f$
       *
       * @param interp_DlnsigmaM interpolating function of
       * \f$\mathrm{d} \ln(\sigma(M)) / \mathrm{d} M\f$
       *
       * @param interp_DN interpolated amplitude of the growing mode
       *
       * @param proxy_pivot mass proxy pivot in the scaling relation
       *
       * @param z_pivot redshift pivot in the scaling relation
       *
       * @param mass_pivot mass pivot in the scaling relation
       *
       * @param log_base logarithmic base used in the scaling relation
       *
       * @param weights_interp weight derived from the selection
       * function, interpolated to include different richness bins
       *
       * @return values of redshift distribution as a function of mass
       * proxy
       */
      double redshift_distribution_proxy (const double redshift, const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const double z_abs_err, const double proxy_rel_err, const double Plambda_a, const double Plambda_b, const double Plambda_c, ProxyFunction fz, const double redshift_min, const double redshift_max, const double proxy_min, const double proxy_max, const std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const double Area, const std::string model_MF, const double Delta, const bool isDelta_critical, const cbl::glob::FuncGrid interp_sigmaM, const cbl::glob::FuncGrid interp_DlnsigmaM, const cbl::glob::FuncGrid interp_DN, const double proxy_pivot, const double z_pivot, const double mass_pivot, const double log_base, cbl::glob::FuncGrid weights_interp);

      /**
       *  @brief the set of multipoles for mixed power spectrum
       *
       *  The function adapts the initial set of multipoles (i.e. give
       *  by the input dataset) for the mixing matrix
       *  convolution. Indeed, the mixed power spectum, in partial sky
       *  survey, is obtained by convolving the full set of
       *  multipoles. The new set of multipoles:
       *
       *  -starts from zero
       *
       * -have few more multipoles at the end, to avoid the loss of
           power due to the truncation of mixing matrix
       *
       *  @param l the scale at which the model is computed
       *
       *  @param bin_width the band-width \f$ Delta_l \f$ in which the
       *  power spectum is averaged
       *
       *  @param limber_approximation true->use the limber
       *  approximation (ell will start from l=0), false->don't use
       *  the limber approximation (ell will start from l=2)
       *
       *  @return the new set of multipoles
       */
      std::vector<double> ell_for_mixing (std::vector<double> l, double bin_width, bool limber_approximation);

      /**
       * @brief limits of the true redshift interval, given the
       * photometric bin width and the photometric error
       *
       * @param zmin zmin
       *
       * @param zmax zmax
       *
       * @param z_abs_err the absolute redshift error, alone or
       * multiplied for (1+z_mean), with z_mean is the average
       * of the photometric redshift bin
       *
       * @return limits of the true redshift interval
       *
       */
      std::vector<double> z_limits (const double zmin, const double zmax, const double z_abs_err);

      /**
       * @brief compute the bias averaged in the redshift bin,
       * using the mass-intrinsic richness scaling relation
       *
       * @param alpha alpha
       *
       * @param beta beta
       *
       * @param gamma gamma
       *
       * @param scatter0 scatter0
       *
       * @param scatterM scatterM
       *
       * @param scatterM_exp scatterM_exp
       *
       * @param scatterz scatterz
       *
       * @param scatterz_exp scatterz_exp
       *
       * @param log_base logarithmic base used in the scaling relation
       *
       * @param mass_pivot mass pivot in the scaling relation
       *
       * @param proxy_pivot mass proxy pivot in the scaling relation
       *
       * @param redshift_pivot redshift pivot in the scaling relation
       *
       * @param scalRel_pars scaling relation parameters
       *
       * @param interp_DN interpolated amplitude of the growing mode
       *
       * @param inputs pointer to the structure that contains the
       * power spectrum angular data model
       *
       * @return values of redshift distribution as a function of redshift
       * and mass proxy
       */
      double bias_from_scaling_relation (const double alpha, const double beta, const double gamma,  const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, double log_base, double mass_pivot, double proxy_pivot, double redshift_pivot, std::vector<double> scalRel_pars, cbl::glob::FuncGrid interp_DN, std::shared_ptr<void> inputs);

      /**
       *  @brief analytic covariance matrix for the angular power spectrum
       *
       *  The function computes the analytic covariance matrix for
       *  angular power spectrum, using a theoretical mass function to
       *  compute dN/dz and theoretical bias
       *
       *  @param l the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the angular power spectrum
       *
       *  @param parameter 1D vector containing the model parameters 
       *
       *  @return the covariance matrix for the angular power spectrum
       */
      std::vector<double> covariance_Cl_bias_MF_theory (const std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       * @brief interpolated bias function
       *
       * @param zz the redshift
       *
       * @param z_center the vector that contains the midpoint values for each bin
       *
       * @param bb the fiducial values for each z_center
       *
       * @return the bias at given redshift
       *
       */      
      double bias_inter(const double zz, const std::vector<double> z_center, const std::vector<double> bb );
      
      
      /**
       * @brief compute the clustering window function for the density fluctuations
       *
       *  \f[ W_i^{g,den}(z) = \frac{H(z)}{c} n_{i}(z) \f]
       *
       * @param zz the redshift
       *
       * @param cosmo the cosmology
       *
       * @param nz the interpolated redshift distribution
       *
       * @return the clustering window function for the density fluctuations at given redshift
       *
       */ 
      double window_function_den (const double zz, const std::shared_ptr<cbl::cosmology::Cosmology> cosmo, const std::shared_ptr<cbl::glob::FuncGrid> nz);
      
      
      /**
       * @brief compute the clustering window function for RSD effects
       *
       *  \f[ W_i^{\rm{g,RSD}}(\ell,z) = \frac{1}{c} \sum_{m=-1}^{+1}
       *  L_m(\ell) \left[H(z_{\rm m})f(z_{\rm m}) \frac{n_{i}(z_{\rm
       *  m})}{b_i(z_{\rm m})} \right] \f]
       *
       *   where \f[z_m \f] gives the redshift corrisponding to the
       *   radial comoving distance r
       *
       *  \f[ z_{\rm m} = z\left(\frac{2\ell+1+4m}{2l+1} r\right) \f]
       *  and  \f[ L_m \f]  is defined as 
       *   \f[ L_{-1} = -\frac{\ell(\ell - 1)}
            {(2 \ell -1)\sqrt{(2 \ell - 3)(2 \ell + 1)}} \f]
       * \f[ L_{0} = \frac{2 \ell^2 + 2 \ell -1}
            {(2 \ell - 1)(2 \ell + 3)} \f]
       * \f[ L_{+1} = -\frac{(\ell + 1)(\ell + 2)}
            {(2\ell + 3)\sqrt{(2\ell + 1)(2 \ell + 5)}} \f]
       *
       * @param zz the redshift
       *
       * @param l the multipole order
       *
       * @param cosmo the cosmology
       *
       * @param nz the the interpolated redshift distribution
       *
       * @param z_r the the interpolated function of the redshift in
       * function of the comoving distance
       *
       * @param bias the interpolated bias
       *
       * @return the clustering window function for the RSD effects at a given redshift
       *
       */ 
      double window_function_rsd (const double zz, const double l, const std::shared_ptr<cbl::cosmology::Cosmology> cosmo, const std::shared_ptr<cbl::glob::FuncGrid> nz, cbl::glob::FuncGrid z_r, cbl::glob::FuncGrid bias);
      
      /**
       * @brief compute the clustering window function for the lensing magnification
       *
       *  \f[ W_{i}^{\mu}(\ell, z) =
            \frac{3}{2}\left ( \frac{H_0}{c}\right )^2
            \frac{\ell (\ell+1)}{(\ell+1/2)^2}
            \Omega_{{\rm m},0} (1 + z)r(z)
            \int_{z}^{z_{max}}{{\rm d}z^{\prime} b_{mag}(z^\prime) n_{i}(z^{\prime})
            \frac{r(z^{\prime}) - r(z)}
            {r(z^{\prime})}} \f]
       *
       * @param zz the redshift
       *
       * @param l the multipole order
       *
       * @param z_max the z_max
       *
       * @param mag_bias the interpolated function of the fiducial
       * value of the magnification bias
       *
       * @param cosmo the cosmology
       *
       * @param nz the the interpolated redshift distribution
       *
       * @param DC the interpolated function of the comoving radial distance 
       *
       * @return the clustering window function for the magnification
       * effects at a given redshift
       *
       */ 
      double window_function_mag (const double zz, const double l, const double z_max, cbl::glob::FuncGrid mag_bias, const std::shared_ptr<cbl::cosmology::Cosmology> cosmo, const std::shared_ptr<cbl::glob::FuncGrid> nz, cbl::glob::FuncGrid DC );
      
      /**
       * @brief compute the angular power spectra in harmonic space 
       *
       *  \f[ C_{ij}(\ell) = c
            \int {\rm d}z
            \frac{W_{i}(z)W_{j}(z)}
            {H(z)r^2(z)} f^2(z)
            P_{\rm \delta\delta}
            \left[ \frac{\ell +1/2}{r(z)}, z \right] \f]
       *
       * @param l the multipole order
       *
       * @param z_min the z_min
       *
       * @param z_max the z_max  
       *
       * @param W_i the interpolated weight function for the first redshift distribution
       *
       * @param W_j the interpolated weight function for the second redshift distribution
       *
       * @param HH the interpolated function of Hubble rate
       *
       * @param DC the interpolated function of the comoving radial distance 
       *
       * @param Pk the 2D interpolated function of Pk in function of (z,k)
       *
       * @param func the interpolated function to multiple the power spectrum (like the bias)
       *
       * @return the angular power spectra in harmonic space 
       *
       */      
      double C_ij_l (const double l,  const double z_min, const double z_max, cbl::glob::FuncGrid W_i,  cbl::glob::FuncGrid W_j, cbl::glob::FuncGrid HH,cbl::glob::FuncGrid DC, cbl::glob::FuncGrid2D Pk,cbl::glob::FuncGrid func );   

      /**
       * @brief compute the angular power spectra in harmonic space 
       *
       *\f[ C_{ij}^{\rm XY}(\ell) = C_{ij}^{\rm xx}(\ell)
            + C_{ij}^{\rm xy}(\ell) + C_{ij}^{yy}(\ell) \f]
       *
       *  \f[ C_{ij}^{ xx}(\ell) = c
            \int {\rm d}z
            \frac{W_{i}^{x}(z)W_{j}^{x}(z)}
            {H(z)r^2(z)} f^2(z)
            P_{\rm \delta \delta }
            \left[ \frac{\ell +1/2}{r(z)}, z \right] \f]
       *
       *  \f[ C_{ij}^{xy}(\ell) = c 
            \int {\rm d}z
            \frac{W_{i}^{y}(z)W_{j}^{x}(z) +
                  W_{i}^{x}(z)W_{j}^{y}(z)}
            {H(z)r^2(z)} f(z)
            P_{\rm \delta \delta}
            \left[\frac{\ell +1/2}{r(z)} , z \right] \f]
       *
       *  \f[  C_{ij}^{yy}(\ell) = c
            \int {\rm d}z
            \frac{W_{i}^{y}(z)W_{j}^{y}(z)}
            {H(z)r^2(z)}
            P_{\delta \delta}
            \left[ \frac{\ell +1/2}{r(z)}, z \right] \f]
       *
       * @param l the multipole order
       *
       * @param z_min the z_min
       *
       * @param z_max the z_max  
       *
       * @param Wx_i the interpolated window function of kind X for
       * the first redshift distribution
       *
       * @param Wx_j the interpolated window function of kind X for
       * the second redshift distribution
       *
       * @param Wy_i the interpolated window function of kind Y for
       * the first redshift distribution
       *
       * @param Wy_j the interpolated window function of kind Y for
       * the second redshift distribution
       *
       * @param HH the interpolated function of Hubble rate
       *
       * @param DC the interpolated function of the comoving radial distance 
       *
       * @param Pk the 2D interpolated function of Pk in function of (z,k)
       *
       * @param func the interpolated function to multiple the power spectrum (like the bias)
       *
       * @return the angular power spectra in harmonic space 
       *
       */ 
      double C_ij_xy (const double l,  const double z_min, const double z_max, cbl::glob::FuncGrid Wx_i,  cbl::glob::FuncGrid Wx_j,cbl::glob::FuncGrid Wy_i,cbl::glob::FuncGrid Wy_j, cbl::glob::FuncGrid HH,cbl::glob::FuncGrid DC, cbl::glob::FuncGrid2D Pk,cbl::glob::FuncGrid func);
      
      /**
       *  @brief the model for the angular power spectrum
       *
       *  the model is the following:
       *
       * \f[ C_{ij}^{\rm gg}(\ell) = c \int {\rm d}z \frac{W_{i}^{\rm
       *     g}(z)W_{j}^{\rm g}(z)} {H(z)r^2(z)} b^2(z) P_{\rm \delta
       *     \delta} \left[ \frac{\ell +1/2}{r(z)}, z \right] \f]
       *
       *  the model has n parameters:
       *    -  N cosmological parameters
       *    
       *  the angular power spectrum is computed using the input
       *  cosmological parameters
       *
       *  @param l the vector of multipoles at which the model is
       *  computed
       *
       *  @param inputs pointer to the structure that contains the
       *  power spectrum angular data model
       *
       *  @param parameter 1D vector containing the linear bias and
       *  the offset and the slope of the normalized dN/dz
       *  distribution
       *  
       *  @return the angular power spectrum
       *
       */
      std::vector<double> C_ij_gg (const std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter); 
      
      /**
       *  @brief the model for the angular power spectrum
       *
       *  the model is the following:
       *
       *  \f[ C_{ij}^{\rm GG}(\ell) = C_{ij}^{\rm gg}(\ell) +
       *  C_{ij}^{\rm g\mu}(\ell) + C_{ij}^{\mu\mu}(\ell) \f]
       *
       *  \f[ C_{ij}^{\rm gg}(\ell) = c \int {\rm d}z \frac{W_{i}^{\rm
       *     g}(z)W_{j}^{\rm g}(z)} {H(z)r^2(z)} b^2(z) P_{\rm \delta
       *     \delta} \left[ \frac{\ell +1/2}{r(z)}, z \right] \f]
       *
       *  \f[ C_{ij}^{\rm g\mu}(\ell) = c \int {\rm d}z
       *     \frac{W_{i}^{\rm \mu}(z)W_{j}^{\rm g}(z) + W_{i}^{\rm
       *     g}(z)W_{j}^{\rm \mu}(z)} {H(z)r^2(z)} b(z) P_{\rm \delta
       *     \delta} \left[\frac{\ell +1/2}{r(z)} , z \right] \f]
       *
       *  \f[ C_{ij}^{\rm \mu\mu}(\ell) = c \int {\rm d}z
       *     \frac{W_{i}^{\rm \mu}(z)W_{j}^{\rm \mu}(z)} {H(z)r^2(z)}
       *     P_{\delta \delta} \left[ \frac{\ell +1/2}{r(z)}, z
       *     \right] \f]
       *
       *  the model has n parameters:
       *    - N cosmological parameters
       *    
       *  the angular power spectrum is computed using the input
       *  cosmological parameters
       *
       *  @param l the vector of multipoles at which the model is
       *  computed
       *
       *  @param inputs pointer to the structure that contains the
       *  power spectrum angular data model
       *
       *  @param parameter 1D vector containing the linear bias and
       *  the offset and the slope of the normalized dN/dz
       *  distribution
       *  
       *  @return the angular power spectrum 
       */
      std::vector<double> C_ij_GG (const std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter); 

    }
  }
}

#endif
