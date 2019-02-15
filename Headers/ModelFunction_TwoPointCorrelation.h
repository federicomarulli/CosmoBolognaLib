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
 *  @file Headers/ModelFunction_TwoPointCorrelation.h
 *
 *  @brief Global functions to model two-point correlation functions
 *  of any type
 *
 *  This file contains all the prototypes of the global functions used
 *  to model two-point correlation functions of any type
 *  
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODFUNCTWOP__
#define __MODFUNCTWOP__

#include "Cosmology.h"


// ============================================================================


namespace cbl {

  namespace modelling {

    namespace twopt {
      
      /**
       *  @struct STR_data_model
       *  @brief the structure STR_data_model
       *
       *  This structure contains the data used for statistical
       *  analyses of the two-point correlation function
       */
      struct STR_data_model {

	/// order of the polynomial that take 
	//systematic effects into account
	int poly_order;

	/// fiducial cosmology
	std::shared_ptr<cosmology::Cosmology> cosmology;

	/// test cosmology
	std::shared_ptr<cosmology::Cosmology> test_cosmology;

	/// redshift
	double redshift;

	/// method to compute the dark matter power spectrum
	std::string method_Pk;
      
	/// the output_dir directory where the output of external codes are written
	std::string output_dir;
      
	/// output root of the parameter file used to compute the dark matter power spectrum
	std::string output_root;

	/// false &rarr; linear power spectrum; true &rarr; non-linear power spectrum
	bool NL;

	/// sigmaNL damping of the wiggles in the linear power spectrum
	double sigmaNL;

	/// sigmaNL damping of the wiggles in the linear power spectrum, perpendicular direction
	double sigmaNL_perp;

	/// sigmaNL damping of the wiggles in the linear power spectrum, parallel direction
	double sigmaNL_par;

	/// 0 &rarr; don't normalize the power spectrum; 1 &rarr; normalize the power spectrum
	int norm;

	/// minimum wave vector module up to which the power spectrum is computed
	double k_min;

	/// maximum wave vector module up to which the power spectrum is computed
	double k_max;

	/// parameter \e a of Eq. 24 of Anderson et al. 2012
	double aa;

	/// 0 &rarr; FFTlog is used; 1 &rarr; the GSL libraries are used
	bool GSL;

	/// accuracy of the GSL integration
	double prec;

	/// name of the parameter file
	std::string file_par;

	/// pointer to a function of func_grid_GSL class, used to interpolate of the two-point correlation function
	std::shared_ptr<glob::FuncGrid> func_xi;
	
	/// barred &xi;(r) as pointer to an interpolation function
	std::shared_ptr<glob::FuncGrid> func_xi_;

	/// double-barred &xi;(r) as pointer to an interpolation function
	std::shared_ptr<glob::FuncGrid> func_xi__;
	
	/// upper limit of integration for the projected correlation function
	double pi_max;
      
	/// minimum separation up to which the binned dark matter correlation function is computed
	double r_min;

	/// maximum separation up to which the binned dark matter correlation function is computed
	double r_max;

	/// number of steps used to compute the binned dark matter correlation function
	int step;
	
	/// the linear growth rate at redshift z
	double linear_growth_rate_z;

	/// &sigma;<SUB>8</SUB> at redshift z
	double sigma8_z;

	/// (1+z)/HH(z)
	double var;

	/// cosmological parameters
	std::vector<cosmology::CosmologicalParameter> Cpar;

	/// FV 0 &rarr; exponential form for f(v); 1 &rarr; Gaussian form for f(v); where f(v) is the velocity distribution function
	int FV;

	/// 0 &rarr; linear bias; 1 &rarr; non-linear bias 
	int bias_nl;
       
	/// non-linear bias parameter
	double bA;
       
	/// 0 &rarr; standard; 1 &rarr; Chuang & Wang model
	int xiType;
       
	/// k<SUB>*</SUB> of the Chuang & Wang model
	double k_star;
       
	/// 0 &rarr; linear two-point correlation function; 1 &rarr; non-linear two-point correlation function
	int xiNL;
       
	/// v_min minimum velocity used in the convolution of the two-point correlation function
	double v_min;
       
	/// v_max maximum velocity used in the convolution of the two-point correlation function
	double v_max;
       
	/// number of steps used in the convolution of the two-point correlation function
	int step_v;
      
	/// index for pre-computed two-point correlation function
	int xi_real_index;

	///  0 &rarr; don't use the pole in the fit; 1 &rarr;  use the pole in the fit
	std::vector<bool> use_pole;

	/// number of (even) multipoles to decompose \f$P(k, \mu)\f$
	std::vector<int> dataset_order;

	/// number of (even) multipoles to decompose \f$P(k, \mu)\f$
	int nmultipoles;

	/// number of two-point correlation function wedges
	int nwedges;

	/// vector of wave vector modules
	std::vector<double> kk;

	/// vector of scales
	std::vector<double> rr;

	/// pointer to a function of FuncGrid class, used to interpolate the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk;

	/// pointer to a function of FuncGrid class, used to interpolate the no-wiggles linear power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_NW;

	/// pointer to a function of FuncGrid class, used to interpolate the power spectrum 1-loop correction
	std::shared_ptr<glob::FuncGrid> func_Pk1loop;

	/// pointer to a vector of FuncGrid objects, used to interpolate the multipoles
	std::vector<std::shared_ptr<glob::FuncGrid>> func_multipoles;

	/// pointer to a vector of FuncGrid objects, used to interpolate the wedges
	std::vector<std::shared_ptr<glob::FuncGrid>> func_wedges;

	/// the \f$P(k,\mu)\f$ model; the possible options are: dispersion_dewiggled, dispersion_modecoupling
	std::string Pk_mu_model;

	/// fiducial bias value
	double bias;

	/// fiducial \f$D_V\f$
	double DVfid;
	
	/// fiducial \f$D_A(z)\f$
	double DAfid;

	/// fiducial \f$H(z)\f$
	double HHfid;

	/// &Delta;: the overdensity, defined as the mean interior density relative to the background
	double Delta;

	/// isDelta_Vir
	bool isDelta_Vir;

	/// Pointer to normal random numbers generator
	std::shared_ptr<cbl::random::NormalRandomNumbers> gau_ran;

	/// cluster masses proxy
	std::shared_ptr<cbl::data::Data> cluster_mass_proxy;

	/// cluster masses proxy standard deviation
	std::vector<double> cluster_mass_proxy_error;	
	
	/// method to estimate the mass function
	std::string model_MF;

	/// method to estimate the bias
	std::string model_bias;

	/// meanType, either the mean bias or the pair mean bias
	std::string meanType;
					       
	///  pointer to a function of FuncGrid class, used to interpolate the \f$\sigma(M)\f$
	std::shared_ptr<glob::FuncGrid> func_sigma;

	/// function to interpolate the effective bias against one cosmological parameter
	std::function<double(const double)> cosmopar_bias_interp_1D;

	/// function to interpolate the effective bias against two cosmological parameters
	std::function<double(const double, const double)> cosmopar_bias_interp_2D;
		
	/// function to interpolate the selection function in mass, at the mean redshift
	std::shared_ptr<glob::FuncGrid> interp_SelectionFunction_cut;
	
	/// function to interpolate the selection function in mass and redshift
	std::shared_ptr<glob::FuncGrid2D> interp_SelectionFunction;

	/// minimum redshift
	double z_min; 

	/// maximum redshift
	double z_max;

	/// numeber of mass steps 
	int mass_step; 

	/// minimum mass 
	double Mass_min; 

	/// maximum mass
	double Mass_max; 

	/// vector containing the masses
	std::vector<double> mass;

	/// cosmology used to measure the cluster masses
	cosmology::Cosmology cosmology_mass;
	
	/// redshift_source vector containing the redshifts of the source galaxies, in case the cluster masses are estimated from weak lensing
	std::vector<double> redshift_source;
	
	/**
	 *  @brief default constructor
	 *  @return object of type STR_data_model
	 */
	STR_data_model () = default;
      };

      
      /**
       *  @brief the power spectrum as a function of k and \f$\mu\f$
       *
       *  this function computes the redshift-space BAO-damped power
       *  spectrum \f$P(k, \mu)\f$ (see e.g. Vargas-Magana et al. 2018
       *  https://arxiv.org/pdf/1610.03506.pdf) :
       *
       *  \f[ P(k, \mu) = \left(1+\beta\mu_p^2 \right)^2 \left(
       *  \frac{1}{1+\left(k_pf\Sigma_S\mu_p\right)^2} \right)^2
       *  P_{NL}(k_p) \f]
       *
       *  where
       *
       *  \f[ P_{NL}(k) = b^2 \left\{ \left[ P_{lin}(k) -
       *  P_{nw}(k) \right] e^{-k^2\Sigma_{NL}^2}+P_{nw}(k) \right\} , \f]
       *
       *  \f[ k_p = \frac{k}{\alpha_\perp} \sqrt{1+\mu^2 \left[ \left(
       *  \frac{\alpha_\parallel}{\alpha_\perp}\right)^{-2}-1 \right]}
       *  \, , \f]
       *
       *  \f[ \mu_p = \mu\frac{\alpha_\parallel}{\alpha_\perp}
       *  \frac{1}{\sqrt{1+\mu^2\left[\left(
       *  \frac{\alpha_\parallel}{\alpha_\perp}
       *  \right)^{-2}-1\right]}} \, , \f]
       *
       *  \f$\mu = k_{\parallel} / k\f$, \f$P_{lin}\f$, \f$P_{nw}\f$
       *  are the linear and the de-wiggled power spectra,
       *  respectively (Eisenstein et al. 1998), and \f$\beta = f/b\f$
       *  with \f$f\f$ the linear growth rate and \f$b\f$ the bias,
       *  and \f$\Sigma_S\f$ is the streaming scale that parameterises 
       *  the Fingers of God effect at small scales. 
       *  The BAO damping is parametrised via \f$\Sigma^2_{NL} = 0.5
       *  (1-\mu_p^2)\Sigma^2_{\perp}+\mu_p^2\Sigma^2_{\parallel} \f$,
       *  where \f$\Sigma_{\perp}\f$ and \f$\Sigma_{\parallel}\f$ are
       *  the damping term in the transverse and parallel directions
       *  to the line of sight, respectively.
       *
       *  @param kk the wave vector module
       *  
       *  @param mu the line of sight cosine
       *
       *  @param alpha_perp the shift transverse to the l.o.s.
       *
       *  @param alpha_par the shift parallel to the l.o.s.
       *
       *  @param sigmaNL_perp the damping in the direction transverse
       *  to the l.o.s.
       *
       *  @param sigmaNL_par the damping in the direction parallel to
       *  the l.o.s.
       *
       *  @param linear_growth_rate the linear growth rate
       *
       *  @param bias the linear bias
       *
       *  @param SigmaS streaming scale
       *
       *  @param Pk linear power spectrum interpolator
       *
       *  @param Pk_NW de-wiggled power spectrum interpolator
       *
       *  @return \f$P(k, \mu)\f$
       */
      double Pkmu_DeWiggled (const double kk, const double mu, const double alpha_perp, const double alpha_par, const double sigmaNL_perp, const double sigmaNL_par, const double linear_growth_rate, const double bias, const double SigmaS, const std::shared_ptr<cbl::glob::FuncGrid> Pk, const std::shared_ptr<cbl::glob::FuncGrid> Pk_NW);
      
      /**
       *  @brief the power spectrum as a function of k and \f$\mu\f$
       *
       *  this function computes the redshift-space power spectrum 
       *  \f$P(k, \mu)\f$  (see e.g. Sanchez et al. 2013
       *  https://arxiv.org/pdf/1312.4854.pdf) :
       *
       *  \f[ P(k, \mu) = \left(1+\beta\mu_p^2 \right)^2 \left(
       *  \frac{1}{1+\left(k_pf\sigma_v\mu_p\right)^2} \right)^2
       *  P_{NL}(k_p) \f]
       *
       *  where
       *
       *  \f[ P_{NL}(k) = b^2 \left\{ P_{L}(k)*e^{-(k\mu\sigma_v)^2}
       *  +A_{MC}P_{MC}(k) \right\} , \f]
       *
       *  \f[ k_p = \frac{k}{\alpha_\perp} \sqrt{1+\mu^2 \left[ \left(
       *  \frac{\alpha_\parallel}{\alpha_\perp}\right)^{-2}-1 \right]}
       *  \, , \f]
       *
       *  \f[ \mu_p = \mu\frac{\alpha_\parallel}{\alpha_\perp}
       *  \frac{1}{\sqrt{1+\mu^2\left[\left(
       *  \frac{\alpha_\parallel}{\alpha_\perp}
       *  \right)^{-2}-1\right]}} \, , \f]
       *
       *  \f$\mu = k_{\parallel} / k\f$, \f$P_{lin}\f$, \f$P_{MC}\f$
       *  are the linear and the linear and 1loop correction power
       *  spectra, and \f$\beta = f/b\f$
       *  with \f$f\f$ the linear growth rate and \f$b\f$ the bias,
       *  \f$\sigma_v\f$ is the streaming scale that parameterises 
       *  the Fingers of God effect at small scales and A_{MC} is the
       *  mode coupling bias. 
       *
       *  @param kk the wave vector module
       *  
       *  @param mu the line of sight cosine
       *
       *  @param alpha_perp the shift transverse to the l.o.s.
       *
       *  @param alpha_par the shift parallel to the l.o.s.
       *
       *  @param linear_growth_rate the linear growth rate
       *
       *  @param bias the linear bias
       *
       *  @param sigmav streaming scale
       *  
       *  @param AMC the mode coupling bias
       *
       *  @param PkLin linear power spectrum interpolator
       *
       *  @param PkMC the 1loop power spectrum correction
       *
       *  @return \f$P(k, \mu)\f$
       */
      double Pkmu_ModeCoupling (const double kk, const double mu, const double alpha_perp, const double alpha_par, const double linear_growth_rate, const double bias, const double sigmav, const double AMC, const std::shared_ptr<cbl::glob::FuncGrid> PkLin, const std::shared_ptr<cbl::glob::FuncGrid> PkMC);

      /**
       *  @brief the power spectrum as a function of k and \f$\mu\f$
       *
       *  this function computes the redshift-space power spectrum 
       *  \f$P(k, \mu)\f$  with one of the following models:
       *
       *  - the dispersion + de-wiggled model, implemented in
       *    cbl::modelling::twopt::Pkmu_DeWiggled
       *
       *  - the dispersion + mode-coupling model, implemented in
       *    cbl::modelling::twopt::Pkmu_DeWiggled
       *
       *  The above models may differ for both the redshit-space
       *  distortions and the non-linear power spectrum
       *  implementation, and may have different numbers of free
       *  parameters
       *
       *  @param kk the wave vector module
       *  
       *  @param mu the line of sight cosine
       *
       *  @param model the \f$P(k,\mu)\f$ model; the possible options
       *  are: dispersion_dewiggled, dispersion_modecoupling
       *
       *  @param parameter vector containing parameter values
       *
       *  @param pk_interp vector containing power spectrum
       *  interpolating functions
       *
       *  @return \f$P(k, \mu)\f$
       */
      double Pkmu (const double kk, const double mu, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp);

      /**
       *  @brief the multipole of order l of the power spectrum
       *
       *  The function computes the legendre polynomial expansion of
       *  the \f$P(k, \mu)\f$:
       *
       *  \f[ P_l(k) = \frac{2l+1}{2} \int_{-1}^{1} \mathrm{d}\mu P(k,
       *  \mu) L_l(\mu) \f]
       *
       *  where \f$l\f$ is the order of the expansion and
       *  \f$L_l(\mu)\f$ is the Legendere polynomial of order \f$l\f$;
       *  \f$P(k, \mu)\f$ is computed by
       *  cbl::modelling::twopt::Pkmu
       *
       *  @param kk the wave vector module
       *
       *  @param l the order of the expansion
       *
       *  @param model the \f$P(k,\mu)\f$ model; the possible options
       *  are: dispersion_dewiggled, dispersion_modecoupling
       *
       *  @param parameter vector containing parameter values
       *
       *  @param pk_interp vector containing power spectrum
       *  interpolating functions
       *
       *  @param prec the integral precision
       *
       *  @return the multipole expansion of \f$P(k, \mu)\f$ at given
       *  \f$k\f$
       */
      double Pk_l (const double kk, const int l, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5);

      /**
       *  @brief the multipole of order l of the power spectrum
       *
       *  The function computes the legendre polynomial expansion of
       *  the \f$P(k, \mu)\f$:
       *
       *  \f[ P_l(k) = \frac{2l+1}{2} \int_{-1}^{1} \mathrm{d}\mu P(k,
       *  \mu) L_l(\mu) \f]
       *
       *  where \f$l\f$ is the order of the expansion and
       *  \f$L_l(\mu)\f$ is the Legendere polynomial of order \f$l\f$;
       *  \f$P(k, \mu)\f$ is computed by
       *  cbl::modelling::twopt::Pkmu
       *
       *  @param kk the wave vector module vector
       *
       *  @param l the order of the expansion
       *
       *  @param model the \f$P(k,\mu)\f$ model; the possible options
       *  are: dispersion_dewiggled, dispersion_modecoupling
       *
       *  @param parameter vector containing parameter values
       *
       *  @param pk_interp vector containing power spectrum
       *  interpolating functions
       *
       *  @param prec the integral precision
       *
       *  @return the multipole expansion of \f$P(k, \mu)\f$ at given
       *  \f$k\f$
       */
      std::vector<double> Pk_l (const std::vector<double> kk, const int l, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5);

      /**
       *  @brief the interpolating function of multipole expansion of the
       *  two-point correlation function at a given order l
       *
       *  The function computes the multipoles of the two-point
       *  correlation function:
       *
       *  \f[ \xi_l(s) = \frac{i^l}{2\pi^2} \int \mathrm{d} k P_l(k)
       *  j_l(ks) \f]
       *
       *  where \f$j_l(ks)\f$ are the bessel functions, and
       *  \f$P_l(k)\f$ is computed by cbl::modelling::twopt::Pk_l
       *
       *  @param kk the wave vector module vector
       *
       *  @param l the order of the expansion
       *
       *  @param model the \f$P(k,\mu)\f$ model; the possible options
       *  are: dispersion_dewiggled, dispersion_modecoupling
       *
       *  @param parameter vector containing parameter values
       *
       *  @param pk_interp vector containing power spectrum
       *  interpolating functions
       *
       *  @param prec the integral precision
       *
       *  @return the interpolation function of the multipole
       *  expansion of two-point correlation function at a given order
       *  l
       */
       cbl::glob::FuncGrid Xil_interp (const std::vector<double> kk, const int l, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5);


      /**
       *  @brief the multipole of order l of the two-point correlation
       *  function
       *
       *  The function computes the multipoles of the two-point
       *  correlation function:
       *
       *  \f[ \xi_l(s) = \frac{i^l}{2\pi^2} \int \mathrm{d} k P_l(k)
       *  j_l(ks) \f]
       *
       *  where \f$j_l(ks)\f$ are the bessel functions, and
       *  \f$P_l(k)\f$ is computed by cbl::modelling::twopt::Pk_l
       *
       *  @param rr vector of scales to compute multipoles
       *  
       *  @param nmultipoles the number of (even) multipoles to
       *  compute
       *
       *  @param model the \f$P(k,\mu)\f$ model; the possible options
       *  are: dispersion_dewiggled, dispersion_modecoupling
       *
       *  @param parameter vector containing parameter values
       *
       *  @param pk_interp vector containing power spectrum
       *  interpolating functions
       *
       *  @param prec the integral precision
       *
       *  @return the multipole of order l of the two-point
       *  correlation function
       */
      std::vector<std::vector<double>> Xi_l (const std::vector<double> rr, const int nmultipoles, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5);

      /**
       *  @brief the multipole of order l of the two-point correlation
       *  function
       *
       *  The function computes the multipoles of the two-point
       *  correlation function:
       *
       *  \f[ \xi_l(s) = \frac{i^l}{2\pi^2} \int \mathrm{d} k P_l(k)
       *  j_l(ks) \f]
       *
       *  where \f$j_l(ks)\f$ are the bessel functions, and
       *  \f$P_l(k)\f$ is computed by cbl::modelling::twopt::Pk_l
       *
       *  @param rr vector of scales to compute multipoles
       *
       *  @param dataset_order vector that specify the multipole
       *  to be computed for each scale 
       *
       *  @param use_pole vector of booleans specifying if a given 
       *  multipole should be computed
       *
       *  @param model the \f$P(k,\mu)\f$ model; the possible options
       *  are: dispersion_dewiggled, dispersion_modecoupling
       *
       *  @param parameter vector containing parameter values
       *
       *  @param pk_interp vector containing power spectrum
       *  interpolating functions
       *
       *  @param prec the integral precision
       *
       *  @return the multipole of order l of the two-point
       *  correlation function
       */
      std::vector<double> Xi_l (const std::vector<double> rr, const std::vector<int> dataset_order, const std::vector<bool> use_pole, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5);

      /**
       *  @brief the wedge of the two-point correlation function
       *
       *  The function computes the wedges of the two-point
       *  correlation function (see Kazin et al. 2013):
       *
       *  \f[ \xi(\delta \mu, s) = \xi_0(s)+\frac{1}{2} \left( \frac{
       *  \mu_{max}^3- \mu_{min}^3}{\mu_{max}-\mu_{min}} -1
       *  \right)\xi_2(s)+ \frac{1}{8} \left( \frac{ 7
       *  \left(\mu_{max}^5-\mu_{min}^5\right)-10 \left(\mu_{max}^3-
       *  \mu_{min}^3 \right)}{\mu_{max}-\mu_{min}}+3 \right)
       *  \xi_4(s).  \f]
       *
       *  where \f$xi_0(s), \xi_2(s), \xi_4(s)\f$ are the two-point
       *  correlation function monopoles
       *
       *  @param rr vector of scales to compute wedges 
       * 
       *  @param nwedges the number of wedges
       *
       *  @param model the \f$P(k,\mu)\f$ model; the possible options
       *  are: dispersion_dewiggled, dispersion_modecoupling
       *
       *  @param parameter vector containing parameter values
       *
       *  @param pk_interp vector containing power spectrum
       *  interpolating functions
       *
       *  @param prec the integral precision
       *
       *  @return the wedges of the two-point correlation function.
       */
      std::vector<std::vector<double>> Xi_wedges (const std::vector<double> rr, const int nwedges, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5);

      /**
       *  @brief the wedge of the two-point correlation function
       *
       *  The function computes the wedges of the two-point
       *  correlation function (see Kazin et al. 2013):
       *
       *  \f[ \xi(\delta \mu, s) = \xi_0(s)+\frac{1}{2} \left( \frac{
       *  \mu_{max}^3- \mu_{min}^3}{\mu_{max}-\mu_{min}} -1
       *  \right)\xi_2(s)+ \frac{1}{8} \left( \frac{ 7
       *  \left(\mu_{max}^5-\mu_{min}^5\right)-10 \left(\mu_{max}^3-
       *  \mu_{min}^3 \right)}{\mu_{max}-\mu_{min}}+3 \right)
       *  \xi_4(s).  \f]
       *
       *  where \f$xi_0(s), \xi_2(s), \xi_4(s)\f$ are the two-point
       *  correlation function monopoles
       *
       *  @param rr vector of scales to compute wedges 
       *
       *  @param dataset_wedge vector that specify the wedges
       *  to be computed for each scale 
       * 
       *  @param nwedges the number of wedges
       *
       *  @param model the \f$P(k,\mu)\f$ model; the possible options
       *  are: dispersion_dewiggled, dispersion_modecoupling
       *
       *  @param parameter vector containing parameter values
       *
       *  @param pk_interp vector containing power spectrum
       *  interpolating functions
       *
       *  @param prec the integral precision
       *
       *  @return the wedges of the two-point correlation function.
       */
      std::vector<double> Xi_wedges (const std::vector<double> rr, const std::vector<int> dataset_wedge, const int nwedges, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5);

      /**
       *  @brief the cartesian two-point correlation function
       *
       *  The function computes the cartesian two-point correlation
       *  function:
       *
       *  \f[ \xi_(r_p, \pi) = \xi_0(s) + \xi_2(s) \mathcal{L}_2(\mu)+
       *   \xi_4(s) \mathcal{L}_4(\mu) \f]
       *
       *  where \f$xi_0(s), \xi_2(s), \xi_4(s)\f$ are the two-point
       *  correlation function monopoles and \f$ \mathcal{L}_l(\mu)\f$
       *  are the Legendre polynomial.
       *
       *  @param rp vector of scales transverse to the line of sight 
       * 
       *  @param pi vector of scales parallel to the line of sight  
       *
       *  @param model the \f$P(k,\mu)\f$ model; the possible options
       *  are: dispersion_dewiggled, dispersion_modecoupling
       *
       *  @param parameter vector containing parameter values
       *
       *  @param pk_interp vector containing power spectrum
       *  interpolating functions    
       *
       *  @param prec the integral precision
       *
       *  @return the cartesian two-point correlation function.
       */
      std::vector<std::vector<double>> Xi_rppi (const std::vector<double> rp, const std::vector<double> pi, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5);

      /**
       *  @brief the projected two-point correlation function
       *
       *  The function computes the projected two-point correlation
       *  function from the 2D two-point correlation function in
       *  Cartesian coordinates:
       *
       *  \f[ w_p(r_p) = \int_0^{\pi_{max}} \mathrm{d}\pi \xi(r_p,
       *   \pi) \f]
       *
       *  where \f$xi(r_p, \pi)\f$ is the Cartesian two-point
       *  correlation function
       *
       *  @param rp vector of scales transverse to the line of sight  
       *
       *  @param pimax the maximum scale of integration
       *
       *  @param model the \f$P(k,\mu)\f$ model; the possible options
       *  are: dispersion_dewiggled, dispersion_modecoupling
       *
       *  @param parameter vector containing parameter values
       *
       *  @param pk_interp vector containing power spectrum
       *  interpolating functions    
       *
       *  @param prec the integral precision
       *
       *  @return the cartesian two-point correlation function.
       */
      std::vector<double> wp_from_Xi_rppi (const std::vector<double> rp, const double pimax, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5);

      /**
       *  @brief the power spectrum terms
       *  obtained integrating the redshift space 2D power spectrum
       *
       *  the function returns the analytic solutions of the integral
       *  of the redshift space 2D power spectrum along \f$\mu\f$:
       *
       *  \f$P(k,\mu) = P_\mathrm{DM}(k) (b+f\mu^2)^2
       *  \exp(-k^2\mu^2\sigma^2)\f$.
       *
       *  Solutions are :
       *
       *  \f[ 
       *   P'(k) = P_\mathrm{DM}(k) \frac{\sqrt{\pi}}{2 k \sigma} \mathrm{erf}(k\sigma) ; \\
       *   P''(k) = \frac{f}{(k\sigma)^3} P_\mathrm{DM}(k) \left[\frac{\sqrt{\pi}}{2}\mathrm{erf}(k\sigma) 
       *    -k\sigma\exp(-k^2\sigma^2)\right] ; \\
       *   P'''(k) = \frac{f^2}{(k\sigma)^5}P_\mathrm{DM}(k) \left\{ \frac{3\sqrt{\pi}}{8}\mathrm{erf}(k\sigma) \right. \\ \left. 
       *   - \frac{k\sigma}{4}\left[2(k\sigma)^2+3\right]\exp(-k^2\sigma^2)\right\} . \\
       *   \f]
       *
       *  @param kk the binned wave vector modules 
       *
       *  @param linear_growth_rate the linear growth rate
       *
       *  @param SigmaS streaming scale
       *
       *  @param PkDM dark matter power spectrum interpolator
       *
       *  @return the damped two-point correlation monopole.
       */
      std::vector<std::vector<double>> damped_Pk_terms (const std::vector<double> kk, const double linear_growth_rate, const double SigmaS, const std::shared_ptr<cbl::glob::FuncGrid> PkDM);

      /**
       *  @brief the damped two-point correlation monopole;
       *  from Sereno et al. 2015
       *
       *  The function computes the  damped two-point correlation 
       *  monopole:
       * 
       *  \f$\xi(s) = b^2 \xi'(s) + b \xi''(s) + \xi'''(s) \, ;\f$
       *
       *  where b is the linear bias and the terms \f$\xi'(s)\f$,
       *  \f$\xi''(s)\f$, \f$\xi'''(s)\f$ are
       *  the Fourier anti-transform of the power spectrum terms
       *  obtained integrating the redshift space 2D power spectrum
       *  along \f$\mu\f$ (see cbl::modelling::twopt.:damped_Pk_terms).
       *
       *  @param ss vector of scales  
       *
       *  @param bias the linear bias
       *
       *  @param linear_growth_rate the linear growth rate
       *
       *  @param SigmaS streaming scale
       *
       *  @param kk the binned wave vector modules
       *
       *  @param PkDM dark matter power spectrum interpolator
       *
       *  @return the damped two-point correlation monopole.
       */
      std::vector<double> damped_Xi (const std::vector<double> ss, const double bias, const double linear_growth_rate, const double SigmaS, const std::vector<double> kk, const std::shared_ptr<cbl::glob::FuncGrid> PkDM);

    }
  }
}

#endif
