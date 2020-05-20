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
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
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

	/// true \f$\rightarrow\f$ the output files created by the Boltmann solver are stored; false \f$\rightarrow\f$ the output files are removed
	bool store_output;

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
	int nWedges;

	/// integral limits used to measure the wedges
	std::vector<std::vector<double>> mu_integral_limits;

	/// vector of wave vector modules
	std::vector<double> kk;

	/// vector of scales
	std::vector<double> rr;

	/// pointer to a function of FuncGrid class, used to interpolate the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk;

	/// pointer to a function of FuncGrid class, used to interpolate the no-wiggles linear power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_NW;

	/// pointer to a function of FuncGrid class, used to interpolate the no-linear power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_nonlin;

	/// pointer to a function of FuncGrid class, used to interpolate the Pk_DeltaDelta power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_DeltaDelta;

	/// pointer to a function of FuncGrid class, used to interpolate the Pk_DeltaTheta power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_DeltaTheta;

	/// pointer to a function of FuncGrid class, used to interpolate the Pk_ThetaTheta power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_ThetaTheta;

	/// pointer to a function of FuncGrid class, used to interpolate the A11 term (TNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_A11;

	/// pointer to a function of FuncGrid class, used to interpolate the A12 term (TNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_A12;

	/// pointer to a function of FuncGrid class, used to interpolate the A22 term (TNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_A22;

	/// pointer to a function of FuncGrid class, used to interpolate the A23 term (TNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_A23;

	/// pointer to a function of FuncGrid class, used to interpolate the A33 term (TNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_A33;

	/// pointer to a function of FuncGrid class, used to interpolate the B12 term (TNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_B12;

	/// pointer to a function of FuncGrid class, used to interpolate the B13 term (TNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_B13;

	/// pointer to a function of FuncGrid class, used to interpolate the B14 term (TNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_B14;

	/// pointer to a function of FuncGrid class, used to interpolate the B22 term (TNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_B22;

	/// pointer to a function of FuncGrid class, used to interpolate the B23 term (TNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_B23;

	/// pointer to a function of FuncGrid class, used to interpolate the B24 term (TNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_B24;

	/// pointer to a function of FuncGrid class, used to interpolate the B33 term (TNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_B33;

	/// pointer to a function of FuncGrid class, used to interpolate the B34 term (TNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_B34;

	/// pointer to a function of FuncGrid class, used to interpolate the B44 term (TNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_B44;

	/// pointer to a function of FuncGrid class, used to interpolate the power spectrum 1-loop correction
	std::shared_ptr<glob::FuncGrid> func_Pk1loop;

	/// pointer to a vector of FuncGrid class, used to interpolate power spectra and power spectra integrals
	std::vector<std::shared_ptr<glob::FuncGrid>> funcs_pk;

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
       * @brief shift kk and mu with Alcock-Paczynski
       * parameters
       *
       *  @param kk the wave vector module
       *
       *  @param mu the line of sight cosine
       *
       *  @param alpha_perp the shift transverse to the l.o.s.
       *
       *  @param alpha_par the shift parallel to the l.o.s.
       *
       *  @return None
       */
      void AP_shift_FourierSpace(double &kk, double &mu, const double alpha_perp, const double alpha_par);

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
      double Pkmu_DeWiggled (const double kk, const double mu, const double sigmaNL_perp, const double sigmaNL_par, const double linear_growth_rate, const double bias, const double SigmaS, const std::shared_ptr<cbl::glob::FuncGrid> Pk, const std::shared_ptr<cbl::glob::FuncGrid> Pk_NW);

      /**
       *  @brief the power spectrum as a function of \f$k\f$ and \f$\mu\f$
       *
       *  this function computes the redshift-space power spectrum
       *  \f$P(k, \mu)\f$ for the so-called dispersion model (see e.g. Pezzotta et al. 2017
       *  https://arxiv.org/abs/1612.05645):
       *
       *  \f[ P(k, \mu) = D_{FoG}(k,\mu,\sigma_{12})\left(1+\frac{f}{b}\mu^2\right)P_{\delta\delta}(k)\f]
       *
       *  @author J.E. Garcia-Farieta
       *  @author joegarciafa@unal.edu.co
       *
       *  @param kk the wave vector module
       *
       *  @param mu the line of sight cosine
       *
       *  @param DFoG the damping factor ("Gaussian" or "Lorentzian")
       *
       *  @param linear_growth_rate the linear growth rate
       *
       *  @param bias the linear bias
       *
       *  @param sigma12 streaming scale
       *
       *  @param Pklin linear power spectrum interpolator
       *
       *  @return \f$P(k, \mu)\f$
       */
      double Pkmu_Dispersion (const double kk, const double mu, const std::string DFoG, const double linear_growth_rate, const double bias, const double sigma12, const std::shared_ptr<cbl::glob::FuncGrid> Pklin);

      /**
       *  @brief the power spectrum as a function of \f$k\f$ and \f$\mu\f$
       *
       *  this function computes the redshift-space power spectrum
       *  \f$P(k, \mu)\f$ ussing fitting functions (see e.g. Pezzotta et al. 2017
       *  https://arxiv.org/abs/1612.05645, and Mohammad et. al. 2018 https://arxiv.org/abs/1807.05999):
       *
       *  \f[ P(k, \mu) = D_{FoG}(k,\mu,\sigma_{12})\left(b^2P_{\delta\delta}(k) + 2fb\mu^2P_{\delta\theta}(k)
       *  + f^2\mu^4P_{\delta\delta}(k)\right)\f]
       *
       *  where
       *
       *  \f[ P_{\delta\theta}(k)=\left(P_{\delta\delta}P^{lin}(k)e^{-k/k*}\right)^{1/2} , \f]
       *
       *  \f[ P_{\theta\theta}(k)=P^{lin}(k)e^{-k/k*}, \f]
       *
       *  @author J.E. Garcia-Farieta
       *  @author joegarciafa@unal.edu.co
       *
       *  @param kk the wave vector module
       *
       *  @param mu the line of sight cosine
       *
       *  @param DFoG the damping factor ("Gaussian" or "Lorentzian")
       *
       *  @param linear_growth_rate the linear growth rate
       *
       *  @param bias the linear bias
       *
       *  @param sigma12 streaming scale
       *
       *  @param kd fitting parameter
       *
       *  @param kt fitting parameter
       *
       *  @param Pklin linear power spectrum interpolator
       *
       *  @param Pknonlin nolinear power spectrum interpolator
       *
       *  @return \f$P(k, \mu)\f$
       */
      double Pkmu_Scoccimarro_fitPezzotta (const double kk, const double mu, const std::string DFoG, const double linear_growth_rate, const double bias, const double sigma12, const double kd, const double kt, const std::shared_ptr<cbl::glob::FuncGrid> Pklin, const std::shared_ptr<cbl::glob::FuncGrid> Pknonlin);

      /**
       *  @brief the power spectrum as a function of \f$k\f$ and \f$\mu\f$
       *
       *  this function computes the redshift-space power spectrum
       *  \f$P(k, \mu)\f$ ussing fitting functions (see e.g. Bel et al. 2019
       *  https://arxiv.org/abs/1809.09338):
       *
       *  \f[ P(k, \mu) = D_{FoG}(k,\mu,\sigma_{12})\left(b^2P_{\delta\delta}(k) + 2fb\mu^2P_{\delta\theta}(k)
       *  + f^2\mu^4P_{\delta\delta}(k)\right)\f]
       *
       *  where
       *
       *  \f[ P_{\delta\theta}(k)=\left(P_{\delta\delta}P^{lin}(k)\right)^{1/2}e^{-k/k_\delta-bk^6} , \f]
       *
       *  \f[ P_{\theta\theta}(k)=P^{lin}(k)e^{-k(a_1+a_2k+a_3k^2)}, \f]
       *
       *  @author J.E. Garcia-Farieta
       *  @author joegarciafa@unal.edu.co
       *
       *  @param kk the wave vector module
       *
       *  @param mu the line of sight cosine
       *
       *  @param DFoG the damping factor ("Gaussian" or "Lorentzian")
       *
       *  @param linear_growth_rate the linear growth rate
       *
       *  @param bias the linear bias
       *
       *  @param sigma12 streaming scale
       *
       *  @param kd fitting parameter
       *
       *  @param bb fitting parameter
       *
       *  @param a1 fitting parameter
       *
       *  @param a2 fitting parameter
       *
       *  @param a3 fitting parameter
       *
       *  @param Pklin linear power spectrum interpolator
       *
       *  @param Pknonlin nolinear power spectrum interpolator
       *
       *  @return \f$P(k, \mu)\f$
       */
      double Pkmu_Scoccimarro_fitBel (const double kk, const double mu, const std::string DFoG, const double linear_growth_rate, const double bias, const double sigma12, const double kd, const double bb, const double a1, const double a2, const double a3, const std::shared_ptr<cbl::glob::FuncGrid> Pklin, const std::shared_ptr<cbl::glob::FuncGrid> Pknonlin);

      /**
       *  @brief the power spectrum as a function of \f$k\f$ and \f$\mu\f$
       *
       *  this function computes the redshift-space power spectrum
       *  \f$P(k, \mu)\f$ in 1-loop aproximation using (standard) Perturbation Theory (see e.g. Scoccimarro 2004
       *  https://arxiv.org/abs/astro-ph/0407214 and Taruya et. al., 2010 https://arxiv.org/abs/1006.0699):
       *
       *  \f[ P(k, \mu) = D_{FoG}(k,\mu,\sigma_{12})\left(b^2P_{\delta\delta}(k) + 2fb\mu^2P_{\delta\theta}(k)
       *  + f^2\mu^4P_{\delta\delta}(k)\right)\f]
       *
       *  where
       *
       * \f[
       * < \delta(k)\delta(k')> = (2\pi)^3\delta(k + k')P_{\delta\delta}(k)
       * \f]
       * \f[
       * < \delta(k)\theta(k')> = (2\pi)^3\delta(k + k')P_{\delta\theta}(k)
       * \f]
       * \f[
       * < \theta(k)\theta(k')> = (2\pi)^3\delta(k + k')P_{\theta\theta}(k)
       * \f]
       *
       *  @author J.E. Garcia-Farieta
       *  @author joegarciafa@unal.edu.co
       *
       *  @param kk the wave vector module
       *
       *  @param mu the line of sight cosine
       *
       *  @param DFoG the damping factor ("Gaussian" or "Lorentzian")
       *
       *  @param linear_growth_rate the linear growth rate
       *
       *  @param bias the linear bias
       *
       *  @param sigma12 streaming scale
       *
       *  @param Pk_DeltaDelta power spectrum interpolator
       *
       *  @param Pk_DeltaTheta power spectrum interpolator
       *
       *  @param Pk_ThetaTheta power spectrum interpolator
       *
       *  @return \f$P(k, \mu)\f$
       */
      double Pkmu_Scoccimarro (const double kk, const double mu, const std::string DFoG, const double linear_growth_rate, const double bias, const double sigma12, const std::shared_ptr<cbl::glob::FuncGrid> Pk_DeltaDelta, const std::shared_ptr<cbl::glob::FuncGrid> Pk_DeltaTheta, const std::shared_ptr<cbl::glob::FuncGrid> Pk_ThetaTheta);

      /**
       *  @brief the TNS power spectrum as a function of \f$k\f$ and \f$\mu\f$
       *
       *  this function computes the TNS redshift-space power spectrum
       *  \f$P(k, \mu)\f$ in 1-loop aproximation using (standard) Perturbation Theory (see e.g. Taruya et. al, 2010
       *  https://arxiv.org/abs/1006.0699 and Taruya et. al., 2013 https://arxiv.org/abs/1301.3624):
       *
       *  \f[ P(k, \mu) = D_{FoG}(k,\mu,\sigma_{12})\left(b^2P_{\delta\delta}(k) + 2fb\mu^2P_{\delta\theta}(k)
       *  + f^2\mu^4P_{\delta\delta}(k) + b^3A(k, \mu, f) + b^4B(k, \mu, f)\right)\f]
       *
       *  where
       *
       * \f[
       * < \delta(k)\delta(k')> = (2\pi)^3\delta(k + k')P_{\delta\delta}(k)
       * \f]
       * \f[
       * < \delta(k)\theta(k')> = (2\pi)^3\delta(k + k')P_{\delta\theta}(k)
       * \f]
       * \f[
       * < \theta(k)\theta(k')> = (2\pi)^3\delta(k + k')P_{\theta\theta}(k)
       * \f]
       *
       * and
       *
       * \f[
       * A(k, \mu ; f)=j_{1} \int d^{3} r e^{i \boldsymbol{k} \cdot \boldsymbol{r}}\left\langle A_{1} A_{2} A_{3}\right\rangle_{c}=k \mu f \int \frac{d^{3} p}{(2 \pi)^{3}} \frac{p_{z}}{p^{2}}\left\{B_{\sigma}(\boldsymbol{p}, \boldsymbol{k}-\boldsymbol{p},-\boldsymbol{k})-B_{\sigma}(\boldsymbol{p}, \boldsymbol{k},-\boldsymbol{k}-\boldsymbol{p})\right\}
       * \f]
       *
       * \f[
       * B(k, \mu ; f)=j_{1}^{2} \int d^{3} r e^{i \boldsymbol{k} \cdot \boldsymbol{r}}\left\langle A_{1} A_{2}\right\rangle_{c}\left\langle A_{1} A_{3}\right\rangle_{c}=(k \mu f)^{2} \int \frac{d^{3} p}{(2 \pi)^{3}} F_{\sigma}(\boldsymbol{p}) F_{\sigma}(\boldsymbol{k}-\boldsymbol{p})
       * \f]
       *
       *
       *  @author J.E. Garcia-Farieta
       *  @author joegarciafa@unal.edu.co
       *
       *  @param kk the wave vector module
       *
       *  @param mu the line of sight cosine
       *
       *  @param DFoG the damping factor ("Gaussian" or "Lorentzian")
       *
       *  @param linear_growth_rate the linear growth rate
       *
       *  @param bias the linear bias
       *
       *  @param sigma12 streaming scale
       *
       *  @param Pk_DeltaDelta power spectrum interpolator
       *
       *  @param Pk_DeltaTheta power spectrum interpolator
       *
       *  @param Pk_ThetaTheta power spectrum interpolator
       *
       *  @param Pk_A11 power spectrum interpolator
       *
       *  @param Pk_A12 power spectrum interpolator
       *
       *  @param Pk_A22 power spectrum interpolator
       *
       *  @param Pk_A23 power spectrum interpolator
       *
       *  @param Pk_A33 power spectrum interpolator
       *
       *  @param Pk_B12 power spectrum interpolator
       *
       *  @param Pk_B13 power spectrum interpolator
       *
       *  @param Pk_B14 power spectrum interpolator
       *
       *  @param Pk_B22 power spectrum interpolator
       *
       *  @param Pk_B23 power spectrum interpolator
       *
       *  @param Pk_B24 power spectrum interpolator
       *
       *  @param Pk_B33 power spectrum interpolator
       *
       *  @param Pk_B34 power spectrum interpolator
       *
       *  @param Pk_B44 power spectrum interpolator
       *
       *  @return \f$P(k, \mu)\f$
       */
      double Pkmu_TNS (const double kk, const double mu, const std::string DFoG, const double linear_growth_rate, const double bias, const double sigma12, const std::shared_ptr<cbl::glob::FuncGrid> Pk_DeltaDelta, const std::shared_ptr<cbl::glob::FuncGrid> Pk_DeltaTheta, const std::shared_ptr<cbl::glob::FuncGrid> Pk_ThetaTheta, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A11, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A12, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A22, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A23, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A33, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B12, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B13, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B14, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B22, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B23, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B24, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B33, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B34, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B44);

      /**
       *  @brief the eTNS power spectrum as a function of \f$k\f$ and \f$\mu\f$
       *
       *  this function computes the extended TNS (eTNS) redshift-space power spectrum
       *  \f$P(k, \mu)\f$ in 1-loop aproximation using (standard) Perturbation Theory (see e.g. Taruya et. al, 2010
       *  https://arxiv.org/abs/1006.0699, Taruya et. al., 2013 https://arxiv.org/abs/1301.3624) and Beutler et. al., https://arxiv.org/abs/1312.4611:
       *
       *  \f[ P(k, \mu) = D_{FoG}\left\{-\left(f k \mu \sigma_{v}\right)^{2}\right\}\left[P_{\mathrm{g}, \delta \delta}(k) +2 f \mu^{2} P_{\mathrm{g}, \delta \theta}(k)+f^{2} \mu^{4} P_{\theta \theta}(k)+b_{1}^{3} A(k, \mu, \beta)+b_{1}^{4} B(k, \mu, \beta)\right]\f]
       *
       *  where
       *
       * \f[
       * P_{\mathrm{g}, \delta \delta}(k)= b_{1}^{2} P_{\delta \delta}(k)+2 b_{2} b_{1} P_{b 2, \delta}(k)+2 b_{s 2} b_{1} P_{b s 2, \delta}(k) \\ +2 b_{3 \mathrm{nl}} b_{1} \sigma_{3}^{2}(k) P_{\mathrm{m}}^{\mathrm{L}}(k)+b_{2}^{2} P_{b 22}(k) \\ +2 b_{2} b_{s 2} P_{b 2 s 2}(k)+b_{s 2}^{2} P_{b s 22}(k)+N
       * \f]
       * \f[
       * P_{\mathrm{g}, \delta \theta}(k)= b_{1} P_{\delta \theta}(k)+b_{2} P_{b 2, \theta}(k)+b_{s 2} P_{b s 2, \theta}(k) \\ +b_{3 \mathrm{n}}\left[\rho_{3}^{2}(k) P_{\mathrm{m}}^{\mathrm{lin}}(k)\right.
       * \f]
       *
       * and
       *
       * \f[
       * A(k, \mu ; f)=j_{1} \int d^{3} r e^{i \boldsymbol{k} \cdot \boldsymbol{r}}\left\langle A_{1} A_{2} A_{3}\right\rangle_{c}=k \mu f \int \frac{d^{3} p}{(2 \pi)^{3}} \frac{p_{z}}{p^{2}}\left\{B_{\sigma}(\boldsymbol{p}, \boldsymbol{k}-\boldsymbol{p},-\boldsymbol{k})-B_{\sigma}(\boldsymbol{p}, \boldsymbol{k},-\boldsymbol{k}-\boldsymbol{p})\right\}
       * \f]
       *
       * \f[
       * B(k, \mu ; f)=j_{1}^{2} \int d^{3} r e^{i \boldsymbol{k} \cdot \boldsymbol{r}}\left\langle A_{1} A_{2}\right\rangle_{c}\left\langle A_{1} A_{3}\right\rangle_{c}=(k \mu f)^{2} \int \frac{d^{3} p}{(2 \pi)^{3}} F_{\sigma}(\boldsymbol{p}) F_{\sigma}(\boldsymbol{k}-\boldsymbol{p})
       * \f]
       *
       *
       *  @author J.E. Garcia-Farieta
       *  @author joegarciafa@unal.edu.co
       *
       *  @param kk the wave vector module
       *
       *  @param mu the line of sight cosine
       *
       *  @param DFoG the damping factor ("Gaussian" or "Lorentzian")
       *
       *  @param linear_growth_rate the linear growth rate
       *
       *  @param bias the linear bias
       *
       *  @param bias2 the second order local bias
       *
       *  @param sigma12 streaming scale
       *
       *  @param Ncorr constant stochasticity term
       *
       *  @param Pk_DeltaDelta power spectrum interpolator
       *
       *  @param Pk_DeltaTheta power spectrum interpolator
       *
       *  @param Pk_ThetaTheta power spectrum interpolator
       *
       *  @param Pk_A11 power spectrum interpolator
       *
       *  @param Pk_A12 power spectrum interpolator
       *
       *  @param Pk_A22 power spectrum interpolator
       *
       *  @param Pk_A23 power spectrum interpolator
       *
       *  @param Pk_A33 power spectrum interpolator
       *
       *  @param Pk_B12 power spectrum interpolator
       *
       *  @param Pk_B13 power spectrum interpolator
       *
       *  @param Pk_B14 power spectrum interpolator
       *
       *  @param Pk_B22 power spectrum interpolator
       *
       *  @param Pk_B23 power spectrum interpolator
       *
       *  @param Pk_B24 power spectrum interpolator
       *
       *  @param Pk_B33 power spectrum interpolator
       *
       *  @param Pk_B34 power spectrum interpolator
       *
       *  @param Pk_B44 power spectrum interpolator
       *
       *  @param Pk_b2d power spectrum interpolator
       *
       *  @param Pk_b2v power spectrum interpolator
       *
       *  @param Pk_b22 power spectrum interpolator
       *
       *  @param Pk_bs2d power spectrum interpolator
       *
       *  @param Pk_bs2v power spectrum interpolator
       *
       *  @param Pk_b2s2 power spectrum interpolator
       *
       *  @param Pk_bs22 power spectrum interpolator
       *
       *  @param sigma32Pklin power spectrum interpolator
       *
       *  @return \f$P(k, \mu)\f$
       */
      double Pkmu_eTNS (const double kk, const double mu, const std::string DFoG, const double linear_growth_rate, const double bias, const double bias2, const double sigma12, const double Ncorr, const std::shared_ptr<cbl::glob::FuncGrid> Pk_DeltaDelta, const std::shared_ptr<cbl::glob::FuncGrid> Pk_DeltaTheta, const std::shared_ptr<cbl::glob::FuncGrid> Pk_ThetaTheta, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A11, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A12, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A22, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A23, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A33, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B12, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B13, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B14, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B22, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B23, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B24, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B33, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B34, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B44, const std::shared_ptr<cbl::glob::FuncGrid> Pk_b2d, const std::shared_ptr<cbl::glob::FuncGrid> Pk_b2v, const std::shared_ptr<cbl::glob::FuncGrid> Pk_b22, const std::shared_ptr<cbl::glob::FuncGrid> Pk_bs2d, const std::shared_ptr<cbl::glob::FuncGrid> Pk_bs2v, const std::shared_ptr<cbl::glob::FuncGrid> Pk_b2s2, const std::shared_ptr<cbl::glob::FuncGrid> Pk_bs22, const std::shared_ptr<cbl::glob::FuncGrid> sigma32Pklin);

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
      double Pkmu_ModeCoupling (const double kk, const double mu, const double linear_growth_rate, const double bias, const double sigmav, const double AMC, const std::shared_ptr<cbl::glob::FuncGrid> PkLin, const std::shared_ptr<cbl::glob::FuncGrid> PkMC);

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
       *    cbl::modelling::twopt::Pkmu_ModeCoupling
       *
       *  Alcock-Paczynski effect has been introduced following Ballinger
       *  et al. 1998 (https://arxiv.org/pdf/astro-ph/9605017.pdf),
       *  Beutler et al. 2016, sec 5.2 (https://arxiv.org/pdf/1607.03150.pdf)
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
       *  @param alpha_perp the shift transverse to the l.o.s.
       *
       *  @param alpha_par the shift parallel to the l.o.s.
       *
       *  @return \f$P(k, \mu)\f$
       */
      double Pkmu (const double kk, const double mu, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double alpha_perp = 1., const double alpha_par = 1.);

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
       *  @param alpha_perp the shift transverse to the l.o.s.
       *
       *  @param alpha_par the shift parallel to the l.o.s.
       *
       *  @return the multipole expansion of \f$P(k, \mu)\f$ at given
       *  \f$k\f$
       */
      double Pk_l (const double kk, const int l, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5, const double alpha_perp = 1., const double alpha_par = 1.);

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
       *  @param alpha_perp the shift transverse to the l.o.s.
       *
       *  @param alpha_par the shift parallel to the l.o.s.
       *
       *  @return the multipole expansion of \f$P(k, \mu)\f$ at given
       *  \f$k\f$
       */
      std::vector<double> Pk_l (const std::vector<double> kk, const int l, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5, const double alpha_perp = 1., const double alpha_par = 1.);

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
       *  @param alpha_perp the shift transverse to the l.o.s.
       *
       *  @param alpha_par the shift parallel to the l.o.s.
       *
       *  @return the interpolation function of the multipole
       *  expansion of two-point correlation function at a given order
       *  l
       */
       cbl::glob::FuncGrid Xil_interp (const std::vector<double> kk, const int l, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5, const double alpha_perp = 1., const double alpha_par = 1.);


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
       *  @param alpha_perp the shift transverse to the l.o.s.
       *
       *  @param alpha_par the shift parallel to the l.o.s.
       *
       *  @return the multipole of order l of the two-point
       *  correlation function
       */
      std::vector<std::vector<double>> Xi_l (const std::vector<double> rr, const int nmultipoles, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5, const double alpha_perp = 1., const double alpha_par = 1.);

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
       *  @param alpha_perp the shift transverse to the l.o.s.
       *
       *  @param alpha_par the shift parallel to the l.o.s.
       *
       *  @return the multipole of order l of the two-point
       *  correlation function
       */
      std::vector<double> Xi_l (const std::vector<double> rr, const std::vector<int> dataset_order, const std::vector<bool> use_pole, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5, const double alpha_perp = 1., const double alpha_par = 1.);

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
       *  @param alpha_perp the shift transverse to the l.o.s.
       *
       *  @param alpha_par the shift parallel to the l.o.s.
       *
       *  @return the cartesian two-point correlation function.
       */
      std::vector<std::vector<double>> Xi_rppi (const std::vector<double> rp, const std::vector<double> pi, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5, const double alpha_perp = 1., const double alpha_par = 1.);

      /**
       *  @brief the cartesian two-point correlation function
       *
       *  The function computes the polar two-point correlation
       *  function from its multipoles as expressed in Kazin et
       *  al. 2013 (https://arxiv.org/pdf/1303.4391.pdf, appendix A)
       *
       *  \f[ \xi_(s_{true}, \mu_{true}) = \sum \xi_l(s_{true}
       *  		     			(s_{fid}, \mu_{fid},
       *  		     			\alpha_{\perp},
       *  		     			\alpha_{\par}) )
       *  		     			L_l(\mu_{true}(
       *  		     			\mu_{fid},
       *  		     			\alpha_{\perp},
       *  		     			\alpha_{\par}))) \f]
       *
       *  where \f$xi_l(s)\f$ are the two-point correlation function
       *  monopoles up to l=4, and \f$ \mathcal{L}_l(\mu)\f$ are the
       *  Legendre polynomial.
       *
       *  The relations between fiducial and true quantities are:
       *
       *  \f[ &s_{\mathrm{true}}=s_{\mathrm{fid}} \cdot
       *     \sqrt{\alpha_{\|}^{2}
       *     \mu_{\mathrm{fid}}}^{2}+\alpha_{\perp}^{2}\left(1-\mu_{\mathrm{fid}}^{2}\right)}\\
       *     &\mu_{\mathrm{true}}=\mu_{\mathrm{fid}}
       *     \frac{\alpha_{\|}}{\sqrt{\alpha_{\|}^{2}
       *     \mu_{\mathrm{f}}^{2}+\alpha_{\perp}^{2}\left(1-\mu_{\mathrm{fid}}^{2}\right)}}
       *     \f]
       *
       *  @param rad_fid fiducial separation
       *
       *  @param mu_fid fiducial \f$\mu\f$
       *
       *  @param alpha_perpendicular Alcock-Paczynski perpendicular
       *  parameter
       *
       *  @param alpha_parallel Alcock-Paczynski perpendicular
       *  parameter
       *
       *  @param xi_multipoles vector containing two-point correlation
       *  function multipoles interpolating functions
       *
       *  @return the polar two-point correlation function.
       */
      double Xi_polar (const double rad_fid, const double mu_fid, const double alpha_perpendicular, const double alpha_parallel, const std::vector<std::shared_ptr<cbl::glob::FuncGrid>> xi_multipoles);

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
       *  @param alpha_perp the shift transverse to the l.o.s.
       *
       *  @param alpha_par the shift parallel to the l.o.s.
       *
       *  @return the cartesian two-point correlation function.
       */
      std::vector<double> wp_from_Xi_rppi (const std::vector<double> rp, const double pimax, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5, const double alpha_perp = 1., const double alpha_par = 1.);

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
