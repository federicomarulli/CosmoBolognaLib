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

#include "Modelling_MassObservableRelation.h"
#include "Modelling_PowerSpectrum_Angular.h"

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

	/// order of the polynomial that take systematic effects into account
	int poly_order;

	/// fiducial cosmology
	std::shared_ptr<cosmology::Cosmology> cosmology;

	/// test cosmology
	std::shared_ptr<cosmology::Cosmology> test_cosmology;
	
	/// mass-observable scaling relation
	std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation;
	
	/// angular power spectrum
	std::shared_ptr<modelling::angularpk::Modelling_PowerSpectrum_angular> angularpk;
	
	/// typical absolute error on redshift
	double z_abs_err = -1;
	
	/// typical relative error on mass proxy
	double proxy_rel_err = -1;

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
	std::vector<std::string> Cpar;

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

	/// pointer to a function of FuncGrid class, used to interpolate the Pk_b2d term (eTNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_b2d;

	/// pointer to a function of FuncGrid class, used to interpolate the Pk_b2v term (eTNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_b2v;

	/// pointer to a function of FuncGrid class, used to interpolate the Pk_b22 term (eTNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_b22;

	/// pointer to a function of FuncGrid class, used to interpolate the Pk_bs2d term (eTNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_bs2d;

	/// pointer to a function of FuncGrid class, used to interpolate the Pk_bs2v term (eTNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_bs2v;

	/// pointer to a function of FuncGrid class, used to interpolate the Pk_b2s2 term (eTNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_b2s2;

	/// pointer to a function of FuncGrid class, used to interpolate the Pk_bs22 term (eTNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_Pk_bs22;

	/// pointer to a function of FuncGrid class, used to interpolate the sigma32Pklin term (eTNS) of the power spectrum
	std::shared_ptr<glob::FuncGrid> func_sigma32Pklin;

	/// pointer to a function of FuncGrid class, used to interpolate the power spectrum 1-loop correction
	std::shared_ptr<glob::FuncGrid> func_Pk1loop;

	/// pointer to a vector of FuncGrid class, used to interpolate power spectra and power spectra integrals
	std::vector<std::shared_ptr<glob::FuncGrid>> funcs_pk;

	/// pointer to a vector of FuncGrid objects, used to interpolate the multipoles
	std::vector<std::shared_ptr<glob::FuncGrid>> func_multipoles;

	/// pointer to a vector of FuncGrid objects, used to interpolate the wedges
	std::vector<std::shared_ptr<glob::FuncGrid>> func_wedges;

	/// the \f$P(k,\mu)\f$ model
	std::string Pk_mu_model;

	/// fiducial bias value
	double bias;
	
	/// if true the bias is estimated directly using mass and redshift of the halo catalogue
	bool bias_from_sample = false;
	
	/// if true the quadrupole of the 3D correlation function is included in the angular correlation function model. It works only in redshift space
	bool include_quadrupole = false;

	/// fiducial \f$D_V\f$
	double DVfid;

	/// fiducial \f$D_A(z)\f$
	double DAfid;

	/// fiducial \f$H(z)\f$
	double HHfid;
	
	/// fiducial \f$r_s\f$
	double rsfid;

	/// &Delta;, the overdensity
	double Delta;
	
	/// the input overdensity
	double Delta_input;

	/// isDelta_critical
	bool isDelta_critical;

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
	
	/// cluster redshift
	std::vector<double> cluster_redshift;

	/// cosmology used to measure the cluster masses
	std::shared_ptr<cosmology::Cosmology> cosmology_mass;

	/// redshift_source vector containing the redshifts of the source galaxies, in case the cluster masses are estimated from weak lensing
	std::vector<double> redshift_source;
	
	/// bestfit parameters for the redshift distribution of_sources
	std::vector<double> dN_par;
	
	/// the distribution type of dN/dz
	std::string distribution_type;

	/// interpoolation type 
	std::string interpType;

	/// specify if you are in real or redshift space (for function xi0_linear_bias_cosmology) 
	bool redshift_space;
	
	/// the area covered by the survey, in square degrees
	double area_degrees;
	
	/// specify angular units (for angular 2pt correlation function) 
	cbl::CoordinateUnits theta_units;

	/// the redshift evolution function in the scaling relation
	std::function<double(const double, const double, const std::shared_ptr<void>)> fz;

	/// weights for completeness and purity of number counts in lambda* redshift bins 
	std::vector<double> SF_weights;

	///minimum proxy in a given proxy bin
	double proxy_min;

	///maximum proxy in a given proxy bin
	double proxy_max;
		
	/**
	 *  @brief default constructor
	 */
	STR_data_model () = default;
      };

      /**
       *  @brief the interpolating function of multipole expansion of
       *  the two-point correlation function at a given order l
       *
       *  The function computes the multipoles of the two-point
       *  correlation function:
       *
       *  \f[ \xi_l(s) = i^l \int \frac{\mathrm{d} k}{2\pi^2} k^2
       *  P_l(k) j_l(ks) \f]
       *
       *  where \f$j_l(ks)\f$ are the Bessel functions, and
       *  \f$P_l(k)\f$ is computed by cbl::modelling::powspec::Pk_l
       *
       *  @param kk the wave vector module vector
       *
       *  @param l the order of the expansion
       *
       *  @param model the \f$P(k,\mu)\f$ model
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
       *  \f[ \xi_l(s) = i^l \int \frac{\mathrm{d} k}{2\pi^2} k^2
       *  P_l(k) j_l(ks) \f]
       *
       *  where \f$j_l(ks)\f$ are the Bessel functions, and
       *  \f$P_l(k)\f$ is computed by cbl::modelling::powspec::Pk_l
       *
       *  @param rr vector of scales to compute multipoles
       *
       *  @param nmultipoles the number of (even) multipoles to
       *  compute
       *
       *  @param model the \f$P(k,\mu)\f$ model
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
       *  \f[ \xi_l(s) = i^l \int \frac{\mathrm{d} k}{2\pi^2} k^2
       *  P_l(k) j_l(ks) \f]
       *
       *  where \f$j_l(ks)\f$ are the Bessel functions, and
       *  \f$P_l(k)\f$ is computed by cbl::modelling::powspec::Pk_l
       *
       *  @param rr vector of scales to compute multipoles
       *
       *  @param dataset_order vector that specify the multipole
       *  to be computed for each scale
       *
       *  @param use_pole vector of booleans specifying if a given
       *  multipole should be computed
       *
       *  @param model the \f$P(k,\mu)\f$ model
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
       *  @param model the \f$P(k,\mu)\f$ model
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
       *  @brief the polar two-point correlation function
       *
       *  The function computes the polar two-point correlation
       *  function from its multipoles as expressed in Kazin et
       *  al. 2013 (https://arxiv.org/pdf/1303.4391.pdf, appendix A)
       *
       *  \f[ \xi(s_{\mathrm{true}}, \mu_{\mathrm{true}}) =
       *  \sum\xi_l(s_{\mathrm{true}}(s_{\mathrm{fid}},
       *  \mu_{\mathrm{fid}}, \alpha_{\perp}, \alpha_{\parallel}))
       *  L_l(\mu_{\mathrm{true}}(\mu_{\mathrm{fid}}, \alpha_{\perp},
       *  \alpha_{\parallel})) \f]
       *
       *  where \f$\xi_l(s)\f$ are the two-point correlation function
       *  monopoles up to l=4, and \f$ \mathcal{L}_l(\mu)\f$ are the
       *  Legendre polynomial.
       *
       *  The relations between fiducial and true quantities are:
       *
       *  \f[ s_{\mathrm{true}} = s_{\mathrm{fid}}
       *     \sqrt{\alpha_{\|}^{2}
       *     \mu_{\mathrm{fid}}^{2}+\alpha_{\perp}^2\left(1-\mu_{\mathrm{fid}}^2\right)}
       *     \f]
       *
       *  \f[ \mu_{\mathrm{true}} = \mu_{\mathrm{fid}}
       *     \frac{\alpha_{\|}}{\sqrt{\alpha_{\|}^{2}
       *     \mu_{\mathrm{fid}}^{2}+\alpha_{\perp}^2\left(1-\mu_{\mathrm{fid}}^2\right)}}
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
       *  @param model the \f$P(k,\mu)\f$ model
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
       *  along \f$\mu\f$ (see cbl::modelling::powspec::damped_Pk_terms).
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
       * @param inputs pointer to the structure that contains the power spectrum angular data model
       *
       * @return values of redshift distribution as a function of redshift
       * and mass proxy
       *
       */
      double bias_from_scaling_relation(const double alpha, const double beta, const double gamma,  const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, double log_base, double mass_pivot, double proxy_pivot, double redshift_pivot, std::vector<double> scalRel_pars, cbl::glob::FuncGrid interp_DN, std::shared_ptr<void> inputs);
    }
  }
}

#endif
