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
 *  @file Headers/Lib/ModelFunction_TwoPointCorrelation1D_monopole.h
 *
 *  @brief Global functions to model the monopole of the two-point
 *  correlation function
 *
 *  This file contains all the prototypes of the functions used
 *  to model the monopole of the two-point correlation function
 *  
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODFUNCTWOP1DMON__
#define __MODFUNCTWOP1DMON__

#include "ModelFunction_TwoPointCorrelation1D.h"


// ============================================================================


namespace cosmobl {

  namespace modelling {

    namespace twopt {

      /**
       *  @brief model for the BAO signal in the monopole of the 
       *  two-point correlation function 
       *
       *  the function computes:
       *
       *  \f[\xi_0(s) = B^2 \cdot \xi_{\rm DM}(\alpha\cdot s, \Sigma_{NL}) 
       *  + A_0 + A_1/s + A_2/s^2\f]
       *
       *  the model has 6 parameters: 
       *    - \f$\Sigma_{NL}\f$
       *    - \f$\alpha\f$
       *    - \f$B\f$
       *    - \f$A_0\f$
       *    - \f$A_1\f$
       *    - \f$A_2\f$ 
       *
       *  the dark matter two-point correlation function is fixed and
       *  provided in input
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  dark matter power spectrum
       *
       *  @param parameter 6D vector containing the input parameters
       *
       *  @return the monopole of the two-point correlation function
       */
      vector<double> xi0_BAO_sigmaNL (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief model for the monopole of the two-point correlation
       *  function
       *
       *  the function computes:
       *
       *  \f[\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3} f\sigma_8
       *  \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2 \right] \cdot
       *  \xi_{\rm DM}(\alpha\cdot s)/\sigma_8^2 + A_0 + A_1/s +
       *  A_2/s^2\f]
       *
       *  the model has 6 parameters: 
       *    - \f$\alpha\f$
       *    - \f$f(z)\sigma_8(z)\f$
       *    - \f$b(z)\sigma_8(z)\f$
       *    - \f$A_0\f$
       *    - \f$A_1\f$
       *    - \f$A_2\f$ 
       *
       *  the dark matter two-point correlation function is fixed and
       *  provided in input
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  dark matter two-point correlation function and
       *  \f$\sigma_8(z)\f$, computed at a given (fixed) cosmology
       *
       *  @param parameter 6D vector containing the input parameters
       *
       *  @return the monopole of the two-point correlation function
       */
      vector<double> xi0_linear (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);
      
      /**
       *  @brief model for the monopole of the two-point correlation
       *  function
       *
       *  the function computes:
       *
       *  \f$\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3} f\sigma_8
       *  \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2 \right] \cdot
       *  \xi_{\rm DM}(\alpha\cdot s)/\sigma_8^2 + A_0 + A_1/s +
       *  A_2/s^2\f$
       *
       *  the model has 6 parameters: 
       *    - \f$\alpha\f$
       *    - \f$f(z)\sigma_8(z)\f$
       *    - \f$b(z)\sigma_8(z)\f$
       *    - \f$A_0\f$
       *    - \f$A_1\f$
       *    - \f$A_2\f$ 
       *
       *  the dark matter two-point correlation function is fixed and
       *  provided in input
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  dark matter two-point correlation function and
       *  \f$\sigma_8(z)\f$, computed at a given (fixed) cosmology
       *
       *  @param parameter 6D vector containing the input parameters
       *
       *  @return the monopole of the two-point correlation function
       */
      vector<double> xi0_linear_LinearPoint (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief model for the monopole of the two-point correlation
       *  function
       *
       *  the function computes the monopole as a polynomial of
       *  used-defined order
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  dark matter two-point correlation function and
       *  \f$\sigma_8(z)\f$, computed at a given (fixed) cosmology
       *
       *  @param parameter 6D vector containing the input parameters
       *
       *  @return the monopole of the two-point correlation function
       */
      vector<double> xi0_polynomial_LinearPoint (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief model for the monopole of the two-point correlation
       *  function in redshift space
       *
       *  the function computes:
       *
       *  \f[\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3}
       *  f(z)\sigma_8(z) \cdot b(z)\sigma_8(z) +
       *  \frac{1}{5}(f(z)^2\sigma_8(z)^2) \right] \cdot \xi_{\rm
       *  DM}(s)\frac{\sigma_8^2}{\sigma_8^2(z)} \f]
       *
       *  the model has 2 parameters: 
       *    - \f$\sigma_8(z)\f$
       *    - \f$b(z)\f$
       *
       *  the dark matter two-point correlation function is computed
       *  using the input cosmological parameters
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector containing the linear bias
       *
       *  @return the monopole of the two-point correlation function
       */
      vector<double> xi0_linear_sigma8_bias (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);

      
      /**
       *  @brief model for the monopole of the two-point correlation
       *  function
       *
       *  the function computes:
       *
       *  \f[\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3} f\sigma_8
       *  \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2 \right] \cdot
       *  \xi_{\rm DM}(\alpha\cdot s)/\sigma_8^2 + A_0 + A_1/s +
       *  A_2/s^2\f]
       *
       *  the model has 6 parameters: 
       *    - \f$\alpha\f$
       *    - \f$f(z)\sigma_8(z)\f$
       *    - \f$b(z)\sigma_8(z)\f$
       *    - \f$A_0\f$
       *    - \f$A_1\f$
       *    - \f$A_2\f$ 
       *
       *  the dark matter two-point correlation function is computed
       *  using the input cosmological parameters
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector containing the linear bias
       *
       *  @return the monopole of the two-point correlation function
       */
      vector<double> xi0_linear_cosmology (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief the damped two-point correlation monopole;
       *  from Sereno et al. 2015
       *
       *  The function computes the damped two-point correlation 
       *  monopole:
       * 
       *  \f$\xi(s) = b^2 \xi'(s) + b \xi''(s) + \xi'''(s) \, ;\f$
       *
       *  where b is the linear bias and the terms \f$\xi'(s)\f$,
       *  \f$\xi''(s)\f$, \f$\xi'''(s)\f$ are
       *  the Fourier anti-transform of the power spectrum terms
       *  obtained integrating the redshift space 2D power spectrum
       *  along \f$\mu\f$ (see cosmobl::modelling::twopt.:damped_Pk_terms,
       *  see cosmobl::modelling::twopt.:damped_Xi).
       *       
       *  the model has 2 parameters:
       *      - bias the linear bias
       *      - \f$\sigma_z$\f$ the redshift error
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector containing the linear bias
       *
       *  @return the damped two-point correlation monopole.
       */
      vector<double> xi0_damped_bias_sigmaz (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief the damped two-point correlation monopole;
       *  from Sereno et al. 2015
       *
       *  The function computes the damped two-point correlation 
       *  monopole:
       * 
       *  \f$\xi(s) = b^2 \xi'(s) + b \xi''(s) + \xi'''(s) \, ;\f$
       *
       *  where b is the linear bias and the terms \f$\xi'(s)\f$,
       *  \f$\xi''(s)\f$, \f$\xi'''(s)\f$ are
       *  the Fourier anti-transform of the power spectrum terms
       *  obtained integrating the redshift space 2D power spectrum
       *  along \f$\mu\f$ (see cosmobl::modelling::twopt.:damped_Pk_terms,
       *  see cosmobl::modelling::twopt.:damped_Xi).
       *       
       *  the model has 4 parameters:
       *      - M0 the intercept of the scaling relation
       *      - slope the slope of the scaling relation
       *      - scatter the scattr of the scaling relation
       *      - \f$\sigma_z$\f$ the redshift error
       *
       *  the linear bias is computed and used in the modelling.
       *  It is provided as an output parameter
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector containing the linear bias
       *
       *  @return the damped two-point correlation monopole.
       */
      vector<double> xi0_damped_scaling_relation_sigmaz (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief model for the monopole of the two-point correlation
       *  function, the bias is computed by the input cluster masses,
       *  with only \f$sigma_8\f$ as a free parameter
       *
       *  the function computes:
       *
       *  \f[\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3} f\sigma_8
       *  \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2 \right] \cdot
       *  \xi_{\rm DM}(s)/\sigma_8^2\f]
       *
       *  the model has 1 parameter: \f$\sigma_8\f$
       *
       *  the dark matter two-point correlation function and the
       *  linear effective bias are computed using the input
       *  cosmological parameters
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector containing the linear bias
       *
       *  @return the monopole of the two-point correlation function
       */
      vector<double> xi0_linear_sigma8_clusters (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief model for the monopole of the two-point correlation
       *  function, the bias is computed by the input cluster masses,
       *  with only \f$sigma_8\f$ as a free parameter
       *
       *  the function computes:
       *
       *  \f[\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3} f\sigma_8
       *  \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2 \right] \cdot
       *  \xi_{\rm DM}(s)/\sigma_8^2\f]
       *
       *  the model has 1 cosmological parameter
       *
       *  the dark matter two-point correlation function and the
       *  linear effective bias are computed using the input
       *  cosmological parameters
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector containing the cosmological
       *  parameter
       *
       *  @return the monopole of the two-point correlation function
       */
      vector<double> xi0_linear_one_cosmo_par_clusters (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief model for the monopole of the two-point correlation
       *  function, the bias is computed by the input cluster masses,
       *  with only \f$sigma_8\f$ as a free parameter
       *
       *  the function computes:
       *
       *  \f[\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3} f\sigma_8
       *  \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2 \right] \cdot
       *  \xi_{\rm DM}(s)/\sigma_8^2\f]
       *
       *  the model has 2 cosmological parameters
       *
       *  the dark matter two-point correlation function and the
       *  linear effective bias are computed using the input
       *  cosmological parameters
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 2D vector containing the cosmological
       *  parameters
       *
       *  @return the monopole of the two-point correlation function
       */
      vector<double> xi0_linear_two_cosmo_pars_clusters (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);
      
      /**
       *  @brief model for the monopole of the two-point correlation
       *  function, the bias is computed by the input cluster masses
       *
       *  the function computes:
       *
       *  \f[\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3} f\sigma_8
       *  \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2 \right] \cdot
       *  \xi_{\rm DM}(s)/\sigma_8^2\f]
       *
       *  the model has N cosmological parameters
       *
       *  the dark matter two-point correlation function and the
       *  linear effective bias are computed using the input
       *  cosmological parameters
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector containing the linear bias
       *
       *  @return the monopole of the two-point correlation function
       */
      vector<double> xi0_linear_cosmology_clusters (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);
      
      /**
       *  @brief model for the monopole of the two-point correlation
       *  function
       *
       *  the function computes:
       *
       *  \f[\xi_0(s) = \left[ (b)^2 + \frac{2}{3}f \cdot b + \frac{1}{5}(f)^2 \right] 
       *  \cdot \xi_{\rm DM}(\frac{D_V(redshift)}{D_V^{fid}(redshift}\cdot s) \f]
       *
       *  the model has 1+n parameters: 
       *    - \f$b\f$
       *    - cosmological paramters
       *
       *  the dark matter two-point correlation function is computed
       *  using the input cosmological parameters
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector containing the linear bias
       *
       *  @return the monopole of the two-point correlation function
       */
      vector<double> xi0_linear_bias_cosmology (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @name Functions of the Halo Occupation Distribution (HOD) models
       */
      ///@{

      /**
       *  @brief the average number of central galaxies hosted in a
       *  dark matter halo of a given mass
       *
       *  this function computes the average number of central
       *  galaxies hosted in a dark matter halo of mass \f$M_h\f$ as
       *  follows (Harikane et al. 2017; see also e.g. Zheng et
       *  al. 2005 for a similar formalism):
       *
       *  \f[N_{cen}(M_h)\equiv<N_{cen}|M_h>=\frac{1}{2}\left[1+{\rm
       *  erf}\left(\frac{\log M_h-\log M_{min}}{\sqrt{2}\sigma_{\log
       *  M_h}}\right)\right]\f]
       *     
       *  @param Mass the mass of the hosting dark matter halo
       *
       *  @param Mmin \f$M_{min}\f$: the mass scale at which 50% of
       *  haloes host a central galaxy
       *
       *  @param sigmalgM \f$\sigma_{\log M_h}\f$: transition width
       *  reflecting the scatter in the luminosity-halo mass relation
       *
       *  @return \f$N_{cen}(M_h)\f$
       */
      double Ncen (const double Mass, const double Mmin, const double sigmalgM);

      /**
       *  @brief the average number of satellite galaxies hosted in a
       *  dark matter halo of a given mass
       *
       *  this function computes the average number of satellite
       *  galaxies hosted in a dark matter halo of mass \f$M_h\f$ as
       *  follows (Harikane et al. 2017; see also e.g. Zheng et
       *  al. 2005 for a similar formalism):
       *
       *  \f[N_{sat}(M_h)\equiv<N_{sat}|M_h> =
       *  N_{cen}(M_h)\left(\frac{M_h-M_0}{M_1}\right)^\alpha\f]
       *     
       *  where \f$N_{cen}\f$ is computed by
       *  cosmobl::modelling::twopt::Ncen
       *
       *  @param Mass the mass of the hosting dark matter halo
       *
       *  @param Mmin \f$M_{min}\f$: the mass scale at which 50% of
       *  haloes host a satellite galaxy
       *
       *  @param sigmalgM \f$\sigma_{\log M_h}\f$: transition width
       *  reflecting the scatter in the luminosity-halo mass relation
       *
       *  @param M0 \f$M_0\f$: the cutoff mass 
       *
       *  @param M1 \f$M_1\f$: the amplitude of the power law
       *
       *  @param alpha \f$\alpha\f$: the slope of the power law
       *
       *  @return \f$N_{sat}(M_h)\f$
       */
      double Nsat (const double Mass, const double Mmin, const double sigmalgM, const double M0, const double M1, const double alpha);
      
      /**
       *  @brief the average number of galaxies hosted in a dark
       *  matter halo of a given mass
       *
       *  the function computes the average number of galaxies hosted
       *  in a dark matter halo of mass \f$M_h\f$ as follows:
       *
       *  \f[N_{gal}(M_h)\equiv<N_{gal}|M_h> = <N_{cen}|M_h> +
       *  <N_{sat}|M_h>\f]
       *
       *  where \f$<N_{cen}|M_h>\f$ and \f$<N_{sat}|M_h>\f$ are
       *  computed by cosmobl::modelling::twopt::Ncen and
       *  cosmobl::modelling::twopt::Nsat, respectively
       *
       *  @param Mass the mass of the hosting dark matter halo
       *
       *  @param Mmin \f$M_{min}\f$: the mass scale at which 50% of
       *  haloes host a satellite galaxy
       *
       *  @param sigmalgM \f$\sigma_{\log M_h}\f$: transition width
       *  reflecting the scatter in the luminosity-halo mass relation
       *
       *  @param M0 \f$M_0\f$: the cutoff mass 
       *
       *  @param M1 \f$M_1\f$: the amplitude of the power law
       *
       *  @param alpha \f$\alpha\f$: the slope of the power law
       *
       *  @return \f$N_{gal}(M_h)\f$
       */
      double Navg (const double Mass, const double Mmin, const double sigmalgM, const double M0, const double M1, const double alpha);

      /**
       *  @brief the integrand function to compute the galaxy number
       *  density
       *
       *  this function is the integrand to compute the galaxy number
       *  density; it returns the following quantity:
       *
       *  \f[n_h(M_h, z)N_{gal}(M_h)\f]
       *
       *  where the halo mass function, \f$n_h(M_h, z)=dn/dM_h\f$, is
       *  computed by cosmology::Cosmology::mass_function, and the
       *  average number of galaxies hosted in a dark matter halo of a
       *  given mass, \f$N_{gal}(M_h)\f$, is computed by
       *  cosmobl::modelling::twopt::Navg
       *
       *  @param mass \f$M_{h}\f$: the halo mass
       *
       *  @param Mmin \f$M_{min}\f$: the mass scale at which 50% of
       *  haloes host a satellite galaxy
       *
       *  @param sigmalgM \f$\sigma_{\log M_h}\f$: transition width
       *  reflecting the scatter in the luminosity-halo mass relation
       *
       *  @param M0 \f$M_0\f$: the cutoff mass 
       *
       *  @param M1 \f$M_1\f$: the amplitude of the power law
       *
       *  @param alpha \f$\alpha\f$: the slope of the power law
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @return the integrand function to compute galaxy number
       *  density
       *
       *  @warning in the current implementation, the integral is
       *  actually computed in the range \f$10^{10}-10^{16}\f$, using 
       */
      double ng_integrand (const double mass, const double Mmin, const double sigmalgM, const double M0, const double M1, const double alpha, const shared_ptr<void> inputs);

      /**
       *  @brief the galaxy number density
       *
       *  this function computes the galaxy number density as follows:
       *
       *  \f[n_{gal}(z) = \int_0^{M_{max}}n_h(M_h, z)N_{gal}(M_h)\,{\rm
       *  d}M_h\f]
       *
       *  where the halo mass function, \f$n_h(M_h, z)=dn/dM_h\f$, is
       *  computed by cosmology::Cosmology::mass_function, and the
       *  average number of galaxies hosted in a dark matter halo of a
       *  given mass, \f$N_{gal}(M_h)\f$, is computed by
       *  cosmobl::modelling::twopt::Navg; the integrand function is
       *  computed by cosmobl::modelling::twopt::ng_integrand
       *
       *  @param Mmin \f$M_{min}\f$: the mass scale at which 50% of
       *  haloes host a satellite galaxy
       *
       *  @param sigmalgM \f$\sigma_{\log M_h}\f$: transition width
       *  reflecting the scatter in the luminosity-halo mass relation
       *
       *  @param M0 \f$M_0\f$: the cutoff mass 
       *
       *  @param M1 \f$M_1\f$: the amplitude of the power law
       *
       *  @param alpha \f$\alpha\f$: the slope of the power law
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @return the galaxy number density
       *
       *  @warning in the current implementation, the integral is
       *  actually computed in the range \f$10^{10}-10^{16}\f$, using 
       */
      double ng (const double Mmin, const double sigmalgM, const double M0, const double M1, const double alpha, const shared_ptr<void> inputs);
      
      /**
       *  @brief the mean galaxy bias
       *
       *  this function computes the mean galaxy bias (e.g. Berlind &
       *  Weinberg 2002, van den Bosch et al. 2012):
       *
       *  \f[\bar{b}(z) =
       *  \frac{1}{\bar{n}_{gal}(z)}\int_{M_{min}}^{M_{max}}
       *  <N_{gal}|M>\, b_{halo}(M, z)\, n_{halo}(M, z) \,{\rm d}M\f]
       *
       *  where \f$\bar{n}_{gal}\f$ is the mean number density of
       *  galaxies, \f$<N_{gal}|M>\f$ is the mean number of galaxies
       *  hosted in haloes of mass M, \f$b_{halo}\f$ is the linear
       *  halo bias and the halo mass function, \f$n_h(M_h,
       *  z)=dn/dM_h\f$, is computed by
       *  cosmology::Cosmology::mass_function
       *
       *  @param Mmin \f$M_{min}\f$: the mass scale at which 50% of
       *  haloes host a satellite galaxy
       *
       *  @param sigmalgM \f$\sigma_{\log M_h}\f$: transition width
       *  reflecting the scatter in the luminosity-halo mass relation
       *
       *  @param M0 \f$M_0\f$: the cutoff mass 
       *
       *  @param M1 \f$M_1\f$: the amplitude of the power law
       *
       *  @param alpha \f$\alpha\f$: the slope of the power law
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @return the mean galaxy bias
       */
      double bias (const double Mmin, const double sigmalgM, const double M0, const double M1, const double alpha, const shared_ptr<void> inputs);

      /**
       *  @brief the mean number of central-satellite galaxy pairs
       *
       *  this function computes the mean number of central-satellite
       *  galaxy pairs, \f$<N_{cen}N_{sat}>(M_h)\f$, as follows:
       *
       *  \f[<N_{cen}N_{sat}>(M_h)=N_{cen}(M_h)N_{sat}(M_h)\f]
       * 
       *  where the average number of central and satellite galaxies,
       *  \f$N_{cen}\f$ and \f$N_{sat}\f$ are computed by
       *  cosmobl::modelling::twopt::Ncen and
       *  cosmobl::modelling::twopt::Nsat, respectively
       *
       *  @param Mass the mass of the hosting dark matter halo
       *
       *  @param Mmin \f$M_{min}\f$: the mass scale at which 50% of
       *  haloes host a satellite galaxy
       *
       *  @param sigmalgM \f$\sigma_{\log M_h}\f$: transition width
       *  reflecting the scatter in the luminosity-halo mass relation
       *
       *  @param M0 \f$M_0\f$: the cutoff mass 
       *
       *  @param M1 \f$M_1\f$: the amplitude of the power law
       *
       *  @param alpha \f$\alpha\f$: the slope of the power law
       *
       *  @return \f$<N_{cen}N_{sat}>(M_h)\f$
       *
       *  @warning The current implementation is valid only for
       *  Poisson distribution of the satellite galaxy's distribution
       *  (see e.g. Harikane et al 2016)
       */
      double NcNs (const double Mass, const double Mmin, const double sigmalgM, const double M0, const double M1, const double alpha);
      
      /**
       *  @brief the mean number of satellite-satellite galaxy pairs
       *
       *  this function computes the mean number of
       *  satellite-satellite galaxy pairs,
       *  \f$<N_{sat}(N_{sat}-1)>(M_h)\f$, as follows:
       *
       *  \f[<N_{sat}(N_{sat}-1)>(M_h)=N_{sat}^2(M_h)\f]
       *
       *  where the average number of satellite galaxies,
       *  \f$N_{sat}\f$ are computed by
       *  cosmobl::modelling::twopt::Nsat
       *
       *  @param Mass the mass of the hosting dark matter halo
       *
       *  @param Mmin \f$M_{min}\f$: the mass scale at which 50% of
       *  haloes host a satellite galaxy
       *
       *  @param sigmalgM \f$\sigma_{\log M_h}\f$: transition width
       *  reflecting the scatter in the luminosity-halo mass relation
       *
       *  @param M0 \f$M_0\f$: the cutoff mass 
       *
       *  @param M1 \f$M_1\f$: the amplitude of the power law
       *
       *  @param alpha \f$\alpha\f$: the slope of the power law
       *
       *  @return \f$<N_{sat}(N_{sat}-1)>(M_h)\f$
       *
       *  @warning The current implementation is valid only for
       *  Poisson distribution of the satellite galaxy's distribution
       *  (see e.g. Harikane et al 2016)
       */
      double NsNs1 (const double Mass, const double Mmin, const double sigmalgM, const double M0, const double M1, const double alpha);
      
      /**
       *  @brief the integrand function to compute the numerator of
       *  the central-satellite part of the 1-halo term of the power
       *  spectrum
       *
       *  this function is the integrand to compute the
       *  central-satellite part of the 1-halo term of power spectrum
       *  (e.g. Harikane et al. 2016); it returns the following
       *  quantity:
       *
       *  \f[<N_{cen}N_{sat}>(M_h)\,n_h(M_h, z)\,\tilde{u}_h(k, M_h,
       *  z)\f]
       *
       *  where the galaxy number density, \f$n_{gal}(z)\f$ is
       *  computed by cosmobl::modelling::twopt::ng,
       *  \f$<N_{cen}N_{sat}>(M_h)\f$ is computed by
       *  cosmobl::modelling::twopt::NcNs, the Fourier transform of
       *  the halo density profile, \f$\tilde{u}_h(k, M_h, z)\f$ if
       *  computed by
       *  cosmobl::cosmology::Cosmology::density_profile_FourierSpace
       *  and the halo mass function, \f$n_h(M_h, z)=dn/dM_h\f$, is
       *  computed by cosmology::Cosmology::mass_function
       *
       *  @param mass \f$M_{h}\f$: the halo mass
       *
       *  @param kk the wave vector module at which the model is
       *  computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the integrand function to compute central-satellite
       *  part of the 1-halo term of the power spectrum
       */
      double Pk_cs_numerator_integrand (const double mass, const double kk, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief model for the central-satellite part of the 1-halo
       *  term of the power spectrum
       *
       *  this function computes the central-satellite part of the
       *  1-halo term of power spectrum as follows (e.g. Harikane et
       *  al. 2016):
       *
       *  \f[P_{cs}(k, z) = \frac{2}{n_{gal}^2(z)}
       *  \int_{M_{min}}^{M_{max}} <N_{cen}N_{sat}>(M_h)\,n_h(M_h,
       *  z)\,\tilde{u}_h(k, M_h, z)\,{\rm d} M_h\f]
       *
       *  where the galaxy number density, \f$n_{gal}(z)\f$ is
       *  computed by cosmobl::modelling::twopt::ng,
       *  \f$<N_{cen}N_{sat}>(M_h)\f$ is computed by
       *  cosmobl::modelling::twopt::NcNs, the Fourier transform of
       *  the halo density profile, \f$\tilde{u}_h(k, M_h, z)\f$ if
       *  computed by
       *  cosmobl::cosmology::Cosmology::density_profile_FourierSpace
       *  and the halo mass function, \f$n_h(M_h, z)=dn/dM_h\f$, is
       *  computed by cosmology::Cosmology::mass_function; the
       *  integrand function of the numerator is computed by
       *  cosmobl::modelling::twopt::Pk_cs_numerator_integrand
       *
       *  @param kk the wave vector module at which the model is
       *  computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the central-satellite part of the 1-halo term of the
       *  power spectrum
       */
      double Pk_cs (const double kk, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief the integrand function to compute the numerator of
       *  the satellite-satellite part of the 1-halo term of the power
       *  spectrum
       *
       *  this function is the integrand to compute the
       *  satellite-satellite part of the 1-halo term of power
       *  spectrum as follows (e.g. Harikane et al. 2016); it returns
       *  the following quantity:
       *
       *  \f[<N_{sat}(N_{sat}-1)>(M_h)\,n_h(M_h, z)\,\tilde{u}_h^2(k,
       *  M_h, z)\f]
       *
       *  where the galaxy number density, \f$n_{gal}(z)\f$ is
       *  computed by cosmobl::modelling::twopt::ng,
       *  \f$<N_{sat}(N_{sat}-1)>(M_h)\f$ is computed by
       *  cosmobl::modelling::twopt::NsNs1, the Fourier transform of
       *  the halo density profile, \f$\tilde{u}_h(k, M_h, z)\f$ if
       *  computed by
       *  cosmobl::cosmology::Cosmology::density_profile_FourierSpace
       *  and the halo mass function, \f$n_h(M_h, z)=dn/dM_h\f$, is
       *  computed by cosmology::Cosmology::mass_function
       *
       *  @param mass \f$M_{h}\f$: the halo mass
       *
       *  @param kk the wave vector module at which the model is
       *  computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the integrand function to compute the
       *  satellite-satellite part of the 1-halo term of the power
       *  spectrum
       */
      double Pk_ss_numerator_integrand (const double mass, const double kk, const shared_ptr<void> inputs, vector<double> &parameter);
      
      /**
       *  @brief model for the satellite-satellite part of the 1-halo
       *  term of the power spectrum
       *
       *  this function computes the satellite-satellite part of the
       *  1-halo term of power spectrum as follows (e.g. Harikane et
       *  al. 2016):
       *
       *  \f[P_{ss}(k, z) = \frac{1}{n_{gal}^2(z)}
       *  \int_{M_{min}}^{M_{max}} <N_{sat}(N_{sat}-1)>(M_h)\,n_h(M_h,
       *  z)\,\tilde{u}_h^2(k, M_h, z)\,{\rm d} M_h\f]
       *
       *  where the galaxy number density, \f$n_{gal}(z)\f$ is
       *  computed by cosmobl::modelling::twopt::ng,
       *  \f$<N_{sat}(N_{sat}-1)>(M_h)\f$ is computed by
       *  cosmobl::modelling::twopt::NsNs1, the Fourier transform of
       *  the halo density profile, \f$\tilde{u}_h(k, M_h, z)\f$ is
       *  computed by
       *  cosmobl::cosmology::Cosmology::density_profile_FourierSpace
       *  and the halo mass function, \f$n_h(M_h, z)=dn/dM_h\f$, is
       *  computed by cosmology::Cosmology::mass_function; the
       *  integrand function of the numerator is computed by
       *  cosmobl::modelling::twopt::Pk_ss_numerator_integrand
       *
       *  @param kk the wave vector module at which the model is
       *  computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the satellite-satellite part of the 1-halo term of
       *  the power spectrum
       */
      double Pk_ss (const double kk, const shared_ptr<void> inputs, vector<double> &parameter);
      
      /**
       *  @brief model for the 1-halo term of the
       *  power spectrum
       *
       *  this function computes the 1-halo term of power spectrum as
       *  follows:
       *
       *  \f[P_{1halo}(k, z) = P_{cs}(k, z)+P_{ss}(k, z)\f]
       *
       *  where the central-satellite term, \f$P_{cs}(k, z)\f$, and
       *  the satellite-satellite term, \f$P_{ss}(k, z)\f$ are
       *  computed by cosmobl::modelling::twopt::Pk_cs and
       *  cosmobl::modelling::twopt::Pk_ss, respectively 
       *
       *  @param kk the wave vector module at which the model is
       *  computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the 1-halo term of the power spectrum
       */
      double Pk_1halo (const double kk, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief model for the 2-halo term of the power spectrum
       *
       *  this function computes the 2-halo term of the power spectrum
       *  as follows (e.g. Harikane et al. 2016):
       *
       *  \f[P_{2halo}(k, z) = P_m(k, z) \left[\frac{1}{n_{gal}(z)}
       *  \int_{M_{min}}^{M_{max}} N_{gal}(M_h)\,n_h(M_h, z)\,b_h(M,
       *  z)\,\tilde{u}_h(k, M_h, z)\,{\rm d} M_h\right]^2\f]
       *
       *  where the matter power spectrum, \f$P_m(k, z)\f$, is
       *  computed by cosmobl::cosmology::Cosmology::Pk, the average
       *  number of galaxies hosted in a dark matter halo of a given
       *  mass, \f$N_{gal}(M_h)\f$, is computed by
       *  cosmobl::modelling::twopt::Navg, the halo mass function,
       *  \f$n_h(M_h, z)=dn/dM_h\f$ is computed by
       *  cosmology::Cosmology::mass_function, the halo bias,
       *  \f$b(M_h, z)\f$, is computed by
       *  cosmology::Cosmology::bias_halo, and the Fourier transform
       *  of the density profile, \f$\tilde{u}_h(k, M_h, z)\f$, is
       *  computed by
       *  cosmobl::cosmology::Cosmology::density_profile_FourierSpace
       *
       *  @param kk the wave vector module at which the model is
       *  computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the 2-halo term of the power spectrum
       */
      double Pk_2halo (const double kk, const shared_ptr<void> inputs, vector<double> &parameter);
	
      /**
       *  @brief HOD model of the power spectrum
       *
       *  this function computes the power spectrum as follows:
       *
       *  \f[P(k, z) = P_{1halo}(k, z)+P_{2halo}(k, z)\f]
       *
       *  where the 1-halo and 2-halo terms of the power spectrum,
       *  \f$P_{1halo}(k, z)\f$ and \f$P_{2halo}(k, z)\f$ are computed
       *  by cosmobl::modelling::twopt::Pk_1halo and
       *  cosmobl::modelling::twopt::Pk_2halo, respectively
       *
       *  @param kk the wave vector module at which the model is
       *  computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the 1-halo term of the power spectrum
       */
      double Pk_HOD (const double kk, const shared_ptr<void> inputs, vector<double> &parameter);
      
      /**
       *  @brief model for the 1-halo term of the monopole of the
       *  two-point correlation function
       *
       *  this function computes the 1-halo term of the two-point
       *  correlation function by Fourier transforming the 1-halo term
       *  of the power spectrum:
       *
       *  \f[\xi_{1halo}(r, z) =
       *  \frac{1}{2\pi^2}\int_{0}^{k_{max}}{\rm d}k\, k^2P_{1halo}(k,
       *  z)\frac{\sin (kr)}{kr}\f]
       *
       *  where the 1-halo term of the power spectrum \f$P_{1halo}(k,
       *  z)\f$, is computed by cosmobl::modelling::twopt::Pk_1halo
       *
       *  so the integral that is actually computed is the following:
       *
       *  \f[\xi_{1halo}(r, z) =
       *  \frac{1}{2\pi^2}\int_{0}^{k_{max}}{\rm d}k\, \left[P_{cs}(k,
       *  z)+P_{ss}(k, z)\right]\frac{k\sin (kr)}{r} = \f]
       *
       *  \f[= \frac{1}{2\pi^2n_{gal}^2(z)}\int_{0}^{k_{max}}{\rm
       *  d}k\, \int_{M_{min}}^{M_{max}}{\rm d} M_h\, \left[
       *  2<N_{cen}N_{sat}>(M_h)\,n_h(M_h, z)\,\tilde{u}_h(k, M_h, z)
       *  + <N_{sat}(N_{sat}-1)>(M_h)\,n_h(M_h, z)\,\tilde{u}_h^2(k,
       *  M_h, z)\right]\frac{k\sin (kr)}{r} \f]
       *
       *  where the galaxy number density, \f$n_{gal}(z)\f$ is
       *  computed by cosmobl::modelling::twopt::ng,
       *  \f$<N_{cen}N_{sat}>(M_h)\f$ is computed by
       *  cosmobl::modelling::twopt::NcNs,
       *  \f$<N_{sat}(N_{sat}-1)>(M_h)\f$ is computed by
       *  cosmobl::modelling::twopt::NsNs1, the Fourier transform of
       *  the halo density profile, \f$\tilde{u}_h(k, M_h, z)\f$ is
       *  computed by
       *  cosmobl::cosmology::Cosmology::density_profile_FourierSpace
       *  and the halo mass function, \f$n_h(M_h, z)=dn/dM_h\f$, is
       *  computed by cosmology::Cosmology::mass_function
       *
       *  the multidimensional integral is implemented with the CUBA
       *  libraries
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the 1-halo term of the monopole of the two-point
       *  correlation function
       */
      vector<double> xi_1halo (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief auxiliary function used to speed up the computation
       *  of the multidimensional integral of the 2-halo term of the
       *  monopole of the two-point correlation function
       *
       *  this function provides a shared pointer to the following
       *  function defined on a grid:
       *
       *  \f[ N_{gal}(M_h)\,n_h(M_h, z)\,b_h(M, z)\,\tilde{u}_h^2(k,
       *  M_h, z) \f]
       *
       *  where the matter power spectrum, \f$P_m(k, z)\f$, is
       *  computed by cosmobl::cosmology::Cosmology::Pk, the average
       *  number of galaxies hosted in a dark matter halo of a given
       *  mass, \f$N_{gal}(M_h)\f$, is computed by
       *  cosmobl::modelling::twopt::Navg, the halo mass function,
       *  \f$n_h(M_h, z)=dn/dM_h\f$ is computed by
       *  cosmology::Cosmology::mass_function, the halo bias,
       *  \f$b(M_h, z)\f$, is computed by
       *  cosmology::Cosmology::bias_halo, and the Fourier transform
       *  of the density profile, \f$\tilde{u}_h^2(k, M_h, z)\f$, is
       *  computed by
       *  cosmobl::cosmology::Cosmology::density_profile_FourierSpace
       *
       *  @param kk vector containting the wave vector modules at which
       *  the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return a shared pointer to the function specified above
       */
      shared_ptr<glob::FuncGrid> func_2halo (const vector<double> kk, const shared_ptr<void> inputs, vector<double> &parameter);
      
      /**
       *  @brief model for the 2-halo term of the monopole of the
       *  two-point correlation function
       *
       *  this function computes the 2-halo term of the two-point
       *  correlation function by Fourier transforming the 2-halo term
       *  of the power spectrum:
       *
       *  \f[\xi_{2halo}(r, z) =
       *  \frac{1}{2\pi^2}\int_{0}^{k_{max}}k^2P_{2halo}(k,
       *  z)\frac{\sin (kr)}{kr}\,{\rm d}k\f]
       *
       *  where the 2-halo term of the power spectrum \f$P_{2halo}(k,
       *  z)\f$, is computed by cosmobl::modelling::twopt::Pk_2halo
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the 2-halo term of the monopole of the two-point
       *  correlation function
       */
      vector<double> xi_2halo (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);
      
      /**
       *  @brief HOD model of the monopole of the two-point
       *  correlation function
       *
       *  the function computes:
       *
       *  \f[\xi(r, z) = \xi_{1halo}(r, z)+\xi_{2halo}(r, z)\f]
       *
       *  where \f$\xi_{1halo}(r, z)\f$ and \f$\xi_{2halo}(r, z)\f$
       *  are computed by cosmobl::modelling::twopt::xi_1halo and
       *  cosmobl::modelling::twopt::xi_2halo, respectively
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the HOD monopole of the two-point correlation
       *  function
       */
      vector<double> xi_HOD (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief function used to compute the redshift-space monopole
       *  of the two-point correlation function
       *
       *  this function is used to compute the redshift-space
       *  spherically averaged correlation function with a modified
       *  version of the Kaiser model (sec. 2.3 of van den Bosch et
       *  al. 2012); it is used to derive the model projected
       *  correlation function accounting for residual redshift-space
       *  distortions caused by finite integration range
       *
       *  \f[\xi(r_p, \pi, z) = \sum_{l=0}^2\xi_{2l}(r,
       *  z)\mathcal{P}_{2l}(\mu)\f]
       *
       *  where
       *
       *  \f[\xi_0(r, z) = \left(1 + \frac{2}{3}\beta +
       *  \frac{1}{5}\beta^2\right)\xi(r, z)\,,\f]
       *
       *  \f[\xi_2(r, z) = \left(\frac{2}{3}\beta +
       *  \frac{4}{7}\beta^2\right)\left[\xi(r, z)-3J_3(r,
       *  z)\right]\,,\f]
       *
       *  \f[\xi_4(r, z) = \frac{8}{35}\beta^2\left[\xi(r,
       *  z)+\frac{15}{2}J_3(r, z)-\frac{35}{2}J_5(r, z)\right]\,,\f]
       *
       *  \f[J_n(r, z) = \frac{1}{r^n}\int_0^r\xi(y,
       *  z)y^{n-1}{\rm d}y\f]
       *
       *  the real-space galaxy correlation function \f$\xi(r, z)\f$
       *  is either the 1-halo, the 2-halo, or the full-shape
       *  real-space correlation function, computed by either
       *  cosmobl::modelling::twopt::xi_1halo,
       *  cosmobl::modelling::twopt::xi_2halo, or
       *  cosmobl::modelling::twopt::xi_HOD, respectively
       *
       *  @param func the two-point correlation function that will be
       *  integrated (it can be either the 1-halo, the 2-halo, or the
       *  full-shape correlation function)
       *
       *  @param rp the scale perpendicular to the line of sight at
       *  which the model is computed
       *
       *  @param pi the scale parallel to the line of sight at which
       *  the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the redshift-space monopole of the two-point
       *  correlation function
       */
      double xi_zspace (function<vector<double>(const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter)> func, const double rp, const double pi, const shared_ptr<void> inputs, vector<double> &parameter);
      
      /**
       *  @brief model for the 1-halo term of the redshift-space
       *  monopole of the two-point correlation function
       *
       *  this function computes the 1-halo term of redshift-space
       *  spherically averaged correlation function with a modified
       *  version of the Kaiser model (sec. 2.3 of van den Bosch et
       *  al. 2012); it is used to derive the model of the 1-halo term
       *  of the projected correlation function accounting for
       *  residual redshift-space distortions caused by finite
       *  integration range
       *
       *  \f[\xi_{1halo}(r_p, \pi, z) = \sum_{l=0}^2\xi_{2l}(r,
       *  z)\mathcal{P}_{2l}(\mu)\f]
       *
       *  where
       *
       *  \f[\xi_0(r, z) = \left(1 + \frac{2}{3}\beta +
       *  \frac{1}{5}\beta^2\right)\xi_{1halo}(r, z)\,,\f]
       *
       *  \f[\xi_2(r, z) = \left(\frac{4}{3}\beta +
       *  \frac{4}{7}\beta^2\right)\left[\xi_{1halo}(r, z)-3J_3(r,
       *  z)\right]\,,\f]
       *
       *  \f[\xi_4(r, z) = \frac{8}{35}\beta^2\left[\xi_{1halo}(r,
       *  z)+\frac{15}{2}J_3(r, z)-\frac{35}{2}J_5(r, z)\right]\,,\f]
       *
       *  \f[J_n(r, z) = \frac{1}{r^n}\int_0^r\xi_{1halo}(y,
       *  z)y^{n-1}{\rm d}y\f]
       *
       *  the real-space galaxy correlation function \f$\xi_{g,
       *  1halo}(r, z)\f$ is computed by
       *  cosmobl::modelling::twopt::xi_1halo
       *
       *  @param rp the scale perpendicular to the line of sight at
       *  which the model is computed
       *
       *  @param pi the scale parallel to the line of sight at which
       *  the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the 1-halo term of the redshift-space monopole of
       *  the two-point correlation function
       */
      double xi_1halo_zspace (const double rp, const double pi, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief model for the 2-halo term of the redshift-space
       *  monopole of the two-point correlation function
       *
       *  this function computes the 2-halo term of redshift-space
       *  spherically averaged correlation function with a modified
       *  version of the Kaiser model (sec. 2.3 of van den Bosch et
       *  al. 2012); it is used to derive the model of the 2-halo term
       *  of the projected correlation function accounting for
       *  residual redshift-space distortions caused by finite
       *  integration range
       *
       *  \f[\xi_{2halo}(r_p, \pi, z) = \sum_{l=0}^2\xi_{2l}(r,
       *  z)\mathcal{P}_{2l}(\mu)\f]
       *
       *  where
       *
       *  \f[\xi_0(r, z) = \left(1 + \frac{2}{3}\beta +
       *  \frac{1}{5}\beta^2\right)\xi_{2halo}(r, z)\,,\f]
       *
       *  \f[\xi_2(r, z) = \left(\frac{4}{3}\beta +
       *  \frac{4}{7}\beta^2\right)\left[\xi_{2halo}(r, z)-3J_3(r,
       *  z)\right]\,,\f]
       *
       *  \f[\xi_4(r, z) = \frac{8}{35}\beta^2\left[\xi_{2halo}(r,
       *  z)+\frac{15}{2}J_3(r, z)-\frac{35}{2}J_5(r, z)\right]\,,\f]
       *
       *  \f[J_n(r, z) = \frac{1}{r^n}\int_0^r\xi_{2halo}(y,
       *  z)y^{n-1}{\rm d}y\f]
       *
       *  the real-space galaxy correlation function \f$\xi_{g,
       *  2halo}(r, z)\f$ is computed by
       *  cosmobl::modelling::twopt::xi_2halo
       *
       *  @param rp the scale perpendicular to the line of sight at
       *  which the model is computed
       *
       *  @param pi the scale parallel to the line of sight at which
       *  the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the 2-halo term of redshift-space monopole of the
       *  two-point correlation function
       */
      double xi_2halo_zspace (const double rp, const double pi, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief HOD model of the redshift-space monopole of the
       *  two-point correlation function
       *
       *  this function computes the redshift-space spherically
       *  averaged correlation function with a modified version of the
       *  Kaiser model (sec. 2.3 of van den Bosch et al. 2012); it is
       *  used to derive the model of the projected correlation
       *  function accounting for residual redshift-space distortions
       *  caused by finite integration range
       *
       *  \f[\xi(r_p, \pi, z) = \xi_{1halo}(r_p, \pi,
       *  z)+\xi_{2halo}(r_p, \pi, z)\f]
       *
       *  where \f$\xi_{1halo}(r_p, \pi, z)\f$ and
       *  \f$\xi_{2halo}(r_p, \pi, z)\f$ are computed by
       *  cosmobl::modelling::twopt::xi_1halo_zspace and
       *  cosmobl::modelling::twopt::xi_2halo_zspace, respectively
       *
       *  @param rp the scale perpendicular to the line of sight at
       *  which the model is computed
       *
       *  @param pi the scale parallel to the line of sight at which
       *  the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the HOD redshift-space monopole of the two-point
       *  correlation function
       */
      double xi_HOD_zspace (const double rp, const double pi, const shared_ptr<void> inputs, vector<double> &parameter);
	
      ///@}

      vector<double> xi0_linear_cosmology_clusters_selection_function (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter);
    }
  }
}

#endif
