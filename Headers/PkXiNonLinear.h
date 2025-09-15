/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *******************************************************************/

/**
 *  @file Headers/PkXiNonLinear.h
 *
 *  @brief The class PkXiNonLinear
 *
 *  This file defines the interface of the class PkXiNonLinear, used
 *  to model nonlinearities in two-point correlation function and
 *  power spectrum
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __PkXiNonLinear__
#define __PkXiNonLinear__

#include "PkXi.h"


// ===================================================================================================


namespace cbl {
  
  namespace cosmology {

    /**
     *  @class PkXiNonLinear PkXiNonLinear.h "Headers/PkXiNonLinear.h"
     *
     *  @brief The class PkXiNonLinear
     *
     *  This class is used to handle objects of type <EM>
     *  PkXiNonLinear </EM>.
     */
    class PkXiNonLinear {

    private:

      /// pointer to the input cosmology
      std::shared_ptr<cosmology::Cosmology> m_cosmology = NULL;
      
    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      PkXiNonLinear () = default;
      
      /**
       *  @brief constructor
       *
       *  @param cosmology pointer to an object of class Cosmology
       */
      PkXiNonLinear (std::shared_ptr<cbl::cosmology::Cosmology> cosmology)
      : m_cosmology(std::move(cosmology)) {}
      
      /**
       *  @brief default destructor
       */
      ~PkXiNonLinear () = default;

      ///@}

      
      /**
       *  @name Functions to get the private members of the class
       */
      ///@{
      
      /**
       *  @brief Get the private member m_cosmology
       *
       *  @return the cosmological model
       */
      std::shared_ptr<cosmology::Cosmology> cosmology ()
      { return move(m_cosmology); }

      ///@}

      
      /**
       *  @name Functions used to model nonlinearities in two-point
       *  correlation function and power spectrum
       */
      ///@{

      /**
       * @brief function used to estimate the non-linear power
       * spectrum
       *
       * this function computes the following integral, which is used
       * to estimate the non-linear power spectrum estimate the
       * non-linear power spectrum following the MPTbreeze scheme
       * (Crocce et al. 2012, https://arxiv.org/pdf/1207.1465.pdf):
       *
       * \f[ f(k) = \int_{q_{min}}^{q_{max}} \mathrm{d}q\,
       * q^2P(k)\frac{1}{504k^3q^5}\left[
       * 6k^7q-79k^5q^3+50k^3q^5-21kq^7+\frac{3}{4}\left(k^2-q^2\right)^3
       * \left(2k^2+7q^2\right)\ln\left(\frac{(k-q)^2}{(k+q)^2}\right)
       * \right] \f]
       *
       * @param k the wavevector module
       *
       * @param PkLin linear power spectrum
       *
       * @param qmin
       *
       * @param qmax
       *
       * @param prec the integral precision
       *
       * @return \f$f(k)\f$
       */
      double f_k (const double k, const std::shared_ptr<cbl::glob::FuncGrid> PkLin, const double qmin, const double qmax, const double prec=1.e-3);

      /**
       * @brief function used to estimate the non-linear power
       * spectrum
       *
       * this function computes the following integral, which is used
       * to estimate the non-linear power spectrum estimate the
       * non-linear power spectrum following the MPTbreeze scheme
       * (Crocce et al. 2012, https://arxiv.org/pdf/1207.1465.pdf):
       *
       * \f[ f(k) = \int_{q_{min}}^{q_{max}} \mathrm{d}q\,
       * q^2P(k)\frac{1}{168k^3q^5}\left[
       * 6k^7q-41k^5q^3+2k^3q^5-3kq^7+\frac{3}{4}\left(k^2-q^2\right)^3
       * \left(2k^2+q^2\right)\ln\left(\frac{(k-q)^2}{(k+q)^2}\right)
       * \right] \f]
       *
       * @param k the wavevector module
       *
       * @param PkLin linear power spectrum
       *
       * @param qmin
       *
       * @param qmax
       *
       * @param prec the integral precision
       *
       * @return \f$g(k)\f$
       */
      double g_k (const double k, const std::shared_ptr<cbl::glob::FuncGrid> PkLin, const double qmin, const double qmax, const double prec=1.e-3);

      /**
       * @brief function used to estimate the non-linear power
       * spectrum
       *
       * this function computes the following equation, which is used
       * to estimate the non-linear power spectrum estimate the
       * non-linear power spectrum following the MPTbreeze scheme
       * (Crocce et al. 2012, https://arxiv.org/pdf/1207.1465.pdf):
       *
       * \f[ F_2 = \frac{5}{7} + \frac{kq}{2}\left(\frac{k}{q} +
       * \frac{q}{k}\right) + \frac{2}{7}(kq)^2 \f]
       *
       * @param k the wavevector module
       *
       * @param q
       *
       * @param kq
       *
       * @return \f$F_2\f$
       */
      double F2 (const double k, const double q, const double kq);

      /**
       * @brief function used to estimate the non-linear power
       * spectrum
       *
       * this function computes the following equation, which is used
       * to estimate the non-linear power spectrum estimate the
       * non-linear power spectrum following the MPTbreeze scheme
       * (Crocce et al. 2012, https://arxiv.org/pdf/1207.1465.pdf):
       *
       * \f[ G_2 = \frac{3}{7} + \frac{kq}{2}\left(\frac{k}{q} +
       * \frac{q}{k}\right) + \frac{4}{7}(kq)^2 \f]
       *
       * @param k the wavevector module
       *
       * @param q
       *
       * @param kq
       *
       * @return \f$ G_2 \f$
       */
      double G2 (const double k, const double q, const double kq);
      
      /**
       * @brief the one-loop power spectrum
       *
       *  this function computes the one-loop power spectrum
       *  corrections (Crocce et al. 2012,
       *  https://arxiv.org/pdf/1207.1465.pdf):
       *
       *  \f[ P^{1loop}(k; f_1, f_2) = 4\pi\int_{q_{min}}^{q_{max}}
       *  \mathrm{d} q\, q^2P(q) \int_{-1}^{1} \mathrm{d} x\, f_1(kq,
       *  q, akq) \, f_2(kq, q, akq) \, P(kq) \f]
       *
       *  where 
       *
       *  \f[ kq = \sqrt{k^2+q^2-2kqx} \f]
       * 
       *  \f[ akq = \frac{kx-q}{kq} \f]
       *
       *  and \f$f_1\f$, \f$f_2\f$ are computed by either
       *  cbl::cosmology::PkXiNonLinear::F2 or
       *  cbl::cosmology::PkXiNonLinear::G2
       *
       *  @param kk the wavevector module
       *
       *  @param PkLin linear power spectrum
       *
       *  @param corrtype 0 \f$\rightarrow\f$ \f$f_1=f_2=F_2\f$; 1
       *  \f$\rightarrow\f$ \f$f_1=F_2\f$, \f$f_2=G_2\f$; 2
       *  \f$\rightarrow\f$ \f$f_1=f_2=G_2\f$
       *
       *  @param qmin minimum q value in the integration
       *
       *  @param qmax maximum q value in the integration
       *
       *  @param prec the integral precision
       *
       *  @return one-loop power spectrum
       */
      double Pk_1loop (const double kk, const std::shared_ptr<cbl::glob::FuncGrid> PkLin, const int corrtype, const double qmin, const double qmax, const double prec=1.e-3);

      /**
       *  @brief the real-space matter non-linear power spectrum
       *  \f$P_{\delta\delta}(k)\f$, computed at 1-loop
       *
       *  this function computes the real-space non-linear matter
       *  power spectrum \f$P_{\delta\delta}(k)\f$ at 1-loop,
       *  following the MPTbreeze scheme (Crocce et al. 2012,
       *  https://arxiv.org/pdf/1207.1465.pdf):
       *
       *  \f[ P_{\delta\delta}(k) =
       *     \left(2\pi\right)^3\left(\exp(f(k)\right)^2
       *     (P^{lin}(k)+P^{1loop}(k; F_2, F_2)) \f]
       *
       *  where \f$P^{lin}(k)\f$ is the linear power spectrum,
       *  \f$f(k)\f$ is the second-order correction of the non-linear
       *  propagator computed by cbl::cosmology::PkXiNonLinear::f_k, and
       *  \f$P^{1loop}\f$ is the one-loop power spectrum correction,
       *  computed by cbl::cosmology::PkXiNonLinear::Pk_1loop
       *
       *  @param kk the wavevector module
       *
       *  @param Pk pointer to a FuncGrid object to interpolate the
       *  linear power spectrum
       *
       *  @param qmin the lower integration limit
       *
       *  @param qmax the upper integration limit
       *
       *  @param prec the integral precision
       *
       *  @return \f$P_{\delta\delta}(k)\f$
       */
      double Pk_DeltaDelta (const double kk, const std::shared_ptr<cbl::glob::FuncGrid> Pk, const double qmin, const double qmax, const double prec=1.e-3);

      /**
       *  @brief the real-space matter non-linear power spectrum
       *  \f$P_{\delta\delta}(k)\f$, computed at 1-loop
       *
       *  this function computes the real-space non-linear matter
       *  power spectrum \f$P_{\delta\delta}(k)\f$ at 1-loop,
       *  following the MPTbreeze scheme (Crocce et al. 2012,
       *  https://arxiv.org/pdf/1207.1465.pdf):
       *
       *  \f[ P_{\delta\delta}(k) =
       *     \left(2\pi\right)^3\left(\exp(f(k)\right)^2
       *     (P^{lin}(k)+P^{1loop}(k; F_2, F_2)) \f]
       *
       *  where \f$P^{lin}(k)\f$ is the linear power spectrum,
       *  \f$f(k)\f$ is the second-order correction of the non-linear
       *  propagator computed by cbl::cosmology::PkXiNonLinear::f_k, and
       *  \f$P^{1loop}\f$ is the one-loop power spectrum correction,
       *  computed by cbl::cosmology::PkXiNonLinear::Pk_1loop
       *
       *  @param kk the wavevector module
       *
       *  @param redshift the redshift
       * 
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If this parameter is different from the default value
       *  it will be used also in the output directory name
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
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will use be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return \f$P_{\delta\delta}(k)\f$
       */
      std::vector<double> Pk_DeltaDelta (const std::vector<double> kk, const double redshift, const std::string method_Pk, const bool store_output=true, const std::string output_root="test", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString, const bool unit1=false);
      
      /**
       *  @brief the real-space matter non-linear power spectrum
       *  \f$P_{\delta\theta}(k)\f$, computed at 1-loop
       *
       *  this function computes the real-space non-linear matter
       *  power spectrum \f$P_{\delta\delta}(k)\f$ at 1-loop,
       *  following the MPTbreeze scheme (Crocce et al. 2012,
       *  https://arxiv.org/pdf/1207.1465.pdf):
       *
       *  \f[ P_{\delta\theta}(k) =
       *     \left(2\pi\right)^3\left(\exp(f(k)\right)^2
       *     (P^{lin}(k)+P^{1loop}(k; F_2, G_2)) \f]
       *
       *  where \f$P^{lin}(k)\f$ is the linear power spectrum,
       *  \f$f(k)\f$ is the second-order correction of the non-linear
       *  propagator computed by cbl::cosmology::PkXiNonLinear::f_k, and
       *  \f$P^{1loop}\f$ is the one-loop power spectrum correction,
       *  computed by cbl::cosmology::PkXiNonLinear::Pk_1loop
       *
       *  @param kk the wavevector module
       *
       *  @param Pk pointer to a FuncGrid object to interpolate
       *  the linear power spectrum
       *
       *  @param qmin the lower integration limit
       *
       *  @param qmax the upper integration limit
       *
       *  @param prec the integral precision
       *
       *  @return \f$P_{\delta\theta}(k)\f$
       */
      double Pk_DeltaTheta (const double kk, const std::shared_ptr<cbl::glob::FuncGrid> Pk, const double qmin, const double qmax, const double prec=1.e-3);

      /**
       *  @brief the real-space matter non-linear power spectrum
       *  \f$P_{\delta\theta}(k)\f$, computed at 1-loop
       *
       *  this function computes the real-space non-linear matter
       *  power spectrum \f$P_{\delta\delta}(k)\f$ at 1-loop,
       *  following the MPTbreeze scheme (Crocce et al. 2012,
       *  https://arxiv.org/pdf/1207.1465.pdf):
       *
       *  \f[ P_{\delta\theta}(k) =
       *     \left(2\pi\right)^3\left(\exp(f(k)\right)^2
       *     (P^{lin}(k)+P^{1loop}(k; F_2, G_2)) \f]
       *
       *  where \f$P^{lin}(k)\f$ is the linear power spectrum,
       *  \f$f(k)\f$ is the second-order correction of the non-linear
       *  propagator computed by cbl::cosmology::PkXiNonLinear::f_k, and
       *  \f$P^{1loop}\f$ is the one-loop power spectrum correction,
       *  computed by cbl::cosmology::PkXiNonLinear::Pk_1loop
       *
       *  @param kk the wavevector module
       *
       *  @param redshift the redshift
       * 
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If this parameter is different from the default value
       *  it will be used also in the output directory name
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
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will use be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return \f$P_{\delta\theta}(k)\f$
       */
      std::vector<double> Pk_DeltaTheta (const std::vector<double> kk, const double redshift, const std::string method_Pk, const bool store_output=true, const std::string output_root="test", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString, const bool unit1=false);
      
      /**
       *  @brief the real-space matter non-linear power spectrum
       *  \f$P_{\theta\theta}(k)\f$, computed at 1-loop
       *
       *  this function computes the real-space non-linear matter
       *  power spectrum \f$P_{\delta\delta}(k)\f$ at 1-loop,
       *  following the MPTbreeze scheme (Crocce et al. 2012,
       *  https://arxiv.org/pdf/1207.1465.pdf):
       *
       *  \f[ P_{\theta\theta} (k) = \left(\exp(g(k)\right)^2
       *     (P^{lin}(k)+P^{1loop}(k; G_2, G_2)) \f]
       *
       *  where \f$P^{lin}(k)\f$ is the linear power spectrum,
       *  \f$f(k)\f$ is the second-order correction of the non-linear
       *  propagator computed by cbl::cosmology::PkXiNonLinear::g_k, and
       *  \f$P^{1loop}\f$ is the one-loop power spectrum correction,
       *  computed by cbl::cosmology::PkXiNonLinear::Pk_1loop 
       *
       *  @param kk the wavevector module
       *
       *  @param Pk pointer to a FuncGrid object to interpolate the
       *  linear power spectrum
       *
       *  @param qmin the lower integration limit
       * 
       *  @param qmax the upper integration limit
       * 
       *  @param prec the integral precision
       * 
       *  @return \f$P_{\theta\theta}(k)\f$
       */
      double Pk_ThetaTheta (const double kk, const std::shared_ptr<cbl::glob::FuncGrid> Pk, const double qmin, const double qmax, const double prec=1.e-3);

      /**
       *  @brief the real-space matter non-linear power spectrum
       *  \f$P_{\theta\theta}(k)\f$, computed at 1-loop
       *
       *  this function computes the real-space non-linear matter
       *  power spectrum \f$P_{\delta\delta}(k)\f$ at 1-loop,
       *  following the MPTbreeze scheme (Crocce et al. 2012,
       *  https://arxiv.org/pdf/1207.1465.pdf):
       *
       *  \f[ P_{\theta \theta}(k) =
       *     \left(2\pi\right)^3\left(\exp(g(k)\right)^2
       *     (P^{lin}(k)+P^{1loop}(k; G_2, G_2)) \f]
       *
       *  where \f$P^{lin}(k)\f$ is the linear power spectrum,
       *  \f$g(k)\f$ is the second-order correction of the non-linear
       *  propagator computed by cbl::cosmology::PkXiNonLinear::g_k,
       *  and \f$P^{1loop}\f$ is the one-loop power spectrum
       *  correction, computed by
       *  cbl::cosmology::PkXiNonLinear::Pk_1loop
       *
       *  @param kk vector of wavevector modules
       *
       *  @param redshift the redshift
       * 
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If this parameter is different from the default value
       *  it will be used also in the output directory name
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
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will use be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return \f$P_{\theta\theta}(k)\f$
       */
      std::vector<double> Pk_ThetaTheta (const std::vector<double> kk, const double redshift, const std::string method_Pk, const bool store_output=true, const std::string output_root="test", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString, const bool unit1=false);
      
      /**
       *  @brief the multipoles of the A and B correction terms for
       *  the TNS model 
       *
       *  the multipoles of the A and B correction terms for the TNS
       *  model are computed at 1-loop using the Standard
       *  Perturbation Theory implemented in the CPT Library
       *  [http://www2.yukawa.kyoto-u.ac.jp/~atsushi.taruya/cpt_pack.html].
       *  Details can be found in Taruya et al. (2010)
       *  [https://arxiv.org/abs/1006.0699].
       *
       *  \f[ A(k,\mu) = A_0(k)L_0(\mu) + A_2(k)L_2(\mu) +
       *  A_4(k)L_4(\mu) \f]
       *
       *  \f[ B(k,\mu) = B_0(k)L_0(\mu) + B_2(k)L_2(\mu) +
       *  B_4(k)L_4(\mu) \f]
       *
       *  where \f$L_l\f$ represents the Legendre polynomial of order
       *  \f$l\f$.
       *
       *  @author J.E. Garcia-Farieta
       *  @author joegarciafa@unal.edu.co
       *
       *  @param kk the wave vector module
       *
       *  @param method method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  CLASS [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created
       *  by CAMB are stored; if false the output files created by
       *  CAMB are removed
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If this parameter is different from the default value
       *  it will be used also in the output directory name
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
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param prec accuracy of the integration
       *
       *  @return multipoles of A and B terms for TNS model (A0, A2, A4, B0, B2, B4)
       */
      std::vector<std::vector<double>> Pk_TNS_AB_multipoles (std::vector<double> kk, const std::string method, const double redshift, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec);

      /**
       *  @brief the A and B correction terms for the TNS model at
       *  1-loop from the multipole expansion
       *
       *  the A and B correction terms for the TNS model at 1-loop
       *  from the multipole expansion are computed using the
       *  Standard Perturbation Theory implemented in the CPT Library
       *  [http://www2.yukawa.kyoto-u.ac.jp/~atsushi.taruya/cpt_pack.html].
       *  Details can be found in Taruya et al. (2010)
       *  [https://arxiv.org/abs/1006.0699].
       *
       *  \f[ A(k,\mu) = A_0(k)L_0(\mu) + A_2(k)L_2(\mu) +
       *  A_4(k)L_4(\mu) \f]
       *
       *  \f[ B(k,\mu) = B_0(k)L_0(\mu) + B_2(k)L_2(\mu) +
       *  B_4(k)L_4(\mu) \f]
       *
       *  where \f$L_l\f$ represents the Legendre polynomial of order
       *  \f$l\f$.
       *
       *  @author J.E. Garcia-Farieta
       *  @author joegarciafa@unal.edu.co
       *
       *  @param kk the wave vector module
       *
       *  @param mu the cosine of the angle between the separation
       *  and the line of sight
       *
       *  @param method method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  CLASS [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created
       *  by CAMB are stored; if false the output files created by
       *  CAMB are removed
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If this parameter is different from the default value
       *  it will be used also in the output directory name
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
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration
       *
       *  @return A and B terms (total contribution) from multipoles
       */
      std::vector<std::vector<double>> Pk_TNS_AB_1loop (std::vector<double> kk, const double mu, const std::string method, const double redshift, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec);

      /**
       *  @brief the expanded A and B correction terms for the TNS
       *  model
       *
       *  the expanded A and B correction terms for the TNS model are
       *  computed at 1-loop using the Standard Perturbation Theory
       *  implemented in the CPT Library
       *  [http://www2.yukawa.kyoto-u.ac.jp/~atsushi.taruya/cpt_pack.html].
       *  Details can be found in Taruya et al. (2010)
       *  [https://arxiv.org/abs/1006.0699].
       *
       *  \f[ A(k,\mu,\beta) = \beta\mu^2A_{11} + \beta^2(\mu^2A_{12}
       *  + \mu^4A_{22}) + \beta^3(\mu^4A_{23} + \mu^6A_{33}) \f] 
       *
       *  \f[ B(k,\mu,\beta) = \mu^2(\beta^2B_{12} + \beta^3B_{13} +
       *  \beta^4B_{14}) + \mu^4(\beta^2B_{22} + \beta^3B_{23} +
       *  \beta^4B_{24}) + \mu^6(\beta^3B_{33} + \beta^4B_{34}) +
       *  \mu^8\beta^4B_{44} \f]
       *
       *  @author J.E. Garcia-Farieta
       *  @author joegarciafa@unal.edu.co
       *
       *  @param kk the wave vector module
       *
       *  @param method method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  CLASS [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created
       *  by CAMB are stored; if false the output files created by
       *  CAMB are removed
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can
       *  be any name. If this parameter is different from the default value
       *  it will be used also in the output directory name
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
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @return A and B terms (expanded) for TNS model: A11, A12,
       *  A22, A23, A33, B12, B13, B14, B22, B23, B24, B33, B34, B44
       */
      std::vector<std::vector<double>> Pk_TNS_AB_terms_1loop (std::vector<double> kk, const std::string method, const double redshift, const bool store_output, const std::string output_root, const int norm, const double k_min=0.001, const double k_max=100., const double prec=1.e-2);

      /**
       *  @brief The expanded correction terms for the extended TNS
       *  model (eTNS)
       *
       *  the expanded correction terms for the eTNS model are
       *  computed at 1-loop using the Standard Perturbation Theory
       *  implemented in the CAMB-integrated version by Shun Saito.
       *  [http://www2.yukawa.kyoto-u.ac.jp/~atsushi.taruya/cpt_pack.html].
       *  Details can be found in Saito et al. (2014) and Beutler et
       *  al. (2014) [https://arxiv.org/abs/1405.1447 and
       *  https://arxiv.org/abs/1312.4611v2].
       *
       *  \f[ P(k, \mu) = D_{FoG}(k', \mu', f,
       *  \sigma_v)\left[P_{\mathrm{g}, \delta \delta}(k) +2 f
       *  \mu^{2} P_{\mathrm{g}, \delta \theta}(k)+f^{2} \mu^{4}
       *  P_{\theta \theta}(k) + b_{1}^{3} A(k, \mu, \beta)+b_{1}^{4}
       *  B(k, \mu, \beta)\right] \f]
       *
       *  where
       *
       *  \f[ P_{\mathrm{g}, \delta \delta}(k) = b_{1}^{2} P_{\delta
       *  \delta}(k)+2 b_{2} b_{1} P_{b 2, \delta}(k)+2 b_{s 2} b_{1}
       *  P_{b s 2, \delta}(k) \\ +2 b_{3 \mathrm{nl}} b_{1}
       *  \sigma_{3}^{2}(k) P_{\mathrm{m}}^{\mathrm{lin}}(k)+b_{2}^{2}
       *  P_{b 22}(k) \\ +2 b_{2} b_{s 2} P_{b 2 s 2}(k)+b_{s 2}^{2}
       *  P_{b s 22}(k)+N \f]
       *
       *  \f[ P_{\mathrm{g}, \delta \theta}(k) = b_{1} P_{\delta
       *  \theta}(k)+b_{2} P_{b 2, \theta}(k)+b_{s 2} P_{b s 2,
       *  \theta}(k) \\ +b_{3 \mathrm{nl}} \sigma_{3}^{2}(k)
       *  P_{\mathrm{m}}^{\mathrm{lin}}(k) \f]
       *
       *  @author J.E. Garcia-Farieta
       *  @author joegarciafa@unal.edu.co
       *
       *  @param kk the wave vector module
       *
       *  @param method method used to compute the linear power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  CLASS [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created
       *  by CAMB are stored; if false the output files created by
       *  CAMB are removed
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can
       *  be any name. If this parameter is different from the default value
       *  it will be used also in the output directory name
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
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @return expanded terms of the eTNS model: Pdd, Pdv, Pvv,
       *  Pb2d, Pb2v, Pb22, Pbs2d, Pbs2v, Pb2s2, Pbs22, sigma32Pklin,
       *  Bb1, Bb2, Bbs2
       */
      std::vector<std::vector<double>> Pk_eTNS_terms_1loop (std::vector<double> kk, const std::string method, const double redshift, const bool store_output, const std::string output_root, const int norm, const double k_min=0.001, const double k_max=100., const double prec=1.e-2);

      /**
       *  @brief the expanded A and B correction terms for the TNS
       *  model
       *
       *  the expanded A and B correction terms are computed at
       *  1-loop using (Standard) Perturbation Theory implemented in
       *  the CPT Library
       *  [http://www2.yukawa.kyoto-u.ac.jp/~atsushi.taruya/cpt_pack.html].
       *  Details can be found in (Taruya et.al, 2010)
       *  [https://arxiv.org/abs/1006.0699].
       *
       *  \f[ A(k,\mu,\beta) = \beta\mu^2A_{11} + \beta^2(\mu^2A_{12}
       *  + \mu^4A_{22}) + \beta^3(\mu^4A_{23} + \mu^6A_{33}) \f]
       *
       *  \f[ B(k,\mu,\beta) = \mu^2(\beta^2B_{12} + \beta^3B_{13} +
       *  \beta^4B_{14}) + \mu^4(\beta^2B_{22} + \beta^3B_{23} +
       *  \beta^4B_{24}) + \mu^6(\beta^3B_{33} + \beta^4B_{34}) +
       *  \mu^8\beta^4B_{44} \f]
       *
       *  @author J.E. Garcia-Farieta
       *  @author joegarciafa@unal.edu.co
       *
       *  @param kk the wave vector module
       *
       *  @param mu the cosine of the angle between the separation and the line of sight
       *
       *  @param linear_growth_rate the linear growth rate
       *
       *  @param bias the bias
       *
       *  @param method method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  CLASS [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created
       *  by CAMB are stored; if false the output files created by
       *  CAMB are removed
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If this parameter is different from the default value
       *  it will be used also in the output directory name
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
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param prec accuracy of the integration
       *
       *  @return A and B terms (total contribution) for TNS model
       */ 
      std::vector<std::vector<double>> Pk_TNS_AB_1loop (std::vector<double> kk, const double mu, const double linear_growth_rate, const double bias, const std::string method, const double redshift, const bool store_output, const std::string output_root, const int norm, const double k_min=0.001, const double k_max=100., const double prec=1.e-2);

      /**
       *  @brief the non-linear \f$\delta-\delta\f$,
       *  \f$\delta-\theta\f$, \f$\theta-\theta\f$ matter power
       *  spectra
       *
       *  the non-linear \f$\delta-\delta\f$, \f$\delta-\theta\f$,
       *  \f$\theta-\theta\f$ matter power spectra are computed at
       *  1-loop using (Standard) Perturbation Theory as implemented
       *  in the CPT Library
       *  [http://www2.yukawa.kyoto-u.ac.jp/~atsushi.taruya/cpt_pack.html].
       *  Details can be found in Taruya et al. (2010)
       *  [https://arxiv.org/abs/1006.0699]. The density,
       *  density-velocity divergence and velocity divergence-velocity
       *  divergence power spectra are defined as:
       *
       *  \f[ < \delta(k)\delta(k')> = (2\pi)^3\delta(k +
       *  k')P_{\delta\delta}(k) \f]
       *
       *  \f[ < \delta(k)\theta(k')> = (2\pi)^3\delta(k +
       *  k')P_{\delta\theta}(k) \f]
       *
       *  \f[ < \theta(k)\theta(k')> = (2\pi)^3\delta(k +
       *  k')P_{\theta\theta}(k) \f]
       *
       *  @author J.E. Garcia-Farieta
       *  @author joegarciafa@unal.edu.co
       *
       *  @param kk the wave vector module
       *
       *  @param method method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  CLASS [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created
       *  by CAMB are stored; if false the output files created by
       *  CAMB are removed
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can
       *  be any name. If this parameter is different from the default value
       *  it will be used also in the output directory name
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
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @return a vector of vectors containing
       *  \f$\left[P_{\delta\delta}(k), P_{\delta\theta}(k),
       *  P_{\theta\theta}(k)\right]\f$
       */
      std::vector<std::vector<double>> Pk_TNS_dd_dt_tt (std::vector<double> kk, const std::string method, const double redshift, const bool store_output, const std::string output_root, const int norm, const double k_min=0.001, const double k_max=100., const double prec=1.e-2);

      ///@}

    };
  }
}

#endif
    
