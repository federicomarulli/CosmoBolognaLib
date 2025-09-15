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
 *  @file Headers/PkXizSpace.h
 *
 *  @brief The class PkXizSpace
 *
 *  This file defines the interface of the class PkXizSpace, used to
 *  model the two-point correlation function and power spectrum in
 *  redshift space
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __PkXizSpace__
#define __PkXizSpace__

#include "PkXi.h"


// ===================================================================================================


namespace cbl {
  
  namespace cosmology {

    /**
     *  @class PkXizSpace PkXizSpace.h "Headers/PkXizSpace.h"
     *
     *  @brief The class PkXizSpace
     *
     *  This class is used to handle objects of type <EM>
     *  PkXizSpace </EM>.
     */
    class PkXizSpace {

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
      PkXizSpace () = default;
      
      /**
       *  @brief constructor
       *
       *  @param cosmology object of class Cosmology
       */
      PkXizSpace (std::shared_ptr<cbl::cosmology::Cosmology> cosmology)
      : m_cosmology(std::move(cosmology)) {}
      
      /**
       *  @brief default destructor
       */
      ~PkXizSpace () = default;

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
       *  @name Functions used to model two-point statistics in
       *  redshift-space
       */
      ///@{

      /**
       *  @brief monopole of the redshift-space two-point correlation
       *  function in the Kaiser limit
       *
       *  this function provides the monopole of the two-point
       *  correlation function of biased tracers, predicted at large
       *  scales in the Kaiser limit:
       *
       *  \f[ \xi_0(s) =
       *  \xi_{DM}(s)\left(\frac{b\sigma_8}{\sigma_8}\right)^2 \left[
       *  1 + \frac{2}{3}\frac{f\sigma_8}{b\sigma_8} +
       *  \frac{1}{5}\left(\frac{f\sigma_8}{b\sigma_8}\right)^2
       *  \right] \f]
       *
       *  The Kaiser factor \f$ \left[ 1 +
       *  \frac{2}{3}\frac{f\sigma_8}{b\sigma_8} +
       *  \frac{1}{5}\left(\frac{f\sigma_8}{b\sigma_8}\right)^2
       *  \right] \f$ is computed by cbl::xi_ratio
       *
       *  @param rad comoving separation
       *
       *  @param f_sigma8 f*&sigma;<SUB>8</SUB>
       *
       *  @param bias_sigma8 b*&sigma;<SUB>8</SUB>
       *
       *  @param method_Pk method used to compute the power spectrum
       *  and &sigma;(mass) (i.e. the Boltzmann solver); valid choices
       *  for method_Pk are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param xiType 0 \f$\rightarrow\f$ standard; 1
       *  \f$\rightarrow\f$ Chuang & Wang model
       *
       *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
       *
       *  @param NL false \f$\rightarrow\f$ linear power spectrum;
       *  true \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param r_min minimum separation up to which the
       *  correlation function is computed
       *
       *  @param r_max maximum separation up to which the
       *  correlation function is computed
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param aa parameter \e a of Eq. 24 of Anderson et
       *  al. 2012
       *
       *  @param GSL true \f$\rightarrow\f$ the GSL libraries are
       *  used
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return &xi;<SUB>0</SUB>
       *
       *  @warning working for \f$z>0\f$
       *
       */
      double xi0_Kaiser (const double rad, const double f_sigma8, const double bias_sigma8, const std::string method_Pk, const double redshift, const bool store_output=true, const std::string output_root="test", const bool xiType=0, const double k_star=-1., const bool NL=false, const int norm=-1, const double r_min=0.1, const double r_max=150., const double k_min=0.001, const double k_max=100., const double aa=0., const bool GSL=false, const double prec=1.e-2, const std::string file_par=par::defaultString);

       /**
       *  @brief monopole of the redshift-space two-point correlation
       *  function in the Kaiser limit
       *
       *  this function provides the monopole of the two-point
       *  correlation function of biased tracers, predicted at large
       *  scales in the Kaiser limit:
       *
       *  \f[ \xi_0(s) = \xi_{\rm DM}(s)b^2\left[1 + \frac{2\beta}{3} +
       *  \frac{\beta^2}{5}\right] \f]
       * 
       *  where \f$ \beta=\beta(z)=\frac{f(z)}{b(z)} \f$. The Kaiser
       *  factor \f$ \left[1 + \frac{2\beta}{3} +
       *  \frac{\beta^2}{5}\right] \f$ is computed by cbl::xi_ratio,
       *  while the linear growth rate \f$ f(z) \f$ is computed by
       *  cbl::cosmology::Cosmology::linear_growth_rate
       *
       *  @param rad comoving separations
       *
       *  @param bias the bias
       *
       *  @param method_Pk method used to compute the power spectrum
       *  and &sigma;(mass) (i.e. the Boltzmann solver); valid choices
       *  for method_Pk are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param NL false \f$\rightarrow\f$ linear power spectrum;
       *  true \f$\rightarrow\f$ non-linear power spectrum
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
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return &xi;<SUB>0</SUB>
       *
       */
      std::vector<double> xi0_Kaiser (const std::vector<double> rad, const double bias, const std::string method_Pk, const double redshift, const bool store_output=true, const std::string output_root="test", const bool NL=false, const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString);

      /**
       *  @brief 2D correlation function, &xi;(r<SUB>p</SUB>,&pi;),
       *  predicted by the dispersion model
       *
       *  @param rp r<SUB>p</SUB>: the comoving separation perpendicular
       *  to the line-of-sight
       *
       *  @param pi &pi;: the comoving separation parallel to the
       *  line-of-sight
       *
       *  @param f_sigma8 f*&sigma;<SUB>8</SUB>
       *
       *  @param bias_sigma8 b*&sigma;<SUB>8</SUB>
       *
       *  @param sigmav &sigma;<SUB>12</SUB>: pairwise peculiar
       *  velocity dispersion
       *
       *  @param method_Pk method used to compute the power spectrum
       *  and &sigma;(mass) (i.e. the Boltzmann solver); valid choices
       *  for method_Pk are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param FV 0 \f$\rightarrow\f$ exponential form for f(v); 1
       *  \f$\rightarrow\f$ Gaussian form for f(v); where f(v) is the
       *  velocity distribution function
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1
       *  \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param rr vector of r, the module of the comoving
       *  separation
       *
       *  @param Xi vector of &xi;(r), the two-point correlation
       *  function of dark matter
       *
       *  @param Xi_ vector of barred &xi;(r),
       *
       *  @param Xi__ vector of double-barred &xi;(r)
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param index internal parameter used when minimizing the
       *  &chi;<SUB>2</SUB>
       *
       *  @param bias_nl 0 \f$\rightarrow\f$ linear bias; 1
       *  \f$\rightarrow\f$ non-linear bias
       *
       *  @param bA b<SUB>a</SUB> non-linear bias parameter
       *
       *  @param xiType 0 \f$\rightarrow\f$ standard; 1
       *  \f$\rightarrow\f$ Chuang & Wang model
       *
       *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
       *
       *  @param xiNL 0 \f$\rightarrow\f$ linear power spectrum; 1
       *  \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param v_min minimum velocity used in the convolution of the
       *  correlation function
       *
       *  @param v_max maximum velocity used in the convolution of the
       *  correlation function
       *
       *  @param step_v number of steps used in the convolution of the
       *  correlation function
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param r_min minimum separation up to which the
       *  correlation function is computed
       *
       *  @param r_max maximum separation up to which the
       *  correlation function is computed
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param aa parameter \e a of Eq. 24 of Anderson et
       *  al. 2012
       *
       *  @param GSL true \f$\rightarrow\f$ the GSL libraries are
       *  used
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return &xi;(r<SUB>p</SUB>,&pi;)
       */
      double xi2D_dispersionModel (const double rp, const double pi, const double f_sigma8, const double bias_sigma8, const double sigmav, const std::string method_Pk, const double redshift, const int FV, const bool NL, std::vector<double> rr, std::vector<double> &Xi, std::vector<double> &Xi_, std::vector<double> &Xi__, const bool store_output=true, const std::string output_root="test", const int index=-1, const bool bias_nl=0, const double bA=-1., const bool xiType=0, const double k_star=-1., const bool xiNL=0, const double v_min=-3000., const double v_max=3000., const int step_v=500, const int norm=-1, const double r_min=0.1, const double r_max=150., const double k_min=0.001, const double k_max=100., const double aa=0., const bool GSL=false, const double prec=1.e-2, const std::string file_par=par::defaultString);
  
      /**
       *  @brief the function &xi;<SUB>g,nw</SUB>(s) of the Chuang &
       *  Wang 2012 model
       *
       *  see Chuang & Wang 2012, 1209.0210
       *
       *  @param rp r<SUB>p</SUB>: the comoving separation perpendicular
       *  to the line-of-sight
       *
       *  @param pi &pi;: the comoving separation parallel to the
       *  line-of-sight
       *      
       *  @param f_sigma8 f*&sigma;<SUB>8</SUB>
       *
       *  @param bias_sigma8 b*&sigma;<SUB>8</SUB>
       *
       *  @param bA b<SUB>a</SUB> non-linear bias parameter
       *
       *  @param redshift the redshift
       *    
       *  @param rr vector of r, the module of the comoving
       *  separation
       *
       *  @param Xi vector of &xi;(r), the two-point correlation
       *  function of dark matter
       *
       *  @param Xi_ vector of barred &xi;(r),
       *
       *  @param Xi__ vector of double-barred &xi;(r)
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @return &xi;<SUB>g,nw</SUB>(s)
       */
      double xisnl_gnw (const double rp, const double pi, const double f_sigma8, const double bias_sigma8, const double bA, const double redshift, std::vector<double> rr, std::vector<double> Xi, std::vector<double> &Xi_, std::vector<double> &Xi__, const bool store_output=true, const std::string output_root="test");
 
      /**
       *  @brief the function &xi;<SUB>g,BAO</SUB>(s) of the Chuang &
       *  Wang 2012 model
       *
       *  see Chuang & Wang 2012, 1209.0210
       *
       *  @param rp r<SUB>p</SUB>: the comoving separation perpendicular
       *  to the line-of-sight
       *
       *  @param pi &pi;: the comoving separation parallel to the
       *  line-of-sight
       *      
       *  @param f_sigma8 f*&sigma;<SUB>8</SUB>
       *
       *  @param bias_sigma8 b*&sigma;<SUB>8</SUB>
       *
       *  @param redshift the redshift
       *
       *  @param rr vector of r, the module of the comoving
       *  separation
       *
       *  @param Xi vector of &xi;(r), the two-point correlation
       *  function of dark matter
       *
       *  @param Xi_ vector of barred &xi;(r),
       *
       *  @param Xi__ vector of double-barred &xi;(r)
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
       *
       *  @param x_min minimum velocity used in the integral of the
       *  Chuang & Wang model
       *
       *  @param x_max maximum velocity used in the integral of the
       *  Chuang & Wang model
       *
       *  @param step_x number of steps in the integral of the Chuang
       *  & Wang model
       *
       *  @return &xi;<SUB>g,BAO</SUB>(s)
       */
      double xis_gBAO (const double rp, const double pi, const double f_sigma8, const double bias_sigma8, const double redshift, std::vector<double> rr, std::vector<double> Xi, std::vector<double> &Xi_, std::vector<double> &Xi__, const bool store_output=true, const std::string output_root="test", const double k_star=-1., const double x_min=-3000., const double x_max=3000., const int step_x=500);
 
      /**
       *  @brief 2D correlation function, &xi;(r<SUB>p</SUB>,&pi;),
       *  predicted by the Chuang & Wang model
       *
       *  @param rp r<SUB>p</SUB>: the comoving separation
       *  perpendicular to the line-of-sight
       *
       *  @param pi &pi;: the comoving separation parallel to the
       *  line-of-sight
       *
       *  @param beta &beta;=f/b, where f is the linear growth rate
       *  and b is the bias
       *
       *  @param bias_lin linear bias
       *
       *  @param bA b<SUB>a</SUB> non-linear bias parameter
       *
       *  @param sigmav0 &sigma;<SUB>0</SUB>(v): parameter of the
       *  velocity distribution function, f(v)
       *
       *  @param cmu parameter of the velocity distribution function,
       *  f(v)
       *
       *  @param cs1 parameter of the velocity distribution function,
       *  f(v)
       *
       *  @param cs2 parameter of the velocity distribution function,
       *  f(v)
       *
       *  @param redshift the redshift
       *
       *  @param rr1 vector of r, the module of the comoving separation
       *
       *  @param Xi1 vector of &xi;(r), the two-point correlation
       *  function of dark matter
       *
       *  @param rr2 vector of r, the module of the comoving separation
       *
       *  @param Xi2 vector of &xi;(r), the two-point correlation
       *  function of dark matter
       *
       *  @param Xi1_ vector of barred &xi;(r),
       *
       *  @param Xi1__ vector of double-barred &xi;(r)
       *   
       *  @param Xi2_ vector of barred &xi;(r),
       *
       *  @param Xi2__ vector of double-barred &xi;(r)
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param BAO 0 \f$\rightarrow\f$ no BAO convolution; 1
       *  \f$\rightarrow\f$ BAO convolution
       *
       *  @param xiType 0 \f$\rightarrow\f$ standard; 1
       *  \f$\rightarrow\f$ Chuang & Wang model
       *
       *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
       *
       *  @param xiNL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
       *  non-linear power spectrum
       *
       *  @param r_min minimum separation up to which the
       *  correlation function is computed
       *
       *  @param r_max maximum separation up to which the
       *  correlation function is computed
       *
       *  @param v_min minimum velocity used in the convolution of the
       *  correlation function
       *
       *  @param v_max maximum velocity used in the convolution of the
       *  correlation function
       *
       *  @param step_v number of steps used in the convolution of the
       *  correlation function
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param x_min minimum velocity used in the integral of the
       *  Chuang & Wang model
       *
       *  @param x_max maximum velocity used in the integral of the
       *  Chuang & Wang model
       *
       *  @param step_x number of steps in the integral of the Chuang
       *  & Wang model
       *
       *  @param aa parameter \e a of Eq. 24 of Anderson et
       *  al. 2012
       *
       *  @param GSL true \f$\rightarrow\f$ the GSL libraries are
       *  used
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return &xi;(r<SUB>p</SUB>,&pi;)
       */
      double xi2D_CW (const double rp, const double pi, const double beta, const double bias_lin, const double bA, const double sigmav0, const double cmu, const double cs1, const double cs2, const double redshift, std::vector<double> rr1, std::vector<double> Xi1, std::vector<double> rr2, std::vector<double> Xi2, std::vector<double> &Xi1_, std::vector<double> &Xi1__, std::vector<double> &Xi2_, std::vector<double> &Xi2__, const bool store_output=true, const std::string output_root="test", const bool BAO=1, const bool xiType=0, const double k_star=-1, const bool xiNL=0, const double r_min=0.1, const double r_max=150., const double v_min=-3000., const double v_max=3000., const int step_v=500, const double k_min=0.001, const double k_max=100., const double x_min=-3000., const double x_max=3000., const int step_x=500, const double aa=0., const bool GSL=false, const double prec=1.e-2, const std::string file_par=par::defaultString);
      
      ///@}

    };
  }
}

#endif
    
