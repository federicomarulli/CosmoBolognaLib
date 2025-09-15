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
 *  @file Headers/Sigma.h
 *
 *  @brief The class Sigma
 *
 *  This file defines the interface of the class Sigma, used to model
 *  the power spectrum and two-point correlation function of
 *  extra-galactic sources
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __SIGMA__
#define __SIGMA__

#include "PkXi.h"


// ===================================================================================================


namespace cbl {
  
  namespace cosmology {

    /**
     *  @class Sigma Sigma.h "Headers/Sigma.h"
     *
     *  @brief The class Sigma
     *
     *  This class is used to handle objects of type <EM> Sigma
     *  </EM>. It is used to model the variance of the cosmic density
     *  field
     */
    class Sigma {

    private:

      /// pointer to the input cosmology
      std::shared_ptr<cosmology::Cosmology> m_cosmology = NULL;
      
      /**
       *  @brief Private member function to compute a vector of 
       *  mass variances and their derivatives
       *
       *  This function computes the mass variances
       *  and their derivatives as follows:
       *
       *  \f[ \sigma^2(R) = \frac{1}{2\pi^2}\int_0^\infty {\rm d}k\,
       *  k^2 P_{lin}(k, z) F^2(k, R)\f]
       *
       *  where \f$F(x)\f$ is a generic filter
       *
       *  @param pkInterp the interpolated matter power spectrum
       * 
       *  @param filter the filter
       * 
       *  @param radius the vector of radii used for filtering
       * 
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       * 
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return the funciton to compute the sigma mass variances
       */
      std::vector<double> m_func_sigma (const std::function<double(double)> pkInterp, std::function<double(double,int)> filter, const std::vector<double> radius, const double k_max=100., const bool unit1=false) const;

      /**
       *  @brief The mass variance, \f$\sigma^2(R)\f$
       *
       *  This private member function computes the variance of the
       *  linear density field as follows:
       *
       *  \f[ \sigma^2(R) = \frac{1}{2\pi^2}\int_0^\infty {\rm d}k\,
       *  k^2 P_{lin}(k, z) W^2(k, R)\f]
       *
       *  where \f$W(x)=(3/x)^3(\sin x-x\cos x)\f$ and
       *  \f$R=(3M/4\pi\rho_m)^{1/3}\f$
       *
       *  @param radius the vector of radii, \f$R\f$
       * 
       *  @param pkInterp the interpolated matter power spectrum
       * 
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *  
       *  @param fact the normalization factor 
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return \f$\sigma^2(R)\f$
       */
      std::vector<double> m_sigma2R (const std::vector<double> radius, const std::function<double(double)> pkInterp, const double k_max=100., const double fact=1., const bool unit1=false) const;

      /**
       *  @brief The mass variance, \f$\sigma^2(M)\f$
       *
       *  This private member function computes the variance of the
       *  linear density field as follows:
       *
       *  \f[ \sigma^2(M) = \frac{1}{2\pi^2}\int_0^\infty {\rm d}k\,
       *  k^2 P_{lin}(k, z) W^2(k, R)\f]
       *
       *  where \f$W(x)=(3/x)^3(\sin x-x\cos x)\f$ and
       *  \f$R=(3M/4\pi\rho_m)^{1/3}\f$
       *
       *  @param mass the vector of masses
       * 
       *  @param pkInterp the interpolated matter power spectrum
       * 
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       * 
       *  @param fact the normalization factor
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return the not-yet-normalised \f$\sigma^2(M)\f$
       */
      std::vector<double> m_sigma2M (const std::vector<double> mass, const std::function<double(double)> pkInterp, const double k_max=100., const double fact=1., const bool unit1=false) const;

      
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Sigma () = default;
      
      /**
       *  @brief constructor
       *
       *  @param cosmology pointer to an object of class Cosmology
       */
      Sigma (std::shared_ptr<cbl::cosmology::Cosmology> cosmology)
	: m_cosmology(std::move(cosmology)) {}
      
      /**
       *  @brief default destructor
       */
      ~Sigma () = default;

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
       *  @name Functions to estimate the variance of the cosmic
       *  density field
       */
      ///@{ 

      /**
       *  @brief The mass variance, \f$\sigma^2(R)\f$
       *
       *  This function computes the variance of the linear density
       *  field:
       *
       *  \f[ \sigma^2(R)=\frac{1}{2\pi^2}\int_0^\infty {\rm d}k\, k^2
       *  P_{lin}(k, z) W^2(k, R)\f]
       *
       *  where \f$W(x)=(3/x)^3(\sin x-x\cos x)\f$
       *
       *  @param radius the vector of radii, \f$R\f$
       *
       *  @param pk0Interp the interpolated matter power spectrum at
       *  \f$z=0\f$
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return \f$\sigma^2(R)\f$
       */
      std::vector<double> sigma2R (const std::vector<double> radius, const std::function<double(double)> pk0Interp, const double k_max=100., const bool unit1=false) const;
      
      /**
       *  @brief The mass variance, \f$\sigma^2(R)\f$
       *
       *  This function computes the variance of the linear density
       *  field:
       *
       *  \f[ \sigma^2(R)=\frac{1}{2\pi^2}\int_0^\infty {\rm d}k\, k^2
       *  P_{lin}(k, z) W^2(k, R)\f]
       *
       *  where \f$W(x)=(3/x)^3(\sin x-x\cos x)\f$
       *
       *  @param radius the vector of radii, \f$R\f$
       *
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
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
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return \f$\sigma^2(R)\f$
       */
      std::vector<double> sigma2R (const std::vector<double> radius, const std::string method_Pk, const double redshift, const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true, const bool unit1=false) const;

      /**
       *  @brief The nth-order derivative of the mass variance,
       *  \f${\rm d}^n\sigma^2(R)/{\rm d}R^n\f$
       *
       *  @param nd the derivative order, \f$n\f$
       *
       *  @param radius the vector of radii, \f$R\f$
       *
       *  @param pk0Interp the interpolated matter power spectrum at
       *  \f$z=0\f$
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return the nth-order derivative of the mass variance
       *
       *  @warning the current implementation computes the derivative
       *  using the simplest numerical approximation, with fixed
       *  incremental step; it is computationally efficient, but the
       *  accuracy might be lowt
       */
      std::vector<double> dnsigma2R (const int nd, const std::vector<double> radius, std::function<double(double)> pk0Interp, const double k_max=100., const bool unit1=false) const;

      /**
       *  @brief The nth-order derivative of the mass variance,
       *  \f${\rm d}^n\sigma^2(R)/{\rm d}R^n\f$
       *
       *  @param nd the derivative order, \f$n\f$
       *
       *  @param radius the vector of radii, \f$R\f$
       *
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
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
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return the nth-order derivative of the mass variance
       *
       *  @warning the current implementation computes the derivative
       *  using the simplest numerical approximation, with fixed
       *  incremental step; it is computationally efficient, but the
       *  accuracy might be lowt
       */
      std::vector<double> dnsigma2R (const int nd, const std::vector<double> radius, const std::string method_Pk, const double redshift, const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true, const bool unit1=false) const;

      /**
       *  @brief The mass variance, \f$\sigma^2(M)\f$
       *
       *  This function computes the variance of the linear density
       *  field:
       *
       *  \f[ \sigma^2(M) = \frac{1}{2\pi^2}\int_0^\infty {\rm d}k\,
       *  k^2 P_{lin}(k, z) W^2(k, R)\f]
       *
       *  where \f$W(x)=(3/x)^3(\sin x-x\cos x)\f$ and
       *  \f$R=(3M/4\pi\rho_m)^{1/3}\f$
       *
       *  @param mass the vector of masses 
       *
       *  @param pk0Interp the interpolated matter power spectrum at
       *  \f$z=0\f$
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return \f$\sigma^2(M)\f$
       */
      std::vector<double> sigma2M (const std::vector<double> mass, const std::function<double(double)> pk0Interp, const double k_max=100., const bool unit1=false) const;

      /**
       *  @brief The mass variance, \f$\sigma^2(M)\f$
       *
       *  This function computes the variance of the linear density
       *  field:
       *
       *  \f[ \sigma^2(M) = \frac{1}{2\pi^2}\int_0^\infty {\rm d}k\,
       *  k^2 P_{lin}(k, z) W^2(k, R)\f]
       *
       *  where \f$W(x)=(3/x)^3(\sin x-x\cos x)\f$ and
       *  \f$R=(3M/4\pi\rho_m)^{1/3}\f$
       *
       *  @param mass the vector of masses 
       *
       *  @param method_Pk the method used to compute the power
       *  spectrum; valid choices for method_Pk are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
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
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return \f$\sigma^2(M)\f$
       */
      std::vector<double> sigma2M (const std::vector<double> mass, const std::string method_Pk, const double redshift, const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true, const bool unit1=false) const;
      
      /**
       *  @brief The first derivative of the mass variance, \f${\rm
       *  d}^n\sigma^2(M)/{\rm d}M^n\f$
       * 
       *  @param nd the derivative order, \f$n\f$
       *
       *  @param mass the vector of masses, \f$M\f$
       *
       *  @param pk0Interp the interpolated matter power spectrum at \f$z=0\f$
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return the first derivative of the mass variance
       *
       *  @warning the current implementation computes the derivative
       *  using the simplest numerical approximation, with fixed
       *  incremental step; it is computationally efficient, but the
       *  accuracy might be low
       */
      std::vector<double> dnsigma2M (const int nd, const std::vector<double> mass, const std::function<double(double)> pk0Interp, const double k_max=100., const bool unit1=false) const;

      /**
       *  @brief The first derivative of the mass variance, \f${\rm
       *  d}^n\sigma^2(M)/{\rm d}M^n\f$
       * 
       *  @param nd the derivative order, \f$n\f$
       *
       *  @param mass the vector of masses, \f$M\f$
       *
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
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
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return the first derivative of the mass variance
       *
       *  @warning the current implementation computes the derivative
       *  using the simplest numerical approximation, with fixed
       *  incremental step; it is computationally efficient, but the
       *  accuracy might be low
       */
      std::vector<double> dnsigma2M (const int nd, const std::vector<double> mass, const std::string method_Pk, const double redshift, const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true, const bool unit1=false) const;

      /**
       *  @brief Auxiliary function to create a grid file with
       *  &sigma;(M)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
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
       *  @param interpType method to interpolate the power spectrum
       *     
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return file_grid name of the file where the grid is stored
       */
      std::string create_grid_sigmaM (const std::string method_SS, const double redshift, const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const; 

      /**
       *  @brief The hierarchical moments S<SUB>n</SUB>
       *
       *  this function provides the hierarchical moments S<SUB>n</SUB>
       *  given by the perturbation theory (see e.g. Juszkiewicz et
       *  al. 1993, Bernardeau 1994, Wolk 2013)
       *
       *  @param nn order of the moment
       *
       *  @param RR comoving separation
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
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
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return the hierarchical moments, S<SUB>n</SUB>, given by the
       *  perturbation theory
       */
      double Sn_PT (const int nn, const double RR, const std::string method_SS, const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;
  
      /**
       *  @brief The deprojected hierarchical moments
       *  &Sigma;<SUB>n</SUB>
       *
       *  this function provides the deprojected hierarchical moments
       *  &Sigma;<SUB>n</SUB> given by the perturbation theory (see
       *  e.g. Juszkiewicz et al. 1993, Bernardeau 1994, Wolk 2013)
       *
       *  @param nn order of the moment
       *
       *  @param RR comoving separation
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
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
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return the deprojected hierarchical moments,
       *  &Sigma;<SUB>n</SUB>, given by the perturbation theory
       */
      double Sigman_PT (const int nn, const double RR, const std::string method_SS, const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;

      ///@}
      
      
    };
  }
  
}

#endif
