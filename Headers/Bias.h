/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/Bias.h
 *
 *  @brief The class Bias
 *
 *  This file defines the interface of the class Bias, used to model
 *  bias of dark matter haloes
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __BIAS__
#define __BIAS__

#include "MassFunction.h"


// ===================================================================================================


namespace cbl {
  
  namespace cosmology {
    
    /**
     *  @class Bias Bias.h "Headers/Bias.h"
     *
     *  @brief The class Bias
     *
     *  This class is used to handle objects of type <EM> Bias
     *  </EM>. It is used to model the bias of dark matter haloes
     */
    class Bias {

    private:
      
      /// pointer to the input cosmology
      std::shared_ptr<cosmology::Cosmology> m_cosmology = NULL;
      
      /// object of class MassFunction
      cbl::cosmology::MassFunction m_mass_function;
      
      /**
       *  @brief auxiliary function to compute the halo bias
       *
       *  @param Sigma &sigma;(mass, z=0): the mass variance at z=0
       *
       *  @param redshift the redshift
       *
       *  @param author author(s) who proposed the bias; valid authors
       *  are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen
       *  2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction
       *  of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param Delta \f$\Delta\f$, the overdensity
       *
       *  @param dlnsigma_dlnR \f$\frac{\mathrm{d log} \sigma}{\mathrm{d log} R}\f$ the derivative of the variance
       *
       *  @return the halo bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ to obtain the background overdensity
       */
      double m_bias_halo_generator (const double Sigma, const double redshift, const std::string author, const double Delta=200., const double dlnsigma_dlnR=0.) const;    
          
      /**
       *  @brief auxiliary function to compute the halo bias
       *
       *  @param Sigma &sigma;(mass, z=0): the mass variance at z=0
       *
       *  @param redshift the redshift
       *
       *  @param D_N the amplitude of the growing mode
       *
       *  @param author author(s) who proposed the bias; valid authors
       *  are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen
       *  2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction
       *  of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param Delta \f$\Delta\f$, the overdensity
       *
       *  @param dlnsigma_dlnR \f$\frac{\mathrm{d log} \sigma}{\mathrm{d log} R}\f$ the derivative of the variance
       *
       *  @return the halo bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ to obtain the background overdensity
       */
      double m_bias_halo_generator (const double Sigma, const double redshift, const double D_N, const std::string author, const double Delta=200., const double dlnsigma_dlnR=0.) const;  

    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Bias () = default;
      
      /**
       *  @brief constructor
       *
       *  @param cosmology pointer to an object of class Cosmology
       */
      Bias (std::shared_ptr<cbl::cosmology::Cosmology> cosmology)
      : m_cosmology(std::move(cosmology)) { cbl::cosmology::MassFunction MF(m_cosmology); m_mass_function = MF; }
      
      /**
       *  @brief default destructor
       */
      ~Bias () = default;
      
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
       *  @name Functions to estimate the bias of the dark matter
       *  haloes
       */
      ///@{
      
      /**
       *  @brief bias of dark matter haloes
       *
       *  @param Mass halo mass
       *
       *  @param redshift the redshift
       *
       *  @param author author(s) who proposed the bias; valid authors
       *  are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen
       *  2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction
       *  of Warren 2004), Tinker (Tinker et al. 2010)
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
       *  @param Delta \f$\Delta\f$, the overdensity
       *  
       *  @param kk wave vector module
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
       *  @return b<SUB>halo</SUB>: the dark matter bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ to obtain the background overdensity
       */
      double bias_halo (const double Mass, const double redshift, const std::string author, const std::string method_SS, const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double Delta=200., const double kk=-1., const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true); 
          
      /**
       *  @brief bias of dark matter haloes
       *
       *  @param Mass halo mass
       *
       *  @param Sigma &sigma;(mass, z=0): the mass variance at z=0
       *
       *  @param redshift the redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010),
       *  Castro (Castro et al. 2023)
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
       *  @param Delta \f$\Delta\f$, the overdensity
       *  
       *  @param kk wave vector module
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
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
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
       *  @return b<SUB>halo</SUB>: the dark matter bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ to obtain the background overdensity
       */
      double bias_halo (const double Mass, const double Sigma, const double redshift, const std::string model_bias, const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double Delta=200., const double kk=-1., const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string method_SS="CAMB", const std::string input_file=par::defaultString, const bool is_parameter_file=true);
                
      /**
       *  @brief bias of dark matter haloes
       *
       *  @param Mass halo mass
       *
       *  @param Sigma &sigma;(mass, z=0): the mass variance at z=0
       *
       *  @param redshift the redshift
       *
       *  @param DN normalised amplitude of the growing mode at a
       *  given redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010),
       *  Castro (Castro et al. 2023)
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
       *  @param Delta \f$\Delta\f$, the overdensity
       *  
       *  @param kk wave vector module
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
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
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
       *  @return b<SUB>halo</SUB>: the dark matter bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ to obtain the background overdensity
       */
      double bias_halo (const double Mass, const double Sigma, const double redshift, const double DN, const std::string model_bias, const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double Delta=200., const double kk=-1., const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string method_SS="CAMB", const std::string input_file=par::defaultString, const bool is_parameter_file=true);
          
      /**
       *  @brief the effective bias of dark matter haloes, with masses
       *  in a given range and at a given mean redshift
       *
       *  this function computes the effective bias of dark matter
       *  haloes:
       *
       *  \f[ b_{eff}(z) = \frac{\int_{M_{min}}^{M_{max}} {\rm d}M\,
       *  b(M, z) \Phi(M, z)}{\int_{M_{min}}^{M_{max}} {\rm
       *  d}M\,\Phi(M, z)} \f]
       *
       *  in the current implementation, the integral is actually
       *  replaced by the Riemann sum, as follows:
       *
       *  \f[ b_{eff}(z) \simeq \frac{\sum_{i_{min}}^{i_{max}} b(M, z)
       *  \Phi(M, z) (M_{i+1}-M_i)}{\sum_{i_{min}}^{i_{max}}
       *  \Phi(M, z) (M_{i+1}-M_i)} \f]
       *
       *  where the halo mass function, \f$\Phi(M, z)\f$, is computed
       *  by cbl::cosmology::MassFunction::mass_function and the
       *  linear bias, \f$b(M, z)\f$, is computed by
       *  cbl::cosmology::Bias::bias_halo
       *
       *  @param Mass_min minimum halo mass
       *
       *  @param Mass_max maximum halo mass
       *
       *  @param redshift the redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010),
       *  Castro (Castro et al. 2023)
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press & Schechter), ST (Sheth &
       *  Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       *  al. 2006), ShenH (halo MF, Shen et al. 2006), ShenF
       *  (filament MF, Shen et al. 2006), ShenS (sheet MF, Shen et
       *  al. 2006), Pan (Pan 2007), Peacock (Peacock at al. 2007),
       *  Reed (Reed et al. 2007), Tinker (Tinker et al. 2008), Crocce
       *  (Crocce et al. 2010), Manera (Manera et al. 2010), Courtin
       *  (Courtin et al. 2010), Bhattacharya (Bhattacharya et
       *  al. 2011), Angulo_FOF (FoF MF, Angulo et al. 2012),
       *  Angulo_Sub (SUBFIND MF, Angulo et al. 2012), Watson_FOF (FoF
       *  MF, Watson et al. 2012), Watson_SOH (Spherical Overdensity
       *  halo MF, Watson et al. 2012), Despali_Z0, Despali_AllZ,
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016)
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  CLASS [http://class-code.net/], EisensteinHu
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
       *  @param Delta \f$\Delta\f$, the overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
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
       *  @return b<SUB>eff</SUB>: the effective dark matter bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ to obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double bias_eff (const double Mass_min, const double Mass_max, const double redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const bool store_output=true, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
          
      /**
       *  @brief effective bias of dark matter haloes, computed by
       *  averaging over the bias of a given set of haloes
       *
       *  this function computes the effective bias of dark matter
       *  haloes as follows:
       *
       *  \f[ b_{eff}(z) = \frac{\int_{M_{min}}^{M_{max}} {\rm d}M\,
       *  b(M, z) \Phi(M, z)}{\int_{M_{min}}^{M_{max}} {\rm
       *  d}M\,\Phi(M, z)} \f]
       *
       *  in the current implementation, the integral is actually
       *  replaced by the Riemann sum, as follows:
       *
       *  \f[ b_{eff}(z) \simeq \frac{\sum_{i_{min}}^{i_{max}} b(M, z)
       *  \Phi(M, z) (M_{i+1}-M_i)}{\sum_{i_{min}}^{i_{max}}
       *  \Phi(M, z) (M_{i+1}-M_i)} \f]
       *
       *  where the halo mass function is provided in input and the
       *  linear bias, \f$b(M, z)\f$, is computed by
       *  cbl::cosmology::Bias::bias_halo
       *
       *  @param MM vector of halo masses
       *
       *  @param MF vector of mass function values, d&Phi;/dM=dn(M)/dM
       *
       *  @param redshift the redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
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
       *  @param Delta \f$\Delta\f$, the overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
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
       *  @return b<SUB>eff</SUB>: the effective dark matter bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ to obtain the background overdensity
       */
      double bias_eff (const std::vector<double> MM, const std::vector<double> MF, const double redshift, const std::string model_bias, const std::string method_SS, const bool store_output=true, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
          
      /**
       *  @brief effective bias of dark matter haloes, computed by
       *  averaging the bias of a set of haloes, with the mass variance
       *  estimated from a grid
       *
       *  this function computes the effective bias of dark matter
       *  haloes by either averaging the bias of a set of haloes with
       *  a given mass:
       *
       *  \f[b_{eff}(z) = \frac{1}{N_{halo}}\sum_{i=1}^{N_{halo}}
       *  b(M_i, z_i) \; , \; (1)\f]
       *
       *  or by averaging over halo pairs:
       *
       *  \f[b_{eff}(z) = \sqrt{ \frac{2}{N_{halo}(N_{halo}-1)}
       *  \sum_{i=1}^{N_{halo}}\sum_{j=i+1}^{N_{halo}} b(M_i,
       *  z_i)b(M_j, z_j)} \; , \; (2)\f]
       *
       *  where the linear bias of the \f$i\f$-th halo, \f$b^{i}(M,
       *  z)\f$, is computed by cbl::cosmology::Bias::bias_halo; the
       *  mass variance is estimated from a grid by
       *  cbl::cosmology::Sigma::create_grid_sigmaM
       *
       *  @param MM vector containing the halo masses
       *
       *  @param redshift vector containing the redshifts; if it has
       *  size=1, it will be considered as the main redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param meanType meanType="mean_bias" \f$\rightarrow\f$ the
       *  effective bias is computed with Eq.(1);
       *  meanType="mean_pair_bias" \f$\rightarrow\f$ the effective
       *  bias is computed with Eq.(2)
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
       *  @param Delta_crit \f$\Delta_{crit}\f$: the critical
       *  overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
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
       *  @return a vector containing the mean and standard deviation
       *  of the effective dark matter bias
       */
      std::vector<double> bias_eff_mass_grid (const std::vector<double> MM, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType="mean_bias", const bool store_output=true, const std::string output_root="test", const double Delta_crit=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
          
      /**
       *  @brief effective bias of dark matter haloes, computed by
       *  averaging the bias of a set of haloes
       *
       *  this function computes the effective bias of dark matter
       *  haloes by either averaging the bias of a set of haloes with
       *  a given mass:
       *
       *  \f[b_{eff}(z) = \frac{1}{N_{halo}}\sum_{i=1}^{N_{halo}}
       *  b(M_i, z_i) \; , \; (1)\f]
       *
       *  or by averaging over halo pairs:
       *
       *  \f[b_{eff}(z) = \sqrt{ \frac{2}{N_{halo}(N_{halo}-1)}
       *  \sum_{i=1}^{N_{halo}}\sum_{j=i+1}^{N_{halo}} b(M_i,
       *  z_i)b(M_j, z_j)} \; , \; (2)\f]
       *
       *  where the linear bias of the \f$i\f$-th halo, \f$b^{i}(M,
       *  z)\f$, is computed by cbl::cosmology::Bias::bias_halo; the
       *  mass variance is computed by cbl::cosmology::Sigma::sigma2M
       *
       *  @param MM vector containing the halo masses
       *
       *  @param redshift vector containing the redshifts; if it has
       *  size=1, it will be considered as the main redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param meanType meanType="mean_bias" \f$\rightarrow\f$ the
       *  effective bias is computed with Eq.(1);
       *  meanType="mean_pair_bias" \f$\rightarrow\f$ the effective
       *  bias is computed with Eq.(2)
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
       *  @param Delta \f$\Delta\f$: the background overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
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
       *  @return a vector containing the mean and standard deviation
       *  of the effective dark matter bias
       */
      std::vector<double> bias_eff_mass (const std::vector<double> MM, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType="mean_bias", const bool store_output=true, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
                
      /**
       *  @brief effective bias of dark matter haloes, computed by
       *  averaging the bias of a set of haloes, interpolating the
       *  mass variance on a grid
       *
       *  this function computes the effective bias of dark matter
       *  haloes by either averaging the bias of a set of haloes with
       *  a given mass:
       *
       *  \f[b_{eff}(z) = \frac{1}{N_{halo}}\sum_{i=1}^{N_{halo}}
       *  b(M_i, z_i) \; , \; (1)\f]
       *
       *  or by averaging over halo pairs:
       *
       *  \f[b_{eff}(z) = \sqrt{ \frac{2}{N_{halo}(N_{halo}-1)}
       *  \sum_{i=1}^{N_{halo}}\sum_{j=i+1}^{N_{halo}} b(M_i,
       *  z_i)b(M_j, z_j)} \; , \; (2)\f]
       *
       *  where the linear bias of the \f$i\f$-th halo, \f$b^{i}(M,
       *  z)\f$, is computed by cbl::cosmology::Bias::bias_halo; the
       *  mass variance is computed by cbl::cosmology::Sigma::sigma2M
       *
       *  @param mass vector containing the halo masses
       *
       *  @param mass_grid vector containing the halo masses on the
       *  grid used to interpolate the mass variance
       *
       *  @param redshift vector containing the redshifts; if it has
       *  size=1, it will be considered as the main redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param meanType meanType="mean_bias" \f$\rightarrow\f$ the
       *  effective bias is computed with Eq.(1);
       *  meanType="mean_pair_bias" \f$\rightarrow\f$ the effective
       *  bias is computed with Eq.(2)
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
       *  @param Delta \f$\Delta\f$: the background overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
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
       *  @return a vector containing the mean and standard deviation
       *  of the effective dark matter bias
       */
      std::vector<double> bias_eff_mass (const std::vector<double> mass, const std::vector<double> mass_grid, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType="mean_bias", const bool store_output=true, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
                
      /**
       *  @brief effective bias of dark matter haloes, computed using
       *  a given selection function; &sigma;(mass) and dln&sigma;/dM
       *  are provided in input
       *
       *  this function computes the effective bias of dark matter
       *  haloes:
       *
       *  \f[ b_{eff}(z) = \frac{\int_{M_{min}}^{M_{max}} {\rm d}M\,
       *  b(M, z) \Phi(M, z) f(M, z)}{\int_{M_{min}}^{M_{max}} {\rm
       *  d}M\,\Phi(M, z) f(M, z)} \f]
       *
       *  in the current implementation, the integral is actually
       *  replaced by the Riemann sum, as follows:
       *
       *  \f[ b_{eff}(z) \simeq \frac{\sum_{i_{min}}^{i_{max}} b(M, z)
       *  \Phi(M, z) f(M, z) (M_{i+1}-M_i)}{\sum_{i_{min}}^{i_{max}}
       *  \Phi(M, z) f(M, z) (M_{i+1}-M_i)} \f]
       *
       *  where the halo mass function, \f$\Phi(M, z)\f$, is computed
       *  by cbl::cosmology::MassFunction::mass_function, the linear
       *  bias, \f$b(M, z)\f$, is computed by
       *  cbl::cosmology::Bias::bias_halo, and f(M, z) is the
       *  selection function
       *
       *  @param interp_sigma FuncGrid object containing the values
       *  of &sigma;(mass) computed on a grid, used by interpolation
       *
       *  @param interp_DnSigma FuncGrid object containing the values
       *  of dln&sigma;/dM computed on a grid, used by interpolation
       *
       *  @param interp_SF FuncGrid object containing the values of
       *  the selection function computed on a grid in mass, at the
       *  mean redshift, used by interpolation
       *
       *  @param Mass_min minimum halo mass
       *
       *  @param Mass_max maximum halo mass
       *
       *  @param redshift vector containing the input redshifts
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press & Schechter), ST (Sheth &
       *  Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       *  al. 2006), ShenH (halo MF, Shen et al. 2006), ShenF
       *  (filament MF, Shen et al. 2006), ShenS (sheet MF, Shen et
       *  al. 2006), Pan (Pan 2007), Peacock (Peacock at al. 2007),
       *  Reed (Reed et al. 2007), Tinker (Tinker et al. 2008), Crocce
       *  (Crocce et al. 2010), Manera (Manera et al. 2010), Courtin
       *  (Courtin et al. 2010), Bhattacharya (Bhattacharya et
       *  al. 2011), Angulo_FOF (FoF MF, Angulo et al. 2012),
       *  Angulo_Sub (SUBFIND MF, Angulo et al. 2012), Watson_FOF (FoF
       *  MF, Watson et al. 2012), Watson_SOH (Spherical Overdensity
       *  halo MF, Watson et al. 2012), Despali_Z0, Despali_AllZ,
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param alpha the \f$\alpha\f$ parameter of the cluster mass
       *  scaling relation
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
       *  @param Delta_crit \f$\Delta_{crit}\f$: the critical
       *  overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
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
       *  @param is_parameter_file true \f$\rightarrow\f$ the
       *  input_file is a parameter file, used to compute the power
       *  spectrum with the method specified by method_Pk; false
       *  \f$\rightarrow\f$ the input_file is a file containing the
       *  power spectrum
       *
       *  @return b<SUB>eff</SUB>: the effective dark matter bias
       *
       *  @warning interp_sigma and interp_DnSigma have to be
       *  computed at the same cosmology of the object. They are
       *  provided as an input just to improve the performances in
       *  some applications (e.g. MCMC) where these quantities can be
       *  computed once
       */
      std::vector<double> bias_eff_selection_function (const glob::FuncGrid interp_sigma, const glob::FuncGrid interp_DnSigma, const glob::FuncGrid interp_SF, const double Mass_min, const double Mass_max, const std::vector<double> redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const double alpha=1., const bool store_output=true, const std::string output_root="test", const double Delta_crit=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
                 
      /**
       *  @brief effective bias of dark matter haloes, computed using
       *  a given selection function; &sigma;(mass) and dln&sigma;/dM
       *  are provided in input
       *
       *  this function computes the effective bias of dark matter
       *  haloes:
       *
       *  \f[ b_{eff}(z) = \frac{\int_{M_{min}}^{M_{max}} {\rm d}M\,
       *  b(M, z) \Phi(M, z) f(M, z)}{\int_{M_{min}}^{M_{max}} {\rm
       *  d}M\,\Phi(M, z) f(M, z)} \f]
       *
       *  in the current implementation, the integral is actually
       *  replaced by the Riemann sum, as follows:
       *
       *  \f[ b_{eff}(z) \simeq \frac{\sum_{i_{min}}^{i_{max}} b(M, z)
       *  \Phi(M, z) f(M, z) (M_{i+1}-M_i)}{\sum_{i_{min}}^{i_{max}}
       *  \Phi(M, z) f(M, z) (M_{i+1}-M_i)} \f]
       *
       *  where the halo mass function, \f$\Phi(M, z)\f$, is computed
       *  by cbl::cosmology::MassFunction::mass_function, the linear
       *  bias, \f$b(M, z)\f$, is computed by
       *  cbl::cosmology::Bias::bias_halo, and f(M, z) is the
       *  selection function
       *
       *  @param interp_sigma FuncGrid object containing the values
       *  of &sigma;(mass) computed on a grid, used by interpolation
       *
       *  @param interp_DnSigma FuncGrid object containing the values
       *  of dln&sigma;/dM computed on a grid, used by interpolation
       *
       *  @param interp_SF FuncGrid2D object containing the values of
       *  the selection function computed on a grid in mass and
       *  redshift, used by interpolation
       *
       *  @param Mass_min minimum halo mass
       *
       *  @param Mass_max maximum halo mass
       *
       *  @param redshift vector containing the input redshifts
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press & Schechter), ST (Sheth &
       *  Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       *  al. 2006), ShenH (halo MF, Shen et al. 2006), ShenF
       *  (filament MF, Shen et al. 2006), ShenS (sheet MF, Shen et
       *  al. 2006), Pan (Pan 2007), Peacock (Peacock at al. 2007),
       *  Reed (Reed et al. 2007), Tinker (Tinker et al. 2008), Crocce
       *  (Crocce et al. 2010), Manera (Manera et al. 2010), Courtin
       *  (Courtin et al. 2010), Bhattacharya (Bhattacharya et
       *  al. 2011), Angulo_FOF (FoF MF, Angulo et al. 2012),
       *  Angulo_Sub (SUBFIND MF, Angulo et al. 2012), Watson_FOF (FoF
       *  MF, Watson et al. 2012), Watson_SOH (Spherical Overdensity
       *  halo MF, Watson et al. 2012), Despali_Z0, Despali_AllZ,
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param alpha the \f$\alpha\f$ parameter of the cluster mass
       *  scaling relation
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
       *  @param Delta_crit \f$\Delta_{crit}\f$: the critical
       *  overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
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
       *  @param is_parameter_file true \f$\rightarrow\f$ the
       *  input_file is a parameter file, used to compute the power
       *  spectrum with the method specified by method_Pk; false
       *  \f$\rightarrow\f$ the input_file is a file containing the
       *  power spectrum
       *
       *  @return b<SUB>eff</SUB>: the effective dark matter bias
       *
       *  @warning interp_sigma and interp_DnSigma have to be
       *  computed at the same cosmology of the object. They are
       *  provided as an input just to improve the performances in
       *  some applications (e.g. MCMC) where these quantities can be
       *  computed once
       */
      std::vector<double> bias_eff_selection_function (const glob::FuncGrid interp_sigma, const glob::FuncGrid interp_DnSigma, const glob::FuncGrid2D interp_SF, const double Mass_min, const double Mass_max, const std::vector<double> redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const double alpha=1., const bool store_output=true, const std::string output_root="test", const double Delta_crit=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
          
      /**
       *  @brief effective bias of dark matter haloes, computed using
       *  a given selection function
       *
       *  this function computes the effective bias of dark matter
       *  haloes:
       *
       *  \f[ b_{eff}(z) = \frac{\int_{M_{min}}^{M_{max}} {\rm d}M\,
       *  b(M, z) \Phi(M, z) f(M, z)}{\int_{M_{min}}^{M_{max}} {\rm
       *  d}M\,\Phi(M, z) f(M, z)} \f]
       *
       *  in the current implementation, the integral is actually
       *  replaced by the Riemann sum, as follows:
       *
       *  \f[ b_{eff}(z) \simeq \frac{\sum_{i_{min}}^{i_{max}} b(M, z)
       *  \Phi(M, z) f(M, z) (M_{i+1}-M_i)}{\sum_{i_{min}}^{i_{max}}
       *  \Phi(M, z) f(M, z) (M_{i+1}-M_i)} \f]
       *
       *  where the halo mass function, \f$\Phi(M, z)\f$, is computed
       *  by cbl::cosmology::MassFunction::mass_function, the linear
       *  bias, \f$b(M, z)\f$, is computed by
       *  cbl::cosmology::Bias::bias_halo, and f(M, z) is the
       *  selection function
       *
       *  @param Mass_min minimum halo mass
       *
       *  @param Mass_max maximum halo mass
       *
       *  @param redshift vector containing the input redshifts
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press & Schechter), ST (Sheth &
       *  Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       *  al. 2006), ShenH (halo MF, Shen et al. 2006), ShenF
       *  (filament MF, Shen et al. 2006), ShenS (sheet MF, Shen et
       *  al. 2006), Pan (Pan 2007), Peacock (Peacock at al. 2007),
       *  Reed (Reed et al. 2007), Tinker (Tinker et al. 2008), Crocce
       *  (Crocce et al. 2010), Manera (Manera et al. 2010), Courtin
       *  (Courtin et al. 2010), Bhattacharya (Bhattacharya et
       *  al. 2011), Angulo_FOF (FoF MF, Angulo et al. 2012),
       *  Angulo_Sub (SUBFIND MF, Angulo et al. 2012), Watson_FOF (FoF
       *  MF, Watson et al. 2012), Watson_SOH (Spherical Overdensity
       *  halo MF, Watson et al. 2012), Despali_Z0, Despali_AllZ,
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param selection_function_file input file with the selection
       *  functon
       *
       *  @param column vector containing the three columns of the
       *  selection function file to be read
       *
       *  @param alpha the \f$\alpha\f$ parameter of the cluster mass
       *  scaling relation
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
       *  @param Delta_crit \f$\Delta_{crit}\f$: the critical
       *  overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
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
       *  @param is_parameter_file true \f$\rightarrow\f$ the
       *  input_file is a parameter file, used to compute the power
       *  spectrum with the method specified by method_Pk; false
       *  \f$\rightarrow\f$ the input_file is a file containing the
       *  power spectrum
       *
       *  @return b<SUB>eff</SUB>: the effective dark matter bias
       */
      std::vector<double> bias_eff_selection_function (const double Mass_min, const double Mass_max, const std::vector<double> redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> column={}, const double alpha=1., const bool store_output=true, const std::string output_root="test", const double Delta_crit=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
                 
      /**
       *  @brief compute the effective bias of dark matter haloes, by
       *  averaging the bias of a set of haloes, interpolating the
       *  mass variance on a grid of masses and of one input
       *  cosmological parameter; this function is used when modelling
       *  the two-point correlation function
       *
       *  this function computes the effective bias of dark matter
       *  haloes by either averaging the bias of a set of haloes with
       *  a given mass:
       *
       *  \f[b_{eff}(z) = \frac{1}{N_{halo}}\sum_{i=1}^{N_{halo}}
       *  b(M_i, z_i) \; , \; (1)\f]
       *
       *  or by averaging over halo pairs:
       *
       *  \f[b_{eff}(z) = \sqrt{ \frac{2}{N_{halo}(N_{halo}-1)}
       *  \sum_{i=1}^{N_{halo}}\sum_{j=i+1}^{N_{halo}} b(M_i,
       *  z_i)b(M_j, z_j)} \; , \; (2)\f]
       *
       *  where the linear bias of the \f$i\f$-th halo, \f$b^{i}(M,
       *  z)\f$, is computed by cbl::cosmology::Bias::bias_halo; the
       *  mass variance is computed by cbl::cosmology::Sigma::sigma2M
       *
       *  @param parameter vector containing the grid of the
       *  cosmological parameters on which the effective bias grid is
       *  computed
       *
       *  @param bias_eff vector containing the effective bias grid
       *
       *  @param dir_output the directory where the effective bias
       *  grid is stored
       *
       *  @param file_bias_eff_grid the file there the effective bias
       *  grid is stored
       *
       *  @param cosmoPar the cosmological parameter for which the
       *  effective bias grid is computed
       *
       *  @param min_par the minimum value for the
       *  parameter where the effective bias is computed
       *  
       *  @param max_par the maximum value for the
       *  parameter where the effective bias is computed
       *
       *  @param nbin_par the number of points for the
       *  parameter where the effective bias is computed
       *
       *  @param mass vector containing the halo masses
       *
       *  @param mass_grid vector containing the halo masses on the
       *  grid used to interpolate the mass variance
       *
       *  @param redshift vector containing the redshifts; if it has
       *  size=1, it will be considered as the main redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param meanType meanType="mean_bias" \f$\rightarrow\f$ the
       *  effective bias is computed with Eq.(1);
       *  meanType="mean_pair_bias" \f$\rightarrow\f$ the effective
       *  bias is computed with Eq.(2)
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
       *  @param Delta_crit \f$\Delta_{crit}\f$: the critical
       *  overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
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
       *  @param cosmology_mass cosmology used to measure the cluster
       *  masses
       *
       *  @param redshift_source vector containing the redshifts of
       *  the source galaxies, in case the cluster masses are
       *  estimated from weak lensing
       */
      void generate_bias_eff_grid_one_cosmopar (std::vector<double> &parameter, std::vector<double> &bias_eff, const std::string dir_output, const std::string file_bias_eff_grid, const std::string cosmoPar, const double min_par, const double max_par, const int nbin_par, const std::vector<double> mass, const std::vector<double> mass_grid, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType="mean_bias", const bool store_output=true, const std::string output_root="test", const double Delta_crit=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true, const std::shared_ptr<cosmology::Cosmology> cosmology_mass=NULL, const std::vector<double> redshift_source={});
                
      /**
       *  @brief effective bias of dark matter haloes, computed by
       *  weighting on the selection function on a grid of one input
       *  cosmological parameter; this function is used when modelling
       *  the two-point correlation function
       *
       *  this function computes the effective bias of dark matter
       *  haloes:
       *
       *  \f[ b_{eff}(z) = \frac{\int_{M_{min}}^{M_{max}} {\rm d}M\,
       *  b(M, z) \Phi(M, z) f(M, z)}{\int_{M_{min}}^{M_{max}} {\rm
       *  d}M\,\Phi(M, z) f(M, z)} \f]
       *
       *  where the linear bias of the \f$i\f$-th halo, \f$b^{i}(M,
       *  z)\f$, is computed by cbl::cosmology::Bias::bias_halo; the
       *  mass variance is computed by cbl::cosmology::Sigma::sigma2M
       *
       *  @param parameter vector containing the grid of the
       *  cosmological parameters on which the effective bias grid is
       *  computed
       *
       *  @param bias_eff vector containing the effective bias grid
       *
       *  @param dir_output the directory where the effective bias
       *  grid is stored
       *
       *  @param file_bias_eff_grid the file there the effective bias
       *  grid is stored
       *
       *  @param cosmoPar the cosmological parameter for which the
       *  effective bias grid is computed
       *
       *  @param min_par the minimum value for the
       *  parameter where the effective bias is computed
       *  
       *  @param max_par the maximum value for the
       *  parameter where the effective bias is computed
       *
       *  @param nbin_par the number of points for the
       *  parameter where the effective bias is computed
       *
       *  @param redshift vector containing the redshifts; if it has
       *  size=1, it will be considered as the main redshift
       *
       *  @param Mass_min minimum cluster mass
       *
       *  @param Mass_max maximum cluster mass
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       * 
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press & Schechter), ST (Sheth &
       *  Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       *  al. 2006), ShenH (halo MF, Shen et al. 2006), ShenF
       *  (filament MF, Shen et al. 2006), ShenS (sheet MF, Shen et
       *  al. 2006), Pan (Pan 2007), Peacock (Peacock at al. 2007),
       *  Reed (Reed et al. 2007), Tinker (Tinker et al. 2008), Crocce
       *  (Crocce et al. 2010), Manera (Manera et al. 2010), Courtin
       *  (Courtin et al. 2010), Bhattacharya (Bhattacharya et
       *  al. 2011), Angulo_FOF (FoF MF, Angulo et al. 2012),
       *  Angulo_Sub (SUBFIND MF, Angulo et al. 2012), Watson_FOF (FoF
       *  MF, Watson et al. 2012), Watson_SOH (Spherical Overdensity
       *  halo MF, Watson et al. 2012), Despali_Z0, Despali_AllZ,
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param alpha the \f$\alpha\f$ parameter of the cluster mass
       *  scaling relation
       *
       *  @param selection_function_file the input selection function
       *  file
       *
       *  @param column vector containing the columns with {mass,
       *  redshift, selection function}
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
       *  @param Delta_crit \f$\Delta_{crit}\f$: the critical
       *  overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
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
       */
      void generate_bias_eff_grid_one_cosmopar (std::vector<double> &parameter, std::vector<double> &bias_eff, const std::string dir_output, const std::string file_bias_eff_grid, const std::string cosmoPar, const double min_par, const double max_par, const int nbin_par, const double redshift, const double Mass_min, const double Mass_max, const std::string model_bias, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> column={}, const double alpha=1., const bool store_output=true, const std::string output_root="test", const double Delta_crit=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
                
      /**
       *  @brief effective bias of dark matter haloes, computed by
       *  averaging the bias of a set of haloes, interpolating the
       *  mass variance on a grid of masses and two input cosmological
       *  parameters; this function is used when modelling the
       *  two-point correlation function
       *
       *  this function computes the effective bias of dark matter
       *  haloes by either averaging the bias of a set of haloes with
       *  a given mass:
       *
       *  \f[b_{eff}(z) = \frac{1}{N_{halo}}\sum_{i=1}^{N_{halo}}
       *  b(M_i, z_i) \; , \; (1)\f]
       *
       *  or by averaging over halo pairs:
       *
       *  \f[b_{eff}(z) = \sqrt{ \frac{2}{N_{halo}(N_{halo}-1)}
       *  \sum_{i=1}^{N_{halo}}\sum_{j=i+1}^{N_{halo}} b(M_i,
       *  z_i)b(M_j, z_j)} \; , \; (2)\f]
       *
       *  where the linear bias of the \f$i\f$-th halo, \f$b^{i}(M,
       *  z)\f$, is computed by cbl::cosmology::Bias::bias_halo; the
       *  mass variance is computed by cbl::cosmology::Sigma::sigma2M
       *
       *  @param parameter1 vector containing the grid of the first
       *  cosmological parameters on which the effective bias grid is
       *  computed
       *
       *  @param parameter2 vector containing the grid of the second
       *  cosmological parameters on which the effective bias grid is
       *  computed
       *
       *  @param bias_eff vector containing the effective bias grid
       *
       *  @param dir_output the directory where the effective bias
       *  grid is stored
       *
       *  @param file_bias_eff_grid the file there the effective bias
       *  grid is stored
       *
       *  @param cosmoPar1 the first cosmological parameter for which
       *  the effective bias grid is computed
       *
       *  @param min_par1 the minimum value for the first
       *  parameter where the effective bias is computed
       *  
       *  @param max_par1 the maximum value for the first
       *  parameter where the effective bias is computed
       *
       *  @param nbin_par1 the number of points for the first
       *  parameter where the effective bias is computed
       *
       *  @param cosmoPar2 the second cosmological parameter for which
       *  the effective bias grid is computed
       *
       *  @param min_par2 the minimum value for the second
       *  parameter where the effective bias is computed
       *  
       *  @param max_par2 the maximum value for the second
       *  parameter where the effective bias is computed
       *
       *  @param nbin_par2 the number of points for the second
       *  parameter where the effective bias is computed
       *
       *  @param mass vector containing the halo masses
       *
       *  @param mass_grid vector containing the halo masses on the
       *  grid used to interpolate the mass variance
       *
       *  @param redshift vector containing the redshifts; if it has
       *  size=1, it will be considered as the main redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param meanType meanType="mean_bias" \f$\rightarrow\f$ the
       *  effective bias is computed with Eq.(1);
       *  meanType="mean_pair_bias" \f$\rightarrow\f$ the effective
       *  bias is computed with Eq.(2)
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
       *  @param Delta \f$\Delta\f$, the overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
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
       *  @param cosmology_mass cosmology used to measure the cluster
       *  masses
       *
       *  @param redshift_source vector containing the redshifts of
       *  the source galaxies, in case the cluster masses are estimated
       *  from weak lensing
       */
      void generate_bias_eff_grid_two_cosmopars (std::vector<double> &parameter1, std::vector<double> &parameter2, std::vector<std::vector<double>> &bias_eff, const std::string dir_output, const std::string file_bias_eff_grid, const std::string cosmoPar1, const double min_par1, const double max_par1, const int nbin_par1, const std::string cosmoPar2, const double min_par2, const double max_par2, const int nbin_par2, const std::vector<double> mass, const std::vector<double> mass_grid, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType="mean_bias", const bool store_output=true, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true, const std::shared_ptr<cosmology::Cosmology> cosmology_mass=NULL, const std::vector<double> redshift_source={});
          
      /**
       *  @brief bias of dark matter haloes, for a vector of masses
       *
       *  @param Mass halo mass
       *
       *  @param Sigma &sigma;(mass, z=0): the mass variance at z=0
       *
       *  @param redshift the redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
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
       *  @param Delta \f$\Delta\f$, the overdensity
       *  
       *  @param kk wave vector module
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
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
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
       *  @return b<SUB>halo</SUB>: the dark matter bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ to obtain the background overdensity
       */
      std::vector<double> bias_halo (const std::vector<double> Mass, const std::vector<double> Sigma, const double redshift, const std::string model_bias, const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double Delta=200., const double kk=-1., const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string method_SS="CAMB", const std::string input_file=par::defaultString, const bool is_parameter_file=true);
      
      ///@}

    };
  }
}

#endif
