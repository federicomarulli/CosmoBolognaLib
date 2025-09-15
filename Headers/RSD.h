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
 *  @file Headers/RSD.h
 *
 *  @brief The class RSD
 *
 *  This file defines the interface of the class RSD, used to model
 *  redshift-space distortions
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __RSD__
#define __RSD__

#include "Bias.h"


// ===================================================================================================


namespace cbl {
  
  namespace cosmology {

    /**
     *  @class RSD RSD.h "Headers/RSD.h"
     *
     *  @brief The class RSD
     *
     *  This class is used to handle objects of type <EM> RSD
     *  </EM>. It is used to model redshift-space distortions
     */
    class RSD {

    private:

      /// pointer to the input cosmology
      std::shared_ptr<cosmology::Cosmology> m_cosmology = NULL;
      
      /// object of class bias
      cbl::cosmology::Bias m_bias;
            
    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      RSD () = default;
      
      /**
       *  @brief constructor
       *
       *  @param cosmology pointer to an object of class Cosmology
       */
      RSD (std::shared_ptr<cbl::cosmology::Cosmology> cosmology)
      : m_cosmology(std::move(cosmology)) {}
      
      /**
       *  @brief default destructor
       */
      ~RSD () = default;
      
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
       *  @name Functions to model redshift-space distortions
       */
      ///@{
      
      /**
       *  @brief f*&sigma;<SUB>8</SUB>: the linear growth rate times
       *  the dark matter rms mass fluctuation within 8 Mpc/h
       *
       *  @param redshift the redshift
       *
       *  @param method_Pk method used to compute the power spectrum
       *  and &sigma;(mass) (i.e. the Boltzmann solver); valid choices
       *  for method_Pk are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
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
       *  @param NL false \f$\rightarrow\f$ linear power spectrum;
       *  false \f$\rightarrow\f$ non-linear power spectrum
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
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return f*&sigma;<SUB>8</SUB>
       */
      double fsigma8 (const double redshift, const std::string method_Pk, const bool store_output=true, const std::string output_root="test", const bool NL=false, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString) const;
          
      /**
       *  @brief the specific growth rate &beta;
       *  @param redshift the redshift
       *  @param bias bias
       *  @return &beta;=f/b, where f is the linear growth rate and b is
       *  the bias
       */
      double beta (const double redshift, const double bias) const;
          
      /**
       *  @brief the error on the specific growth rate &beta;
       *  @param redshift the redshift
       *  @param bias bias
       *  @param err_bias error on the bias
       *  @return error on &beta;=f/b, where f is the linear growth rate
       *  and b is the bias
       */
      double error_beta (const double redshift, const double bias, const double err_bias) const;
          
      /**
       *  @brief the error on the specific growth rate &beta;
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
       *  @return &beta;=f/b, where f is the linear growth rate and b is
       *  the bias
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
      double beta (const double Mass_min, const double Mass_max, const double redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const bool store_output=true, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
          
      /**
       *  @brief the specific growth rate &beta;
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
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  CLASS [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param err_bias error on the bias
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
       *  @return error on &beta;=f/b, where f is the linear growth
       *  rate and b is the bias
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
      double error_beta (const double Mass_min, const double Mass_max, const double redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const double err_bias, const bool store_output=true, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true); 
          
      /**
       *  @brief the specific growth rate &beta;
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
       *  @return &beta;=f/b, where f is the linear growth rate and b is
       *  the bias
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
      double beta (const std::vector<double> MM, const std::vector<double> MF, const double redshift, const std::string model_bias, const std::string method_SS, const bool store_output=true, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
          
      /**
       *  @brief the error on the specific growth rate &beta;
       *
       *  @param MM vector of halo masses
       *
       *  @param MF vector of mass function values, d&Phi;/dM=dn(M)/dM
       *
       *  @param redshift the redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid authors
       *  are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen
       *  2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction
       *  of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  CLASS [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param err_bias error on the bias
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
       *  @return error on &beta;=f/b, where f is the linear growth
       *  rate and b is the bias
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
      double error_beta (const std::vector<double> MM, const std::vector<double> MF, const double redshift, const std::string model_bias, const std::string method_SS, const double err_bias, const bool store_output=true, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
          
      /**
       *  @brief the error on the specific growth rate &beta; from
       * Bianchi et al. 2012
       *
       *  @param Volume comoving volume 
       *
       *  @param density comoving density
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
       *  @return error on &beta;=f/b, where f is the linear growth rate
       *  and b is the bias
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
      double error_beta_measured (const double Volume, const double density, const double Mass_min, const double Mass_max, const double redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const bool store_output=true, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true); 
          
      /**
       *  @brief the normalised quadrupole Q
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
       *  @return Q: the normalised quadrupole
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
      double quadrupole (const double Mass_min, const double Mass_max, const double redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const bool store_output=true, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true); 
          
      /**
       *  @brief the normalised quadrupole Q
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
       *  @return Q: the normalised quadrupole
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ to obtain the background overdensity
       */
      double quadrupole (const std::vector<double> MM, const std::vector<double> MF, const double redshift, const std::string model_bias, const std::string method_SS, const bool store_output=true,  const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
          
      /**
       *  @brief the non-linear dark matter power spectrum using
       *  fitting functions given by Bel et. al (2019)
       *
       *  the non-linear dark matter power spectrum,
       *  P<SUB>&delta;&delta;</SUB>, is computed using the fitting
       *  functions given by Bel et. al (2019)
       *  [https://arxiv.org/abs/1809.09338].
       *
       *  @author J.E. Garcia-Farieta
       *  @author joegarciafa@unal.edu.co
       *
       *  @param kk the wave vector module
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
       *  @param author author(s) who proposed the fitting functions; valid
       *  authors are: Pezzotta (Pezzotta, et.al, 2017), Bel (Bel et.al., 2019)
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
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return P<SUB>&delta;&delta;</SUB>: the non linear power
       *  spectrum using fitting functions
       */
      double Pk_DeltaDelta_fitting_function (const double kk, const std::string method_Pk, const double redshift, const std::string author, const bool store_output, const std::string output_root, const int norm, double k_min, double k_max, const double prec, const std::string file_par, const bool unit1);
          
      /**
       *  @brief the dark matter cross power spectrum 
       *
       *  the dark matter cross power spectrum,
       *  P<SUB>&delta;&theta;</SUB>, is computed using the fitting
       *  functions given by Bel et. al (2019)
       *  [https://arxiv.org/abs/1809.09338].
       *
       *  \f[ P_{\delta\theta} (k) =
       *  \left(P_{\delta\delta}(k)P_{\theta\theta}^{Lin}(k)\right)^{1/2}e^{-k/k_\delta-bk^6}
       *  \f]
       *
       *  with \f$b = 0.091 + 0.702*\sigma8^2\f$
       *
       *  @author J.E. Garcia-Farieta
       *  @author joegarciafa@unal.edu.co
       *
       *  @param kk the wave vector module
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
       *  @param author author(s) who proposed the fitting functions; valid
       *  authors are: Pezzotta (Pezzotta, et.al, 2017), Bel (Bel et.al., 2019)
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
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$ non-linear
       *  power spectrum
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
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return P<SUB>&delta;&theta;</SUB>: the cross power spectrum using fitting functions
       */
      double Pk_DeltaTheta_fitting_function (const double kk, const std::string method_Pk, const double redshift, const std::string author, const bool store_output, const std::string output_root, const bool NL, const int norm, double k_min, double k_max, const double prec, const std::string file_par, const bool unit1);
          
      /**
       *  @brief the dark matter velocity divergence power spectrum 
       *
       *  the dark matter velocity divergence power spectrum,
       *  P<SUB>&theta;&theta;</SUB>, is computed using the fitting
       *  functions given by Bel et. al (2019)
       *  [https://arxiv.org/abs/1809.09338].
       *
       *  \f[ P_{\theta\theta}(k) =
       *  P^{Lin}_{\theta\theta}(k)e^{-k(a_1+a_2k+a_3k^2)} \f]
       *
       *  with
       *
       *  \f[a_1 = -0.817 + 3.198*\sigma_8\f] \f[a_2 = 0.877 -
       *  4.191*\sigma_8\f] \f[a_3 = -1.199 + 4.629*\sigma_8\f]
       *
       *  @author J.E. Garcia-Farieta
       *  @author joegarciafa@unal.edu.co
       *
       *  @param kk the wave vector module
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
       *  @param author author(s) who proposed the fitting functions;
       *  valid authors are: Pezzotta (Pezzotta, et.al, 2017), Bel
       *  (Bel et.al., 2019)
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
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return P<SUB>&theta;&theta;</SUB>: the velocity divergence
       *  power spectrum using fitting functions
       */
      double Pk_ThetaTheta_fitting_function (const double kk, const std::string method_Pk, const double redshift, const std::string author, const bool store_output, const std::string output_root, const bool NL, const int norm, double k_min, double k_max, const double prec, const std::string file_par, const bool unit1);
          
      /**
       *  @brief the linear-order one-dimensional pairwise velocity
       *  dispersion, \f$\sigma_{\mathrm{v}, \mathrm{lin}}\f$
       *
       *  This function computes the linear-order one-dimensional
       *  pairwise velocity dispersion, \f$\sigma_{\mathrm{v},
       *  \mathrm{lin}}\f$, as defined in Taruya et al. (2010) (Eq. 7
       *  of https://arxiv.org/abs/1006.0699)
       *
       *  \f[ \sigma_{\mathrm{v}, \mathrm{lin}}^{2} = \frac{1}{3} \int
       *  \frac{d^{3} \boldsymbol{q}}{(2 \pi)^{3}}
       *  \frac{P_{\mathrm{lin}}(q, z)}{q^{2}} \f]
       *
       *  where \f$P_{\mathrm{lin}}(q, z)\f$ is the linear power
       *  spectrum at a given redshift
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
       *  power spectrum is computed to estimate the power spectrum
       *  normalisation; this parameter is used only if norm=1
       *
       *  @param bin_k number of wave vector modules used for the
       *  integration
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return \f$\sigma_{\mathrm{v}, \mathrm{lin}}\f$
       */
      double sigma_v (const double redshift=0., const std::string method_Pk="CAMB", const bool store_output=true, const std::string output_root="test", const int norm=-1, const double k_min=0.001, const double k_max=100., const int bin_k=512, const double prec=1.e-2, const std::string file_par=par::defaultString, const bool unit1=false);
      
      ///@}

    };
  }
}

#endif
