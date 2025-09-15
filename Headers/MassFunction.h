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
 *  @file Headers/MassFunction.h
 *
 *  @brief The class MassFunction
 *
 *  This file defines the interface of the class MassFunction, used to
 *  model the mass function of dark matter haloes
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __MASSFUNCTION__
#define __MASSFUNCTION__

#include "NG.h"
#include "Sigma.h"


// ===================================================================================================


namespace cbl {

  namespace cosmology {

    /**
     *  @class MassFunction MassFunction.h "Headers/MassFunction.h"
     *
     *  @brief The class MassFunction
     *
     *  This class is used to handle objects of type <EM> MassFunction
     *  </EM>. It is used to model the mass function of dark matter
     *  haloes
     */
    class MassFunction {

    private:
      
      /// pointer to the input cosmology
      std::shared_ptr<cosmology::Cosmology> m_cosmology = NULL;
      
      /**
       *  @brief auxiliary function to compute the mass function of
       *  dark matter haloes (filaments and sheets)
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass mass
       *
       *  @param mass_function_params function to a container of the
       *  mass function parameters
       *
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       */
      double m_mass_function (const double Mass, std::shared_ptr<void> mass_function_params) const;
          
      /**
       *  @brief auxiliary function to compute the mass function
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Sigma &sigma;(mass): the mass variance
       *
       *  @param Dln_Sigma dln&sigma;/dM: the derivative of the mass
       *  variance
       *
       *  @param redshift the redshift
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
       *  Castro (Castro et al. 2023)
       *
       *  @param Delta \f$\Delta\f$, the overdensity
       *
       *  @param default_delta true = using function
       *  cbl::cosmology::deltac; false = using delta_t*growth
       *  factor
       *  
       *  @param delta_t user defined density contrast at \f$z = 0\f$
       *
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double m_MF_generator (const double Sigma, const double Dln_Sigma, const double redshift, const std::string model_MF, const double Delta=200., const bool default_delta=true, const double delta_t=1.686)const ; 
          
      /**
       *  @brief auxiliary function to compute the mass function
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Sigma &sigma;(mass): the mass variance
       *
       *  @param Dln_Sigma dln&sigma;/dM: the derivative of the mass
       *  variance
       *
       *  @param redshift the redshift
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
       *  @param D_N the growth factor, precomputed. 
       *
       *  @param Delta \f$\Delta\f$, the overdensity
       *
       *  @param default_delta true = using function
       *  cbl::cosmology::deltac; false = using delta_t*growth
       *  factor
       *  
       *  @param delta_t user defined density contrast at \f$z = 0\f$
       *
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double m_MF_generator (const double Sigma, const double Dln_Sigma, const double redshift, const double D_N, const std::string model_MF, const double Delta = 200., const bool default_delta = true, const double delta_t = 1.686)const ;

    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      MassFunction () = default;

      /**
       *  @brief constructor
       *
       *  @param cosmology pointer to an object of class Cosmology
       */
      MassFunction (std::shared_ptr<cbl::cosmology::Cosmology> cosmology)
      : m_cosmology(std::move(cosmology)) {}
      
      /**
       *  @brief default destructor
       */
      ~MassFunction () = default;

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
       *  @name Functions to estimate the mass function of the dark
       *  matter haloes
       */
      ///@{

      /**
       *  @brief the mass function of dark matter haloes (filaments
       *  and sheets)
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass mass
       *
       *  @param redshift the redshift
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
       *  Castro (Castro et al. 2023)
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
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
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
       *  @param default_delta true = using function
       *  cbl::cosmology::deltac; false = using delta_t*growth
       *  factor
       *  
       *  @param delta_t user defined density contrast at \f$z = 0\f$
       *
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double mass_function (const double Mass, const double redshift, const std::string model_MF, const std::string method_SS, const bool store_output = true, const std::string output_root = "test", const double Delta = 200., const std::string interpType = "Linear", const int norm = -1, const double k_min = 0.001, const double k_max = 100., const double prec = 1.e-2, const std::string input_file = par::defaultString, const bool is_parameter_file = true, const bool default_delta = true, const double delta_t = 1.686)const ;

      /**
       *  @brief the mass function of dark matter haloes (filaments
       *  and sheets)
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass the vector of masses
       *
       *  @param redshift the redshift
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
       *  Castro (Castro et al. 2023)
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
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
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
       *  @param default_delta true = using function
       *  cbl::cosmology::deltac; false = using delta_t*growth
       *  factor
       *  
       *  @param delta_t user defined density contrast at \f$z = 0\f$
       *
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      std::vector<double> mass_function (const std::vector<double> Mass, const double redshift, const std::string model_MF, const std::string method_SS, const bool store_output=true, const std::string output_root="test", const double Delta=200., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true, const bool default_delta=true, const double delta_t=1.686) const;

      /**
       *  @brief the mass function of dark matter haloes in f(R)
       *  cosmologies (see Hu & Sawicki 2007) computed with the
       *  Boltzmann solver MGCAMB
       *
       *  @author Leonardo Gabriele Coppola and Sofia Contarini
       *
       *  @author leonardo.coppola@studio.unibo.it,
       *  sofia.contarini3@unibo.it
       *
       *  @param Mass the vector of masses
       *
       *  @param redshift the redshift
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
       *  Castro (Castro et al. 2023)
       *
       *  @param f_R0 value of the parameter \f$f_\mathrm{R0}\f$
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param Delta \f$\Delta\f$, the overdensity
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
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
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
       *  @param default_delta true = using function
       *  cbl::cosmology::deltac; false = using delta_t*growth
       *  factor
       *  
       *  @param delta_t user defined density contrast at \f$z = 0\f$
       *
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      std::vector<double> mass_function_fR (const std::vector<double> Mass, const double redshift, const std::string model_MF, const double f_R0=0., const bool store_output=true, const double Delta=200., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true, const bool default_delta=true, const double delta_t=1.686) const;
                
      /**
       *  @brief the mass function of dark matter haloes (filaments and
       *  sheets) computed quickly using a grid
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass mass
       *
       *  @param redshift the redshift
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
       *  Castro (Castro et al. 2023)
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
       *  name. . If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param Delta \f$\Delta\f$, the overdensity
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
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
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
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double mass_function_fast (const double Mass, const double redshift, const std::string model_MF, const std::string method_SS, const bool store_output=true, const std::string output_root="test", const double Delta=200., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true) const; 
          
      /**
       *  @brief the mass function of dark matter haloes (filaments and
       *  sheets) computed quickly passing directly the mass variance
       *  and its derivative as inputs
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass mass
       *
       *  @param Sigma &sigma;(mass): the mass variance
       *
       *  @param Dln_Sigma dln&sigma;/dM: the derivative of the mass
       *  variance
       *
       *  @param redshift the redshift
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
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
       *  @param Delta \f$\Delta\f$, the overdensity
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed to estimate the power spectrum
       *  normalisation; this parameter is used only if norm=1
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
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
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double mass_function (const double Mass, const double Sigma, const double Dln_Sigma, const double redshift, const std::string model_MF, const bool store_output=true, const std::string output_root="test", const double Delta=200., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string method_SS="CAMB", const std::string input_file=par::defaultString, const bool is_parameter_file=true) const; 
                
      /**
       *  @brief the mass function of dark matter haloes (filaments and
       *  sheets) computed quickly passing directly the mass variance
       *  and its derivative as inputs
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass mass
       *
       *  @param Sigma &sigma;(mass): the mass variance
       *
       *  @param Dln_Sigma dln&sigma;/dM: the derivative of the mass
       *  variance
       *
       *  @param redshift the redshift
       *
       *  @param D_N the amplitude of the growing mode
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
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
       *  @param Delta \f$\Delta\f$, the overdensity
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed to estimate the power spectrum
       *  normalisation; this parameter is used only if norm=1
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
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
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double mass_function (const double Mass, const double Sigma, const double Dln_Sigma, const double redshift, const double D_N, const std::string model_MF, const bool store_output=true, const std::string output_root="test", const double Delta=200., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string method_SS="CAMB", const std::string input_file=par::defaultString, const bool is_parameter_file=true) const; 
          
      /**
       *  @brief number of dark matter haloes per steradian or square
       *  degree, for a given redshift range
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass_min minimum mass
       *
       *  @param Mass_max maximum mass
       *
       *  @param z_min minimum redshift
       *
       *  @param z_max maximum redshift
       *
       *  @param angle_rad 0 \f$\rightarrow\f$ &Omega; in square
       *  degrees; 1 \f$\rightarrow\f$ &Omega; in steradians
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
       *  Castro (Castro et al. 2023)
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
       *  @return n<SUB>haloes</SUB>: the number density of dark matter
       *  haloes (per steradian or square degree)
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double n_haloes (const double Mass_min, const double Mass_max, const double z_min, const double z_max, const bool angle_rad, const std::string model_MF, const std::string method_SS, const bool store_output=true, const std::string output_root="test", const double Delta=200, const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;
                
      /**
       *  @brief number of dark matter haloes per volume at fixed
       *  redshift 
       *
       *  this function computes the number of dark matter haloes per
       *  volume at fixed redshift as follows:
       *
       *  \f[ N_h = \int_{M_{min}}^{M_{max}} d M \Phi(M)\f]
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass_min minimum mass
       *
       *  @param Mass_max maximum mass
       *
       *  @param Volume the volume
       *
       *  @param redshift the redshift
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
       *  Castro (Castro et al. 2023)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param nbin_mass number of bin for the mass function
       *  computation
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
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed to estimate the power spectrum
       *  normalisation; this parameter is used only if norm=1
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
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
       *  @param default_delta true = using function
       *  cbl::cosmology::deltac; false = using delta_t*growth
       *  factor
       *  
       *  @param delta_t user defined density contrast at \f$z = 0\f$
       *
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double n_haloes (const double Mass_min, const double Mass_max, const double Volume, const double redshift, const std::string model_MF, const std::string method_SS, const int nbin_mass=0, const bool store_output=true, const std::string output_root="test", const double Delta=200., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true, const bool default_delta=true, const double delta_t=1.686) const;
                
      /**
       *  @brief minimum halo mass, given the number of haloes in a
       *  given region of sky
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param n_halo number density of dark matter haloes
       *
       *  @param Area sky area
       *
       *  @param angle_rad 0 \f$\rightarrow\f$ &Omega; in square
       *  degrees; 1 \f$\rightarrow\f$ &Omega; in steradians
       *
       *  @param z_min minimum redshift
       * 
       *  @param z_max maximum redshift
       *
       *  @param Mmax maximum mass
       *  
       *  @param lgM1_guess logarithm of the minimum mass used by the
       *  root finder
       *
       *  @param lgM2_guess logarithm of the maximum mass used by the
       *  root finder
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
       *  Castro (Castro et al. 2023)
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
       *  @return minimum halo mass
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double MhaloMin (const int n_halo, const double Area, const bool angle_rad, const double z_min, const double z_max, const double Mmax, const double lgM1_guess, const double lgM2_guess, const std::string model_MF, const std::string method_SS, const bool store_output=true, const std::string output_root="test", const double Delta=200, const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;
          
      /**
       *  @brief the unevolved mass function
       *  @author Carlo Giocoli
       *  @author cgiocoli@gmail.com
       *  @param mass_accr mass accreted
       *  @return the unevolved mass function
       */
      double unevolved_mass_function (const double mass_accr) const; 
          
      /**
       *  @brief number of dark matter haloes per steradian or square
       *  degree, for a given redshift range and with selection function
       *  defined on a grid
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass_min minimum mass
       *
       *  @param Mass_max maximum mass
       *
       *  @param z_min minimum redshift
       *
       *  @param z_max maximum redshift
       *
       *  @param angle_rad 0 \f$\rightarrow\f$ &Omega; in square
       *  degrees; 1 \f$\rightarrow\f$ &Omega; in steradians
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
       *  Castro (Castro et al. 2023)
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  CLASS [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param selection_function_file input file where the selection
       *  function is stored
       *
       *  @param column the columns to be read
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
       *  @param isDelta_critical \f$\rightarrow\f$ \f$\Delta\f$ is the
       *  overdensity defined with respect to the critical density
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
       *  @return n<SUB>haloes</SUB>: the number density of dark matter
       *  haloes (per steradian or square degree)
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double n_haloes_selection_function (const double Mass_min, const double Mass_max, const double z_min, const double z_max, const bool angle_rad, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> column={}, const bool store_output=true, const std::string output_root="test", const double Delta=200, const bool isDelta_critical=false, const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;
          
      /**
       *  @brief mass function for a range of masses
       *
       *  @author Alfonso Veropalumbo
       *
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param mass vector of mass
       *
       *  @param z_min minimum redshift
       *
       *  @param z_max maximum redshift
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
       *  Castro (Castro et al. 2023)
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
       *  @param isDelta_critical \f$\rightarrow\f$ \f$\Delta\f$ is the
       *  overdensity defined with respect to the critical density
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
       *  @return n<SUB>haloes</SUB>: the number density of dark matter
       *  haloes (per steradian or square degree)
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      std::vector<double> mass_function (const std::vector<double> mass, const double z_min, const double z_max, const std::string model_MF, const std::string method_SS, const bool store_output=true, const std::string output_root="test", const double Delta=200, const bool isDelta_critical=false, const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;
          
      /**
       *  @brief mass function given a selection function
       *
       *  @author Alfonso Veropalumbo
       *
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param mass vector of mass
       *
       *  @param z_min minimum redshift
       *
       *  @param z_max maximum redshift
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
       *  Castro (Castro et al. 2023)
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  CLASS [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param selection_function_file input file where the selection
       *  function is stored
       *
       *  @param column the columns to be read
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
       *  @param isDelta_critical \f$\rightarrow\f$ \f$\Delta\f$ is the
       *  overdensity defined with respect to the critical density
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
       *  @return n<SUB>haloes</SUB>: the number density of dark matter
       *  haloes (per steradian or square degree)
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      std::vector<double> mass_function_selection_function_vector (const std::vector<double> mass, const double z_min, const double z_max, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> column={}, const bool store_output=true, const std::string output_root="test", const double Delta=200, const bool isDelta_critical=false, const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;
          
      /**
       *  @brief redshift distribution of dark matter haloes
       *
       *  @author Alfonso Veropalumbo
       *
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param z_min minimum redshift
       *
       *  @param z_max maximum redshift
       *
       *  @param step_z redshift step
       *
       *  @param Area_degrees the survey area, in degrees
       *
       *  @param Mass_min minimum halo mass
       *
       *  @param Mass_max maximum halo mass
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
       *  Castro (Castro et al. 2023)
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
       *  @param isDelta_critical \f$\rightarrow\f$ \f$\Delta\f$ is the
       *  overdensity defined with respect to the critical density
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
       *  @return the redshift distribution of dark matter haloes
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      std::vector<double> redshift_distribution_haloes (const double z_min, const double z_max, const int step_z, const double Area_degrees, const double Mass_min, const double Mass_max, const std::string model_MF, const std::string method_SS, const bool store_output=true, const std::string output_root="test", const double Delta=200, const bool isDelta_critical=false, const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;
          
      /**
       *  @brief redshift distribution of dark matter haloes, given a
       *  selection function
       *
       *  @author Alfonso Veropalumbo
       *
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param redshift vector containing the redshift at which the
       *  halo distribution will be computed
       *
       *  @param Area_degrees the survey area, in degrees
       *
       *  @param Mass_min minimum halo mass
       *
       *  @param Mass_max maximum halo mass
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
       *  Castro (Castro et al. 2023)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param selection_function_file input file where the selection
       *  function is stored
       *
       *  @param column the columns to be read
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
       *  @param isDelta_critical \f$\rightarrow\f$ \f$\Delta\f$ is the
       *  overdensity defined with respect to the critical density
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
       *  @return the redshift distribution of dark matter haloes
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      std::vector<double> redshift_distribution_haloes_selection_function (const std::vector<double> redshift, const double Area_degrees, const double Mass_min, const double Mass_max, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> column={}, const bool store_output=true, const std::string output_root="test", const double Delta=200, const bool isDelta_critical=false, const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;
          
      /**
       *  @brief the mean redshift of a dark matter haloe sample,
       *  given a selection function
       *
       *  @author Alfonso Veropalumbo
       *
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param z_min minimum redshift
       *
       *  @param z_max maximum redshift
       *
       *  @param Mass_min minimum halo mass
       *
       *  @param Mass_max maximum halo mass
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
       *  Castro (Castro et al. 2023)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param selection_function_file input file where the selection
       *  function is stored
       *
       *  @param column the columns to be read
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
       *  @param isDelta_critical \f$\rightarrow\f$ \f$\Delta\f$ is the
       *  overdensity defined with respect to the critical density
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
       *  @return the mean redshift of dark matter haloes
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double mean_redshift_haloes_selection_function (const double z_min, const double z_max, const double Mass_min, const double Mass_max, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> column={}, const bool store_output=true, const std::string output_root="test", const double Delta=200, const bool isDelta_critical=false, const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;
          
      /**
       *  @brief convert a cluster mass estimated in a different
       *  cosmology
       *
       *  this function converts a cluster mass estimated assuming a
       *  different cosmology, following Eq.C4 of Sereno & Ettori 2015
       *  (https://arxiv.org/abs/1407.7868)
       *
       *  \f[ M^{(2)} = M^{(1)} \left( \frac{D_{ds}^{(1)}}{D_s^{(1)}}
       *  \right)^{3/2} H^{(1)} \left(
       *  \frac{D_{ds}^{(2)}}{D_s^{(2)}}\right)^{-3/2} H^{(2) -1} \f]
       *
       *  where and \f$D_s\f$ and \f$D_{ds}\f$ are the source and the
       *  lens-source angular diameter distances, respectively
       *  
       *  @param mass the cluster mass to be converted (since
       *  estimated in a different cosmology), \f$M^{(1)}\f$
       *
       *  @param assumed_cosmology pointer to the cosmology assumed to
       *  measure the cluster mass
       * 
       *  @param redshift redshift of the cluster
       *
       *  @param redshift_source redshift of the source, if the
       *  cluster mass is estimated from weak lensing, -1 otherwise
       *
       *  @return the cluster mass converted in this cosmology,
       *  \f$M^{(2)}\f$
       */
      double converted_mass (const double mass, const std::shared_ptr<cbl::cosmology::Cosmology> assumed_cosmology, const double redshift, const double redshift_source=-1.) const;
                
      /**
       *  @brief the mass function of dark matter haloes (filaments and
       *  sheets) computed quickly passing directly the mass variance
       *  and its derivative as inputs. This function takes vector in input
       *  and compute the mass function for the input masses
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass mass
       *
       *  @param Sigma &sigma;(mass): the mass variance
       *
       *  @param Dln_Sigma dln&sigma;/dM: the derivative of the mass
       *  variance
       *
       *  @param redshift the redshift
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
       *  Despali_AllZAllCosmo, Despali_HighM (Despali et al. 2016),
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
       *  @param Delta \f$\Delta\f$, the overdensity
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed to estimate the power spectrum
       *  normalisation; this parameter is used only if norm=1
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
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
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; divide \f$\Delta_{crit}\f$ by
       *  \f$\Omega_M(z)\f$ (see cbl::cosmology::Cosmology::OmegaM) to
       *  obtain the background overdensity
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      std::vector<double> mass_function (const std::vector<double> Mass, const std::vector<double> Sigma, const std::vector<double> Dln_Sigma, const double redshift, const std::string model_MF, const bool store_output=true, const std::string output_root="test", const double Delta=200., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string method_SS="CAMB", const std::string input_file=par::defaultString, const bool is_parameter_file=true) const; 

      ///@}
      
    };
  }

  namespace glob {

    struct STR_MF
    {
      double redshift;
      std::string model_MF;
      std::string method_SS;
      bool store_output;
      std::string output_root;
      double Delta;
      std::string interpType;
      int norm; 
      double k_min; 
      double k_max;
      double prec;
      std::string input_file;
      bool is_parameter_file;
      bool default_delta;
      double delta_t;
    };

  } 
}

#endif
