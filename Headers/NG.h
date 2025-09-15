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
 *  @file Headers/NG.h
 *
 *  @brief The class NF
 *
 *  This file defines the interface of the class NG, used to model
 *  correction factors to be used in cosmologies with primordial non-Gaussianity
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __NG__
#define __NG__

#include "Sigma.h"


// ===================================================================================================


namespace cbl {

  namespace cosmology {

    /**
     *  @class NG NG.h "Headers/NG.h"
     *
     *  @brief The class NG
     *
     *  This class is used to handle objects of type <EM> NG </EM>. It
     *  is used to model correction factors to be used in cosmologies
     *  with primordial non-Gaussianity
     */
    class NG {

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
      NG () = default;
      
      /**
       *  @brief constructor
       *
       *  @param cosmology pointer to an object of class Cosmology
       */
      NG (std::shared_ptr<cbl::cosmology::Cosmology> cosmology)
      : m_cosmology(std::move(cosmology)) {}
      
      /**
       *  @brief default destructor
       */
      ~NG () = default;
          
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
       *  @name Functions to estimate the correction factors to be
       *  used in cosmologies with primordial non-Gaussianity
       */
      ///@{
      
      /**
       *  @brief the amplitude of the matter power spectrum
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
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
       *  @return A<SUB>m</SUB>
       */
      double Am (const std::string method_Pk, const bool store_output=true, const std::string output_root="test", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString); 
          
      /**
       *  @brief the potential spectral amplitude 
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
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
       *  @return the potential spectral amplitude
       */
      double potential_spectral_amplitude (const std::string method_Pk, const bool store_output=true, const std::string output_root="test", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString);
          
      /**
       *  @brief the bispectrum
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
       *
       *  @param kk wave vector module
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
       *  @return the potential spectral amplitude
       */
      double bispectrum (const std::vector<double> kk, const std::string method_Pk, const bool store_output=true, const std::string output_root="test", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString);
              
      /**
       *  @brief auxiliary function to estimate cosmological quantities
       *  in non-Gaussian cosmologies
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
       *
       *  @param kk wave vector module
       *
       *  @param mass halo mass
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
       *  @return mrk
       */
      double mrk (const double kk, const double mass, const std::string method_Pk, const bool store_output=true, const std::string output_root="test", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString);
          
      /**
       *  @brief auxiliary function to estimate cosmological quantities
       *  in non-Gaussian cosmologies
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
       *
       *  @param kk wave vector module
       *
       *  @param mass halo mass
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
       *  @return frk
       */
      double frk (const double kk, const double mass, const std::string method_Pk, const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
          
      /**
       *  @brief correction to the halo bias in non-Gaussian cosmologies
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
       *
       *  @param kk wave vector module
       *
       *  @param mass halo mass
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
       *  @return bias correction
       */
      double bias_correction (const double kk, const double mass, const std::string method_Pk, const bool store_output=true, const std::string  output_root="test", const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
          
      /**
       *  @brief the skewness
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
       *
       *  @param mass halo mass
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
       *  @return skewness
       */
      double skewness (const double mass, const std::string method_Pk, const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
          
      /**
       *  @brief the derivative of the skewness, ds/dM
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
       *
       *  @param mass halo mass
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
       *  @return derivative of the skewness
       */
      double dskewnessdM (const double mass, const std::string method_Pk, const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
          
      /**
       *  @brief correction to the halo mass in non-Gaussian cosmologies
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
       *
       *  @param mass the halo mass
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
       *  @param interpType method to interpolate the power spectrum
       *   
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum
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
       *  @return bias correction
       */
      double MF_correction (const double mass, const double redshift, const std::string method_Pk, const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      ///@}
      
    };
  

    /// @cond glob

    double bias_kernel (double, void *); // non-Gaussian cosmologies

    double skewness_kernel (double *, size_t, void *); // non-Gaussian cosmologies

    /// @endcond

    
    struct GSL_f_pars {
      double kt;
      double mass;
      std::string method_Pk;
      bool store_output;
      std::string output_root;
      int norm;
      double k_min;
      double k_max;
      bool GSL;
      double prec;
      std::string file_par;
      cosmology::Cosmology *pt_Cosmology;
    };
    
    struct STR_NG
    {
      double Omega_matter;
      double Omega_baryon;
      double Omega_neutrinos;
      double massless_neutrinos;
      double massive_neutrinos;
      double Omega_DE;
      double Omega_radiation;
      double hh;
      double scalar_amp;
      double scalar_pivot;
      double n_spec;
      double w0;
      double wa;
      double fNL;
      int type_NG;
      double tau;
      std::string model;
      bool unit;
      double kt;
      double mass;
      std::string method_Pk;
      bool store_output;
      std::string output_root;
      int norm;
      double k_min;
      double k_max;
      bool GSL;
      double prec;
      std::string file_par;
    };

  }
}

#endif
