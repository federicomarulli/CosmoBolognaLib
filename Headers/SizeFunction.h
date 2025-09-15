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
 *  @file Headers/SizeFunction.h
 *
 *  @brief The class SizeFunction
 *
 *  This file defines the interface of the class SizeFunction, used to
 *  model the size function of cosmic voids
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __SIZEFUNCTION__
#define __SIZEFUNCTION__

#include "MassFunction.h"


// ===================================================================================================


namespace cbl {
  
  namespace cosmology {
    
    /**
     *  @class SizeFunction SizeFunction.h "Headers/SizeFunction.h"
     *
     *  @brief The class SizeFunction
     *
     *  This class is used to handle objects of type <EM> SizeFunction
     *  </EM>. It is used to model the size function of cosmic voids
     */
    class SizeFunction {

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
      SizeFunction () = default;
      
      /**
       *  @brief constructor
       *
       *  @param cosmology pointer to an object of class Cosmology
       */
      SizeFunction (std::shared_ptr<cbl::cosmology::Cosmology> cosmology)
      : m_cosmology(std::move(cosmology)) {}
      
      /**
       *  @brief default destructor
       */
      ~SizeFunction () = default;

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
       *  @name Functions to estimate the size function of cosmic
       *  voids
       */
      ///@{
      
      /**
       *  @brief Linear (under)density contrast
       *
       *  @author Tommaso Ronconi
       *  @author tommaso.ronconi@studio.unibo.it
       *
       *  @param deltav_NL the non linear density contrast:
       *  \f$\rho_v/\rho_m\f$ (default value set to \f$-0.795\f$)
       *
       *  @param b_eff the effective bias of the sample
       *
       *  @param slope first coefficent to convert the effective bias
       *  (default value set to \f$0.854\f$)
       *
       *  @param offset second coefficent to convert the effective
       *  bias (default value set to \f$0.420\f$)
       *
       *  @return The linear density contrast used as second barrier
       *  in the excursion set formalism for voids, based on the fit
       *  by Bernardeu (1994): \f$\delta_v^L \equiv \frac{\rho_v -
       *  \rho_m}{\rho_m} \approx C [1 - (\rho_v/\rho_m)^{- 1/C}]\f$
       *  where \f$\rho_v =\ \f$ average void density, \f$\rho_m =\
       *  \f$ average density of the surrounding Universe and \f$C =
       *  1.594\f$, a costant.
       */
      double deltav_L (const double deltav_NL, const double b_eff, double slope=0.854, double offset=0.420) const;

      /**
       *  @brief Non-Linear (under)density contrast
       *
       *  @author Tommaso Ronconi
       *  @author tommaso.ronconi@studio.unibo.it
       *
       *  @param deltav the linear density contrast: \f$\delta_v\f$
       *  (default value set to \f$-2.71\f$)
       *  
       *  @return The non linear density contrast to be used in
       *  samples of tracers, independently of the bias value, based
       *  on the fit by Bernardeu (1994): \f$ \delta_v^{NL} \equiv
       *  \frac{\rho_v}{\rho_m} - 1 \approx \bigl(1 -
       *  C^{-1}\delta_v^L\bigr)^{-C} - 1\f$ where \f$\rho_v =\ \f$
       *  average void density, \f$\rho_m =\ \f$ average density of
       *  the surrounding Universe and \f$C = 1.594\f$, a costant.
       */
      double deltav_NL (const double deltav=-2.71) const;

      /**
       *  @brief expansion factor
       *
       *  @author Tommaso Ronconi
       *  @author tommaso.ronconi@studio.unibo.it
       *
       *  @param deltav the linear density contrast: \f$\delta_v^L\f$
       *  (default value set to \f$-2.71\f$)
       *  
       *  @return the expansion factor: \f$\frac{r}{r_L} = \bigl((1 -
       *  C^{-1}\delta_v^L\bigr)^{C/3}\f$ where \f$C = 1.594\f$
       */
      double r_rL (const double deltav = -2.71) const;

      /**
       *  @brief \f$f_{\ln \sigma}(\sigma)\f$ (approximation)
       *
       *  @author Tommaso Ronconi
       *  @author tommaso.ronconi@studio.unibo.it
       *
       *  @param SS variance of the linear density field
       *  (\f$\sigma^2(R)\f$)
       *
       *  @param del_v linear density contrast defining a void
       *
       *  @param del_c critical value of the linear density field
       *  
       *  @return the fraction of trajectories that evolve into voids,
       *  as given in equation (8) of Jennings et al. (2013)
       */
      double f_nu (const double SS, const double del_v, const double del_c) const;

      /**
       *  @brief Supplementary function to compute a correction factor
       *  to apply to the void size function, to predict the measured
       *  abundance of voids affected by geometrical distortions
       *  (Alcock-Paczyński change of volume). These are caused by the
       *  assumption of a fiducial cosmology different from the true
       *  one, see Correa et al. 2020
       *  (https://arxiv.org/pdf/2007.12064.pdf)
       *
       *  @author Sofia Contarini
       *  @author sofia.contarini3@unibo.it
       *
       *  @param cosm_true the true cosmology of the catalogue
       *
       *  @param redshift the vector of redshifts of different samples
       *  of voids
       *
       *  @return a vector with the correction factors relative to the
       *  void samples at different redshifts
       */
      std::vector<double> AP_corr (const std::shared_ptr<cosmology::Cosmology> cosm_true, const std::vector<double> redshift);

      /**
       *  @brief the void size function
       *
       *  @author Tommaso Ronconi
       *  @author tommaso.ronconi@studio.unibo.it
       *
       *  @param RV radius
       *
       *  @param redshift the redshift
       *
       *  @param model size function model name; valid choices for
       *  model name are SvdW (Sheth and van de Weygaert, 2004),
       *  linear and Vdn (Jennings et al., 2013)
       *
       *  @param b_eff the effective bias of the sample
       *
       *  @param slope first coefficent to convert the effective bias
       *  (default value set to \f$0.854\f$)
       *
       *  @param offset second coefficent to convert the effective
       *  bias (default value set to \f$0.420\f$)
       *
       *  @param deltav_NL the non linear density contrast:
       *  \f$\rho_v/\rho_m\f$ (default value set to \f$-0.795\f$)
       *
       *  @param del_c critical value of the linear density field
       *  (default value set to \f$1.06\f$)
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
       *  @return the number density of voids as a function of radius.
       *  Volume Conserving Model, equation (17) from Jennings et
       *  al.(2013)
       */
      double size_function (const double RV, const double redshift, const std::string model, const double b_eff, double slope=0.854, double offset=0.420, const double deltav_NL=-0.795, const double del_c=1.69, const std::string method_Pk="EisensteinHu", const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;

      /**
       *  @brief the void size function
       *
       *  @author Sofia Contarini
       *  @author sofia.contarini3@unibo.it
       *
       *  @param RV vector of radii
       *
       *  @param redshift the redshift
       *
       *  @param model size function model name; valid choices for
       *  model name are SvdW (Sheth and van de Weygaert, 2004),
       *  linear and Vdn (Jennings et al., 2013)
       *
       *  @param b_eff the effective bias of the sample
       *
       *  @param slope first coefficent to convert the effective bias
       *  (default value set to \f$0.854\f$)
       *
       *  @param offset second coefficent to convert the effective
       *  bias (default value set to \f$0.420\f$)
       *
       *  @param deltav_NL the non linear density contrast:
       *  \f$\rho_v/\rho_m\f$ (default value set to \f$-0.795\f$)
       *
       *  @param del_c critical value of the linear density field
       *  (default value set to \f$1.06\f$)
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
       *  @return the number density of voids as a function of radius.
       *  Volume Conserving Model, equation (17) from Jennings et
       *  al.(2013) for each radius
       */
      std::vector<double> size_function (const std::vector<double> RV, const double redshift, const std::string model, const double b_eff, double slope=0.854, double offset=0.420, const double deltav_NL=-0.795, const double del_c=1.69, const std::string method_Pk="EisensteinHu", const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;

      /**
       *  @brief number of voids computed from the void size function
       *  model for bins of radii spaced in log scale and for a
       *  specified survey/simulation volume
       *
       *  @author Sofia Contarini
       *  @author sofia.contarini3@unibo.it
       *
       *  @param min_r the minimum void radius
       *
       *  @param max_r the maximum void radius
       *
       *  @param num_bins number of bins of void radius
       *
       *  @param mean_z the mean redshift of the sample
       *
       *  @param Volume the volume of the survey/simulation in units
       *  of \f$(Mpc/h)^3\f$
       *
       *  @param model size function model name; valid choices for
       *  model name are SvdW (Sheth and van de Weygaert, 2004),
       *  linear and Vdn (Jennings et al., 2013)
       *
       *  @param b_eff the effective bias of the sample
       *
       *  @param slope first coefficent to convert the effective bias
       *  (default value set to \f$0.854\f$)
       *
       *  @param offset second coefficent to convert the effective
       *  bias (default value set to \f$0.420\f$)
       *
       *  @param deltav_NL the non linear density contrast:
       *  \f$\rho_v/\rho_m\f$ (default value set to \f$-0.795\f$)
       *
       *  @param del_c critical value of the linear density field
       *  (default value set to \f$1.06\f$)
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
       *  @return a vector of vectors with dimensions 2xnum_bins. The
       *  first component contains the void radii (centres of the bins
       *  computed in log scale), while the second the predicted void
       *  counts
       */
      std::vector<std::vector<double>> Nvoids (const double min_r, const double max_r, const int num_bins, const double mean_z, const double Volume, const std::string model, const double b_eff, double slope=0.854, double offset=0.420, const double deltav_NL=-0.795, const double del_c=1.69, const std::string method_Pk="EisensteinHu", const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;

      /**
       *  @brief number of voids computed from the void size function
       *  model for bins of radii spaced in log scale. The considered
       *  volume is computed from a (fraction of) a sphere shell,
       *  included in between two different redshifts
       *
       *  @author Sofia Contarini
       *  @author sofia.contarini3@unibo.it
       *
       *  @param min_r the minimum void radius
       *
       *  @param max_r the maximum void radius
       *
       *  @param num_bins number of bins of void radius
       *
       *  @param min_z the minimum redshift of the shell
       *
       *  @param max_z the maximum redshift of the shell
       *
       *  @param mean_z the mean redshift of the sample
       *
       *  @param Area sky area in units of squares degrees
       *
       *  @param model size function model name; valid choices for
       *  model name are SvdW (Sheth and van de Weygaert, 2004),
       *  linear and Vdn (Jennings et al., 2013)
       *
       *  @param b_eff the effective bias of the sample
       *
       *  @param slope first coefficent to convert the effective bias
       *  (default value set to \f$0.854\f$)
       *
       *  @param offset second coefficent to convert the effective
       *  bias (default value set to \f$0.420\f$)
       *
       *  @param deltav_NL the non linear density contrast:
       *  \f$\rho_v/\rho_m\f$ (default value set to \f$-0.795\f$)
       *
       *  @param del_c critical value of the linear density field
       *  (default value set to \f$1.06\f$)
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
       *  @return a vector of vectors with dimensions 2xnum_bins. The
       *  first component contains the void radii (centres of the bins
       *  computed in log scale), while the second the predicted void
       *  counts
       */
      std::vector<std::vector<double>> Nvoids (const double min_r, const double max_r, const int num_bins, const double min_z, const double max_z, const double mean_z, const double Area, const std::string model, const double b_eff, double slope=0.854, double offset=0.420, const double deltav_NL=-0.795, const double del_c=1.69, const std::string method_Pk="EisensteinHu", const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;

      /**
       *  @brief number of voids computed from the void size function
       *  model for bins of radii spaced in log scale. The considered
       *  volume is computed from a (fraction of) a sphere shell,
       *  included in between two different redshifts
       *
       *  @author Sofia Contarini
       *  @author sofia.contarini3@unibo.it
       *
       *  @param min_r the minimum void radius
       *
       *  @param max_r the maximum void radius
       *
       *  @param num_bins number of bins of void radius
       *
       *  @param min_z the minimum redshift of the shell
       *
       *  @param max_z the maximum redshift of the shell
       *
       *  @param mean_z the mean redshift of the sample
       *
       *  @param Area sky area in units of squares degrees
       *
       *  @param Sigma2R a 2D vector containing a vector of radii
       *  \f$R\f$ and the corresponding vector \f$\sigma^2(R)\f$
       * 
       *  @param model size function model name; valid choices for
       *  model name are SvdW (Sheth and van de Weygaert, 2004),
       *  linear and Vdn (Jennings et al., 2013)
       *
       *  @param b_eff the effective bias of the sample
       *
       *  @param slope first coefficent to convert the effective bias
       *  (default value set to \f$0.854\f$)
       *
       *  @param offset second coefficent to convert the effective
       *  bias (default value set to \f$0.420\f$)
       *
       *  @param deltav_NL the non linear density contrast:
       *  \f$\rho_v/\rho_m\f$ (default value set to \f$-0.795\f$)
       *
       *  @param del_c critical value of the linear density field
       *  (default value set to \f$1.06\f$)
       *
       *  @return a vector of vectors with dimensions 2xnum_bins. The
       *  first component contains the void radii (centres of the bins
       *  computed in log scale), while the second the predicted void
       *  counts
       */
      std::vector<std::vector<double>> Nvoids (const double min_r, const double max_r, const int num_bins, const double min_z, const double max_z, const double mean_z, const double Area, const std::vector<std::vector<double>> Sigma2R, const std::string model, const double b_eff, double slope=0.854, double offset=0.420, const double deltav_NL=-0.795, const double del_c=1.69) const;

      /**
       *  @brief the void size function
       *
       *  @author Tommaso Ronconi
       *  @author tommaso.ronconi@studio.unibo.it
       *
       *  @param RV radius
       *
       *  @param redshift the redshift
       *
       *  @param model_mf author(s) who proposed the mass function;
       *  valid authors are: PS (Press & Schechter), ST (Sheth &
       *  Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       *  al. 2006), Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH
       *  (halo MF by Shen et al. 2006), ShenF (filaments MF by Shen
       *  et al. 2006), ShenS (sheets MF by Shen et al. 2006), Tinker
       *  (Tinker et al. 2008), Crocce (Crocce et al. 2010),
       *  Angulo_FOF (FOF MF by Angulo et al. 2012), Angulo_Sub
       *  (SUBFIND MF by Angulo et al. 2012), Watson_FOF(FOF MF by
       *  Watson et al. 2012), Watson_SOH (MF for Spherical Overdensity
       *  Haloes by Watson et al. 2012), Manera (Manera et al. 2010),
       *  Bhattacharya (Bhattacharya et al. 2011), Courtin (Courtin 
       *  et al. 2010), Peacock (by Peacock at al. 2007)
       *
       *  @param del_v linear density contrast defining a void
       *
       *  @param model_sf size function model name; valid choices for
       *  model name are SvdW (Sheth and van de Weygaert, 2004),
       *  linear and Vdn (Jennings et al., 2013)
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
       *  @param is_parameter_file true \f$\rightarrow\f$ the
       *  input_file is a parameter file, used to compute the power
       *  spectrum with the method specified by method_Pk; false
       *  \f$\rightarrow\f$ the input_file is a file containing the
       *  power spectrum
       *  
       *  @return the number density of voids as a function of radius.
       *  Volume Conserving Model, equation (17) from Jennings et
       *  al.(2013)
       */
      double size_function (const double RV, const double redshift, const std::string model_mf, const double del_v, const std::string model_sf, const std::string method_Pk="EisensteinHu", const bool store_output=true, const std::string output_root="test", const double Delta=200., const std::string interpType="Linear", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
      
      ///@}

    };
  }
}

#endif
