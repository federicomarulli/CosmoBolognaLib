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
 *  @file Headers/ModelFunction_NumberCounts.h
 *
 *  @brief Global functions to model number counts
 *  of any type
 *
 *  This file contains all the prototypes of the global functions used
 *  to model number counts of any type
 *  
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODFUNCNC__
#define __MODFUNCNC__

#include "Cosmology.h"


// ============================================================================


namespace cbl {

  namespace modelling {

    namespace numbercounts {
      
      /**
       *  @struct STR_NC_data_model
       *  @brief the structure STR_NC_data_model
       *
       *  This structure contains the data used for statistical
       *  analyses of the numbercounts
       */
      struct STR_NC_data_model {

	/// false &rarr; data not from a simulation snapshot; true &rarr; data from a simulation snapshot 
	bool isSnapshot;

	/// fiducial cosmology
	std::shared_ptr<cosmology::Cosmology> cosmology;

	/// cosmological parameters
	std::vector<cosmology::CosmologicalParameter> Cpar;

	/// redshift
	double redshift;

	/// method to compute the dark matter power spectrum
	std::string method_Pk;
      
	/// false &rarr; linear power spectrum; true &rarr; non-linear power spectrum
	bool NL;

	/// minimum wave vector module up to which the power spectrum is computed
	double k_min;

	/// maximum wave vector module up to which the power spectrum is computed
	double k_max;

	/// number of steps used to compute the binned dark matter correlation function
	int step;

	/// vector of wave vector modules
	std::vector<double> kk;

	/// the output_dir directory where the output of external codes are written
	std::string output_dir;

	/// true \f$\rightarrow\f$ the output files created by the Boltzmann solver are stored; false \f$\rightarrow\f$ the output files are removed
	bool store_output;
      
	/// output root of the parameter file used to compute the dark matter power spectrum
	std::string output_root;

	/// 0 &rarr; don't normalize the power spectrum; 1 &rarr; normalize the power spectrum
	int norm;

	/// name of the parameter file
	std::string file_par;

	/// accuracy of the GSL integration
	double prec;

	/// &Delta;: the overdensity, defined as the mean interior density relative to the background
	double Delta;

	/// isDelta_Vir
	bool isDelta_Vir;

	/// author(s) who proposed the mass function
	std::string model_MF;

        /// minimum redshift
	double z_min; 

	/// maximum redshift
	double z_max;

	/// number of redshift steps used to compute the binned mass function
	int z_step;

	/// vector of redshifts
	std::vector<double> z_vector;

	/// minimum mass 
	double Mass_min; 

	/// maximum mass
	double Mass_max; 

	/// number of mass steps used to compute the binned mass function
	int Mass_step;

	/// vector of masses
	std::vector<double> Mass_vector;

	/// the survey aperture Area
	double area_rad;

	/// the survey volume
	double Volume;

	/// false &rarr; don't use the selection function; true &rarr; use the selection function
        bool use_SF;

	/// true &rarr; sigma8 is a free parameter; false &rarr; sigma8 can be considered a derived parameter
	bool is_sigma8_free;
	
	/// function to interpolate the selection function in mass and redshift
	std::shared_ptr<glob::FuncGrid2D> interp_SelectionFunction;

	/**
	 *  @brief default constructor
	 *  @return object of type STR_data_model
	 */
	STR_NC_data_model () = default;
      };


   
      /**
       *  @struct STR_NCSF_data_model
       *  @brief the structure STR_NCSF_data_model
       *
       *  This structure contains the data used for statistical
       *  analyses of the void numbercounts
       */
      struct STR_NCSF_data_model {

	/// vector of radii
	std::vector<double> radii;

	/// fiducial cosmology
	std::shared_ptr<cosmology::Cosmology> cosmology;

	/// cosmological parameters
	std::vector<cosmology::CosmologicalParameter> Cpar;

	/// redshift
	double redshift;

	/// the size function model
	std::string model_SF;

	/// the effective bias of the sample
	double b_eff;

	/// first coefficent to convert the effective bias
	double b_slope;

	/// second coefficent to convert the effective bias
	double b_offset;

	/// the non linear density contrast:
	double deltav_NL;

	/// critical value of the linear density field
	double delta_c;
	
	/// method to compute the dark matter power spectrum
	std::string method_Pk;

	/// true \f$\rightarrow\f$ the output files created by the Boltzmann solver are stored; false \f$\rightarrow\f$ the output files are removedt
	bool store_output;

	/// output root of the parameter file used to compute the dark matter power spectrum
	std::string output_root;

	/// interpType method to interpolate the power spectrum
	std::string interpType;

	/// maximum wave vector module up to which the power spectrum is computed
	double k_max;

	/// either the parameter file or the power spectrum file
	std::string input_file;

	/// if the input_file is a file containing the power spectrum
	/// or a parameter file used to compute it
	bool is_parameter_file;
      
	/**
	 *  @brief default constructor
	 *  @return object of type STR_data_model
	 */
	STR_NCSF_data_model () = default;
      };
      

      /**
       * @brief the filter to compute \f$\sigma(R)\f$
       *
       * @param kk the wave vector module
       *
       * @param radius the radius
       *
       * @return the value of the top-hat filter
       */
      double Filter_sigmaR (const double kk, const double radius); 

      /**
       * @brief the filter to compute \f$\mathrm{d} \sigma(R) / \mathrm{d} R\f$
       *
       * @param kk the wave vector module
       *
       * @param radius radius of the filter
       *
       * @return the value of the filter to compute 
       * \f$\mathrm{d} \sigma(R) / \mathrm{d} R\f$
       */
      double Filter_dsigmaR (const double kk, const double radius); 

      /**
       * @brief compute \f$ \sigma(M), \mathrm{d} \ln(\sigma(M)) / \mathrm{d} M \f$
       *
       * @param sigmaM vector containing the values of \f$\sigma(M)\f$
       *
       * @param dlnsigmaM vector containing the values of \f$\mathrm{d} \ln(\sigma(M)) / \mathrm{d} M\f$
       *
       * @param mass the mass value
       *
       * @param interp_Pk function to interpolate the power spectrum
       *
       * @param kmax the maximum value of the wavevector module
       *
       * @param rho the mean density of the Universe
       *
       * @return none
       */
      void sigmaM_dlnsigmaM (double &sigmaM, double &dlnsigmaM, const double mass, const cbl::glob::FuncGrid interp_Pk, const double kmax, const double rho);

      /**
       * @brief compute \f$ \sigma(M), \mathrm{d} \ln(\sigma(M)) / \mathrm{d} M \f$
       *
       * @param sigmaM vector containing the values of \f$\sigma(M)\f$
       *
       * @param dlnsigmaM vector containing the values of \f$\mathrm{d} \ln(\sigma(M)) / \mathrm{d} M\f$
       *
       * @param mass vector of mass used for computation
       *
       * @param kk vector of wavevector modules
       *
       * @param Pk vector containing the dark matter power spectrum
       *
       * @param interpType the kind of interpolation
       *
       * @param kmax the maximum value of the wavevector module
       *
       * @param rho the mean density of the Universe
       *
       * @return none
       */
      void sigmaM_dlnsigmaM (std::vector<double> &sigmaM, std::vector<double> &dlnsigmaM, const std::vector<double> mass, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax, const double rho);

      /**
       * @brief compute \f$ \sigma(M), \mathrm{d} \ln(\sigma(M)) / \mathrm{d} M \f$
       * and return them as interpolating function
       *
       * @param mass vector of mass used for computation
       *
       * @param cosmology the cosmology 
       *
       * @param kk vector of wavevector modules
       *
       * @param Pk vector containing the dark matter power spectrum
       *
       * @param interpType the kind of interpolation
       *
       * @param kmax the maximum value of the wavevector module
       *
       * @return vector containing the interpolation function for
       * \f$ \sigma(M), \mathrm{d} \ln(\sigma(M)) / \mathrm{d} M \f$
       */
      std::vector<cbl::glob::FuncGrid> sigmaM_dlnsigmaM (const std::vector<double> mass, cosmology::Cosmology cosmology, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax);

      /**
       * @brief compute the mass function
       *
       * @param mass the mass
       *
       * @param cosmology the cosmology 
       *
       * @param redshift the redshift
       *
       * @param model_MF author(s) who proposed the mass function;
       * valid authors are: PS (Press & Schechter), ST (Sheth &
       * Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       * al. 2006), Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH
       * (halo MF by Shen et al. 2006), ShenF (filament MF by Shen et
       * al. 2006), ShenS (sheet MF by Shen et al. 2006), Tinker
       * (Tinker et al. 2008), Crocce (Crocce et al. 2010),
       * Angulo_FOF (FoF MF by Angulo et al. 2012), Angulo_Sub
       * (SUBFIND MF by Angulo et al. 2012), Watson_FOF (FoF MF by
       * Watson et al. 2012), Watson_SOH (Spherical Overdensity halo
       * MF by Watson et al. 2012), Manera (Manera et al. 2010),
       * Bhattacharya (Bhattacharya et al. 2011), Courtin (Courtin et
       * al. 2010), Peacock (by Peacock at al. 2007)
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       * @param Delta \f$\Delta\f$: the overdensity, defined as the
       * mean interior density relative to the background
       *
       * @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
       * virial overdensity
       *
       * @param interp_Pk function to interpolate the power spectrum
       *
       * @param kmax the maximum value of the wavevector module
       *
       * @return value of the mass function
       */
      double mass_function (const double mass, cosmology::Cosmology cosmology, const double redshift, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_vir, const cbl::glob::FuncGrid interp_Pk, const double kmax);

      /**
       * @brief compute the mass function
       *
       * @param mass vector of masses
       *
       * @param cosmology the cosmology 
       *
       * @param redshift the redshift
       *
       * @param model_MF author(s) who proposed the mass function;
       * valid authors are: PS (Press & Schechter), ST (Sheth &
       * Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       * al. 2006), Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH
       * (halo MF by Shen et al. 2006), ShenF (filament MF by Shen et
       * al. 2006), ShenS (sheet MF by Shen et al. 2006), Tinker
       * (Tinker et al. 2008), Crocce (Crocce et al. 2010),
       * Angulo_FOF (FoF MF by Angulo et al. 2012), Angulo_Sub
       * (SUBFIND MF by Angulo et al. 2012), Watson_FOF (FoF MF by
       * Watson et al. 2012), Watson_SOH (Spherical Overdensity halo
       * MF by Watson et al. 2012), Manera (Manera et al. 2010),
       * Bhattacharya (Bhattacharya et al. 2011), Courtin (Courtin et
       * al. 2010), Peacock (by Peacock at al. 2007)
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       * @param Delta \f$\Delta\f$: the overdensity, defined as the
       * mean interior density relative to the background
       *
       * @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
       * virial overdensity
       *
       * @param kk vector of wavevector modules
       *
       * @param Pk vector containing the dark matter power spectrum
       *
       * @param interpType the kind of interpolation
       *
       * @param kmax the maximum value of the wavevector module
       *
       * @return values of the mass function
       */
      std::vector<double> mass_function (const std::vector<double> mass, cosmology::Cosmology cosmology, const double redshift, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_vir, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax);

      /**
       * @brief compute the mass function as function
       * of mass and redshift
       *
       * @param redshift vector of redshift
       *
       * @param mass vector of masses
       *
       * @param cosmology the cosmology 
       *
       * @param model_MF author(s) who proposed the mass function;
       * valid authors are: PS (Press & Schechter), ST (Sheth &
       * Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       * al. 2006), Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH
       * (halo MF by Shen et al. 2006), ShenF (filament MF by Shen et
       * al. 2006), ShenS (sheet MF by Shen et al. 2006), Tinker
       * (Tinker et al. 2008), Crocce (Crocce et al. 2010),
       * Angulo_FOF (FoF MF by Angulo et al. 2012), Angulo_Sub
       * (SUBFIND MF by Angulo et al. 2012), Watson_FOF (FoF MF by
       * Watson et al. 2012), Watson_SOH (Spherical Overdensity halo
       * MF by Watson et al. 2012), Manera (Manera et al. 2010),
       * Bhattacharya (Bhattacharya et al. 2011), Courtin (Courtin et
       * al. 2010), Peacock (by Peacock at al. 2007)
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       * @param Delta \f$\Delta\f$: the overdensity, defined as the
       * mean interior density relative to the background
       *
       * @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
       * virial overdensity
       *
       * @param kk vector of wavevector modules
       *
       * @param Pk vector containing the dark matter power spectrum
       *
       * @param interpType the kind of interpolation
       *
       * @param kmax the maximum value of the wavevector module
       *
       * @return values of the mass function as a function of redshift
       * and mass
       */
      std::vector<std::vector<double>> mass_function (const std::vector<double> redshift, const std::vector<double> mass, cosmology::Cosmology cosmology, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_vir, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax);

      /**
       *  @brief the void size function
       *
       *  @param cosmology the cosmology 
       *
       *  @param radii the void radii
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
       *  name
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
      std::vector<double> size_function (cosmology::Cosmology cosmology, const std::vector<double> radii, const double redshift, const std::string model, const double b_eff, double slope=0.854, double offset=0.420, const double deltav_NL=-0.795, const double del_c=1.69, const std::string method_Pk="CAMB", const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      /**
       * @brief compute the number counts as function
       * of mass and redshift
       *
       * @param redshift_min minimum redshift
       *
       * @param redshift_max maximum redshift
       *
       * @param Mass_min minimum mass
       *
       * @param Mass_max maximum mass
       *
       * @param cosmology the cosmology 
       *
       * @param Area the area in degrees
       *
       * @param model_MF author(s) who proposed the mass function;
       * valid authors are: PS (Press & Schechter), ST (Sheth &
       * Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       * al. 2006), Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH
       * (halo MF by Shen et al. 2006), ShenF (filament MF by Shen et
       * al. 2006), ShenS (sheet MF by Shen et al. 2006), Tinker
       * (Tinker et al. 2008), Crocce (Crocce et al. 2010),
       * Angulo_FOF (FoF MF by Angulo et al. 2012), Angulo_Sub
       * (SUBFIND MF by Angulo et al. 2012), Watson_FOF (FoF MF by
       * Watson et al. 2012), Watson_SOH (Spherical Overdensity halo
       * MF by Watson et al. 2012), Manera (Manera et al. 2010),
       * Bhattacharya (Bhattacharya et al. 2011), Courtin (Courtin et
       * al. 2010), Peacock (by Peacock at al. 2007)
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       * @param Delta \f$\Delta\f$: the overdensity, defined as the
       * mean interior density relative to the background
       *
       * @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
       * virial overdensity
       *
       * @param interp_sigmaM interpolating function of \f$
       * \sigma(M)\f$
       *
       * @param interp_DlnsigmaM interpolating function of \f$
       * \mathrm{d} \ln(\sigma(M)) / \mathrm{d} M \f$
       *
       * @param npt_redshift number of redshifts bins to compute the
       * integral
       *
       * @param npt_mass number of mass bins to compute the integral
       *
       * @return values of the mass function as a function of redshift
       * and mass
       */
      double number_counts(const double redshift_min, const double redshift_max, const double Mass_min, const double Mass_max, cosmology::Cosmology cosmology, const double Area, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_vir, const glob::FuncGrid interp_sigmaM, const  glob::FuncGrid interp_DlnsigmaM, const int npt_redshift=10, const int npt_mass=10);

    }
  }
}

#endif
