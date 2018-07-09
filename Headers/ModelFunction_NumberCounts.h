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
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
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

	/// false &rarr; don't use selection function; true &rarr; use selection function
        bool use_SF;
	
	/// function to interpolate the selection function in mass and redshift
	std::shared_ptr<glob::FuncGrid2D> interp_SelectionFunction;

	/**
	 *  @brief default constructor
	 *  @return object of type STR_data_model
	 */
	STR_NC_data_model () = default;
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
      double mass_function (const double mass, cosmology::Cosmology cosmology, const double redshift, const std::string model_MF, const double Delta, const bool isDelta_vir, const cbl::glob::FuncGrid interp_Pk, const double kmax);

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
      std::vector<double> mass_function (const std::vector<double> mass, cosmology::Cosmology cosmology, const double redshift, const std::string model_MF, const double Delta, const bool isDelta_vir, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax);

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
      std::vector<std::vector<double>> mass_function (const std::vector<double> redshift, const std::vector<double> mass, cosmology::Cosmology cosmology, const std::string model_MF, const double Delta, const bool isDelta_vir, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax);

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
       * @param Delta \f$\Delta\f$: the overdensity, defined as the
       * mean interior density relative to the background
       *
       * @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
       * virial overdensity
       *
       * @param interp_sigmaM interpolating function of \f$ \sigma(M)\f$
       *
       * @param interp_DlnsigmaM interpolating function of \f$ \mathrm{d} \ln(\sigma(M)) / \mathrm{d} M \f$
       *
       * @param npt_redshift number of redshifts bins to compute the integral
       *
       * @param npt_mass number of mass bins to compute the integral
       *
       * @return values of the mass function as a function of redshift
       * and mass
       */
      double number_counts(const double redshift_min, const double redshift_max, const double Mass_min, const double Mass_max, cosmology::Cosmology cosmology, const double Area, const std::string model_MF, const double Delta, const bool isDelta_vir, const glob::FuncGrid interp_sigmaM, const  glob::FuncGrid interp_DlnsigmaM, const int npt_redshift=10, const int npt_mass=10);

    }
  }
}

#endif
