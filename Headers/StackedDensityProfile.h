/********************************************************************
 *  Copyright (C) 2020 by Giorgio Lesci and Federico Marulli        *
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
 *  @file Headers/StackedDensityProfile.h
 *
 *  @brief The class StackedDensityProfile
 *
 *  This file defines the interface of the class StackedDensityProfile,
 *  used to measure the stacked density profile of galaxy clusters
 *  from lensing data
 *
 *  @author Giorgio Lesci (and Fabio Bellagamba, Federico Marulli)
 *
 *  @author giorgio.lesci2@unibo.it (and fabiobg83@gmail.com, federico.marulli3@unibo.it)
 */

#ifndef __STACKPROFILE__
#define __STACKPROFILE__

#include "Catalogue.h"
#include "Cosmology.h"
#include "Measure.h"


// ===================================================================================================


namespace cbl {

  namespace measure {

    /**
     *  @brief The namespace of the <B> stacked density profile of clusters
     *  </B>
     *  
     *  The \e measure::stackprofile namespace contains all the functions and
     *  classes to measure the stacked density profile of galaxy clusters
     */
    namespace stackprofile {

      /**
       *  @class StackedDensityProfile StackedDensityProfile.h
       *  "Headers/StackedDensityProfile.h"
       *
       *  @brief The class StackedDensityProfile
       *
       *  Class used to measure the stacked surface density profiles of galaxy clusters,
       *  i.e. \f$\Delta\Sigma(r)\f$ [\f$h\f$ M\f$_\odot\f$/pc\f$^2\f$].
       *  Cosmological units are forced.
       *
       */
      class StackedDensityProfile : public Measure {

      protected :	  

	/**
	 *  @name Protected member functions and members to measure the cluster stacked density profile
	 *  and the relative covariance matrix
	 */
	///@{

	
	/**
	 *  @brief set the multiplicative shear calibration, m
	 *
	 */
	void m_set_multiplicative_calib();
	
	/**
	 *  @brief colour selection function in the colour-colour
	 *  space (\f$r-i\f$) vs (\f$g-r\f$)
	 *
	 *  @param i_gal galaxy index
	 *
	 *  @param clu_zbin index of the redshift bin of the lensing cluster
	 *
	 *  @return true if the galaxy satisfies the conditions,
	 *  otherwise it returns false
	 *
	 */
	bool m_colourSelection_gri (const int i_gal, const int clu_zbin);
	
	/**
	 *  @brief redshift selection function
	 *
	 *  @param i_clu cluster index
	 *
	 *  @param i_gal galaxy index
	 *
	 *  @param clu_zbin index of the redshift bin of the lensing cluster
	 *
	 *  @return true if the galaxy satisfies the conditions,
	 *  otherwise it returns false
	 *
	 */
	bool m_photzSelection (const int i_clu, const int i_gal, const int clu_zbin);
	
	/**
	 *  @brief check if all the necessary variables in the
	 *  galaxy and cluster catalogues are set
	 *
	 */
	void m_check_catalogue_variables();
	
	/**
	 *  @brief resize the private member arrays
	 *
	 *  @param rad_min minimum distance from the
	 *  cluster centre considered
	 *
	 *  @param rad_max maximum distance from the
	 *  cluster centre considered
	 *
	 *  @param nRad number of cluster radial bins
	 *
	 *  @param log_rad if true radial bins logarithmically
	 *  spaced
	 *
	 */
	void m_resize(const double rad_min, const double rad_max, const int nRad, const bool log_rad);
	
	/**
	 *  @brief linked galaxy list
	 *
	 */
	void m_linked_list();
	
	/**
	 *  @brief add galaxies to profiles
	 *
	 *  @param i_gal galaxy index
	 *
	 *  @param clu_index cluster index
	 *
	 *  @param coscludec cos of cluster Dec
	 *
	 *  @param clu_dist cluster angular diameter distance
	 *
	 */
	void m_add_galaxy(const int i_gal, const int clu_index, const double coscludec, const double clu_dist);
	
	/**
	 *  @brief compute a single cluster profile
	 *
	 *  @param clu_index cluster index
	 * 
	 */
	void m_profile(const int clu_index);
	
	/**
	 *  @brief perform the stacking
	 * 
	 */
	void m_stacker();

	/**
	 *  @brief bootstrap resampling to retrieve
	 *  the covariance matrix of the stacked signal,
	 *  providing the stacked profile measurement in output
	 *
	 *  @param z_proxy_bin vector containing the redshift
	 *  and proxy bins to be stored, respectively
	 * 
	 *  @return a shared pointer to Data
	 *
	 */
	std::shared_ptr<data::Data> m_make_bootstrap(const std::vector<int> z_proxy_bin);
	
	/**
	 *  @brief check if the results obtained from the stacking
	 *  have already been written to a file
	 *
	 *  @param checked_file the file to be checked
	 *
	 *  @param z_proxy_bin vector containing the redshift
	 *  and proxy bins to be stored, respectively
	 *
	 *  @return true if the stacking file already exists
	 *
	 */
	bool m_check_file(const std::string checked_file, const std::vector<int> z_proxy_bin);
	
	/**
	 *  @brief write the results obtained from the stacking on file
	 *
	 *  @param output_dir output directory
	 *
	 *  @param output_file output file
	 *
	 */
	void m_write(const std::string output_dir, const std::string output_file);
	
	/// input cosmology
	std::shared_ptr<cosmology::Cosmology> m_cosm;
	
	/// input galaxy catalogue
	std::shared_ptr<catalogue::Catalogue> m_galData;
	
	/// input cluster catalogue
	std::shared_ptr<catalogue::Catalogue> m_cluData;
	
	/// vector of pointers to the colour selection function s
	std::vector<bool (StackedDensityProfile::*)(const int, const int)> m_colourSel; // std::vector<std::function<bool(const int)>> m_colourSel;
	
	/// vector of pointers to the photo-z selection functions
	std::vector<bool (StackedDensityProfile::*)(const int, const int, const int)> m_photzSel;
	
	/// logic operator between colour and photo-z selection
	std::vector<std::function<bool(const std::vector<bool>)>> m_logicSel;
	
	/// vector for checking if the colour selection is set in all the z bins
	std::vector<bool> m_isSet_colourSel;
	
	/// vector for checking if the photo-z selection is set in all the z bins
	std::vector<bool> m_isSet_photzSel;
	
	/// vector for checking if the logic selection linking colour and z selections is set in all the z bins
	std::vector<bool> m_isSet_logicSel;
	
	/// vector of indices of the background galaxies selected through the redshift selection
	std::vector<std::vector<int>> m_background_idx_z;
	
	/// vector of indices of the background galaxies selected through the colour selection
	std::vector<std::vector<int>> m_background_idx_colour;
	
	/// bool stating if the measure is read from file
	bool m_measure_is_read;
	
	/// minimum interval between the cluster and the source redshifts
	double m_delta_redshift;
	
	/// minimum signal-to-noise
	double m_SN_min;
	
	/// pixel size (in deg)
	double m_pix_size;
	
	/// the multiplicative calibration bias mean and standard deviation in the z bins
	std::vector<std::vector<double>> m_m_calib;
	
	/// the redshift edges within which the multiplicative calibration bias is evaluated
	std::vector<std::vector<double>> m_m_calib_zEdges;
	
	/// number of regions for the resampling
	double m_n_resampling;
	
	/// vector of the inputs stored as strings
	std::vector<std::string> m_inputs_to_str;
	
	/// slope for the observable weighted mean
	double m_rad_alpha;
	
	/// gamma slope for lensing-weighted observables
	double m_obs_gamma;
	
	/// the parameters for the colour selection
	std::vector<std::vector<double>> m_colour_sel_pars;
	
	/// the parameters for the redshift selection
	std::vector<std::vector<double>> m_zphot_sel_pars;
	
	/// the logic selection linking redshift and colour selections
	std::vector<std::vector<std::string>> m_logic_sel_par;
	
	/// radius bin edges
	std::vector<double> m_rad_arr;
	
	/// phot-z bin edges
	std::vector<double> m_z_binEdges;
	
	/// mass proxy bin edges
	std::vector<std::vector<double>> m_proxy_binEdges;
	
	/// a "first" for each pixel
	std::vector<std::vector<int>> m_first;
	
	/// a "last" for each pixel
	std::vector<std::vector<int>> m_last;
	
	/// a "next" for each galaxy
	std::vector<int> m_next;
	
	/// number of clusters in a bin of z and mass proxy
	std::vector<std::vector<double>> m_nClu_inBin;
	
	/// number of galaxies for annulus
	std::vector<int> m_ngal_arr;
	
	/// number of galaxies for annulus for cluster
	std::vector<std::vector<std::vector<std::vector<int>>>> m_single_ngal_arr;
	
	/// stacked number of galaxies for annulus
	std::vector<std::vector<std::vector<int>>> m_stacked_ngal_arr;
	
	/// tangential surface density
	std::vector<double> m_deltasigma_t;
	
	/// tangential surface density, stored for each single cluster
	std::vector<std::vector<std::vector<std::vector<double>>>> m_single_deltasigma_t;
	
	/// tangential stacked surface density
	std::vector<std::vector<std::vector<double>>> m_stacked_deltasigma_t;
	
	/// tangential stacked surface density error
	std::vector<std::vector<std::vector<double>>> m_stacked_deltasigma_t_err;
	
	/// cross surface density
	std::vector<double> m_deltasigma_x;
	
	/// cross surface density, stored for each single cluster
	std::vector<std::vector<std::vector<std::vector<double>>>> m_single_deltasigma_x;
	
	/// cross stacked surface density
	std::vector<std::vector<std::vector<double>>> m_stacked_deltasigma_x;
	
	/// cross stacked surface density error
	std::vector<std::vector<std::vector<double>>> m_stacked_deltasigma_x_err;
	
	/// error on the surface density
	std::vector<double> m_deltasigma_err;
	
	/// error on the surface density, stored for each single cluster
	std::vector<std::vector<std::vector<std::vector<double>>>> m_single_deltasigma_err;
	
	/// error on the stacked surface density
	std::vector<std::vector<std::vector<double>>> m_stacked_deltasigma_err;
	
	/// wetasquareSum
	std::vector<double> m_wetasquareSum;
	
	/// stacked wetasquareSum
	std::vector<std::vector<std::vector<double>>> m_stacked_wetasquareSum;
	
	/// wSum
	std::vector<double> m_wSum;
	
	/// stacked wSum
	std::vector<std::vector<std::vector<double>>> m_stacked_wSum;
	
	/// wetaSum
	std::vector<double> m_wetaSum;
	
	/// stacked wetaSum
	std::vector<std::vector<std::vector<double>>> m_stacked_wetaSum;
	
	/// deltasigmaSum
	std::vector<double> m_deltasigmaSum;
	
	/// stacked deltasigmaSum
	std::vector<std::vector<std::vector<double>>> m_stacked_deltasigmaSum;
	
	/// stacked deltasigma_wei
	std::vector<std::vector<std::vector<double>>> m_stacked_deltasigma_wei;
	
	/// effective radii
	std::vector<double> m_rad_eff_arr; 
	
	/// stacked effective radii
	std::vector<std::vector<std::vector<double>>> m_stacked_rad_eff_arr;
	
	/// effective radii errors
	std::vector<double> m_rad_sigma_arr;
	
	/// stacked effective radii errors
	std::vector<std::vector<std::vector<double>>> m_stacked_rad_sigma_arr;
	
	/// constant factor for lensing quantities (c*c/(4piG) in M_sun/Mpc)
	double m_sigma_fac;
	
	/// starting RA
	double m_ra_start;
	
	/// starting Dec
	double m_dec_start;
	
	/// number of pixels in RA and Dec
	std::vector<int> m_nPix;
	
	/// effective mass proxy linked to the stacked signal
	std::vector<std::vector<double>> m_proxy_eff;
	
	/// error on the effective mass proxy linked to the stacked signal
	std::vector<std::vector<double>> m_proxy_sigma;

	/// effective redshift linked to the stacked signal
	std::vector<std::vector<double>> m_z_eff;
	
	/// error on the effective redshift linked to the stacked signal
	std::vector<std::vector<double>> m_z_sigma;
	
	/// Bootstrap covariance matrix for the stacked signal
	std::vector<std::vector<std::vector<std::vector<double>>>> m_deltasigma_cov_matr;
	
	/// 
	
	/// 

	///@}

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *
	 *  
	 */
	StackedDensityProfile () = default;

	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~StackedDensityProfile () = default;


	/**
	   *  @brief constructor performing the stacking 
	   *  by reading a .fits galaxy file and a cluster file
	   *
	   *  @param cosm cosmological model
	   *
	   *  @param gal_cat catalogue of galaxies
	   *
	   *  @param clu_cat catalogue of galaxy clusters
	   *
	   *  @param delta_redshift minimum interval between 
	   *  the cluster and the source redshifts. That is,
	   *  it is \f$\Delta_z\f$ in the equation 
	   *  \f$z_g>z_c+\Delta_z\f$, where \f$z_g\f$ and 
	   *  \f$z_c\f$ are the galaxy and cluster 
	   *  mean redshifts, respectively.
	   *  This condition must be satisfied along with
	   *  \f$p(z)\f$ or/and color conditions.
	   *
	   *  @param z_binEdges redshift bin edges
	   *
	   *  @param proxy_binEdges proxy bin edges
	   *
	   *  @param rad_min minimum distance from the
	   *  cluster centre considered (in Mpc/h)
	   *
	   *  @param rad_max maximum distance from the
	   *  cluster centre considered (in Mpc/h)
	   *
	   *  @param nRad number of cluster radial bins
	   *
	   *  @param log_rad if true, the radial bins are logarithmically
	   *  spaced. Otherwise, a linear binning is used
	   *
	   *  @param SN_min minimum signal-to-noise considered 
	   *  for the clusters
	   *
	   *  @param pix_size pixel size in deg
	   *
	   *  @param multiplicative_calibration_stats a vector of vectors containing
	   *  mean and standard deviation of the multiplicative shear calibration
	   *  parameter, usually denoted as m, in all the redshift bins where it is evaluated. 
	   *  If not provided, the galaxy-by-galaxy values of m are used
	   *
	   *  @param multiplicative_calibration_zEdges vector of vectors containing
	   *  the lower and upper edge of the redshift bins for each estimate of the 
	   *  multiplicative shear calibration parameter, m
	   *
	   *  @param rad_alpha slope for the observable weighted mean
	   *
	   *  @param obs_gamma gamma slope for lensing-weighted observables
	   *
	   *  Cluster stacked density profile
	   */
	StackedDensityProfile (cosmology::Cosmology cosm, const catalogue::Catalogue gal_cat, const catalogue::Catalogue clu_cat, const double delta_redshift, std::vector<double> z_binEdges, std::vector<std::vector<double>> proxy_binEdges, const double rad_min, const double rad_max, const int nRad, const bool log_rad, const double SN_min, const double pix_size, const std::vector<std::vector<double>> multiplicative_calibration_stats={}, const std::vector<std::vector<double>> multiplicative_calibration_zEdges={}, const double rad_alpha=1., const double obs_gamma=1.);
	///@}
	
	/**
	 *  @name Member functions to measure the stacked profile
	 */
	///@{
	
	/**
	 *  @brief Set the colour selection in the i-th cluster redshift bin.
	 *
	 *  @param colour_space the colour space; possibilities are: 
	 *  "ri_gr" (\f$r-i\f$ vs \f$g-r\f$)
	 *
	 *  @param z_bin index of the cluster redshift bin
	 *
	 *  @param C1 for example, for the colour space \f$(r-i)\f$ vs \f$(g-r)\f$,
	 *  it is \f$C_1\f$ in \f$ (g-r)<C_1 \f$
	 *
	 *  @param C2 for example, for the colour space \f$(r-i)\f$ vs \f$(g-r)\f$,
	 *  it is \f$C_2\f$ in \f$ (r-i)>C_2 \f$
	 *
	 *  @param C3 for example, for the colour space \f$(r-i)\f$ vs \f$(g-r)\f$,
	 *  it is \f$C_3\f$ in \f$ (r-i)>C_3*(g-r)+C_4 \f$
	 *
	 *  @param C4 for example, for the colour space \f$(r-i)\f$ vs \f$(g-r)\f$,
	 *  it is \f$C_4\f$ in \f$ (r-i)>C_3*(g-r)+C_4 \f$
	 *
	 *  @param z_min minimum redshift of the selected galaxies
	 *
	 */
	void set_colour_selection(const std::string colour_space, const int z_bin, const double C1, const double C2, const double C3, const double C4, const double z_min);
	
	/**
	 *  @brief Set the colour selection in all cluster redshift bins.
	 *
	 *  @param colour_space the colour space; possibilities are: 
	 *  "ri_gr" (\f$r-i\f$ vs \f$g-r\f$)
	 *
	 *  @param C1 for example, for the colour space \f$(r-i)\f$ vs \f$(g-r)\f$,
	 *  it is \f$C_1\f$ in \f$ (g-r)<C_1 \f$
	 *
	 *  @param C2 for example, for the colour space \f$(r-i)\f$ vs \f$(g-r)\f$,
	 *  it is \f$C_2\f$ in \f$ (r-i)>C_2 \f$
	 *
	 *  @param C3 for example, for the colour space \f$(r-i)\f$ vs \f$(g-r)\f$,
	 *  it is \f$C_3\f$ in \f$ (r-i)>C_3*(g-r)+C_4 \f$
	 *
	 *  @param C4 for example, for the colour space \f$(r-i)\f$ vs \f$(g-r)\f$,
	 *  it is \f$C_4\f$ in \f$ (r-i)>C_3*(g-r)+C_4 \f$
	 *
	 *  @param z_min minimum redshift of the selected galaxies
	 *
	 */
	void set_colour_selection(const std::string colour_space, const double C1, const double C2, const double C3, const double C4, const double z_min);
	
	/**
	 *  @brief Set the redshift selection in the i-th cluster redshift bin
	 *
	 *  @param z_bin index of the cluster redshift bin
	 *
	 *  @param deltaz additive term, \f$\Delta z\f$, to
	 *  the cluster redshift, i.e. how much the minimum galaxy redshift,
	 *  \f$z_{\rm g,\,min}\f$, (in the tail of the relative posterior) must be higher 
	 *  than the cluster redshift, \f$z_c\f$, that is \f$z_{\rm g,\,min} > z_c+\Delta z\f$
	 *
	 *  @param zgal_min minimum value for the galaxy mean redshift
	 *
	 *  @param zgal_max maximum value for the galaxy mean redshift
	 *
	 *  @param ODDS_min minimum value for the ODDS parameter
	 *
	 */
	 
	void set_zphot_selection(const int z_bin, const double deltaz, const double zgal_min, const double zgal_max, const double ODDS_min);
	
	/**
	 *  @brief Set the redshift selection in all cluster redshift bins
	 *
	 *  @param deltaz additive term, \f$\Delta z\f$, to
	 *  the cluster redshift, i.e. how much the minimum galaxy redshift,
	 *  \f$z_{\rm g,\,min}\f$, (in the tail of the relative posterior) must be higher 
	 *  than the cluster redshift, \f$z_c\f$, that is \f$z_{\rm g,\,min} > z_c+\Delta z\f$
	 *
	 *  @param zgal_min minimum value for the galaxy mean redshift
	 *
	 *  @param zgal_max maximum value for the galaxy mean redshift
	 *
	 *  @param ODDS_min minimum value for the ODDS parameter
	 *
	 */
	 
	void set_zphot_selection(const double deltaz, const double zgal_min, const double zgal_max, const double ODDS_min);
	
	/**
	 *  @brief Set the logic operator linking redshift and colour
	 *  selections in the i-th cluster redshift bin
	 *
	 *  @param z_bin index of the cluster redshift bin
	 *
	 *  @param sel "and" or "or"
	 *
	 */
	 
	void set_logic_selection(const int z_bin, const std::string sel);
	
	/**
	 *  @brief Set the logic operator linking redshift and colour
	 *  selections in the all the cluster redshift bins
	 *
	 *  @param sel "and" or "or"
	 *
	 */
	 
	void set_logic_selection(const std::string sel);
	
	/**
	 *  @brief measure the stacked profiles in all the bins of redshift and 
	 *  mass proxy, providing in output the stacked profile in the redshift 
	 *  and proxy bins chosen through the parameter z_proxy_bin.
	 *
	 *  Note that with this function the stacking is performed in all the
	 *  redshift and proxy bins, and the results are written on file. 
	 *  In the first line of the header of such file, all the 
	 *  parameters used for the stacking (colour and redshift selections,
	 *  binnings, ...) are written, as well as the cosmological parameters. 
	 *  If such a file has already been written, the code reads it instead
	 *  of performing again the stacking procedure.
	 *
	 *  In addition, this function creates the following folders:
	 *
	 *  - covariance/: it contains \f$N\f$ files, where \f$N\f$ is
	 *  the number of redshift-proxy bins. Each file name contains
	 *  a string of two integers (e.g. "00"), where the first integer
	 *  corresponds to the index of the redshift bin, while the second
	 *  integer represents the index of the proxy bin.
	 *
	 *  - background_galaxies/: it contains sub-folders, one for each 
	 *  redshift bin. In such folders, two files are stored: one containing
	 *  the indices of the galaxies selected with the redshift selection,
	 *  while the second file contains the indices of those galaxies selected
	 *  through the colour selection. Such indices correspond to those
	 *  of the input catalogue of galaxies, and start from 0.
	 *
	 *  @param z_proxy_bin vector containing the indices of the redshift
	 *  and proxy bins to be stored, respectively
	 *
	 *  @param output_dir output directory for the output_file
	 *
	 *  @param output_file_root root name of the output file, which contains the
	 *  stacked profiles in all the redshift and mass 
	 *  proxy bins
	 *
	 *  @param errorType the type of error assigned to the 
	 *  density profile (only the diagonal of the covariance
	 *  matrix is considered)
	 *
	 *  @param n_resampling number of resampling regions for the bootstrap
	 *  procedure used to evaluate the uncertainty on \f$\Delta\Sigma(r)\f$
	 *
	 */
	void measure(const std::vector<int> z_proxy_bin, const std::string output_dir, const std::string output_file_root, const ErrorType errorType=ErrorType::_Bootstrap_, const int n_resampling=10000);
	
	///@}
	
	/**
	 *  @name input/output member functions (customized in all the derived classes)
	 */
	///@{
	
	/**
	 *  @brief write on file the measure in a given
	 *  bin of redshift and mass proxy
	 *
	 *  @param dir output directory
	 *
	 *  @param file name of the output file
	 *
	 */
	void write(const std::string dir, const std::string file);
	
	///@}
	
      };
    }
  }
}

#endif
