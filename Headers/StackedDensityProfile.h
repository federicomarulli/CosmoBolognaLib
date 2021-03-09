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
	 *  @brief set the colour and photo-z selection, along with the
	 *  logic operator between them
	 *
	 */
	void m_set_selections();
	
	/**
	 *  @brief colour selection by Oguri et al. 2012
	 *
	 *  @param i_gal galaxy index
	 *
	 *  @return true if the conditions are satisfied
	 *
	 */	
	bool m_colourSelection_Oguri (const int i_gal);
	
	/**
	 *  @brief photo-z selection accounting for the odds
	 *
	 *  @param x cluster and galaxy indices, respctively
	 *
	 *  @return true if the conditions are satisfied
	 *
	 */
	bool m_zphotSelection_ODDS (const std::vector<int> x);
	
	/**
	 *  @brief photo-z selection not accounting for the odds
	 *
	 *  @param x cluster and galaxy indices, respctively
	 *
	 *  @return true if the conditions are satisfied
	 *
	 */
	bool m_zphotSelection_noODDS (const std::vector<int> x);
	
	/**
	 *  @brief set the logic operator between color and phot-z selections
	 *
	 *  @param sel the logic selection operator ("and" or "or")
	 *
	 */
	void m_set_logicSelection(const std::string sel);
	
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
	
	/// pointer to the colour selection function
	bool (StackedDensityProfile::*m_colourSel)(int);
	
	/// pointer to the photo-z selection function
	bool (StackedDensityProfile::*m_photzSel)(std::vector<int>);
	
	/// logic operator between colour and photo-z selection
	std::function<bool(const std::vector<bool>)> m_logicSelection;
	
	/// the colour selection
	std::string m_colour_sel;
	
	/// the photo-z selection
	std::string m_zphot_sel;
	
	/// the logic operator between colour and photo-z selections
	std::string m_logic_sel;
	
	/// minimum signal-to-noise
	double m_SN_min;
	
	/// pixel size (in deg)
	double m_pix_size;
	
	/// number of regions for the resampling
	double m_n_resampling;
	
	/// vector of the inputs stored as strings
	std::vector<std::string> m_inputs_to_str;
	
	/// slope for the observable weighted mean
	double m_rad_alpha;
	
	/// gamma slope for lensing-weighted observables
	double m_obs_gamma;
	
	/// the parameters for the redshift selection
	std::vector<double> m_zphot_sel_pars;
	
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
	   *  @param colour_sel colour selection criterion. 
	   *  The possible choices are : "Oguri" (Oguri et al. 2012).
	   *
	   *  @param zphot_sel photo-z selection criterion. 
	   *  The possible choices are : "Odds", "NoOdds".
	   *  If "NoOdds" is chosen, the selection is not based on the odds.
	   *
	   *  @param zphot_sel_pars array of 4 parameters for the 
	   *  phot-z selection, following this order : 1) additive term, \f$\Delta z\f$, to
	   *  the cluster redshift, i.e. how much the minimum galaxy redshift,
	   *  \f$z_g\f$, (in the tail of the relative posterior) must be higher 
	   *  than the cluster redshift, \f$z_c\f$, that is \f$z_g > z_c+\Delta z\f$
	   *  , 2) minimum value for the galaxy mean redshift, 
	   *  3) maximum value for the galaxy mean redshift, 4) minimum value
	   *  for the ODDS parameter (if the ODDS are considered).
	   *
	   *  @param logic_sel logic operator between
	   *  colour and phot-z selection criteria ("and" or "or")
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
	   *  @param n_resampling number of resampling regions for the bootstrap
	   *  procedure used to evaluate the uncertainty on \f$\Delta\Sigma(r)\f$
	   *
	   *  @param rad_alpha slope for the observable weighted mean
	   *
	   *  @param obs_gamma gamma slope for lensing-weighted observables
	   *
	   *  Cluster stacked density profile
	   */
	StackedDensityProfile (cosmology::Cosmology cosm, const catalogue::Catalogue gal_cat, const catalogue::Catalogue clu_cat, const std::string colour_sel, const std::string zphot_sel, std::vector<double> zphot_sel_pars, const std::string logic_sel, std::vector<double> z_binEdges, std::vector<std::vector<double>> proxy_binEdges, const double rad_min, const double rad_max, const int nRad, const bool log_rad, const double SN_min, const double pix_size, const int n_resampling=10000, const double rad_alpha=1., const double obs_gamma=1.);
	///@}
	
	/**
	 *  @name Member functions to measure the stacked profile
	 */
	///@{
	
	/**
	 *  @brief measure the stacked profiles in all the bins of redshift and 
	 *  mass proxy, providing in output the stacked profile in the redshift 
	 *  and proxy bins chosen through the parameter z_proxy_bin.
	 *
	 *  Note that with this function the stacking is performed in all the
	 *  redshift and proxy bins, and the results are written on file. 
	 *
	 *  In the first line of the header of such file, all the 
	 *  parameters used for the stacking (colour and redshift selections,
	 *  binnings, ...) are written, as well as the cosmological parameters. 
	 *
	 *  If such a file has already been written, the code reads it instead
	 *  of performing again the stacking procedure.
	 *
	 *  @param z_proxy_bin vector containing the indices of the redshift
	 *  and proxy bins to be stored, respectively
	 *
	 *  @param output_dir output directory for the output_file
	 *
	 *  @param output_file output file, containing the
	 *  stacked profiles in all the redshift and mass 
	 *  proxy bins
	 *
	 *  @param errorType the type of error assigned to the 
	 *  density profile (only the diagonal of the covariance
	 *  matrix is considered)
	 *
	 */
	void measure(const std::vector<int> z_proxy_bin, const std::string output_dir, const std::string output_file, const ErrorType errorType=ErrorType::_Bootstrap_);
	
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
