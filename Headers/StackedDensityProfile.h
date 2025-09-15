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
       *  i.e. \f$\Delta\Sigma(r)\f$ [\f$h\f$ M\f$_\odot\f$/pc\f$^2\f$], or the reduced shear.
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
	 *  @brief null colour selection
	 *
	 *  @param i_gal unused parameter
	 *
	 *  @param z_c unused parameter
	 *
	 *  @return false
	 *
	 */
	bool m_colourSelection_null (const int i_gal, const double z_c) {(void)i_gal; (void)z_c; return false;};
	
	/**
	 *  @brief colour selection by Oguri+12, based on gri photometry, 
	 *  valid for lenses at \f$z<0.7\f$
	 *
	 *  @param i_gal galaxy index
	 *
	 *  @param z_c cluster redshift
	 *
	 *  @return true if the galaxy satisfies the conditions,
	 *  otherwise it returns false
	 *
	 */
	bool m_colourSelection_Oguri12 (const int i_gal, const double z_c);
	
	/**
	 *  @brief colour selection by Lesci+23 based on griz photometry, 
	 *  valid in the lens redshift range \f$z\in[0.1,0.8]\f$
	 *
	 *  @param i_gal galaxy index
	 *
	 *  @param z_c cluster redshift
	 *
	 *  @return true if the galaxy satisfies the conditions,
	 *  otherwise it returns false
	 *
	 */
	bool m_colourSelection_Lesci23_griz (const int i_gal, const double z_c);
	
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
	
	/// if true, create the folder containing the indices of the background galaxies
	bool m_write_background;
	
	/// if true, create a file containing the signal of individual clusters
	bool m_write_individual_signal;
	
	/// "DeltaSigma" \f$\rightarrow\f$ measure the excess surface density; "g" \f$\rightarrow\f$ measure the reduced shear
	std::string m_probe;
	
	/// input cosmology
	std::shared_ptr<cosmology::Cosmology> m_cosmology;
	
	/// input galaxy catalogue
	std::shared_ptr<catalogue::Catalogue> m_galData;
	
	/// input cluster catalogue
	std::shared_ptr<catalogue::Catalogue> m_cluData;
	
	/// magnitude fail value
	double m_mag_fail = -99;
	
	/// vector of pointers to the colour selection functions
	std::vector<bool (StackedDensityProfile::*)(const int, const double)> m_colourSel;
	
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
	
	/// the colour selections authors
	std::vector<std::string> m_CCsel_author;
	
	/// boolean indicating whether the CC sel. includes conditions on the galaxy photo-z
	std::vector<bool> m_CCsel_conservative;
	
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
	
	/// indices of the clusters in bins of z and mass proxy
	std::vector<std::vector<std::vector<double>>> m_idxClu_inBin;
	
	/// number of galaxies for annulus
	std::vector<int> m_ngal_arr;
	
	/// vector containing the source redshifts (related to a single cluster)
	std::vector<double> m_z_source_single;
	
	/// vector of vectors containing the number of sources as a function of redshift for each stack
	std::vector<std::vector<double>> m_Nz_source_stack;
	
	/// vector of vectors containing the number of sources as a function of redshift for each cluster
	std::vector<std::vector<double>> m_Nz_source_single;
	
	/// source redshift bin edges for computing m_Nz_source_stack and m_Nz_source_single
	std::vector<double> m_z_source_bin_edges;
	
	/// number of galaxies for annulus for cluster
	std::vector<std::vector<std::vector<std::vector<int>>>> m_single_ngal_arr;
	
	/// stacked number of galaxies for annulus
	std::vector<std::vector<std::vector<int>>> m_stacked_ngal_arr;
	
	/// tangential surface density
	std::vector<double> m_deltasigma_t;
	
	/// effective radius, stored for each single cluster
	std::vector<std::vector<std::vector<std::vector<double>>>> m_single_rad_eff;
	
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
	
	/// stacked wSum
	std::vector<std::vector<std::vector<double>>> m_stacked_wSum;
	
	/// stacked wetaSum
	std::vector<std::vector<std::vector<double>>> m_stacked_wetaSum;
	
	/// term including the multiplicative shear bias
	std::vector<double> m_K;
	
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
	   *  @brief Constructor used to perform the signal stacking.
	   *  Note that the function
	   *  cbl::measure::stackprofile::StackedDensityProfile::measure
	   *  can also write a file containing the signal of single
	   *  clusters and the associated background source
	   *  distributions, if requested.
	   *
	   *  By default, the excess surface density is measured.
	   *  There is the possibility to compute the reduced shear.
	   *  See the description of the variable named "probe".
	   *
	   *  The observed ellipticity components, \f$e_1\f$ and \f$e_2\f$,
	   *  are separated into a tangential component, \f$e_+\f$,
	   *  and a cross component, \f$e_\times\f$, as follows
	   *  (see, e.g., Viola et al. 2015)
	   *
	   *  \f$ e_+ = -e_1\cos(2\varphi) -e_2\sin(2\varphi)\, \f$
	   *
	   *  \f$ e_\times = e_1\sin(2\varphi) -e_2\cos(2\varphi)\, \f$
	   *
	   *  where \f$\varphi\f$ is the position angle of the source 
	   *  with respect to the lens centre. The estimator of the
	   *  tangential excess surface density, \f$\Delta\Sigma_+\f$,
	   *  is expressed as follows
	   *
	   *  \f$ \Delta\Sigma_+(R_j)=\Sigma_{\rm cr}\,g_+(R_j) =
	   *  \left( \frac{\sum_{i\in j}\;(w_i\,\Sigma_{{\rm cr},i}^{-2})\,
	   *  e_{+,i}\,\Sigma_{{\rm cr},i}}{\sum_{i\in j}\;
	   *  (w_i\,\Sigma_{{\rm cr},i}^{-2})} \right)\frac{1}{1+K(R_j)}\,, \f$
	   *
	   *  where \f$\Sigma_{\rm cr}\f$ is the critical density, \f$g_+\f$
	   *  is the tangential reduced shear,
	   *  \f$j\f$ is the radial annulus index, with mean radius \f$R_j\f$, 
	   *  and \f$w_i\f$ is the statistical weight assigned to the measure of the 
	   *  source ellipticity of the \f$i\f$th galaxy (Sheldon et al. 2004).
	   *  In addition, \f$K(R_j)\f$ is the average correction due to the 
	   *  multiplicative noise bias in the shear estimate, defined as
	   *
	   *  \f$ K(R_j)=\frac{\sum_{i\in j}\;(w_i\,\Sigma_{{\rm cr},i}^{-2})\,m_i}
	   *  {\sum_{i\in j}\;(w_i\,\Sigma_{{\rm cr},i}^{-2})}, \f$
	   *
	   *  where \f$m_i\f$ is the multiplicative shear bias of the \f$i\f$th 
	   *  galaxy. The stacked excess surface density for the \f$K\f$th 
	   *  cluster bin is expressed as
	   *
	   *  \f$ \Delta\Sigma_{+,K}(R_j) = \frac{\sum_{k\in K}W_{k,j}\,
	   *  \Delta\Sigma_{+,k}(R_j)}{\sum_{k\in K}W_{k,j}}, \f$
	   *
	   *  where \f$k\f$ runs over the clusters in the bin, and \f$W_{k,j}\f$ is the 
	   *  total weight for the \f$j\f$th radial bin of the \f$k\f$th cluster, estimated as
	   *
	   *  \f$ W_{k,j} = \sum_i w_i\,\Sigma_{{\rm cr},i}^{-2}, \f$
	   *
	   *  where \f$i\f$ runs over the background galaxies in the \f$j\f$th radial bin. 
	   *  For each stack, we derive the effective cluster mass proxy, 
	   *  \f$\mathcal{O}_{K,\rm eff}\f$, as follows (Umetsu et al. 2014)
	   *
	   *  \f$ O_{K,\rm eff} = \frac{\sum_{k\in K}\,W_kO_k}{\sum_{k\in K}\,W_k}\,, \f$
	   *
	   *  where \f$W_k\f$ does not account for the subdivision in radial bins.
	   *  The effective redshift and radius of the stack are derived analogously.
	   *
	   *  The reduced shear estimator is expressed as follows
	   *
	   *  \f$ g_+(R_j) = \left( \frac{\sum_{i\in j}\,w_i\,e_{+,i}}
	   *  {\sum_{i\in j}w_i} \right)\frac{1}{1+(\sum_{i\in j}\,w_im_i)
	   *  \,/\,(\sum_{i\in j}\,w_i)}\,, \f$
	   *
	   *  and the stacked reduced shear for the \f$K\f$th 
	   *  cluster bin is expressed as
	   *
	   *  \f$ g_{+,K}(R_j) = \frac{\sum_{k\in K}\tilde{W}_{k,j}\,
	   *  g_{+,k}(R_j)}{\sum_{k\in K}\tilde{W}_{k,j}}, \f$
	   *
	   *  where
	   *
	   *  \f$ \tilde{W}_{k,j} = \sum_i w_i\,. \f$
	   *
	   *  In this case, the effective mass proxy is a weighted sum analogous
	   *  to that presented above, where \f$W_{k,j}\f$ is replaced
	   *  by \f$\tilde{W}_{k,j}\f$. The same holds for the effective
	   *  redshift and radius.
	   *
	   *  @param cosmology cosmological model
	   *
	   *  @param gal_cat catalogue of galaxies
	   *
	   *  @param clu_cat catalogue of galaxy clusters
	   *
	   *  @param delta_redshift lens redshift buffer considered in both 
	   *  photo-\f$z\f$ and colour selections. For photo-\f$z\f$ selections,
	   *  it is \f$\Delta_z\f$ in the equation 
	   *  \f$z_{\rm g}>z_{\rm l}+\Delta_z\f$, where \f$z_{\rm g}\f$ and 
	   *  \f$z_{\rm l}\f$ are the galaxy and cluster 
	   *  mean redshifts, respectively. For colour selections,
	   *  the considered lens redshift is \f$z_{\rm l}'=z_{\rm l}+\Delta z\f$.
	   *  For details,
	   *  see cbl::measure::stackprofile::StackedDensityProfile::set_zphot_selection
	   *  and cbl::measure::stackprofile::StackedDensityProfile::set_colour_selection
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
	   *  @param pix_size linked list pixel size in deg. 1 degree
	   *  should be optimal
	   *
	   *  @param multiplicative_calibration_stats a vector of
	   *  vectors containing mean and standard deviation of the
	   *  multiplicative shear calibration parameter, usually
	   *  denoted as m, in all the redshift bins where it is
	   *  evaluated.  If not provided, the galaxy-by-galaxy values
	   *  of m are used
	   *
	   *  @param multiplicative_calibration_zEdges vector of
	   *  vectors containing the lower and upper edge of the
	   *  redshift bins for each estimate of the multiplicative
	   *  shear calibration parameter, m
	   *
	   *  @param rad_alpha slope for the observable weighted mean
	   *
	   *  @param obs_gamma gamma slope for lensing-weighted observables
	   *
	   *  @param probe "DeltaSigma" \f$\rightarrow\f$ measure the excess
	   *  surface density; "g" \f$\rightarrow\f$ measure the reduced shear
	   *
	   *  @warning Note that cbl::catalogue::Var::\_Shear1\_ and
	   *  cbl::catalogue::Var::\_Shear2\_, in the input galaxy
	   *  catalogue, here are interpreted as ellipticity
	   *  components.
	   *
	   *  @warning This constructor does not set any colour or
	   *  redshift background selections.  By default, colour and
	   *  redshift selections are linked by the "or" logic
	   *  operator. To change these settings, see
	   *  cbl::measure::stackprofile::StackedDensityProfile::set_colour_selection,
	   *  cbl::measure::stackprofile::StackedDensityProfile::set_zphot_selection,
	   *  cbl::measure::stackprofile::StackedDensityProfile::set_logic_selection
	   *
	   *  @warning If the probe is \f$\Delta\Sigma\f$, only the galaxies with mean
	   *  redshift larger than the lens mean redshift are considered.
	   */
	StackedDensityProfile (const std::shared_ptr<cosmology::Cosmology> cosmology, std::shared_ptr<catalogue::Catalogue> gal_cat, std::shared_ptr<catalogue::Catalogue> clu_cat, const double delta_redshift, std::vector<double> z_binEdges, std::vector<std::vector<double>> proxy_binEdges, const double rad_min, const double rad_max, const int nRad, const bool log_rad, const double SN_min, const double pix_size=1., const std::vector<std::vector<double>> multiplicative_calibration_stats={}, const std::vector<std::vector<double>> multiplicative_calibration_zEdges={}, const double rad_alpha=1., const double obs_gamma=1., const std::string probe="DeltaSigma");
	
	///@}

	
	/**
	 *  @name Member functions to measure the stacked profile
	 */
	///@{
	
	/**
	 *  @brief Set the colour selection in the i-th cluster
	 *  redshift bin. By default, the class constructor sets no
	 *  colour selections.
	 *
	 *  The lens redshift considered in the colour selections is
	 *  \f$z_{\rm l}' = z_{\rm l} + \Delta z\f$, where \f$\Delta z\f$
	 *  is the delta_redshift parameter set through the constructor 
	 *  cbl::measure::stackprofile::StackedDensityProfile
	 *
	 *  @param author the colour selection author; possibilities
	 *  are: "Oguri12" (Eq.s (4)-(6) in
	 *  https://ui.adsabs.harvard.edu/abs/2012MNRAS.420.3213O/abstract),
	 *  valid for lenses at \f$z<0.7\f$; "Lesci23_griz"
	 *  (Table A.1 in https://ui.adsabs.harvard.edu/abs/2023arXiv231116239E/abstract),
	 *  valid in the lens redshift range \f$z\in[0.1,0.8]\f$
	 *
	 *  @param z_bin index of the cluster redshift bin
	 *
	 *  @param mag_fail value corresponding to magnitude measurement fails.
	 *  It is the same for all redshift bins. Fails
	 *  are converted into NaNs
	 *
	 *  @param conservative_selection if true, the candidate background galaxies
	 *  must satisfy also the condition \f$z_{\rm g} > z_{\rm l} + \Delta z\f$, 
	 *  where \f$z_{\rm g}\f$ is the mean galaxy redshift
	 */
	void set_colour_selection(const std::string author, const int z_bin, const double mag_fail=-99, const bool conservative_selection=false);
	
	/**
	 *  @brief Set the colour selection in all cluster redshift bins.
	 *  By default, the class constructor sets no colour selections.
	 *
	 *  The lens redshift considered in the colour selections is
	 *  \f$z_{\rm l}' = z_{\rm l} + \Delta z\f$, where \f$\Delta z\f$
	 *  is the delta_redshift parameter set through the constructor 
	 *  cbl::measure::stackprofile::StackedDensityProfile
	 *
	 *  @param author the colour selection author; possibilities are: 
	 *  "Oguri12" (Eq.s (4)-(6) in
	 *  https://ui.adsabs.harvard.edu/abs/2012MNRAS.420.3213O/abstract), 
	 *  valid for lenses at \f$z<0.7\f$; 
	 *  "Lesci23_griz"
	 *  (Table A.1 in https://ui.adsabs.harvard.edu/abs/2023arXiv231116239E/abstract),
	 *  valid in the lens redshift range \f$z\in[0.1,0.8]\f$
	 *
	 *  @param mag_fail value corresponding to magnitude measurement fails. Fails
	 *  are converted into NaNs
	 *
	 *  @param conservative_selection if true, the candidate background galaxies
	 *  must satisfy also the condition \f$z_{\rm g} > z_{\rm l} + \Delta z\f$, 
	 *  where \f$z_{\rm g}\f$ is the mean galaxy redshift
	 */
	void set_colour_selection(const std::string author, const double mag_fail=-99, const bool conservative_selection=false);
	
	/**
	 *  @brief Set the redshift selection in the i-th cluster redshift bin.
	 *  By default, the class constructor sets no redshift selections.
	 *
	 *  The redshift selection is expressed as
	 *
	 *  \f$ z_{\rm g,\,min} > z_{\rm l} + (\Delta z)' \,\,\land \f$
	 *
	 *  \f$ {\rm ODDS} > {\rm ODDS}_{\rm min} \,\,\land \f$
	 *
	 *  \f$ z_{\rm g} \in [z_{\rm g}^{\rm low},z_{\rm g}^{\rm up}]\,\,\land \f$
	 *  
	 *  \f$ z_{\rm g} > z_{\rm l} + \Delta z \,, \f$
	 *
	 *  where \f$ z_{\rm g,\,min}\f$ is the minimum galaxy
	 *  redshift excluding the left PDZ tail, \f$(\Delta z)'\f$ is
	 *  the cluster redshift buffer, while \f$ z_{\rm g} \f$ is the
	 *  galaxy most probable redshift value. \f$\Delta z\f$ is the delta_redshift
	 *  parameter set through the constructor 
	 *  cbl::measure::stackprofile::StackedDensityProfile
	 *
	 *  @param z_bin index of the cluster redshift bin
	 *
	 *  @param deltaz additive term, \f$(\Delta z)'\f$, to
	 *  the cluster redshift, i.e. how much the minimum galaxy redshift,
	 *  \f$z_{\rm g,\,min}\f$, (in the tail of the relative posterior) must be higher 
	 *  than the cluster redshift, \f$z_{\rm l}\f$, that is \f$z_{\rm g,\,min} > z_{\rm l}+(\Delta z)'\f$
	 *
	 *  @param zgal_min minimum value for the galaxy most probable redshift, \f$ z_{\rm g}^{\rm low} \f$
	 *
	 *  @param zgal_max maximum value for the galaxy most probable redshift, \f$ z_{\rm g}^{\rm up} \f$
	 *
	 *  @param ODDS_min minimum value for the ODDS parameter
	 *
	 */
	 
	void set_zphot_selection(const int z_bin, const double deltaz, const double zgal_min, const double zgal_max, const double ODDS_min);
	
	/**
	 *  @brief Set the redshift selection in all cluster redshift bins.
	 *  By default, the class constructor sets no redshift selections.
	 *
	 *  The redshift selection is expressed as
	 *
	 *  \f$ z_{\rm g,\,min} > z_{\rm l} + (\Delta z)' \,\,\land \f$
	 *
	 *  \f$ {\rm ODDS} > {\rm ODDS}_{\rm min} \,\,\land \f$
	 *
	 *  \f$ z_{\rm g} \in [z_{\rm g}^{\rm low},z_{\rm g}^{\rm up}]\,\,\land \f$
	 *  
	 *  \f$ z_{\rm g} > z_{\rm l} + \Delta z \,, \f$
	 *
	 *  where \f$ z_{\rm g,\,min}\f$ is the minimum galaxy
	 *  redshift excluding the left PDZ tail, \f$(\Delta z)'\f$ is
	 *  the cluster redshift buffer, while \f$ z_{\rm g} \f$ is the
	 *  galaxy most probable redshift value. \f$\Delta z\f$ is the delta_redshift
	 *  parameter set through the constructor 
	 *  cbl::measure::stackprofile::StackedDensityProfile
	 *
	 *  @param deltaz additive term, \f$(\Delta z)'\f$, to
	 *  the cluster redshift, i.e. how much the minimum galaxy redshift,
	 *  \f$z_{\rm g,\,min}\f$, (in the tail of the relative posterior) must be higher 
	 *  than the cluster redshift, \f$z_{\rm l}\f$, that is \f$z_{\rm g,\,min} > z_{\rm l}+(\Delta z)'\f$
	 *
	 *  @param zgal_min minimum value for the galaxy most probable redshift, \f$ z_{\rm g}^{\rm low} \f$
	 *
	 *  @param zgal_max maximum value for the galaxy most probable redshift, \f$ z_{\rm g}^{\rm up} \f$
	 *
	 *  @param ODDS_min minimum value for the ODDS parameter
	 *
	 */
	 
	void set_zphot_selection(const double deltaz, const double zgal_min, const double zgal_max, const double ODDS_min);
	
	/**
	 *  @brief Set the logic operator linking redshift and colour
	 *  selections in the i-th cluster redshift bin.
	 *  By default, the class constructor sets the "or" operator
	 *
	 *  @param z_bin index of the cluster redshift bin
	 *
	 *  @param sel "and" or "or"
	 *
	 */
	 
	void set_logic_selection(const int z_bin, const std::string sel);
	
	/**
	 *  @brief Set the logic operator linking redshift and colour
	 *  selections in the all the cluster redshift bins.
	 *  By default, the class constructor sets the "or" operator
	 *
	 *  @param sel "and" or "or"
	 *
	 */
	 
	void set_logic_selection(const std::string sel);
	
	/**
	 *  @brief The files produced through this function may be very large.
	 *  Be careful!
	 *
	 *  This function measures the stacked profiles in all the bins of redshift and 
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
	 *  In addition, this function creates the following outputs:
	 *
	 *  - covariance/: folder containing \f$N\f$ files, where \f$N\f$ is
	 *  the number of redshift-proxy bins. Each file name contains
	 *  a string of two integers (e.g. "00"), where the first integer
	 *  corresponds to the index of the redshift bin, while the second
	 *  integer represents the index of the proxy bin.
	 *
	 *  - clusters_information/: folder containing \f$N\f$ files following the same
	 *  naming convention described above. These files contain information on
	 *  the clusters contributing to the lensing signal in each stack.
	 *  The clusters without background galaxies are not included. Such files
	 *  are useful when 
	 *  cbl::modelling::densityprofile::Modelling_DensityProfile::set_data_model
	 *  is used.
	 *
	 *  - background_galaxies/: folder created if write_background is set to true.
	 *  This folder contains sub-folders, one for each 
	 *  redshift bin. In such folders, two files are stored: one containing
	 *  the indices of the galaxies selected with the redshift selection,
	 *  while the second file contains the indices of those galaxies selected
	 *  through the colour selection. Such indices correspond to those
	 *  of the input catalogue of galaxies, and start from 0.
	 *
	 *  - individual_signal.fits: fits file containing the signal of individual
	 *  clusters, created only if write_individual_signal is true.
	 *  Only the clusters having at least one background galaxy are included.
	 *  The first entries represent R.A., Dec, redshift, mass proxy, and
	 *  all the effective radii, in this order, followed by
	 *  the entries representing the tangential signal, the cross signal, and the
	 *  statistical uncertainty, in this order. Of course, the number of entries
	 *  depends on the requested number of radial bins.
	 *
	 *  - individual_source_distribution.fits: fits file containing the 
	 *  background source distribution of individual
	 *  clusters, created only if write_individual_signal is true.
	 *  Only the clusters having at least one background galaxy are included.
	 *  The first column contains the central values of the source redshift
	 *  in the bin. Each subsequent column corresponds to the number of
	 *  sources as a function of redshift for a single cluster. Then, the
	 *  number of columns depends on the number of input clusters. The 
	 *  \f$i\f$th column in this file, with \f$i>1\f$ to exclude the first
	 *  column containing the central source redshift values, corresponds
	 *  to the \f$(i-1)\f$th row in individual_signal.fits.
	 *
	 *  @param z_proxy_bin vector containing the indices of the redshift
	 *  and proxy bins to be stored, respectively
	 *
	 *  @param output_dir output directory for the output_file
	 *
	 *  @param output_file_root root name of the output file, which contains the
	 *  stacked profiles in all the redshift and mass 
	 *  proxy bins. The same root is used to produce a file including the
	 *  source redshift distribution for each cluster redshift bin adopted
	 *  for the stacking.
	 *
	 *  @param errorType the type of error assigned to the 
	 *  density profile (only the diagonal of the covariance
	 *  matrix is considered)
	 *
	 *  @param n_resampling number of resampling regions for the bootstrap
	 *  procedure used to evaluate the uncertainty on \f$\Delta\Sigma(r)\f$
	 *
	 *  @param write_background if true, create the background_galaxies/ folder.
	 *  WARNING: this folder may be very large
	 *
	 *  @param write_individual_signal if true, create the file containing 
	 *  the signal of individual clusters, named individual_signal.fits,
	 *  and another one containing the distribution of sources for each cluster,
	 *  named individual_source_distribution.fits. 
	 *  Only the clusters having at least one background galaxy are included.
	 *  WARNING: these files may be very large
	 *
	 */
	void measure(const std::vector<int> z_proxy_bin, const std::string output_dir, const std::string output_file_root, const ErrorType errorType=ErrorType::_Bootstrap_, const int n_resampling=10000, const bool write_background=false, const bool write_individual_signal=false);
	
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
	
	/**
	 *  @brief test the chosen colour selection on all the
	 *  input source galaxies. This function does not account
	 *  for the lens dimension on the sky
	 *
	 *  @param z_lim redshift lower limit
	 *
	 *  @return vector containing completeness, purity, and
	 *  foreground fraction
	 *
	 */
	std::vector<double> check_colour_selection(const double z_lim);
	
	///@}
	
      };
    }
  }
}

#endif
