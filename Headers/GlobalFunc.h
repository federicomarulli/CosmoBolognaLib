/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli                         *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file Headers/GlobalFunc.h
 *
 *  @brief Generic functions that use one or more classes of the
 *  CosmoBolognaLib
 *
 *  This file contains the prototypes of a set of generic functions
 *  that use one or more classes of the CosmoBolognaLib
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __GLOBALFUNC__
#define __GLOBALFUNC__


#include "ThreePointCorrelation.h"


// ============================================================================================


namespace cbl {

  /**
   *  @name Generic functions that use the class Catalogue
   */
  ///@{

  /**
   *  @brief overloading of the + operator, to sum two catalogues
   *  @param c1 object of class Catalogue
   *  @param c2 object of class Catalogue
   *  @return object of class catalogue
   */
  inline catalogue::Catalogue operator + (const catalogue::Catalogue &c1, const catalogue::Catalogue &c2)
    {
      catalogue::Catalogue ctemp = c1;
      return ctemp += c2;
    }
  
  /**
   *  @name Generic functions that use the class Cosmology
   */
  ///@{

  /**
   *  @brief get a smoothed distribution of comoving distances,
   *  estimated with the V<SUB>max</SUB> method
   *
   *  this function is useful to construct random catalogues
   *
   *  @param [out] dc d<SUB>c</SUB>: vector containing the output
   *  values of the binned comoving distances 
   *
   *  @param [out] nObj vector containing the smoothed number of
   *  objects at each d<SUB>c</SUB>, estimated with the
   *  V<SUB>max</SUB> method
   *
   *  @param [in] D_C vector containing the comoving distances of the
   *  objects
   *
   *  @param [in] zobj_min minimum redshift of the objects
   *
   *  @param [in] zobj_max maximum redshift of the objects
   *
   *  @param [in] z_min minimum redshift used to enlarge the range,
   *  useful to smooth the redshift distributions
   *
   *  @param [in] z_max maximum redshift used to enlarge the range,
   *  useful to smooth the redshift distributions
   *
   *  @param [in] zbin_min minimum redshift of the output binning
   *
   *  @param [in] zbin_max maximum redshift of the output binning
   *
   *  @param [in] cosm object of class Cosmology
   *
   *  @param [in] Area area of the survey
   *
   *  @param [in] nObjRan number of random objects used to assess the
   *  smoothed distribution
   *
   *  @param [in] norm 0 &rarr; don't normalize; 1 &rarr; normalize
   *  the distribution to the number of objects
   *
   *  @param [in] file_Vmax the output file used to store the
   *  smoothed distribution
   *
   *  @param [in] delta_dc_Vmax &Delta;d<SUB>c</SUB>: the bin size of
   *  the output smoothed distribution
   *
   *  @param [in] seed the random seed
   *
   *  @return none
   */
  void Vmax_DC_distribution (std::vector<double> &dc, std::vector<double> &nObj, const std::vector<double> D_C, const std::vector<double> zobj_min, const std::vector<double> zobj_max, const double z_min, const double z_max, const double zbin_min, const double zbin_max, cosmology::Cosmology &cosm, const double Area, const int nObjRan, const bool norm=1, const std::string file_Vmax=par::defaultString, const double delta_dc_Vmax=100., const int seed=3213);

  /**
   *  @brief the Alcock-Pacinski factor used to shift comoving
   *  distances
   *
   *  this function is used to model geometric distortions
   *
   *  @param redshift the redshift
   *  @param cosm1 object of class Cosmology
   *  @param cosm2 object of class Cosmology
   *
   *  @return D<SUB>V</SUB>[cosm2]/D<SUB>V</SUB>[cosm1], where
   *  D<SUB>V</SUB> is Cosmology::D_V
   */
  double AP_shift_r (const double redshift, const cosmology::Cosmology &cosm1, const cosmology::Cosmology &cosm2);

  /**
   *  @brief the Alcock-Pacinski factor used to shift comoving
   *  distances perpendicular to the line-of-sight, r<SUB>p</SUB>
   *
   *  this function is used to model geometric distortions
   *
   *  @param redshift the redshift
   *  @param cosm1 object of class Cosmology
   *  @param cosm2 object of class Cosmology
   *
   *  @return D<SUB>A</SUB>[cosm1]/D<SUB>A</SUB>[cosm2], where
   *  D<SUB>A</SUB> is Cosmology::D_A
   */
  double AP_shift_rp (const double redshift, const cosmology::Cosmology &cosm1, const cosmology::Cosmology &cosm2);

  /**
   *  @brief the Alcock-Pacinski factor used to shift comoving
   *  distances parallel to the line-of-sight, &pi;
   *
   *  this function is used to model geometric distortions
   *
   *  @param redshift the redshift
   *  @param cosm1 object of class Cosmology
   *  @param cosm2 object of class Cosmology
   *
   *  @return H[cosm2]/H[cosm1], where H is Cosmology::HH
   */
  double AP_shift_pi (const double redshift, const cosmology::Cosmology &cosm1, const cosmology::Cosmology &cosm2);

  /**
   *  @brief the maximum comoving separations to be used for the AP
   *  test, for a given set of test cosmologies
   *
   *  @param [in] Rp_max the maximum value of the comoving distance
   *  perpendicular to the line-of-sight, r<SUB>p,MAX</SUB>, in the
   *  cosmology cosm1
   *
   *  @param [in] Pi_max the maximum value of the comoving distance
   *  parallel to the line-of-sight, &pi;<SUB>MAX</SUB>, in the
   *  cosmology cosm1
   *
   *  @param [in] redshift the redshift
   *
   *  @param [in] cosm1 object of class Cosmology: the assumed (or real)
   *  cosmology
   *
   *  @param [in] cosm2 a vector of objects of class Cosmology: the test
   *  cosmologies
   *
   *  @param [out] rpM_AP the maximum value of the comoving distance
   *  perpendicular to the line-of-sight, r<SUB>p,MAX</SUB>, over all
   *  the test cosmologies cosm2
   *
   *  @param [out] piM_AP the maximum value of the comoving distance
   *  parallel to the line-of-sight, &pi;<SUB>MAX</SUB>, over all the
   *  test cosmologies cosm2
   *
   *  @param [out] rM_AP the maximum value of the comoving distance
   *  over all the test cosmologies cosm2
   * 
   *  @return none
   */
  void max_separations_AP (const double Rp_max, const double Pi_max, const double redshift, const cosmology::Cosmology &cosm1, const std::vector<cosmology::Cosmology> &cosm2, double &rpM_AP, double &piM_AP, double &rM_AP);

  /**
   *  @brief the 1D two-point correlation function converted from one
   *  cosmology to another one
   *  
   *  @param RR the comoving separation, R
   *
   *  @param redshift the redshift
   *
   *  @param rr vector containing the comoving separations, r
   *
   *  @param Xi vector containing the two-point correlation function,
   *  &xi;(r)
   *  
   *  @param cosm1 object of class Cosmology
   *
   *  @param cosm2 object of class Cosmology 
   *
   *  @param direction 0 &rarr; cosm2 \f$ \rightarrow \f$ cosm1; 1 &rarr; cosm1
   *  \f$ \rightarrow \f$ cosm2;
   *
   *  @return the converted two-point correlation function, &xi;(R)
   */
  double converted_xi (const double RR, const double redshift, const std::vector<double> rr, const std::vector<double> Xi, const cosmology::Cosmology &cosm1, const cosmology::Cosmology &cosm2, const bool direction);

  /**
   *  @brief the 2D two-point correlation function converted from one
   *  cosmology to another one
   *  
   *  @param RP the comoving separation perpendicular to the
   *  line-of-sight, R<SUB>p</SUB>
   *
   *  @param PI the comoving separation parallel to the line-of-sight,
   *  &Pi;
   *
   *  @param redshift the redshift
   *
   *  @param rp vector containing the comoving separations
   *  perpendicular to the line-of-sight, r<SUB>p</SUB>
   *
   *  @param pi vector containing the comoving separations
   *  parallel to the line-of-sight, &pi;
   *
   *  @param Xi matrix containing the two-point correlation function,
   *  &xi;(r<SUB>p</SUB>,&pi;)
   *  
   *  @param cosm1 object of class Cosmology
   *
   *  @param cosm2 object of class Cosmology 
   *
   *  @param direction 0 &rarr; cosm2 \f$ \rightarrow \f$ cosm1; 1 &rarr; cosm1
   *  \f$ \rightarrow \f$ cosm2;
   *
   *  @return the converted two-point correlation function, &xi;(R<SUB>p</SUB>,&Pi;)
   */
  double converted_xi (const double RP, const double PI, const double redshift, const std::vector<double> rp, const std::vector<double> pi, const std::vector<std::vector<double> > Xi, const cosmology::Cosmology &cosm1, const cosmology::Cosmology &cosm2, const bool direction); 

  ///@}

  
  /**
   *  @name Generic functions that use one or more classes of the CosmoBolognaLib
   */

  ///@{

  /**
   *  @brief compute the redsfhit range of a simulation box centered
   *  at z=mean_redshift
   *
   *  @param [in] mean_redshift the mean redshift
   *  @param [in] boxSide the box side
   *  @param [in] real_cosm an object of class Cosmology
   *  @param [out] redshift_min the minimum redshift
   *  @param [out] redshift_max the maximum redshift
   *  @return none
   */
  void redshift_range (const double mean_redshift, const double boxSide, cosmology::Cosmology &real_cosm, double &redshift_min, double &redshift_max); 

  /**
   *  @brief get the volume of a simulation box
   *
   *  @param boxSize the box side
   *
   *  @param frac the side fraction (if the input box is a sub-box of
   *  a box with side equal to boxSize)
   *
   *  @param Bord the redshift interval that is cutted at the bords of
   *  the box
   *
   *  @param mean_redshift the mean redshift
   *
   *  @param real_cosm an object of class Cosmology
   *
   *  @return the volume of the simulation
   */
  double volume (const double boxSize, const int frac, const double Bord, const double mean_redshift, cosmology::Cosmology &real_cosm);

  /**
   *  @brief convert a set of coordinates from real-space to
   *  redshift-space
   *
   *  @param[in,out] ra vector containing the Right Ascensions 
   *
   *  @param[in,out] dec vector containing the Declinations
   *
   *  @param[in,out] redshift vector containing the redshifts
   *
   *  @param[in,out] xx vector containing the x coordinates
   *
   *  @param[in,out] yy vector containing the y coordinates
   *
   *  @param[in,out] zz vector containing the z coordinates
   *
   *  @param[in] vx vector containing the peculiar velocities
   *  along the x direction
   *
   *  @param[in] vy vector containing the peculiar velocities
   *  along the y direction
   *
   *  @param[in] vz vector containing the peculiar velocities
   *  along the z direction
   *
   *  @param[in] sigmaV the error on the peculiar velocities
   *
   *  @param[in] real_cosm an object of class Cosmology
   *
   *  @param[in] mean_redshift the mean redshift
   *
   *  @param[in] redshift_min the minimum redshift
   *
   *  @param[in] redshift_max the maximum redshift
   *
   *  @param[in] seed the random seed
   *
   *  @return none
   */
  void coord_zSpace (std::vector<double> &ra, std::vector<double> &dec, std::vector<double> &redshift, std::vector<double> &xx, std::vector<double> &yy, std::vector<double> &zz, const std::vector<double> vx, const std::vector<double> vy, const std::vector<double> vz, const double sigmaV, cosmology::Cosmology &real_cosm, const double mean_redshift, const double redshift_min, const double redshift_max, const int seed=3213);

  /**
   *  @brief create a mock catalogue, subdividing a box into sub-boxes
   *  and recentering
   *
   *  @param [in] xx vector containing the x coordinates
   *
   *  @param [in] yy vector containing the y coordinates
   *
   *  @param [in] zz vector containing the z coordinates
   *
   *  @param [in] vx vector containing the peculiar velocities along
   *  the x direction
   *
   *  @param [in] vy vector containing the peculiar velocities along
   *  the y direction
   *
   *  @param [in] vz vector containing the peculiar velocities along
   *  the z direction
   *
   *  @param [in] var1 vector containing a generic quantity
   *  (e.g. masses or luminosities)
   *
   *  @param [in] var2 vector containing a generic quantity
   *  (e.g. masses or luminosities)
   *
   *  @param [in] var3 vector containing a generic quantity
   *  (e.g. masses or luminosities)
   *
   *  @param [in] output_dir name of directory used to store the outputs
   *
   *  @param [in] boxSize the box side
   *
   *  @param [in] frac the side fraction (if the input box is a
   *  sub-box of a box with side equal to boxSize)
   *
   *  @param [in] Bord the redshift interval that is cutted at the
   *  bords of the box
   *
   *  @param [in] mean_redshift the mean redshift
   *
   *  @param [in] real_cosm an object of class Cosmology
   *
   *  @param [in] REAL 0 &rarr; redshift-space; 1 &rarr; real-space
   *
   *  @param [in] sigmaV the error on the peculiar velocities
   *
   *  @param [in] idum the random seed
   *
   *  @param [out] Volume the mock volume
   *
   *  @return none
   */
  void create_mocks (const std::vector<double> xx, const std::vector<double> yy, const std::vector<double> zz, const std::vector<double> vx, const std::vector<double> vy, const std::vector<double> vz, const std::vector<double> var1, const std::vector<double> var2, const std::vector<double> var3, const std::string output_dir, const double boxSize, const int frac, const double Bord, const double mean_redshift, cosmology::Cosmology &real_cosm, const int REAL, const double sigmaV, const int idum, double &Volume);

  /**
   *  @brief set the object region in sub-boxes
   *  @param data input data catalogue
   *  @param nx side fraction used to divide the box in the x direction 
   *  @param ny side fraction used to divide the box in the y direction 
   *  @param nz side fraction used to divide the box in the z direction 
   *  @return none
   */
  void set_ObjectRegion_SubBoxes (catalogue::Catalogue &data, const int nx, const int ny, const int nz);

  /**
   *  @brief set the object region in angular SubBoxes
   *  @param data input data catalogue
   *  @param Cell_size size of the cell in degrees
   *  @return none
   */
  void set_ObjectRegion_RaDec (catalogue::Catalogue &data, const double Cell_size);

  /**
   *  @brief set the object region in sub-regions using mangle
   *  @param data input data catalogue
   *  @param nSamples number of sub-regions
   *  @param polygonfile name of the input file with polygons
   *  @return none
   */
  void set_ObjectRegion_mangle (catalogue::Catalogue &data, const int nSamples, const std::string polygonfile);

  /**
   *  @brief set the object region in sub-boxes
   *  @param data input data catalogue
   *  @param random random catalogue
   *  @param nx side fraction used to divide the box in the x direction 
   *  @param ny side fraction used to divide the box in the y direction 
   *  @param nz side fraction used to divide the box in the z direction 
   *  @return none
   */
  void set_ObjectRegion_SubBoxes (catalogue::Catalogue &data, catalogue::Catalogue &random, const int nx, const int ny, const int nz);

  /**
   *  @brief set the object region in angular SubBoxes
   *  @param data input data catalogue
   *  @param random random catalogue
   *  @param Cell_size size of the cell in degrees
   *  @return none
   */
  void set_ObjectRegion_RaDec (catalogue::Catalogue &data, catalogue::Catalogue &random, const double Cell_size);

  /**
   *  @brief set the object region in sub-regions using mangle
   *  @param data input data catalogue
   *  @param random random catalogue
   *  @param nSamples number of sub-regions
   *  @param polygonfile name of the input file with polygons
   *  @return none
   */
  void set_ObjectRegion_mangle (catalogue::Catalogue &data, catalogue::Catalogue &random, const int nSamples, const std::string polygonfile);

  /**
   *  @brief set the object region in SDSS stripes
   *  @param data input data catalogue
   *  @param random random catalogue
   *  @return none
   */
  void set_ObjectRegion_SDSS_stripes (catalogue::Catalogue &data, catalogue::Catalogue &random);

  /**
   *  @brief check if the subdivision process produced the correct results
   *  @param data input data catalogue
   *  @param random random catalogue
   *  @return none
   */
  void check_regions (catalogue::Catalogue &data, catalogue::Catalogue &random);

  ///@}

  /**
   *  @name Generic functions for density field reconstruction
   */

  ///@{

  /**
   * @brief compute the non linear displacements
   * of the density field
   *
   * @param data the data catalogue
   *
   * @param random the random catalogue
   *
   * @param random_RSD true \f$ \rightarrow \f$ RSD displacements for the 
   * random sample; false no-RSD displacements for the random sample
   *
   * @param cosmology the cosmology
   *
   * @param redshift the redshift
   *
   * @param bias the bias
   *
   * @param cell_size the cell size for density field
   * computation
   *
   * @param smoothing_radius the smoothing radius for density field
   * computation
   *
   * @param interpolation_type 0 \f$ \rightarrow \f$ compute density 
   * field using Nearest Grid Point (NGP) method; 0 \f$ \rightarrow \f$ 
   * compute density field using Cloud-in-cell (CIC) method;
   *
   * @return none
   */
  void reconstruction_fourier_space(const catalogue::Catalogue data, const catalogue::Catalogue random, const bool random_RSD, const cosmology::Cosmology cosmology, const double redshift, const double bias, const double cell_size, const double smoothing_radius, const int interpolation_type=0);

  /**
   * @brief return a sample with objects displaced, according to the
   * internal variables m_x_displacement, m_y_displacement, m_z_displacement 
   * @param input_catalogue input catalogue
   * @return the displaced catalogue
   */
  catalogue::Catalogue displaced_catalogue (const catalogue::Catalogue input_catalogue);

  ///@}

  /**
   *  @name Functions to provide 2PCF/3PCF forecasts
   */

  ///@{

  /**
   * @brief fit the input covariance matrix using
   * the gaussian model, varying the number of objects
   * and the volume of the sample
   *
   * @param mean the 2PCF mean value
   *
   * @param mock_xi0 the 2CPF of the mocks
   *
   * @param doJK 0 &rarr; normalize to 1/(n-1); 1 &rarr; normalize
   *  to n-1/n (for Jackknife)
   *
   * @param cosmology the cosmology
   *
   * @param nObjects the number of objects in the sample
   *
   * @param Volume the volume of the sample
   *
   * @param bias the bias of the sample
   *
   * @param redshift the redshift of the sample
   *
   * @param rMin the minimum scale
   *
   * @param rMax the maximum scale
   *
   * @param nbins the number of bins
   *
   * @param bin_type the binning type
   *
   * @param method_Pk method used to compute the power spectrum;
   * valid choices for method_Pk are: CAMB [http://camb.info/],
   * classgal_v1 [http://class-code.net/], MPTbreeze-v1
   * [http://arxiv.org/abs/1207.1465], EisensteinHu
   * [http://background.uchicago.edu/~whu/transfer/transferpage.html]
   *
   * @param sigma_NL the BAO damping parameter
   *
   * @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1
   * \f$\rightarrow\f$ non-linear power spectrum
   *
   * @return the best-fit values 
   */
  std::vector<double> fit_covariance_matrix_2PCF_monopole (const std::vector<double> mean, const std::vector<std::vector<double>> mock_xi0, const bool doJK, const cbl::cosmology::Cosmology cosmology, const double nObjects, const double Volume, const double bias, const double redshift, const double rMin, const double rMax, const int nbins, const cbl::BinType bin_type, const std::string method_Pk="CAMB", const double sigma_NL=0., const bool NL=true);

  /**
   * @brief generate mock measurementes 
   * of the 2PCF monopole from gaussian
   * covariance matrix
   *
   * @param cosmology the cosmology
   *
   * @param bias the bias of the sample
   *
   * @param nObjects the number of objects in the sample
   *
   * @param Volume the volume of the sample
   *
   * @param redshift the redshift of the sample
   *
   * @param rMin the minimum scale
   *
   * @param rMax the maximum scale
   *
   * @param nbins the number of bins
   *
   * @param bin_type the binning type
   *
   * @param method_Pk method used to compute the power spectrum;
   * valid choices for method_Pk are: CAMB [http://camb.info/],
   * classgal_v1 [http://class-code.net/], MPTbreeze-v1
   * [http://arxiv.org/abs/1207.1465], EisensteinHu
   * [http://background.uchicago.edu/~whu/transfer/transferpage.html]
   *
   * @param sigma_NL the BAO damping parameter
   *
   * @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1
   * \f$\rightarrow\f$ non-linear power spectrum
   *
   * @return the best-fit values 
   */
  std::shared_ptr<cbl::data::Data> generate_mock_2PCF_monopole (const cbl::cosmology::Cosmology cosmology, const double bias, const double nObjects, const double Volume, const double redshift, const double rMin, const double rMax, const int nbins, const cbl::BinType bin_type, const std::string method_Pk="CAMB", const double sigma_NL=0., const bool NL=true);

  /**
   * @brief generate mock measurementes 
   * of the 2PCF multipoles from gaussian
   * covariance matrix
   *
   * @param cosmology the cosmology
   *
   * @param bias the bias of the sample
   *
   * @param nObjects the number of objects in the sample
   *
   * @param Volume the volume of the sample
   *
   * @param redshift the redshift of the sample
   *
   * @param rMin the minimum scale
   *
   * @param rMax the maximum scale
   *
   * @param nbins the number of bins
   *
   * @param bin_type the binning type
   *
   * @param method_Pk method used to compute the power spectrum;
   * valid choices for method_Pk are: CAMB [http://camb.info/],
   * classgal_v1 [http://class-code.net/], MPTbreeze-v1
   * [http://arxiv.org/abs/1207.1465], EisensteinHu
   * [http://background.uchicago.edu/~whu/transfer/transferpage.html]
   *
   * @param sigma_NL the BAO damping parameter
   *
   * @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1
   * \f$\rightarrow\f$ non-linear power spectrum
   *
   * @return the best-fit values 
   */
  std::shared_ptr<cbl::data::Data> generate_mock_2PCF_multipoles (const cbl::cosmology::Cosmology cosmology, const double bias, const double nObjects, const double Volume, const double redshift, const double rMin, const double rMax, const int nbins, const cbl::BinType bin_type, const std::string method_Pk="CAMB", const double sigma_NL=0., const bool NL=true);

  ///@}

  namespace glob {

    /**
     *  @class spherical_harmonics_coeff GlobalFunc.h "Headers/GlobalFunc.h"
     *
     *  @brief The class spherical_harmonics_coeff
     *
     *  This class is used to handle objects of type 
     *  <EM> spherical_harmonics_coeff </EM>. 
     *  It is used to perform the Slepian&Eisenstein algorithm to
     *  compute the three-point correlation function
     *  It contains all methods to compute spherical harmonics expansion \f$a_{lm}\f$ 
     *  of a collection of points in separation bins, and the product of \f$a_{lm}\f$
     *  for two separation bins.
     */
    class spherical_harmonics_coeff {

      protected:

	/// the number of separation bins
	int m_nbins;

	/// the number of multipoles, \f$ l_{max}+1 \f$
	int m_norder;

	/// the maximum multipole \f$ l_{max} \f$
	int m_lmax;

	/// the total number of spherical harmonics
	int m_n_sph;

	/// the number of spherical harmonics for a given choice of \f$l\f$
	std::vector<int> m_n_sph_l;

	/// the spherical harmonics expansion coefficients in separation bins
	std::vector<std::vector<std::complex<double>>> m_alm;

      public:
	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 * @brief default constructor
	 *
	 * @return object of type spherical_harmonics_coeff
	 */
	spherical_harmonics_coeff () {}

	/**
	 * @brief default constructor
	 *
	 * @param norder the number of multipoles, \f$ l_{max}+1 \f$
	 *
	 * @param nbins the number of separation bins
	 *
	 * @return object of type spherical_harmonics_coeff
	 */
	spherical_harmonics_coeff (const int norder, const int nbins=1) { initialize(norder, nbins); }

	/**
	 * @brief default descructor
	 *
	 * @return none
	 */
	~spherical_harmonics_coeff () {}

	///@}

	/**
	 * @brief return the real part of the n-th coefficient
	 * of the expansion for a given separation bin
	 *
	 * @param n the n-th coefficient of the spherical harmonics
	 * expansion
	 *
	 * @param bin the separation bin
	 *
	 * @return  the real part of the n-th coefficient
	 * of the expansion for a given separation bin
	 */
	double real (const int n, const int bin=0) { return m_alm[bin][n].real();} 

	/**
	 * @brief return the imaginary part of the n-th coefficient
	 * of the expansion for a given separation bin
	 *
	 * @param n the n-th coefficient of the spherical harmonics
	 * expansion
	 *
	 * @param bin the separation bin
	 *
	 * @return  the imaginary part of the n-th coefficient
	 * of the expansion for a given separation bin
	 */
	double imag (const int n, const int bin=0) { return m_alm[bin][n].imag();} 

	/**
	 * @brief initialize the internal quantities
	 *
	 * @param norder the number of multipoles, \f$ l_{max}+1 \f$
	 *
	 * @param nbins the number of separation bins 
	 *
	 * @return none
	 */
	void initialize (const int norder, const int nbins=1);

	/**
	 * @brief reset the internal quantities
	 *
	 * @return none
	 */
	void reset ();

	/**
	 * @brief compute the \f$ a_{lm}\f$ for the normalized
	 * coordinates \f$ \lbrace x, y, z \rbrace \f$
	 *
	 * @param xx the x coordinate
	 *
	 * @param yy the y coordinate
	 *
	 * @param zz the z coordinate
	 *
	 * @return vector containing the \f$ a_{lm}\f$ for the 
	 * normalized coordinates \f$ \lbrace x, y, z \rbrace \f$
	 */
	std::vector<std::complex<double>> alm(const double xx, const double yy, const double zz);

	/**
	 * @brief add the \f$ a_{lm}\f$ to a specific separation
	 * bin with a weight
	 *
	 * @param alm the spherical harmonics expansion coefficients
	 *
	 * @param ww the weight
	 *
	 * @param bin the separation bin
	 *
	 * @return none
	 */
	void add (const std::vector<std::complex<double>> alm, const double ww, const int bin=0);

	/**
	 * @brief compute add the \f$ a_{lm}\f$ for the normalized
	 * coordinates \f$ \lbrace x, y, z \rbrace \f$ to a specific 
	 * separation bin with a weight
	 *
	 * @param xx the x coordinate
	 *
	 * @param yy the y coordinate
	 *
	 * @param zz the z coordinate
	 *
	 * @param ww the weight
	 *
	 * @param bin the separation bin
	 *
	 * @return none
	 */
	void add (const double xx, const double yy, const double zz, const double ww, const int bin=0);

	/**
	 * @brief compute the product of the 
	 * \f$ a_{lm} \f$ in two separation bin.
	 *
	 * This function computes the product of the 
	 * \f$ a_{lm} \f$ in two separation bin:
	 *
	 * \f[
	 *   \zeta_l(r_1, r_2) = \sum_{m=-l}^l a_{lm}(r_1) a^*_{lm} (r_2)
	 * \f]
	 *
	 * @param l the coefficient order
	 *
	 * @param bin1 the first separation bin
	 *
	 * @param bin2 the second separation bin
	 *
	 * @return none
	 */
	double power (const int l, const int bin1, const int bin2);
    };


    /**
     * @brief compute the triples using
     * Slepian, Eisenstein 2015 approach
     *
     * @param r12_min the minimum separation of the first shell
     *
     * @param r12_max the maximum separation of the first shell
     * 
     * @param r13_min the minimum separation of the second shell
     * 
     * @param r13_max the maximum separation of the second shell
     * 
     * @param norders the number of multipoles, \f$ l_{max}+1 \f$
     * 
     * @param catalogue the catalogue
     *
     * @return the triplets
     */
    std::vector<double> count_triplets_SphericalHarmonics (const double r12_min, const double r12_max, const double r13_min, const double r13_max, const int norders, const catalogue::Catalogue catalogue);

    /**
     * @brief compute the triples using
     * Slepian, Eisenstein 2015 approach
     *
     * @param pairs the pairs in the separation bins
     *
     * @param triplets the multipoles expansion of the triplets
     * for the separation bins
     *
     * @param rmin the minimum separation
     *
     * @param rmax the maximum separation
     * 
     * @param nbins the number of separation bins 
     * 
     * @param norders the number of multipoles, \f$ l_{max}+1 \f$
     * 
     * @param catalogue the catalogue
     *
     * @return none
     */
    void count_triplets_SphericalHarmonics (std::vector<double> &pairs, std::vector<std::vector<std::vector<double>>> &triplets, const double rmin, const double rmax, const int nbins, const int norders, const catalogue::Catalogue catalogue);

    /**
     * @brief compute the reconstructed triplets from
     * the multipoles expansion obtained from
     * Slepian, Eisenstein 2015 approach in 
     * cbl::count_triplets_SphericalHarmonics
     *
     * @param rmin the minimum separation
     *
     * @param rmax the maximum separation
     * 
     * @param nbins the number of separation bins 
     * 
     * @param norders the number of multipoles, \f$ l_{max}+1 \f$
     * 
     * @param catalogue the catalogue
     *
     * @param output_dir the output directory
     *
     * @param output_file_pairs the file to store pairs
     *
     * @param output_file_triplets the file to store triplets
     *
     * @return none
     */
    void count_triplets_SphericalHarmonics (const double rmin, const double rmax, const int nbins, const int norders, const catalogue::Catalogue catalogue, const std::string output_dir, const std::string output_file_pairs, const std::string output_file_triplets);

    /**
     * @brief compute the reconstructed triplets from
     * the multipoles expansion obtained from
     * Slepian, Eisenstein 2015 approach in 
     * cbl::count_triplets_SphericalHarmonics
     * 
     * @param nbins the number of \f$ \theta \f$ bins to reconstruct
     * triplet counts
     *
     * @param side_s the size of r<SUB>12</SUB>
     *
     * @param side_u the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
     * 
     * @param perc_increase define the shell size as side_s*perc_increase
     * and side_s*side_u*perc_increase
     *
     * @param norders the number of multipoles, \f$ l_{max}+1 \f$
     * 
     * @param catalogue the catalogue
     *
     * @param random_catalogue the random catalogue
     *
     * @param output_dir output directory
     *
     * @param output_file output_file to write three-point
     * correlation function
     *
     * @return the three-point correlation function
     */
    std::vector<double> zeta_SphericalHarmonics (const int nbins, const double side_s, const double side_u, const double perc_increase, const int norders, const catalogue::Catalogue catalogue, const catalogue::Catalogue random_catalogue, const std::string output_dir, const std::string output_file);

    /**
     * @brief compute the reconstructed triplets from
     * the multipoles expansion obtained from
     * Slepian, Eisenstein 2015 approach in 
     * cbl::count_triplets_SphericalHarmonics
     * 
     * @param nbins the number of \f$ \theta \f$ bins to reconstruct
     * triplet counts
     *
     * @param r12_min the minimum separation of the first shell
     *
     * @param r12_max the maximum separation of the first shell
     * 
     * @param r13_min the minimum separation of the second shell
     * 
     * @param r13_max the maximum separation of the second shell
     * 
     * @param norders the number of multipoles, \f$ l_{max}+1 \f$
     * 
     * @param catalogue the catalogue
     *
     * @param random_catalogue the random catalogue
     *
     * @param output_dir output directory
     *
     * @param output_file output_file to write three-point
     * correlation function
     *
     * @return the three-point correlation function
     */
    std::vector<double> zeta_SphericalHarmonics (const int nbins, const double r12_min, const double r12_max, const double r13_min, const double r13_max, const int norders, const catalogue::Catalogue catalogue, const catalogue::Catalogue random_catalogue, const std::string output_dir, const std::string output_file);

    /**
     * @brief compute the triples using
     * Slepian, Eisenstein 2015 approach
     *
     * This function computes the multipoles expansion
     * of the binned triplet counts from a the joined data+random 
     * catalogue. The random catalogue objects are weighted by
     * \f[ \hat{w}_i = -\frac{N_D}{N_R} w_i; \f]
     *
     * with \f$ N_D, N_R \f$ the weighted number of data and random objects
     * and \f$ w_i \f$ the original weight of the random object.
     *
     * The multipoles expansion of triplets in the joined catalogue corresponds
     * to the numerator of the Szapudi-Szalay estimator of the three-point
     * correlation function.
     *
     * By looking at negative weights we can compute the multipoles
     * expansion 
     *
     * @param NNN the multipoles expansion of the joined catalogue
     * triplets
     *
     * @param RRR the multipoles expansion of the random catalogue
     * triplets
     *
     * @param r12_min the minimum separation of the first shell
     *
     * @param r12_max the maximum separation of the first shell
     * 
     * @param r13_min the minimum separation of the second shell
     * 
     * @param r13_max the maximum separation of the second shell
     * 
     * @param norders the number of multipoles, \f$ l_{max}+1 \f$
     * 
     * @param catalogue the catalogue
     *
     * @return the triplets
     */
    void  count_triplets_SphericalHarmonics (std::vector<double> &NNN, std::vector<double> &RRR, const double r12_min, const double r12_max, const double r13_min, const double r13_max, const int norders, const catalogue::Catalogue catalogue);

    /**
     * @brief compute the reconstructed triplets from
     * the multipoles expansion obtained from
     * Slepian, Eisenstein 2015 approach in 
     * cbl::count_triplets_SphericalHarmonics
     * 
     * @param nbins the number of \f$ \theta \f$ bins to reconstruct
     * triplet counts
     *
     * @param side_s the size of r<SUB>12</SUB>
     *
     * @param side_u the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
     * 
     * @param perc_increase define the shell size as side_s*perc_increase
     * and side_s*side_u*perc_increase
     *
     * @param norders the number of multipoles, \f$ l_{max}+1 \f$
     * 
     * @param catalogue the catalogue
     *
     * @param random_catalogue the random catalogue
     *
     * @param output_dir output directory
     *
     * @param output_file output_file to write three-point
     * correlation function
     *
     * @param count_triplets  1 &rarr; count the multipoles
     * expansion of the triplets 0 &rarr; read the multipoles
     * expansion of the triplets
     *
     * @param dir_triplets the triplets directory
     *
     * @return the three-point correlation function
     */
    std::vector<double> zeta_SphericalHarmonics_AllInOne (const int nbins, const double side_s, const double side_u, const double perc_increase, const int norders, const catalogue::Catalogue catalogue, const catalogue::Catalogue random_catalogue, const std::string output_dir, const std::string output_file, const bool count_triplets=true, const std::string dir_triplets=cbl::par::defaultString);

    /**
     * @brief compute the reconstructed triplets from
     * the multipoles expansion obtained from
     * Slepian, Eisenstein 2015 approach in 
     * cbl::count_triplets_SphericalHarmonics
     * 
     * @param nbins the number of \f$ \theta \f$ bins to reconstruct
     * triplet counts
     *
     * @param r12_min the minimum separation of the first shell
     *
     * @param r12_max the maximum separation of the first shell
     * 
     * @param r13_min the minimum separation of the second shell
     * 
     * @param r13_max the maximum separation of the second shell
     *
     * @param norders the number of multipoles, \f$ l_{max}+1 \f$
     * 
     * @param catalogue the catalogue
     *
     * @param random_catalogue the random catalogue
     *
     * @param output_dir output directory
     *
     * @param output_file output_file to write three-point
     * correlation function
     *
     * @param count_triplets  1 &rarr; count the multipoles
     * expansion of the triplets 0 &rarr; read the multipoles
     * expansion of the triplets
     *
     * @param dir_triplets the triplets directory
     *
     * @return the three-point correlation function
     */
    std::vector<double> zeta_SphericalHarmonics_AllInOne (const int nbins, const double r12_min, const double r12_max, const double r13_min, const double r13_max, const int norders, const catalogue::Catalogue catalogue, const catalogue::Catalogue random_catalogue, const std::string output_dir, const std::string output_file, const bool count_triplets=true, const std::string dir_triplets=cbl::par::defaultString);

    /**
     * @brief compute binned triplets using
     * the direct triplet-finding method
     *
     * @param r12_min the minimum separation of the first shell
     *
     * @param r12_max the maximum separation of the first shell
     * 
     * @param r13_min the minimum separation of the second shell
     * 
     * @param r13_max the maximum separation of the second shell
     *
     * @param nbins the number of bins to bin triplet
     *
     * @param catalogue the catalogue
     *
     * @param tripletType the triplet type
     *
     * @return the binned triplets counts with 
     * two sides fixed
     */
    std::vector<double> count_triplets_classic (const double r12_min, const double r12_max, const double r13_min, const double r13_max, const int nbins, const cbl::catalogue::Catalogue catalogue, cbl::triplets::TripletType tripletType);

    /**
     * @brief compute edge-correction
     * and return legendre coefficients of the 3PCF
     *
     * @param NNN the Data-Random triplet legendre coefficients
     * @param RRR the random triplet legendre coefficients
     * @param normalization the normalization factor: 	
     *        \f$ n_R*(n_R-1)*(n_R-2)/(n_G*(n_G-1)*(n_G-2))\f$
     * @return the legendre coefficients of the 3PCF
     */
    std::vector<double> zeta_SphericalHarmonics_edgeCorrection (const std::vector<double> NNN, const std::vector<double> RRR, const double normalization=1.);
  }

}


#endif
