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

}


#endif
