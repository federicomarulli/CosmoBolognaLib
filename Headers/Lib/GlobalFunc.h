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
 *  @file Headers/Lib/GlobalFunc.h
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


namespace cosmobl {

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
   *  @param [in] idum the random seed
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
   *  @return none
   */
  void Vmax_DC_distribution (vector<double> &, vector<double> &, const vector<double>, const vector<double>, const vector<double>, const double, const double, const double, const double, Cosmology &, const double, const int, const int, const bool norm=1, const string file_Vmax=par::defaultString, const double delta_dc_Vmax=100.);

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
  double AP_shift_r (const double, const Cosmology &, const Cosmology &);

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
  double AP_shift_rp (const double, const Cosmology &, const Cosmology &);

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
  double AP_shift_pi (const double, const Cosmology &, const Cosmology &);

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
  void max_separations_AP (const double, const double, const double, const Cosmology &, const vector<Cosmology> &, double &, double &, double &);

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
   *  @param direction 0 &rarr; cosm2 &rArr; cosm1; 1 &rarr; cosm1
   *  &rArr; cosm2;
   *
   *  @return the converted two-point correlation function, &xi;(R)
   */
  double converted_xi (const double, const double, const vector<double>, const vector<double>, const Cosmology &, const Cosmology &, const bool);

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
   *  @param redshift
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
   *  @param direction 0 &rarr; cosm2 &rArr; cosm1; 1 &rarr; cosm1
   *  &rArr; cosm2;
   *
   *  @return the converted two-point correlation function, &xi;(R<SUB>p</SUB>,&Pi;)
   */
  double converted_xi (const double, const double, const double, const vector<double>, const vector<double>, const vector<vector<double> >, const Cosmology &, const Cosmology &, const bool); 

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
  void redshift_range (const double, const double, Cosmology &, double &, double &); 

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
  double volume (const double, const int, const double, const double, Cosmology &);

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
   *  @param[in] idum the random seed
   *
   *  @return none
   */
  void coord_zSpace (vector<double> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &, const vector<double>, const vector<double>, const vector<double>, const double, Cosmology &, const double, const double, const double, const int);

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
  void create_mocks (const vector<double>, const vector<double>, const vector<double>, const vector<double>, const vector<double>, const vector<double>, const vector<double>, const vector<double>, const vector<double>, const string, const double, const int, const double, const double, Cosmology &, const int, const double, const int, double &);

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
   *  @param dir output directories
   *  @return none
   */
  void set_ObjectRegion_mangle (catalogue::Catalogue &data, catalogue::Catalogue &random, const int nSamples, const string polygonfile, const string dir);

  /**
   *  @brief check if the subdivision process produced the correct results
   *  @param data input data catalogue
   *  @param random random catalogue
   *  @return none
   */
  void check_regions (catalogue::Catalogue &data, catalogue::Catalogue &random);

  ///@}
}


#endif
