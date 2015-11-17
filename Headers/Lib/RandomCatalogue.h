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
 *  @file Headers/Lib/RandomCatalogue.h
 *
 *  @brief Functions for random catalogues
 *
 *  This file contains the prototypes of a set of functions to create
 *  random catalogues
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __RANDOMCAT__
#define __RANDOMCAT__

#include "LogNormal.h"


// ============================================================================


namespace cosmobl {


  /**
   *  @name Functions to create or read random catalogues
   */
  ///@{

  /**
   *  @brief read a random catalogue
   *  @param file vector containing the files where the random
   *  catalogues are stored
   *  @param nSub the fracton of objects randomly selected (nSub=1
   *  &rArr; all objects are selected)
   *  @return an object of class Catalogue
   */
  shared_ptr<Catalogue> random_catalogue_fromFile (const vector<string>, const double nSub=1.1); 

  /**
   *  @brief read a random catalogue with polar coordinates [ra, dec,
   redshift]
   *  @param file_in the name of the input file
   *  @param z_min the minimum redshift of the catalogue
   *  @param z_max the maximum redshift of the catalogue
   *  @param cosm object of class Cosmology
   *  @param nSub the fracton of objects randomly selected (nSub=1
   *  &rArr; all objects are selected)
   *  @param fact conversion factor
   *  @return an object of class Catalogue
   */
  shared_ptr<Catalogue> random_catalogue_radecred_fromFile (const string, const double, const double, const Cosmology &, const double nSub=1.1, const double fact=1.); 
  
  /**
   *  @brief create a random catalogue in a box
   *  @param catalogue object of class Catalogue
   *  @param nRandom number of random objects
   *  @return an object of class Catalogue
   */
  shared_ptr<Catalogue> random_catalogue_box (const shared_ptr<Catalogue>, const int);

  /**
   *  @brief create a random catalogue in a box
   *
   *  this function reads a cubic random catalogue from a file,
   *  generated in a given cosmology, and trasforms it into a new one
   *  in a different cosmology
   *
   *  @param real_cosm object of class Cosmology representing the \e
   *  real (or \e assumed) cosmology
   *
   *  @param test_cosm object of class Cosmology representing the \e
   *  test cosmology
   *
   *  @param dir_in the input directory where the original random
   *  catalogue is stored
   *
   *  @param dir_out the output directory where the new random
   *  catalogue will be stored
   *     
   *  @param Zguess_min minimum redshift used to search the redshift
   *
   *  @param Zguess_max maximum redshift used to search the redshift
   *
   *  @return an object of class Catalogue
   */
  shared_ptr<Catalogue> warped_random_catalogue (const Cosmology &, const Cosmology &, const string, const string, const double, const double);

  /**
   *  @brief create a random catalogue in a cone
   *
   *  @param [in] catalogue object of class Catalogue
   *
   *  @param [in] nRandom the number of random objects
   *
   *  @param [in] cosm object of class Cosmology 
   *
   *  @param [in] Angle angle of the cone 
   *
   *  @param [in] step_redshift the number of steps in redshift used to
   *  redshift distribution of the random object; if step_redshift=0
   *  the redshift distribution is estimated from the convolvolution
   *  of N(D<SUB>C</SUB>)
   *
   *  @param [in] redshift vector containing the redshift of the object in
   *  the real catalogue
   *
   *  @param [out] dc vector containing the central values of the binned comoving distances
   *  of the random objects
   *
   *  @param [out] convol vector containing the central values of the
   *  binned smoothed distribution of comoving distances of the random
   *  objects
   *
   *  @param [in] idum the random seed
   *  @return an object of class Catalogue
   */
  shared_ptr<Catalogue> random_catalogue_cone (const shared_ptr<Catalogue>, const int, const Cosmology &, const double, const int, const vector<double>, vector<double> &, vector<double> &, const int idum=13);

  /**
   *  @brief create a random catalogue for a mock sample (with polar
   * coordinates [ra, dec, redshift])
   *  
   *  @param [in] catalogue object of class Catalogue
   *
   *  @param [in] nRandom the number of random objects
   *
   *  @param [in] cosm object of class Cosmology 
   *
   *  @param [in] dir the directory where the random catalogue is stored
   *
   *  @param [in] step_redshift the number of steps in redshift used to
   *  redshift distribution of the random object; if step_redshift=0
   *  the redshift distribution is estimated from the convolvolution
   *  of N(D<SUB>C</SUB>)
   *
   *  @param [in] redshift vector containing the redshift of the object in
   *  the real catalogue
   *
   *  @param [out] dc vector containing the central values of the binned comoving distances
   *  of the random objects
   *
   *  @param [out] convol vector containing the central values of the
   *  binned smoothed distribution of comoving distances of the random
   *  objects
   *
   *  @param [in] idum the random seed
   *  @return an object of class Catalogue
   */
  shared_ptr<Catalogue> random_catalogue_mock (shared_ptr<Catalogue>, int &, Cosmology &, string &, int &, vector<double>, vector<double> &, vector<double> &, int idum=13);
  
  /**
   *  @brief compute the redshift distribution of a random catalogue
   *
   *  @author Alfonso Veropalumbo
   *  @author alfonso.veropalumbo@unibo.it
   *
   *  @param random object of class Catalogue
   *
   *  @param data object of class Catalogue
   *
   *  @param dir_random the directory where the random catalogue is stored
   *
   *  @param nbin number of redshift bin
   *
   *  @param convolution 0 &rarr; don't convolve the distribution; 1
   *  &rarr; convolve the distribution with a gaussian function
   *
   *  @param sigma &sigma;: the standard deviation of the gaussian
   *  function used to convolve the distribution
   *  
   *  @return none
   */
  void random_redshift_distribution (const shared_ptr<Catalogue> random, const shared_ptr<Catalogue> data, const string dir_random, const int nbin, const bool convolution, const double sigma);

  /// @cond extrandom
  shared_ptr<Catalogue> random_catalogue_mock_cone (const shared_ptr<Catalogue>, const int, const Cosmology &, const string, const int, vector<double>);

  shared_ptr<Catalogue> random_catalogue_VIPERS (const int, const Cosmology &, const string, const int, const vector<double>, const vector<double>, const vector<double>, const vector<double>, const bool, string, const string, const string, const int idum=13); 

  shared_ptr<Catalogue> random_sdss_angular_distribution (const int, const string, const string, const bool veto=0); 
  /// @endcond
  
  ///@}

}

#endif
