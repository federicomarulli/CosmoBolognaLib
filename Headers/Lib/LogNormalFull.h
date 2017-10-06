/*******************************************************************
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
 *  @file Headers/Lib/LogNormalFull.h
 *
 *  @brief Implementation of the lognormal data structure
 *
 *  This file defines the interface of the class LogNormalFull
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __LOGNORMALF__
#define __LOGNORMALF__

#include "Catalogue.h"
#include <fftw3.h>

namespace cosmobl {
    
  /**
   *  @brief The namespace of the functions and classes used to
   *  construct <B> log-normal mocks </B> 
   *  
   *  The \e lognormal namespace contains all the functions and
   *  classes used to construct log-normal mocks
   */
  namespace lognormal {

    /**
     *  @class LogNormalFull LogNormalFull.h
     *  "Headers/Lib/LogNormalFull.h"
     *
     *  @brief The class LogNormalFull
     *
     *  This class is used to handle objects of type <EM> lognormalFull
     *  </EM>
     */
    class LogNormalFull {

    protected:

      /// pointer to the fiducial cosmology
      shared_ptr<cosmology::Cosmology> m_cosmology;

      /// string containing the author for the model power spectrum
      string m_author;

      /// 0 \f$\rightarrow\f$ compute the linear power spectrum; 1 \f$\rightarrow\f$ compute the non-linear power spectrum
      bool m_NL;

      /// generator for random numbers
      default_random_engine m_generator;

      /// the clustering signal 
      shared_ptr<data::Field3D> m_clustering_signal;

      /// approximate cell size of the density field
      double m_rmin;

      /// the number of cells along x-axis
      int m_nx;

      /// the number of cells along y-axis
      int m_ny;

      /// the number of cells along z-axis
      int m_nz;

      /// the number of cells along z-axis, Fourier Space
      int m_nzF;

      /// the minimum x-value
      double m_xMin;

      /// the maximum x-value
      double m_xMax;

      /// the minimum y-value
      double m_yMin;

      /// the maximum y-value
      double m_yMax;

      /// the minimum z-value
      double m_zMin;

      /// the maximum z-value
      double m_zMax;

      /// the effective cell numbers
      int m_nCells_eff;

      /// the number of lognormal realization to produce/read
      int m_nLN;

      /// smoothing radius
      double m_rsmooth;

      /// 1 \f$\rightarrow\f$ use visibility, 0 \f$\rightarrow\f$ don't use visibility
      bool m_use_visibility;

      /// pointer to the input datasets
      vector<shared_ptr<catalogue::Catalogue>> m_data;

      /// pointer to the random catalogues
      vector<shared_ptr<catalogue::Catalogue>> m_random;

      /// vector containing pointers to the LogNormal realizations
      shared_ptr<data::Field3D> m_density;  

      /// vector containing pointers to the LogNormal realizations
      shared_ptr<data::Field3D> m_densityG;  

      /// density field variance
      double m_sigma2G;

      /// vector containing pointers to the LogNormal realizations
      shared_ptr<data::Field3D> m_rsd_displacement;  

      /// vector containing pointers to the LogNormal realizations
      shared_ptr<data::Field3D> m_displacement; 

      /// vector containing pointers to the LogNormal realizations
      shared_ptr<data::Field3D> m_potential;  

      /// vector containing pointers to the LogNormal realizations
      shared_ptr<data::Field3D> m_velocity; 

      /// vector containing pointers to the LogNormal realizations
      shared_ptr<data::Field3D> m_los_velocity; 

      /// vector containing pointers to the LogNormal realizations
      shared_ptr<data::Field3D> m_visibility; 

      /// vector containing pointers to the LogNormal realizations
      vector<shared_ptr<data::Field3D>> m_visibility_random; 

      /// vector containing pointers to the LogNormal realizations
      vector<vector<shared_ptr<catalogue::Catalogue>>> m_LNCat;  

      /// interpolator \f$z(D_c)\f$ 
      shared_ptr<glob::FuncGrid> m_func_redshift;

      /// interpolator \f$D_c(z)\f$ 
      shared_ptr<glob::FuncGrid> m_func_DC;

      /// interpolator \f$H(z)\f$ 
      shared_ptr<glob::FuncGrid> m_func_HH;

      /// interpolator \f$b(z)\f$ 
      vector<shared_ptr<glob::FuncGrid> > m_func_bias;

      /// interpolator \f$gg(z)\f$ 
      shared_ptr<glob::FuncGrid> m_func_growth_factor;

      /// interpolator \f$f(z)\f$ 
      shared_ptr<glob::FuncGrid> m_func_growth_rate;

      /// interpolator \f$P(k)\f$ 
      shared_ptr<glob::FuncGrid> m_func_pk;

      /**
       * @brief set the grid parameters
       * @return none
       */
      void set_grid_parameters ();

      /**
       * @brief set the fields
       *
       *  @param use_random true \f$\rightarrow\f$ use random for
       *  visibility mask; false \f$\rightarrow\f$ do not use random
       *
       *  @param doRSD true \f$\rightarrow\f$do redshift-space
       *  mocks; false \f$\rightarrow\f$ do real-space mocks
       *
       *  @return none
       */
      void set_fields (const bool use_random, const bool doRSD);

      /**
       *  @brief reset the fields
       *
       *  @param use_random true \f$\rightarrow\f$ use random for
       *  visibility mask, false \f$\rightarrow\f$ do not use random
       *
       *  @param doRSD true \f$\rightarrow\f$do redshift-space
       *  mocks, false \f$\rightarrow\f$ do real-space mocks
       *
       *  @return none
       */
      void reset_fields (const bool use_random, const bool doRSD);

      /**
       *  @brief set the target dark matter clustering signal
       *
       *  @return none
       */
      void set_clustering_signal ();

      /**
       *  @brief set the density field
       *
       *  @param smoothing_radius the gaussian kernel size
       *
       *  @return none
       */
      void set_density_field (const double smoothing_radius);

      /**
       *  @brief set the gravitational potential field
       *  @return none
       */
      void set_potential ();

      /**
       *  @brief set the radial velocity field
       *  @return none
       */
      void set_radial_velocity ();

      /**
       *  @brief set the visibility
       *  @return none
       */
      void set_visibility ();

      /**
       *  @brief set the visibility from random catalogue
       *
       *  @return none
       */
      void set_visibility_from_random ();

      /**
       *  @brief set the visibility for redshift space mocks
       *
       *  @return none
       */
      void set_visibility_from_random_RSD ();

      /**
       *  @brief extract points from lognormal density fields
       *
       *  @param nObjects number of generated points
       *
       *  @param doRSD true \f$\rightarrow\f$do redshift-space
       *  mocks, false \f$\rightarrow\f$ do real-space mocks
       *
       *  @param redshift vector containing the redshift
       *
       *  @param bias vector containing the bias in function of redshift
       *
       *  @param visibility pointer to the visibility
       *
       *  @param file_out out file for the extracted sample
       *
       *  @return none
       */
      void extract_points_lognormal_field (const double nObjects, const bool doRSD, const vector<double> redshift, const vector<double> bias, const shared_ptr<data::Field3D> visibility, const string file_out);

	
    public:

      /**
       *  @brief default constructor
       *  @return object of class LogNormal
       */
      LogNormalFull () = default;

      /**
       *  @brief constructor
       *  @param cosmology the input cosmology
       *  @param redshift_min the minimum redshift
       *  @param redshift_max the maximum redshift
       *  @param nredshift the number or redshift bins
       *  @param author the linear power spectrum method
       *  @return object of class LogNormal
       */
      LogNormalFull (const cosmology::Cosmology cosmology, const double redshift_min=0., const double redshift_max=10., const int nredshift=500, const string author="CAMB");

      /**
       *  @brief constructor
       *  @param rmin cell size
       *  @param xMin the minimum x-value
       *  @param xMax the minimum x-value
       *  @param yMin the maximum y-value
       *  @param yMax the maximum x-value
       *  @param zMin the minimum z-value
       *  @param zMax the maximum z-value
       *  @param cosmology the input cosmology
       *  @param redshift_min the minimum redshift
       *  @param redshift_max the maximum redshift
       *  @param nredshift the number or redshift bins
       *  @param author the linear power spectrum
       *  @return object of class LogNormal
       */
      LogNormalFull (const double rmin, const double xMin, const double xMax, const double yMin, const double yMax, const double zMin, const double zMax, const cosmology::Cosmology cosmology, const double redshift_min=0., const double redshift_max=10., const int nredshift=500, const string author="CAMB");

      /**
       *  @brief constructor
       *  @param rmin cell size
       *  @param random vector containing random samples
       *  @param pad the size of padding area around grid
       *  @param cosmology the input cosmology
       *  @param redshift_min the minimum redshift
       *  @param redshift_max the maximum redshift
       *  @param nredshift the number or redshift bins
       *  @param author the linear power spectrum
       *  @return object of class LogNormal
       */
      LogNormalFull (const double rmin, const vector<shared_ptr<catalogue::Catalogue>> random, const double pad, const cosmology::Cosmology cosmology, const double redshift_min=0., const double redshift_max=10., const int nredshift=500, const string author="CAMB");

      /**
       *  @brief default destructor
       *  @return none
       */
      ~LogNormalFull () = default;

      /**
       *  @brief set cosmological functions
       *  @param cosmology the input cosmology
       *  @param redshift_min the minimum redshift
       *  @param redshift_max the maximum redshift
       *  @param nredshift the number or redshift bins
       *  @param author the linear power spectrum
       *  @return none
       */
      void set_cosmo_function (const cosmology::Cosmology cosmology, const double redshift_min=0., const double redshift_max=10., const int nredshift=500, const string author="CAMB");

      /**
       *  @brief set grid parameters
       *
       *  @param rmin cell size
       *  @param xMin the minimum x-value
       *  @param xMax the minimum x-value
       *  @param yMin the maximum y-value
       *  @param yMax the maximum x-value
       *  @param zMin the minimum z-value
       *  @param zMax the maximum z-value
       *
       *  @return none
       */
      void set_grid_parameters (const double rmin, const double xMin, const double xMax, const double yMin, const double yMax, const double zMin, const double zMax);

      /**
       *  @brief set grid parameters
       *  @param rmin cell size
       *  @param random vector containing random samples
       *  @param pad the size of padding area around grid
       *  @return none
       */
      void set_grid_parameters (const double rmin, const vector<shared_ptr<catalogue::Catalogue>> random, const double pad);

      /**
       *  @brief generate lognormal realizations
       *         
       *  @param start index of the first realization
       *
       *  @param stop index of the last realization
       *
       *  @param doRSD true \f$\rightarrow\f$ do redshift-space mocks,
       *  false \f$\rightarrow\f$ do real-space mocks
       *
       *  @param smoothing_radius the smoothing radius
       *
       *  @param nObjects number of generated points
       *
       *  @param redshift vector containing the redshift
       *
       *  @param bias vector containing the bias in function of
       *  redshift
       *
       *  @param dir output directory
       *
       *  @param filename output filename
       *
       *  @param seed seed for random generator
       *
       *  @param setfields true \f$\rightarrow\f$ set the fields,
       *  false \f$\rightarrow\f$ don't set the fields
       *
       *  @param use_random true \f$\rightarrow\f$ use random for
       *  visibility mask, false \f$\rightarrow\f$ don't use random
       *  for visibility mask
       *
       *  @return none
       */
      void generate_lognormal (const int start, const int stop, const bool doRSD, const double smoothing_radius, const vector<double> nObjects, const vector<vector<double>> redshift, const vector<vector<double> > bias, const string dir, const string filename, const int seed, const bool setfields=1, const bool use_random=true);

    };
  }
}

#endif
