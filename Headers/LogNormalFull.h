/********************************************************************
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
 *  @file Headers/LogNormalFull.h
 *
 *  @brief Implementation of the log-normal data structure
 *
 *  This file defines the interface of the class LogNormalFull
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __LOGNORMALF__
#define __LOGNORMALF__

#include "Catalogue.h"
#include <fftw3.h>

namespace cbl {
    
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
     *  "Headers/LogNormalFull.h"
     *
     *  @brief The class LogNormalFull
     *
     *  This class is used to handle objects of type <EM>
     *  LogNormalFull </EM>. It is an improved version of the
     *  LogNormal class, aimed at providing log-normal mocks in
     *  redshift-space with reliable clustering multipoles. Still to
     *  be tested!
     */
    class LogNormalFull {

    protected:
      
      /// pointer to the input datasets
      std::vector<std::shared_ptr<catalogue::Catalogue>> m_data;

      /// pointer to the random catalogues
      std::vector<std::shared_ptr<catalogue::Catalogue>> m_random;
      
      /// pointer to the fiducial cosmology
      std::shared_ptr<cosmology::Cosmology> m_cosmology;

      /// std::string containing the author of the model power spectrum (i.e. the Boltzmann solver)
      std::string m_author;

      /// false \f$\rightarrow\f$ compute the linear power spectrum; true \f$\rightarrow\f$ compute the non-linear power spectrum
      bool m_NL;

      /// generator for random numbers
      std::default_random_engine m_generator;

      /// the clustering signal 
      std::shared_ptr<data::Field3D> m_clustering_signal;

      /// approximate cell size of the density field
      double m_rmin;

      /// the number of cells along x-axis
      int m_nx;

      /// the number of cells along y-axis
      int m_ny;

      /// the number of cells along z-axis
      int m_nz;

      /// the number of cells along z-axis, Fourier space
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

      /// the number of log-normal realizations to produce/read
      int m_nLN;

      /// smoothing radius
      double m_rsmooth;

      /// 1 \f$\rightarrow\f$ use the visibility mask, 0 \f$\rightarrow\f$ don't use the visibility mask
      bool m_use_visibility;

      /// vector containing pointers to the LogNormal realizations
      std::shared_ptr<data::Field3D> m_density;  

      /// vector containing pointers to the LogNormal realizations
      std::shared_ptr<data::Field3D> m_densityG;  

      /// density field variance
      double m_sigma2G;

      /// vector containing pointers to the LogNormal realizations
      std::shared_ptr<data::Field3D> m_rsd_displacement;  

      /// vector containing pointers to the LogNormal realizations
      std::shared_ptr<data::Field3D> m_displacement; 

      /// vector containing pointers to the LogNormal realizations
      std::shared_ptr<data::Field3D> m_potential;  

      /// vector containing pointers to the LogNormal realizations
      std::shared_ptr<data::Field3D> m_velocity; 

      /// vector containing pointers to the LogNormal realizations
      std::shared_ptr<data::Field3D> m_los_velocity; 

      /// vector containing pointers to the LogNormal realizations
      std::shared_ptr<data::Field3D> m_visibility; 

      /// vector containing pointers to the LogNormal realizations
      std::vector<std::shared_ptr<data::Field3D>> m_visibility_random; 

      /// vector containing pointers to the LogNormal realizations
      std::vector<std::vector<std::shared_ptr<catalogue::Catalogue>>> m_LNCat;  

      /// interpolator for \f$z(D_c)\f$ 
      std::shared_ptr<glob::FuncGrid> m_func_redshift;

      /// interpolator for \f$D_c(z)\f$ 
      std::shared_ptr<glob::FuncGrid> m_func_DC;

      /// interpolator for \f$H(z)\f$ 
      std::shared_ptr<glob::FuncGrid> m_func_HH;

      /// interpolator for \f$b(z)\f$ 
      std::vector<std::shared_ptr<glob::FuncGrid> > m_func_bias;

      /// interpolator for \f$g(z)\f$ 
      std::shared_ptr<glob::FuncGrid> m_func_growth_factor;

      /// interpolator for \f$f(z)\f$ 
      std::shared_ptr<glob::FuncGrid> m_func_growth_rate;

      /// interpolator for \f$P(k)\f$ 
      std::shared_ptr<glob::FuncGrid> m_func_pk;

      /**
       * @brief set the grid parameters
       */
      void m_set_grid_parameters ();

      /**
       *  @brief set the fields
       *
       *  @param use_random true \f$\rightarrow\f$ use an input random
       *  catalogue for the visibility mask; false \f$\rightarrow\f$
       *  do not use random
       *
       *  @param doRSD true \f$\rightarrow\f$ construct redshift-space
       *  mocks; false \f$\rightarrow\f$ construct real-space mocks
       */
      void m_set_fields (const bool use_random, const bool doRSD);

      /**
       *  @brief reset the fields
       *
       *  @param use_random true \f$\rightarrow\f$ use an input random
       *  catalogue for visibility mask, false \f$\rightarrow\f$ do
       *  not use random
       *
       *  @param doRSD true \f$\rightarrow\f$do redshift-space
       *  mocks, false \f$\rightarrow\f$ do real-space mocks
       */
      void m_reset_fields (const bool use_random, const bool doRSD);

      /**
       *  @brief set the target dark matter clustering signal
       */
      void m_set_clustering_signal ();

      /**
       *  @brief set the density field
       *
       *  @param smoothing_radius the gaussian kernel size
       */
      void m_set_density_field (const double smoothing_radius);

      /**
       *  @brief set the gravitational potential field
       */
      void m_set_potential ();

      /**
       *  @brief set the radial velocity field
       */
      void m_set_radial_velocity ();

      /**
       *  @brief set the visibility
       */
      void m_set_visibility ();

      /**
       *  @brief set the visibility from random catalogue
       */
      void m_set_visibility_from_random ();

      /**
       *  @brief set the visibility for redshift space mocks
       */
      void m_set_visibility_from_random_RSD ();

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
       */
      void m_extract_points_lognormal_field (const double nObjects, const bool doRSD, const std::vector<double> redshift, const std::vector<double> bias, const std::shared_ptr<data::Field3D> visibility, const std::string file_out);

	
    public:

      /**
       *  @brief default constructor
       */
      LogNormalFull () = default;

      /**
       *  @brief constructor
       *  @param cosmology the input cosmology
       *  @param redshift_min the minimum redshift
       *  @param redshift_max the maximum redshift
       *  @param n_redshift_bins the number or redshift bins
       *  @param author the linear power spectrum method
       */
      LogNormalFull (const cosmology::Cosmology cosmology, const double redshift_min=0., const double redshift_max=10., const int n_redshift_bins=500, const std::string author="CAMB");

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
       *  @param n_redshift_bins the number or redshift bins
       *  @param author the linear power spectrum
       */
      LogNormalFull (const double rmin, const double xMin, const double xMax, const double yMin, const double yMax, const double zMin, const double zMax, const cosmology::Cosmology cosmology, const double redshift_min=0., const double redshift_max=10., const int n_redshift_bins=500, const std::string author="CAMB");

      /**
       *  @brief constructor
       *  @param rmin cell size
       *  @param random vector containing random samples
       *  @param pad the size of padding area around grid
       *  @param cosmology the input cosmology
       *  @param redshift_min the minimum redshift
       *  @param redshift_max the maximum redshift
       *  @param n_redshift_bins the number or redshift bins
       *  @param author the linear power spectrum
       */
      LogNormalFull (const double rmin, const std::vector<std::shared_ptr<catalogue::Catalogue>> random, const double pad, const cosmology::Cosmology cosmology, const double redshift_min=0., const double redshift_max=10., const int n_redshift_bins=500, const std::string author="CAMB");

      /**
       *  @brief default destructor
       */
      ~LogNormalFull () = default;

      /**
       *  @brief set cosmological functions
       *  @param cosmology the input cosmology
       *  @param redshift_min the minimum redshift
       *  @param redshift_max the maximum redshift
       *  @param n_redshift_bins the number or redshift bins
       *  @param author the linear power spectrum
       */
      void set_cosmo_function (const cosmology::Cosmology cosmology, const double redshift_min=0., const double redshift_max=10., const int n_redshift_bins=500, const std::string author="CAMB");

      /**
       *  @brief set grid parameters
       *  @param rmin cell size
       *  @param xMin the minimum x-value
       *  @param xMax the minimum x-value
       *  @param yMin the maximum y-value
       *  @param yMax the maximum x-value
       *  @param zMin the minimum z-value
       *  @param zMax the maximum z-value
       */
      void set_grid_parameters (const double rmin, const double xMin, const double xMax, const double yMin, const double yMax, const double zMin, const double zMax);

      /**
       *  @brief set grid parameters
       *  @param rmin cell size
       *  @param random vector containing random samples
       *  @param pad the size of padding area around grid
       */
      void set_grid_parameters (const double rmin, const std::vector<std::shared_ptr<catalogue::Catalogue>> random, const double pad);

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
       *  @param set_fields true \f$\rightarrow\f$ set the fields,
       *  false \f$\rightarrow\f$ don't set the fields
       *
       *  @param use_random true \f$\rightarrow\f$ use an input random
       *  catalogue for the visibility mask; false \f$\rightarrow\f$
       *  don't use random
       */
      void generate_lognormal (const int start, const int stop, const bool doRSD, const double smoothing_radius, const std::vector<double> nObjects, const std::vector<std::vector<double>> redshift, const std::vector<std::vector<double> > bias, const std::string dir, const std::string filename, const int seed, const bool set_fields=true, const bool use_random=true);

    };
  }
}

#endif
