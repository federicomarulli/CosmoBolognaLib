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
 *  @file Headers/Lib/LogNormal.h
 *
 *  @brief Implementation of the lognormal data structure
 *
 *  This file defines the interface of the class LogNormal
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __LOGNORMAL__
#define __LOGNORMAL__

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
     *  @class LogNormal LogNormal.h "Headers/Lib/LogNormal.h"
     *
     *  @brief The class LogNormal
     *
     *  This class is used to handle objects of type <EM> lognormal
     *  </EM>
     */
    class LogNormal_full {
    
    protected:
    
      /// pointer to the fiducial cosmology
      shared_ptr<cosmology::Cosmology> m_cosmology;
    
      /// string containing the author for the model power spectrum
      string m_author;
    
      /// 0 &rarr; compute the linear power spectrum; 1 &rarr; compute the non-linear power spectrum
      bool m_NL;

      default_random_engine m_generator;

      /// the clustering signal 
      shared_ptr<data::Field3D> m_clustering_signal;

      /// approximate cell size of the density field
      double m_rmin;

      int m_nx;

      int m_ny;

      int m_nz;

      int m_nzF;

      double m_xMin;

      double m_xMax;

      double m_yMin;

      double m_yMax;

      double m_zMin;

      double m_zMax;

      int m_nCells_eff;
      
      /// the number of lognormal realization to produce/read
      int m_nLN;
    
      /// smoothing radius
      double m_rsmooth;

      bool m_use_visibility;

      /// pointer to the input datasets
      vector<shared_ptr<catalogue::Catalogue>> m_data;
    
      /// pointer to the random catalogues
      vector<shared_ptr<catalogue::Catalogue>> m_random;

      /// vector containing pointers to the LogNormal realizations
      shared_ptr<data::Field3D> m_density;  

      /// vector containing pointers to the LogNormal realizations
      shared_ptr<data::Field3D> m_densityG;  

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

      shared_ptr<glob::FuncGrid> m_func_redshift;

      shared_ptr<glob::FuncGrid> m_func_DC;

      shared_ptr<glob::FuncGrid> m_func_HH;

      vector<shared_ptr<glob::FuncGrid> > m_func_bias;

      shared_ptr<glob::FuncGrid> m_func_growth_factor;

      shared_ptr<glob::FuncGrid> m_func_growth_rate;

      shared_ptr<glob::FuncGrid> m_func_pk;

      void set_grid_parameters();

      void set_fields(const bool use_random, const bool doRSD);

      void reset_fields(const bool use_random, const bool doRSD);

      void set_clustering_signal();
      
      void set_density_field(const double smoothing_radius);
      
   //   void set_displacement_field();
      
   //   void set_rsd_displacement_field();

      void set_potential();

      void set_radial_velocity();

      void set_visibility();

      void set_visibility_from_random();
      
      void set_visibility_from_random_RSD();

      void extract_points_lognormal_field(const double nObjects, const bool doRSD, const vector<double> redshift, const vector<double> bias, const shared_ptr<data::Field3D> visibility, const string file_out);

    public:

      /**
       *  @brief default constructor
       *  @return object of class LogNormal
       */
      LogNormal_full() {};

      LogNormal_full(const cosmology::Cosmology cosmology, const double redshift_min=0., const double redshift_max=10., const int nredshift=500, const string author="CAMB");

      LogNormal_full(const double rmin, const double xMin, const double xMax, const double yMin, const double yMax, const double zMin, const double zMax, const cosmology::Cosmology cosmology, const double redshift_min=0., const double redshift_max=10., const int nredshift=500, const string author="CAMB");

      LogNormal_full(const double rmin, const vector<shared_ptr<catalogue::Catalogue>> random, const double pad, const cosmology::Cosmology cosmology, const double redshift_min=0., const double redshift_max=10., const int nredshift=500, const string author="CAMB");

      /**
       *  @brief default destructor
       *  @return none
       */
      ~LogNormal_full () = default;

      void set_cosmo_function(const cosmology::Cosmology cosmology, const double redshift_min=0., const double redshift_max=10., const int nredshift=500, const string author="CAMB");

      void set_grid_parameters(const double rmin, const double xMin, const double xMax, const double yMin, const double yMax, const double zMin, const double zMax);

      void set_grid_parameters(const double rmin, const vector<shared_ptr<catalogue::Catalogue>> random, const double pad);

      void generate_lognormal(const int start, const int stop, const bool doRSD, const double smoothing_radius, const vector<double> nObjects, const vector<vector<double>> redshift, const vector<vector<double> > bias, const string dir, const string filename, const int seed, const bool setfields=1, const bool use_random = 1);

    };

  }

}

#endif
