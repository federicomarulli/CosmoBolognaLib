/********************************************************************
 *  Copyright (C) 2022 by Sofia Contarini                           *
 *  sofia.contarini3@unibo.it                                       *
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
 *  @file Headers/Path.h
 *
 *  @brief The class Path used to handle the Cosmobolognalib paths
 *
 *  This file defines the interface of the class Path, used to handle
 *  the Cosmobolognalib paths
 *
 *  @author Sofia Contarini
 *
 *  @author sofia.contarini3@unibo.it
 */

#ifndef __PATH__
#define __PATH__

#include <iostream>
#include <sstream>
#include <unordered_map>
#include <fstream>
#include <cmath>
#include <complex>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <memory>
#include <numeric>
#include <functional>
#include <stdlib.h>
#include <unistd.h>
#include <random>
#include <map>
#include <omp.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <ctype.h>
#include <sys/stat.h>
#include <errno.h>

#ifdef LINUX
#include "sys/types.h"
#include "sys/sysinfo.h"
#endif

// Save compiler switches
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wparentheses"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wunused-variable"

/// @cond GSLinc
#include <gsl/gsl_errno.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
/// @endcond

/// @cond FFTWinc
#include <fftw3.h>
/// @endcond

/// @cond BOOSTinc
#include <boost/numeric/odeint.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/beta.hpp>
/// @endcond

/// @cond FFTWinc
#include <Eigen/Dense>
/// @endcond

// Restore compiler switches
#pragma GCC diagnostic pop



/**
 *  @brief The global namespace of the <B> \e CosmoBolognaLib </B>
 *  
 *  The \e cbl namespace contains all the main functions and
 *  classes of the CosmoBolognaLib
 */
namespace cbl {

  /**
   *  @class Path
   *
   *  @brief The class Path
   *
   *  This class defines the interface of the class Path, used to handle
   *  the Cosmobolognalib main paths
   *
   */
   class Path
  {
      
    /// directory where the CosmoBolognaLib are stored
    std::string m_DirCosmo;
  
    /// local directory of the main code
    std::string m_DirLoc;

      
  public:
    
    /**
     *  @name Constructors/destructors
     */
    ///@{

    /**
     *  @brief default constructor
     *  
     */
    Path ();
    
    /**
     *  @brief default destructor
     */
    ~Path () = default;
 
    /**
     *  @brief get the directory where the CosmoBolognaLbi are stored
     *
     *  @return m_DirCosmo, that is the directory where the
     *  CosmoBolognaLbi are stored
     */
    inline std::string DirCosmo () { return m_DirCosmo; }

    /**
     *  @brief get the local directory
     *
     *  @return m_DirLoc, that is the local directory
     */
    inline std::string DirLoc () { return m_DirLoc; }

    /**
     *  @brief set the default directories
     *
     *  @param input_DirCosmo directory where the CosmoBolognaLib are
     *  stored
     *
     *  @param input_DirLoc local directory of the main code
     */
    inline void SetDirs (const std::string input_DirCosmo, const std::string input_DirLoc="./")
    { m_DirCosmo = input_DirCosmo; m_DirLoc = input_DirLoc; }
       
    /**
     *  @brief substitute ~ with the full path
     *
     *  @param path the relative path
     *
     *  @param isDir true \f$\rightarrow\f$ directory path, false
     *  \f$\rightarrow\f$ otherwise
     *
     *  @return std::string containing the full path
     */
    std::string fullpath (std::string path, const bool isDir=true);
      
  };
  
}


#endif
