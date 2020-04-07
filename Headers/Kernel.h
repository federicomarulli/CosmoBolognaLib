/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Headers/Kernel.h
 *
 *  @brief Useful generic functions
 *
 *  This file contains the prototypes of the kernel functions of the
 *  CosmoBolognaLib
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __KERNEL__
#define __KERNEL__

#include <iostream>
#include <sstream>
#include <unordered_map>
#include <fstream>
#include <cmath>
#include <complex>
#include <iomanip>
#include <vector>
#include <limits>
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

/// @cond FFTWinc
#include <Eigen/Dense>
/// @endcond

#include "Constants.h"
#include "EnumCast.h"
#include "Exception.h"


// ============================================================================================


/**
 *  @example vectors.cpp 
 *
 *  This example shows how to print vectors and matrices, and how to
 *  remove elements
 */
/**
 *  @example eigen.cpp 
 *
 *  This example shows how to use eigen vectorization
 */
/**
 *  @example randomNumbers.cpp 
 *
 *  This example shows how to generate random numbers extracted
 *  from a normal distribution
 */
/**
 *  @example randomNumbers_custom.cpp 
 *
 *  This example shows how to generate random numbers extracted from a
 *  generic probability distribution
 */
/**
 *  @example histogram.cpp 
 *
 *  This example shows how use the histogram class
 */
/**
 *  @example correlated_samples.cpp 
 *
 *  This example shows how to generate correlated random samples using
 *  the Cholesky decomposition
 */
/**
 *  @example distances.cpp  
 *
 *  This example shows how to convert redshifts into comoving
 *  distances
 */
/**
 *  @example integration_cuba.cpp 
 *
 *  This example shows how to use the wrapper for CUBA multidimensiona
 *  integration functions
 */
/**
 *  @example integration_gsl.cpp  
 *
 *  This example shows how to integrate a function using the GSL
 *  libraries
 */
/**
 *  @example minimisation_gsl.cpp
 *
 *  This example shows how to minimize a function using the GSL
 *  libraries
 */
/**
 *  @example fits.cpp
 *
 *  This example shows how to read/write a fits file
 */
/**
 *  @example covsample.cpp  
 *
 *  This example shows how to generate correlated samples 
 */
/**
 *  @example cosmology.cpp
 *
 *  This example shows how to set a cosmological model
 */
/**
 *  @example fsigma8.cpp
 *
 *  This example shows how to estimate f*sigma8(z=1)
 */
/**
 *  @example model_cosmology.cpp
 *
 *  This example shows how to estimate f*sigma8(z=1)
 */
/**
 *  @example prior.cpp
 *
 *  This example shows how to menage priors
 */
/**
*  @example data1D.cpp
*
*  this example shows how to construct an object of class Data1D, used
*  to handle 1D datasets of any type
*/
/**
 *  @example fit.cpp
 *
 *  This example shows how to fit a data set with a generic model
 */
/**
 *  @example sampler.cpp
 *
 *  This example shows how to fit a function with a generic model
 */
/**
 *  @example catalogue.cpp
 *
 *  This example shows how to construct a catalogue of extragalactic
 *  objects
 */
/**
 * @example 2pt_monopole.cpp 
 *
 * This example shows how to measure the monopole of the two-point
 * correlation function
 */
/**
 * @example 2pt_monopole_errors.cpp 
 *
 * This example shows how to measure the monopole of the two-point
 * correlation function and estimate the errors with different methods
 */
/**
 * @example 2pt_multipoles.cpp 
 *
 * This example shows how to measure the multipoles of the two-point
 * correlation function using the "direct" and "integrated" method.
 */
/**
 * @example 2pt_2D.cpp 
 *
 * This example shows how to measure the 2D two-point correlation
 * function
 */
/**
 * @example 2pt_projected.cpp 
 *
 * This example shows how to measure the projected two-point
 * correlation function
 */
/**
 * @example 2pt_angular.cpp 
 *
 * This example shows how to measure the angular two-point correlation
 * function
 */
/**
 * @example 3pt.cpp 
 *
 * This example shows how to measure the three-point correlation
 * function
 */
/**
 * @example 3pt_multipoles.cpp 
 *
 * This example shows how to measure the three-point correlation
 * function legendre coefficients
 */
/**
 * @example model_2pt_monopole_BAO.cpp
 *
 * This example shows how to model baryon acoustic oscillations in the
 * monopole of the two-point correlation function
 */
/**
 * @example model_2pt_monopole_RSD.cpp
 *
 * This example shows how to model redshift-space distortions in the
 * monopole of the two-point correlation function
 */
/**
 * @example model_2pt_projected.cpp
 *
 * This example shows how to model the projected two-point correlation
 * function to constrain the linear bias
 */
/**
 * @example model_2pt_2D.cpp
 *
 * This example shows how to model the 2D two-point correlation
 * function in redshift space
 */
/**
 * @example model_2pt_multipoles.cpp
 *
 * This example shows how to model the multipoles of the two-point
 * correlation function in redshift space
 */
/**
 * @example model_3pt.cpp
 *
 * This example shows how to model the reduced 
 * three-point correlation function
 */
/**
 * @example readParameterFile.cpp
 *
 * This example shows how to read parameters from a standard *.ini
 * file
 */
/**
 * @example numberCounts.cpp
 *
 * This example shows how to how to measure the number counts of a catalogue
 */
/**
 * @example numberCounts_errors.cpp
 *
 * This example shows how to how to measure the number counts of a
 * catalogue, computing Poissonian errors
 */
/**
 * @example sizeFunction.cpp
 *
 * This example shows how to compute the theoretical size function of
 * cosmic voids
 */
/**
 * @example cleanVoidCatalogue.cpp
 *
 * This example shows how to clean a cosmic void catalogue, in order
 * to extract cosmological constraints from void counting
 */
/**
 * @example modelling_VoidAbundances.cpp
 *
 * This example shows how to measure and model the void size function,
 * extracting constraints on the cosmological parameters of the model
 */
/** @example table.py
 *
 *  This example shows how read and write files in a Table object 
 */
/** @example distances.py
 *
 *  This example shows how to convert redshifts into comoving
 *  distances 
 */
/**
 *  @example funcgrid_bspline.py
 *
 *  This example shows how to interpolate a function using b-spline
 */
/**
 *  @example catalogue.py
 *
 *  This example shows how to construct a catalogue of extragalactic
 *  objects
 */
/**
 *  @example divide_catalogue.py
 *
 *  This example shows how to divide the a catalogue in sub regions
 */
/**
 *  @example mask_catalogue.py
 *
 *  This example shows how to extract subsamples from a catalogue
 */
/**
 *  @example fft_fftlog.py
 *
 *  This example shows how to computes the discrete Fourier
 * of a logarithmically spaced periodic sequence using the FFTlog libraries
 */
/**
 * @example sizeFunction.py
 *
 * This example shows how to compute the theoretical size function of
 * cosmic voids
 */
/**
 *  @example prior.py
 *
 *  This example shows how to manage priors
 */
/**
 *  @example fit.py
 *
 *  This example shows how to fit a data set with a generic model
 */
/**
 *  @example 2pt_model.py
 *
 *  This example shows how to compute power spectrum and two-point
 *  correlation function models
 */
/**
 *  @example 2pt_model_zErrors.py
 *
 *  This example shows how to compute the two-point correlation
 *  function model taking into account the effect of redshift errors
 */
/**
 *  @example 2pt_monopole.py
 *
 *  This example shows how to compute the monopole of the two-point
 *  correlation function
 */
/**
 * @example 3pt.py
 *
 * This example shows how to measure the three-point correlation
 * function
 */
/**
 * @example cleanVoidCatalogue.py
 *
 *  This example shows how to clean a cosmic void catalogue, making
 *  use of a parameter file, in order to extract cosmological
 *  constraints from void counting
 */
/**
 *  @example covarianceMatrix.ipynb
 *  
 *  This \b notebook shows basic functionalities of the
 *  CovarianceMatrix and TaperedCovarianceMatrix classes
 *
 *  To see the notebook, click here: <a
 *  href="https://github.com/federicomarulli/CosmoBolognaLib/blob/master/Examples/data/covarianceMatrix.ipynb">
 *  notebook</a>
 *
 */
/**
 *  @example analyzeChains.ipynb 
 *  
 *  This \b notebook is aimed to help in post-processing Markov chain
 *  Monte Carlo outputs
 *
 *  To see the notebook, click here: <a
 *  href="https://github.com/federicomarulli/CosmoBolognaLib/blob/master/Examples/statistics/codes/analyzeChains.ipynb">
 *  notebook</a>
 */
/**
 *  @example CLeaning_Algorithm_for_Void_Abundances.ipynb
 *  
 *  This \b notebook explains how to clean cosmic void catalogues and
 *  extract cosmological constraints from void statistics.
 *
 *  To see the notebook, click here: <a
 *  href="https://github.com/federicomarulli/CosmoBolognaLib/blob/master/Examples/cosmicVoids/codes/CLeaning_Algorithm_for_Void_Abundances.ipynb">
 *  notebook</a>
 */
/**
 *  @example 2pt_monopole.ipynb 
 *  
 *  This \b notebook explains how to compute the monopole of the
 *  two-point correlation function
 *
 *  To see the notebook, click here: <a
 *  href="https://colab.research.google.com/drive/1tpQ1Tj06RRwOfZtxYXPC1Sr_Hdb679w8">
 *  notebook</a>
 */
/**
 *  @example no_wiggles_pk.ipynb 
 *  
 *  This \b notebook explains how to compute the de-wiggled power
 *  spectrum model
 */
/**
 *  @example BAO_primer.ipynb 
 *  
 *  This \b notebook explains how to extract cosmological constraints
 *  from Baryon Acoustic Oscillations
 *  
 *  To see the notebook, click here: <a
 *  href="https://colab.research.google.com/drive/1VTn7b0dy8XPDR7at11HecATsXDbu5vzP">
 *  notebook</a>
 */

/**
 *  @brief The global namespace of the <B> \e CosmoBolognaLib </B>
 *  
 *  The \e cbl namespace contains all the main functions and
 *  classes of the CosmoBolognaLib
 */
namespace cbl {

  /**
   *  @enum Dim
   *  @brief the dimension, used e.g. for pair and triplet vectors
   */
  enum class Dim {
    
    /// 1D, used e.g. for 1D pairs, in angular or comoving separations
    _1D_,
    
    /// 2D pair, used e.g. for 2D pairs, in Cartesian or polar coordinates
    _2D_
    
  };

  /**
   * @brief return a vector containing the
   * Dim names
   * @return a vector containing the
   * Dim names
   */
  inline std::vector<std::string> DimNames () { return {"1D", "2D"}; }
  
  /**
   *  @enum BinType
   *  @brief the binning type
   */
  enum class BinType { 

    /// linear binning
    _linear_,
      
    /// logarithmic binning
    _logarithmic_
      
  };

  /**
   * @brief return a vector containing the
   * BinType names
   * @return a vector containing the
   * BinType names
   */
  inline std::vector<std::string> BinTypeNames () { return {"linear", "logarithmic"}; }

  /**
   * @brief cast an enum of type BinType
   * from its index
   * @param binTypeIndex the binType index
   * @return object of class BinType
   */
  inline BinType BinTypeCast (const int binTypeIndex) { return castFromValue<BinType>(binTypeIndex); }

  /**
   * @brief cast an enum of type BinType
   * from its name
   * @param binTypeName the binType name
   * @return object of class BinType
   */
  inline BinType BinTypeCast (const std::string binTypeName) { return castFromName<BinType>(binTypeName, BinTypeNames()); }

  /**
   * @brief cast an enum of type BinType
   * from indeces
   * @param binTypeIndeces the binType indeces
   * @return object of class BinType
   */
  inline std::vector<BinType> BinTypeCast (const std::vector<int> binTypeIndeces) { return castFromValues<BinType>(binTypeIndeces); } 

  /**
   * @brief cast an enum of type BinType
   * from thier names
   * @param binTypeNames the binType names
   * @return objects of class BinType
   */
  inline std::vector<BinType> BinTypeCast (const std::vector<std::string> binTypeNames) { return castFromNames<BinType>(binTypeNames, BinTypeNames()); }

  /**
   *  @enum CoordinateUnits
   *  @brief the coordinate units
   */
  enum class CoordinateUnits {

    /// angle in radians
    _radians_,
    
    /// angle in degrees
    _degrees_,

    /// angle in arcseconds
    _arcseconds_,

    /// angle in arcminutes
    _arcminutes_
    
  };
  

  /**
   * @brief return a std::vector containing the
   * CoordinateUnits names
   * @return a std::vector containing the
   * CoordinateUnits names
   */
  inline std::vector<std::string> CoordinateUnitsNames () { return {"radians", "degrees", "arcseconds", "arcminutes"}; }

  /**
   * @brief cast an enum of type CoordinateUnits
   * from its index
   * @param coordinateUnitsIndex the coordinateUnits index
   * @return object of class CoordinateUnits
   */
  inline CoordinateUnits CoordinateUnitsCast (const int coordinateUnitsIndex) { return castFromValue<CoordinateUnits>(coordinateUnitsIndex); }

  /**
   * @brief cast an enum of type CoordinateUnits
   * from its name
   * @param coordinateUnitsName the coordinateUnits name
   * @return object of class CoordinateUnits
   */
  inline CoordinateUnits CoordinateUnitsCast (const std::string coordinateUnitsName) { return castFromName<CoordinateUnits>(coordinateUnitsName, CoordinateUnitsNames()); }

  /**
   * @brief cast an enum of type CoordinateUnits
   * from indeces
   * @param coordinateUnitsIndeces the coordinateUnits indeces
   * @return object of class CoordinateUnits
   */
  inline std::vector<CoordinateUnits> CoordinateUnitsCast (const std::vector<int> coordinateUnitsIndeces) { return castFromValues<CoordinateUnits>(coordinateUnitsIndeces); } 

  /**
   * @brief cast an enum of type CoordinateUnits
   * from thier names
   * @param coordinateUnitsNames the coordinateUnits names
   * @return objects of class CoordinateUnits
   */
  inline std::vector<CoordinateUnits> CoordinateUnitsCast (const std::vector<std::string> coordinateUnitsNames) { return castFromNames<CoordinateUnits>(coordinateUnitsNames, CoordinateUnitsNames()); }

  
  /**
   *  @enum CoordinateType
   *  @brief the coordinate type
   */
  enum class CoordinateType {

    /// comoving coordinates (x, y, z)
    _comoving_,
    
    /// observed coordinates (R.A., Dec, redshift)
    _observed_
    
  };

  /**
   * @brief return a std::vector containing the
   * CoordinateType names
   * @return a std::vector containing the
   * CoordinateType names
   */
  inline std::vector<std::string> CoordinateTypeNames () { return {"comoving", "observed"}; }

  /**
   * @brief cast an enum of type CoordinateType
   * from its index
   * @param coordinateTypeIndex the coordinateType index
   * @return object of class CoordinateType
   */
  inline CoordinateType CoordinateTypeCast (const int coordinateTypeIndex) { return castFromValue<CoordinateType>(coordinateTypeIndex); }

  /**
   * @brief cast an enum of type CoordinateType
   * from its name
   * @param coordinateTypeName the coordinateType name
   * @return object of class CoordinateType
   */
  inline CoordinateType CoordinateTypeCast (const std::string coordinateTypeName) { return castFromName<CoordinateType>(coordinateTypeName, CoordinateTypeNames()); }

  /**
   * @brief cast an enum of type CoordinateType
   * from indeces
   * @param coordinateTypeIndeces the coordinateType indeces
   * @return object of class CoordinateType
   */
  inline std::vector<CoordinateType> CoordinateTypeCast (const std::vector<int> coordinateTypeIndeces) { return castFromValues<CoordinateType>(coordinateTypeIndeces); } 

  /**
   * @brief cast an enum of type CoordinateType
   * from thier names
   * @param coordinateTypeNames the coordinateType names
   * @return objects of class CoordinateType
   */
  inline std::vector<CoordinateType> CoordinateTypeCast (const std::vector<std::string> coordinateTypeNames) { return castFromNames<CoordinateType>(coordinateTypeNames, CoordinateTypeNames()); }

  struct comovingCoordinates { double xx; double yy; double zz; };
  struct observedCoordinates { double ra; double dec; double redshift; };

  /// Eigen 3D std::vector
  typedef Eigen::Matrix<double, 3, 1> Vector3D;

  /// Eigen 4D matrix
  typedef Eigen::Matrix<double, 4, 1> Vector4D;

  /// Eigen complex std::vector
  typedef Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic> VectorComplex;

  /// typedef of a function returning a double with a double in input
  typedef std::function<double(double)> FunctionDoubleDouble;

  /// typedef of a function returning a double with a vector in input
  typedef std::function<double(std::vector<double>)> FunctionDoubleVector;

  /// typedef of a function returning a double with a vector reference in input
  typedef std::function<double(std::vector<double> &)> FunctionDoubleVectorRef;

  /// typedef of a function returning a double with a double, a pointer and a vector reference in input
  typedef std::function< double(double, std::shared_ptr<void>, std::vector<double> &)> FunctionDoubleDoublePtrVectorRef;

  /// typedef of a function returning a double with two double, a pointer and a vector reference in input
  typedef std::function< double(double, double, std::shared_ptr<void>, std::vector<double> &)> FunctionDoubleDoubleDoublePtrVectorRef;

  /// typedef of a function returning a double with a vector, a pointer and a vector reference in input
  typedef std::function< double(std::vector<double>, std::shared_ptr<void>, std::vector<double> &)> FunctionDoubleVectorPtrVectorRef;

  /// typedef of a function returning a vector with a vector, a pointer and a vector reference in input
  typedef std::function< std::vector<double>(std::vector<double>, std::shared_ptr<void>, std::vector<double> &)> FunctionVectorVectorPtrVectorRef;

  /// typedef of a 3D Tensor of int
  typedef std::vector<std::vector<std::vector<int>>> Tensor3Di;
    
  /// typedef of a 3D Tensor of double
  typedef std::vector<std::vector<std::vector<double>>> Tensor3Dd;
    
  /// typedef of a 4D Tensor of int
  typedef std::vector<std::vector<std::vector<std::vector<int>>>> Tensor4Di;
         

  /**
   *  @name Functions of generic use  
   */
  ///@{

  /**
   *  @brief provide the header for all internal messages
   *  @param stream an std::ostream object
   *  @return the header for internal messages
   */
  inline std::ostream &headerCBL (std::ostream &stream)
  {
    stream << par::col_blue << "CBL > " << par::col_default;
    return stream;
  }
  
#define coutCBL std::cout << headerCBL

  
  /**
   *  @brief set the default directories
   *
   *  @param input_DirCosmo directory where the CosmoBolognaLib are
   *  stored
   *
   *  @param input_DirLoc local directory of the main code
   *
   *  @return none
   */
  inline void SetDirs (const std::string input_DirCosmo, const std::string input_DirLoc)
  { par::DirCosmo = input_DirCosmo; par::DirLoc = input_DirLoc; }
  
  /**
   *  @brief internal CBL warning message
   *
   *  @param msg std::string containing the warning message
   *
   *  @param functionCBL the CBL function in which the warning message
   *  is called
   *
   *  @param fileCBL the CBL file containing the function in which the
   *  warning message is called
   *
   *  @return none
   */
  inline void WarningMsgCBL (const std::string msg, const std::string functionCBL, const std::string fileCBL)
  { std::cerr << std::endl << par::col_bred << "CBL > Warning in the CBL function " << cbl::par::col_yellow << functionCBL << cbl::par::col_bred << " of " << fileCBL << ": " << cbl::par::col_default << msg << std::endl << std::endl; }

  /**
   *  @brief throw an exception
   *  @param msg the message describing the exception
   *  @param exitCode the exit status
   *  @param header header of the error message
   *  @return none
   */
  inline int Error (const std::string msg, const cbl::glob::ExitCode exitCode=cbl::glob::ExitCode::_error_, const std::string header="\n")
  { throw cbl::glob::Exception(msg, exitCode, header); }

  /**
   *  @brief throw an exception: it is used for handling exceptions
   *  inside the CosmoBolognaLib
   *
   *  @param msg the message describing the exception
   *
   *  @param functionCBL the CBL function where the exception is
   *  raised
   *
   *  @param fileCBL the CBL file containing the function where the
   *  exception is raised
   *
   *  @param exitCode the exit status
   *
   *  @return none
   */
  inline int ErrorCBL (const std::string msg, const std::string functionCBL, const std::string fileCBL, const cbl::glob::ExitCode exitCode=cbl::glob::ExitCode::_error_)
  { throw cbl::glob::Exception(msg, exitCode, cbl::par::ErrorMsg, functionCBL, fileCBL); }
  
  /**
   *  @brief produce a beep using the software totem
   *  @return none
   */
  inline void Beep ()
  { if (system("say beep")) {} }

  /**
   *  @brief check if the value of a [string] variable has already
   *  been set
   *
   *  @param var a string variable
   *
   *  @return if var is different than par::defaultString \f$ \rightarrow \f$ true;
   *  else \f$ \rightarrow \f$ false
   */
  inline bool isSet (const std::string var) 
  { return (var!=cbl::par::defaultString) ? true : false; }

  /**
   *  @brief check if the value of a [int] variable has already been
   *  set
   *
   *  @param var a int variable
   *
   *  @return if var>par::defaultInt \f$ \rightarrow \f$ true; else \f$ \rightarrow \f$ false
   */
  inline bool isSet (const int var) 
  { return (var>par::defaultInt) ? true : false; }

  /**
   *  @brief check if the value of a [long] variable has already been
   *  set
   *
   *  @param var a long variable
   *
   *  @return if var>par::defaultLong \f$ \rightarrow \f$ true; else \f$ \rightarrow \f$ false
   */
  inline bool isSet (const long var) 
  { return (var>par::defaultLong) ? true : false; }
  
  /**
   *  @brief check if the value of a [double] variable has already
   *  been set
   *
   *  @param var a double variable
   *
   *  @return if var>par::defaultDouble \f$ \rightarrow \f$ true; else \f$ \rightarrow \f$ false
   */
  inline bool isSet (const double var) 
  { return (var>par::defaultDouble) ? true : false; }
  
  /**
   *  @brief check if the values of a [double] std::vector have already
   *  been set
   *
   *  @param vect a vactor of double values
   *
   *  @return if vect[i]<par::defaultDouble \f$\forall\f$ i \f$ \rightarrow \f$ false;
   *  else \f$ \rightarrow \f$ true
   */
  inline bool isSet (const std::vector<double> vect) 
  {
    bool is = true;
    size_t ind = 0;
    while (is && ind<vect.size()) 
      if (vect[ind++]<par::defaultDouble*1.000001) is = false;
    return is;
  }

  /**
   *  @brief convert a number to a std::string
   *  @param val number of any type
   *  @param fact output format 
   *  @return a std::string containing T
   */
  template <typename T> std::string conv (const T val, const char *fact)
    {
      char VAL[20]; sprintf(VAL, fact, val); 
      return std::string(VAL);
    }
  
  /**
   *  @brief the nearest integer
   *  @param val a number
   *  @return the integer value nearest to val
   */
  template <typename T> 
    int nint (const T val) 
    { return (val<0) ? val-0.5 : val+0.5; }
  
  /**
   *  @brief common logarithm (i.e. logarithm to base 10)
   *  @param val a number
   *  @return if val>0 \f$ \rightarrow \f$ log10(val); else \f$ \rightarrow \f$
   *  par::defaultDouble
   */
  template <typename T> 
    T Log (const T val) 
    { return (val>0) ? log10(val) : par::defaultDouble; }

  /**
   *  @brief given a number x, return the closest of two values a, b
   *  @param x the starting value
   *  @param a the first number to test
   *  @param b the second number to test
   *  @return a if x is closer to a, b if x is closer to b
   */
  template <typename T>
    T closest (T x, T a, T b)
    { 
      if (a>b) ErrorCBL("the input parameter a must be <= than the input parameter b!", "closest", "Kernel.h");
      else if (a==b) return a;
      else return (fabs(x-a) < fabs(x-b)) ? a : b;
      return 1;
    }

  /**
   *  @brief given a number x, return the index of the closest element
   *  to x in vv
   *
   *  @param x the value
   *  @param vv the std::vector
   *  @return the index of the closest element to x in vv
   */
  template <typename T>
    T index_closest (T x, std::vector<T> vv)
    { 
      if (vv.size()==0) ErrorCBL("vv is an empty std::vector!", "index_closest", "Kernel.h");
      std::vector<double>::iterator low, up;
      low = lower_bound(vv.begin(), vv.end(), x);
      up = upper_bound(vv.begin(), vv.end(), x);
      int index = (closest(x, *low, *up)==*low) ? low-vv.begin() : up-vv.begin();
      return index;
    }

  /**
   *  @brief given a number x, return the closest value in a
   *  std::vector
   *
   *  @param x the starting value
   *  @param values std::vector of values
   *  @return the closest value in the std::vector
   */
  template <typename T>
    T closest (T x, std::vector<T> values)
    { return values[index_closest(x, values)]; }
  
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
  
  /**
   *  @brief endian conversion of a short variable
   *  @param s a short variable
   *  @return the converted variable s
   */
  short ShortSwap (const short s);

  /**
   *  @brief endian conversion of an integer variable
   *  @param i an integer variable
   *  @return the converted variable i
   */
  int IntSwap (const int i);

  /**
   *  @brief endian conversion of a long integer variable
   *  @param i a long integer variable
   *  @return the converted variable i
   */
  long long LongSwap (const long long i);

  /**
   *  @brief endian conversion of a float variable
   *  @param f a flot variable
   *  @return the converted variable f
   */
  float FloatSwap (const float f);

  /**
   *  @brief endian conversion of a double variable
   *  @param d a double variable
   *  @return the converted variable d
   */
  double DoubleSwap (const double d);

  /**
   *  @brief reduce the digit figures of an input double
   *
   *  e.g. round(0.2363, 2) = 0.24, round(13.24, 1) = 10
   *
   *  @param num the input double number
   *
   *  @param ndigits the number of digit figures
   *
   *  @return the input number with the required digit figures
   */
  double round_to_digits (const double num, const int ndigits);

  /**
   *  @brief reduce the precision of an input double
   *
   *  e.g. round(0.2363, 2) = 0.23, round(13.24, 1) = 13.2
   *
   *  @param num the input double number
   *
   *  @param ndigits the number of digit figures
   *
   *  @return the input number with the required digit figures
   */
  double round_to_precision (const double num, const int ndigits);

  /**
   *  @brief check if an input file can be opened
   *  @param fin std::ifstream object
   *  @param file the file name
   *  @return none
   */
  void checkIO (const std::ifstream &fin, const std::string file="NULL");

  /**
   *  @brief check if an output file can be opened
   *  @param fout std::ofstream object
   *  @param file the file name
   *  @return none
   */
  void checkIO (const std::ofstream &fout, const std::string file="NULL");

  /**
   *  @brief set evironment variables
   *  @param Var std::vector containing the evironment variables to be set
   *  @return none
   */
  void set_EnvVar (const std::vector<std::string> Var);

  /**
   *  @brief check if an environment variable exists
   *  @param Var the evironment variable to be checked
   *  @return none
   */
  void check_EnvVar (const std::string Var); 

  /**
   *  @brief get the memory used by current process in kB
   *
   *  @warning this function works only on Linux systems
   *
   *  @param type 1 \f$\rightarrow\f$ Physical Memory (RAM); 2 \f$\rightarrow\f$ Virtual
   *  Memory
   *
   *  @return the Physical (RAM) or Virtual Memory used by current
   *  process in kB
   */
  int used_memory (const int type);

  /**
   *  @brief check if the memory used by current process is larger
   *  than a given fraction of the available memory
   *
   *  @warning this function works only on Linux systems
   *
   *  @param frac the fraction of the available memory that is allowed
   *
   *  @param func a std::string that should contain the name of the
   *  function from which check_memory is called; it is used when
   *  printing the error message
   *
   *  @param exit 0 \f$\rightarrow\f$ warning message; 1
   *  \f$\rightarrow\f$ error message; (and exit)
   *
   *  @param type 1 \f$\rightarrow\f$ Physical Memory (RAM); 2
   *  \f$\rightarrow\f$ Virtual Memory
   *
   *  @return 0 \f$\rightarrow\f$ memory problems; 1 \f$\rightarrow\f$
   *  no memory problems
   */
  int check_memory (const double frac, const bool exit=true, const std::string func="", const int type=1);

  ///@}

  
  // ============================================================================================
  

  /**
   *  @name Functions to manipulate std::vectors and matrices
   */
  ///@{

  /**
   *  @brief function to print values with a proper homegenised format
   *  
   *  @param value the value print
   *
   *  @param prec decimal precision
   *
   *  @param ww number of characters to be used as field width
   *
   *  @param header string that is added at the beginning of the line
   *
   *  @param end string that is added at the end of the line
   *
   *  @param use_coutCBL if true, coutCBL is used instead of std::cout
   *
   *  @param stream object of class std::ostream 
   *
   *  @param colour output colour 
   *
   *  @return none
   */
  template <typename T> 
    void Print (const T value, const int prec, const int ww, const std::string header="", const std::string end="\n", const bool use_coutCBL=true, std::ostream &stream=std::cout, const std::string colour=cbl::par::col_default) 
    {
      const int bp = std::cout.precision(); 
      if (fabs(value)<pow(10, -prec) || fabs(value)>pow(10, prec)) {
	if (use_coutCBL)
	  coutCBL << header << colour << std::scientific << std::setprecision(prec) << std::setw(ww) << std::right << value << par::col_default << end;
	else
	  stream << header << std::scientific << std::setprecision(prec) << std::setw(ww) << std::right << value << end;
      }
      else {
	if (use_coutCBL) 
	  coutCBL << header << colour << std::fixed << std::setprecision(prec) << std::setw(ww) << std::right << value << par::col_default << end;
	else 
	  stream << header << std::fixed << std::setprecision(prec) << std::setw(ww) << std::right << value << end;
      }
      std::cout.precision(bp);
    }
  
  /**
   *  @brief print the elements of a std::vector of non string values
   *  on the screen
   *  
   *  @param vect a std::vector
   *  @param prec decimal precision
   *  @param ww number of characters to be used as field width
   *  @return none
   */
  template <typename T> 
    void Print (const std::vector<T> vect, const int prec=4, const int ww=8) 
    {
      for (auto &&vi : vect)
	Print(vi, prec, ww);
    }
  
  /**
   *  @brief print the elements of a std::vector of string values
   *  on the screen
   *  
   *  @param vect a std::vector
   *  @return none
   */
  inline void Print (const std::vector<std::string> vect)
  {
    for (auto &&vi : vect) 
      coutCBL << vi << std::endl;
  }
  

  /**
   *  @brief print the elements of a two std::vectors of non string
   *  values on the screen
   *
   *  @param vect1 a std::vector
   *  @param vect2 a std::vector
   *  @param prec decimal precision
   *  @param ww number of characters to be used as field width
   *  @return none
   */
  template <typename T> 
    void Print (const std::vector<T> vect1, const std::vector<T> vect2, const int prec=4, const int ww=8) 
    {
      if (vect1.size()!=vect2.size())
	ErrorCBL("the two input vectors to be printed must have the same dimension!", "Print", "Kernel.h");
      
      for (size_t i=0; i<vect1.size(); i++) {
	Print(vect1[i], prec, ww, "", "  ");
	Print(vect2[i], prec, ww);
      }
    }


  /**
   *  @brief print the elements of two std::vectors of string values
   *  on the screen
   *
   *  @param vect1 a std::vector
   *
   *  @param vect2 a std::vector
   *
   *  @return none
   */
  inline void Print (const std::vector<std::string> vect1, const std::vector<std::string> vect2) 
    {
      if (vect1.size()!=vect2.size())
	ErrorCBL("the two input vectors to be printed must have the same dimension!", "Print", "Kernel.h");
      
      for (size_t i=0; i<vect1.size(); i++) 
	coutCBL << vect1[i] << "   " << vect2[i] << std::endl;
    }


  /**
   *  @brief print the elements of a three std::vectors of non string
   *  values on the screen
   *
   *  @param vect1 a std::vector
   *
   *  @param vect2 a std::vector
   *
   *  @param vect3 a std::vector
   *
   *  @param prec decimal precision
   *
   *  @param ww number of characters to be used as field width
   *
   *  @return none
   */
  template <typename T> 
    void Print (const std::vector<T> vect1, const std::vector<T> vect2, const std::vector<T> vect3, const int prec=4, const int ww=8) 
    {
      if (vect1.size()!=vect2.size() || vect1.size()!=vect3.size())
	ErrorCBL("the three input vectors to be printed must have the same dimension!", "Print", "Kernel.h");
      
      for (size_t i=0; i<vect1.size(); i++) {
	Print(vect1[i], prec, ww, "", "  ");
	Print(vect2[i], prec, ww, "", "  ");
	Print(vect3[i], prec, ww);	
      }
    }


  /**
   *  @brief print the elements of two std::vectors of string values
   *  on the screen
   *
   *  @param vect1 a std::vector
   *
   *  @param vect2 a std::vector
   *
   *  @param vect3 a std::vector
   *
   *  @return none
   */
  inline void Print (const std::vector<std::string> vect1, const std::vector<std::string> vect2, const std::vector<std::string> vect3) 
    {
      if (vect1.size()!=vect2.size() || vect1.size()!=vect3.size())
	ErrorCBL("the three input vectors to be printed must have the same dimension!", "Print", "Kernel.h");
      
      for (size_t i=0; i<vect1.size(); i++) 
	coutCBL << vect1[i] << "   " << vect2[i] << "   " << vect3[i] << std::endl;
    }
  
  
  /**
   *  @brief print the elements of a matrix of non string values on
   *  the screen
   *
   *  @param mat a matrix (i.e. a std::vector of std::vectors)
   *
   *  @param prec decimal precision
   *
   *  @param ww number of characters to be used as field width
   *
   *  @return none
   */
  template <typename T> 
    void Print (const std::vector<std::vector<T>> mat, const int prec=4, const int ww=8) 
    {
      for (size_t i=0; i<mat.size(); i++) {
	for (size_t j=0; j<mat[i].size(); j++)
	  if (j==0)
	    Print(mat[i][j], prec, ww, "", "  ");
	  else
	    Print(mat[i][j], prec, ww, "", "  "); 
	std::cout << std::endl;
      }
    }

   /**
   *  @brief print the elements of a matrix of string values on the
   *  screen
   *
   *  @param mat a matrix (i.e. a std::vector of std::vectors)
   *
   *  @return none
   */ 
  inline void Print (const std::vector<std::vector<std::string>> mat) 
    {
      for (size_t i=0; i<mat.size(); i++) {
	for (size_t j=0; j<mat[i].size(); j++)
	  if (j==0) 
	    coutCBL << mat[i][j] << "   ";
	  else 
	    std::cout << mat[i][j] << "   ";
	std::cout << std::endl;
      }
    }
  

  /**
   *  @brief minimum element of a std::vector
   *  @param vect a std::vector
   *  @return the minimum element of the std::vector vect
   */
  template <typename T> 
    T Min (const std::vector<T> vect) 
    {
      if (vect.size()==0) ErrorCBL("vect.size=0!", "Min", "Kernel.h");
      return *min_element(vect.begin(), vect.end());
    }

  /**
   *  @brief maximum element of a std::vector
   *  @param vect a std::vector
   *  @return the maximum element of the std::vector vect
   */
  template <typename T> 
    T Max (const std::vector<T> vect) 
    {
      if (vect.size()==0) ErrorCBL("vect.size=0!", "Max", "Kernel.h");
      return *max_element(vect.begin(), vect.end());
    }

  /**
   *  @brief get the unique elements of a std::vector
   *  @param [in] vect_input the input std::vector
   *  @return std::vector containing the unique elements of the input
   *  std::vector
   */
  template <typename T> 
    std::vector<T> different_elements (const std::vector<T> vect_input) 
    {
      std::vector<T> vect = vect_input;
      sort(vect.begin(), vect.end());
      typename std::vector<T>::iterator it = unique(vect.begin(), vect.end()); 
      vect.resize(it-vect.begin());    
      return vect;
    }

  /**
   *  @brief get the number of unique elements of a std::vector
   *  @param vect_input the input std::vector
   *  @return the number of unique elements of the input std::vector
   */
  template <typename T> 
    int N_different_elements (const std::vector<T> vect_input) 
    {
      std::vector<T> vect = different_elements<T>(vect_input);
      return vect.size();
    }

  /**
   *  @brief erase all the equal elements of the input std::vector
   *  @param [in,out] vv a std::vector of integer values
   *  @return none
   */
  void unique_unsorted (std::vector<int> & vv);

  /**
   *  @brief erase all the equal elements of the input std::vector
   *  @param [in,out] vv a std::vector of double values
   *  @return none
   */
  void unique_unsorted (std::vector<double> & vv);

  /**
   *  @brief erase all the equal elements of the input std::vector
   *  @param [in,out] vv a std::vector of integer values
   *  @return none
   */
  void unique_unsorted (std::vector<std::string> & vv);

  /**
   *  @brief erase some elements of a std::vector
   *  @param [in,out] vv a std::vector
   *  @param [in] ind a std::vector containing the elements of the input
   *  std::vector vv to be erased
   *  @return none
   */
  template <typename T> 
    void Erase (std::vector<T> &vv, std::vector<int> ind) 
    {
      for (auto &&i : ind) 
	if (i>=int(vv.size())) ErrorCBL("the input value of ind is too large!", "Erase", "Kernel.h");

      unique_unsorted(ind);
      int tt = 0;
      for (auto &&i : ind) 
	vv.erase(vv.begin()+i-(tt++));
    }

  /**
   *  @brief erase some lines of a matrix
   *  @param [in,out] Mat a matrix (i.e. a std::vector of std::vectors)
   *  @param [in] ll a std::vector containing the lines of the input matrix
   *  Mat to be erased
   *  @return none
   */
  template <typename T> 
    void Erase_lines (std::vector<std::vector<T> > &Mat, std::vector<int> ll) 
    {
      for (auto &&i : ll)
	if (i>=int(Mat.size())) ErrorCBL("the dimension of the input vector ll is too large!", "Erase_lines", "Kernel.h");

      unique_unsorted(ll);
      int tt = 0;
      for (auto &&i : ll)
	Mat.erase(Mat.begin()+i-(tt++));
    }

  /**
   *  @brief erase some columns of a matrix
   *  @param [in,out] Mat a matrix (i.e. a std::vector of std::vectors)
   *  @param [in] col a std::vector containing the columns of the input
   *  matrix Mat to be erased
   *  @return none
   */
  template <typename T> 
    void Erase_columns (std::vector<std::vector<T> > &Mat, std::vector<int> col) 
    {
      for (auto &&i : col)
	for (auto &&j : Mat)
	  if (i>=int(j.size())) ErrorCBL("the dimension of the input vector col is too large!", "Erase_columns", "Kernel.h");

      unique_unsorted(col);
      int tt = 0;
      for (auto &&i : col) {
	for (auto &&j : Mat)
	  j.erase(j.begin()+i-tt);
	tt ++;
      }
    }
 
  /**
   *  @brief select a submatrix containing lines and columns with all
   *  elements major than \e val
   *
   *  @warning this function works only with rectangular matrices and
   *  only for a small set of specific cases: use it with care and
   *  check carefully the results!
   *
   *  @param [in,out] xx the std::vector x 
   *  @param [in,out] yy the std::vector y
   *  @param [in,out] Mat the matrix Mat(x,y) (i.e. a std::vector of
   *  std::vectors)
   *  @param [in] val a number 
   *  @return none
   */
  template <typename T> 
    void SubMatrix (std::vector<T> &xx, std::vector<T> &yy, std::vector<std::vector<T> > &Mat, T val) 
    { 
      std::vector<int> line, column;

      for (unsigned int i=0; i<xx.size(); i++) {
	if (i>=Mat.size()) ErrorCBL("the dimension of the input vector xx is too large!", "SubMatrix", "Kernel.h");
	bool ll = 0;

	for (unsigned int j=0; j<yy.size(); j++) {
	  if (j>=Mat[i].size()) ErrorCBL("the dimension of the input vector yy is too large!", "SubMatrix", "Kernel.h");
	  if (Mat[i][j]<val) {
	    if (j<int(yy.size()*0.5)) {line.push_back(i); ll = 1;}
	    else if (ll==0) {column.push_back(j);}
	  }
	}
      }
      
      Erase(xx, line);
      Erase_lines(Mat, line);
      
      Erase(yy, column);
      Erase_columns(Mat, column);
    }

  /**
   *  @brief check if the dimensions of two std::vectors are equal
   *  @param vect1 a std::vector
   *  @param vect2 a std::vector
   *  @return 0 \f$\rightarrow\f$ the dimensions are different; 1 \f$\rightarrow\f$ the
   *  dimensions are equal
   */
  template <typename T> 
    bool isDimEqual (const std::vector<T> vect1, const std::vector<T> vect2) 
    {
      return (vect1.size()==vect2.size()) ? 1 : 0;
    }

  /**
   *  @brief check if the dimensions of two matrices are equal
   *  @param mat1 a matrix
   *  @param mat2 a matrix
   *  @return 0 \f$\rightarrow\f$ the dimensions are different; 1 \f$\rightarrow\f$ the
   *  dimensions are equal
   */
  template <typename T> 
    bool isDimEqual (const std::vector<std::vector<T> > mat1, const std::vector<std::vector<T> > mat2) 
    {
      bool is = (mat1.size()==mat2.size()) ? 1 : 0;
      if (is) 
	for (unsigned int i=0; i<mat1.size(); i++) {
	  if (mat1[i].size()!=mat2[i].size()) is = 0;
	}
      return is;
    }

  /**
   *  @brief check if the dimension of a std::vector is equal/lower
   *  than an input value
   *
   *  @param vect a std::vector
   *   
   *  @param val the input value
   *
   *  @param vector the name of the std::vector (used only to write
   *  the error message)
   *
   *  @param equal true \f$\rightarrow\f$ check if the dimension is
   *  equal to val; false \f$\rightarrow\f$ check if the dimension is
   *  lower than val
   *
   *  @return none
   */
  template <typename T> 
    void checkDim (const std::vector<T> vect, const int val, const std::string vector, bool equal=true) 
    {
      if (equal) {
	if ((int)vect.size()!=val) 
	  ErrorCBL("the dimension of " + vector + " is: " + conv(vect.size(), par::fINT) + " ( != " + conv(val, par::fINT) + " )", "checkDim", "Kernel.h");
      }
      else { 
	if ((int)vect.size()<val)
	  ErrorCBL("the dimension of " + vector + " is: " + conv(vect.size(), par::fINT) + " ( < " + conv(val, par::fINT) + " )", "checkDim", "Kernel.h");
      }
    }
  
  /**
   *  @brief check if the dimensions of a matrix are higher than two
   *  input values
   *  @param mat a matrix
   *  @param val_i an input value
   *  @param val_j an input value
   *  @param matrix the name of the matrix (using only to write the
   *  error message)
   *  @param equal true \f$\rightarrow\f$ check if the dimension is equal to val;
   *  false \f$\rightarrow\f$ check if the dimension is lower than val
   *  @return none
   */
  template <typename T> 
    void checkDim (const std::vector<T> mat, const int val_i, const int val_j, const std::string matrix, const bool equal=true) 
    {
      if (equal) {
	if (int(mat.size())!=val_i) 
	  ErrorCBL("the dimension of: " + matrix + " is:" + conv(mat.size(), par::fINT) + " != " + conv(val_i, par::fINT) + "!", "checkDim", "Kernel.h");
	else 
	  for (size_t k=0; k<mat.size(); k++)
	    if (int(mat[k].size())!=val_j) 
	      ErrorCBL("the dimension of: " + matrix + " is:" + conv(mat[k].size(), par::fINT) + " != " + conv(val_j, par::fINT) + "!", "checkDim", "Kernel.h");
      }
      else {
	if (int(mat.size())<val_i) 
	  ErrorCBL("the dimension of: " + matrix + " is:" + conv(mat.size(), par::fINT) + " < " + conv(val_i, par::fINT) + "!", "checkDim", "Kernel.h");
	else 
	  for (size_t k=0; k<mat.size(); k++)
	    if (int(mat[k].size())<val_j) 
	      ErrorCBL("the dimension of: " + matrix + " is:" + conv(mat[k].size(), par::fINT) + " < " + conv(val_j, par::fINT) + "!", "checkDim", "Kernel.h");
      }
    }

  /**
   *  @brief check if two std::vectors are equal
   *
   *  @param vect1 a std::vector
   *
   *  @param vect2 a std::vector
   *
   *  @return error if the given std::vectors are different; nothing
   *  if they are equal
   */
  template <typename T> 
    void checkEqual (const std::vector<T> vect1, const std::vector<T> vect2) 
    {
      checkDim(vect2, vect1.size(), "vect2");
      for (size_t i=0; i<vect1.size(); i++)
	if (vect1[i]!=vect2[i])
	  ErrorCBL("vect1 and vect2 are different!", "checkEqual", "Kernel.h");
    }
  
  /**
   *  @brief fill a std::vector with linearly spaced values
   *  @param [in] nn the number of steps, i.e. the final dimension of
   *  vv
   *  @param [in] min the minimum value of the range of values
   *  @param [in] max the maximum value of the range of values
   *  @return none
   */
  template <typename T> 
    std::vector<T> linear_bin_vector (const size_t nn, const T min, const T max)
    {
      std::vector<T> vv(nn);
      for (size_t i = 0; i<nn; i++)
	vv[i] = min+(max-min)*T(i)/T(nn-1);
      return vv;
    }

  /**
   *  @brief fill a std::vector with logarithmically spaced values
   *  @param [in] nn the number of steps, i.e. the final dimension of
   *  vv
   *  @param [in] min the minimum value of the range of values
   *  @param [in] max the maximum value of the range of values
   *  @return none
   */
  template <typename T> 
    std::vector<T> logarithmic_bin_vector (const size_t nn, const T min, const T max)
    {
      std::vector<T> vv(nn);
      for (size_t i=0; i<nn; i++)
	vv[i] = exp(log(min)+(log(max)-log(min))*T(i)/T(nn-1));
      return vv;
    }

  /**
   *  @brief locate a value in a given std::vector
   *  @author Carlo Giocoli
   *  @author cgiocoli@gmail.com
   *  @param vv a std::vector of generic values
   *  @param xx a generic number
   *  @return the std::vector index i such that vv[i]~xx
   */
  template <typename T> 
    int locate (const std::vector<T> &vv, const T xx) 
    {
      size_t nn = vv.size ();
      int jl = -1;
      int ju = nn;
      bool as = (vv[nn-1] >= vv[0]);
      while (ju-jl > 1)
	{
	  int jm = (ju+jl)*0.5;
	  if ((xx >= vv[jm]) == as)
	    jl = jm;
	  else
	    ju = jm;
	}
      if (xx == vv[0])
	return 0;
      else if (xx == vv[nn-1])
	return nn-2;
      else
	return jl;
    }

  /**
   *  @brief extract elements from a given vector
   *
   *  @param vec the input vector
   *  @param index vector containing the index of values
   *  to extract
   *  @return the vector with requested elements
   */
  template <typename T>
    std::vector<T> extract_elements (std::vector<T> vec, std::vector<unsigned int> index)
    {
      std::vector<T> vv;
      for (unsigned int i=0; i< index.size(); i++)
	vv.push_back(vec[index[i]]);
      return vv;
    }
  
  /**
   *  @brief flatten a \f$ left( n\times m \right)  \f$ matrix
   *  in a vector of size \f$x\times m\f$.
   *
   *  @param matrix the input matrix
   *  
   *  @return the vectorized matrix
   */
  template <typename T>
    std::vector<T> flatten(std::vector<std::vector<T>> matrix)
    {
      std::vector<T> flatten;
      for (size_t i=0; i<matrix.size(); i++)
	for (size_t j=0; j<matrix[i].size(); j++)
	  flatten.push_back(matrix[i][j]);

      return flatten;
    }

  /**
   *  @brief reshape a vector into a matrix of
   *  given number of rows and columns
   *
   *  @param vec the input vector
   *  @param size1 the number of rows
   *  @param size2 the number of columns
   *  
   *  @return the reshaped matrix
   */
  template <typename T>
    std::vector<std::vector<T>> reshape (std::vector<T> vec, const int size1, const int size2)
    {
      if (size1*size2!=int(vec.size()))
	ErrorCBL("sizes does not match! "+conv(size1*size2, par::fINT)+" should be equal to "+conv(int(vec.size()), par::fINT), "reshape", "Kernel.h");

      std::vector<std::vector<T>> matrix(size1, std::vector<T> (size2, 0));

      for (int i=0; i<size1; i++)
	for (int j=0; j<size2; j++)
	  matrix[i][j] = vec[j+i*size2];

      return matrix;
    }
  
  /**
   *  @brief transpose a matrix
   *
   *  @param matrix the input matrix
   *  
   *  @return the transposed matrix
   */
  template <typename T>
    std::vector<std::vector<T>> transpose (std::vector<std::vector<T>> matrix)
    {
      const int size1 = matrix.size();
      const int size2 = matrix[0].size();
      
      std::vector<std::vector<T>> TRmatrix(size2, std::vector<T> (size1, 0));

      for (int i=0; i<size1; i++)
	for (int j=0; j<size2; j++)
	  TRmatrix[j][i] = matrix[i][j];

      return TRmatrix;
    }

  /**
   *  @brief return false \f$ \rightarrow \f$ value
   *  outside the range; true \f$ \rightarrow \f$ value
   *  inside the range
   *
   *  @param value the value
   *  @param min the lower limit
   *  @param max the upper limit
   *  @param include_limits true \f$ \rightarrow \f$ include limits;
   *  false \f$ \rightarrow \f$ exclude limits.
   *  
      @return false \f$ \rightarrow \f$ 
   *  values outside the range; 
   *  true \f$ \rightarrow \f$ value inside the range
   */
  template <typename T>
    bool inRange (T value, T min, T max, bool include_limits = true) 
    {
      if (include_limits){
	if (value>=min && max>=value)
	  return true;
      }
      else{
	if (value>min && max>value)
	  return true;
      }

      return false;
    }

  /**
   *  @brief return false \f$ \rightarrow \f$ if values are
   *  outside the range; true \f$ \rightarrow \f$ if values
   *  are inside the range
   *
   *  @param value vector containing the values
   *  @param min vector containing the lower limits
   *  @param max vector containing the upper limits
   *  @param include_limits true \f$ \rightarrow \f$ include limits;
   *  false \f$ \rightarrow \f$ exclude limits.
   *  
   *  @return vector of booleans: false \f$ \rightarrow \f$ 
   *  values outside the range; 
   *  true \f$ \rightarrow \f$ values inside the range
   */
  template <typename T>
    bool inRange (std::vector<T> value, std::vector<T> min, std::vector<T> max, bool include_limits = true) 
    {
      bool in_range = true;

      for (size_t i=0; i<value.size(); i++)
	in_range *= inRange(value[i], min[i], max[i], include_limits);

      return in_range;
    }

  /**
   *  @brief return false \f$ \rightarrow \f$ values
   *  outside the range; true \f$ \rightarrow \f$ values
   *  inside the range
   *
   *  @param value vector containing the values
   *  @param ranges vector containing the ranges
   *  @param include_limits true \f$ \rightarrow \f$ include limits;
   *  false \f$ \rightarrow \f$ exclude limits.
   *  
   *  @return vector of booleans: false \f$ \rightarrow \f$ 
   *  if values outside the range; 
   *  true \f$ \rightarrow \f$ values inside the range
   */
  template <typename T>
    bool inRange (std::vector<T> value, std::vector<std::vector<T>> ranges, bool include_limits = true) 
    {
      bool in_range = true;

      for (size_t i=0; i<value.size(); i++)
	in_range *= inRange(value[i], ranges[i][0], ranges[i][1], include_limits);

      return in_range;
    }

  /**
   *  @brief return the value of 
   *  \f[ \vec{x} M \vec{x}^T \f]
   *
   *  @param [in] vv the std::vector v
   *  @param [in] MM the matrix M
   *  @return result
   */
  template <typename T> 
    T v_M_vt (const std::vector<T> vv, const std::vector<std::vector<T>> MM)
    {
      const int size = vv.size();

      std::vector<double> ivv(size, 0);
      for(int i=0; i<size; i++)
	for(int j=0; j<size; j++)
	  ivv[i] += vv[j]*MM[i][j];

      double res=0;
      for(int i=0; i<size; i++)
	res += vv[i]*ivv[i];

      return res;
    }

  // sort two or more std::vectors at the same time
  namespace glob {
    class CL {
    public:
      std::vector<double> VV;
      CL (std::vector<double> vv) {VV = vv;};
    };
    /// @cond glob
    bool operator<(const CL &, const CL &);
    /// @endcond
  }
  
  /**
   *  @brief sort the elements of a std::vectors, and the elements of
   *  a second std::vector according to the first sorting
   *
   *  @param p1 iterator to the first std::vector
   *  @param p2 iterator to the second std::vector
   *  @param dim dimension of the two std::vectors 
   *  @return none
   */
  void sort_2vectors (std::vector<double>::iterator p1, std::vector<double>::iterator p2, const int dim);

  /**
   *  @brief sort the elements of a std::vectors, and the elements of
   *  two other std::vectors according to the first sorting
   *
   *  @param p1 iterator to the first std::vector
   *  @param p2 iterator to the second std::vector
   *  @param p3 iterator to the third std::vector
   *  @param dim dimension of the three std::vectors 
   *  @return none
   */
  void sort_3vectors (std::vector<double>::iterator p1, std::vector<double>::iterator p2, std::vector<double>::iterator p3, const int dim);

  /**
   *  @brief sort the elements of a std::vectors, and the elements of three
   *  other std::vectors according to the first sorting
   *
   *  @param p1 iterator to the first std::vector
   *  @param p2 iterator to the second std::vector
   *  @param p3 iterator to the third std::vector
   *  @param p4 iterator to the four std::vector
   *  @param dim dimension of the four std::vectors 
   *  @return none
   */
  void sort_4vectors (std::vector<double>::iterator p1, std::vector<double>::iterator p2, std::vector<double>::iterator p3, std::vector<double>::iterator p4, const int dim);

  /**
   *  @brief function to create multiple directories 
   *
   *  http://mylinuxtechcorner.blogspot.com/2012/09/c-version-for-mkdir-p.html
   *
   *  @param path the name of the directory (or directories) to be
   *  created recursively (with parents)
   *
   *  @param rootPath the path to the root directory
   *
   *  @param mode the permissions for the directory
   *
   *  @param verbose if true it shows a warning message when the
   *  directory to be created already exists
   *
   *  @return 0 if no error occurs
   */
  int makeDir (std::string path, const std::string rootPath=".", const mode_t mode=0777, const bool verbose=false);
  
  /**
   *  @brief matrix multiplication
   *
   *  overloading of the * operator used to multiplicate two matrices
   *
   *  @param Mat1 STL std::vector of std::vectors, i.e. a matrix
   *  @param Mat2 STL std::vector of std::vectors, i.e. a matrix
   *  @return Mat1*Mat2
   */
  inline std::vector<std::vector<double> > operator * (const std::vector<std::vector<double> > &Mat1, const std::vector<std::vector<double> > &Mat2)
  {   
    std::vector<std::vector<double> > MatP(Mat1.size(), std::vector<double>(Mat2[0].size(),0.));
  
    for (unsigned int i=0; i<Mat1.size(); i++) 
      for (unsigned int j=0; j<Mat2[0].size(); j++) {
	double temp = 0.;
	for (unsigned int k=0; k<Mat1[0].size(); k++) 
	  temp += Mat1[i][k]*Mat2[k][j];
	MatP[i][j] = temp;
      }

    return MatP;
  }

  /**
   *  @brief slice a std::vector from start to stop
   *
   *  @param v original std::vector
   *  @param start starting position
   *  @param end ending position
   *
   *  @return the sliced std::vector
   */
  template <typename T>
    std::vector<T> slice (const std::vector<T> v, const int start=0, const  int end=-1)
    {
      int oldlen = v.size();
      int newlen;

      if (end==-1 || end>=oldlen)
	newlen = oldlen-start;
      else 
	newlen = end-start;
      
      std::vector<T> nv(newlen);

      for (int i=0; i<newlen; i++) 
	nv[i] = v[start+i];
      
      return nv;
    }
  
  ///@}
  
}


#endif
