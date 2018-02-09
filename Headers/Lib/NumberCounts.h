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
 *  @file Headers/Lib/NumberCounts.h
 *
 *  @brief The class NumberCounts
 *
 *  This file defines the interface of the class NumberCounts,
 *  used to measure the number counts
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __NCOUNTS__
#define __NCOUNTS__


#include "Catalogue.h"
#include "Measure.h"


// ===================================================================================================


namespace cosmobl {

  namespace measure {

    /**
     *  @brief The namespace of the <B> number counts
     *  </B>
     *  
     *  The \e measure::numbercounts namespace contains all the functions and
     *  classes to measure the number counts
     */
    namespace numbercounts {

      /**
       * @enum NCType
       * @brief the number counts type
       */
      enum NCType { 

	_dn_dV_,

	_dn_dlogV_,

	_N_V_,

	_n_V_
      };


      /**
       *  @class NumberCounts NumberCounts.h
       *  "Headers/Lib/NumberCounts.h"
       *
       *  @brief The class NumberCounts
       *
       *  This is the base class used to measure the 
       *  number counts
       *
       */
      class NumberCounts : public Measure {

	protected :

	  /**
	   *  @name Input catalogue
	   */
	  ///@{

	  /// input data catalogue
	  shared_ptr<catalogue::Catalogue> m_data;

	  ///@}

	  /**
	   *  @name number counts variable
	   */
	  ///@{

	  /// the catalogue variable
	  catalogue::Var m_variable;

	  /// the minimum of the variable
	  double m_MinVar;
	  
	  /// the maximum of the variable
	  double m_MaxVar;

	  ///@}

	  /**
	   *  @name Number Counts function data
	   */
	  ///@{

	  /// number counts type
	  NCType m_ncType;

	  ///@}
	  ///@}
	  /**
	   *  @name Protected member functions to measure the number counts
	   */
	  ///@{

	  /**
	   *  @brief measure the number counts with Poisson 
	   *  errors
	   *
	   *  @param nbin the number of bins
	   *
	   *  @param bintype the binning type
	   *  
	   *  @param Volume the volume of the catalogue
	   *
	   *  @return none
	   */
	  virtual shared_ptr<data::Data> m_measurePoisson (const int nbin=10, const binType bintype=binType::_linear_, const double Volume=1.) const;

	  /**
	   *  @brief measure the number counts with Jackknife 
	   *  covariance matrix
	   *
	   *  @param nbin the number of bins
	   *
	   *  @param bintype the binning type
	   *  
	   *  @param Volume the volume of the catalogue
	   *
	   *  @param dir_output_resample output directory of the
	   *  resampled correlation function
	   *
	   *  @return none
	   */
	  virtual shared_ptr<data::Data> m_measureJackknife (const int nbin=10, const binType bintype=binType::_linear_, const double Volume=1., const string dir_output_resample=par::defaultString) const;

	  /**
	   *  @brief measure the number counts with
	   *  Bootstrap covariance matrix
	   *
	   *  @param nbin the number of bins
	   *
	   *  @param bintype the binning type
	   *  
	   *  @param Volume the volume of the catalogue
	   *
	   *  @param dir_output_resample output directory of the
	   *  resampled correlation function
	   *
	   *  @param nMocks number of resampling used for bootstrap
	   *
	   *  @param seed the seed for random number generation
	   *
	   *  @return none
	   */
	  virtual shared_ptr<data::Data> m_measureBootstrap (const int nbin=10, const binType bintype=binType::_linear_, const double Volume=1., const string dir_output_resample=par::defaultString, const int nMocks=0, const int seed=3213) const;

	  /**
	   *  @brief measure the number counts with Jackknife 
	   *  covariance matrix
	   *
	   *  @param nbin the number of bins
	   *
	   *  @param bintype the binning type
	   *  
	   *  @param Volume the volume of the catalogue
	   *
	   *  @param dir_output_resample output directory of the
	   *  resampled correlation function
	   *
	   *  @return none
	   */
	  virtual shared_ptr<data::Data> m_measureJackknifeObjects (const int nbin=10, const binType bintype=binType::_linear_, const double Volume=1., const string dir_output_resample=par::defaultString) const;

          /**
           *  @brief measure the number counts with
           *  Bootstrap covariance matrix
           *
           *  @param nbin the number of bins
           *
           *  @param bintype the binning type
           *  
           *  @param Volume the volume of the catalogue
           *
           *  @param dir_output_resample output directory of the
           *  resampled correlation function
           *
           *  @param nMocks number of resampling used for bootstrap
           *
           *  @param seed the seed for random number generation
           *
           *  @return none
           */
          virtual shared_ptr<data::Data> m_measureBootstrapObjects (const int nbin=10, const binType bintype=binType::_linear_, const double Volume=1., const string dir_output_resample=par::defaultString, const int nMocks=0, const int seed=3213) const;

	  ///@}

	public:

	  /**
	   *  @name Constructors/destructors
	   */
	  ///@{

	  /**
	   *  @brief default constructor
	   *
	   *  @return object of class NumberCounts
	   */
	  NumberCounts () {}

	  /**
	   *  @brief constructor
	   *
	   *  @param data object of class Catalogue containing the input
	   *  catalogue
	   *
	   *  @param variable the number counts variable
	   *
	   *  @param MinVar the minimum of the variable
	   *  
	   *  @param MaxVar the maximum of the variable
	   *
	   *  @param ncType the type of number counts
	   *
	   *  @return object of class NumberCounts
	   */
	  NumberCounts (const catalogue::Catalogue data, const catalogue::Var variable, const double MinVar=par::defaultDouble, const double MaxVar=par::defaultDouble, const NCType ncType=NCType::_dn_dV_) : 
	    m_data(make_shared<catalogue::Catalogue>(catalogue::Catalogue(move(data)))), m_variable(variable), m_MinVar(MinVar), m_MaxVar(MaxVar), m_ncType(ncType) {}

	  /**
	   *  @brief default destructor
	   *  @return none
	   */
	  virtual ~NumberCounts () = default;

	  ///@}

	  /**
	   *  @name Member functions to set protected members
	   */
	  ///@{

	  /**
	   *  @brief add a data catalogue
	   *  @param data object of class Catalogue 
	   *  @return none
	   */
	  void set_data (const catalogue::Catalogue data)
	  { m_data = make_shared<catalogue::Catalogue>(catalogue::Catalogue(move(data))); }

	  /**
	   *  @brief set the number counts variable
	   *  @param variable  the number counts variable 
	   *  @param MinVar minimum value of the variable
	   *  @param MaxVar var maximum value of the variable
	   *  @return none
	   */
	  void set_variable (const catalogue::Var variable, const double MinVar=par::defaultDouble, const double MaxVar=par::defaultDouble)
	  { m_variable = variable; m_MinVar = MinVar; m_MaxVar=MaxVar; }

	  /**
	   *  @brief set the number counts type
	   *  
	   *  @param ncType the type of number counts
	   *
	   *  @return none
	   */
	  void set_NumberCounts_Type (const NCType ncType) { m_ncType = ncType; }

	  ///@}

	  /**
	   *  @name Member functions to get protected members
	   */
	  ///@{

	  /**
	   *  @brief function to get the protected member m_data
	   *  @return return the protected member m_data
	   */
	  catalogue::Catalogue catalogue () { return m_data; }

	  /**
	   *  @brief get the number counts variable
	   *  @return the protected member m_variable
	   */
	  catalogue::Var variable () { return m_variable;}

	  /**
	   *  @brief get the number counts type
	   *  
	   *  @return the protected member m_ncType;
	   */
	  NCType set_NumberCounts_Type () { return m_ncType; }

	  ///@}
	  /**
	   *  @name Member functions to measure the number counts
	   */
	  ///@{

	  /**
	   *  @brief measure the number counts
	   *
	   *  @param errorType type of error
	   *
	   *  @param nbin the number of bins
	   *
	   *  @param bintype the binning type
	   *  
	   *  @param Volume the volume of the catalogue
	   *
	   *  @param dir_output_resample output directory of the
	   *  resampled correlation function
	   *
	   *  @param nMocks number of resampling used for bootstrap
	   *
	   *  @param seed the seed for random number generation
	   *
	   *  @return none
	   */
	  virtual void measure (const ErrorType errorType=ErrorType::_Poisson_, const int nbin=10, const binType bintype=binType::_linear_, const double Volume=1., const string dir_output_resample=par::defaultString, const int nMocks=0, const int seed=3213);

	  ///@}

	  /**
	   *  @name Input/Output member functions (customized in all the derived classes)
	   */
	  ///@{

	  /**
	   *  @brief read the measured two-point correlation
	   *  @param dir input directory
	   *  @param file input file
	   *  @return none
	   */
	  virtual void read (const string dir, const string file);

	  /**
	   *  @brief write the measured two-point correlation
	   *  @param dir output directory
	   *  @param file output file
	   *  @param rank cpu index (for MPI usage)
	   *  @return none
	   */
	  virtual void write (const string dir=par::defaultString, const string file=par::defaultString, const int rank=0);

	  ///@}

	  /**
	   *  @name Member functions to estimate the errors and covariance matrices
	   */
	  ///@{ 

	  /**
	   *  @brief read the measured covariance matrix
	   *  @param dir input directory
	   *  @param file input file
	   *  @return none
	   */
	  virtual void read_covariance (const string dir, const string file);

	  /**
	   *  @brief write the measured two-point correlation
	   *  @param dir output directory
	   *  @param file output file
	   *  @return none
	   */
	  virtual void write_covariance (const string dir, const string file) const;

	  /**
	   *  @brief compute the covariance matrix
	   *  @param xi vector containing the measure correlation
	   *  functions used to compute the covariance matrix
	   *  @param JK true &rarr; compute the jackknife covariance
	   *  matrix; false compute the standard covariance matrix
	   *  @return none
	   */
	  virtual void compute_covariance (const vector<shared_ptr<data::Data>> xi, const bool JK);

	  /**
	   *  @brief compute the covariance matrix
	   *  @param file vector containing the input files with the
	   *  measured correlation functions used to compute the
	   *  covariance matrix
	   *  @param JK true &rarr; compute the jackknife covariance
	   *  matrix; false compute the standard covariance matrix
	   *  @return none
	   */
	  virtual void compute_covariance (const vector<string> file, const bool JK);

	  ///@}

      };
    }
  }
}

#endif
