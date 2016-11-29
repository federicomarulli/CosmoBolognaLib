/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli                          *
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
 *  @file Headers/Lib/Data1D_extra.h
 *
 *  @brief The class Data1D_extra
 *
 *  This file defines the interface of the class Data1D_extra
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unbo.it
 */

#ifndef __DATA1D_EXTRA__
#define __DATA1D_EXTRA__

#include "Data1D.h"

namespace cosmobl {

  namespace data {
    
    /**
     *  @class Data1D_extra Data1D_extra.h
     *  "Headers/Lib/Data1D_extra.h"
     *
     *  @brief The class Data1D_extra
     *
     *  This is the base class used to manage 1D data with extra
     *  information
     */
    class Data1D_extra : public Data1D
    {
      
    protected:

      /**
       *  @name Data input
       */
      ///@{
      
      /// extra information
      vector<vector<double>> m_extra_info;

      ///@}

      
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Data1D_extra
       */
      Data1D_extra () { set_dataType(DataType::_1D_data_extra_); }
      
      /**
       *  @brief constructor
       *  @param x vector containing x points 
       *  @param fx vector containing f(x) 
       *  @param extra_info vector containing vectors of extra generic
       *  information
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @param dataType the data type
       *  @return an object of class Data1D_extra
       */
      Data1D_extra (const vector<double> x, const vector<double> fx, const vector<vector<double>> extra_info, const double xmin=par::defaultDouble, const double xmax=-par::defaultDouble, const DataType dataType=DataType::_1D_data_extra_)
	: Data1D(x, fx, xmin, xmax, dataType), m_extra_info(extra_info) {}

      /**
       *  @brief constructor
       *  @param x vector containing x points 
       *  @param fx vector containing f(x) 
       *  @param error_fx vector containing error on f(x) 
       *  @param extra_info vector containing vectors of extra generic
       *  information
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @param dataType the data type
       *  @return an object of class Data1D_extra
       */
      Data1D_extra (const vector<double> x, const vector<double> fx, const vector<double> error_fx, const vector<vector<double>> extra_info={}, const double xmin=par::defaultDouble, const double xmax=-par::defaultDouble, const DataType dataType=DataType::_1D_data_extra_)
	: Data1D(x, fx, error_fx, xmin, xmax, dataType), m_extra_info(extra_info) {}
      
      /**
       *  @brief constructor
       *  @param x vector containing x points 
       *  @param fx vector containing f(x) 
       *  @param covariance vector containing f(x) covariance matrix 
       *  @param extra_info vector containing vectors of extra generic
       *  information
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used
       *  @param dataType the data type
       *  @return an object of class Data1D_extra
       */
      Data1D_extra (const vector<double> x, const vector<double> fx, const vector<vector<double>> covariance, const vector<vector<double>> extra_info, const double xmin=par::defaultDouble, const double xmax=-par::defaultDouble, const DataType dataType=DataType::_1D_data_extra_)
	: Data1D(x, fx, covariance, xmin, xmax, dataType), m_extra_info(extra_info) {}

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Data1D_extra () = default;

      ///@}


      /**
       *  @name Member functions used to get the private members
       */
      ///@{ 

      /**
       *  @brief return the value of the extra information at index i,j
       *  @param i index
       *  @param j index
       *  @return the value of the extra_info at position i,j
       */
      double extra_info (const int i, const int j) const override { return m_extra_info[i][j]; }

      /**
       *  @brief return the m_exta_info vector
       *  @return vector containing vectors with extra information
       */
      vector<vector<double>> extra_info () const override { return m_extra_info; }
      
      ///@}

      
      /**
       *  @name Member functions used to set the private members
       */
      ///@{
      
      /**
       *  @brief set interval variable m_error_fx
       *  @param extra_info vector containing vectors with extra
       *  information
       *  @return none
       */
      void set_extra_info (const vector<vector<double>> extra_info) override { m_extra_info = extra_info; }

      ///@}

      
      /**
       *  @name Member functions for Input/Output 
       */
      ///@{

      /**
       *  @brief read the data
       *  @param input_file input data file
       *  @param skip_nlines the header lines to be skipped
       *  @return none
       */
      virtual void read (const string input_file, const int skip_nlines=0) override;

      /**
       *  @brief write the data
       *  @param dir output directory
       *  @param file output file
       *  @param header text with the variable names to be written at
       *  the first line of the output file
       *  @param rank cpu index (for MPI usage)
       *  @return none
       */
      virtual void write (const string dir, const string file, const string header, const int rank=0) const override;

      ///@}
      
    };
    
  }
}

#endif
