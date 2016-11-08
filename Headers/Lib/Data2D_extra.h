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
 *  @file Headers/Lib/Data2D_extra.h
 *
 *  @brief The class Data2D_extra
 *
 *  This file defines the interface of the class Data2D_extra
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unbo.it
 */

#ifndef __DATA2D_EXTRA__
#define __DATA2D_EXTRA__

#include "Data2D.h"

namespace cosmobl {

  namespace data {
    
    /**
     *  @class Data2D_extra Data2D_extra.h
     *  "Headers/Lib/Data2D_extra.h"
     *
     *  @brief The class Data2D_extra
     *
     *  This is the base class used to manage 2D data with extra
     *  information
     */
    class Data2D_extra : public Data2D
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
       *  @return an object of class Data2D_extra
       */
      Data2D_extra () { set_dataType(DataType::_2D_data_extra_); }

      /**
       *  @brief constructor
       *  @param x vector containing x points 
       *  @param y vector containing y points 
       *  @param fxy vector containing f(x,y) values
       *  @param extra_info vector containing vectors of extra generic
       *  information
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @param ymin maximun value of y to be used 
       *  @param ymax maximun value of y to be used 
       *  @param dataType the data type
       *  @return an object of class Data2D_extra
       */
      Data2D_extra (const vector<double> x, const vector<double> y, const vector<vector<double>> fxy, const vector<vector<double>> extra_info, const double xmin=par::defaultDouble, const double xmax=-par::defaultDouble, const double ymin=par::defaultDouble, const double ymax=-par::defaultDouble, const DataType dataType=DataType::_2D_data_extra_)
	: Data2D(x, y, fxy, xmin, xmax, ymin, ymax, dataType), m_extra_info(extra_info) {}

      /**
       *  @brief constructor
       *  @param x vector containing x points 
       *  @param y vector containing y points 
       *  @param fxy vector containing f(x,y) values
       *  @param error_fxy vector containing error on f(x,y)
       *  @param extra_info vector containing vectors of extra generic
       *  information
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @param ymin maximun value of y to be used 
       *  @param ymax maximun value of y to be used 
       *  @param dataType the data type
       *  @return shared pointer to an object of class Data2D_extra
       */
      Data2D_extra (const vector<double> x, const vector<double> y, const vector<vector<double>> fxy, const vector<vector<double>> error_fxy, const vector<vector<double>> extra_info, const double xmin=par::defaultDouble, const double xmax=-par::defaultDouble, const double ymin=par::defaultDouble, const double ymax=-par::defaultDouble, const DataType dataType=DataType::_2D_data_extra_)
	: Data2D(x, y, fxy, error_fxy, xmin, xmax, ymin, ymax, dataType), m_extra_info(extra_info) {}

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Data2D_extra () = default;

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
      double extra_info (const int i, const int j) const { return m_extra_info[i][j]; }

      /**
       *  @brief return the m_exta_info vector
       *  @return vector containing the extra information
       */
      vector<vector<double>> extra_info () const { return m_extra_info; }
      
      ///@}

      
      /**
       *  @name Member functions used to set the private members
       */
      ///@{
      
      /**
       *  @brief set interval variable m_error_fx
       *  @param extra_info vector containing vectors with extra info
       *  @return none
       */
      void set_extra_info (const vector<vector<double>> extra_info) { m_extra_info = extra_info; }

      ///@}


      /**
       *  @name Input/Output member functions
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
       *  @param full false &rarr; simply store the data; true &rarr;
       *  duplicate the data in the other three quadrands (usefull
       *  e.g. when storing the 2D correlation function)
       *  @param rank cpu index (for MPI usage)
       *  @return none
       */
      virtual void write (const string dir, const string file, const string header, const bool full, const int rank=0) const override;
      
      ///@}
      
    };
    
  }
}

#endif
