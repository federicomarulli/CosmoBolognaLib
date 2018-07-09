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
 *  @file Headers/Data2D_extra.h
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

namespace cbl {

  namespace data {
    
    /**
     *  @class Data2D_extra Data2D_extra.h
     *  "Headers/Data2D_extra.h"
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
      std::vector<std::vector<double>> m_extra_info;

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
      Data2D_extra () { set_dataType(DataType::_2D_extra_); }

      /**
       *  @brief constructor
       *  @param x vector containing x points 
       *  @param y vector containing y points 
       *  @param data vector containing f(x,y) values
       *  @param extra_info vector containing vectors of extra generic
       *  information
       *  @return an object of class Data2D_extra
       */
      Data2D_extra (const std::vector<double> x, const std::vector<double> y, const std::vector<std::vector<double>> data, const std::vector<std::vector<double>> extra_info)
	: Data2D(x, y, data), m_extra_info(extra_info) { set_dataType(DataType::_2D_extra_); }

      /**
       *  @brief constructor
       *  @param x vector containing x points 
       *  @param y vector containing y points 
       *  @param data vector containing data values
       *  @param covariance vector containing covariance matrix
       *  @param extra_info vector containing vectors of extra generic
       *  information
       *  @return an object of class Data2D_extra
       */
      Data2D_extra (const std::vector<double> x, const std::vector<double> y, const std::vector<std::vector<double>> data, const std::vector<std::vector<double>> covariance, const std::vector<std::vector<double>> extra_info) : Data2D(x, y, data, covariance), m_extra_info(extra_info) { set_dataType(DataType::_2D_extra_); }

      /**
       *  @brief constructor
       *  @param x vector containing x points 
       *  @param y vector containing y points 
       *  @param data vector containing data values
       *  @param error vector containing error values
       *  @param extra_info vector containing vectors of extra generic
       *  information
       *  @return shared pointer to an object of class Data2D_extra
       */
      Data2D_extra (const std::vector<double> x, const std::vector<double> y, const std::vector<double> data, const std::vector<double> error, const std::vector<std::vector<double>> extra_info) : Data2D(x, y, data, error), m_extra_info(extra_info) { set_dataType(DataType::_2D_extra_); }

      /**
       *  @brief constructor
       *  @param x vector containing x points 
       *  @param y vector containing y points 
       *  @param data vector containing data values
       *  @param covariance vector containing covariance matrix
       *  @param extra_info vector containing vectors of extra generic
       *  information
       *  @return shared pointer to an object of class Data2D_extra
       */
      Data2D_extra (const std::vector<double> x, const std::vector<double> y, const std::vector<double> data, const std::vector<std::vector<double>> covariance, const std::vector<std::vector<double>> extra_info) : Data2D(x, y, data, covariance), m_extra_info(extra_info) { set_dataType(DataType::_2D_extra_); }

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Data2D_extra () = default;

      /**
       *  @brief static factory used to construct objects of class
       *  Data1D
       *  @return a shared pointer to an object of class Data
       */
      std::shared_ptr<Data> as_factory () { return move(std::unique_ptr<Data2D_extra>(this)); }

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
      std::vector<std::vector<double>> extra_info () const { return m_extra_info; }

      /**
       *  @brief get the independet variable, to be used 
       *  in model computation
       *  @param i index of the extra_info containing
       *  the first independent variable
       *  @param j index of the extra_info containing
       *  the second independent variable
       *  @return the independent variable
       */
      std::vector<std::vector<double>> IndipendentVariable(const int i=-1, const int j=-1) const 
      {
	std::vector<std::vector<double>> iv;
	iv.push_back(((i>0) ? m_extra_info[i] : m_x));
	iv.push_back(((j>0) ? m_extra_info[i] : m_y));
	return iv;
      }

      ///@}

      
      /**
       *  @name Member functions used to set the private members
       */
      ///@{
      
      /**
       *  @brief set the extra info
       *  @param extra_info vector containing vectors with extra info
       *  @return none
       */
      void set_extra_info (const std::vector<std::vector<double>> extra_info) { m_extra_info = extra_info; }

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
      virtual void read (const std::string input_file, const int skip_nlines=0) override;

      /**
       *  @brief write the data
       *  @param dir output directory
       *  @param file output file
       *  @param header text with the variable names to be written at
       *  the first line of the output file
       *  @param full false &rarr; simply store the data; true &rarr;
       *  duplicate the data in the other three quadrands (usefull
       *  e.g. when storing the 2D correlation function)
       *  @param precision the floating point precision for the output
       *  file
       *  @param rank cpu index (for MPI usage)
       *  @return none
       */
      virtual void write (const std::string dir, const std::string file, const std::string header, const bool full, const int precision=4, const int rank=0) const override;
      
      ///@}
      
    };
    
  }
}

#endif
