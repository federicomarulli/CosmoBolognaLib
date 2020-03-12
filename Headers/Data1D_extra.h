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
 *  @file Headers/Data1D_extra.h
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

namespace cbl {

  namespace data {
    
    /**
     *  @class Data1D_extra Data1D_extra.h
     *  "Headers/Data1D_extra.h"
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
      std::vector<std::vector<double>> m_extra_info;

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
      Data1D_extra ()
	: Data1D() { set_dataType(DataType::_1D_extra_); }

      /**
       *  @brief constructor which reads the data from file
       *
       *  @param n_extra_info number of extra info
       *
       *  @param input_file the input data file
       *
       *  @param skip_nlines the header lines to be skipped
       *
       *  @param column_x the column of x values in the input file; if
       *  it is not provided, the first column will be used by default
       *
       *  @param column_data the column of data values in the input
       *  file; the size of column_data is the number of data to be
       *  read (e.g. the size should be 3 in the case of the 3
       *  multipole moments of the two-point correlation function); if
       *  the size of column_data is larger than 1, more than 1 data
       *  vectors are read and then added one after the other in a
       *  single data object; if column_data is not provided, the
       *  first column after column_x will be used by default,
       *  assuming that only 1 data vector has to be read
       *
       *  @param column_errors the column of error values in the input
       *  file; the size of column_error must be equal to the size of
       *  column_data; if the size of column_error is larger than 1,
       *  more than 1 error vectors are read and then added one after
       *  the other in a single data object; if column_random is not
       *  provided, the second column after column_x will be used by
       *  default, assuming that only 1 random vector has to be read;
       *  if the input file has only 2 columns, the errors will be set
       *  to 1
       *
       *  @return object of class Data1D
       */
      Data1D_extra (const int n_extra_info, const std::string input_file, const int skip_nlines=0, const int column_x=0, const std::vector<int> column_data={}, const std::vector<int> column_errors={})
	: Data1D()
	{ m_extra_info.resize(n_extra_info); read(input_file, skip_nlines, column_x, column_data, column_errors); set_dataType(DataType::_1D_extra_); }
      
      /**
       *  @brief constructor which gets the data from input vectors
       *
       *  @param x vector containing the x values 
       *
       *  @param data vector containing the data values 
       *
       *  @param extra_info vector containing vectors of extra generic
       *  information
       *
       *  @return an object of class Data1D_extra
       */
      Data1D_extra (const std::vector<double> x, const std::vector<double> data, const std::vector<std::vector<double>> extra_info)
	: Data1D(x, data), m_extra_info(extra_info)
      { set_dataType(DataType::_1D_extra_); }

      /**
       *  @brief constructor which gets both the data and the errors
       *  from input vectors
       *   
       *  @param x vector containing the x values 
       *
       *  @param data vector containing the data
       *
       *  @param error vector containing the errors
       *
       *  @param extra_info vector containing vectors of extra generic
       *  information
       *
       *  @return an object of class Data1D_extra
       */
      Data1D_extra (const std::vector<double> x, const std::vector<double> data, const std::vector<double> error, const std::vector<std::vector<double>> extra_info={})
	: Data1D(x, data, error), m_extra_info(extra_info)
      { set_dataType(DataType::_1D_extra_); }
      
      /**
       *  @brief constructor which gets both the data and the
       *  covariance matrix from input vectors
       *
       *  @param x vector containing the x values 
       *
       *  @param data vector containing the data
       *
       *  @param covariance matrix containing the covariance 
       * 
       *  @param extra_info vector containing vectors of extra generic
       *  information
       *
       *  @return an object of class Data1D_extra
       */
      Data1D_extra (const std::vector<double> x, const std::vector<double> data, const std::vector<std::vector<double>> covariance, const std::vector<std::vector<double>> extra_info)
	: Data1D(x, data, covariance), m_extra_info(extra_info)
      { set_dataType(DataType::_1D_extra_); }

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Data1D_extra () = default;

      /**
       *  @brief static factory used to construct objects of class
       *  Data1D_extra
       *  @return a shared pointer to an object of class Data
       */
      std::shared_ptr<Data> as_factory () { return move(std::unique_ptr<Data1D_extra>(this)); }

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
      std::vector<std::vector<double>> extra_info () const override { return m_extra_info; }

      /**
       *  @brief get the independet variable, to be used 
       *  in model computation
       *  @param i index of the extra_info containing
       *  the independent variable
       *  @param j index of the second extra_info, not used
       *  @return the independent variable
       */
      std::vector<std::vector<double>> IndipendentVariable (const int i=-1, const int j=-1) const 
      {
	(void)j;
	std::vector<std::vector<double>> iv;
	iv.push_back(((i>0) ? m_extra_info[i] : m_x));
	return iv;
      }

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
      void set_extra_info (const std::vector<std::vector<double>> extra_info) override { m_extra_info = extra_info; }

      ///@}

      
      /**
       *  @name Member functions for Input/Output 
       */
      ///@{

      /**
       *  @brief read the data
       *
       *  @param input_file the input data file
       *
       *  @param skip_nlines the header lines to be skipped
       *
       *  @param column_x the column of x values in the input file; if
       *  it is not provided, the first column will be used by default
       *
       *  @param column_data the column of data values in the input
       *  file; the size of column_data is the number of data to be
       *  read (e.g. the size should be 3 in the case of the 3
       *  multipole moments of the two-point correlation function); if
       *  the size of column_data is larger than 1, more than 1 data
       *  vectors are read and then added one after the other in a
       *  single data object; if column_data is not provided, the
       *  first column after column_x will be used by default,
       *  assuming that only 1 data vector has to be read
       *
       *  @param column_errors the column of error values in the input
       *  file; the size of column_error must be equal to the size of
       *  column_data; if the size of column_error is larger than 1,
       *  more than 1 error vectors are read and then added one after
       *  the other in a single data object; if column_random is not
       *  provided, the second column after column_x will be used by
       *  default, assuming that only 1 random vector has to be read;
       *  if the input file has only 2 columns, the errors will be set
       *  to 1
       *
       *  @return none
       */
      virtual void read (const std::string input_file, const int skip_nlines=0, const int column_x=0, const std::vector<int> column_data={}, const std::vector<int> column_errors={}) override;

      /**
       *  @brief print the data on screen
       *
       *  @param precision the float precision
       *
       *  @return none
       */
      virtual void Print (const int precision=4) const override;
      
      /**
       *  @brief write the data
       *  @param dir output directory
       *  @param file output file
       *  @param header text with the variable names to be written at
       *  the first line of the output file
       *  @param prec the floating point precision
       *  @param ww number of characters to be used as field width
       *  @param rank cpu index (for MPI usage)
       *  @return none
       */
      virtual void write (const std::string dir, const std::string file, const std::string header, const int prec=4, const int ww=8, const int rank=0) const override;

      ///@}


      /**
       *  @name Member functions for data cut
       */
      
      ///@{
      
      /**
       * @brief cut the data
       * @param xmin minumum value for the independent variable x
       * @param xmax maximum value for the independent variable x
       * @return pointer to an object of type Data1D_extra
       */
      std::shared_ptr<Data> cut (const double xmin, const double xmax) const override;

      /**
       * @brief cut the data, for Data1D_extra
       * @param [in] mask vector containing values to be masked
       * @return pointer to an object of type Data1D_extra
       */
      std::shared_ptr<Data> cut (const std::vector<bool> mask) const;

      ///@}
      
    };
    
  }
}

#endif
