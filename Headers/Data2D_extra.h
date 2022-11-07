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
 *  @authors federico.marulli3@unibo.it
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
       */
      Data2D_extra () { set_dataType(DataType::_2D_extra_); }
      
      /**
       *  @brief constructor which reads the data from file
       *
       *  @param n_extra_info number of extra info
       *
       *  @param input_file the input data file
       *
       *  @param skip_nlines the header lines to be skipped
       *
       *  @param column vector containing the column(s) of x (and y)
       *  values in the input file; if it is not provided, the first
       *  column will be used by default
       *
       *  @param column_data the column of data values in the input
       *  file; the size of column_data is the number of data to be
       *  read (e.g. the size should be 3 in the case of the 3
       *  multipole moments of the two-point correlation function); if
       *  the size of column_data is larger than 1, more than 1 data
       *  vectors are read and then added one after the other in a
       *  single data object; if column_data is not provided, the
       *  first column after column will be used by default, assuming
       *  that only 1 data vector has to be read
       *
       *  @param column_errors the column of error values in the input
       *  file; the size of column_error must be equal to the size of
       *  column_data; if the size of column_error is larger than 1,
       *  more than 1 error vectors are read and then added one after
       *  the other in a single data object; if column_random is not
       *  provided, the second column after column will be used by
       *  default, assuming that only 1 random vector has to be read;
       *  if the input file has less than 4 columns, the errors will
       *  be set to 1
       *
       *  @param column_edges vector containing the columns of x and y
       *  bin edge values in the input file; if it is not provided,
       *  the third and four columns after the column of x values will
       *  be used; if these columns do no exist the edges are not read
       */
      Data2D_extra (const int n_extra_info, const std::string input_file, const int skip_nlines=0, const std::vector<int> column={1, 2}, const std::vector<int> column_data={}, const std::vector<int> column_errors={}, const std::vector<int> column_edges={})
	: Data2D()
      { m_extra_info.resize(n_extra_info); read(input_file, skip_nlines, column, column_data, column_errors, column_edges); set_dataType(DataType::_2D_extra_); }
      
      /**
       *  @brief constructor which gets the data from an input matrix
       *
       *  @param x vector containing the x values
       * 
       *  @param y vector containing the y values
       * 
       *  @param data matrix containing the f(x,y) values
       *
       *  @param extra_info vector containing vectors of extra generic
       *  information
       *
       *  @param bin_edges_x the x variable bin edges
       *
       *  @param bin_edges_y the y variable bin edges
       */
      Data2D_extra (const std::vector<double> x, const std::vector<double> y, const std::vector<std::vector<double>> data, const std::vector<std::vector<double>> extra_info, const std::vector<double> bin_edges_x={}, const std::vector<double> bin_edges_y={})
	: Data2D(x, y, data, bin_edges_x, bin_edges_y), m_extra_info(extra_info) { set_dataType(DataType::_2D_extra_); }

      /**
       *  @brief constructor which gets both the data and the errors
       *  from input matrices
       *
       *  @param x vector containing the x values
       * 
       *  @param y vector containing the y values 
       *
       *  @param data vector containing the data
       *
       *  @param error vector containing the errors
       *
       *  @param extra_info vector containing vectors of extra generic
       *  information
       *
       *  @param bin_edges_x the x variable bin edges
       *
       *  @param bin_edges_y the y variable bin edges
       */
      Data2D_extra (const std::vector<double> x, const std::vector<double> y, const std::vector<double> data, const std::vector<double> error, const std::vector<std::vector<double>> extra_info, const std::vector<double> bin_edges_x={}, const std::vector<double> bin_edges_y={}) : Data2D(x, y, data, error, bin_edges_x, bin_edges_y), m_extra_info(extra_info) { set_dataType(DataType::_2D_extra_); }
      
      /**
       *  @brief constructor which gets both the data and the
       *  covariance matrix from input matrices
       *
       *  @param x vector containing the x values 
       *
       *  @param y vector containing the y values
       *
       *  @param data matrix containing the data
       *
       *  @param covariance matrix containing the covariance
       *
       *  @param extra_info vector containing vectors of extra generic
       *  information
       *
       *  @param bin_edges_x the x variable bin edges
       *
       *  @param bin_edges_y the y variable bin edges
       */
      Data2D_extra (const std::vector<double> x, const std::vector<double> y, const std::vector<std::vector<double>> data, const std::vector<std::vector<double>> covariance, const std::vector<std::vector<double>> extra_info, const std::vector<double> bin_edges_x={}, const std::vector<double> bin_edges_y={}) : Data2D(x, y, data, covariance, bin_edges_x, bin_edges_y), m_extra_info(extra_info) { set_dataType(DataType::_2D_extra_); }
      
      /**
       *  @brief constructor which gets both the data and the
       *  covariance matrix from input matrices
       *
       *  @param x vector containing the x values 
       *
       *  @param y vector containing the y values 
       *
       *  @param data matrix containing the data
       *
       *  @param covariance matrix containing the covariance
       *
       *  @param extra_info vector containing vectors of extra generic
       *  information
       *
       *  @param bin_edges_x the x variable bin edges
       *
       *  @param bin_edges_y the y variable bin edges
       */
      Data2D_extra (const std::vector<double> x, const std::vector<double> y, const std::vector<double> data, const std::vector<std::vector<double>> covariance, const std::vector<std::vector<double>> extra_info, const std::vector<double> bin_edges_x={}, const std::vector<double> bin_edges_y={}) : Data2D(x, y, data, covariance, bin_edges_x, bin_edges_y), m_extra_info(extra_info) { set_dataType(DataType::_2D_extra_); }

      /**
       *  @brief default destructor
       */
      virtual ~Data2D_extra () = default;

      /**
       *  @brief static factory used to construct objects of class
       *  Data2D_extra
       *
       *  @return a shared pointer to an object of class Data2D_extra
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
       *  @brief get the independet variable, to be used in model
       *  computation
       *
       *  @param i index of the extra_info containing the first
       *  independent variable
       *
       *  @param j index of the extra_info containing the second
       *  independent variable
       *
       *  @return the independent variable
       */
      std::vector<std::vector<double>> IndipendentVariable (const int i=-1, const int j=-1) const 
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
       */
      void set_extra_info (const std::vector<std::vector<double>> extra_info) { m_extra_info = extra_info; }

      ///@}
      
      
      /**
       *  @name Input/Output member functions
       */
      ///@{
      
      /**
       *  @brief read the data
       *
       *  @param input_file input data file
       *
       *  @param skip_nlines the header lines to be skipped
       * 
       *  @param column vector containing the column(s) of x (and y)
       *  values in the input file; if it is not provided, the first
       *  column will be used by default
       *
       *  @param column_data the column of data values in the input
       *  file; the size of column_data is the number of data to be
       *  read (e.g. the size should be 3 in the case of the 3
       *  multipole moments of the two-point correlation function); if
       *  the size of column_data is larger than 1, more than 1 data
       *  vectors are read and then added one after the other in a
       *  single data object; if column_data is not provided, the
       *  first column after column will be used by default, assuming
       *  that only 1 data vector has to be read
       *
       *  @param column_errors the column of error values in the input
       *  file; the size of column_error must be equal to the size of
       *  column_data; if the size of column_error is larger than 1,
       *  more than 1 error vectors are read and then added one after
       *  the other in a single data object; if column_random is not
       *  provided, the second column after column will be used by
       *  default, assuming that only 1 random vector has to be read;
       *  if the input file has less than 4 columns, the errors will
       *  be set to 1
       *
       *  @param column_edges vector containing the columns of x and y
       *  bin edge values in the input file; if it is not provided,
       *  the third and four columns after the column of x values will
       *  be used; if these columns do no exist the edges are not read
       */
      virtual void read (const std::string input_file, const int skip_nlines=0, const std::vector<int> column={1, 2}, const std::vector<int> column_data={}, const std::vector<int> column_errors={}, const std::vector<int> column_edges={}) override;

      /**
       *  @brief print the data on screen
       *
       *  @param precision the float precision
       */
      virtual void Print (const int precision=4) const override;
      
      /**
       *  @brief write the data
       *  @param dir output directory
       *  @param file output file
       *  @param header text with the variable names to be written at
       *  the first line of the output file
       *  @param full false &rarr; simply store the data; true &rarr;
       *  duplicate the data in the other three quadrands (usefull
       *  e.g. when storing the 2D correlation function)
       *  @param prec the floating point precision for the output
       *  file
       *  @param ww number of characters to be used as field width
       *  @param rank cpu index (for MPI usage)
       */
      virtual void write (const std::string dir, const std::string file, const std::string header, const bool full, const int prec=4, const int ww=8, const int rank=0) const override;
      
      ///@}

      
      /**
       *  @name Member functions for data cut
       */

      ///@{

      /**
       * @brief cut the data
       * @param xmin minumum value for the independent variable x
       * @param xmax maximum value for the independent variable x
       * @param ymin minumum value for the independent variable y
       * @param ymax maximum value for the independent variable y
       * @return pointer to an object of type Data2D
       */
      std::shared_ptr<Data> cut (const double xmin, const double xmax, const double ymin, const double ymax) const;

      ///@}
      
    };
    
  }
}

#endif
