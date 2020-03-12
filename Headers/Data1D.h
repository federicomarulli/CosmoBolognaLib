/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/Data1D.h
 *
 *  @brief The class Data1D
 *
 *  This file defines the interface of the class Data1D
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __DATA1D__
#define __DATA1D__

#include "Data.h"

namespace cbl {

  namespace data {
    
    /**
     *  @class Data1D Data1D.h
     *  "Headers/Data1D.h"
     *
     *  @brief The class Data1D
     *
     *  This is the base class used to manage 1D data
     */
    class Data1D : public Data
    {
      
    protected:

      /**
       *  @name Data input
       */
      ///@{

      /// ordered x axis points
      std::vector<double> m_x;

      ///@}
     
      /// number of points along x
      int m_xsize;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Data1D
       */
      Data1D ()
	: Data(DataType::_1D_) {}

      /**
       *  @brief constructor which reads the data from file
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
      Data1D (const std::string input_file, const int skip_nlines=0, const int column_x=1, const std::vector<int> column_data={}, const std::vector<int> column_errors={})
	: Data(cbl::data::DataType::_1D_)
	{ read(input_file, skip_nlines, column_x, column_data, column_errors); }

      /**
       *  @brief constructor which gets the data from input vectors
       *
       *  @param x vector containing the x values 
       *
       *  @param data vector containing data
       *
       *  @return object of class Data1D
       */
      Data1D (const std::vector<double> x, const std::vector<double> data)
	: Data(cbl::data::DataType::_1D_, data)
	{ set_xx(x); }
      
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
       *  @return object of class Data1D
       */
      Data1D (const std::vector<double> x, const std::vector<double> data, const std::vector<double> error) 
	: Data(cbl::data::DataType::_1D_, data, error)
	{ set_xx(x); }
      
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
       *  @return object of class Data1D
       */
      Data1D (const std::vector<double> x, const std::vector<double> data, const std::vector<std::vector<double>> covariance)
	: Data(cbl::data::DataType::_1D_, data, covariance)
	{ set_xx(x); }
      
      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Data1D () = default;

      /**
       *  @brief static factory used to construct objects of class
       *  Data1D
       *  @return a shared pointer to an object of class Data
       */
      std::shared_ptr<Data> as_factory ();

      ///@}

      /**
       *  @name Member functions used to get the private members
       */
      ///@{

      /**
       *  @brief get the number of points along x
       *  @return the number of points along x
       */
      int xsize () const
      { return m_xsize; }

      /**
       *  @brief get the value of x at the i-th bin
       *  @param i index
       *  @return value of x at the i-th bin
       */
      double xx (const int i) const override
      { return m_x[i]; }  

      /**
       *  @brief get the x std::vector
       *  @return the x std::vector
       */
      std::vector<double> xx () const override
	{ return m_x; }  

      /**
       *  @brief get the independet variable, to be used 
       *  in model computation
       *  @param i index of the extra_info containing
       *  the independent variable
       *  @param j index of the extra_info containing the
       *  independent variable, not used
       *  @return the independent variable
       */
      std::vector<std::vector<double>> IndipendentVariable (const int i=-1, const int j=-1) const override
	{ (void)i; (void)j; return {m_x}; }

      /**
       *  @brief get data for Data1D
       *  @param [out] data vector containing the dataset
       *  @return none
       */
      void get_data (std::vector<double> &data) const override
      { data = m_data; }

      /**
       *  @brief get standard deviation for Data1D
       *  @param [out] error vector containing the standard deviation
       *  @return none
       */
      void get_error (std::vector<double> &error) const override
      { error = m_error; }

      ///@}

      
      /**
       *  @name Member functions used to set the private members
       */
      ///@{

      /**
       *  @brief set interval variable m_x
       *  @param x vector containing x points
       *  @return none
       */
      void set_xx (const std::vector<double> x);

      ///@}

      
      /**
       *  @name Member functions for Input/Output 
       */
      ///@{

      /**
       *  @brief read the data
       *
       *  @param input_file input data file
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
      void read (const std::string input_file, const int skip_nlines=0, const int column_x=1, const std::vector<int> column_data={}, const std::vector<int> column_errors={}) override;

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
       *
       *  @param dir output directory
       *
       *  @param file output file
       *
       *  @param header text with the variable names to be written at
       *  the first line of the output file
       *
       *  @param prec the float precision
       *
       *  @param ww number of characters to be used as field width
       *
       *  @param rank cpu index (for MPI usage)
       *
       *  @return none
       */
      void write (const std::string dir, const std::string file, const std::string header, const int prec=4, const int ww=8, const int rank=0) const override;
      
      /**
       *  @brief write the covariance
       *
       *  @param dir the output directory
       *
       *  @param file the output file
       *
       *  @param precision the float precision
       *
       *  @return none
       */
      void write_covariance (const std::string dir, const std::string file, const int precision=10) const override;
      
      ///@}
      

      /**
       *  @name Member functions for data cut
       */

      ///@{
      
      /**
       * @brief cut the data
       * @param xmin minumum value for the independent variable x
       * @param xmax maximum value for the independent variable x
       * @return pointer to an object of type Data1D
       */
      std::shared_ptr<Data> cut (const double xmin, const double xmax) const override;

      /**
       * @brief cut the data, for Data1D
       * @param [in] mask vector containing values to be masked
       * @return pointer to an object of type Data1D
       */
      std::shared_ptr<Data> cut (const std::vector<bool> mask) const;

      ///@}

    };
    
  }
}

#endif
