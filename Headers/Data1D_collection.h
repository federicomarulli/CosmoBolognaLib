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
 *  @file Headers/Data1D_collection.h
 *
 *  @brief The class Data1D_collection
 *
 *  This file defines the interface of the class Data1D_collection
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __DATA1DC__
#define __DATA1DC__

#include "Data1D.h"

namespace cbl {

  namespace data {
  
    /**
     *  @class Data1D_collection Data1D_collection.h
     *  "Headers/Data1D_collection.h"
     *
     *  @brief The class Data1D_collection
     *
     *  This is the base class used to manage collections of 1D data
     */
    class Data1D_collection : public Data
    {
    protected:

      /**
       *  @name Data input
       */
      ///@{
	
      /// The number of datasets
      int m_ndataset;

      /// vector containing the number of data in each dataset
      std::vector<int> m_xsize;

      /// vector containing indipendent variables
      std::vector<std::vector<double>> m_x;

      /// vector containing matrix-to-vector indexes
      std::vector<std::vector<int>> m_index;

      ///@}
      
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *
       *  @return object of class Data1D
       */
      Data1D_collection () : Data() {set_dataType(DataType::_1D_collection_);}

      /**
       *  @brief constructor which reads the data from file
       *
       *  @param input_file file with datasets
       *
       *  @param skip_header number of lines to be skipped 
       *
       *  @return object of class Data1D_collection
       */
      Data1D_collection (const std::string input_file, const int skip_header=0);

      /**
       *  @brief constructor which reads the data from a vector of
       *  input files
       *
       *  @param input_files files with datasets
       *
       *  @param skip_header number of lines to be skipped 
       *
       *  @return object of class Data1D_collection
       */
      Data1D_collection (const std::vector<std::string> input_files, const int skip_header=0);

      /**
       *  @brief constructor which gets the data from input matrices
       *
       *  @param xx vector containing the independent variable xx for
       *  each dataset
       *
       *  @param data vector containing the data for each dataset
       *
       *  @return object of class Data1D_collection
       */
      Data1D_collection (const std::vector<std::vector<double>> xx, const std::vector<std::vector<double>> data);

      /**
       *  @brief constructor which gets both the data and the
       *  covariance matrix from input matrices
       *
       *  @param xx vector containing the independent variable xx 
       *  for each dataset
       *
       *  @param data vector containing the data for each dataset
       *
       *  @param covariance vector containing the data covariance
       *
       *  @return object of class Data1D_collection
       */
      Data1D_collection (const std::vector<std::vector<double>> xx, const std::vector<std::vector<double>> data, const std::vector<std::vector<double>> covariance);

      /**
       *  @brief constructor which gets both the data and the
       *  errors from input vectors
       *
       *  @param xx vector containing the independent variable xx 
       *  for each dataset
       *
       *  @param data vector containing the data for each dataset
       *
       *  @param error vector containing the data error
       *
       *  @return object of class Data1D_collection
       */
      Data1D_collection (const std::vector<std::vector<double>> xx, const std::vector<double> data, const std::vector<double> error);

      /**
       *  @brief constructor which gets both the data and the
       *  covariance matrix from input vectors
       *
       *  @param xx vector containing the independent variable xx 
       *  for each dataset
       *
       *  @param data vector containing the data for each dataset
       *
       *  @param covariance vector containing the data covariance
       *
       *  @return object of class Data1D_collection
       */
      Data1D_collection (const std::vector<std::vector<double>> xx, const std::vector<double> data, const std::vector<std::vector<double>> covariance);

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      virtual ~Data1D_collection () = default;

      /**
       *  @brief static factory used to construct objects of class
       *  Data1D
       *
       *  @return a shared pointer to an object of class Data
       */
      std::shared_ptr<Data> as_factory () { return move(std::unique_ptr<Data1D_collection>(this)); }

      ///@}

      
      /**
       *  @name Member functions used to get the private members
       */
      ///@{
    
      int xsize (const int i) const
      { return m_xsize[i]; }

      /**
       *  @brief get value of x at index j in the i-th dataset
       *
       *  @param i the i-th dataset
       *
       *  @param j index
       *
       *  @return the value of the x vector at position j in the i-th dataset
       */
      double xx (const int i, const int j) const { return m_x[i][j]; }  

      /**
       *  @brief get x values
       *
       *  @param x [out] the x values
       *
       *  @return the value of the m_x vector for the i-th dataset
       */
      void xx (std::vector<std::vector<double> > x) const { x = m_x; }

      /**
       *  @brief get the independet variable, to be used 
       *  in model computation
       *
       *  @param i index (not used)
       *
       *  @param j index (not used)
       *
       *  @return the independent variable
       */
      std::vector<std::vector<double>> IndipendentVariable(const int i=-1, const int j=-1) const {(void)i; (void)j; return m_x;}

      /**
       *  @brief get data at index i,j
       *
       *  @param i index
       *
       *  @param j index
       *
       *  @return the value of the m_data vector at position i,j
       */
      double data (const int i, const int j) const {return m_data[m_index[i][j]];} 

      /**
       *  @brief get data
       *
       *  @param [out] data vector containing the dataset
       *
       *  @return none
       */
      void data (std::vector<std::vector<double>> &data) const;

      /**
       *  @brief get value of f(x) error at index i,j
       *
       *  @param i index
       *
       *  @param j index
       *
       *  @return the value of the m_error vector at position i,j
       */
      double error (const int i, const int j) const { return m_error[m_index[i][j]]; } 

      /**
       *  @brief get standard deviation
       *
       *  @param [out] error vector containing the error
       *
       *  @return none
       */
      void error (std::vector<std::vector<double>> &error) const;

      ///@}
      
      /**
       *  @name Member functions used to set the private members
       */
      ///@{

      /**
       *  @brief set interval variable m_x in the i-th dataset
       *
       *  @param i index to the i-th dataset
       *
       *  @param x vector containing x points
       *
       *  @return none
       */
      void set_xx (const int i, const std::vector<double> x) { m_x[i] = x; } 

      ///@}
      
      /**
       *  @name Member functions to compute data properties
       */
      ///@{

      /**
       * @brief function that returns total number of datasets
       * @return total number of dataset
       */
      int ndataset () const { return m_ndataset; }  
      
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
       *  second column will be used by default, assuming that only 1
       *  data vector has to be read
       *
       *  @param column_errors the column of error values in the input
       *  file; the size of column_error must be equal to the size of
       *  column_data; if the size of column_error is larger than 1,
       *  more than 1 error vectors are read and then added one after
       *  the other in a single data object; if column_random is not
       *  provided, the third column will be used by default, assuming
       *  that only 1 random vector has to be read; if the input file
       *  has only 2 columns, the errors will be set to 1
       *
       *  @return none
       *
       *  @warning column_x, column_data, column_errors are not used
       *  in the current implementation: work in progress!
       */
      void read (const std::string input_file, const int skip_nlines=0, const int column_x=0, const std::vector<int> column_data={}, const std::vector<int> column_errors={}) override;
      
      /**
       *  @brief read the data
       *
       *  @param input_files input data files
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
       *  second column will be used by default, assuming that only 1
       *  data vector has to be read
       *
       *  @param column_errors the column of error values in the input
       *  file; the size of column_error must be equal to the size of
       *  column_data; if the size of column_error is larger than 1,
       *  more than 1 error vectors are read and then added one after
       *  the other in a single data object; if column_random is not
       *  provided, the third column will be used by default, assuming
       *  that only 1 random vector has to be read; if the input file
       *  has only 2 columns, the errors will be set to 1
       *
       *  @return none
       *
       *  @warning column_x, column_data, column_errors are not used
       *  in the current implementation: work in progress!
       */
      void read (const std::vector<std::string> input_files, const int skip_nlines=0, const int column_x=0, const std::vector<int> column_data={}, const std::vector<int> column_errors={}) override;

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
       *  @brief write the data
       *
       *  @param dir output directory
       *
       *  @param files output files
       *
       *  @param header text with the variable names to be written at
       *  the first line of the output file
       *
       *  @param precision the float precision
       *
       *  @param ww number of characters to be used as field width 
       *
       *  @param rank cpu index (for MPI usage)
       *
       *  @return none
       */
      void write (const std::string dir, const std::vector<std::string> files, const std::string header, const int precision=10, const int ww=8, const int rank=0) const override;

      /**
       *  @brief write the interval variable m_covariance on a file,
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
       * @brief cut the data, for Data1D_collection
       * @param dataset the i-th dataset
       * @param xmin minumum value for the independet variable x
       * @param xmax maximum value for the independent variable x
       * @return pointer to an object of type Data1D
       */
      std::shared_ptr<Data> cut(const int dataset, const double xmin, const double xmax) const override;

      /**
       * @brief cut the data, for Data1D_collection type
       * @param xmin vector containing minumum values for the independet variable x
       * @param xmax vector containing maximum values for the independent variable x
       * @return pointer to an object of type Data1D_collection
       */
      std::shared_ptr<Data> cut(const std::vector<double> xmin, const std::vector<double> xmax) const override;

      ///@}

    };
   
  }
}

#endif
