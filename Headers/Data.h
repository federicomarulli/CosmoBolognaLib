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
 *  @file Headers/Data.h
 *
 *  @brief The class Data
 *
 *  This file defines the interface of the class Data
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __DATA__
#define __DATA__

#include "Func.h"

namespace cbl {

  /**
   *  @brief The namespace of the functions and classes used to handle
   *  <B> data </B> of any kind
   *  
   *  The \e data namespace contains all the main functions and
   *  classes to handle data of any kind
   */
  namespace data {

    /**
     *  @enum DataType 
     *  @brief the data type
     */
    enum class DataType { 

      /// 1D dataset
      _1D_,

      /// 2D dataset
      _2D_,
     
      /// collection of 1D datasets
      _1D_collection_,

      /// 1D dataset with extra information
      _1D_extra_,

      /// 2D dataset with extra information
      _2D_extra_,

    };

    /**
     * @brief return a vector containing the
     * DataType names
     * @return a vector containing the
     * DataType names
     */
    inline std::vector<std::string> DataTypeNames () { return {"1D", "2D", "1D_collection", "1D_extra", "2D_extra"}; }
  
    /**
     *  @class Data Data.h
     *  "Headers/Data.h"
     *
     *  @brief The class Data
     *
     *  This is the base class used to manage 1D-2D data
     */
    class Data
    {
      
    protected:

      /// type of data
      DataType m_dataType;

      /// number of data
      int m_ndata;

      /// data values
      std::vector<double> m_data;
      
      /// standard deviations
      std::vector<double> m_error;
      
      /// bin edges for the x variable
      std::vector<double> m_edges_xx;
      
      /// bin edges for the y variable
      std::vector<double> m_edges_yy;
      
      /// covariance matrix
      std::vector<std::vector<double>> m_covariance;
      
      /// inverse covariance matrix
      std::vector<std::vector<double>> m_inverse_covariance;

      /**
       *  @brief set the data type 
       *  @param dataType the data type
       */
      void set_dataType (const DataType dataType) { m_dataType = dataType; } 

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *
       *  
       */
      Data () = default;
      
      /**
       *  @brief constructor which sets only the data type
       *
       *  @param dataType the data type
       *
       *  
       */
      Data (const DataType dataType) : m_dataType(dataType) {}

      /**
       *  @brief constructor which sets both the data type and the
       *  number of data
       *
       *  @param dataType the data type
       *
       *  @param ndata the number of data
       *
       *  
       */
      Data (const DataType dataType, const int ndata)
	: m_dataType(dataType)
      { reset(ndata); }

      /**
       *  @brief constructor which gets the data from an input vector
       *
       *  @param dataType the data type
       *
       *  @param data vector containing the data
       *
       *  
       */
      Data (const DataType dataType, const std::vector<double> data); 

      /**
       *  @brief constructor which gets both the data and the errors
       *  from input vectors
       *
       *  @param dataType the data type
       *
       *  @param data vector containing the data
       *
       *  @param error vector containing the errors
       *
       *  
       */
      Data (const DataType dataType, const std::vector<double> data, const std::vector<double> error);

      /**
       *  @brief constructor which gets both the data and the
       *  covariance matrix from input vectors
       *
       *  @param dataType the data type
       *
       *  @param data vector containing the data
       *
       *  @param covariance matrix containing the covariance
       *
       *  
       */
      Data (const DataType dataType, const std::vector<double> data, const std::vector<std::vector<double>> covariance);

      /**
       *  @brief default destructor
       */
      virtual ~Data () = default;

      /**
       *  @brief static factory used to construct objects of class
       *  Data1D
       *  @param dataType the data type
       *  @return a shared pointer to an object of class Data
       */
      static std::shared_ptr<Data> Create (const DataType dataType); 

      /**
       *  @brief static factory used to construct objects of class
       *  Data1D
       *  @return a shared pointer to an object of class Data
       */
      virtual std::shared_ptr<Data> as_factory () 
      { ErrorCBL("", "as_factory", "Data.h"); return NULL; }

      ///@}

      
      /**
       *  @name Non-virtual member functions
       */
      ///@{
      
      /**
       *  @brief get the value of the data correlation at index i,j
       *  @param i index
       *  @param j index
       *  @return the value of the correlation
       * \f$ Corr_{i,j} = \frac{Cov_{i,j}}{\sqrt{Cov_{i,i} \cdot Cov{j,j}}} \f$
       */
      double correlation (const int i, const int j) const { return m_covariance[i][j]/sqrt(m_covariance[i][i]*m_covariance[j][j]); }

      /**
       *  @brief get the value of the data correlation at index i,j
       *  @return the value of the correlation
       * \f$ Corr_{i,j} = \frac{Cov_{i,j}}{\sqrt{Cov_{i,i} \cdot Cov{j,j}}} \f$
       */
      std::vector<std::vector<double>> correlation () const;
      
      /**
       *  @brief get the value of data inverse_covariance at index i,j
       *  @param i index
       *  @param j index
       *  @return the value of the m_inverse_covariance matrix at position i,j
       */
      double inverse_covariance (const int i, const int j) const { return m_inverse_covariance[i][j]; }

      /**
       *  @brief get the m_inverse_covariance std::vector
       *  @return the std::vector containing the inverse convariance matrix
       */
      std::vector<std::vector<double>> inverse_covariance () const { return m_inverse_covariance; }

      /**
       *  @brief reset data object with new empty arrays large enough
       *  to store ndata data
       *
       *  @param ndata the new number of data
       */
      void reset (const int ndata);

      /**
       *  @brief invert the covariance matrix
       *
       *  @param prec the precision required in the inversion of the
       *  covariance matrix
       *
       *  @param Nres \f$N_{res}\f$, the number of catalogue
       *  resamplings used to estimate the covariance matrix;
       *  \f$N_{res}=-1\f$ if the covariance matrix has not been
       *  estimated with resampling methods
       *
       *  @return none
       */
      virtual void invert_covariance (const double prec, const int Nres) { invert_matrix(m_covariance, m_inverse_covariance, prec, Nres); }

      ///@}
      

      /**
       *  @name Member functions to get the private/protected members
       */
      ///@{

      /**
       *  @brief get value of x at index i
       *  @param i index
       *  @return the value of the m_x std::vector at position i
       */
      virtual double xx (const int i) const
      { (void)i; ErrorCBL("", "xx", "Data.h"); return 0; } 

      /**
       *  @brief get the x std::vector
       *  @return the x std::vector
       */
      virtual std::vector<double> xx () const 
      { ErrorCBL("", "xx", "Data.h"); return {}; } 

      /**
       *  @brief get value of x at position i,j, for Data1D_collection
       *  @param i index
       *  @param j index
       *  @return the value of the m_x std::vector at position i,j
       */
      virtual double xx (const int i, const int j) const
      { (void)i; (void)j; ErrorCBL("", "xx", "Data.h"); return 0.; }

      /**
       *  @brief get value of y at index i, for Data2D
       *  @param i index
       *  @return the value of the m_y std::vector at position i
       */
      virtual double yy (const int i) const 
      { (void)i; ErrorCBL("", "yy", "Data.h"); return 0.; }

      /**
       *  @brief get the y std::vector
       *  @return the y std::vector
       */
      virtual std::vector<double> yy () const 
      { ErrorCBL("", "yy", "Data.h"); return {}; } 
      
      /**
       *  @brief get the independet variable, to be used 
       *  in model computation
       *  @param i first indipendent variable index
       *  @param j second indipendent variable index
       *  @return the independent variable
       */
      virtual std::vector<std::vector<double>> IndipendentVariable(const int i=-1, const int j=-1) const
      { (void)i; (void)j; ErrorCBL("", "IndipendentVariable", "Data.h"); std::vector<std::vector<double>> pp; return {pp}; }

      /**
       *  @brief get data at index i,j for Data1D_collection, Data2D
       *  @param i index
       *  @param j index
       *  @return the value of the m_data std::vector at position i,j
       */
      virtual double data (const int i, const int j) const
      { (void)i; (void)j; ErrorCBL("", "data", "Data.h"); return 0.; }

      /**
       *  @brief get data for Data1D
       *  @param [out] data std::vector containing the dataset
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void get_data (std::vector<double> &data) const
      { (void)data; ErrorCBL("", "data", "Data.h"); }
      
      /**
       *  @brief get the data type
       *  @return the data type
       */
      virtual DataType dataType () const
      { return m_dataType; }

      /**
       *   @brief the total number of data
       *   @return total data number
       */
      virtual int ndata () const { return m_ndata; }

      /**
       *  @brief get data at index i
       *  @param i index
       *  @return the value of the data vector at position i
       */
      virtual double data (const int i) const { return m_data[i]; }

      /**
       *  @brief get data
       *  @return the dataset
       */
      virtual std::vector<double> data () const { return m_data; }

      /**
       *  @brief get value of data standard deviation at index i
       *  @param i index
       *  @return the value of the error vector at position i
       */
      virtual double error (const int i) const { return m_error[i]; }

      /**
       *  @brief get standard deviation
       *  @return the standard deviation
       */
      virtual std::vector<double> error () const { return m_error; }

      /**
       *  @brief get the value of the data covariance at index i,j
       *  @param i index
       *  @param j index
       *  @return the value of the m_covariance matrix at position i,j
       */
      virtual double covariance (const int i, const int j) const { return m_covariance[i][j]; }

      /**
       *  @brief get the m_covariance vector
       *  @return the vector containing the covariance matrix
       */
      virtual std::vector<std::vector<double>> covariance () const { return m_covariance; }
      
      /**
       *  @brief get value of x variable bin edge at index i
       *  @param i index
       *  @return the value of the x variable bin edge vector at position i
       */
      virtual double edges_xx (const int i) const { return m_edges_xx[i]; }

      /**
       *  @brief get x variable bin edges
       *  @return the x variable bin edges
       */
      virtual std::vector<double> edges_xx () const { return m_edges_xx; }
      
      /**
       *  @brief get value of y variable bin edge at index i
       *  @param i index
       *  @return the value of the y variable bin edge vector at position i
       */
      virtual double edges_yy (const int i) const { return m_edges_yy[i]; }

      /**
       *  @brief get y variable bin edges
       *  @return the y variable bin edges
       */
      virtual std::vector<double> edges_yy () const { return m_edges_yy; }

      /**
       *  @brief get data for Data1D_collection, Data2D
       *  @param [out] data std::vector containing the dataset
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void get_data (std::vector<std::vector<double>> &data) const
      { (void)data; ErrorCBL("", "data", "Data.h"); }

      /**
       *  @brief get value of f(x) error at index i,j for Data1D_collection, Data2D
       *  @param i index
       *  @param j index
       *  @return the value of the m_error std::vector at position i,j
       */
      virtual double error (const int i, const int j) const
      { (void)i; (void)j; ErrorCBL("", "error", "Data.h"); return 0.; } 

      /**
       *  @brief get standard deviation for Data1D
       *  @param [out] error std::vector containing the staandard deviation
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void get_error (std::vector<double> &error) const
      { (void)error; ErrorCBL("", "error", "Data.h"); }

      /**
       *  @brief get standard deviation for Data1D_Collection, Data2D
       *  @param [out] error std::vector containing the staandard deviation
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void get_error (std::vector<std::vector<double>> &error) const
      { (void)error; ErrorCBL("", "error", "Data.h"); }

       /**
       *  @brief return the value of the extra information at index i,j
       *  @param i index
       *  @param j index
       *  @return the value of the extra_info at position i,j
       */
      virtual double extra_info (const int i, const int j) const
      { (void)i; (void)j; ErrorCBL("", "extra_info", "Data.h"); return 0.; }

      /**
       *  @brief return the m_exta_info std::vector
       *  @return std::vector containing std::vectors with extra information
       */
      virtual std::vector<std::vector<double>> extra_info () const
      { ErrorCBL("", "extra_info", "Data.h"); std::vector<std::vector<double>> x; return x; }
      
      ///@}

      /**
       *  @name Member functions to set the private/protected members
       */
      ///@{

      /**
       *  @brief set interval variable m_x
       *  @param x std::vector containing x points
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_xx (const std::vector<double> x)
      { (void)x; ErrorCBL("", "set_xx", "Data.h"); }

      /**
       *  @brief set interval variable m_y, for Data2D
       *  @param y std::vector containing y points
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_yy (const std::vector<double> y)
      { (void)y; ErrorCBL("", "set_yy", "Data.h"); }

      /**
       *  @brief set interval variable m_x in the i-th dataset,
       *  for Data1D_collection
       *  @param i index to the i-th dataset
       *  @param x std::vector containing x points
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_xx (const int i, const std::vector<double> x) 
      { (void)i; (void)x; ErrorCBL("", "set_xx", "Data.h"); }

      /**
       *  @brief set interval variable m_x, for Data1D_collection
       *  @param x std::vector containing x points
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_xx (const std::vector<std::vector<double>> x)
      { (void)x; ErrorCBL("", "set_xx", "Data.h"); }
      
      /**
       *  @brief set interval variable m_edges_xx
       *  @param edges std::vector containing the x bin edges
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_edges_xx (const std::vector<double> edges)
      { (void)edges; ErrorCBL("", "set_edges_xx", "Data.h"); }
      
      /**
       *  @brief set interval variable m_edges_yy
       *  @param edges std::vector containing the y bin edges
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_edges_yy (const std::vector<double> edges)
      { (void)edges; ErrorCBL("", "set_edges_yy", "Data.h"); }

      /**
       *  @brief set interval variable m_data, for Data1D_collection,
       *  Data2D
       *  @param data std::vector containing data points 
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_data (const std::vector<std::vector<double>> data) 
      { (void)data; ErrorCBL("", "set_data", "Data.h"); }

      /**
       *  @brief set interval variable m_error_fx
       *  @param extra_info std::vector containing std::vectors with extra
       *  information
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_extra_info (const std::vector<std::vector<double>> extra_info)
      { (void)extra_info; ErrorCBL("", "set_extra_info", "Data.h"); }

      /**
       *  @brief set interval variable data
       *  @param data std::vector containing data points 
       */
      void set_data (const std::vector<double> data);

      /**
       *  @brief set interval variable m_error_fx
       *  @param error std::vector containing data standard deviation
       */
      void set_error (const std::vector<double> error);

      /**
       *  @brief set interval variable m_error_fx
       *  @param covariance std::vector containing the covariance matrix
       */
      void set_error (const std::vector<std::vector<double>> covariance);

      /**
       *  @brief set the interval variable m_covariance, reading from
       *  an input file
       *
       *  @param filename file containing the covariance matrix in the
       *  format: column 0 \f$ \rightarrow \f$ x<SUB>i</SUB>, column 1
       *  \f$ \rightarrow \f$ x<SUB>j</SUB>, column cov_col &rarr;
       *  cov(x<SUB>i</SUB>,x<SUB>j</SUB>)
       *
       *  @param cov_col covariance matrix column, starting from 0
       *
       *  @param skipped_lines comment lines to be skipped
       */
      void set_covariance (const std::string filename, const int cov_col=2, const int skipped_lines=0);

      /**
       *  @brief set interval the variable m_covariance
       *  @param value covariance matrix value
       *  @param i the first index 
       *  @param j the second index 
       */
      void set_covariance (const double value, const int i, const int j) { m_covariance[i][j] = value; }

      /**
       *  @brief set interval the variable m_covariance
       *  @param covariance std::vector containing the covariance matrix
       */
      void set_covariance (const std::vector<std::vector<double>> covariance);

      /**
       *  @brief set interval the variable m_covariance
       *  @param error std::vector containing the data standard deviation
       */
      void set_covariance (const std::vector<double> error);
      
      ///@}

      /**
       *  @name Member functions to compute data properties
       */
      ///@{

      /**
       * @brief function that returns number of data for one 
       * dataset, for Data1D_collection
       * @param i index to the i-th dataset
       * @return total number of data
       */
      virtual int ndata (const int i) const 
      { (void)i; return ErrorCBL("", "ndata", "Data.h"); }

      /**
       * @brief function that returns total number of datasets
       * @return total number of dataset
       */
      virtual int ndataset () const
      { return ErrorCBL("", "ndataset", "Data.h"); }

      /**
       *  @brief get the number of points along x
       *  @return the number of points along x
       */
      virtual int xsize () const
      { return ErrorCBL("", "xsize", "Data.h"); }

      /**
       *  @brief get the number of points along x
       *  @param i index to the i-th dataset
       *  @return the number of points along x
       */
      virtual int xsize (const int i) const
      { (void)i; return ErrorCBL("", "xsize", "Data.h"); }

      /**
       *  @brief get the number of points along y
       *  @return the number of points along y
       */
      virtual int ysize () const
      { return ErrorCBL("", "ysize", "Data.h"); }

      ///@}


      /**
       *  @name Member functions for Input/Output
       */
      ///@{
      
      /**
       *  @brief read the data from one input file
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
       *  first column after the column of x values will be used by
       *  default, assuming that only 1 data vector has to be read
       *
       *  @param column_errors the column of error values in the input
       *  file; the size of column_error must be equal to the size of
       *  column_data; if the size of column_error is larger than 1,
       *  more than 1 error vectors are read and then added one after
       *  the other in a single data object; if column_random is not
       *  provided, the second column after the column of x values
       *  will be used by default, assuming that only 1 random vector
       *  has to be read; if the input file has only 2 columns, the
       *  errors will be set to 1
       *
       *  @param column_edges vector containing the columns of x (and
       *  y) bin edge values in the input file; if it is not provided,
       *  the third and four columns after the column of x values will
       *  be used; if these columns do no exist the edges are not read
       *
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void read (const std::string input_file, const int skip_nlines=0, const std::vector<int> column={0}, const std::vector<int> column_data={}, const std::vector<int> column_errors={}, const std::vector<int> column_edges={})
      { (void)input_file; (void)skip_nlines; (void)column, (void)column_data, (void)column_errors, (void)column_edges, ErrorCBL("", "read", "Data.h"); }

      /**
       *  @brief read the data from a set of input files
       *
       *  @param input_file vector containing the input data files
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
       *  first column after the column of x values will be used by
       *  default, assuming that only 1 data vector has to be read
       *
       *  @param column_errors the column of error values in the input
       *  file; the size of column_error must be equal to the size of
       *  column_data; if the size of column_error is larger than 1,
       *  more than 1 error vectors are read and then added one after
       *  the other in a single data object; if column_random is not
       *  provided, the second column after the column of x values
       *  will be used by default, assuming that only 1 random vector
       *  has to be read; if the input file has only 2 columns, the
       *  errors will be set to 1
       *
       *  @param column_edges vector containing the columns of x (and
       *  y) bin edge values in the input file; if it is not provided,
       *  the third and four columns after the column of x values will
       *  be used; if these columns do no exist the edges are not read
       *
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void read (const std::vector<std::string> input_file, const int skip_nlines=0, const std::vector<int> column={0}, const std::vector<int> column_data={}, const std::vector<int> column_errors={}, const std::vector<int> column_edges={})
      { (void)input_file; (void)skip_nlines; (void)column, (void)column_data, (void)column_errors, (void)column_edges, ErrorCBL("", "read", "Data.h"); }
      
      /**
       *  @brief print the data on screen
       *
       *  @param prec the float precision
       *
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void Print (const int prec=4) const 
      { (void)prec; ErrorCBL("", "Print", "Data.h"); }
      
      /**
       *  @brief write the data
       *  @param dir output directory
       *  @param file output file
       *  @param header text with the variable names to be written at
       *  the first line of the output file
       *  @param prec the float precision
       *  @param ww number of characters to be used as field width
       *  @param rank cpu index (for MPI usage)
       *
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void write (const std::string dir, const std::string file, const std::string header, const int prec=4, const int ww=8, const int rank=0) const 
      { (void)dir; (void)file; (void)header; (void)prec; (void)ww; (void)rank; ErrorCBL("", "write", "Data.h"); }

      /**
       *  @brief write the data
       *  @param dir output directory
       *  @param file output file
       *  @param header text with the variable names to be written at
       *  the first line of the output file
       *  @param full false &rarr; simply store the data; true &rarr;
       *  duplicate the data in the other three quadrands (usefull
       *  e.g. when storing the 2D correlation function)
       *  @param prec the float precision
       *  @param ww number of characters to be used as field width
       *  @param rank cpu index (for MPI usage)
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void write (const std::string dir, const std::string file, const std::string header, const bool full, const int prec=10, const int ww=8, const int rank=0) const
      { (void)dir; (void)file; (void)header; (void)full; (void)prec; (void)ww; (void)rank; ErrorCBL("", "write", "Data.h"); }

      /**
       *  @brief write the data
       *  @param dir output directory
       *  @param files output file
       *  @param header text with the variable names to be written at
       *  the first line of the output file
       *  @param prec the float precision
       *  @param ww number of characters to be used as field width
       *  @param rank cpu index (for MPI usage)
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void write (const std::string dir, const std::vector<std::string> files, const std::string header, const int prec=10, const int ww=8, const int rank=0) const
      { (void)dir; (void)files; (void)header; (void)prec; (void)ww; (void)rank; ErrorCBL("", "write", "Data.h"); }
      
      /**
       *  @brief write the interval variable m_covariance on a file,
       *  @param dir the output directory
       *  @param file the output file
       *  @param precision the float precision
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void write_covariance (const std::string dir, const std::string file, const int precision=10) const
      { (void)dir; (void)file; (void)precision; ErrorCBL("", "write_covariance", "Data.h"); }

      ///@}

      
      /**
       *  @name Member functions for data cut
       */

      ///@{

      /**
       *  @brief cut the dataset using a mask
       *  @param [in] mask std::vector containing values to be masked
       *  @param [out] data std::vector containing data
       *  @param [out] error std::vector containing data standard deviations
       *  @param [out] covariance_matrix std::vector containing data covariance matrix
       */
      void cut (const std::vector<bool> mask, std::vector<double> &data, std::vector<double> &error, std::vector<std::vector<double>> &covariance_matrix) const;

      /**
       * @brief cut the data, for Data1D
       * @param [in] mask std::vector containing values to be masked
       * @return pointer to an object of type Data1D
       */
      virtual std::shared_ptr<Data> cut (const std::vector<bool> mask) const
      { (void)mask; ErrorCBL("", "cut", "Data.h"); std::shared_ptr<Data> dd; return dd; } 

      /**
       * @brief cut the data, for Data1D
       * @param xmin minumum value for the independet variable x
       * @param xmax maximum value for the independent variable x
       * @return pointer to an object of type Data1D
       */
      virtual std::shared_ptr<Data> cut (const double xmin, const double xmax) const
      { (void)xmin; (void)xmax; ErrorCBL("", "cut", "Data.h"); std::shared_ptr<Data> dd; return dd; } 

      /**
       * @brief cut the data, for Data2D
       * @param xmin minumum value for the independent variable x
       * @param xmax maximum value for the independent variable x
       * @param ymin minumum value for the independent variable y
       * @param ymax maximum value for the independent variable y
       * @return pointer to an object of type Data2D
       */
      virtual std::shared_ptr<Data> cut (const double xmin, const double xmax, const double ymin, const double ymax) const
      { (void)xmin; (void)xmax; (void)ymin; (void)ymax; ErrorCBL("", "cut", "Data.h"); std::shared_ptr<Data> dd; return dd; } 

      /**
       * @brief cut the data, for Data1D_collection
       * @param dataset the dataset index
       * @param xmin minumum value for the independet variable x
       * @param xmax maximum value for the independent variable x
       * @return pointer to an object of type Data1D
       */
      virtual std::shared_ptr<Data> cut (const int dataset, const double xmin, const double xmax) const
      { (void)dataset; (void)xmin; (void)xmax; ErrorCBL("", "cut", "Data.h"); std::shared_ptr<Data> dd; return dd; } 

      /**
       * @brief cut the data, for Data1D_collection type
       * @param xmin std::vector containing minumum values for the independet variable x
       * @param xmax std::vector containing maximum values for the independent variable x
       * @return pointer to an object of type Data1D_collection
       */
      virtual std::shared_ptr<Data> cut (const std::vector<double> xmin, const std::vector<double> xmax) const
      { (void)xmin; (void)xmax; ErrorCBL("", "cut", "Data.h"); std::shared_ptr<Data> dd; return dd; } 

      ///@}
      
    };
    

    /**
     *  @brief merge dataset (only work for one dataset type)
     *  @param dataset std::vector containing the dataset to merge
     *  @return pointer to an object of class Data
     */
    std::shared_ptr<data::Data> join_dataset (std::vector<std::shared_ptr<data::Data>> dataset);	

    /**
     *  @brief merge dataset of type _1D_
     *  @param dataset std::vector containing the dataset to merge
     *  @return pointer to an object of class Data
     */
    std::shared_ptr<data::Data> join_dataset_1D (std::vector<std::shared_ptr<data::Data>> dataset);	

    /**
     *  @brief merge dataset of type _1D_extra_
     *  @param dataset std::vector containing the dataset to merge
     *  @return pointer to an object of class Data
     */
    std::shared_ptr<data::Data> join_dataset_1D_extra (std::vector<std::shared_ptr<data::Data>> dataset);	

  }
}

#endif
