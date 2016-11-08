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
 *  @file Headers/Lib/Data.h
 *
 *  @brief The class Data
 *
 *  This file defines the interface of the class Data
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __DATA__
#define __DATA__

#include "GSLfunction.h"


namespace cosmobl {

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
    enum DataType { 

      /// 1D dataset
      _1D_data_,

      /// 2D dataset
      _2D_data_,
     
      /// collection of 1D datasets
      _1D_data_collection_,

      /// 1D dataset with extra information
      _1D_data_extra_,

      /// 2D dataset with extra information
      _2D_data_extra_,

    };

  
    /**
     *  @class Data Data.h
     *  "Headers/Lib/Data.h"
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

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return an object of class Data
       */
      Data () = default;

      /**
       *  @brief default constructor
       *  @param dataType the data type
       *  @return an object of class Data
       */
      Data (const DataType dataType) : m_dataType(dataType) {}

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Data () = default;

      /**
       *  @brief static factory used to construct objects of class
       *  Data1D
       *  @param dataType the data type
       *  @return a shared pointer to an object of class Data
       */
      static shared_ptr<Data> Create (const DataType dataType); 

      /**
       *  @brief static factory used to construct objects of class
       *  Data1D
       *  @param input_file input data file 
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @return a shared pointer to an object of class Data
       */
      static shared_ptr<Data> Create (const string input_file, const double xmin=par::defaultDouble, const double xmax=-par::defaultDouble); 

      /**
       *  @brief static factory used to construct objects of class
       *  Data1D
       *  @param x vector containing x points 
       *  @param fx vector containing f(x) 
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @return a shared pointer to an object of class Data
       */
      static shared_ptr<Data> Create (const vector<double> x, const vector<double> fx, const double xmin=par::defaultDouble, const double xmax=-par::defaultDouble); 

      /**
       *  @brief static factory used to construct objects of class
       *  Data1D
       *  @param x vector containing x points 
       *  @param fx vector containing f(x) 
       *  @param error_fx vector containing error on f(x) 
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @return a shared pointer to an object of class Data
       */
      static shared_ptr<Data> Create (const vector<double> x, const vector<double> fx, const vector<double> error_fx, const double xmin=par::defaultDouble, const double xmax=-par::defaultDouble); 

      /**
       *  @brief static factory used to construct objects of class
       *  Data1D
       *  @param x vector containing x points 
       *  @param fx vector containing f(x) 
       *  @param covariance vector containing f(x) covariance matrix 
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @return a shared pointer to an object of class Data
       */
      static shared_ptr<Data> Create (const vector<double> x, const vector<double> fx, const vector<vector<double>> covariance, const double xmin=par::defaultDouble, const double xmax=-par::defaultDouble);

      /**
       *  @brief static factory used to construct objects of class
       *  Data2D
       *  @param x vector containing x points 
       *  @param y vector containing y points 
       *  @param fxy vector containing f(x,y) values
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @param ymin maximun value of y to be used 
       *  @param ymax maximun value of y to be used 
       *  @return a shared pointer to an object of class Data
       */
      static shared_ptr<Data> Create (const vector<double> x, const vector<double> y, const vector<vector<double>> fxy, const double xmin, const double xmax, const double ymin, const double ymax); 

      /**
       *  @brief static factory used to construct objects of class
       *  Data2D
       *  @param x vector containing x points 
       *  @param y vector containing y points 
       *  @param fxy vector containing f(x,y) values
       *  @param error_fxy vector containing error on f(x,y)
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @param ymin maximun value of y to be used 
       *  @param ymax maximun value of y to be used 
       *  @return a shared pointer to an object of class Data
       */
      static shared_ptr<Data> Create (const vector<double> x, const vector<double> y, const vector<vector<double>> fxy, const vector<vector<double>> error_fxy, const double xmin=par::defaultDouble, const double xmax=-par::defaultDouble, const double ymin=par::defaultDouble, const double ymax=-par::defaultDouble);
    
      ///@}


      /**
       *  @name Member functions to get the private/protected members
       */
      ///@{
      
      /**
       *  @brief get the data type
       *  @return the data type
       */
      DataType dataType () const
      { return m_dataType; }

      /**
       *  @brief get the index of the first x used
       *  @return the index of the first x used
       */
      virtual int x_down () const
      { ErrorCBL("Error in x_down of Data.h!"); return 0; }

      /**
       *  @brief get the index of the last x used
       *  @return the index of the last x used
       */
      virtual int x_up () const
      { ErrorCBL("Error in x_up of Data.h!"); return 0; }

      /**
       *  @brief get the index of the first y used
       *  @return the index of the first y used
       */
      virtual int y_down () const
      { ErrorCBL("Error in y_down of Data.h!"); return 0; }

      /**
       *  @brief get the index of the last y used
       *  @return the index of the last y used
       */
      virtual int y_up () const
      { ErrorCBL("Error in y_up of Data.h!"); return 0; }

      /**
       *  @brief get value of x at index i
       *  @param i index
       *  @return the value of the m_x vector at position i
       */
      virtual double xx (const int i) const
      { (void)i; ErrorCBL("Error in xx of Data.h!"); return 0.; }

      /**
       *  @brief get value of y at index i
       *  @param i index
       *  @return the value of the m_y vector at position i
       */
      virtual double yy (const int i) const 
      { (void)i; ErrorCBL("Error in yy of Data.h!"); return 0.; }

      /**
       *  @brief get f(x) at index i
       *  @param i index
       *  @return the value of the m_fx vector at position i
       */
      virtual double fx (const int i) const
      { (void)i; ErrorCBL("Error in fx of Data.h!"); return 0.; }

      /**
       *  @brief get value of f(x) error at index i
       *  @param i index
       *  @return the value of the m_error_fx vector at position i
       */
      virtual double error_fx (const int i) const
      { (void)i; ErrorCBL("Error in error_fx of Data.h!"); return 0.; }
      
      /**
       *  @brief get the value of the f(x) covariance at index i,j
       *  @param i index
       *  @param j index
       *  @return the value of the m_covariance matrix at position i,j
       */
      virtual double covariance (const int i, const int j) const
      { (void)i; (void)j; ErrorCBL("Error in covariance of Data.h!"); return 0.; }

      /**
       *  @brief get covariance on f(x,y) at index i1, i2, j1, j2
       *  @param i1 index
       *  @param i2 index
       *  @param j1 index
       *  @param j2 index
       *  @return the value of the m_error_fxy vector at position (i1, i2), (j1, j2)
       */
      virtual double covariance (const int i1, const int i2, const int j1, const int j2) const
      { (void)i1; (void)j1; (void)i2; (void)j2; ErrorCBL("Error in covariance_fxy of Data.h!"); return 0.; }

      /**
       *  @brief get the value of f(x) inverse_covariance at index i,j
       *  @param i index
       *  @param j index
       *  @return the value of the m_inverse_covariance matrxi at position i,j
       */
      virtual double inverse_covariance (const int i, const int j) const
      { (void)i; (void)j; ErrorCBL("Error in inverse_covariance of Data.h!"); return 0.; }

      /**
       *  @brief get value of f(x) inverted covariance at index i,j
       *  @param d the d-th dataset
       *  @param i the i-th x element of the covariance matrix 
       *  @param j the j-th x element of the covariance matrix 
       *  @return the value of the inverted covariance matrix for at position i,j
       */
      virtual double inverse_covariance (const int d, const int i, const int j) const
      { (void)d; (void)i; (void)j; ErrorCBL("Error in inverse_covariance of Data.h!"); return 0.; }

      /**
       *  @brief get value of f(x,y) at index i,j
       *  @param i index
       *  @param j index
       *  @return the value of the m_fxy vector at position i,j
       */
      virtual double fxy (const int i, const int j) const
      { (void)i; (void)j; ErrorCBL("Error in fxy of Data.h!"); return 0.; }

      /**
       *  @brief get error on f(x,y) at index i,j
       *  @param i index
       *  @param j index
       *  @return the value of the m_error_fxy vector at position i,j
       */
      virtual double error_fxy (const int i, const int j) const
      { (void)i; (void)j; ErrorCBL("Error in error_fxy of Data.h!"); return 0.; }

      /**
       *  @brief the x vector
       *  @return the vector containing the x values
       */
      virtual vector<double> xx () const  
      { ErrorCBL("Error in xx of Data.h!"); vector<double> x; return x; }

      /**
       *  @brief get the m_fx vector
       *  @return the vector containing the fx values
       */
      virtual vector<double> fx () const 
      { ErrorCBL("Error in fx of Data.h!"); vector<double> x; return x; }

      /**
       *  @brief get the m_error_fx vector
       *  @return the vector containing the values of fx error
       */
      virtual vector<double> error_fx () const 
      { ErrorCBL("Error in error_fx of Data.h!"); vector<double> x; return x; }

      /**
       **  @brief get the m_covariance vector
       *  @return the vector<containing the covariance matrix
       */
      virtual vector<vector<double>> covariance () const
      { ErrorCBL("Error in covariance of Data.h!"); vector<vector<double>> x; return x; }

      /**
       *  @brief get the m_inverse_covariance vector
       *  @return the vector containing the inverse convariance matrix
       */
      virtual vector<vector<double>> inverse_covariance () const
      { ErrorCBL("Error in inverse_covariance of Data.h!"); vector<vector<double>> x; return x; }

      /**
       *  @brief get the y vector
       *  @return the vector containing the y values
       */
      virtual vector<double> yy () const  
      { ErrorCBL("Error in xx of Data.h!"); vector<double> x; return x; }

      /**
       *  @brief get the m_fx vector
       *  @return the vector containing the fx values
       */
      virtual vector<vector<double>> fxy () const 
      { ErrorCBL("Error in fxy of Data.h!"); vector<vector<double>> x; return x; }

      /**
       *  @brief get the m_error_fx vector
       *  @return the vector containing the values of fx errors
       */
      virtual vector<vector<double>> error_fxy () const 
      { ErrorCBL("Error in error_fxy of Data.h!"); vector<vector<double>> x; return x; }
    
      /**
       *  @brief index of the first x used in the i-th dataset
       *  @param i the i-th dataset
       *  @return the index of the first x used
       */
      virtual int x_down (const int i) const 
      { (void)i; ErrorCBL("Error in x_down of Data.h!"); return 0; }

      /**
       *  @brief index of the last x used in the i-th dataset
       *  @param i the i-th dataset
       *  @return the index of the last x used
       */
      virtual int x_up (const int i) const 
      { (void)i; ErrorCBL("Error in x_up of Data.h!"); return 0; }

      /**
       *  @brief the value of x at index j in the i-th dataset
       *  @param i the i-th dataset
       *  @param j index
       *  @return the value of the x[j] in the i-th dataset
       */
      virtual double xx (const int i, const int j) const  
      { (void)i; (void)j; ErrorCBL("Error in xx of Data.h!"); return 0.; }

      /**
       *  @brief the function f(x[j]) of the i-th dataset
       *  @param i the i-th dataset
       *  @param j index
       *  @return f(x[j]) in the i-th dataset
       */
      virtual double fx (const int i, const int j) const 
      { (void)i; (void)j; ErrorCBL("Error in fx of Data.h!"); return 0.; }

      /**
       *  @brief the error on f(x[j]) of the i-th dataset
       *  @param i the i-th dataset
       *  @param j index
       *  @return the error of f(x[j]) of the i-th dataset
       */
      virtual double error_fx (const int i, const int j) const 
      { (void)i; (void)j; ErrorCBL("Error in error_fx of Data.h!"); return 0.; }

       /**
       *  @brief return the value of the extra information at index i,j
       *  @param i index
       *  @param j index
       *  @return the value of the extra_info at position i,j
       */
      virtual double extra_info (const int i, const int j) const
      { (void)i; (void)j; ErrorCBL("Error in extra_info of Data.h!"); return 0.; }

      /**
       *  @brief return the m_exta_info vector
       *  @return vector containing vectors with extra information
       */
      virtual vector<vector<double>> extra_info () const
      { ErrorCBL("Error in extra_info of Data.h!"); vector<vector<double>> x; return x; }
      
      ///@}


      /**
       *  @name Member functions to set the private/protected members
       */
      ///@{

      /**
       *  @brief set the data type 
       *  @param dataType the data type
       *  @return none
       */
      void set_dataType (const DataType dataType) { m_dataType = dataType; } 
      
      /**
       *  @brief set interval variables
       *  @param min maximun value to be used 
       *  @param max maximun value to be used 
       *  @param axis the asis used
       *  @return none
       */
      virtual void set_limits (const double min, const double max, const bool axis) 
      { (void)min; (void)max; (void)axis; ErrorCBL("Error in set_limits of Data.h!"); }

      /**
       *  @brief set interval variables for x range
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @param ymin maximun value of y to be used 
       *  @param ymax maximun value of y to be used 
       *  @return none
       */
      virtual void set_limits (const double xmin, const double xmax, const double ymin, const double ymax) 
      { (void)xmin; (void)xmax; (void)ymin; (void)ymax; ErrorCBL("Error in set_limits of Data.h!"); }     

      /**
       *  @brief set interval variables for x range
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @return none
       */
      virtual void set_limits (const double xmin, const double xmax) 
      { (void)xmin; (void)xmax; ErrorCBL("Error in set_limits of Data.h!"); }

      /**
       *  @brief set interval variable m_x
       *  @param x vector containing x points
       *  @return none
       */
      virtual void set_xx (const vector<double> x)
      { (void)x; ErrorCBL("Error in set_xx of Data.h!"); }

      /**
       *  @brief set interval variable m_y
       *  @param y vector containing y points
       *  @return none
       */
      virtual void set_yy (const vector<double> y)
      { (void)y; ErrorCBL("Error in set_yy of Data.h!"); }

      /**
       *  @brief set interval variable m_fx
       *  @param fx vector containing f(x) values 
       *  @return none
       */
      virtual void set_fx (const vector<double> fx)
      { (void)fx; ErrorCBL("Error in set_fx of Data.h!"); }

      /**
       *  @brief set interval variable m_error_fx
       *  @param error_fx vector containing error on f(x)
       *  @return none
       */
      virtual void set_error_fx (const vector<double> error_fx)
      { (void)error_fx; ErrorCBL("Error in set_error_fx of Data.h!"); }

      /**
       *  @brief set interval variable m_covariance, reading from an input file
       *  @param filename file containing the covariance matrix in the format:
       * column 0 &rarr x<SUB>i</SUB>, column 1 &rarr x<SUB>j</SUB>, column 2 &rarr cov(x<SUB>i</SUB>,x<SUB>j</SUB>)
       *  @return none
       */
      virtual void set_covariance (const string filename)
      { (void)filename; ErrorCBL("Error in set_covariance of Data.h!"); }

      /**
       *  @brief set interval variable m_covariance
       *  @param covariance vector containing f(x) covariance matrix 
       *  @return none
       */
      virtual void set_covariance (const vector<vector<double>> covariance) 
      { (void)covariance; ErrorCBL("Error in set_covariance of Data.h!"); }

      /**
       *  @brief set interval variable m_fxy
       *  @param fxy vector containing f(x,y) 
       *  @return none
       */
      virtual void set_fxy (const vector<vector<double>> fxy) 
      { (void)fxy; ErrorCBL("Error in set_fxy of Data.h!"); }

      /**
       *  @brief set interval variable m_error_fxy
       *  @param error_fxy vector containing errors on f(x,y)
       *  @return none
       */ 
      virtual void set_error_fxy (const vector<vector<double>> error_fxy)
      { (void)error_fxy; ErrorCBL("Error in set_error_fxy of Data.h!"); }

      /**
       *  @brief set interval variables for x range in the i-th dataset
       *  @param i index to the i-th dataset
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @return none
       */
      virtual void set_limits (const int i, const double xmin, const double xmax)
      { (void)i; (void)xmin; (void)xmax; ErrorCBL("Error in set_limits of Data.h"); }

      /**
       *  @brief set interval variable m_x in the i-th dataset
       *  @param i index to the i-th dataset
       *  @param x vector containing x points
       *  @return none
       */
      virtual void set_xx (const int i, const vector<double> x) 
      { (void)i; (void)x; ErrorCBL("Error in set_xx of Data.h"); }

      /**
       *  @brief set interval variable m_fx in the i-th dataset
       *  @param i index to the i-th dataset
       *  @param fx vector containing f(x) values 
       *  @return none
       */
      virtual void set_fx (const int i, const vector<double> fx) 
      { (void)i; (void)fx; ErrorCBL("Error in set_fx of Data.h"); }

      /**
       *  @brief set interval variable m_error_fx in the i-th dataset
       *  @param i index to the i-th dataset
       *  @param error_fx vector containing error on f(x)
       *  @return none
       */
      virtual void set_error_fx (const int i, const vector<double> error_fx)
      { (void)i; (void)error_fx; ErrorCBL("Error in set_error_fx of Data.h"); }

      /**
       *  @brief set interval variable m_covariance, in the i-th
       *  dataset reading from an input file; also compute inverted
       *  covariance matrix
       *  @param i index to the i-th dataset
       *  @param filename file containing the covariance matrix in the
       *  format: column 0 &rarr x<SUB>i</SUB>, column 1 &rarr
       *  x<SUB>j</SUB>, column 2 &rarr
       *  cov(x<SUB>i</SUB>,x<SUB>j</SUB>)
       *  @return none
       */
      virtual void set_covariance (const int i, const string filename)
      { (void)i; (void)filename; ErrorCBL("Error in set_covariance of Data.h"); }

      /**
       *  @brief set interval variable m_covariance, in the i-th
       *  dataset reading from an input file; also compute inverted
       *  covariance matrix
       *  @param i index to the i-th dataset
       *  @param covariance vector containing f(x) covariance matrix 
       *  @return none
       */
      virtual void set_covariance (const int i, const vector<vector<double>> covariance) 
      { (void)i; (void)covariance; ErrorCBL("Error in set_covariance of Data.h"); }

      /**
       *  @brief set cross-correlation between i-th and j-th datasets
       *  @param i index to the i-th dataset
       *  @param j index to the j-th dataset
       *  @param covariance vector containing f(x) cross-covariance matrix 
       *  @return none
       */
      virtual void set_covariance (const int i, const int j, const vector<vector<double>> covariance)
      { (void)i; (void)j; (void)covariance; ErrorCBL("Error in set_covariance of Data.h"); }

      /**
       *  @brief set interval variable m_covariance_matrix, from covariance matrix
       *  of datasets, the result is a block covariance matrix
       *  @return none
       */
      virtual void set_covariance ()
      { ErrorCBL("Error in set_covariance of Data.h"); }

      /**
       *  @brief set interval variable m_error_fx
       *  @param extra_info vector containing vectors with extra
       *  information
       *  @return none
       */
      virtual void set_extra_info (const vector<vector<double>> extra_info)
      { (void)extra_info; ErrorCBL("Error in set_extra_info of Data.h"); }

      ///@}


      /**
       *  @name Member functions to compute data properties
       */
      ///@{

      /**
       *   @brief the effective number of data 
       *   @return effective number of data
       */
      virtual int ndata_eff () const
      { ErrorCBL("Error in ndata of Data.h!"); return 0; }

      /**
       *   @brief the total number of data
       *   @return total number of data
       */
      virtual int ndata () const
      { ErrorCBL("Error in ndata of Data.h!"); return 0; }

      /**
       * @brief function that returns effective number of data between defined limits
       * @param i index to the i-th dataset
       * @return effective number of data between defined limits
       */
      virtual int ndata_eff (const int i) const
      { (void)i; ErrorCBL("Error in ndata of Data.h!"); return 0; }

      /**
       * @brief function that returns total number of data
       * @param i index to the i-th dataset
       * @return total number of data
       */
      virtual int ndata (const int i) const 
      { (void)i; ErrorCBL("Error in ndata of Data.h!"); return 0; }

      /**
       * @brief function that returns total number of datasets
       * @return total number of dataset
       */
      virtual int ndataset () const
      { ErrorCBL("Error in ndataset of Data.h!"); return 0; }

      /**
       *  @brief invert the covariance matrix
       *  @return none
       */
      virtual void invert_covariance ()
      { ErrorCBL("Error in invert_covariance of Data.h!"); }
      
      ///@}


      /**
       *  @name Member functions for Input/Output
       */
      ///@{

      /**
       *  @brief read the data
       *  @param input_file file containing input data in 4 columns:
       *  first column &rarr x points, second column y points, thrid column &rarr f(x,y), fourth column &rarr 
       *  f(x,y) error
       *  @param skip_nlines the header lines to be skipped
       *  @return none
       */
      virtual void read (const string input_file, const int skip_nlines=0)
      { (void)input_file; (void)skip_nlines; ErrorCBL("Error in read of Data.h!"); }

      /**
       *  @brief write the data
       *  @param dir output directory
       *  @param file output file
       *  @param header text with the variable names to be written at
       *  the first line of the output file
       *  @param rank cpu index (for MPI usage)
       *  @return none
       */
      virtual void write (const string dir, const string file, const string header, const int rank=0) const 
      { (void)dir; (void)file; (void)header; (void)rank; ErrorCBL("Error in write of Data.h!"); }

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
      virtual void write (const string dir, const string file, const string header, const bool full, const int rank=0) const
      { (void)dir; (void)file; (void)header; (void)full; (void)rank; ErrorCBL("Error in write of Data.h!"); }
      
      /**
       *  @brief write the interval variable m_covariance on a file,
       *  @param dir the output directory
       *  @param file the output file
       *  @param xname name for the x variable
       *  @param fxname name for the f(x) variable
       *  @return none
       */
      virtual void write_covariance (const string dir, const string file, const string xname="x", const string fxname="fx") const
      { (void)dir; (void)file; (void)xname; (void)fxname; ErrorCBL("Error in write_covariance of Data.h!"); }
      
      ///@}
      
    };
  }
}

#endif
