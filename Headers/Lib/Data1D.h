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
 *  @file Headers/Lib/Data1D.h
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

namespace cosmobl {

  namespace data {
    
    /**
     *  @class Data1D Data1D.h
     *  "Headers/Lib/Data1D.h"
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
      vector<double> m_x;

      /// f(x) values
      vector<double> m_fx;
      
      /// errors at f(x)
      vector<double> m_error_fx;
      
      /// covariance matrix of f(x)
      vector<vector<double>> m_covariance;
      
      /// inverse covariance matrix of f(x)
      vector<vector<double>> m_inverse_covariance;
      
      ///@}
     
      /// index of the first x data to be used
      int m_x_down;
      
      /// index of the last x data to be used
      int m_x_up;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Data1D
       */
      Data1D () { set_dataType(DataType::_1D_data_); }

      /**
       *  @brief constructor, reading data from an input file
       *  @param input_file the input data file
       *  @param skipped_lines the header lines to be skipped
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @param dataType the data type
       *  @return object of class Data1D
       */
      Data1D (const string input_file, const int skipped_lines=0, const double xmin=par::defaultDouble, const double xmax=-par::defaultDouble, const DataType dataType=DataType::_1D_data_); 

      /**
       *  @brief constructor
       *  @param x vector containing x points 
       *  @param fx vector containing f(x) 
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used
       *  @param dataType the data type
       *  @return object of class Data1D
       */
      Data1D (const vector<double> x, const vector<double> fx, const double xmin=par::defaultDouble, const double xmax=-par::defaultDouble, const DataType dataType=DataType::_1D_data_); 

      /**
       *  @brief constructor
       *  @param x vector containing x points 
       *  @param fx vector containing f(x) 
       *  @param error_fx vector containing error on f(x) 
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @param dataType the data type
       *  @return object of class Data1D
       */
      Data1D (const vector<double> x, const vector<double> fx, const vector<double> error_fx, const double xmin=par::defaultDouble, const double xmax=-par::defaultDouble, const DataType dataType=DataType::_1D_data_); 

      /**
       *  @brief Constructor
       *  @param x vector containing x points 
       *  @param fx vector containing f(x) 
       *  @param covariance vector containing f(x) covariance matrix 
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @param dataType the data type
       *  @return object of class Data1D
       */
      Data1D (const vector<double> x, const vector<double> fx, const vector<vector<double> > covariance, const double xmin=par::defaultDouble, const double xmax=-par::defaultDouble, const DataType dataType=DataType::_1D_data_);

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Data1D () = default;

      ///@}


      /**
       *  @name Member functions used to get the private members
       */
      ///@{

      /**
       *  @brief get index of the first x used
       *  @return int containing the index of the first x used
       */
      int x_down () const override { return m_x_down; } 

      /**
       *  @brief get index of the last x used
       *  @return int containing the index of the last x used
       */
      int x_up () const override { return m_x_up; } 

      /**
       *  @brief get the value of x at the i-th bin
       *  @param i index
       *  @return value of x at the i-th bin
       */
      double xx (const int i) const override { return m_x[i]; }  

      /**
       *  @brief get the value of f(x) at the i-th bin
       *  @param i i-th bin
       *  @return value of f(x) at the i-th bin
       */
      double fx (const int i) const override { return m_fx[i]; } 

      /**
       *  @brief get the value of the f(x) error at the i-th bin
       *  @param i index
       *  @return value of the f(x) error at the i-th bin
       */
      double error_fx (const int i) const override { return m_error_fx[i]; } 

      /**
       *  @brief get the value of f(x) covariance at index i,j
       *  @param i index
       *  @param j index
       *  @return value of the m_covariance matrix at position i,j
       */
      double covariance (const int i, const int j) const override { return m_covariance[i][j]; }

      /**
       *  @brief get the value of f(x) inverse_covariance at index
       *  i,j
       *  @param i index
       *  @param j index
       *  @return value of the m_inverse_covariance matrix at position i,j
       */
      double inverse_covariance (const int i, const int j) const override  { return m_inverse_covariance[i][j]; }
      
      /**
       *  @brief get the x vector
       *  @return the vector containing the x values
       */
      vector<double> xx () const override;   

      /**
       *  @brief get the m_fx vector
       *  @return the vector containing the fx values
       */
      vector<double> fx () const override;  

      /**
       *  @brief get the m_error_fx vector
       *  @return the vector containing the values of fx errors
       */
      vector<double> error_fx () const override;  

      /**
       *  @brief get the m_covariance vector
       *  @return the vector containing the covariance matrix
       */
      vector<vector<double>> covariance () const override;

      /**
       *  @brief get the m_inverse_covariance vector
       *  @return the vector containing the inverse convariance matrix
       */
      vector<vector<double>> inverse_covariance () const override;

      ///@}

      
      /**
       *  @name Member functions used to set the private members
       */
      ///@{
      
      /**
       *  @brief set interval variables for x range
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @return none
       */
      void set_limits (const double xmin, const double xmax) override;

      /**
       *  @brief set interval variable m_x
       *  @param x vector containing x points
       *  @return none
       */
      void set_xx (const vector<double> x) override { m_x = x; } 

      /**
       *  @brief set interval variable m_fx
       *  @param fx vector containing f(x) values 
       *  @return none
       */
      void set_fx (const vector<double> fx) override { m_fx = fx; }

      /**
       *  @brief set interval variable m_error_fx
       *  @param error_fx vector containing error on f(x)
       *  @return none
       */
      void set_error_fx (const vector<double> error_fx) override { m_error_fx = error_fx; }

      /**
       *  @brief set the interval variable m_covariance, reading from an input file
       *
       *  @param filename file containing the covariance matrix in the
       *  format: column 0 &rarr x<SUB>i</SUB>, column 1 &rarr
       *  x<SUB>j</SUB>, column 2 &rarr
       *  cov(x<SUB>i</SUB>,x<SUB>j</SUB>)
       *
       *  @param skipped_lines comment lines to be skipped
       *
       *  @return none
       */
      void set_covariance (const string filename, const int skipped_lines=0) override;

      /**
       *  @brief set interval the variable m_covariance
       *  @param covariance vector containing f(x) covariance matrix
       *  @return none
       */
      void set_covariance (const vector<vector<double>> covariance) override; 

      ///@}
      

      /**
       *  @name Member functions to compute data properties
       */
      ///@{
      
      /**
       *  @brief function that returns effective number of data between
       *  defined limits
       *  @return effective number of data between defined limits
       */
      int ndata_eff () const override { return m_x_up-m_x_down; }

      /**
       *  @brief function that returns total number of data
       *  @return total number of data
       */
      int ndata () const override { return m_x.size(); }
      
      /**
       *  @brief invert the covariance matrix
       *  @return none
       */
      void invert_covariance () override;

      ///@}

      
      /**
       *  @name Member functions for Input/Output 
       */
      ///@{

      /**
       *  @brief read the data
       *  @param input_file input data file
       *  @param skipped_lines the header lines to be skipped
       *  @return none
       */
      virtual void read (const string input_file, const int skipped_lines=0) override;

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
      
      /**
       *  @brief write the covariance
       *  @param dir the output directory
       *  @param file the output file
       *  @param xname name for the x variable
       *  @param fxname name for the f(x) variable
       *  @return none
       */
      virtual void write_covariance (const string dir, const string file, const string xname="x", const string fxname="fx") const override;
      
      ///@}
      
    };
    
  }
}

#endif
