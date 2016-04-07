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
      vector<vector<double> > m_covariance_fx;
      
      /// inverse covariance matrix of f(x)
      vector<vector<double> > m_inverse_covariance_fx;
      
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
      Data1D () {}

      /**
       *  @brief constructor of Data1D
       *  @param input_file file containing input data in 3 columns:
       *  first column &rarr x points, second column &rarr f(x), third column &rarr 
       *  f(x) error
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @return object of class Data1D
       */
      Data1D (const string input_file , const double xmin=-par::defaultDouble, const double xmax=par::defaultDouble); 

      /**
       *  @brief constructor of Data1D
       *  @param x vector containing x points 
       *  @param fx vector containing f(x) 
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @return object of class Data1D
       */
      Data1D (const vector<double> x, const vector<double> fx, const double xmin=-par::defaultDouble, const double xmax=par::defaultDouble); 

      /**
       *  @brief Constructor of Data1D
       *  @param x vector containing x points 
       *  @param fx vector containing f(x) 
       *  @param error_fx vector containing error on f(x) 
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @return object of class Data1D
       */
      Data1D (const vector<double> x, const vector<double> fx, const vector<double> error_fx, const double xmin=-par::defaultDouble, const double xmax=par::defaultDouble); 

      /**
       *  @brief Constructor of Data1D
       *  @param x vector containing x points 
       *  @param fx vector containing f(x) 
       *  @param covariance_fx vector containing f(x) covariance matrix 
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @return object of class Data1D
       */
      Data1D (const vector<double> x, const vector<double> fx, const vector<vector<double> > covariance_fx, const double xmin=-par::defaultDouble, const double xmax=par::defaultDouble);

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Data1D () {}

      ///@}

      /**
       *  @brief return index of the first x used
       *  @return int containing the index of the first x used
       */
      int x_down () const override { return m_x_down; } 

      /**
       *  @brief return index of the last x used
       *  @return int containing the index of the last x used
       */
      int x_up () const override { return m_x_up; } 

      /**
       *  @brief return value of x at index i
       *  @param i index
       *  @return value of the m_x vector at position i
       */
      double xx (const int i) const override { return m_x[i]; }  

      /**
       *  @brief return f(x) at index i
       *  @param i index
       *  @return value of the m_fx vector at position i
       */
      double fx (const int i) const override { return m_fx[i]; } 

      /**
       *  @brief return value of f(x) error at index i
       *  @param i index
       *  @return value of the m_error_fx vector at position i
       */
      double error_fx (const int i) const override { return m_error_fx[i]; } 

      /**
       *  @brief return value of f(x) covariance at index i,j
       *  @param i index
       *  @param j index
       *  @return value of the m_covariance_fx vector at position i,j
       */
      double covariance_fx (const int i, const int j) const override { return m_covariance_fx[i][j]; }

      /**
       *  @brief return value of f(x) inverse_covariance at index i,j
       *  @param i index
       *  @param j index
       *  @return value of the m_inverse_covariance_fx vector at position i,j
       */
      double inverse_covariance_fx (const int i, const int j) const override  { return m_inverse_covariance_fx[i][j]; }

      /**
       *  @brief return the x vector
       *  @return vector containing the x values
       */
      vector<double> xx () const override;   

      /**
       *  @brief return the m_fx vector
       *  @return vector containing the fx values
       */
      vector<double> fx () const override;  

      /**
       *  @brief return the m_error_fx vector
       *  @return vector containing the values of fx error
       */
      vector<double> error_fx () const override;  

      /**
       *  @brief return the m_covariance vector
       *  @return vector<containing the covariance matrix
       */
      vector<vector<double> > covariance_fx () const override;

      /**
       *  @brief return the m_inverse_covariance vector
       *  @return vector containing the inverse convariance matrix
       */
      vector<vector<double> > inverse_covariance_fx () const override;

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
      void set_fx (const vector<double> fx) override { m_fx=fx; }

      /**
       *  @brief set interval variable m_error_fx
       *  @param error_fx vector containing error on f(x)
       *  @return none
       */
      void set_error_fx (const vector<double> error_fx) override {m_error_fx=error_fx;}

      /**
       *  @brief set interval variable m_covariance_fx, reading from an input file;
       *  also compute inverted covariance matrix
       *  @param filename file containing the covariance matrix in the format:
       *  column 0 &rarr x<SUB>i</SUB>, column 1 &rarr x<SUB>j</SUB>, column 2 &rarr cov(x<SUB>i</SUB>,x<SUB>j</SUB>)
       *  @return none
       */
      void set_covariance_fx (const string filename) override;

      /**
       *  @brief set interval variable m_covariance_fx,
       *  also compute inverted covariance matrix
       *  @param covariance_fx vector containing f(x) covariance matrix 
       *  @return none
       */
      void set_covariance_fx (const vector<vector<double> > covariance_fx) override; 

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
       *  @brief read data from file
       *  @param input_file file containing input data in 3 columns:
       *  first column &rarr x points, second column &rarr f(x), third
       *  column &rarr f(x) error
       *  @return none
       */
      virtual void read (const string input_file=par::defaultString) override;

      /**
       *  @brief write the measured two-point correlation
       *  @param dir output directory
       *  @param file output file
       *  @param xname name for the x variable
       *  @param fxname name for the f(x)
       *  @param rank cpu index (for MPI usage)
       *  @return none
       */
      virtual void write (const string dir=par::defaultString, const string file=par::defaultString, const string xname="x", const string fxname="f(x)", const int rank=0) const override;

  };
}
#endif
