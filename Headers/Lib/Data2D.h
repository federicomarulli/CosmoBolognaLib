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
 *  @file Headers/Lib/Data2D.h
 *
 *  @brief The class Data2D
 *
 *  This file defines the interface of the class Data2D
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __DATA2D__
#define __DATA2D__

#include "Data.h"

namespace cosmobl {

  /**
   *  @class Data2D Data2D.h
   *  "Headers/Lib/Data2D.h"
   *
   *  @brief The class Data2D
   *
   *  This is the base class used to manage 2D data
   */
  class Data2D : public Data
  {
    protected:

      /**
       *  @name Data input
       */
      ///@{

      /// ordered x axis points
      vector<double> m_x;

      /// ordered y axis points
      vector<double> m_y;

      /// f(x,y) values
      vector<vector<double> > m_fxy;

      /// errors of f(x,y)
      vector<vector<double> > m_error_fxy;

      ///@}

      /// index of the first x data to be used
      int m_x_down;

      /// index of the last x data to be used
      int m_x_up;

      /// index of the first y data to be used
      int m_y_down;

      /// index of the last y data to be used
      int m_y_up;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Data
       */
      Data2D(){}
   
      /**
       *  @brief static factory used to construct Data2D
       *  @param x vector containing x points 
       *  @param y vector containing y points 
       *  @param fxy vector containing f(x,y) values
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @param ymin maximun value of y to be used 
       *  @param ymax maximun value of y to be used 
       *  @return shared pointer to an object of class Data
       */
      Data2D (const vector<double> x, const vector<double> y, const vector< vector<double> > fxy, const double xmin=-1.e30, const double xmax=1.e30, const double ymin=-1.e30, const double ymax=1.e30); 

      /**
       *  @brief static factory used to construct Data2D
       *  @param x vector containing x points 
       *  @param y vector containing y points 
       *  @param fxy vector containing f(x,y) values
       *  @param error_fxy vector containing error on f(x,y)
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @param ymin maximun value of y to be used 
       *  @param ymax maximun value of y to be used 
       *  @return shared pointer to an object of class Data
       */
      Data2D (const vector<double> x, const vector<double> y, const vector< vector<double> > fxy, const vector< vector<double> > error_fxy, const double xmin=-1.e30, const double xmax=1.e30, const double ymin=-1.e30, const double ymax=1.e30); 

      /**
       *  @brief default destructor
       *  @return none
       */
      ~Data2D () {}

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
       *  @brief return index of the first y used
       *  @return int containing the index of the first y used
       */
      int y_down () const override { return m_y_down; }

      /**
       *  @brief return index of the last y used
       *  @return int containing the index of the last y used
       */
      int y_up () const override { return m_y_up; }

      /**
       *  @brief return value of x at index i
       *  @param i index
       *  @return value of the m_x vector at position i
       */
      double xx (const int i) const override { return m_x[i]; }

      /**
       *  @brief return value of y at index i
       *  @param i index
       *  @return value of the m_y vector at position i
       */
      double yy (const int i) const override { return m_y[i]; }

      /**
       *  @brief return value of f(x,y) at index i,j
       *  @param i index
       *  @param j index
       *  @return value of the m_fxy vector at position i,j
       */
      double fxy (const int i, const int j) const override { return m_fxy[i][j]; }

      /**
       *  @brief return error on f(x,y) at index i,j
       *  @param i index
       *  @param j index
       *  @return value of the m_error_fxy vector at position i,j
       */
      double error_fxy (const int i, const int j) const override { return m_error_fxy[i][j]; }

      /**
       *  @brief return the x vector
       *  @return vector containing the x values
       */
      vector<double> xx () const override { return m_x; }

      /**
       *  @brief return the y vector
       *  @return vector containing the x values
       */
      vector<double> yy () const override { return m_y; }

      /**
       *  @brief return the m_fx vector
       *  @return vector containing the fx values
       */
      vector<vector<double> > fxy () const override { return m_fxy; }

      /**
       *  @brief return the m_error_fx vector
       *  @return vector containing the values of fx error
       */
      vector<vector<double> > error_fxy () const override { return m_error_fxy; }

      /**
       *  @brief set interval variables for x range
       *  @param min maximun value of x to be used 
       *  @param max maximun value of x to be used 
       *  @param axis 0 &rarr; acts on x axis ; 1 &rarr;
       *  acts on y axis
       *  @return none
       */
      void set_limits (const double min, const double max, bool axis) override;

      /**
       *  @brief set interval variables for x range
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @param ymin maximun value of y to be used 
       *  @param ymax maximun value of y to be used 
       *  @return none
       */
      void set_limits (const double xmin, const double xmax, const double ymin, const double ymax) override;

      /**
       *  @brief set interval variable m_x
       *  @param x vector containing x points
       *  @return none
       */
      void set_xx (const vector<double> x) override { m_x = x; }

      /**
       *  @brief set interval variable m_y
       *  @param y vector containing y points
       *  @return none
       */
      void set_yy (const vector<double> y) override { m_y = y; }

      /**
       *  @brief set interval variable m_fxy
       *  @param fxy vector containing f(x,y) 
       *  @return none
       */
      void set_fxy (const vector<vector<double> > fxy) override { m_fxy = fxy; }

      /**
       *  @brief set interval variable m_error_fxy
       *  @param error_fxy vector containing errors on f(x,y)
       *  @return none
       */ 
      void set_error_fxy (const vector<vector<double> > error_fxy) override { m_error_fxy = error_fxy; }

      /**
       * @brief function that returns effective number of data between
       * defined limits
       * @return effective number of data between defined limits
       */
      int ndata_eff () const override { return (m_x_up-m_x_down)*(m_y_up-m_y_down); }

      /**
       * @brief function that returns total number of data
       * @return total number of data
       */
      int ndata () const override { return m_x.size()*m_y.size(); }

      /**
       *  @brief read data from file
       *  @param input_file file containing input data in 3 columns:
       *  first column &rarr x points, second column &rarr f(x), third column &rarr 
       *  f(x) error
       *  @return none
       */
      virtual void read (const string input_file=par::defaultString) override;

      /**
       *  @brief write the measured two-point correlation
       *  @param dir output directory
       *  @param file output file
       *  @param xname name for the x variable
       *  @param yname name for the y variable
       *  @param fxyname name for the f(x,y)
       *  @param full 0 &rarr; simply store the data; 1 &rarr; duplicate
       *  the data in the other three quadrands (usefull when storing
       *  the 2D clustering)
       *  @param rank cpu index (for MPI usage)
       *  @return none
       */
      virtual void write (const string dir=par::defaultString, const string file=par::defaultString, const string xname="x", const string yname="y",const string fxyname="f(x,y)", const bool full=0, const int rank=0) const override;
      
  };

}
#endif
