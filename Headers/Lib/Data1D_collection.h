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
 *  @file Headers/Lib/Data1D_collection.h
 *
 *  @brief The class Data1D_collection
 *
 *  This file defines the interface of the class Data1D
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __DATA1DC__
#define __DATA1DC__

#include "Data.h"
#include "Data1D.h"

namespace cosmobl {

  /**
   *  @class Data1D_collection Data1D_collection.h
   *  "Headers/Lib/Data1D_collection.h"
   *
   *  @brief The class Data1D_collection
   *
   *  This is the base class used to manage collection
   *  of 1D data
   */
  class Data1D_collection : public Data
  {
    protected:

      /**
       *  @name Data input
       */
      ///@{

      /// ordered x axis points
      vector<Data1D> m_data;

      /// covariance_matrix
     vector<vector<double> > m_covariance_matrix;
     
      /// inverse covariance_matrix
     vector<vector<double> > m_inverse_covariance_matrix;

      /// members of the collection
      int m_n_data;

      ///@}
      
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Data1D
       */
      Data1D_collection () {}

      /**
       *  @brief constructor of class Data1D_collection
       *  @param n_data the number of dataset
       *  @return object of class Data1D_collection
       */
      Data1D_collection (const int n_data);

      /**
       *  @brief constructor of class Data1D_collection
       *  @param data vector containing the datasets
       *  @param x_min vector containing miminum values for the datasets 
       *  @param x_max vector containing maximum values for the datasets 
       *  @return object of class Data1D_collection
       */
      Data1D_collection(const vector<Data1D> data, const vector<double> x_min={}, const vector<double> x_max={});

      /**
       *  @brief constructor of class Data1D_collection
       *  @param data vector containing the datasets
       *  @param covariance_matrix vector containing covariance_matrix values of the datasets
       *  @param x_min vector containing miminum values for the datasets 
       *  @param x_max vector containing maximum values for the datasets 
       *  @return object of class Data1D_collection
       */
      Data1D_collection (const vector<Data1D> data, const vector<vector<double > > covariance_matrix, const vector<double> x_min={}, const vector<double> x_max={});

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Data1D_collection () {}

      ///@}
      

      /**
       *  @brief return index of the first x used in the i-th dataset
       *  @param i the i-th dataset
       *  @return int containing the index of the first x used
       */
      int x_down (const int i) const 
      { return m_data[i].x_down(); } 

      /**
       *  @brief return index of the last x usedin the i-th dataset
       *  @param i the i-th dataset
       *  @return int containing the index of the last x used
       */
      int x_up (const int i) const 
      { return m_data[i].x_up(); } 

      /**
       *  @brief return value of x at index j in the i-th dataset
       *  @param i the i-th dataset
       *  @param j index
       *  @return value of the x vector at position j in the i-th dataset
       */
      double xx (const int i, const int j) const 
      { return m_data[i].xx(j); }  

      /**
       *  @brief return the x vector for all datasets, concatenated
       *  @return vector containing the x values for all datasets, concatenated
       */
      vector<double> xx () const override;

      /**
       *  @brief return f(x) at index j in the i-th dataset
       *  @param i the i-th dataset
       *  @param j index
       *  @return value of the fx vector at position j in the i-th dataset
       */
      double fx (const int i, const int j) const
      { return m_data[i].fx(j); } 

      /**
       *  @brief return the fx vector for all datasets, concatenated
       *  @return vector containing the fx values for all datasets, concatenated
       */
      vector<double> fx () const override;

      /**
       *  @brief return error on f(x) at index j in the i-th dataset
       *  @param i the i-th dataset
       *  @param j index
       *  @return value of the error_fx vector at position j in the i-th dataset
       */
      double error_fx (const int i, const int j) const
      { return m_data[i].error_fx(j); } 

      /**
       *  @brief return the error_fx vector for all datasets, concatenated
       *  @return vector containing the error_fx values for all datasets, concatenated
       */
      vector<double> error_fx () const override;

      /**
       *  @brief return value of f(x) covariance at index i,j
       *  @param i the i-th x element of the covariance matrix
       *  @param j the j-th x element of the covariance matrix
       *  @return value of the covariance matrix for datasets i,j at position ax1,ax2
       */
      double covariance (const int i, const int j) const
      {return m_covariance_matrix[i][j];}

      /**
       *  @brief invert the covariance matrix
       *  @return none
       */  
      void invert_covariance() override;

      /**
       *  @brief return value of f(x) inverted covariance at index i,j
       *  @param i the i-th x element of the covariance matrix 
       *  @param j the j-th x element of the covariance matrix 
       *  @return value of the inverted covariance matrix at position i,j
       */
      double inverse_covariance (const int i, const int j) const
      {return m_inverse_covariance_matrix[i][j];}

      /**
       *  @brief return value of f(x) inverted covariance at index i,j
       *  @param d the d-th dataset
       *  @param i the i-th x element of the covariance matrix 
       *  @param j the j-th x element of the covariance matrix 
       *  @return value of the inverted covariance matrix for at position i,j
       */
      double inverse_covariance (const int d, const int i, const int j) const
      {return m_data[d].inverse_covariance(i,j);}
      
      /**
       *  @brief set interval variables for x range in the i-th dataset
       *  @param i index to the i-th dataset
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @return none
       */
      void set_limits (const int i, const double xmin, const double xmax) 
      { m_data[i].set_limits(xmin, xmax);}

      /**
       *  @brief set interval variable m_x in the i-th dataset
       *  @param i index to the i-th dataset
       *  @param x vector containing x points
       *  @return none
       */
      void set_xx (const int i, const vector<double> x) 
      { m_data[i].set_xx(x); } 

      /**
       *  @brief set interval variable m_fx in the i-th dataset
       *  @param i index to the i-th dataset
       *  @param fx vector containing f(x) values 
       *  @return none
       */
      void set_fx (const int i, const vector<double> fx) 
      { m_data[i].set_fx(fx); }

      /**
       *  @brief set interval variable m_error_fx in the i-th dataset
       *  @param i index to the i-th dataset
       *  @param error_fx vector containing error on f(x)
       *  @return none
       */
      void set_error_fx (const int i, const vector<double> error_fx)
      { m_data[i].set_error_fx(error_fx); }

      /**
       *  @brief set interval variable m_covariance reading from an input file; 
       *  @param filename file containing the covariance matrix in the
       *  format: column 0 &rarr x<SUB>i</SUB>, column 1 &rarr
       *  x<SUB>j</SUB>, column 2 &rarr
       *  cov(x<SUB>i</SUB>,x<SUB>j</SUB>)
       *  @return none
       */
      void set_covariance (const string filename) override;

      /**
       *  @brief set the covariance matrix for the i-th dataset
       *  @param i index for the i-th dataset 
       *  @param covariance vector containing f(x) covariance matrix 
       *  @return none
       */
      void set_covariance (const int i, const vector<vector<double> > covariance)
      { m_data[i].set_covariance(covariance); }

      /**
       *  @brief set interval variable m_covariance_matrix, from covariance matrix
       *  of datasets, the result is a block covariance matrix
       *  @return none
       */
      void set_covariance () override;

      /**
       *  @brief set interval variable m_covariance
       *  @param covariance vector containing f(x) covariance matrix 
       *  @return none
       */
      void set_covariance (const vector<vector<double> > covariance)
      { m_covariance_matrix = covariance; }

      /**
       * @brief function that returns effective number of data between defined limits
       * @return effective number of data between defined limits
       */
      int ndata_eff () const override; 

      /**
       * @brief function that returns total number of data
       * @return total number of data
       */
      int ndata () const override;  

      /**
       * @brief function that returns effective number of data between defined limits
       * @param i index to the i-th dataset
       * @return effective number of data between defined limits
       */
      int ndata_eff (const int i) const
      { return m_data[i].ndata_eff();}

      /**
       * @brief function that returns total number of data
       * @param i index to the i-th dataset
       * @return total number of data
       */
      int ndata (const int i) const 
      { return m_data[i].ndata();} 

      /**
       * @brief function that returns total number of datasets
       * @return total number of dataset
       */
      int ndataset () const 
      { return m_data.size();}  
  };
}
#endif
