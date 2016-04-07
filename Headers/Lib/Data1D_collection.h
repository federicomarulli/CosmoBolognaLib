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
       *  @return object of class Data1D_collection
       */
      Data1D_collection (const int n_data);

      /**
       *  @brief constructor of class Data1D_collection
       *  @param x1 vector containing first data set x values
       *  @param fx1 vector containing first data set f(x) values
       *  @param x2 vector containing second data set x values
       *  @param fx2 vector containing second data set f(x) values
       *  @param x1_min miminum value of x1 to be used 
       *  @param x1_max maximun value of x1 to be used 
       *  @param x2_min minimum value of x2 to be used 
       *  @param x2_max maximun value of x2 to be used 
       *  @return object of class Data1D_collection
       */
      Data1D_collection (const vector<double> x1, const vector<double> fx1, const vector<double> x2, const vector<double> fx2, const double x1_min=-1.e30, const double x1_max=1.e30, const double x2_min=-1.e30, const double x2_max=1.e30 );

      /**
       *  @brief constructor of class Data1D_collection
       *  @param x1 vector containing first data set x values
       *  @param fx1 vector containing first data set f(x) values
       *  @param x2 vector containing second data set x values
       *  @param fx2 vector containing second data set f(x) values
       *  @param x3 vector containing third data set x values
       *  @param fx3 vector containing third data set f(x) values
       *  @param x1_min miminum value of x1 to be used 
       *  @param x1_max maximun value of x1 to be used 
       *  @param x2_min minimum value of x2 to be used 
       *  @param x2_max maximun value of x2 to be used 
       *  @param x3_min minimum value of x3 to be used 
       *  @param x3_max maximun value of x3 to be used 
       *  @return object of class Data1D_collection
       */
      Data1D_collection(const vector<double> x1, const vector<double> fx1, const vector<double> x2, const vector<double> fx2, const vector<double> x3, const vector<double> fx3, const double x1_min=-1.e30, const double x1_max=1.e30, const double x2_min=-1.e30, const double x2_max=1.e30 , const double x3_min=-1.e30, const double x3_max=1.e30);

      /**
       *  @brief constructor of class Data1D_collection
       *  @param x vector containing x values of the datasets
       *  @param fx vector containing fx values of the datasets
       *  @param x_min vector containing miminum values for the datasets 
       *  @param x_max vector containing maximum values for the datasets 
       *  @return object of class Data1D_collection
       */
      Data1D_collection(const vector< vector< double > > x, const vector< vector< double > > fx, const vector<double> x_min={}, const vector<double> x_max={});

      /**
       *  @brief constructor of class Data1D_collection
       *  @param x1 vector containing first data set x values
       *  @param fx1 vector containing first data set f(x) values
       *  @param error_fx1 vector containing first data set error on f(x)
       *  @param x2 vector containing second data set x values
       *  @param fx2 vector containing second data set f(x) values
       *  @param error_fx2 vector containing second data set error on f(x)
       *  @param x1_min miminum value of x1 to be used 
       *  @param x1_max maximun value of x1 to be used 
       *  @param x2_min minimum value of x2 to be used 
       *  @param x2_max maximun value of x2 to be used 
       *  @return object of class Data1D_collection
       */
      Data1D_collection(const vector<double> x1, const vector<double> fx1, const vector<double> error_fx1, const vector<double> x2, const vector<double> fx2, const vector<double> error_fx2, const double x1_min=-1.e30, const double x1_max=1.e30, const double x2_min=-1.e30, const double x2_max=1.e30 );

      /**
       *  @brief constructor of class Data1D_collection
       *  @param x1 vector containing first data set x values
       *  @param fx1 vector containing first data set f(x) values
       *  @param error_fx1 vector containing first data set error on f(x)
       *  @param x2 vector containing second data set x values
       *  @param fx2 vector containing second data set f(x) values
       *  @param error_fx2 vector containing second data set error on f(x)
       *  @param x3 vector containing third data set x values
       *  @param fx3 vector containing third data set f(x) values
       *  @param error_fx3 vector containing third data set error on f(x)
       *  @param x1_min miminum value of x1 to be used 
       *  @param x1_max maximun value of x1 to be used 
       *  @param x2_min minimum value of x2 to be used 
       *  @param x2_max maximun value of x2 to be used 
       *  @param x3_min minimum value of x3 to be used 
       *  @param x3_max maximun value of x3 to be used 
       *  @return object of class Data1D_collection
       */
      Data1D_collection(const vector<double> x1, const vector<double> fx1, const vector<double> error_fx1, const vector<double> x2, const vector<double> fx2, const vector<double> error_fx2, const vector<double> x3, const vector<double> fx3, const vector<double> error_fx3, const double x1_min=-1.e30, const double x1_max=1.e30, const double x2_min=-1.e30, const double x2_max=1.e30, const double x3_min=-1.e30, const double x3_max=1.e30);

      /**
       *  @brief constructor of class Data1D_collection
       *  @param x vector containing x values of the datasets
       *  @param fx vector containing f(x) values of the datasets
       *  @param error_fx vector containing error_fx values of the datasets
       *  @param x_min vector containing miminum values for the datasets 
       *  @param x_max vector containing maximum values for the datasets 
       *  @return object of class Data1D_collection
       */
      Data1D_collection(const vector<vector<double > > x, const vector<vector<double > > fx, const vector<vector<double > > error_fx, const vector<double> x_min={}, const vector<double> x_max={});

      /**
       *  @brief constructor of class Data1D_collection
       *  @return object of class Data1D_collection
       */
      Data1D_collection(const vector<double> x1, const vector<double> fx1, const vector<vector<double>> covariance_fx1, const vector<double> x2, const vector<double> fx2, const vector<vector<double>> covariance_fx2, const vector<vector<double> > covariance_12 = {}, const double x1_min=-1.e30, const double x1_max=1.e30, const double x2_min=-1.e30, const double x2_max=1.e30 );

      /**
       *  @brief constructor of class Data1D_collection
       *  @return object of class Data1D_collection
       */
      Data1D_collection(const vector<double> x1, const vector<double> fx1, const vector<vector<double>> covariance_fx1, const vector<double> x2, const vector<double> fx2, const vector<vector<double>> covariance_fx2, const vector<double> x3, const vector<double> fx3, const vector<vector<double>> covariance_fx3, const vector<vector<double> > covariance_12 = {}, const vector<vector<double> > covariance_13 = {}, const vector<vector<double> > covariance_23 = {}, double x1_min=-1.e30, const double x1_max=1.e30, const double x2_min=-1.e30, const double x2_max=1.e30, const double x3_min=-1.e30, const double x3_max=1.e30);

      /**
       *  @brief constructor of class Data1D_collection
       *  @param data vector containing the datasets
       *  @return object of class Data1D_collection
       */
      Data1D_collection (const vector<Data1D> data);

      /**
       *  @brief constructor of class Data1D_collection
       *  @param input_file vector containing input files for the datasets 
       *  in the format x,f(x),error
       *  @param x_min vector containing miminum values for the datasets 
       *  @param x_max vector containing maximum values for the datasets 
       *  @return object of class Data1D_collection
       */
      Data1D_collection (const vector<string> input_file, const vector<double> x_min={}, const vector<double> x_max={});

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
      int x_down (const int i) const override { return m_data[i].x_down(); } 

      /**
       *  @brief return index of the last x usedin the i-th dataset
       *  @param i the i-th dataset
       *  @return int containing the index of the last x used
       */
      int x_up (const int i) const override { return m_data[i].x_up(); } 

      /**
       *  @brief return value of x at index j in the i-th dataset
       *  @param i the i-th dataset
       *  @param j index
       *  @return value of the x vector at position j in the i-th dataset
       */
      double xx (const int i, const int j) const override { return m_data[i].xx(j); }  

      /**
       *  @brief return f(x) at index j in the i-th dataset
       *  @param i the i-th dataset
       *  @param j index
       *  @return value of the fx vector at position j in the i-th dataset
       */
      double fx (const int i, const int j) const override { return m_data[i].fx(j); } 

      /**
       *  @brief return error on f(x) at index j in the i-th dataset
       *  @param i the i-th dataset
       *  @param j index
       *  @return value of the error_fx vector at position j in the i-th dataset
       */
      double error_fx (const int i, const int j) const override { return m_data[i].error_fx(j); } 

      /**
       *  @brief return value of f(x) covariance at index i,j
       *  @param i index of the i-th dataset
       *  @param j index of the i-th dataset
       *  @param ax1 index to the ax1-th x element for the dataset i 
       *  @param ax2 index to the ax2-th x element for the dataset j
       *  @return value of the covariance matrix for datasets i,j at position ax1,ax2
       */
      double covariance_fx (const int i, const int j, const int ax1, const int ax2) const override;

      /**
       *  @brief return value of f(x) inverted covariance at index i,j
       *  @param i index of the i-th dataset
       *  @param j index of the i-th dataset
       *  @param ax1 index to the ax1-th x element for the dataset i 
       *  @param ax2 index to the ax2-th x element for the dataset j
       *  @return value of the inverted covariance matrix for datasets i,j at position ax1,ax2
       */
      double inverse_covariance_fx (const int i, const int j, const int ax1, const int ax2) const override;

      /**
       *  @brief set interval variables for x range in the i-th dataset
       *  @param i index to the i-th dataset
       *  @param xmin maximun value of x to be used 
       *  @param xmax maximun value of x to be used 
       *  @return none
       */
      void set_limits (const int i, const double xmin, const double xmax) override;

      /**
       *  @brief set interval variable m_x in the i-th dataset
       *  @param i index to the i-th dataset
       *  @param x vector containing x points
       *  @return none
       */
      void set_xx (const int i, const vector<double> x) override { m_data[i].set_xx(x); } 

      /**
       *  @brief set interval variable m_fx in the i-th dataset
       *  @param i index to the i-th dataset
       *  @param fx vector containing f(x) values 
       *  @return none
       */
      void set_fx (const int i,const vector<double> fx) override { m_data[i].set_fx(fx); }

      /**
       *  @brief set interval variable m_error_fx in the i-th dataset
       *  @param i index to the i-th dataset
       *  @param error_fx vector containing error on f(x)
       *  @return none
       */
      void set_error_fx (const int i, const vector<double> error_fx) override { m_data[i].set_error_fx(error_fx); }

      /**
       *  @brief set interval variable m_covariance_fx,  in the i-th dataset
       * reading from an input file; also compute inverted covariance matrix
       *  @param i index to the i-th dataset
       *  @param filename file containing the covariance matrix in the format:
       * column 0 &rarr x<SUB>i</SUB>, column 1 &rarr x<SUB>j</SUB>, column 2 &rarr cov(x<SUB>i</SUB>,x<SUB>j</SUB>)
       *  @return none
       */
      void set_covariance_fx (const int i, const string filename) override;

      /**
       *  @brief set interval variable m_covariance_fx,  in the i-th dataset
       * reading from an input file; also compute inverted covariance matrix
       *  @param i index to the i-th dataset
       *  @param covariance_fx vector containing f(x) covariance matrix 
       *  @return none
       */
      void set_covariance_fx (const int i, const vector<vector<double> > covariance_fx) override;

      /**
       *  @brief set cross-correlation between i-th and j-th datasets
       *  @param i index to the i-th dataset
       *  @param j index to the j-th dataset
       *  @param covariance_fx vector containing f(x) cross-covariance matrix 
       *  @return none
       */
      void set_covariance_fx (const int i, const int j, const vector<vector<double> > covariance_fx) override; 

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
  };
}
#endif
