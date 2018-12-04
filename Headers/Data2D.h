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
 *  @file Headers/Data2D.h
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

namespace cbl {

  namespace data {
  
    /**
     *  @class Data2D Data2D.h
     *  "Headers/Data2D.h"
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
      std::vector<double> m_x;

      /// ordered y axis points
      std::vector<double> m_y;

      /// number of points along x
      int m_xsize;

      /// number of points along y
      int m_ysize;

      ///@}

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return an object of class Data2D
       */
      Data2D ()
	: Data(DataType::_2D_) {}

       /**
       *  @brief constructor used 
       *  @param input_file input file 
       *  @param skip_nlines skip header lines
       *  @return an object of class Data2D
       */
      Data2D (const std::string input_file, const int skip_nlines=0)    
	: Data(DataType::_2D_)
	{ read(input_file, skip_nlines); }
      
      /**
       *  @brief constructor used 
       *  @param x vector containing x points 
       *  @param y vector containing y points 
       *  @param data vector containing f(x,y) values
       *  @return an object of class Data2D
       */
      Data2D (const std::vector<double> x, const std::vector<double> y, const std::vector<std::vector<double>> data);

      /**
       *  @brief constructor used 
       *  @param x vector containing x points 
       *  @param y vector containing y points 
       *  @param data vector containing data points
       *  @param error vector containing error in data points
       *  @return an object of class Data2D
       */
      Data2D (const std::vector<double> x, const std::vector<double> y, const std::vector<std::vector<double>> data, const std::vector<std::vector<double>> error); 

      /**
       *  @brief constructor used
       *  @param x vector containing x points 
       *  @param y vector containing y points 
       *  @param data vector containing data values
       *  @param error vector containing error on data points
       *  @return an object of class Data2D
       */
      Data2D (const std::vector<double> x, const std::vector<double> y, const std::vector<double> data, const std::vector<double> error); 

      /**
       *  @brief constructor used
       *  @param x vector containing x points 
       *  @param y vector containing y points 
       *  @param data vector containing data values
       *  @param covariance vector containing error on data points
       *  @return an object of class Data2D
       */
      Data2D (const std::vector<double> x, const std::vector<double> y, const std::vector<double> data, const std::vector<std::vector<double>> covariance); 

      /**
       *  @brief default destructor
       *  @return none
       */
      ~Data2D () = default;

      /**
       *  @brief static factory used to construct objects of class
       *  Data2D
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
      int xsize () const override
      { return m_xsize; }

      /**
       *  @brief get the number of points along y
       *  @return the number of points along y
       */
      int ysize () const override
      { return m_ysize; }

      /**
       *  @brief get the value of x at index i
       *  @param i index
       *  @return the value of the m_x vector at position i
       */
      double xx (const int i) const override
      { return m_x[i]; }

      /**
       *  @brief get the value of y at index i
       *  @param i index
       *  @return the value of the m_y vector at position i
       */
      double yy (const int i) const override
      { return m_y[i]; }

      /**
       *  @brief get the independet variable, to be used 
       *  in model computation
       *  @param i index of the extra_info containing
       *  the first independent variable
       *  @param j index of the extra_info containing
       *  the second independent variable
       *  @return the independent variable
       */
      std::vector<std::vector<double>> IndipendentVariable (const int i=-1, const int j=-1) const override
	{ (void)i; (void)j;  return {m_x, m_y}; } 

      /**
       *  @brief get data at index i,j
       *  @param i index
       *  @param j index
       *  @return the value of the m_data vector at position i,j
       */
      double data (const int i, const int j) const override
      { return m_data[j+i*m_ysize]; } 

      /**
       *  @brief get data
       *  @param [out] data vector containing the dataset
       *  @return none
       */
      void get_data (std::vector<std::vector<double>> &data) const override;

      /**
       *  @brief get error at index i,j
       *  @param i index
       *  @param j index
       *  @return the value of the m_error vector at position i,j
       */
      double error (const int i, const int j) const override
      { return m_error[j+i*m_ysize]; } 

      /**
       *  @brief get error
       *  @param [out] error vector containing the error
       *  @return none
       */
      void get_error (std::vector<std::vector<double>> &error) const override;

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
      void set_xx (const std::vector<double> x) override { m_x = x; m_xsize=m_x.size();}

      /**
       *  @brief set interval variable m_y
       *  @param y vector containing y points
       *  @return none
       */
      void set_yy (const std::vector<double> y) override { m_y = y; m_ysize=m_y.size(); }

      ///@}

      
      /**
       *  @name Input/Output member functions
       */
      ///@{
      
      /**
       *  @brief read the data
       *  @param input_file input data file
       *  @param skip_nlines the header lines to be skipped
       *  @return none
       */
      virtual void read (const std::string input_file, const int skip_nlines=0) override;

      /**
       *  @brief write the data
       *  @param dir output directory
       *  @param file output file
       *  @param header text with the variable names to be written at
       *  the first line of the output file
       *  @param full false &rarr; simply store the data; true &rarr;
       *  duplicate the data in the other three quadrands (usefull
       *  e.g. when storing the 2D correlation function)
       *  @param precision the floating point precision
       *  @param rank cpu index (for MPI usage)
       *  @return none
       */
      void write (const std::string dir, const std::string file, const std::string header, const bool full, const int precision=4, const int rank=0) const;
          
      /**
       *  @brief write the covariance
       *  @param dir the output directory
       *  @param file the output file
       *  @param precision the floating point precision
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
       * @param ymin minumum value for the independent variable y
       * @param ymax maximum value for the independent variable y
       * @return pointer to an object of type Data2D
       */
      std::shared_ptr<Data> cut(const double xmin, const double xmax, const double ymin, const double ymax) const;

      ///@}

    };

  }
}

#endif
