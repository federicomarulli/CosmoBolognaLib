/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file Headers/Lib/FuncGrid.h
 *
 *  @brief Class used to handle functions stored on a grid
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __FUNCGRID__
#define __FUNCGRID__ 


// =====================================================================================


namespace cosmobl {

  namespace glob {
    
    /**
     *  @class FuncGrid FuncGrid.h "Headers/Lib/FuncGrid.h"
     *
     *  @brief The class FuncGrid
     *
     *  This class is used to handle functions stored on a
     *  grid. Specifically, it contains member functions to
     *  interpolate, find minima, compute derivatives and integrals.
     */
    class FuncGrid
    {
      
    private:

      /// x values
      vector<double> m_x;

      /// y values
      vector<double> m_y;
      
      /// size of the x,y vectors, i.e. the grid size
      size_t m_size;
      
      /// method used to interpolate
      string m_interpType;

      /// GSL object used to interpolate 
      gsl_spline *m_spline;

      /// GSL object used to set the interpolation type 
      const gsl_interp_type *m_type;

      /// GSL accelerator object
      gsl_interp_accel *m_acc;

      /// minimum x value
      double m_xmin;

      /// maximum x value
      double m_xmax;

      
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{
      
      /**
       *  @brief default constructor
       *  @return object of class FuncGrid
       */
      FuncGrid () = default;

      /**
       *  @brief constructor
       *
       *  @param x vector containing the x values
       *  @param y vector containing the y values
       *  @param interpType interpolation method
       *
       *  @return object of class FuncGrid
       */
      FuncGrid (const vector<double> x, const vector<double> y, const string interpType);

      /**
       *  @brief default destructor
       *  @return none
       */
      ~FuncGrid () = default;

      ///@}

      
      /**
       *  @name Functions to get the private members of the class
       */
      ///@{

      /**
       *  @brief get the private member FuncGrid::m_x
       *  @return the vector containing the x values
       */
      vector<double> x () const { return m_x; }
      
      /**
       *  @brief get the i-th element of the private member
       *  FuncGrid::m_x
       *  @param i the i-th index
       *  @return the i-th value of x
       */
      double x (const int i) const { return m_x[i]; }
      
      /**
       *  @brief get the private member FuncGrid::m_y
       *  @return the vector containing the y values
       */
      vector<double> y () const { return m_y; }
      
      /**
       *  @brief get the i-th element of the private member
       *  FuncGrid::m_y
       *  @param i the i-th index
       *  @return the i-th value of y
       */
      double y (const int i) const { return m_y[i]; }

      /**
       *  @brief get the private member FuncGrid::m_size
       *  @return the size of the x,y vectors, i.e. the grid size
       */
      double size () const { return m_size; }

      /**
       *  @brief get the private member FuncGrid::m_xmin
       *  @return minimum x value
       */
      double xmin () const { return m_xmin; }

      /**
       *  @brief get the private member FuncGrid::m_xmax
       *  @return maximum x value
       */                            
      double xmax () const { return m_xmax; }
			 
      ///@}

      
      /**
       *  @name Functions of generic usage
       */
      ///@{

      /**
       *  @brief free the GSL objects
       *  @return none
       */   
      void free ();

      /**
       *  @brief overloading of the () operator
       *  @param xx the value at which the function will be
       *  evaluated
       *  @return the function evaluated at xx
       */   
      double operator () (const double xx) const;

      /**
       *  @brief evaluate the function at the xx points
       *  @param xx the values at which the function will be
       *  evaluated
       *  @return the function evaluated at the xx points
       */   
      vector<double> eval_func (const vector<double> xx) const;

      /**
       *  @brief compute the first derivative at xx
       *  @param xx the value at which the derivative is computed
       *  @return the first derivative
       */   
      double D1v (const double xx) const;
      
      /**
       *  @brief compute the second derivative at xx
       *  @param xx the value at which the derivative is computed
       *  @return the second derivative
       */  
      double D2v (const double xx) const;
   
      /**
       *  @brief compute the definite integral with GSL qag method
       *  @param a the lower limit of the integral
       *  @param b the upper limit of the integral
       *  @param prec the relative error tolerance
       *  @param limit_size the maximum size of workspace
       *  @param rule the rule of integration
       *  @return the definite integral of the function
       */  
      double integrate_qag (const double a, const double b, const double prec=1.e-2, const int limit_size=6, const int rule=6);

      /**
       *  @brief compute the definite integral with GSL qaws method
       *  @param a the lower limit of the integral
       *  @param b the upper limit of the integral
       *  @param alpha &alpha;
       *  @param beta &beta;
       *  @param mu &mu;
       *  @param nu &nu;
       *  @param prec the relative error tolerance
       *  @param limit_size the maximum size of workspace
       *  @return the definite integral of the function
       */  
      double integrate_qaws (const double a, const double b, const double alpha=0, const double beta=0, const int mu=0, const int nu=0, const double prec=1.e-2, const int limit_size=6);
      
      /**
       *  @brief find roots with GSL brent method 
       *  @param x_low the lower limit 
       *  @param x_up the upper limit
       *  @param fx0 fx0
       *  @param prec the relative error tolerance
       *  @return the root
       */
      double root (const double x_low, const double x_up, const double fx0=0, const double prec=1.e-2);

      /**
       *  @brief find roots with GSL brent method 
       *  for the first derivative
       *  @param x_low the lower limit 
       *  @param x_up the upper limit
       *  @param fx0 fx0
       *  @param prec the relative error tolerance
       *  @return the root
       */
      double root_D1v (const double x_low, const double x_up, const double fx0=0, const double prec=1.e-2); 
 
      /**
       *  @brief find roots with GSL brent method 
       *  for the second derivative
       *  @param x_low the lower limit 
       *  @param x_up the upper limit
       *  @param fx0 fx0
       *  @param prec the relative error tolerance
       *  @return the root
       */
      double root_D2v (const double x_low, const double x_up, const double fx0=0, const double prec=1.e-2);
    };

  }
}
    
#endif
