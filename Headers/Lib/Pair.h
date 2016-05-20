/*******************************************************************
 *  Copyright (C) 2015 by Federico Marulli                         *
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
 *  @file Headers/Lib/Pair.h
 *
 *  @brief The class Pair
 *
 *  This file defines the interface of the base class Pair, used to
 *  handle pairs of objects of any kind
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unbo.it
 */

#ifndef __PAIR__
#define __PAIR__


#include "Catalogue.h"


// ===================================================================================================


namespace cosmobl {
  
  /**
   *  @brief The namespace of the pairs 
   *  
   *  The \e pairs namespace contains all the functions and classes to
   *  handle pairs of objects
   */
  namespace pairs {
    
    /**
     * @enum PairType
     * @brief the pair type
     */
    enum PairType { 
      
      /// 1D pair in angular coordinates in linear bins
      _angular_lin_,
    
      /// 1D pair in angular coordinates in logarithmic bins
      _angular_log_,

      /// 1D pair in comoving coordinates in linear bins
      _comoving_lin_,
    
      /// 1D pair in comoving coordinates in logarithmic bins
      _comoving_log_,
    
      /// 2D pair in comoving Cartesian coordinates (r<SUB>p</SUB>, &pi;) in linear-linear bins
      _comovingCartesian_linlin_,

      /// 2D pair in comoving Cartesian coordinates (r<SUB>p</SUB>, &pi;) in linear-logarithmic bins
      _comovingCartesian_linlog_,

      /// 2D pair in comoving Cartesian coordinates (r<SUB>p</SUB>, &pi;) in logarithmic-linear bins
      _comovingCartesian_loglin_,

      /// 2D pair in comoving Cartesian coordinates (r<SUB>p</SUB>, &pi;) in logarithmic-logarithmic bins
      _comovingCartesian_loglog_,

      /// 2D pair in comoving polar coordinates (r, &mu;) in linear-linear bins
      _comovingPolar_linlin_,

      /// 2D pair in comoving polar coordinates (r, &mu;) in linear-logarithmic bins
      _comovingPolar_linlog_,

      /// 2D pair in comoving polar coordinates (r, &mu;) in logarithmic-linear bins
      _comovingPolar_loglin_,

      /// 2D pair in comoving polar coordinates (r, &mu;) in logarithmic-logarithmic bins
      _comovingPolar_loglog_

    };

  
    /**
     *  @class Pair Pair.h "Headers/Lib/Pair.h"
     *
     *  @brief The class Pair
     *
     *  This class is used to handle objects of type <EM> Pair
     *  </EM>. It contains all virtual member functions implemented in the
     *  derived classes 
     */
    class Pair {

    private:
      
      /**
       *  @name Member functions used to set the binning parameters (customized in all the derived classes) 
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  @return none
       */
      virtual void m_set_parameters_nbins () = 0;
  
      /**
       *  @brief set the binning parameters given the bin size
       *  @return none
       */
      virtual void m_set_parameters_binSize () = 0;
  
      ///@}


    protected:
      
      /// pair dimension
      Dim m_pairDim;
              
      /// pair type
      PairType m_pairType;

      /// angular units
      CoordUnits m_angularUnits;
  
      /// angular weight function
      function<double(double)> m_angularWeight;
      
      /// 

      
    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Pair
       */
      Pair () = default;

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Pair () = default;
    
      /**
       *  @brief static factory used to construct pairs of any type
       *  
       *  @param type the pair type; it can be: _angular_lin_,
       *  _angular_log_, _comoving_lin_, _comoving_log_
       *
       *  @param Min minimum value of the separation (comoving or
       *  angular) used to count the pairs
       *  @param Max maximum value of the separation (comoving or
       *  angular) used to count the pairs
       *  @param nbins number of bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *
       *  @return a pointer to an object of class Pair of a given type
       */
      static shared_ptr<Pair> Create (const PairType type, const double Min, const double Max, const int nbins, const double shift, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr);

      /**
       *  @brief static factory used to construct pairs of any type
       *
       *  @param type the pair type; it can be: _angular_lin_,
       *  _angular_log_, _comoving_lin_, _comoving_log_
       *
       *  @param Min minimum value of the separation (comoving or
       *  angular) used to count the pairs
       *  @param Max maximum value of the separation (comoving or
       *  angular) used to count the pairs
       *  @param binSize the bin size
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *
       *  @return a pointer to an object of class Pair of a given type
       */
      static shared_ptr<Pair> Create (const PairType type, const double Min, const double Max, const double binSize, const double shift, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr);

      /**
       *  @brief static factory used to construct pairs of any type
       *
       *  @param type the pair type; it can be:
       *  D_comovingCartesian_linlin_, D_comovingCartesian_linlog_,
       *  D_comovingCartesian_loglin_, D_comovingCartesian_loglog_,
       *  D_comovingPolar_linlin_, D_comovingPolar_linlog_,
       *  D_comovingPolar_loglin_, D_comovingPolar_loglog_
       *
       *  @param Min_D1 minimum value of the separation (comoving or
       *  angular) in the first direction used to count the pairs
       *  @param Max_D1 maximum value of the separation (comoving or
       *  angular) in the first direction used to count the pairs
       *  @param nbins_D1 number of bins in the first direction
       *  @param shift_D1 shift parameter in the first direction,
       *  i.e. the radial shift is binSize*shift
       *
       *  @param Min_D2 minimum value of the separation (comoving or
       *  angular) in the second direction used to count the pairs
       *  @param Max_D2 maximum value of the separation (comoving or
       *  angular) in the second direction used to count the pairs
       *  @param nbins_D2 number of bins in the second direction
       *  @param shift_D2 shift parameter in the second direction,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *
       *  @return a pointer to an object of class Pair of a given type
       */
      static shared_ptr<Pair> Create (const PairType type, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const double Min_D2, const double Max_D2, const int nbins_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr);

      /**
       *  @brief static factory used to construct pairs of any type
       *
       *  @param type the pair type; it can be:
       *  D_comovingCartesian_linlin_, D_comovingCartesian_linlog_,
       *  D_comovingCartesian_loglin_, D_comovingCartesian_loglog_,
       *  D_comovingPolar_linlin_, D_comovingPolar_linlog_,
       *  D_comovingPolar_loglin_, D_comovingPolar_loglog_
       *
       *  @param Min_D1 minimum value of the separation (comoving or
       *  angular) in the first direction used to count the pairs
       *  @param Max_D1 maximum value of the separation (comoving or
       *  angular) in the first direction used to count the pairs
       *  @param binSize_D1 the bin size in the first direction
       *  @param shift_D1 shift parameter in the first direction,
       *  i.e. the radial shift is binSize*shift
       *
       *  @param Min_D2 minimum value of the separation (comoving or
       *  angular) in the second direction used to count the pairs
       *  @param Max_D2 maximum value of the separation (comoving or
       *  angular) in the second direction used to count the pairs
       *  @param binSize_D2 the bin size in the second direction
       *  @param shift_D2 shift parameter in the second direction,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *
       *  @return a pointer to an object of class Pair of a given type
       */
      static shared_ptr<Pair> Create (const PairType type, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const double Min_D2, const double Max_D2, const double binSize_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr);
      
      ///@}

  
      /**
       *  @name Member functions used to get private/protected parameters
       */
      ///@{
      
      /**
       *  @brief get the dimension of the pair vectors
       *  @return the dimension of the pair vectors
       */
      Dim pairDim () const { return m_pairDim; }

      /**
       *  @brief get the pair type
       *  @return the pair type
       */
      PairType pairType () const { return m_pairType; }

      /**
       *  @brief get the angular units
       *  @return the angular units
       */
      CoordUnits angularUnits () { return m_angularUnits; }

      /**
       *  @brief get the m_angularWeight function
       *  @return the m_angularWeight function 
       */
      function<double(double)> angularWeight () { return m_angularWeight; }
       
      /**
       *  @brief get the member m_scale[i]
       *  @param i the bin index
       *  @return the i-th binned scale
       */
      virtual double scale (const int i) const
      { cosmobl::ErrorMsg("Error in double scale(i) of Pair.h!"); return 0; }

      /**
       *  @brief get the member vector<double> m_scale
       *  @return the vector containing the binned scales
       */
      virtual vector<double> scale () const
      { cosmobl::ErrorMsg("Error in vector<double> scale() of Pair.h!"); vector<double> vv; return vv; }
      
      /**
       *  @brief get the member m_PP1D[i]
       *  @param i the bin index
       *  @return the number of pairs in the i-th bin
       */
      virtual double PP1D (const int i) const
      { cosmobl::ErrorMsg("Error in double PP1D(i) of Pair.h!"); return 0; }

      /**
       *  @brief get the member vector<double> m_PP1D
       *  @return the vector containing the binned number of pairs
       */
      virtual vector<double> PP1D () const
      { cosmobl::ErrorMsg("Error in vector<double> PP1D() of Pair.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the member m_scale_D1[i]
       *  @param i the bin index
       *  @return the i-th binned scale in the first dimension
       */
      virtual double scale_D1 (const int i) const
      { cosmobl::ErrorMsg("Error in double scale(i) of Pair.h!"); return 0; }

      /**
       *  @brief get the member vector<double> m_scale_D1
       *  @return the vector containing the binned scales in the first
       *  dimension
       */
      virtual vector<double> scale_D1 () const
      { cosmobl::ErrorMsg("Error in vector<double> scale() of Pair.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the member m_scale_D2[i]
       *  @param i the bin index
       *  @return the i-th binned scale in the second dimension
       */
      virtual double scale_D2 (const int i) const
      { cosmobl::ErrorMsg("Error in double scale(i) of Pair.h!"); return 0; }

      /**
       *  @brief get the member vector<double> m_scale_D2
       *  @return the vector containing the binned scales in the
       *  second dimension
       */
      virtual vector<double> scale_D2 () const
      { cosmobl::ErrorMsg("Error in vector<double> scale() of Pair.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the member m_PP2D[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the number of pairs in the i-th bin
       */
      virtual double PP2D (const int i, const int j) const
      { cosmobl::ErrorMsg("Error in double PP2D(i,j) of Pair.h!"); return 0;}

      /**
       *  @brief get the member vector<vector<double>> m_PP2D
       *  @return the vector containing the binned number of pairs
       */
      virtual vector<vector<double> > PP2D () const
      { cosmobl::ErrorMsg("Error in double PP2D(i,j) of Pair.h!"); vector<vector<double> > vv; return vv; }
      
      /**
       *  @brief get the member m_binSize_inv
       *  @return the inverse of the bin size
       */
      virtual double binSize_inv () const 
      { cosmobl::ErrorMsg("Error in binSize_inv() of Pair.h!"); return 0; }
  
      /**
       *  @brief get the member m_nbins
       *  @return the number of bins
       */
      virtual int nbins () const
      { cosmobl::ErrorMsg("Error in nbins() of Pair.h!"); return 0; }
    
      /**
       *  @brief get the member m_shift
       *  @return the radial shift used to centre the output bins
       */
      virtual double shift () const
      { cosmobl::ErrorMsg("Error in shift() of Pair.h!"); return 0; }
  
      /**
       *  @brief get the member m_binSize_inv_D1
       *  @return the inverse of the bin size in the first dimension
       */
      virtual double binSize_inv_D1 () const 
      { cosmobl::ErrorMsg("Error in binSize_inv_D1() of Pair.h!"); return 0; }
  
      /**
       *  @brief get the member m_nbins_D1
       *  @return the number of bins in the first dimension
       */
      virtual int nbins_D1 () const 
      { cosmobl::ErrorMsg("Error in nbins_D1() of Pair.h!"); return 0; }
  
      /**
       *  @brief get the member m_shift_D1
       *  @return the radial shift in the first dmension used
       *  to centre the output bins
       */
      virtual double shift_D1 () const 
      { cosmobl::ErrorMsg("Error in shift_D1() of Pair.h!"); return 0; }
  
      /**
       *  @brief get the member m_binSize_inv_D2
       *  @return the inverse of the bin size in the second dimension
       */
      virtual double binSize_inv_D2 () const 
      { cosmobl::ErrorMsg("Error in binSize_inv_D2() of Pair.h!"); return 0; }
  
      /**
       *  @brief get the member m_nbins_D2
       *  @return the number of bins in the second dimension
       */
      virtual int nbins_D2 () const 
      { cosmobl::ErrorMsg("Error in nbins_D2() of Pair.h!"); return 0; }
  
      /**
       *  @brief get the member m_shift_D2
       *  @return the radial shift in the second dimension used to centre
       *  the output bins
       */
      virtual double shift_D2 () const 
      { cosmobl::ErrorMsg("Error in shift_D2() of Pair.h!"); return 0; }

      /**
       *  @brief get the minimum separation
       *  @return the minimum separation used to count the pairs
       */
      virtual double sMin () const
      { cosmobl::ErrorMsg("Error in sMin() of Pair.h!"); return 0; }

      /**
       *  @brief get the the minimum separation
       *  @return the maximum separation used to count the pairs
       */
      virtual double sMax () const
      { cosmobl::ErrorMsg("Error in sMax() of Pair.h!"); return 0; }
      
      /**
       *  @brief get the minimum separation in the first dimension
       *  @return the minimum separation in the first dimension used
       *  to count the pairs
       */
      virtual double sMin_D1 () const
      { cosmobl::ErrorMsg("Error in sMin_D1() of Pair.h!"); return 0; }

      /**
       *  @brief get the the minimum separation in the first dimension
       *  @return the maximum separation in the first dimension used
       *  to count the pairs
       */
      virtual double sMax_D1 () const
      { cosmobl::ErrorMsg("Error in sMax_D1() of Pair.h!"); return 0; }

      /**
       *  @brief get the minimum separation in the second dimension
       *  @return the minimum separation in the second dimension used
       *  to count the pairs
       */
      virtual double sMin_D2 () const
      { cosmobl::ErrorMsg("Error in sMin_D2() of Pair.h!"); return 0; }

      /**
       *  @brief get the the minimum separation in the second dimension
       *  @return the maximum separation in the second dimension used
       *  to count the pairs
       */
      virtual double sMax_D2 () const
      { cosmobl::ErrorMsg("Error in sMax_D2() of Pair.h!"); return 0; }
      
      ///@}
    

      /**
       *  @name Member functions used to set private/protected members
       */
      ///@{

      /**
       *  @brief set the member m_PP1D[i]
       *  @param i the bin index
       *  @param pp the number of pairs in the bin
       *  @return none
       */
      virtual void set_PP1D (const int i, const double pp)
      { cosmobl::ErrorMsg("Error in set_PP1D() of Pair.h!"); }

      /**
       *  @brief set the protected member Pair1D::m_PP1D[i] adding the
       *  number of pairs
       *  @param i the bin index
       *  @param pp the number of pairs in the bin
       *  @return none
       */
      virtual void add_PP1D (const int i, const double pp)
      { cosmobl::ErrorMsg("Error in add_PP1D() of Pair.h!"); }
      
      /**
       *  @brief set the member m_PP2D[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @param pp the number of pairs in the bin
       *  @return none
       */
      virtual void set_PP2D (const int i, const int j, const double pp)
      { cosmobl::ErrorMsg("Error in set_PP2D() of Pair.h!"); }

      /**
       *  @brief set the protected member Pair1D::m_PP2D[i][j] adding
       *  the number of pairs
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @param pp the number of pairs in the bin
       *  @return none
       */
      virtual void add_PP2D (const int i, const int j, const double pp)
      { cosmobl::ErrorMsg("Error in add_PP2D() of Pair.h!"); }
      
      ///@}

  
      /**
       *  @name Member functions used to handle object pairs (customized in all the derived classes) 
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      virtual void put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2) = 0;
  
      /**
       *  @brief sum the number of binned pairs
       *  @param pp an object of class Pair
       *  @param ww the weight
       *  @return none
       */
      virtual void Sum (const shared_ptr<Pair> pp, const double ww=1) = 0;

      ///@}

    };
  }
}

#endif
