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
   *  @brief The namespace of the functions and classes used to handle
   *  <B> pairs of objects </B> 
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
     * @enum PairInfo
     * @brief the information contained in the pairs
     */
    enum PairInfo { 

      /// standard: the object contains only the number of pairs 
      _standard_,
      
      /// extra: the object contains the number of pairs plus extra information, such as the mean scale separation and redshift
      _extra_
      
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

      /// pair info
      PairInfo m_pairInfo;

      /// angular units
      CoordUnits m_angularUnits;
  
      /// angular weight function
      function<double(double)> m_angularWeight;

      
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
       *  @param type the pair type; it can be: \_angular_lin\_,
       *  \_angular_log\_, \_comoving_lin\_, \_comoving_log\_
       *
       *  @param info the pair information; it can be: \_standard\_ or
       *  \_extra\_
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
      static shared_ptr<Pair> Create (const PairType type, const PairInfo info, const double Min, const double Max, const int nbins, const double shift, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr);

      /**
       *  @brief static factory used to construct pairs of any type
       *
       *  @param type the pair type; it can be: \_angular_lin\_,
       *  \_angular_log\_, \_comoving_lin\_, \_comoving_log\_
       *
       *  @param info the pair information; it can be: \_standard\_ or
       *  \_extra\_
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
      static shared_ptr<Pair> Create (const PairType type, const PairInfo info, const double Min, const double Max, const double binSize, const double shift, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr);

      /**
       *  @brief static factory used to construct pairs of any type
       *
       *  @param type the pair type; it can be:
       *  \_comovingCartesian_linlin\_, \_comovingCartesian_linlog\_,
       *  \_comovingCartesian_loglin\_, \_comovingCartesian_loglog\_,
       *  \_comovingPolar_linlin\_, \_comovingPolar_linlog\_,
       *  \_comovingPolar_loglin\_, \_comovingPolar_loglog\_
       *
       *  @param info the pair information; it can be: \_standard\_ or
       *  \_extra\_
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
      static shared_ptr<Pair> Create (const PairType type, const PairInfo info, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const double Min_D2, const double Max_D2, const int nbins_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr);

      /**
       *  @brief static factory used to construct pairs of any type
       *
       *  @param type the pair type; it can be:
       *  \_comovingCartesian_linlin\_, \_comovingCartesian_linlog\_,
       *  \_comovingCartesian_loglin\_, \_comovingCartesian_loglog\_,
       *  \_comovingPolar_linlin\_, \_comovingPolar_linlog\_,
       *  \_comovingPolar_loglin\_, \_comovingPolar_loglog\_
       *
       *  @param info the pair information; it can be: \_standard\_ or
       *  \_extra\_
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
      static shared_ptr<Pair> Create (const PairType type, const PairInfo info, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const double Min_D2, const double Max_D2, const double binSize_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr);
      
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
       *  @brief get the pair information type
       *  @return the pair information type
       */
      PairInfo pairInfo () const { return m_pairInfo; }

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
      { (void)i; cosmobl::ErrorCBL("Error in scale() of Pair.h!"); return 0; }

      /**
       *  @brief get the member vector<double> m_scale
       *  @return the vector containing the binned scales
       */
      virtual vector<double> scale () const
      { cosmobl::ErrorCBL("Error in scale() of Pair.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the member m_scale_mean[i]
       *  @param i the bin index
       *  @return the mean scale in the i-th bin
       */
      virtual double scale_mean (const int i) const
      { (void)i; cosmobl::ErrorCBL("Error in scale_mean() of Pair.h!"); return 0; }

      /**
       *  @brief get the member vector<double> m_scale_mean
       *  @return the vector containing the mean scales 
       */
      virtual vector<double> scale_mean () const
      { cosmobl::ErrorCBL("Error in scale_mean() of Pair.h!"); vector<double> vv; return vv; }
      
      /**
       *  @brief get the member m_scale_sigma[i]
       *  @param i the bin index
       *  @return the standard deviation of the scale distribution in
       *  the i-th bin
       */
      virtual double scale_sigma (const int i) const
      { (void)i; cosmobl::ErrorCBL("Error in scale_sigma() of Pair.h!"); return 0; }

      /**
       *  @brief get the member vector<double> m_scale_sigma
       *  @return the vector containing the standard deviations of the
       *  scale distribution
       */
      virtual vector<double> scale_sigma () const
      { cosmobl::ErrorCBL("Error in scale_sigma() of Pair.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the protected member Pair1D_extra::m_z_mean[i]
       *  @param i the bin index
       *  @return the mean redshift in the i-th bin
       */
      virtual double z_mean (const int i) const
      { (void)i; cosmobl::ErrorCBL("Error in z_mean() of Pair.h!"); return 0; }

      /**
       *  @brief get the protected member Pair1D_extra::m_z_mean
       *  @return the vector containing the mean redshifts 
       */
      virtual vector<double> z_mean () const
      { cosmobl::ErrorCBL("Error in z_mean() of Pair.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the protected member Pair1D_extra::m_z_sigma[i]
       *  @param i the bin index
       *  @return the standard deviation of the redshift distribution
       *  in the i-th bin
       */
      virtual double z_sigma (const int i) const
      { (void)i; cosmobl::ErrorCBL("Error in z_sigma() of Pair.h!"); return 0; }

      /**
       *  @brief get the protected member Pair1D_extra::m_z_sigma
       *  @return the vector containing the standard deviations of the
       *  redshift distribution
       */
      virtual vector<double> z_sigma () const
      { cosmobl::ErrorCBL("Error in z_sigma() of Pair.h!"); vector<double> vv; return vv; }
      
      /**
       *  @brief get the member m_PP1D[i]
       *  @param i the bin index
       *  @return the number of pairs in the i-th bin
       */
      virtual double PP1D (const int i) const
      { (void)i; cosmobl::ErrorCBL("Error in PP1D() of Pair.h!"); return 0; }

      /**
       *  @brief get the member vector<double> m_PP1D
       *  @return the vector containing the binned number of pairs
       */
      virtual vector<double> PP1D () const
      { cosmobl::ErrorCBL("Error in PP1D() of Pair.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the member m_scale_D1[i]
       *  @param i the bin index
       *  @return the i-th binned scale in the first dimension
       */
      virtual double scale_D1 (const int i) const
      { (void)i; cosmobl::ErrorCBL("Error in scale_D1() of Pair.h!"); return 0; }

      /**
       *  @brief get the member vector<double> m_scale_D1
       *  @return the vector containing the binned scales in the first
       *  dimension
       */
      virtual vector<double> scale_D1 () const
      { cosmobl::ErrorCBL("Error in scale_D1() of Pair.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the member m_scale_D2[i]
       *  @param i the bin index
       *  @return the i-th binned scale in the second dimension
       */
      virtual double scale_D2 (const int i) const
      { (void)i; cosmobl::ErrorCBL("Error in scale_D2() of Pair.h!"); return 0; }

      /**
       *  @brief get the member vector<double> m_scale_D2
       *  @return the vector containing the binned scales in the
       *  second dimension
       */
      virtual vector<double> scale_D2 () const
      { cosmobl::ErrorCBL("Error in scale_D2() of Pair.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the protected member \e m_scale_D1[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the mean scale in the first dimension, in the i-j bin
       */
      virtual double scale_D1_mean (const int i, const int j) const
      { (void)i; (void)j; cosmobl::ErrorCBL("Error in scale_D1_mean() of Pair.h!"); return 0; }

      /**
       *  @brief get the protected member \e m_scale_D1
       *  @return the matrix containing the mean scales in the first dimension
       */
      virtual vector<vector<double>> scale_D1_mean () const
      { cosmobl::ErrorCBL("Error in scale_D1_mean() of Pair.h!"); vector<vector<double>> vv; return vv; }

      /**
       *  @brief get the protected member \e m_scale_D2[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the mean scale in the second dimension, in the i-j bin
       */
      virtual double scale_D2_mean (const int i, const int j) const
      { (void)i; (void)j; cosmobl::ErrorCBL("Error in scale_D2_mean() of Pair.h!"); return 0; }

      /**
       *  @brief get the protected member \e m_scale_D2
       *  @return the matrix containing the mean scales in the second
       *  dimension
       */
      virtual vector<vector<double>> scale_D2_mean () const
      { cosmobl::ErrorCBL("Error in scale_D2_mean() of Pair.h!"); vector<vector<double>> vv; return vv; }

      /**
       *  @brief get the member m_scale_D1_sigma[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the standard deviation of the scale distribution in
       *  the i-th bin, in the first dimension
       */
      virtual double scale_D1_sigma (const int i, const int j) const
      { (void)i; (void)j; cosmobl::ErrorCBL("Error in scale_D1_sigma() of Pair.h!"); return 0; }

      /**
       *  @brief get the member vector<double> m_scale_D1
       *  @return the vector containing the standard deviations of the
       *  scale distribution in the first dimension
       */
      virtual vector<vector<double>> scale_D1_sigma () const
      { cosmobl::ErrorCBL("Error in scale_D1_sigma() of Pair.h!"); vector<vector<double>> vv; return vv; }

      /**
       *  @brief get the member m_scale_D2[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the standard deviation of the scale distribution in
       *  the i-th bin, in the second dimension
       */
      virtual double scale_D2_sigma (const int i, const int j) const
      { (void)i; (void)j; cosmobl::ErrorCBL("Error in scale_D2_sigma() of Pair.h!"); return 0; }

      /**
       *  @brief get the member vector<double> m_scale_D2
       *  @return the vector containing the standard deviations of the
       *  scale distribution in the second dimension
       */
      virtual vector<vector<double>> scale_D2_sigma () const
      { cosmobl::ErrorCBL("Error in scale_D2_sigma() of Pair.h!"); vector<vector<double>> vv; return vv; }

      /**
       *  @brief get the protected member \e m_z_mean[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the mean redshift, in the i-j bin
       */
      virtual double z_mean (const int i, const int j) const
      { (void)i; (void)j; cosmobl::ErrorCBL("Error in z_mean() of Pair.h!"); return 0; }

      /**
       *  @brief get the protected member \e m_z_mean
       *  @return the matrix containing the mean redshifts
       */
      virtual vector<vector<double>> z_mean2D () const
      { cosmobl::ErrorCBL("Error in z_mean2D() of Pair.h!"); vector<vector<double>> vv; return vv; }

      /**
       *  @brief get the member m_z_sigma[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the standard deviation of the redshift distribution
       *  in the i-th bin
       */
      virtual double z_sigma (const int i, const int j) const
      { (void)i; (void)j; cosmobl::ErrorCBL("Error in z_sigma() of Pair.h!"); return 0; }

      /**
       *  @brief get the member vector<double> m_z_sigma
       *  @return the vector containing the standard deviations of the
       *  redshift distribution
       */
      virtual vector<vector<double>> z_sigma2D () const 
      { cosmobl::ErrorCBL("Error in z_sigma2D() of Pair.h!"); vector<vector<double>> vv; return vv; }

      /**
       *  @brief get the member m_PP2D[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the number of pairs in the i-th bin
       */
      virtual double PP2D (const int i, const int j) const
      { (void)i; (void)j; cosmobl::ErrorCBL("Error in PP2D() of Pair.h!"); return 0;}

      /**
       *  @brief get the member vector<vector<double>> m_PP2D
       *  @return the vector containing the binned number of pairs
       */
      virtual vector<vector<double>> PP2D () const
      { cosmobl::ErrorCBL("Error in PP2D() of Pair.h!"); vector<vector<double>> vv; return vv; }
      
      /**
       *  @brief get the member m_binSize_inv
       *  @return the inverse of the bin size
       */
      virtual double binSize_inv () const 
      { cosmobl::ErrorCBL("Error in binSize_inv() of Pair.h!"); return 0; }
  
      /**
       *  @brief get the member m_nbins
       *  @return the number of bins
       */
      virtual int nbins () const
      { cosmobl::ErrorCBL("Error in nbins() of Pair.h!"); return 0; }
    
      /**
       *  @brief get the member m_shift
       *  @return the radial shift used to centre the output bins
       */
      virtual double shift () const
      { cosmobl::ErrorCBL("Error in shift() of Pair.h!"); return 0; }
  
      /**
       *  @brief get the member m_binSize_inv_D1
       *  @return the inverse of the bin size in the first dimension
       */
      virtual double binSize_inv_D1 () const 
      { cosmobl::ErrorCBL("Error in binSize_inv_D1() of Pair.h!"); return 0; }
  
      /**
       *  @brief get the member m_nbins_D1
       *  @return the number of bins in the first dimension
       */
      virtual int nbins_D1 () const 
      { cosmobl::ErrorCBL("Error in nbins_D1() of Pair.h!"); return 0; }
  
      /**
       *  @brief get the member m_shift_D1
       *  @return the radial shift in the first dmension used
       *  to centre the output bins
       */
      virtual double shift_D1 () const 
      { cosmobl::ErrorCBL("Error in shift_D1() of Pair.h!"); return 0; }
  
      /**
       *  @brief get the member m_binSize_inv_D2
       *  @return the inverse of the bin size in the second dimension
       */
      virtual double binSize_inv_D2 () const 
      { cosmobl::ErrorCBL("Error in binSize_inv_D2() of Pair.h!"); return 0; }
  
      /**
       *  @brief get the member m_nbins_D2
       *  @return the number of bins in the second dimension
       */
      virtual int nbins_D2 () const 
      { cosmobl::ErrorCBL("Error in nbins_D2() of Pair.h!"); return 0; }
  
      /**
       *  @brief get the member m_shift_D2
       *  @return the radial shift in the second dimension used to centre
       *  the output bins
       */
      virtual double shift_D2 () const 
      { cosmobl::ErrorCBL("Error in shift_D2() of Pair.h!"); return 0; }

      /**
       *  @brief get the minimum separation
       *  @return the minimum separation used to count the pairs
       */
      virtual double sMin () const
      { cosmobl::ErrorCBL("Error in sMin() of Pair.h!"); return 0; }

      /**
       *  @brief get the the minimum separation
       *  @return the maximum separation used to count the pairs
       */
      virtual double sMax () const
      { cosmobl::ErrorCBL("Error in sMax() of Pair.h!"); return 0; }
      
      /**
       *  @brief get the minimum separation in the first dimension
       *  @return the minimum separation in the first dimension used
       *  to count the pairs
       */
      virtual double sMin_D1 () const
      { cosmobl::ErrorCBL("Error in sMin_D1() of Pair.h!"); return 0; }

      /**
       *  @brief get the the minimum separation in the first dimension
       *  @return the maximum separation in the first dimension used
       *  to count the pairs
       */
      virtual double sMax_D1 () const
      { cosmobl::ErrorCBL("Error in sMax_D1() of Pair.h!"); return 0; }

      /**
       *  @brief get the minimum separation in the second dimension
       *  @return the minimum separation in the second dimension used
       *  to count the pairs
       */
      virtual double sMin_D2 () const
      { cosmobl::ErrorCBL("Error in sMin_D2() of Pair.h!"); return 0; }

      /**
       *  @brief get the the minimum separation in the second dimension
       *  @return the maximum separation in the second dimension used
       *  to count the pairs
       */
      virtual double sMax_D2 () const
      { cosmobl::ErrorCBL("Error in sMax_D2() of Pair.h!"); return 0; }
      
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
      { (void)i; (void)pp; cosmobl::ErrorCBL("Error in set_PP1D() of Pair.h!"); }
      
      /**
       *  @brief set the member m_PP2D[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @param pp the number of pairs in the bin
       *  @return none
       */
      virtual void set_PP2D (const int i, const int j, const double pp)
      { (void)i; (void)j; (void)pp; cosmobl::ErrorCBL("Error in set_PP2D() of Pair.h!"); }

       /**
       *  @brief set the protected members by adding new 1D data
       *  @param i the bin index
       *  @param data vector containing the new data to be added
       *  @return none
       */
      virtual void add_data1D (const int i, const vector<double> data)
      { (void)i; (void)data; cosmobl::ErrorCBL("Error in void add_data1D() of Pair.h!"); }

      /**
       *  @brief set the protected members by adding new 1D data
       *  @param i the bin index
       *  @param pair pointer to an object of class Pair
       *  @param ww a multiplicative factor used for bootstrap
       *  @return none
       */
      virtual void add_data1D (const int i, const shared_ptr<pairs::Pair> pair, const double ww=1.)
      { (void)i; (void)pair; (void)ww; cosmobl::ErrorCBL("Error in add_data1D() of Pair.h!"); }
      
      /**
       *  @brief set the protected members by adding new 2D data
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @param data vector containing the new data to be added
       *  @return none
       */
      virtual void add_data2D (const int i, const int j, const vector<double> data)
      { (void)i; (void)j; (void)data; cosmobl::ErrorCBL("Error in add_data2D() of Pair.h!"); }

      /**
       *  @brief set the protected members by adding new 2D data
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @param pair pair pointer to an object of class Pair
       *  @param ww a multiplicative factor used for bootstrap
       *  @return none
       */
      virtual void add_data2D (const int i, const int j, const shared_ptr<pairs::Pair> pair, const double ww=1.)
      { (void)i; (void)j; (void)pair; (void)ww; cosmobl::ErrorCBL("Error in add_data2D() of Pair.h!"); }
      
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

      /**
       *  @brief finalise the computation of the extra information
       *  @return none
       */
      virtual void finalise ()
      { cosmobl::ErrorCBL("Error in finalise() of Pair.h!"); }
      
      ///@}

    };
  }
}

#endif
