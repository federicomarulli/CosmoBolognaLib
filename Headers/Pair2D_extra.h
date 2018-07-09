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
 *  @file Headers/Pair2D_extra.h
 *
 *  @brief The classes Pair2D_extra*
 *
 *  This file defines the interface of all the classes Pair2D_extra*,
 *  used to handle 2D pairs of objects of any kind, with extra
 *  information, such as the mean redshift and separation
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unbo.it
 */

#ifndef __PAIR2D_EXTRA__
#define __PAIR2D_EXTRA__


#include "Pair2D.h"


// ===================================================================================================


namespace cbl {
  
  namespace pairs {

    /**
     *  @class Pair2D_extra Pair2D_extra.h "Headers/Pair2D_extra.h"
     *
     *  @brief The class Pair2D_extra
     *
     *  This class is used to handle objects of type <EM> Pair2D_extra
     *  </EM>.
     */
    class Pair2D_extra : public virtual Pair2D {
      
    protected:

      /// the mean scales in each bin, in the first dimension
      std::vector<std::vector<double>> m_scale_D1_mean;
      
      /// the mean scales in each bin, in the second dimension
      std::vector<std::vector<double>> m_scale_D2_mean;

      /// the square of the standard deviations of the scale distributions in each bin, in the first dimension, multiplied by the total weight
      std::vector<std::vector<double>> m_scale_D1_S;
      
      /// the square of the standard deviations of the scale distributions in each bin, in the second dimension, multiplied by the total weight
      std::vector<std::vector<double>> m_scale_D2_S;
      
      /// the standard deviations of the scale distributions in each bin, in the first dimension
      std::vector<std::vector<double>> m_scale_D1_sigma;
      
      /// the standard deviations of the scale distributions in each bin, in the second dimension
      std::vector<std::vector<double>> m_scale_D2_sigma;

      /// the mean redshifts in each bin
      std::vector<std::vector<double>> m_z_mean;

      /// the square of the standard deviations of the redshift distributions in each bin, multiplied by the total weight
      std::vector<std::vector<double>> m_z_S;

      /// the standard deviations of the redshift distributions in each bin
      std::vector<std::vector<double>> m_z_sigma;

      
    public:
  
      /**
       *  @name Constructors/destructors
       */
      ///@{
  
      /**
       *  @brief default constructor
       *  @return object of class Pair2D_extra
       */
      Pair2D_extra () = default;

      /**
       *  @brief constructor
       *  @param binSize_D1 the bin size in the first dimension
       *  @param nbins_D1 number of bins in the first dimension
       *  @param shift_D1 radial shift used to centre the output bins
       *  in the first dimension
       *  @param binSize_D2 the bin size in the second dimension
       *  @param nbins_D2 number of bins in the second dimension
       *  @param shift_D2 radial shift used to centre the output bins
       *  in the second dimension
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_extra
       */
      Pair2D_extra (const double binSize_D1, const int nbins_D1, const double shift_D1, const double binSize_D2, const int nbins_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(binSize_D1, nbins_D1, shift_D1, binSize_D2, nbins_D2, shift_D2, angularUnits, angularWeight) {}
  
      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Pair2D_extra () = default;
  
      ///@}
  
  
      /**
       *  @name Member functions to get the protected members 
       */
      ///@{

      /**
       *  @brief get the protected member \e m_scale_D1_mean[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the mean scale in the first dimension, in the i-j bin
       */
      double scale_D1_mean (const int i, const int j) const override { return m_scale_D1_mean[i][j]; }

      /**
       *  @brief get the protected member \e m_scale_D1_mean
       *  @return the matrix containing the mean scales in the first
       *  dimension
       */
      std::vector<std::vector<double>> scale_D1_mean () const override { return m_scale_D1_mean; }

      /**
       *  @brief get the protected member \e m_scale_D2_mean[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the mean scale in the second dimension, in the i-j bin
       */
      double scale_D2_mean (const int i, const int j) const override { return m_scale_D2_mean[i][j]; }

      /**
       *  @brief get the protected member \e m_scale_D2_mean
       *  @return the matrix containing the mean scales in the second
       *  dimension
       */
      std::vector<std::vector<double>> scale_D2_mean () const override { return m_scale_D2_mean; }

      /**
       *  @brief get the member m_scale_D1_S[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the square of the standard deviations of the scale
       *  distribution, multiplied by the total weight, in the i-th
       *  bin, in the first dimension
       */
      double scale_D1_S (const int i, const int j) const override { return m_scale_D1_S[i][j]; }

      /**
       *  @brief get the member std::vector<double> m_scale_D1_S
       *  @return the vector the square of the standard deviations of
       *  the scale distribution, multiplied by the total weight, in
       *  the first dimension
       */
      std::vector<std::vector<double>> scale_D1_S () const override { return m_scale_D1_S; }

      /**
       *  @brief get the member m_scale_D2_S[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the square of the standard deviations of the scale
       *  distribution, multiplied by the total weight, in the i-th
       *  bin, in the second dimension
       */
      double scale_D2_S (const int i, const int j) const override { return m_scale_D2_S[i][j]; }

      /**
       *  @brief get the member std::vector<double> m_scale_D2_S
       *  @return the vector the square of the standard deviations of
       *  the scale distribution, multiplied by the total weight, in
       *  the second dimension
       */
      std::vector<std::vector<double>> scale_D2_S () const override { return m_scale_D2_S; }

      /**
       *  @brief get the member m_scale_D1_sigma[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the standard deviation of the scale distribution in
       *  the i-th bin, in the first dimension
       */
      double scale_D1_sigma (const int i, const int j) const override { return m_scale_D1_sigma[i][j]; }

      /**
       *  @brief get the member std::vector<double> m_scale_D1_sigma
       *  @return the vector containing the standard deviations of the
       *  scale distribution in the first dimension
       */
      std::vector<std::vector<double>> scale_D1_sigma () const override { return m_scale_D1_sigma; }

      /**
       *  @brief get the member m_scale_D2_sigma[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the standard deviation of the scale distribution in
       *  the i-th bin, in the second dimension
       */
      double scale_D2_sigma (const int i, const int j) const override { return m_scale_D2_sigma[i][j]; }

      /**
       *  @brief get the member std::vector<double> m_scale_D2_sigma
       *  @return the vector containing the standard deviations of the
       *  scale distribution in the second dimension
       */
      std::vector<std::vector<double>> scale_D2_sigma () const override { return m_scale_D2_sigma; }

      /**
       *  @brief get the protected member \e m_z_mean[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the mean redshift, in the i-j bin
       */
      double z_mean (const int i, const int j) const override { return m_z_mean[i][j]; }

      /**
       *  @brief get the member m_z_S[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the square of the standard deviations of the
       *  redshift distribution, multiplied by the total weight, in
       *  the i-th bin
       */
      double z_S (const int i, const int j) const override { return m_z_S[i][j]; }
      
      /**
       *  @brief get the member m_z_sigma[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the standard deviation of the redshift distribution
       *  in the i-th bin
       */
      double z_sigma (const int i, const int j) const override { return m_z_sigma[i][j]; }

      ///@}

    
      /**
       *  @name Member functions used to set the protected members
       */
      ///@{

      /**
       *  @brief set the protected member Pair2D::m_scale_D1[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @param ss the mean scale in the first dimension
       *  @return none
       */
      void set_scale_D1_mean (const int i, const int j, const double ss) { checkDim(m_scale_D1_mean, i, j, "m_scale_D1_mean"); m_scale_D1_mean[i][j] = ss; }

      /**
       *  @brief set the protected member Pair2D::m_scale_D2[i][j]
       *  @param i the bin index in the second dimension
       *  @param j the bin index in the second dimension
       *  @param ss the mean scale in the second dimension
       *  @return none
       */
      void set_scale_D2_mean (const int i, const int j, const double ss) { checkDim(m_scale_D2_mean, i, j, "m_scale_D2_mean"); m_scale_D2_mean[i][j] = ss; }
      
      /**
       *  @brief set the protected member Pair2D::m_scale_D1[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @param ss the standard deviation of the scale distribution
       *  in the first dimension
       *  @return none
       */
      void set_scale_D1_sigma (const int i, const int j, const double ss) { checkDim(m_scale_D1_sigma, i, j, "m_scale_D1_sigma"); m_scale_D1_sigma[i][j] = ss; }

      /**
       *  @brief set the protected member Pair2D::m_scale_D2[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @param ss the standard deviation of the scale distribution
       *  in the second dimension
       *  @return none
       */
      void set_scale_D2_sigma (const int i, const int j, const double ss) { checkDim(m_scale_D2_sigma, i, j, "m_scale_D2_sigma"); m_scale_D2_sigma[i][j] = ss; }

      /**
       *  @brief set the protected member Pair2D::m_z_mean[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @param ss the mean redshift 
       *  @return none
       */
      void set_z_mean (const int i, const int j, const double ss) { checkDim(m_z_mean, i, j, "m_z_mean"); m_z_mean[i][j] = ss; }
      
      /**
       *  @brief set the protected member Pair2D::m_z_sigma[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @param ss the standard deviation of the redshift
       *  distribution 
       *  @return none
       */
      void set_z_sigma (const int i, const int j, const double ss) { checkDim(m_z_sigma, i, j, "m_z_sigma"); m_z_sigma[i][j] = ss; }

      /**
       *  @brief set the protected members by adding new data
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @param data vector containing the new data to be added
       *  @return none
       */
      void add_data2D (const int i, const int j, const std::vector<double> data) override;
      
      /**
       *  @brief set the protected members by adding new data
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @param pair pair pointer to an object of class Pair
       *  @param ww a multiplicative factor used for bootstrap
       *  @return none
       */
      void add_data2D (const int i, const int j, const std::shared_ptr<pairs::Pair> pair, const double ww=1.) override;
      
      ///@}


      /**
       *  @name Member functions used to handle object pairs
       */
      ///@{
      
      /**
       *  @brief sum the number of binned pairs
       *  @param pair an object of class Pair
       *  @param ww the weight
       *  @return none
       */
      void Sum (const std::shared_ptr<Pair> pair, const double ww=1) override;
      
      ///@}
    
    };


    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingCartesian_extra Pair2D_extra.h
     *  "Headers/Pair2D_extra.h"
     *
     *  @brief The class Pair2D_comovingCartesian_extra
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingCartesian_extra </EM>.
     */
    class Pair2D_comovingCartesian_extra : public virtual Pair2D_extra, public virtual Pair2D_comovingCartesian {

    public:
  
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Pair2D_comovingCartesian_extra
       */
      Pair2D_comovingCartesian_extra () = default;

      /**
       *  @brief constructor   
       *  @param rpMin minimum perpendicular separation used to count
       *  the pairs
       *  @param rpMax maximum perpendicular separation used to count
       *  the pairs
       *  @param nbins_rp number of bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular
       *  separation, i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the
       *  pairs
       *  @param piMax maximum parallel separation used to count the
       *  pairs
       *  @param nbins_pi number of bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian_extra
       */
      Pair2D_comovingCartesian_extra (const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingCartesian(rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight) {}
  
      /**
       *  @brief constructor
       *  @param rpMin minimum perpendicular separation used to count
       *  the pairs
       *  @param rpMax maximum perpendicular separation used to count
       *  the pairs
       *  @param binSize_rp size of the bins in the perpendicular
       *  separation
       *  @param shift_rp shift parameter in the perpendicular
       *  separation, i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param binSize_pi size of the bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian_extra
       */
      Pair2D_comovingCartesian_extra (const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingCartesian(rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight) {}
  
      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Pair2D_comovingCartesian_extra () = default;

      ///@}
    
    };


    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingCartesian_linlin_extra Pair2D_extra.h
     *  "Headers/Pair2D_extra.h"
     *
     *  @brief The class Pair2D_comovingCartesian_linlin_extra
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingCartesian_linlin_extra </EM>.
     */
    class Pair2D_comovingCartesian_linlin_extra : public virtual Pair2D_comovingCartesian_extra, public virtual Pair2D_comovingCartesian_linlin {
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Pair2D_comovingCartesian_linlin_extra
       */
      Pair2D_comovingCartesian_linlin_extra ()
	{
	  m_pairType = PairType::_comovingCartesian_linlin_;
	  m_pairInfo = PairInfo::_extra_;
	} 

      /**
       *  @brief constructor   
       *  @param rpMin minimum perpendicular separation used to count the
       *  pairs
       *  @param rpMax maximum perpendicular separation used to count the
       *  pairs
       *  @param nbins_rp number of bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular separation,
       *  i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param nbins_pi number of bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian_linlin_extra
       */
      Pair2D_comovingCartesian_linlin_extra (const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(1., nbins_rp, shift_rp, 1., nbins_pi, shift_pi, angularUnits, angularWeight), Pair2D_comovingCartesian(rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comovingCartesian_linlin_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_nbins();
	  m_PP2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_PP2D_weighted.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	}
  
      /**
       *  @brief constructor
       *  @param rpMin minimum perpendicular separation used to count the
       *  pairs
       *  @param rpMax maximum perpendicular separation used to count the
       *  pairs
       *  @param binSize_rp size of the bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular separation,
       *  i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param binSize_pi size of the bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian_linlin_extra
       */
      Pair2D_comovingCartesian_linlin_extra (const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(binSize_rp, 50, shift_rp, binSize_pi, 50, shift_pi, angularUnits, angularWeight), Pair2D_comovingCartesian(rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comovingCartesian_linlin_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_binSize();
	  m_PP2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_PP2D_weighted.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	} 
  
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Pair2D_comovingCartesian_linlin_extra () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;
  
      ///@}

    };


    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingCartesian_loglin_extra Pair2D_extra.h
     *  "Headers/Pair2D_extra.h"
     *
     *  @brief The class Pair2D_comovingCartesian_loglin_extra
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingCartesian_loglin_extra </EM>.
     */
    class Pair2D_comovingCartesian_loglin_extra : public virtual Pair2D_comovingCartesian_extra, public virtual Pair2D_comovingCartesian_loglin {
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Pair2D_comovingCartesian_loglin_extra
       */
      Pair2D_comovingCartesian_loglin_extra ()
	{
	  m_pairType = PairType::_comovingCartesian_loglin_;
	  m_pairInfo = PairInfo::_extra_;
	} 

      /**
       *  @brief constructor   
       *  @param rpMin minimum perpendicular separation used to count the
       *  pairs
       *  @param rpMax maximum perpendicular separation used to count the
       *  pairs
       *  @param nbins_rp number of bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular separation,
       *  i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param nbins_pi number of bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian_loglin_extra
       */
      Pair2D_comovingCartesian_loglin_extra (const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(1., nbins_rp, shift_rp, 1., nbins_pi, shift_pi, angularUnits, angularWeight), Pair2D_comovingCartesian(rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comovingCartesian_loglin_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_nbins();
	  m_PP2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_PP2D_weighted.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	}
  
      /**
       *  @brief constructor
       *  @param rpMin minimum perpendicular separation used to count the
       *  pairs
       *  @param rpMax maximum perpendicular separation used to count the
       *  pairs
       *  @param binSize_rp size of the bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular separation,
       *  i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param binSize_pi size of the bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian_loglin_extra
       */
      Pair2D_comovingCartesian_loglin_extra (const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(binSize_rp, 50, shift_rp, binSize_pi, 50, shift_pi, angularUnits, angularWeight), Pair2D_comovingCartesian(rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comovingCartesian_loglin_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_binSize();
	  m_PP2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_PP2D_weighted.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	} 
  
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Pair2D_comovingCartesian_loglin_extra () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;
  
      ///@}

    };

    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingCartesian_linlog_extra Pair2D_extra.h
     *  "Headers/Pair2D_extra.h"
     *
     *  @brief The class Pair2D_comovingCartesian_linlog_extra
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingCartesian_linlog_extra </EM>.
     */
    class Pair2D_comovingCartesian_linlog_extra : public virtual Pair2D_comovingCartesian_extra, public virtual Pair2D_comovingCartesian_linlog {
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Pair2D_comovingCartesian_linlog_extra
       */
      Pair2D_comovingCartesian_linlog_extra ()
	{
	  m_pairType = PairType::_comovingCartesian_linlog_;
	  m_pairInfo = PairInfo::_extra_;
	} 

      /**
       *  @brief constructor   
       *  @param rpMin minimum perpendicular separation used to count the
       *  pairs
       *  @param rpMax maximum perpendicular separation used to count the
       *  pairs
       *  @param nbins_rp number of bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular separation,
       *  i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param nbins_pi number of bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian_linlog_extra
       */
      Pair2D_comovingCartesian_linlog_extra (const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(1., nbins_rp, shift_rp, 1., nbins_pi, shift_pi, angularUnits, angularWeight), Pair2D_comovingCartesian(rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comovingCartesian_linlog_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_nbins();
	  m_PP2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_PP2D_weighted.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	}
  
      /**
       *  @brief constructor
       *  @param rpMin minimum perpendicular separation used to count the
       *  pairs
       *  @param rpMax maximum perpendicular separation used to count the
       *  pairs
       *  @param binSize_rp size of the bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular separation,
       *  i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param binSize_pi size of the bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian_linlog_extra
       */
      Pair2D_comovingCartesian_linlog_extra (const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(binSize_rp, 50, shift_rp, binSize_pi, 50, shift_pi, angularUnits, angularWeight), Pair2D_comovingCartesian(rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comovingCartesian_linlog_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_binSize();
	  m_PP2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_PP2D_weighted.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	} 
  
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Pair2D_comovingCartesian_linlog_extra () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;
  
      ///@}

    };

    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingCartesian_loglog_extra Pair2D_extra.h
     *  "Headers/Pair2D_extra.h"
     *
     *  @brief The class Pair2D_comovingCartesian_loglog_extra
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingCartesian_loglog_extra </EM>.
     */
    class Pair2D_comovingCartesian_loglog_extra : public virtual Pair2D_comovingCartesian_extra, public virtual Pair2D_comovingCartesian_loglog {
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Pair2D_comovingCartesian_loglog_extra
       */
      Pair2D_comovingCartesian_loglog_extra ()
	{
	  m_pairType = PairType::_comovingCartesian_loglog_;
	  m_pairInfo = PairInfo::_extra_;
	} 

      /**
       *  @brief constructor   
       *  @param rpMin minimum perpendicular separation used to count the
       *  pairs
       *  @param rpMax maximum perpendicular separation used to count the
       *  pairs
       *  @param nbins_rp number of bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular separation,
       *  i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param nbins_pi number of bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian_loglog_extra
       */
      Pair2D_comovingCartesian_loglog_extra (const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(1., nbins_rp, shift_rp, 1., nbins_pi, shift_pi, angularUnits, angularWeight), Pair2D_comovingCartesian(rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comovingCartesian_loglog_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_nbins();
	  m_PP2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_PP2D_weighted.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	}
  
      /**
       *  @brief constructor
       *  @param rpMin minimum perpendicular separation used to count the
       *  pairs
       *  @param rpMax maximum perpendicular separation used to count the
       *  pairs
       *  @param binSize_rp size of the bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular separation,
       *  i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param binSize_pi size of the bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian_loglog_extra
       */
      Pair2D_comovingCartesian_loglog_extra (const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(binSize_rp, 50, shift_rp, binSize_pi, 50, shift_pi, angularUnits, angularWeight), Pair2D_comovingCartesian(rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comovingCartesian_loglog_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_binSize();
	  m_PP2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_PP2D_weighted.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	} 
  
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Pair2D_comovingCartesian_loglog_extra () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;
  
      ///@}

    };

  
    // ============================================================================================
    // ============================================================================================

  
    /**
     *  @class Pair2D_comovingPolar_extra Pair2D_extra.h "Headers/Pair2D_extra.h"
     *
     *  @brief The class Pair2D_comovingPolar_extra
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingPolar_extra </EM>.
     */
    class Pair2D_comovingPolar_extra : public virtual Pair2D_extra, public virtual Pair2D_comovingPolar {
   
    public:
  
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Pair2D_comovingPolar_extra
       */
      Pair2D_comovingPolar_extra () = default;

      /**
       *  @brief constructor   
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins_D1 number of bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param nbins_D2 number of bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar_extra
       */
      Pair2D_comovingPolar_extra (const double rMin, const double rMax, const int nbins_D1, const double shift_D1, const double muMin, const double muMax, const int nbins_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingPolar(rMin, rMax, nbins_D1, shift_D1, muMin, muMax, nbins_D2, shift_D2, angularUnits, angularWeight) {}
  
      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize_D1 size of the bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param binSize_D2 size of the bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar_extra
       */
      Pair2D_comovingPolar_extra (const double rMin, const double rMax, const double binSize_D1, const double shift_D1, const double muMin, const double muMax, const double binSize_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingPolar(rMin, rMax, binSize_D1, shift_D1, muMin, muMax, binSize_D2, shift_D2, angularUnits, angularWeight) {}
  
      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Pair2D_comovingPolar_extra () = default;

      ///@}

    };

  
    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingPolar_linlin_extra Pair2D_extra.h
     *  "Headers/Pair2D_extra.h"
     *
     *  @brief The class Pair2D_comovingPolar_linlin_extra
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingPolar_linlin_extra </EM>.
     */
    class Pair2D_comovingPolar_linlin_extra : public virtual Pair2D_comovingPolar_extra, public virtual Pair2D_comovingPolar_linlin {
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Pair2D_comovingPolar_linlin_extra
       */
      Pair2D_comovingPolar_linlin_extra ()
	{
	  m_pairType = PairType::_comovingPolar_linlin_;
	  m_pairInfo = PairInfo::_extra_;
	} 

      /**
       *  @brief constructor   
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins_D1 number of bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param nbins_D2 number of bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar_linlin_extra
       */
      Pair2D_comovingPolar_linlin_extra (const double rMin, const double rMax, const int nbins_D1, const double shift_D1, const double muMin, const double muMax, const int nbins_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(1., nbins_D1, shift_D1, 1., nbins_D2, shift_D2, angularUnits, angularWeight), Pair2D_comovingPolar(rMin, rMax, nbins_D1, shift_D1, muMin, muMax, nbins_D2, shift_D2, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comovingPolar_linlin_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_nbins();
	  m_PP2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_PP2D_weighted.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	}
  
      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize_D1 size of the bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param binSize_D2 size of the bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar_linlin_extra
       */
      Pair2D_comovingPolar_linlin_extra (const double rMin, const double rMax, const double binSize_D1, const double shift_D1, const double muMin, const double muMax, const double binSize_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(binSize_D1, 50, shift_D1, binSize_D2, 50, shift_D2, angularUnits, angularWeight), Pair2D_comovingPolar(rMin, rMax, binSize_D1, shift_D1, muMin, muMax, binSize_D2, shift_D2, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comovingPolar_linlin_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_binSize();
	  m_PP2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_PP2D_weighted.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	} 
  
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Pair2D_comovingPolar_linlin_extra () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;
  
      ///@}

    };

  
    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingPolar_loglin_extra Pair2D_extra.h
     *  "Headers/Pair2D_extra.h"
     *
     *  @brief The class Pair2D_comovingPolar_loglin_extra
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingPolar_loglin_extra </EM>.
     */
    class Pair2D_comovingPolar_loglin_extra : public virtual Pair2D_comovingPolar_extra, public virtual Pair2D_comovingPolar_loglin {
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Pair2D_comovingPolar_loglin_extra
       */
      Pair2D_comovingPolar_loglin_extra ()
	{
	  m_pairType = PairType::_comovingPolar_loglin_;
	  m_pairInfo = PairInfo::_extra_;
	} 

      /**
       *  @brief constructor   
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins_D1 number of bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param nbins_D2 number of bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar_loglin_extra
       */
      Pair2D_comovingPolar_loglin_extra (const double rMin, const double rMax, const int nbins_D1, const double shift_D1, const double muMin, const double muMax, const int nbins_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(1., nbins_D1, shift_D1, 1., nbins_D2, shift_D2, angularUnits, angularWeight), Pair2D_comovingPolar(rMin, rMax, nbins_D1, shift_D1, muMin, muMax, nbins_D2, shift_D2, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comovingPolar_loglin_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_nbins();
	  m_PP2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_PP2D_weighted.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	}  
  
      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize_D1 size of the bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param binSize_D2 size of the bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar_loglin_extra
       */
      Pair2D_comovingPolar_loglin_extra (const double rMin, const double rMax, const double binSize_D1, const double shift_D1, const double muMin, const double muMax, const double binSize_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(binSize_D1, 50, shift_D1, binSize_D2, 50, shift_D2, angularUnits, angularWeight), Pair2D_comovingPolar(rMin, rMax, binSize_D1, shift_D1, muMin, muMax, binSize_D2, shift_D2, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comovingPolar_loglin_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_binSize();
	  m_PP2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_PP2D_weighted.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	}
      
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Pair2D_comovingPolar_loglin_extra () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;
  
      ///@}

    };

  
    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingPolar_linlog_extra Pair2D_extra.h
     *  "Headers/Pair2D_extra.h"
     *
     *  @brief The class Pair2D_comovingPolar_linlog_extra
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingPolar_linlog_extra </EM>.
     */
    class Pair2D_comovingPolar_linlog_extra : public virtual Pair2D_comovingPolar_extra, public virtual Pair2D_comovingPolar_linlog {
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Pair2D_comovingPolar_linlog_extra
       */
      Pair2D_comovingPolar_linlog_extra ()
	{
	  m_pairType = PairType::_comovingPolar_linlog_;
	  m_pairInfo = PairInfo::_extra_;
	} 

      /**
       *  @brief constructor   
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins_D1 number of bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param nbins_D2 number of bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar_linlog_extra
       */
      Pair2D_comovingPolar_linlog_extra (const double rMin, const double rMax, const int nbins_D1, const double shift_D1, const double muMin, const double muMax, const int nbins_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(1., nbins_D1, shift_D1, 1., nbins_D2, shift_D2, angularUnits, angularWeight), Pair2D_comovingPolar(rMin, rMax, nbins_D1, shift_D1, muMin, muMax, nbins_D2, shift_D2, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comovingPolar_linlog_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_nbins();
	  m_PP2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_PP2D_weighted.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	}
  
      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize_D1 size of the bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param binSize_D2 size of the bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar_linlog_extra
       */
      Pair2D_comovingPolar_linlog_extra (const double rMin, const double rMax, const double binSize_D1, const double shift_D1, const double muMin, const double muMax, const double binSize_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(binSize_D1, 50, shift_D1, binSize_D2, 50, shift_D2, angularUnits, angularWeight), Pair2D_comovingPolar(rMin, rMax, binSize_D1, shift_D1, muMin, muMax, binSize_D2, shift_D2, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comovingPolar_linlog_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_binSize();
	  m_PP2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_PP2D_weighted.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	} 
  
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Pair2D_comovingPolar_linlog_extra () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;
  
      ///@}

    };

    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingPolar_loglog_extra Pair2D_extra.h
     *  "Headers/Pair2D_extra.h"
     *
     *  @brief The class Pair2D_comovingPolar_loglog_extra
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingPolar_loglog_extra </EM>.
     */
    class Pair2D_comovingPolar_loglog_extra : public virtual Pair2D_comovingPolar_extra, public virtual Pair2D_comovingPolar_loglog {
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Pair2D_comovingPolar_loglog_extra
       */
      Pair2D_comovingPolar_loglog_extra ()
	{
	  m_pairType = PairType::_comovingPolar_loglog_;
	  m_pairInfo = PairInfo::_extra_;
	} 

      /**
       *  @brief constructor   
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins_D1 number of bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param nbins_D2 number of bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar_loglog_extra
       */
      Pair2D_comovingPolar_loglog_extra (const double rMin, const double rMax, const int nbins_D1, const double shift_D1, const double muMin, const double muMax, const int nbins_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(1., nbins_D1, shift_D1, 1., nbins_D2, shift_D2, angularUnits, angularWeight), Pair2D_comovingPolar(rMin, rMax, nbins_D1, shift_D1, muMin, muMax, nbins_D2, shift_D2, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comovingPolar_loglog_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_nbins();
	  m_PP2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_PP2D_weighted.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	}
  
      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize_D1 size of the bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param binSize_D2 size of the bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar_loglog_extra
       */
      Pair2D_comovingPolar_loglog_extra (const double rMin, const double rMax, const double binSize_D1, const double shift_D1, const double muMin, const double muMax, const double binSize_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair2D(binSize_D1, 50, shift_D1, binSize_D2, 50, shift_D2, angularUnits, angularWeight), Pair2D_comovingPolar(rMin, rMax, binSize_D1, shift_D1, muMin, muMax, binSize_D2, shift_D2, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comovingPolar_loglog_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_binSize();
	  m_PP2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_PP2D_weighted.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D1_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_scale_D2_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_mean.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_S.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	  m_z_sigma.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.));
	} 
  
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Pair2D_comovingPolar_loglog_extra () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;
  
      ///@}

    };  
  }
}

#endif
