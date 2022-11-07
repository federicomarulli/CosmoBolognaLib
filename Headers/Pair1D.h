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
 *  @file Headers/Pair1D.h
 *
 *  @brief The classes Pair1D*
 *
 *  This file defines the interface of all the classes Pair1D*, used to
 *  handle 1D pairs of objects of any kind
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unibo.it
 */

#ifndef __PAIR1D__
#define __PAIR1D__


#include "Pair.h"


// ===================================================================================================


namespace cbl {
  
  namespace pairs {

    /**
     *  @class Pair1D Pair1D.h "Headers/Pair1D.h"
     *
     *  @brief The class Pair1D
     *
     *  This class is used to handle objects of type <EM> Pair1D
     *  </EM>.
     */
    class Pair1D : public virtual Pair {
      
    protected:
      
      /**
       *  @name Member functions used to set the binning parameters (customized in all the derived classes) 
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       */
      virtual void m_set_parameters_nbins () = 0;
  
      /**
       *  @brief set the binning parameters given the bin size
       */
      virtual void m_set_parameters_binSize () = 0;
  
      ///@}


    protected:
      
      /// the binned scales
      std::vector<double> m_scale;
      
      /// the number of binned pairs 
      std::vector<double> m_PP1D;
      
      /// the number of binned weighted pairs 
      std::vector<double> m_PP1D_weighted;

      /**
       *  @name binning parameters
       */
      ///@{
  
      /// the inverse of the bin size
      double m_binSize_inv;
  
      /// number of bins
      int m_nbins;
  
      /// radial shift used to centre the output bins 
      double m_shift;
  
      ///@}

  
    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D () : m_binSize_inv(1.), m_nbins(50), m_shift(0.5)
	{ m_pairDim = Dim::_1D_; m_angularUnits = CoordinateUnits::_radians_; m_angularWeight = nullptr; }

      /**
       *  @brief constructor
       *  @param binSize bin size
       *  @param nbins number of bins
       *  @param shift radial shift used to centre the output bins 
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D (const double binSize, const int nbins, const double shift, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight=nullptr)
	: m_binSize_inv(1./binSize), m_nbins(nbins), m_shift(shift)
      {
	m_pairDim = Dim::_1D_;
	m_PP1D.resize(m_nbins+1, 0.);
	m_PP1D_weighted.resize(m_nbins+1, 0.);
	m_angularUnits = angularUnits;
	m_angularWeight = angularWeight;
      }
    
      /**
       *  @brief default destructor
       *  
       */
      virtual ~Pair1D () = default;

      ///@}

  
      /**
       *  @name Member functions used to get the protected members
       */
      ///@{
      
      /**
       *  @brief get the protected member Pair1D::m_scale[i]
       *  @param i the bin index
       *  @return the i-th binned scale
       */
      double scale (const int i) const override { return m_scale[i]; }

      /**
       *  @brief get the protected member Pair1D::m_scale
       *  @return the vector containing the binned scales
       */
      std::vector<double> scale () const override { return m_scale; }
      
      /**
       *  @brief get the protected member Pair1D::m_PP1D[i]
       *  @param i the bin index
       *  @return the number of pairs in the i-th bin
       */
      double PP1D (const int i) const override { return m_PP1D[i]; }

      /**
       *  @brief get the protected member Pair1D::m_PP1D
       *  @return the vector containing the binned number of pairs
       */
      std::vector<double> PP1D () const override { return m_PP1D; }

      /**
       *  @brief get the protected member Pair1D::m_PP1D_weighted[i]
       *  @param i the bin index
       *  @return the number of weighted pairs in the i-th bin
       */
      double PP1D_weighted (const int i) const override { return m_PP1D_weighted[i]; }

      /**
       *  @brief get the protected member Pair1D::m_PP1D
       *  @return the vector containing the binned number of weighted
       *  pairs
       */
      std::vector<double> PP1D_weighted () const override { return m_PP1D_weighted; }
    
      /**
       *  @brief get the protected member Pair1D::m_binSize_inv
       *  @return the inverse of the bin size
       */
      double binSize_inv () const override { return m_binSize_inv; }

      /**
       *  @brief get the protected member Pair1D::m_nbins
       *  @return the number of bins
       */
      int nbins () const override { return m_nbins; }
    
      /**
       *  @brief get the protected member Pair1D::m_shift
       *  @return the radial shift used to centre the output bins
       */
      double shift () const override { return m_shift; }

      ///@}

  
      /**
       *  @name Member functions used to set the protected members
       */
      ///@{

      /**
       * @brief reset the pair counts
       *
       * @details set to 0 the m_PP1D and 
       * m_PP1D_weighted vector elements
       *
       * 
       */
      void reset () override;
      
      /**
       *  @brief set the protected member Pair1D::m_scale[i]
       *  @param i the bin index
       *  @param ss the scale value
       *  
       */
      void set_scale (const int i, const double ss) { checkDim(m_scale, i, "m_scale", false); m_scale[i] = ss; }
      
      /**
       *  @brief set the protected member Pair1D::m_PP1D[i]
       *  @param i the bin index
       *  @param pp the number of pairs in the bin
       *  
       */
      void set_PP1D (const int i, const double pp) { checkDim(m_PP1D, i, "m_PP1D", false); m_PP1D[i] = pp; }

      /**
       *  @brief set the protected member Pair1D::m_PP1D_weighted[i]
       *  @param i the bin index
       *  @param pp the number of weighted pairs in the bin
       *  
       */
      void set_PP1D_weighted (const int i, const double pp) { checkDim(m_PP1D_weighted, i, "m_PP1D_weighted", false); m_PP1D_weighted[i] = pp; }
      
      /**
       *  @brief set the protected member by adding new data
       *  @param i the bin index
       *  @param data the number of pairs in the bin
       *  
       */
      void add_data1D (const int i, const std::vector<double> data) override;

      /**
       *  @brief set the protected member by adding new data
       *  @param i the bin index
       *  @param pair pair pointer to an object of class Pair
       *  @param ww a multiplicative factor used for bootstrap
       *  
       */
      void add_data1D (const int i, const std::shared_ptr<pairs::Pair> pair, const double ww=1.) override;
      
      ///@}

    
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
    
      /**
       *  @brief sum the number of binned pairs
       *  @param pair an object of class Pair
       *  @param ww the weight
       */
      virtual void Sum (const std::shared_ptr<Pair> pair, const double ww=1) override;

      ///@}
    
    };


    // ============================================================================================
    // ============================================================================================

  
    /**
     *  @class Pair1D_angular Pair1D.h "Headers/Pair1D.h"
     *
     *  @brief The class Pair1D_angular
     *
     *  This class is used to handle objects of type <EM> Pair1D_angular
     *  </EM>.
     */
    class Pair1D_angular : public virtual Pair1D {

    protected:
      
      /**
       *  @name Member functions used to set the binning parameters (customized in all the derived classes) 
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       */
      virtual void m_set_parameters_nbins () = 0;
  
      /**
       *  @brief set the binning parameters given the bin size
       */
      virtual void m_set_parameters_binSize () = 0;
  
      ///@}

      
    protected:
  
      /**
       *  @name binning parameters
       */
      ///@{

      /// the minimum value of the angle &theta; used to count the number of pairs
      double m_thetaMin; 
  
      /// the maximum value of the angle &theta; used to count the number of pairs
      double m_thetaMax;
  
      ///@}

  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_angular () : m_thetaMin(0.1), m_thetaMax(50.) {} 

      /**
       *  @brief constructor
       *  @param thetaMin minimum value of the angle &theta; used to count
       *  the pairs
       *  @param thetaMax maximum value of the angle &theta; used to count
       *  the pairs
       *  @param nbins number of bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_angular (const double thetaMin, const double thetaMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(1., nbins, shift, angularUnits, angularWeight), m_thetaMin(thetaMin), m_thetaMax(thetaMax) {}
  
      /**
       *  @brief constructor
       *
       *  @param thetaMin minimum value of the angle &theta; used to
       *  count the pairs
       *
       *  @param thetaMax maximum value of the angle &theta; used to
       *  count the pairs
       *
       *  @param binSize size of the bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_angular (const double thetaMin, const double thetaMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(binSize, 50, shift, angularUnits, angularWeight), m_thetaMin(thetaMin), m_thetaMax(thetaMax) {}
  
      /**
       *  @brief default destructor
       *  
       */
      virtual ~Pair1D_angular () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to get the protected members
       */
      ///@{

      /**
       *  @brief get the protected member Pair1D_angular::m_thetaMin
       *  @return the minimum value of the angle &theta; used to count
       *  the number of pairs
       */
      double sMin () const override { return m_thetaMin; }

      /**
       *  @brief get the protected member Pair1D_angular::m_thetaMax
       *  @return the maximum value of the angle &theta; used to count
       *  the number of pairs
       */    
      double sMax () const override { return m_thetaMax; }

      ///@}
    
    };

  
    // ============================================================================================
    // ============================================================================================

    /**
     *  @class Pair1D_angular_lin Pair1D.h "Headers/Pair1D.h"
     *
     *  @brief The class Pair1D_angular_lin
     *
     *  This class is used to handle objects of type <EM> Pair1D_angular_lin
     *  </EM>.
     */
    class Pair1D_angular_lin : public virtual Pair1D_angular {

    protected:
    
      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  
       */
      void m_set_parameters_nbins () override;
    
      /**
       *  @brief set the binning parameters given the bin size
       *  
       */
      void m_set_parameters_binSize () override;

      ///@}
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_angular_lin ()
	{
	  m_pairType = PairType::_angular_lin_;
	  m_pairInfo = PairInfo::_standard_;
	} 

      /**
       *  @brief constructor
       *  @param thetaMin minimum value of the angle &theta; used to count
       *  the pairs
       *  @param thetaMax maximum value of the angle &theta; used to count
       *  the pairs
       *  @param nbins number of bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_angular_lin (const double thetaMin, const double thetaMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(1., nbins, shift, angularUnits, angularWeight), Pair1D_angular(thetaMin, thetaMax, nbins, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_angular_lin_;
	  m_pairInfo = PairInfo::_standard_; 
	  m_set_parameters_nbins();
	  m_PP1D.resize(m_nbins+1, 0.);
	  m_PP1D_weighted.resize(m_nbins+1, 0.);
	}
  
      /**
       *  @brief constructor
       *  @param thetaMin minimum value of the angle &theta; used to count
       *  the pairs
       *  @param thetaMax maximum value of the angle &theta; used to count
       *  the pairs
       *  @param binSize size of the bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_angular_lin (const double thetaMin, const double thetaMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(binSize, 50, shift, angularUnits, angularWeight), Pair1D_angular(thetaMin, thetaMax, binSize, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_angular_lin_;
	  m_pairInfo = PairInfo::_standard_; 
	  m_set_parameters_binSize();
	  m_PP1D.resize(m_nbins+1, 0.);
	  m_PP1D_weighted.resize(m_nbins+1, 0.);
	} 
  
      /**
       *  @brief default destructor
       *  
       */
      ~Pair1D_angular_lin () = default;
      
      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief get the pair index and weight
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @param kk index of the pairs
       *  @param wkk weight of the pair
       *  
       */
      void get (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2, int &kk, double &wkk) override;

      /**
       *  @brief set the pair vector 
       *  @param kk index of the pairs
       *  @param wkk weight of the pair
       *  @param weight the weght of the region
       *  
       */
      void set (const int kk, const double wkk, const double weight=1) override;

      /**
       *  @brief estimate the distance between two objects and update
       *  the pair vector accordingly
       *
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;

      ///@}
    
    };

  
    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair1D_angular_log Pair1D.h "Headers/Pair1D.h"
     *
     *  @brief The class Pair1D_angular_log
     *
     *  This class is used to handle objects of type <EM> Pair1D_angular_log
     *  </EM>.
     */
    class Pair1D_angular_log : public virtual Pair1D_angular {

    protected:
    
      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  
       */
      void m_set_parameters_nbins () override;
    
      /**
       *  @brief set the binning parameters given the bin size
       *  
       */
      void m_set_parameters_binSize () override;

      ///@}
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_angular_log ()
	{
	  m_pairType = PairType::_angular_log_;
	  m_pairInfo = PairInfo::_standard_;
	} 

      /**
       *  @brief constructor
       *  @param thetaMin minimum value of the angle &theta; used to count
       *  the pairs
       *  @param thetaMax maximum value of the angle &theta; used to count
       *  the pairs
       *  @param nbins number of bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_angular_log (const double thetaMin, const double thetaMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(1., nbins, shift, angularUnits, angularWeight), Pair1D_angular(thetaMin, thetaMax, nbins, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_angular_log_;
	  m_pairInfo = PairInfo::_standard_; 
	  m_set_parameters_nbins();
	  m_PP1D.resize(m_nbins+1, 0.);
	  m_PP1D_weighted.resize(m_nbins+1, 0.);
	}
  
      /**
       *  @brief constructor
       *  @param thetaMin minimum value of the angle &theta; used to count
       *  the pairs
       *  @param thetaMax maximum value of the angle &theta; used to count
       *  the pairs
       *  @param binSize size of the bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_angular_log (const double thetaMin, const double thetaMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(binSize, 50, shift, angularUnits, angularWeight), Pair1D_angular(thetaMin, thetaMax, binSize, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_angular_log_;
	  m_pairInfo = PairInfo::_standard_; 
	  m_set_parameters_binSize();
	  m_PP1D.resize(m_nbins+1, 0.);
	  m_PP1D_weighted.resize(m_nbins+1, 0.);
	} 
  
      /**
       *  @brief default destructor
       *  
       */
      ~Pair1D_angular_log () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
    
      /**
       *  @brief get the pair index and weight
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @param kk index of the pairs
       *  @param wkk weight of the pair
       *  
       */
      void get (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2, int &kk, double &wkk) override;

      /**
       *  @brief set the pair vector 
       *  @param kk index of the pairs
       *  @param wkk weight of the pair
       *  @param weight the weght of the region
       *  
       */
      void set (const int kk, const double wkk, const double weight=1) override;

      /**
       *  @brief estimate the distance between two objects and update
       *  the pair vector accordingly
       *
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;

      ///@}
    
    };

  
    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair1D_comoving Pair1D.h "Headers/Pair1D.h"
     *
     *  @brief The class Pair1D_comoving
     *
     *  This class is used to handle objects of type <EM> Pair1D_comoving
     *  </EM>.
     */
    class Pair1D_comoving : public virtual Pair1D {

    protected:
      
      /**
       *  @name Member functions used to set the binning parameters (customized in all the derived classes) 
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       */
      virtual void m_set_parameters_nbins () = 0;
  
      /**
       *  @brief set the binning parameters given the bin size
       */
      virtual void m_set_parameters_binSize () = 0;
  
      ///@}
      
      
    protected:
  
      /**
       *  @name binning parameters
       */
      ///@{

      /// minimum separation used to count the pairs
      double m_rMin;
  
      /// maximum separation used to count the pairs
      double m_rMax;

      ///@}

  
    public:
  
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_comoving () : m_rMin(0.1), m_rMax(50.) {} 
    
      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins number of bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_comoving (const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(1., nbins, shift, angularUnits, angularWeight), m_rMin(rMin), m_rMax(rMax) {} 
  
      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize size of the bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_comoving (const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(binSize, 50, shift, angularUnits, angularWeight), m_rMin(rMin), m_rMax(rMax) {} 
  
      /**
       *  @brief default destructor
       *  
       */
      virtual ~Pair1D_comoving () = default;

      ///@}

  
      /**
       *  @name Member functions used to get the protected parameters
       */
      ///@{

      /**
       *  @brief get the protected member Pair1D::m_rMin
       *  @return the minimum separation used to count the pairs
       */
      double sMin () const override { return m_rMin; }

      /**
       *  @brief get the protected member Pair1D::m_rMax
       *  @return the maximum separation used to count the pairs
       */
      double sMax () const override { return m_rMax; }

      ///@}
    
    };

  
    // ============================================================================================
    // ============================================================================================

  
    /**
     *  @class Pair1D_comoving_lin Pair1D.h "Headers/Pair1D.h"
     *
     *  @brief The class Pair1D_comoving_lin
     *
     *  This class is used to handle objects of type <EM> Pair1D_comoving_lin
     *  </EM>.
     */
    class Pair1D_comoving_lin : public virtual Pair1D_comoving {

    protected:
    
      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  
       */
      void m_set_parameters_nbins () override;
    
      /**
       *  @brief set the binning parameters given the bin size
       *  
       */
      void m_set_parameters_binSize () override;

      ///@}
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_comoving_lin ()
	{
	  m_pairType = PairType::_comoving_lin_;
	  m_pairInfo = PairInfo::_standard_;
	} 

      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins number of bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_comoving_lin (const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(1., nbins, shift, angularUnits, angularWeight), Pair1D_comoving(rMin, rMax, nbins, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comoving_lin_;
	  m_pairInfo = PairInfo::_standard_; 
	  m_set_parameters_nbins();
	  m_PP1D.resize(m_nbins+1, 0.);
	  m_PP1D_weighted.resize(m_nbins+1, 0.);
	}
  
      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize size of the bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_comoving_lin (const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(binSize, 50, shift, angularUnits, angularWeight), Pair1D_comoving(rMin, rMax, binSize, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comoving_lin_;
	  m_pairInfo = PairInfo::_standard_; 
	  m_set_parameters_binSize();
	  m_PP1D.resize(m_nbins+1, 0.);
	  m_PP1D_weighted.resize(m_nbins+1, 0.);
	} 
  
      /**
       *  @brief default destructor
       *  
       */
      ~Pair1D_comoving_lin () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{ 
      
      /**
       *  @brief get the pair index and weight
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @param kk index of the pairs
       *  @param wkk weight of the pair
       *  
       */
      void get (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2, int &kk, double &wkk) override;

      /**
       *  @brief set the pair vector 
       *  @param kk index of the pairs
       *  @param wkk weight of the pair
       *  @param weight the weght of the region
       *  
       */
      void set (const int kk, const double wkk, const double weight=1) override;
  
      /**
       *  @brief estimate the distance between two objects and update
       *  the pair vector accordingly
       *
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;
  
      ///@}
    
    };

  
    // ============================================================================================
    // ============================================================================================

  
    /**
     *  @class Pair1D_comoving_log Pair1D.h "Headers/Pair1D.h"
     *
     *  @brief The class Pair1D_comoving_log
     *
     *  This class is used to handle objects of type <EM> Pair1D_comoving_log
     *  </EM>.
     */
    class Pair1D_comoving_log : public virtual Pair1D_comoving {

    protected:
      
      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  
       */
      void m_set_parameters_nbins () override;
    
      /**
       *  @brief set the binning parameters given the bin size
       *  
       */
      void m_set_parameters_binSize () override;

      ///@}
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_comoving_log ()
	{
	  m_pairType = PairType::_comoving_log_;
	  m_pairInfo = PairInfo::_standard_; 
	}
      
      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins number of bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_comoving_log (const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(1., nbins, shift, angularUnits, angularWeight), Pair1D_comoving(rMin, rMax, nbins, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comoving_log_;
	  m_pairInfo = PairInfo::_standard_; 
	  m_set_parameters_nbins();
	  m_PP1D.resize(m_nbins+1, 0.);
	  m_PP1D_weighted.resize(m_nbins+1, 0.);
	}
  
      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize size of the bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_comoving_log (const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(binSize, 50, shift, angularUnits, angularWeight), Pair1D_comoving(rMin, rMax, binSize, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comoving_log_;
	  m_pairInfo = PairInfo::_standard_; 
	  m_set_parameters_binSize();
	  m_PP1D.resize(m_nbins+1, 0.);
	  m_PP1D_weighted.resize(m_nbins+1, 0.);
	} 
  
      /**
       *  @brief default destructor
       *  
       */
      ~Pair1D_comoving_log () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
        
      /**
       *  @brief get the pair index and weight
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @param kk index of the pairs
       *  @param wkk weight of the pair
       *  
       */
      void get (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2, int &kk, double &wkk) override;

      /**
       *  @brief set the pair vector 
       *  @param kk index of the pairs
       *  @param wkk weight of the pair
       *  @param weight the weght of the region
       *  
       */
      void set (const int kk, const double wkk, const double weight=1) override;

      /**
       *  @brief estimate the distance between two objects and update
       *  the pair vector accordingly
       *
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;

      ///@}
    
    };

  
    // ============================================================================================
    // ============================================================================================

     /**
     *  @class Pair1D_comoving_multipoles Pair1D.h "Headers/Pair1D.h"
     *
     *  @brief The class Pair1D_comoving_multipoles
     *
     *  This class is used to handle objects of type <EM> Pair1D_comoving_multipoles
     *  </EM>.
     */
    class Pair1D_comoving_multipoles : public virtual Pair1D {

    protected:
      
      /**
       *  @name Member functions used to set the binning parameters (customized in all the derived classes) 
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       */
      virtual void m_set_parameters_nbins () = 0;
  
      /**
       *  @brief set the binning parameters given the bin size
       */
      virtual void m_set_parameters_binSize () = 0;
  
      ///@}
      
      
    protected:
  
      /**
       *  @name binning parameters
       */
      ///@{

      /// minimum separation used to count the pairs
      double m_rMin;
  
      /// maximum separation used to count the pairs
      double m_rMax;

      ///@}

  
    public:
  
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_comoving_multipoles () : m_rMin(0.1), m_rMax(50.) {} 
    
      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins number of bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_comoving_multipoles (const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(1., nbins, shift, angularUnits, angularWeight), m_rMin(rMin), m_rMax(rMax) 
      {
	m_PP1D.resize(3*(m_nbins+1), 0.);
	m_PP1D_weighted.resize(3*(m_nbins+1), 0.);
      }

      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize size of the bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_comoving_multipoles (const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(binSize, 50, shift, angularUnits, angularWeight), m_rMin(rMin), m_rMax(rMax)
      {
	m_PP1D.resize(3*(m_nbins+1), 0.);
	m_PP1D_weighted.resize(3*(m_nbins+1), 0.);
      }

      /**
       *  @brief default destructor
       *  
       */
      virtual ~Pair1D_comoving_multipoles () = default;

      ///@}

  
      /**
       *  @name Member functions used to get the protected parameters
       */
      ///@{

      /**
       *  @brief get the protected member Pair1D::m_rMin
       *  @return the minimum separation used to count the pairs
       */
      double sMin () const override { return m_rMin; }

      /**
       *  @brief get the protected member Pair1D::m_rMax
       *  @return the maximum separation used to count the pairs
       */
      double sMax () const override { return m_rMax; }

      ///@}
    
    };

  
    // ============================================================================================
    // ============================================================================================
 
    /**
     *  @class Pair1D_comoving_multipoles_lin Pair1D.h "Headers/Pair1D.h"
     *
     *  @brief The class Pair1D_comoving_multipoles_lin
     *
     *  This class is used to handle objects of type <EM> Pair1D_comoving_multipoles_lin
     *  </EM>.
     */
    class Pair1D_comoving_multipoles_lin : public virtual Pair1D_comoving_multipoles {

    protected:
    
      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  
       */
      void m_set_parameters_nbins () override;
    
      /**
       *  @brief set the binning parameters given the bin size
       *  
       */
      void m_set_parameters_binSize () override;

      ///@}
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_comoving_multipoles_lin ()
	{
	  m_pairType = PairType::_comoving_lin_;
	  m_pairInfo = PairInfo::_standard_;
	} 

      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins number of bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_comoving_multipoles_lin (const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(1., nbins, shift, angularUnits, angularWeight), Pair1D_comoving_multipoles(rMin, rMax, nbins, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comoving_multipoles_lin_;
	  m_pairInfo = PairInfo::_standard_; 
	  m_set_parameters_nbins();
	  m_PP1D.resize(3*(m_nbins+1), 0.);
	  m_PP1D_weighted.resize(3*(m_nbins+1), 0.);
	}
  
      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize size of the bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_comoving_multipoles_lin (const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(binSize, 50, shift, angularUnits, angularWeight), Pair1D_comoving_multipoles(rMin, rMax, binSize, shift, angularUnits, angularWeight) 
	{
	  m_pairType = PairType::_comoving_multipoles_lin_;
	  m_pairInfo = PairInfo::_standard_; 
	  m_set_parameters_binSize();
	  m_PP1D.resize(3*(m_nbins+1), 0.);
	  m_PP1D_weighted.resize(3*(m_nbins+1), 0.);
	} 
  
      /**
       *  @brief default destructor
       *  
       */
      ~Pair1D_comoving_multipoles_lin () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{ 
      
      /**
       *  @brief get the pair index and weight
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @param kk index of the pairs
       *  @param cosmu cosine of the angle between objects
       *  @param wkk weight of the pair
       *  
       */
      void get (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2, int &kk, double &cosmu, double &wkk) override;

      /**
       *  @brief set the pair vector 
       *  @param wkk weight of the pair
       *  @param kk index of the pairs
       *  @param cosmu cosine of the angle between objects
       *  @param weight the weght of the region
       *  
       */
      void set (const double cosmu, const int kk, const double wkk, const double weight=1) override;
  
      /**
       *  @brief estimate the distance between two objects and update
       *  the pair vector accordingly
       *
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;
  
      ///@}
    
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{

      /**
       *  @brief sum the number of binned pairs
       *  @param pair an object of class Pair
       *  @param ww the weight
       */
      virtual void Sum (const std::shared_ptr<Pair> pair, const double ww=1) override;

      ///@}
    };

  
    // ============================================================================================
    // ============================================================================================

  
    /**
     *  @class Pair1D_comoving_multipoles_log Pair1D.h "Headers/Pair1D.h"
     *
     *  @brief The class Pair1D_comoving_multipoles_log
     *
     *  This class is used to handle objects of type <EM> Pair1D_comoving_multipoles_log
     *  </EM>.
     */
    class Pair1D_comoving_multipoles_log : public virtual Pair1D_comoving_multipoles {

    protected:
      
      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  
       */
      void m_set_parameters_nbins () override;
    
      /**
       *  @brief set the binning parameters given the bin size
       *  
       */
      void m_set_parameters_binSize () override;

      ///@}
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_comoving_multipoles_log ()
	{
	  m_pairType = PairType::_comoving_log_;
	  m_pairInfo = PairInfo::_standard_; 
	}

      
      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins number of bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_comoving_multipoles_log (const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(1., nbins, shift, angularUnits, angularWeight), Pair1D_comoving_multipoles (rMin, rMax, nbins, shift, angularUnits, angularWeight) 
      {
	m_pairType = PairType::_comoving_multipoles_log_;
	m_pairInfo = PairInfo::_standard_; 
	m_set_parameters_nbins();
	m_PP1D.resize(3*(m_nbins+1), 0.);
	m_PP1D_weighted.resize(3*(m_nbins+1), 0.);
      } 

      /**
       *  @brief constructor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize size of the bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_comoving_multipoles_log (const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(binSize, 50, shift, angularUnits, angularWeight), Pair1D_comoving_multipoles(rMin, rMax, binSize, shift, angularUnits, angularWeight) 
      {
	m_pairType = PairType::_comoving_multipoles_log_;
	m_pairInfo = PairInfo::_standard_; 
	m_set_parameters_binSize();
	m_PP1D.resize(3*(m_nbins+1), 0.);
	m_PP1D_weighted.resize(3*(m_nbins+1), 0.);
      } 
  
      /**
       *  @brief default destructor
       *  
       */
      ~Pair1D_comoving_multipoles_log () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs, 
       *  for direct multipoles computations
       */
      ///@{
        
      /**
       *  @brief get the pair index and weight
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @param kk index of the pairs
       *  @param cosmu cosine of the angle between objects
       *  @param wkk weight of the pair
       *  
       */
      void get (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2, int &kk, double &cosmu, double &wkk) override;

      /**
       *  @brief set the pair vector 
       *  @param cosmu cosine of the angle between objects
       *  @param kk index of the pairs
       *  @param wkk weight of the pair
       *  @param weight the weght of the region
       *  
       */
      void set (const double cosmu, const int kk, const double wkk, const double weight=1) override;

      /**
       *  @brief estimate the distance between two objects and update
       *  the pair vector accordingly
       *
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;

      ///@}

      /**
       *  @name Member functions used to handle pairs
       */
      ///@{

      /**
       *  @brief sum the number of binned pairs
       *  @param pair an object of class Pair
       *  @param ww the weight
       */
      virtual void Sum (const std::shared_ptr<Pair> pair, const double ww=1) override;

      ///@}

    };

  }
}

#endif
