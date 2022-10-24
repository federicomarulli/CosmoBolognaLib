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
 *  @file Headers/Pair1D_extra.h
 *
 *  @brief The classes Pair1D_extra*
 *
 *  This file defines the interface of all the classes Pair1D_extra*,
 *  used to handle 1D pairs of objects of any kind with extra
 *  information, such as the mean redshift and separation
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unibo.it
 */

#ifndef __PAIR1D_EXTRA__
#define __PAIR1D_EXTRA__


#include "Pair1D.h"


// ===================================================================================================


namespace cbl {
  
  namespace pairs {

     /**
     *  @class Pair1D_extra Pair1D_extra.h "Headers/Pair1D_extra.h"
     *
     *  @brief The class Pair1D_extra
     *
     *  This class is used to handle objects of type <EM> Pair1D_extra
     *  </EM>.
     */
    class Pair1D_extra : public virtual Pair1D {
      
    protected:
      
      /// the mean scales in each bin
      std::vector<double> m_scale_mean;

      /// the square of the standard deviations of the scale distributions in each bin, multiplied by the total weight
      std::vector<double> m_scale_S;
      
      /// the standard deviations of the scale distributions in each bin
      std::vector<double> m_scale_sigma;

      /// the mean redshift in each bin
      std::vector<double> m_z_mean;

      /// the square of the standard deviations of the redshift distributions in each bin, multiplied by the total weight
      std::vector<double> m_z_S;
      
      /// the standard deviations of the redshift distributions in each bin
      std::vector<double> m_z_sigma;

      
      /**
       *  @name parameters used for internal computations
       */
      ///@{
      
      /// factor used to combine two standard deviations of the scale distribution
      std::vector<double> m_fact_scale;
      
      /// factor used to combine two standard deviations of the redshift distribution
      std::vector<double> m_fact_z;

      ///@}
      
      
    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_extra () = default;

      /**
       *  @brief constructor
       *  @param binSize bin size
       *  @param nbins number of bins
       *  @param shift radial shift used to centre the output bins 
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_extra (const double binSize, const int nbins, const double shift, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(binSize, nbins, shift, angularUnits, angularWeight) {}
    
      /**
       *  @brief default destructor
       *  
       */
      virtual ~Pair1D_extra () = default;

      ///@}

  
      /**
       *  @name Member functions used to get the protected members
       */
      ///@{
      
      /**
       *  @brief get the protected member Pair1D_extra::m_scale_mean[i]
       *  @param i the bin index
       *  @return the mean scale in the i-th bin
       */
      double scale_mean (const int i) const override { return m_scale_mean[i]; }

      /**
       *  @brief get the protected member Pair1D_extra::m_scale_mean
       *  @return the vector containing the mean scales 
       */
      std::vector<double> scale_mean () const override { return m_scale_mean; }

      /**
       *  @brief get the protected member Pair1D_extra::m_scale_S[i]
       *  @param i the bin index
       *  @return the square of the standard deviations of the scale
       *  distribution, multiplied by the total weight, in the i-th
       *  bin
       */
      double scale_S (const int i) const override { return m_scale_S[i]; }

      /**
       *  @brief get the protected member Pair1D_extra::m_scale_S
       *  @return the vector the square of the standard deviations of
       *  the scale distribution, multiplied by the total weight
       */
      std::vector<double> scale_S () const override { return m_scale_S; }

      /**
       *  @brief get the protected member Pair1D_extra::m_scale_sigma[i]
       *  @param i the bin index
       *  @return the standard deviation of the scale distribution in
       *  the i-th bin
       */
      double scale_sigma (const int i) const override { return m_scale_sigma[i]; }

      /**
       *  @brief get the protected member Pair1D_extra::m_scale_sigma
       *  @return the vector containing the standard deviations of the
       *  scale distribution
       */
      std::vector<double> scale_sigma () const override { return m_scale_sigma; }

      /**
       *  @brief get the protected member Pair1D_extra::m_z_mean[i]
       *  @param i the bin index
       *  @return the mean redshift in the i-th bin
       */
      double z_mean (const int i) const override { return m_z_mean[i]; }

      /**
       *  @brief get the protected member Pair1D_extra::m_z_mean
       *  @return the vector containing the mean redshifts 
       */
      std::vector<double> z_mean () const override { return m_z_mean; }

      /**
       *  @brief get the protected member Pair1D_extra::m_z_S[i]
       *  @param i the bin index
       *  @return the square of the standard deviations of the
       *  redshift distribution, multiplied by the total weight, in
       *  the i-th bin
       */
      double z_S (const int i) const override { return m_z_S[i]; }

      /**
       *  @brief get the protected member Pair1D_extra::m_z_S
       *  @return the vector the square of the standard deviations of
       *  the redshift distribution, multiplied by the total weight
       */
      std::vector<double> z_S () const override { return m_z_S; }
      
      /**
       *  @brief get the protected member Pair1D_extra::m_z_sigma[i]
       *  @param i the bin index
       *  @return the standard deviation of the redshift distribution
       *  in the i-th bin
       */
      double z_sigma (const int i) const override { return m_z_sigma[i]; }

      /**
       *  @brief get the protected member Pair1D_extra::m_z_mean
       *  @return the vector containing the standard deviations of the
       *  redshift distribution
       */
      std::vector<double> z_sigma () const override { return m_z_sigma; }
      
      ///@}


      /**
       *  @name Member functions used to set the protected members
       */
      ///@{

      /**
       *  @brief set the protected member Pair1D_extra::m_scale_mean[i]
       *  @param i the bin index 
       *  @param ss the mean scale
       *  
       */
      void set_scale_mean (const int i, const double ss) { checkDim(m_scale_mean, i, "m_scale_mean"); m_scale_mean[i] = ss; }

      /**
       *  @brief set the protected member Pair1D_extra::m_scale_sigma[i]
       *  @param i the bin index 
       *  @param ss the standard deviation of the scale distribution
       *  
       */
      void set_scale_sigma (const int i, const double ss) { checkDim(m_scale_sigma, i, "m_scale_sigma"); m_scale_sigma[i] = ss; }

      /**
       *  @brief set the protected member Pair1D_extra::m_z_mean[i]
       *  @param i the bin index 
       *  @param ss the mean redshift
       *  
       */
      void set_z_mean (const int i, const double ss) { checkDim(m_z_mean, i, "m_z_mean"); m_z_mean[i] = ss; }

      /**
       *  @brief set the protected member Pair1D_extra::m_z_sigma[i]
       *  @param i the bin index 
       *  @param ss the standard deviation of the redshift distribution
       *  
       */
      void set_z_sigma (const int i, const double ss) { checkDim(m_z_sigma, i, "m_z_sigma"); m_z_sigma[i] = ss; }
      
      /**
       *  @brief set the protected members by adding new data
       *  @param i the bin index
       *  @param data vector containing the new data to be added
       *  
       */
      void add_data1D (const int i, const std::vector<double> data) override;
      
      /**
       *  @brief set the protected members by adding new data
       *  @param i the bin index
       *  @param pair pair pointer to an object of class Pair
       *  @param ww a multiplicative factor used for bootstrap
       *  
       */
      void add_data1D (const int i, const std::shared_ptr<pairs::Pair> pair, const double ww=1.) override;
      
      ///@}

      
      /**
       *  @name Member functions used to handle object pairs
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
     *  @class Pair1D_angular_extra Pair1D_extra.h "Headers/Pair1D_extra.h"
     *
     *  @brief The class Pair1D_angular_extra
     *
     *  This class is used to handle objects of type <EM> Pair1D_angular_extra
     *  </EM>.
     */
    class Pair1D_angular_extra : public virtual Pair1D_extra, public virtual Pair1D_angular {
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_angular_extra () = default;

      /**
       *  @brief constructor
       *
       *  @param thetaMin minimum value of the angle &theta; used to
       *  count the pairs
       *
       *  @param thetaMax maximum value of the angle &theta; used to
       *  count
       *
       *  the pairs
       *  @param nbins number of bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_angular_extra (const double thetaMin, const double thetaMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D_angular(thetaMin, thetaMax, nbins, shift, angularUnits, angularWeight) {}
  
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
      Pair1D_angular_extra (const double thetaMin, const double thetaMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D_angular(thetaMin, thetaMax, binSize, shift, angularUnits, angularWeight) {}
  
      /**
       *  @brief default destructor
       *  
       */
      virtual ~Pair1D_angular_extra () = default;

      ///@}
    
    };
    
  
    // ============================================================================================
    // ============================================================================================

    
    /**
     *  @class Pair1D_angular_lin_extra Pair1D_extra.h "Headers/Pair1D_extra.h"
     *
     *  @brief The class Pair1D_angular_lin_extra
     *
     *  This class is used to handle objects of type <EM> Pair1D_angular_lin_extra
     *  </EM>.
     */
    class Pair1D_angular_lin_extra : public virtual Pair1D_angular_extra, public virtual Pair1D_angular_lin {
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_angular_lin_extra ()
	{
	  m_pairType = PairType::_angular_lin_;
	  m_pairInfo = PairInfo::_extra_;
	} 

      /**
       *  @brief constructor
       *  @param thetaMin minimum value of the angle &theta; used to
       *  count the pairs
       *  @param thetaMax maximum value of the angle &theta; used to
       *  count the pairs
       *  @param nbins number of bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       */
      Pair1D_angular_lin_extra (const double thetaMin, const double thetaMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(1., nbins, shift, angularUnits, angularWeight), Pair1D_angular(thetaMin, thetaMax, nbins, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_angular_lin_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_nbins();
	  m_scale_mean.resize(m_nbins+1, 0.);
	  m_scale_S.resize(m_nbins+1, 0.);
	  m_scale_sigma.resize(m_nbins+1, 0.);
	  m_z_mean.resize(m_nbins+1, 0.);
	  m_z_S.resize(m_nbins+1, 0.);
	  m_z_sigma.resize(m_nbins+1, 0.);
	  m_fact_scale.resize(m_nbins+1, -1.);
	  m_fact_z.resize(m_nbins+1, -1.);
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
      Pair1D_angular_lin_extra (const double thetaMin, const double thetaMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(binSize, 50, shift, angularUnits, angularWeight), Pair1D_angular(thetaMin, thetaMax, binSize, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_angular_lin_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_binSize();
	  m_scale_mean.resize(m_nbins+1, 0.);
	  m_scale_S.resize(m_nbins+1, 0.);
	  m_scale_sigma.resize(m_nbins+1, 0.);
	  m_z_mean.resize(m_nbins+1, 0.);
	  m_z_S.resize(m_nbins+1, 0.);
	  m_z_sigma.resize(m_nbins+1, 0.);
	  m_fact_scale.resize(m_nbins+1, -1.);
	  m_fact_z.resize(m_nbins+1, -1.);
	}
      
      /**
       *  @brief default destructor
       *  
       */
      ~Pair1D_angular_lin_extra () = default;
      
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
       *  
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;
  
      ///@}
      
    };


    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair1D_angular_log_extra Pair1D_extra.h "Headers/Pair1D_extra.h"
     *
     *  @brief The class Pair1D_angular_log_extra
     *
     *  This class is used to handle objects of type <EM> Pair1D_angular_log_extra
     *  </EM>.
     */
    class Pair1D_angular_log_extra : public virtual Pair1D_angular_extra, public virtual Pair1D_angular_log {
      
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_angular_log_extra ()
	{
	  m_pairType = PairType::_angular_log_;
	  m_pairInfo = PairInfo::_extra_;
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
      Pair1D_angular_log_extra (const double thetaMin, const double thetaMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(1., nbins, shift, angularUnits, angularWeight), Pair1D_angular(thetaMin, thetaMax, nbins, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_angular_log_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_nbins();
	  m_scale_mean.resize(m_nbins+1, 0.);
	  m_scale_S.resize(m_nbins+1, 0.);
	  m_scale_sigma.resize(m_nbins+1, 0.);
	  m_z_mean.resize(m_nbins+1, 0.);
	  m_z_S.resize(m_nbins+1, 0.);
	  m_z_sigma.resize(m_nbins+1, 0.);
	  m_fact_scale.resize(m_nbins+1, -1.);
	  m_fact_z.resize(m_nbins+1, -1.);
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
      Pair1D_angular_log_extra (const double thetaMin, const double thetaMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(binSize, 50, shift, angularUnits, angularWeight), Pair1D_angular(thetaMin, thetaMax, binSize, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_angular_log_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_binSize();
	  m_scale_mean.resize(m_nbins+1, 0.);
	  m_scale_S.resize(m_nbins+1, 0.);
	  m_scale_sigma.resize(m_nbins+1, 0.);
	  m_z_mean.resize(m_nbins+1, 0.);
	  m_z_S.resize(m_nbins+1, 0.);
	  m_z_sigma.resize(m_nbins+1, 0.);
	  m_fact_scale.resize(m_nbins+1, -1.);
	  m_fact_z.resize(m_nbins+1, -1.);
	} 
  
      /**
       *  @brief default destructor
       *  
       */
      ~Pair1D_angular_log_extra () = default;

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
       *  
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;

      ///@}
    
    };

  
    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair1D_comoving_extra Pair1D_extra.h "Headers/Pair1D_extra.h"
     *
     *  @brief The class Pair1D_comoving_extra
     *
     *  This class is used to handle objects of type <EM> Pair1D_comoving_extra
     *  </EM>.
     */
    class Pair1D_comoving_extra : public virtual Pair1D_extra, public virtual Pair1D_comoving {
  
    public:
  
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_comoving_extra () = default;
    
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
      Pair1D_comoving_extra (const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D_comoving(rMin, rMax, nbins, shift, angularUnits, angularWeight) {} 
  
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
      Pair1D_comoving_extra (const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D_comoving(rMin, rMax, binSize, shift, angularUnits, angularWeight) {} 
  
      /**
       *  @brief default destructor
       *  
       */
      virtual ~Pair1D_comoving_extra () = default;

      ///@}

    };

  
    // ============================================================================================
    // ============================================================================================

  
    /**
     *  @class Pair1D_comoving_lin_extra Pair1D_extra.h "Headers/Pair1D_extra.h"
     *
     *  @brief The class Pair1D_comoving_lin_extra
     *
     *  This class is used to handle objects of type <EM> Pair1D_comoving_lin_extra
     *  </EM>.
     */
    class Pair1D_comoving_lin_extra : public virtual Pair1D_comoving_extra, public virtual Pair1D_comoving_lin {

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_comoving_lin_extra ()
	{
	  m_pairType = PairType::_comoving_lin_;
	  m_pairInfo = PairInfo::_extra_;
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
      Pair1D_comoving_lin_extra (const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(1., nbins, shift, angularUnits, angularWeight), Pair1D_comoving(rMin, rMax, nbins, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comoving_lin_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_nbins();
	  m_scale_mean.resize(m_nbins+1, 0.);
	  m_scale_S.resize(m_nbins+1, 0.);
	  m_scale_sigma.resize(m_nbins+1, 0.);
	  m_z_mean.resize(m_nbins+1, 0.);
	  m_z_S.resize(m_nbins+1, 0.);
	  m_z_sigma.resize(m_nbins+1, 0.);
	  m_fact_scale.resize(m_nbins+1, -1.);
	  m_fact_z.resize(m_nbins+1, -1.);
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
      Pair1D_comoving_lin_extra (const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(binSize, 50, shift, angularUnits, angularWeight), Pair1D_comoving(rMin, rMax, binSize, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comoving_lin_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_binSize();
	  m_scale_mean.resize(m_nbins+1, 0.);
	  m_scale_S.resize(m_nbins+1, 0.);
	  m_scale_sigma.resize(m_nbins+1, 0.);
	  m_z_mean.resize(m_nbins+1, 0.);
	  m_z_S.resize(m_nbins+1, 0.);
	  m_z_sigma.resize(m_nbins+1, 0.);
	  m_fact_scale.resize(m_nbins+1, -1.);
	  m_fact_z.resize(m_nbins+1, -1.);
	} 
   
  
      /**
       *  @brief default destructor
       *  
       */
      ~Pair1D_comoving_lin_extra () = default;

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
       *  
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;
  
      ///@}
    
    };

 
    // ============================================================================================
    // ============================================================================================

  
    /**
     *  @class Pair1D_comoving_log_extra Pair1D_extra.h "Headers/Pair1D_extra.h"
     *
     *  @brief The class Pair1D_comoving_log_extra
     *
     *  This class is used to handle objects of type <EM> Pair1D_comoving_log_extra
     *  </EM>.
     */
    class Pair1D_comoving_log_extra : public virtual Pair1D_comoving_extra, public virtual Pair1D_comoving_log {
 
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_comoving_log_extra ()
	{
	  m_pairType = PairType::_comoving_log_;
	  m_pairInfo = PairInfo::_extra_;
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
      Pair1D_comoving_log_extra (const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(1., nbins, shift, angularUnits, angularWeight), Pair1D_comoving(rMin, rMax, nbins, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comoving_log_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_nbins();
	  m_scale_mean.resize(m_nbins+1, 0.);
	  m_scale_S.resize(m_nbins+1, 0.);
	  m_scale_sigma.resize(m_nbins+1, 0.);
	  m_z_mean.resize(m_nbins+1, 0.);
	  m_z_S.resize(m_nbins+1, 0.);
	  m_z_sigma.resize(m_nbins+1, 0.);
	  m_fact_scale.resize(m_nbins+1, -1.);
	  m_fact_z.resize(m_nbins+1, -1.);
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
      Pair1D_comoving_log_extra (const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(binSize, 50, shift, angularUnits, angularWeight), Pair1D_comoving(rMin, rMax, binSize, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comoving_log_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_binSize();
	  m_scale_mean.resize(m_nbins+1, 0.);
	  m_scale_S.resize(m_nbins+1, 0.);
	  m_scale_sigma.resize(m_nbins+1, 0.);
	  m_z_mean.resize(m_nbins+1, 0.);
	  m_z_S.resize(m_nbins+1, 0.);
	  m_z_sigma.resize(m_nbins+1, 0.);
	  m_fact_scale.resize(m_nbins+1, -1.);
	  m_fact_z.resize(m_nbins+1, -1.);
	}
  
      /**
       *  @brief default destructor
       *  
       */
      ~Pair1D_comoving_log_extra () = default;

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
       *  
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;

      ///@}
    
    };
  
    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair1D_comoving_multipoles_extra Pair1D_extra.h "Headers/Pair1D_extra.h"
     *
     *  @brief The class Pair1D_comoving_multipoles_extra
     *
     *  This class is used to handle objects of type <EM> Pair1D_comoving_multipoles_extra
     *  </EM>.
     */
    class Pair1D_comoving_multipoles_extra : public virtual Pair1D_extra, public virtual Pair1D_comoving_multipoles {
  
    public:
  
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_comoving_multipoles_extra () = default;
    
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
      Pair1D_comoving_multipoles_extra (const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D_comoving_multipoles(rMin, rMax, nbins, shift, angularUnits, angularWeight) {} 
  
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
      Pair1D_comoving_multipoles_extra (const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D_comoving_multipoles(rMin, rMax, binSize, shift, angularUnits, angularWeight) {} 
  
      /**
       *  @brief default destructor
       *  
       */
      virtual ~Pair1D_comoving_multipoles_extra () = default;

      ///@}

    };

    // ============================================================================================
    // ============================================================================================

    /**
     *  @class Pair1D_comoving_multipoles_lin_extra Pair1D_extra.h "Headers/Pair1D_extra.h"
     *
     *  @brief The class Pair1D_comoving_multipoles_lin_extra
     *
     *  This class is used to handle objects of type <EM> Pair1D_comoving_multipoles_lin_extra
     *  </EM>.
     */
    class Pair1D_comoving_multipoles_lin_extra : public virtual Pair1D_comoving_multipoles_extra, public virtual Pair1D_comoving_multipoles_lin {

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_comoving_multipoles_lin_extra ()
	{
	  m_pairType = PairType::_comoving_multipoles_lin_;
	  m_pairInfo = PairInfo::_extra_;
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
      Pair1D_comoving_multipoles_lin_extra (const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(1., nbins, shift, angularUnits, angularWeight), Pair1D_comoving_multipoles(rMin, rMax, nbins, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comoving_multipoles_lin_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_nbins();
	  m_PP1D.resize(3*(m_nbins+1), 0.);
	  m_PP1D_weighted.resize(3*(m_nbins+1), 0.);
	  m_scale_mean.resize(3*(m_nbins+1), 0.);
	  m_scale_S.resize(3*(m_nbins+1), 0.);
	  m_scale_sigma.resize(3*(m_nbins+1), 0.);
	  m_z_mean.resize(3*(m_nbins+1), 0.);
	  m_z_S.resize(3*(m_nbins+1), 0.);
	  m_z_sigma.resize(3*(m_nbins+1), 0.);
	  m_fact_scale.resize(3*(m_nbins+1), -1.);
	  m_fact_z.resize(3*(m_nbins+1), -1.);
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
      Pair1D_comoving_multipoles_lin_extra (const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(binSize, 50, shift, angularUnits, angularWeight), Pair1D_comoving_multipoles(rMin, rMax, binSize, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comoving_multipoles_lin_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_binSize();
	  m_PP1D.resize(3*(m_nbins+1), 0.);
	  m_PP1D_weighted.resize(3*(m_nbins+1), 0.);
	  m_scale_mean.resize(3*(m_nbins+1), 0.);
	  m_scale_S.resize(3*(m_nbins+1), 0.);
	  m_scale_sigma.resize(3*(m_nbins+1), 0.);
	  m_z_mean.resize(3*(m_nbins+1), 0.);
	  m_z_S.resize(3*(m_nbins+1), 0.);
	  m_z_sigma.resize(3*(m_nbins+1), 0.);
	  m_fact_scale.resize(3*(m_nbins+1), -1.);
	  m_fact_z.resize(3*(m_nbins+1), -1.);
	} 
   
  
      /**
       *  @brief default destructor
       *  
       */
      ~Pair1D_comoving_multipoles_lin_extra () = default;

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
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;

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
     *  @class Pair1D_comoving_multipoles_log_extra Pair1D_extra.h "Headers/Pair1D_extra.h"
     *
     *  @brief The class Pair1D_comoving_multipoles_log_extra
     *
     *  This class is used to handle objects of type <EM> Pair1D_comoving_multipoles_log_extra
     *  </EM>.
     */
    class Pair1D_comoving_multipoles_log_extra : public virtual Pair1D_comoving_multipoles_extra, public virtual Pair1D_comoving_multipoles_log {

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Pair1D_comoving_multipoles_log_extra ()
	{
	  m_pairType = PairType::_comoving_multipoles_log_;
	  m_pairInfo = PairInfo::_extra_;
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
      Pair1D_comoving_multipoles_log_extra (const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(1., nbins, shift, angularUnits, angularWeight), Pair1D_comoving_multipoles(rMin, rMax, nbins, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comoving_multipoles_log_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_nbins();
	  m_PP1D.resize(3*(m_nbins+1), 0.);
	  m_PP1D_weighted.resize(3*(m_nbins+1), 0.);
	  m_scale_mean.resize(3*(m_nbins+1), 0.);
	  m_scale_S.resize(3*(m_nbins+1), 0.);
	  m_scale_sigma.resize(3*(m_nbins+1), 0.);
	  m_z_mean.resize(3*(m_nbins+1), 0.);
	  m_z_S.resize(3*(m_nbins+1), 0.);
	  m_z_sigma.resize(3*(m_nbins+1), 0.);
	  m_fact_scale.resize(3*(m_nbins+1), -1.);
	  m_fact_z.resize(3*(m_nbins+1), -1.);
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
      Pair1D_comoving_multipoles_log_extra (const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr)
	: Pair1D(binSize, 50, shift, angularUnits, angularWeight), Pair1D_comoving_multipoles(rMin, rMax, binSize, shift, angularUnits, angularWeight)
	{
	  m_pairType = PairType::_comoving_multipoles_log_;
	  m_pairInfo = PairInfo::_extra_;
	  m_set_parameters_binSize();
	  m_PP1D.resize(3*(m_nbins+1), 0.);
	  m_PP1D_weighted.resize(3*(m_nbins+1), 0.);
	  m_scale_mean.resize(3*(m_nbins+1), 0.);
	  m_scale_S.resize(3*(m_nbins+1), 0.);
	  m_scale_sigma.resize(3*(m_nbins+1), 0.);
	  m_z_mean.resize(3*(m_nbins+1), 0.);
	  m_z_S.resize(3*(m_nbins+1), 0.);
	  m_z_sigma.resize(3*(m_nbins+1), 0.);
	  m_fact_scale.resize(3*(m_nbins+1), -1.);
	  m_fact_z.resize(3*(m_nbins+1), -1.);
	} 
   
  
      /**
       *  @brief default destructor
       *  
       */
      ~Pair1D_comoving_multipoles_log_extra () = default;

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
       *  
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2) override;

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
