/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/TwoPointCorrelation2D.h
 *
 *  @brief The class TwoPointCorrelation2D
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation2D, used to measure the 2D two-point
 *  correlation function in Cartesian and polar coordinates
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINT2D__
#define __TWOPOINT2D__


#include "TwoPointCorrelation.h"


// ===================================================================================================


namespace cbl {

  namespace measure {

    namespace twopt {
    
      /**
       *  @class TwoPointCorrelation2D TwoPointCorrelation2D.h
       *  "Headers/TwoPointCorrelation2D.h"
       *
       *  @brief The class TwoPointCorrelation2D
       *
       *  This class is used to handle objects of type <EM>
       *  TwoPointCorrelation2D </EM>. It is used to measure the 2D
       *  two-point correlation function, in Cartesian and polar
       *  coordinates, \f$\xi(r_p,\pi)\f$ or \f$\xi(r,\mu)\f$, that is as
       *  a function of perpendicular, \f$r_p\f$, and parallel, \f$\pi\f$,
       *  line-of-sight separations, and as a function of absolute
       *  separation, \f$r=\sqrt{r_p^2+\pi^2}\f$, and the cosine of the
       *  angle between the separation vector and the line of sight,
       *  \f$\mu\equiv\cos\theta=s_\parallel/s\f$, respectively.
       */
      class TwoPointCorrelation2D : public TwoPointCorrelation {

      protected :
    
	/**
	 *  @name Internal input/output methods
	 */
	///@{
    
	/**
	 *  @brief write the number of pairs
	 *  @param PP pointer to an object of class Pair
	 *  @param dir output directory
	 *  @param file output file
	 *  
	 */
	void write_pairs (const std::shared_ptr<pairs::Pair> PP, const std::string dir, const std::string file) const override;

	/**
	 *  @brief read the number of pairs
	 *  @param [out] PP pointer to an object of class Pair
	 *  @param [in] dir vector of input directories
	 *  @param [in] file input file
	 *  
	 */
	void read_pairs (std::shared_ptr<pairs::Pair> PP, const std::vector<std::string> dir, const std::string file) const override;

	/**
	 *  @brief write the number of pairs
	 *  @param PP pointer to a vector of objects of class Pair
	 *  @param dir output directory
	 *  @param file output file
	 *  
	 */
	void write_pairs (const std::vector<std::shared_ptr<pairs::Pair>>  PP, const std::string dir, const std::string file) const override;

	/**
	 *  @brief read the number of pairs
	 *  @param [out] PP pointer to a vector of objects of class Pair
	 *  @param [in] dir vector of input directories
	 *  @param [in] file input file
	 *  
	 */
	void read_pairs (std::vector<std::shared_ptr<pairs::Pair>> PP, const std::vector<std::string> dir, const std::string file) const override;

	///@}

	/**
	 *  @name Member functions to measure the two-point correlation function
	 */
	///@{

	/**
	 *  @brief return a data object with extra info
	 *  
	 *  @param dd pointer to an object of type Pair containing the
	 *  data-data pairs
	 *  @param scale_D1 vector containing the binned scales along
	 *  the first dimension
	 *  @param scale_D2 vector containing the binned scales along
	 *  the second dimension
	 *  @param xi matrix containing the binned 2D two-point correlation function
	 *  @param error matrix containing the errors
	 *
	 *  @return pointer to an object of type Data
	 */
	std::shared_ptr<data::Data> data_with_extra_info (const std::shared_ptr<pairs::Pair> dd, const std::vector<double> scale_D1, const std::vector<double> scale_D2, const std::vector<std::vector<double>> xi, const std::vector<std::vector<double>> error) const;
      
	/**
	 *  @brief get a dataset containing the two-point correlation
	 *  function measured with the natural estimator, and its
	 *  Poisson errors
	 *  
	 *  @param dd pointer to an object of type Pair containing the
	 *  data-data pairs
	 *
	 *  @param rr pointer to an object of type Pair containing the
	 *  random-random pairs
	 *
	 *  @param nData number of objects in the data catalogue
	 *
	 *  @param nData_weighted weighted number of objects in the
	 *  data catalogue
	 *
	 *  @param nRandom number of objects in the random catalogue
	 *
	 *  @param nRandom_weighted weighted number of objects in the
	 *  random catalogue
	 *
	 *  @return pointer to an object of type Data
	 */
	std::shared_ptr<data::Data> correlation_NaturalEstimator (const std::shared_ptr<pairs::Pair> dd, const std::shared_ptr<pairs::Pair> rr, const int nData=0, const double nData_weighted=0., const int nRandom=0, const double nRandom_weighted=0.) override;

	/**
	 *  @brief get a dataset containing the two-point correlation
	 *  function measured with the Landy-Szalay estimator, and its
	 *  Poisson errors
	 *  
	 *  @param dd pointer to an object of type Pair containing the
	 *  data-data pairs
	 *
	 *  @param rr pointer to an object of type Pair containing the
	 *  random-random pairs
	 *
	 *  @param dr pointer to an object of type Pair containing the
	 *  data-random pairs
	 *
	 *  @param nData number of objects in the data catalogue
	 *
	 *  @param nData_weighted weighted number of objects in the data
	 *  catalogue
	 *
	 *  @param nRandom number of objects in the random catalogue
	 *
	 *  @param nRandom_weighted weighted number of objects in the
	 *  random catalogue
	 *
	 *  @return pointer to an object of type Data
	 */
	std::shared_ptr<data::Data> correlation_LandySzalayEstimator (const std::shared_ptr<pairs::Pair> dd, const std::shared_ptr<pairs::Pair> rr, const std::shared_ptr<pairs::Pair> dr, const int nData=0, const double nData_weighted=0., const int nRandom=0, const double nRandom_weighted=0.) override;

	/**
	 *  @brief measure the jackknife resampling of the two-point
	 *  correlation function, &xi;(r)
	 *
	 *  @param dd vector of data-data pairs, divided per regions
	 *
	 *  @param rr vector of random-random pairs, divided per regions
	 *
	 *  @return vector of pointers to objects of type Data
	 */
	std::vector<std::shared_ptr<data::Data>> XiJackknife (const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr) override;

	/**
	 *  @brief measure the jackknife resampling of the two-point correlation
	 *  function, &xi;(r)         
	 *
	 *  @param dd vector of data-data pairs, divided per regions
	 *
	 *  @param rr vector of random-random pairs, divided per regions
	 *
	 *  @param dr vector of random-random pairs, divided per regions   
	 *
	 *  @return vector of pointers to objects of type Data
	 */
	std::vector<std::shared_ptr<data::Data>> XiJackknife (const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr, const std::vector<std::shared_ptr<pairs::Pair>> dr) override;

	/**
	 *  @brief measure the bootstrap resampling of the two-point correlation
	 *  function, &xi;(r)  
	 *
	 *  @param nMocks number of bootstrap resampling
	 *
	 *  @param dd vector of data-data pairs, divided per regions
	 *
	 *  @param rr vector of random-random pairs, divided per
	 *  regions
	 *
	 *  @param seed the seed for random number generation
	 *
	 *  @return vector of pointers to objects of type Data
	 */
	std::vector<std::shared_ptr<data::Data>> XiBootstrap (const int nMocks, const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr, const int seed=3213) override;

	/**
	 *  @brief measure the bootstrap resampling of the two-point
	 *  correlation function, &xi;(r)
	 *
	 *  @param nMocks number of bootstrap resampling
	 *
	 *  @param dd vector of data-data pairs, divided per regions
	 *
	 *  @param rr vector of random-random pairs, divided per
	 *  regions
	 *
	 *  @param dr vector of random-random pairs, divided per
	 *  regions
	 *
	 *  @param seed the seed for random number generation
	 * 
	 *  @return vector of pointers to objects of type Data
	 */
	std::vector<std::shared_ptr<data::Data>> XiBootstrap (const int nMocks, const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr, const std::vector<std::shared_ptr<pairs::Pair>> dr, const int seed=3213) override;
      
	///@}


      public:
    
	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *  2D
	 */
	TwoPointCorrelation2D () { m_dataset = data::Data::Create(data::DataType::_2D_); }

	/**
	 *  @brief constructor
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *  2D
	 */
	TwoPointCorrelation2D (const catalogue::Catalogue data, const catalogue::Catalogue random, const bool compute_extra_info=false, const double random_dilution_fraction=1.) 
	  : TwoPointCorrelation(data, random, compute_extra_info, random_dilution_fraction)
	  { m_dataset = (!compute_extra_info) ? data::Data::Create(data::DataType::_2D_) : data::Data::Create(data::DataType::_2D_extra_); }

	/**
	 *  @brief default destructor
	 *  
	 */
	~TwoPointCorrelation2D () = default;

	///@}

        /**
	 *  @name Member functions to get the private/protected members
	 */
	///@{

	/**
	 *  @brief get the protected member m_x
	 *  @return the x coordinates
	 */
	std::vector<double> xx () const override
	  { return m_dataset->xx(); }

	/**
	 *  @brief get the protected member m_y
	 *  @return the y coordinates
	 */
	std::vector<double> yy () const override
	  { return m_dataset->yy(); }

	/**
	 *  @brief get the protected member m_fxy
	 *  @return the binned correlation function 
	 */
	std::vector<std::vector<double>> xi2D () const override
	  { std::vector<std::vector<double>> vv; m_dataset->get_data(vv); return vv; }

	/**
	 *  @brief get the protected member m_error_fxy
	 *  @return the error on the binned correlation function
	 *  function
	 */
	std::vector<std::vector<double>> error2D () const override
	  { std::vector<std::vector<double> > vv; m_dataset->get_error(vv); return vv; }
      
	///@}


	/**
	 *  @name Member functions to count measure the two-point correlation function
	 */
	///@{

	/**
	 *  @brief measure the two-point correlation function
	 *
	 *  @param errorType type
	 *  
	 *  @param dir_output_pairs output directory used to store the
	 *  number of pairs
	 *
	 *  @param dir_input_pairs vector of input directories used to
	 *  store the number of pairs (if the pairs are read from files)
	 *
	 *  @param dir_output_resample output directory of the
	 *  resampled correlation functions; if an empty string
	 *  (i.e. "" or "NULL") is provided, no output will be stored
	 *
	 *  @param nMocks number of resampling for bootstrap
	 *
	 *  @param count_dd true &rarr; count the number of data-data
	 *  opairs; false &rarr; read the number of data-data pairs from
	 *  file
	 *
	 *  @param count_rr true &rarr; count the number of
	 *  random-random opairs; false &rarr; read the number of
	 *  random-random pairs from file
	 *
	 *  @param count_dr true &rarr; count the number of data-random
	 *  pairs; false &rarr; read the number of data-random pairs
	 *
	 *  @param tcount true &rarr; activate the time counter; false &rarr;
	 *  no time counter
	 *
	 *  @param estimator the estimator used to measure the two-point
	 *  correlation function
	 *
	 *  @param fact factor used to compute the cell size of the
	 *  chain mesh: it is multiplied by the maximum distance
	 *  considered for the couples and can be setted by the user
	 *  to optimize the count of the couples
	 *
	 *  @param seed the seed for random number generation
	 */
	virtual void measure (const ErrorType errorType=ErrorType::_Poisson_, const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample=par::defaultString, const int nMocks=0, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const double fact=0.1, const int seed=3213) = 0;
      
	///@}
      

	/**
	 *  @name Input/Output member functions (customized in all the derived classes)
	 */
	///@{

	/**
	 *  @brief read the measured two-point correlation
	 *  @param dir input directory
	 *  @param file input file
	 */
	virtual void read (const std::string dir, const std::string file)
	{ (void)dir; (void)file; ErrorCBL("", "read", "TwoPointCorrelation2D.h"); }	

	/**
	 *  @brief write the measured two-point correlation
	 *  @param dir output directory
	 *  @param file output file
	 *  @param rank cpu index (for MPI usage)
	 */
	virtual void write (const std::string dir, const std::string file, const int rank=0) const
	{ (void)dir; (void)file; (void)rank; ErrorCBL("", "write", "TwoPointCorrelation2D.h"); }	

	///@}


	/**
	 *  @name Member functions to compute, read and write the covariance matrix (customised in all the derived classes)
	 */
	///@{ 

	/**
	 *  @brief read the measured covariance matrix
	 *  @param dir input directory
	 *  @param file input file
	 */
	virtual void read_covariance (const std::string dir, const std::string file)
	{ (void)dir; (void)file; ErrorCBL("", "read_covariance", "TwoPointCorrelation2D.h"); }

	/**
	 *  @brief write the measured two-point correlation
	 *  @param dir output directory
	 *  @param file output file
	 */
	virtual void write_covariance (const std::string dir, const std::string file) const
	{ (void)dir; (void)file; ErrorCBL("", "write_covariance", "TwoPointCorrelation2D.h"); }
      
	/**
	 *  @brief compute the covariance matrix
	 *  @param xi vector containing the measure correlation
	 *  functions used to compute the covariance matrix
	 *  @param JK true &rarr; compute the jackknife covariance
	 *  matrix; false compute the standard covariance matrix
	 */
	virtual void compute_covariance (const std::vector<std::shared_ptr<data::Data>> xi, const bool JK)
	{ (void)xi; (void)JK; ErrorCBL("", "compute_covariance", "TwoPointCorrelation2D.h"); }
 
	/**
	 *  @brief compute the covariance matrix
	 *  @param file vector containing the input files with the
	 *  measured correlation functions used to compute the
	 *  covariance matrix
	 *  @param JK true &rarr; compute the jackknife covariance
	 *  matrix; false compute the standard covariance matrix
	 */
	virtual void compute_covariance (const std::vector<std::string> file, const bool JK)
	{ (void)file; (void)JK; ErrorCBL("", "compute_covariance", "TwoPointCorrelation2D.h"); }
      
	///@} 

      };
    }
  }
}

#endif
