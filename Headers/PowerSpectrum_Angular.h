/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Headers/PowerSpectrum_Angular.h
 *
 *  @brief The class PowerSpectrum_angular
 *
 *  This file defines the interface of the class
 *  PowerSpectrum_angular, used to measure the angular power spectrum
 */

#ifndef __POWSPECTRUMANG__
#define __POWSPECTRUMANG__

#include "TwoPointCorrelation1D_angular.h"


// ===================================================================================================


namespace cbl {   
  
  namespace measure {

    /**
     *  @brief The namespace of the <B> angular power spectrum </B>
     *  
     *  The \e measure::twopt namespace contains all the functions and
     *  classes to measure the angular power spectrum
     */
    namespace angularpk {
      
      /**
       * @enum cbl::measure::angularpk::Estimator
       *
       * @brief the angular two-point correlation estimator type
       *
       * @var cbl::measure::angularpk::Estimator::\_Fast\_
       *
       * @brief the faset estimator
       *
       * @var cbl::measure::angularpk::Estimator::\_SphericalArmonic\_
       *
       * @brief the spherical armonic estimator
       */
      enum class Estimator { _Fast_, _SphericalArmonic_ };
      
      /**
       *  @class PowerSpectrum_angular
       *
       *  @brief The class PowerSpectrum_angular
       *
       *  This class is used to handle objects of type <EM>
       *  PowerSpectrum_angular </EM>. It is used to measure the
       *  angular power spectrum, \f$C_l\f$
       */
      class PowerSpectrum_angular : public Measure {
	
      private:
	
	///vector containing the angular separation
	std::vector<double> m_theta;
	
	///vector containing the angular correlation function
	std::vector<double> m_w;
	
	///vector containing the angular correlation function error
	std::vector<double> m_w_err;
	
	///vector containing the multipoles at which power spectrum is computed
	std::vector<double> m_ell;
	
	/// the minimum multipole
	double m_ell_min;
	
	///the maximum multipole
	double m_ell_max;
	
	///the number of multipoles at which the angular power spectrum is computed
	double m_Nell;
	
	///vector containing the angular power spectrum
	std::vector<double> m_Celle;
	
	///the directory to read w
	std::string m_dir_correlation_input;
	
	///the file to read w
	std::string m_file_correlation_input;
	
	///the header lines to skip
	int m_n_lines_header;
	
	///the coordinate units of the input file
	cbl::CoordinateUnits m_inputUnits;
	
	///the directory to write w (if computed with cbl)
	std::string m_dir_correlation_output;
	
	///the file to write w (if computed with cbl)
	std::string m_file_correlation_output;
	
	///the coordinate units of the correlation function (if computed with cbl)
	cbl::CoordinateUnits m_angularUnits;
	
	///the bin type of the correlation function (if computed with cbl)
	cbl::BinType m_binType;
	
	///variable that store the constructor for the measure of angular correlation function
	std::shared_ptr<twopt::TwoPointCorrelation1D_angular> m_TwoPointCorrelation1D_angular;
	
	/**
	 *  @brief converts the input angle in radians
	 *
	 *  @param inputUnits input units
	 */
	void m_convert_angular_units (cbl::CoordinateUnits inputUnits);
	
	/**
	 *  @brief compute the dtheta*theta for the integral
	 *
	 *  @param binType the binning type, linear or logarithmic
	 *
	 *  @param i index to select the angle theta[i]
	 *
	 *  @return the dtheta*theta product
	 *  
	 */
	double m_dtheta_theta (const BinType binType, int i);
	
	/**
	 *  @brief measure the angular power spectrum errors 
	 *
	 *  @param binType the binning type, linear or logarithmic
	 *  
	 *  @param w_err angular correlation function vector error
	 *
	 *  @return the error vector
	 */
	std::vector<double> m_error (const BinType binType, std::vector<double> w_err);
	
      public:
	
	/**
	 *  @name Constructors/destructors
	 */
	///@{
	
	/**
	 *  @brief default constructor Power_spectrum_angular
	 */
	PowerSpectrum_angular () = default;
	
	/**
	 *  @brief constructor
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *
	 *  @param random of class Catalogue containing the random
	 *  data catalogue
	 *
	 *  @param ell_min Minimum angular power spectrum multipole
	 *
	 *  @param ell_max Maximum angular power spectrum multipole
	 *
	 *  @param Nell number of multipoles
	 *
	 *  @param binType binning type
	 *
	 *  @param thetaMin minimum angular separation used to count
	 *  the pairs
	 *
	 *  @param thetaMax maximum angular separation used to count
	 *  the pairs
	 *
	 *  @param nbins number of bins
	 *
	 *  @param shift shift parameter, i.e. the radial shift is
	 *  binSize*shift
	 *
	 *  @param angularUnits angular units (of the output angular
	 *  correlation file)
	 *
	 */     
	PowerSpectrum_angular (const catalogue::Catalogue data, const catalogue::Catalogue random, const double ell_min, const double ell_max, const int Nell, const BinType binType, const double thetaMin, const double thetaMax, const int nbins, const double shift, const CoordinateUnits angularUnits);   
	
	/**
	 *  @brief default destructor
	 */
	~PowerSpectrum_angular () = default;
	
	///@}
	
	/**
	 *  @brief generate vector which contains the multipoles at
	 *  which the angular power spectrum is computed
	 *
	 *  @param binType the binning type, linear or logarithmic 
	 */
	void set_ell_output (const BinType binType);
	
	/**
	 *  @brief measure the angular power spectrum
	 *
	 *  @param estimator the power spectrum estimator
	 *
	 *  @param dir_correlation_input Angular correlation function
	 *  input directory
	 *
	 *  @param file_correlation_input Angular correlation function
	 *  input file
	 *
	 *  @param n_lines_header the header lines to skip
	 *
	 *  @param inputUnits angular units in input file
	 *
	 *  @param input_binType the binning type, linear or
	 *  logarithmic
	 *
	 *  @param errorType the error type, _Poisson_, _Jackknife_ or
	 *  _Bootstrap_
	 *
	 *  @param dir_correlation_output Angular correlation function
	 *  ouput directory (if measured by CBL)
	 *
	 *  @param file_correlation_output Angular correlation
	 *  function ouput file (if measured by CBL)
	 *
	 *  @param nMocks number of mocks for bootstrap
	 *
	 */
	void measure (const Estimator estimator=Estimator::_Fast_, const std::string dir_correlation_input="", const std::string file_correlation_input="", const int n_lines_header=1, CoordinateUnits inputUnits=CoordinateUnits::_arcminutes_, BinType input_binType=BinType::_logarithmic_, cbl::measure::ErrorType errorType=ErrorType::_Poisson_, const std::string dir_correlation_output= par::defaultString, const std::string file_correlation_output="xi_angular.dat", const int nMocks=0);
	
	/**
	 * @brief get the private member m_theta
	 * @return the angular separation vector
	 */
	std::vector<double> theta () { return m_theta; }
	
	/**
	 * @brief get the private member Celle
	 * @return the angular power spectrum vector
	 */
	std::vector<double> Celle () { return m_Celle; }
	
	/**
	 * @brief get the private member m_theta
	 * @return the multipoles vector
	 */
	std::vector<double> ell () { return m_ell; }
	
	/**
	 *  @name Input/Output methods
	 */
	///@{
	
	/**
	 *  @brief write the angular power spectrum on file
	 *  @param dir output directory
	 *  @param file output file
	 */
	void write (const std::string dir, const std::string file);
	
	///@}
	
      };
    }
  }
}

#endif
