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
	 *  @brief read mixing matrix from file
	 *
	 *  @param dir input directory
	 *
	 *  @param file input file
	 *
	 *  @param ll vector of multipoles
	 *
	 *  @param matrix mixing_matrix
	 *
	 */
      void read_mixing_matrix (const std::string dir, const std::string file, std::vector<double> &ll, std::vector<std::vector<double>> &matrix);
       
      /**
       * @enum cbl::measure::angularpk::AngularEstimator
       *
       * @brief the angular two-point correlation estimator type
       *
       * @var cbl::measure::angularpk::AngularEstimator::\_Fast\_
       *
       * @brief the fast estimator
       *
       * @var cbl::measure::angularpk::AngularEstimator::\_SphericalArmonic\_
       *
       * @brief the spherical armonic estimator
       */
      enum class AngularEstimator { _Fast_, _SphericalArmonic_ };
      
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

	///vector containing the x coordinate read from file
	std::vector<double> m_x;
	
	///vector containing the y coordinate read from file
	std::vector<double> m_y;
	
	///vector containing the y errors read from file
	std::vector<double> m_y_err;
	
	///vector containing the multipoles at which power spectrum is computed
	std::vector<double> m_l;

	///vector containing the angular power spectrum
	std::vector<double> m_AngularPowerSpectrum;
	
	///vector containing the average of the multipoles
	std::vector<double> m_l_av;

	///vector containing the angular power spectrum average
	std::vector<double> m_AngularPowerSpectrum_av;
	
	///the minimum multipole
	double m_l_min;
	
	///the maximum multipole
	double m_l_max;

	///the number of multipoles at which the angular power spectrum is computed
	double m_Nl;
	
	///the number of multipoles at which the angular power spectrum is computed
	double m_nbins;

	///vector containing the angular power spectrum of the window function
	std::vector<double> m_Wl;

	///the mixing matrix
	std::vector<std::vector<double>> m_RR;
	
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

	///variable that store the constructor for the catalogue
	std::shared_ptr<catalogue::Catalogue> m_catalogue;

	///vector containing the colatitude vector of the mask
	std::vector<double> m_theta_mask;

	///vector containing the min colatitude vector of the mask
	std::vector<double> m_theta_mask_min;
	
	///vector containing the max colatitude vector of the mask
	std::vector<double> m_theta_mask_max;
	
	///vector containing the RA vector of the mask
	std::vector<double> m_RA_mask;
	
	///vector containing the min RA vector of the mask
	std::vector<double> m_RA_mask_min;
	
	///vector containing the max RA vector of the mask
	std::vector<double> m_RA_mask_max;

	///pixel area
	double m_pixel_area;
	
	///the survey area
	double m_survey_area;
	
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

	/**
	 *  @brief set the catalogue for Spherical armonic estimator 
	 *
	 *  @param catalogue the catalogue to set
	 *  
	 */
	
	void set_catalogue (cbl::catalogue::Catalogue catalogue);


	/**
	 *  @brief compute the average of the power spectrum
	 *
	 */
	void m_averagePowerSpectrum ();

	
	/**
	 *  @brief read the angular mask from file
	 *
	 *  @param mask_file the file with pixel mask positions
	 *
	 *  @param mask_type "Healpix" in case of binary mask produced with pixelfunc.pix2ang in healpix. Else other case
	 *
	 *  @param pixel_area the pixel area of healpix mask
	 *
	 *  @param n_lines_header the header lines to skip
	 *
	 */
	void m_read_angular_mask (const std::string mask_file, std::string mask_type, double pixel_area, int n_lines_header=1);
	
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
	 *  @param l_min Minimum angular power spectrum multipole
	 *
	 *  @param l_max Maximum angular power spectrum multipole
	 *
	 *  @param Nl number of multipoles
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
	 *  @param correlation_binType the bin type of the 
	 *  correlation function
	 *
	 */     
	PowerSpectrum_angular (const catalogue::Catalogue data, const catalogue::Catalogue random, const double l_min, const double l_max, const int Nl, const BinType binType, const double thetaMin, const double thetaMax, const int nbins, const double shift, const CoordinateUnits angularUnits, const BinType correlation_binType);   
	
	/**
	 *  @brief constructor
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *
	 *  @param l_min Minimum angular power spectrum multipole
	 *
	 *  @param l_max Maximum angular power spectrum multipole
	 *
	 *  @param bandwidth the bandwidth
	 *
	 *  @param mask_file the file containing the pixel center positions of the mask in colatitude-RA for every pixel (obtained with Healpix pixelfunc.pix2ang), or theta, theta_min, theta_max, ra, ra_min, ra_max and pixel_area, for every pixel, i.e. pixel centers, borders and area. Leave it empty if there is not a mask file, in this case the survey area will be computed as \f$ (cos(\theta_{min})-cos(\theta_{max}))*(RA_{max}-RA_{min})\f$
	 *
	 *  @param mask_type the mask_type. "Healpix" if the mask file is obtained with the healpix function pixelfunc.pix2ang. Anything else if themask is obtained in other way

	 *  @param pixel_area the average pixel area 
	 *
	 *  @param n_lines_header the number of lines to skip when reading the mask file
	 */     
	PowerSpectrum_angular (const catalogue::Catalogue data, const double l_min, const double l_max, const int bandwidth=1, const std::string mask_file="", const std::string mask_type="", const double pixel_area=0, const int n_lines_header=1);   
	
	
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
	void set_l_output (const BinType binType);
	
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
	void measure (const AngularEstimator estimator=AngularEstimator::_Fast_, const std::string dir_correlation_input="", const std::string file_correlation_input="", const int n_lines_header=1, CoordinateUnits inputUnits=CoordinateUnits::_arcminutes_, BinType input_binType=BinType::_logarithmic_, cbl::measure::ErrorType errorType=ErrorType::_Poisson_, const std::string dir_correlation_output= par::defaultString, const std::string file_correlation_output="xi_angular.dat", const int nMocks=0);

	/**
	 *  @brief measure the angular power spectrum with fast estimator
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
	void measureFast (const std::string dir_correlation_input="", const std::string file_correlation_input="", const int n_lines_header=1, CoordinateUnits inputUnits=CoordinateUnits::_arcminutes_, BinType input_binType=BinType::_logarithmic_, cbl::measure::ErrorType errorType=ErrorType::_Poisson_, const std::string dir_correlation_output= par::defaultString, const std::string file_correlation_output="xi_angular.dat", const int nMocks=0);

	/**
	 *  @brief measure the angular power spectrum with spherical 
	 *  armonic estimator
	 *
	 */
	void measureSphericalArmonic ();
	
	/**
	 * @brief get the private member m_theta
	 * @return the angular separation vector
	 */
	std::vector<double> theta () { return m_theta; }

	/**
	 * @brief get the private member Wl
	 * @return the angular power spectrum of the window function vector
	 */
	std::vector<double> Wl () { return m_Wl; }

	/**
	 * @brief get the private member AngularPowerSpectrum
	 * @return the angular power spectrum vector
	 */
	std::vector<double> AngularPowerSpectrum () { return m_AngularPowerSpectrum; }
	
	/**
	 * @brief get the private member m_l
	 * @return the multipoles vector
	 */
	std::vector<double> l () { return m_l; }
	
	/**
	 *  @brief compute the mixing matrix
	 *
	 *  @param dir_window_input input directory harmonic coefficients of the mask 
	 *
	 *  @param file_window_input input file harmonic coefficients of the mask 
	 *
	 */
	void compute_mixing_matrix(std::string dir_window_input="", std::string file_window_input="");

	/**
	 *  @brief include the binary angular mask
	 *
	 *  @param theta the colatitude to evaluate 
	 *
	 *  @param RA the right ascension to evaluate
	 *
	 *  @return 0 if theta, RA is within a masked pixel, 1 otherwise
	 *
	 */
	double angular_mask(double theta, double RA);

	/**
	 *  @name Input/Output methods
	 */
	///@{
	
	/**
	 *  @brief write the angular power spectrum on file
	 *
	 *  @param dir output directory
	 *
	 *  @param file output file
	 *
	 */
	void write (const std::string dir, const std::string file);

	/**
	 *  @brief write the angular power spectrum on file
	 *
	 *  @param dir output directory
	 *
	 *  @param file output file
	 *
	 */
	void write_average (const std::string dir, const std::string file);
	
	/**
	 *  @brief write the angular power spectrum on file
	 *
	 *  @param dir output directory
	 *
	 *  @param file output file
	 *	 
	 *  @param store_window true → store the harmonic coefficients of the mask; false → do not store the harmonic coefficients of the mask 
	 *
	 */
	void write_mixing_matrix (const std::string dir, const std::string file,  bool store_window=false);
 
	/**
	 *  @brief read data from file
	 *
	 *  @param dir output directory
	 *
	 *  @param file output file
	 *
	 *  @param column_x vector of x column
	 *
	 *  @param column_y vector of y column
	 *
	 *  @param column_error vector of error column
	 *
	 *  @param n_lines_header the header lines to skip
	 *
	 */
      void read (const std::string dir, const std::string file, const std::vector<int> column_x={1}, const std::vector<int> column_y={2}, const std::vector<int> column_error={3}, const int n_lines_header=1);
	
	///@}
	
      };
    }
  }
}

#endif
