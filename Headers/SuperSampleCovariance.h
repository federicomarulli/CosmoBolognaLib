/********************************************************************
 *  Copyright (C) 2021 by Federico Marulli and Giorgio Lesci        *
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
 *  @file Headers/SuperSampleCovariance.h
 *
 *  @brief The class SuperSampleCovariance
 *
 *  This file defines the interface of the class SuperSampleCovariance,
 *  used to compute the \f$S_{ij}\f$ matrix for the super-sample covariance.
 *
 *  Given two redshift bins, labelled as \f$i\f$ and \f$j\f$, \f$S_{ij}\f$ is:
 *
 *  \f$ S_{ij} = \frac{1}{\Omega} \frac{1}{2\pi^2} 
 *  \int {\rm d} k\,\, k^2 P(k) \frac{U_i(k)}{I_i} \frac{U_j(k)}{I_j}, \f$
 *  
 *  where \f$\Omega\f$ is the survey area, \f$P(k)\f$ is the power spectrum,
 *  and \f$U_i(k)\f$ and \f$I_i\f$ are expressed as:
 *
 *  \f$ U_i(k) = \int {\rm d} V_i \,\, W^2_i g(z_j) j_0(kr_j), \f$
 *
 *  \f$ I_i = \int {\rm d} V_i \,\, W^2_i, \f$
 *
 *  where \f$V_i\f$ is the comoving volume within the \f$i\f$-th
 *  redshift bin, \f$g\f$ is the growth factor, \f$j_0\f$ the
 *  Bessel spherical function, and \f$W_i\f$ is the window function.
 *
 *  Physical units are forced.
 *
 *  This code is a reimplementation of the Python code presented in Lacasa & Grain 2019.
 *  The original code can be found here: https://github.com/fabienlacasa/PySSC
 *
 *  @author Giorgio Lesci
 *
 *  @author giorgio.lesci2@unibo.it
 */


#ifndef __SSC__
#define __SSC__


#include "Cosmology.h"
#include "Modelling.h"


namespace cbl {

  namespace statistics {

     
    /**
     *  @class SuperSampleCovariance SuperSampleCovariance.h
     *  "Headers/SuperSampleCovariance.h"
     *
     *  @brief The class SuperSampleCovariance
     *
     *  This is the base class used to manage 
     *  the super-sample covariance matrix
     */
    class SuperSampleCovariance
    {
      
    protected:
    
      /// pointer to the response function of the probe
      std::vector<std::shared_ptr<statistics::Model>> m_response_func;
    
      /// pointer to the Cosmology object
      std::shared_ptr<cosmology::Cosmology> m_cosmo;
      
      /// names of the cosmological parameters
      std::vector<cbl::cosmology::CosmologicalParameter> m_cosmo_param;
      
      /// method used to compute the power spectrum
      std::string m_method_Pk;
      
      /// linear o non-linear power-spectrum
      bool m_NL;
      
      /// store Pk output
      bool m_store_output;
      
      /// number of redshift bins
      int m_nbins;
      
      /// number of redshift steps where the S matrix is computed
      int m_nsteps;
      
      /// vector containing the redshift values where the S matrix is computed
      std::vector<double> m_redshifts;
      
      /// survey area in steradians
      double m_area;
      
      /// precision of the \f$log k\f$ array
      double m_precision;
      
      /// window functions in the redshift bins
      std::vector<std::vector<double>> m_windows;

      /**
       * @brief compute the top-hat window functions 
       * in the redshift bins
       *
       * @param delta_z redshift step used to construct the
       *  window function
       *
       * @param redshift_edges redshift bin edges
       *
       */
      void m_compute_topHat_window (const double delta_z, const std::vector<double> redshift_edges);
      
      /**
       * @brief compute the Gaussian window functions 
       * in the redshift bins
       *
       * @param delta_z redshift step used to construct the
       *  window function
       *
       *  @param W_mean vector of mean values for the 
       *  Gaussian window functions
       *
       *  @param W_std vector of standard deviation values
       *  for the Gaussian window functions
       *
       */
      void m_compute_gaussian_window (const double delta_z, const std::vector<double> W_mean, const std::vector<double> W_std);


    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief Default constructor, used to compute the \f$S_{ij}\f$ matrix for the super-sample covariance.
       *
       *  Given two redshift bins, labelled as \f$i\f$ and \f$j\f$, \f$S_{ij}\f$ is:
       *
       *  \f$ S_{ij} = \frac{1}{\Omega} \frac{1}{2\pi^2} 
       *  \int {\rm d} k\,\, k^2 P(k) \frac{U_i(k)}{I_i} \frac{U_j(k)}{I_j}, \f$
       *  
       *  where \f$\Omega\f$ is the survey area, \f$P(k)\f$ is the power spectrum,
       *  and \f$U_i(k)\f$ and \f$I_i\f$ are expressed as:
       *
       *  \f$ U_i(k) = \int {\rm d} V_i \,\, W^2_i g(z_j) j_0(kr_j), \f$
       *
       *  \f$ I_i = \int {\rm d} V_i \,\, W^2_i, \f$
       *
       *  where \f$V_i\f$ is the comoving volume within the \f$i\f$-th
       *  redshift bin, \f$g\f$ is the growth factor, \f$j_0\f$ the
       *  Bessel spherical function, and \f$W_i\f$ is the window function.
       *
       *  Physical units are forced, in set_SSC.
       *
       *  @param modelling pointers to Modelling objects
       *  
       */
      SuperSampleCovariance (std::vector<std::shared_ptr<cbl::modelling::Modelling>> modelling);

      /**
       *  @brief default destructor
       */
      virtual ~SuperSampleCovariance () = default;

      ///@}
      
      
      /**
       *  @name Member functions to set the private/protected members
       */
      ///@{
      
      /**
       *  @brief Set the
       *  super-sample covariance matrix \f$S_{ij}\f$,
       *  assuming top-hat window functions. The top-hat
       *  window function can be used when the redshift errors
       *  are smaller than the width of the redshift bins
       *
       *  @param cosm the Cosmology object
       *
       *  @param cosmo_param the cosmological parameters
       *  set in the modelling (useful only if a modelling
       *  is performed)
       *
       *  @param redshift_edges redshift bin edges
       *
       *  @param area effective area of the survey,
       *  expressed in squared degrees
       *
       *  @param method_Pk the method used to compute
       *  the power spectrum
       *
       *  @param delta_z redshift step used to construct the
       *  window function
       *  
       *  @param precision precision of the \f$\log k\f$ array, 
       *  defined by \f$2^{\rm precision}\f$
       *
       *  @param NL false \f$\rightarrow\f$ linear power spectrum;
       *  true \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param store_output if true the output files created 
       *  by the Boltzmann solver are stored; if false the 
       *  output files are removed
       *
       */
      void set_SSC (cbl::cosmology::Cosmology cosm, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<double> redshift_edges, const double area, const std::string method_Pk="EisensteinHu", const double delta_z=0.001, const double precision=10, const bool NL=false, const bool store_output=false);
      
      
      /**
       *  @brief Set the
       *  super-sample covariance matrix \f$S_{ij}\f$,
       *  assuming Gaussian window functions.
       *
       *  @param cosm the Cosmology object
       *
       *  @param cosmo_param the cosmological parameters
       *  set in the modelling (useful only if a modelling
       *  is performed)
       *
       *  @param area effective area of the survey,
       *  expressed in squared degrees
       *
       *  @param W_mean vector of mean values for the 
       *  Gaussian window functions, corresponding to the
       *  centre of the redshift bins
       *
       *  @param W_std vector of standard deviation values
       *  for the Gaussian window functions
       *
       *  @param method_Pk the method used to compute
       *  the power spectrum
       *
       *  @param delta_z redshift step used to construct the
       *  window function
       *  
       *  @param precision precision of the \f$\log k\f$ array, 
       *  defined by \f$2^{\rm precision}\f$
       *
       *  @param NL false \f$\rightarrow\f$ linear power spectrum;
       *  true \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param store_output if true the output files created 
       *  by the Boltzmann solver are stored; if false the 
       *  output files are removed
       *
       */
      void set_SSC (cbl::cosmology::Cosmology cosm, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const double area, const std::vector<double> W_mean, const std::vector<double> W_std, const std::string method_Pk="EisensteinHu", const double delta_z=0.001, const double precision=10, const bool NL=false, const bool store_output=false);
      
      
      ///@}
      
      /**
       *  @name Member functions to compute \f$S_{ij}\f$
       */
      ///@{
      
      
      /**
       * @brief compute the S_ij matrix
       *
       * @param cosmo the cosmological model
       *
       * @return the covariance matrix
       *
       */
      std::vector<std::vector<double>> compute_Sij (cbl::cosmology::Cosmology cosmo) const;
      
      /**
       *  @brief get \f$S_{ij}\f$
       *
       *  @param parameter the parameters of interest in the covariance matrix
       *
       *  @return the covariance matrix
       */
      std::vector<std::vector<double>> operator () (std::vector<double> &parameter) const;
      
      /**
       *  @brief get the response function of all the probes
       *
       *  @param xx the points where the response is evaluated
       *
       *  @param parameter the parameters of interest in the covariance matrix
       *
       *  @return the values of the response function
       */
      std::vector<std::vector<double>> get_response (std::vector<std::vector<double>> xx, std::vector<double> &parameter) const;
      
      /**
       *  @brief get the response function of the i-th probe
       *
       *  @param i index of the probe
       *
       *  @param xx the points where the response is evaluated
       *
       *  @param parameter the parameters of interest in the covariance matrix
       *
       *  @return the values of the response function
       */
      std::vector<double> get_response (int i, std::vector<double> xx, std::vector<double> &parameter) const;
       
       ///@}
       
       
      /**
       *  @name Member functions to get the private members of the class
       */
      ///@{
      
      
      /**
       * @brief return the window functions
       *
       * @return the window function
       */
       std::vector<std::vector<double>> get_window_function ();
       
       /**
       * @brief return the dimension of the Sij matrix
       *
       * @return the dimension of the Sij matrix
       */
       int Sij_dimension ();
       
       ///@}
       
       
       /**
       *  @name Member functions that write on file the products of the class
       */
      ///@{
      
      
      /**
       * @brief write the window functions on file
       *
       * @param dir output directory
       *
       * @param file output file
       *
       */
       void write_window_function (const std::string dir, const std::string file);
       
       /**
       * @brief Write the \f$S_{ij}\f$ matrix on file.
       *
       * Since the objects of the class SuperSampleCovariance are
       * used in parallelized MCMC computations, \f$S_{ij}\f$ is
       * not set as a member of this class. Therefore, this function
       * calculates \f$S_{ij}\f$
       *
       * @param dir output directory
       *
       * @param file output file
       *
       */
       void write_Sij (const std::string dir, const std::string file);
       
       ///@}
       
      
    };


  }
}

#endif
