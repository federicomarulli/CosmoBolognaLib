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

#include "PkXi.h"


// ===================================================================================================


namespace cbl {

  namespace cosmology {

     
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
    
      /// if true, the window function is set
      bool m_isSet_window = false;
      
      /// if true, the SSC is computed in the full-sky case
      bool m_isFullSky;
      
      /// Cls of the survey mask
      std::vector<double> m_Cl_mask;
      
      /// pointer to the Cosmology object
      std::shared_ptr<cosmology::Cosmology> m_cosmology = NULL;
      
      /// names of the cosmological parameters
      std::vector<std::string> m_cosmo_param;
      
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
      
      /// precision of the \f$log k\f$ array
      double m_precision;
      
      /// window functions in the redshift bins
      std::vector<std::vector<double>> m_windows;
      
      /**
       * @brief compute the S_ij matrix
       *
       * @return the super-sample covariance matrix
       */
      std::vector<std::vector<double>> m_compute_Sij () const;


    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief Constructor used to compute the \f$S_{ij}\f$
       *  matrix in the full-sky case.
       *
       *  Given two redshift bins, labelled as \f$i\f$ and \f$j\f$, \f$S_{ij}\f$ is:
       *
       *  \f$ S_{ij} = \frac{1}{2\pi^2} 
       *  \int {\rm d} k\,\, k^2 P(k) \frac{U_i(k)}{I_i} \frac{U_j(k)}{I_j}, \f$
       *  
       *  where \f$P(k)\f$ is the power spectrum computed at \f$z=0\f$,
       *  and \f$U_i(k)\f$ and \f$I_i\f$ are expressed as:
       *
       *  \f$ U_i(k) = \int {\rm d} V_i \,\, W_i(z) g(z) j_0(kr), \f$
       *
       *  \f$ I_i = \int {\rm d} V_i \,\, W_i(z), \f$
       *
       *  where \f$V_i\f$ is the comoving volume within the \f$i\f$-th
       *  redshift bin, \f${\rm d}V = \frac{{\rm d}V}{{\rm d}z}{\rm d}z\f$,
       *  \f$g\f$ is the growth factor, \f$j_0\f$ the
       *  spherical Bessel function of zeroth order,
       *  \f$r\f$ is the comoving distance,
       *  and \f$W_i\f$ is the window function.
       *
       *  We note that
       *
       *  \f$ S_{ij} = \frac{C_S(\ell=0,i,j)}{4\pi}\,, \f$
       *
       *  with
       *
       *  \f$ C_S(\ell=0,i,j) = \frac{2}{\pi} \int{\rm d}k\,k^2
       *  P(k)\frac{U_i(k)}{I_i} \frac{U_j(k)}{I_j}\,. \f$
       *
       *  Physical units are forced.
       *
       *  @param cosmology the Cosmology object
       *
       *  @param cosmo_param the cosmological parameters
       *  set in the modelling (useful only if a modelling
       *  is performed)
       *
       *  @param method_Pk the method used to compute
       *  the power spectrum
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
      SuperSampleCovariance (std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const std::vector<std::string> cosmo_param, const std::string method_Pk="EisensteinHu", const double precision=10, const bool NL=false, const bool store_output=false);
      
      /**
       *  @brief Constructor used to compute the \f$S_{ij}\f$
       *  matrix in the partial-sky case.
       *
       *  Given two redshift bins, labelled as \f$i\f$ and \f$j\f$, \f$S_{ij}\f$ is:
       *
       *  \f$ S_{ij} = \frac{1}{8\pi^3 f_{\rm sky}^2} \sum_\ell (2\ell+1) \, \mathcal{C}_{\rm m}(\ell)
       *  \int {\rm d} k\,\, k^2 P(k) \frac{U_i(k,\ell)}{I_i} \frac{U_j(k,\ell)}{I_j}, \f$
       *  
       *  where \f$f_{\rm sky}\f$ is the sky fraction covered by the survey,
       *  \f$\mathcal{C}_{\rm m}(\ell)\f$ is the mask angular power spectrum,
       *  \f$P(k)\f$ is the power spectrum computed at \f$z=0\f$,
       *  and \f$U_i(k,\ell)\f$ and \f$I_i\f$ are expressed as:
       *
       *  \f$ U_i(k,\ell) = \int {\rm d} V_i \,\, W_i(z) g(z) j_\ell(kr), \f$
       *
       *  \f$ I_i = \int {\rm d} V_i \,\, W_i(z), \f$
       *
       *  where \f$V_i\f$ is the comoving volume within the \f$i\f$-th
       *  redshift bin, \f${\rm d}V = \frac{{\rm d}V}{{\rm d}z}{\rm d}z\f$,
       *  \f$g\f$ is the growth factor, \f$j_\ell\f$ the
       *  spherical Bessel function,
       *  \f$r\f$ is the comoving distance,
       *  and \f$W_i\f$ is the window function.
       *
       *  We note that
       *
       *  \f$ f_{\rm sky} = \sqrt{\frac{\mathcal{C}_{\rm m}(\ell=0)}{4\pi}}\,, \f$
       *
       *  Physical units are forced.
       *
       *  @param cosmology the Cosmology object
       *
       *  @param Cl_mask power spectrum of the angular survey mask. The code assumes
       *  that Cl_mask[0] corresponds to \f$\ell=0\f$, and that the length of Cl_mask
       *  corresponds to the maximum \f$\ell\f$ considered minus 1
       *
       *  @param cosmo_param the cosmological parameters
       *  set in the modelling (useful only if a modelling
       *  is performed)
       *
       *  @param method_Pk the method used to compute
       *  the power spectrum
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
      SuperSampleCovariance (std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const std::vector<double> Cl_mask, const std::vector<std::string> cosmo_param, const std::string method_Pk="EisensteinHu", const double precision=10, const bool NL=false, const bool store_output=false);

      /**
       *  @brief default destructor
       */
      virtual ~SuperSampleCovariance () = default;

      ///@}
      
      /**
       *  @name Member functions to compute \f$S_{ij}\f$
       */
      ///@{
      
      /**
       * @brief set the top-hat window functions 
       * in the redshift bins. A top-hat
       * window function can be used when the redshift errors are
       * smaller than the width of the redshift bins.
       *
       * @param redshift_edges redshift bin edges
       *
       * @param delta_z redshift step used to construct the
       *  window function
       *
       */
      void set_topHat_window (const std::vector<double> redshift_edges, const double delta_z=0.001);
      
      /**
       * @brief set the Gaussian window functions 
       * in the redshift bins
       *
       *  @param W_mean vector of mean values for the 
       *  Gaussian window functions
       *
       *  @param W_std vector of standard deviation values
       *  for the Gaussian window functions
       *
       * @param delta_z redshift step used to construct the
       *  window function
       *
       */
      void set_gaussian_window (const std::vector<double> W_mean, const std::vector<double> W_std, const double delta_z=0.001);
      
      /**
       *  @brief get \f$S_{ij}\f$
       *
       *  @param parameter the parameters of interest in the covariance matrix
       *
       *  @return the covariance matrix
       */
      std::vector<std::vector<double>> operator () (std::vector<double> &parameter) const;
       
      ///@}
      
      /**
       *  @name Member functions to get the private members of the class
       */
      ///@{
      
      /**
       *  @brief Get the private member m_cosmology
       *
       *  @return the cosmological model
       */
      std::shared_ptr<cosmology::Cosmology> cosmology ()
      { return move(m_cosmology); }

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
