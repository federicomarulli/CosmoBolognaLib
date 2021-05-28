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
 *******************************************************************/

/**
 *  @file Headers/LogNormal.h
 *
 *  @brief Implementation of the log-normal data structure
 *
 *  This file defines the interface of the class LogNormal
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __LOGNORMAL__
#define __LOGNORMAL__

#include "Catalogue.h"
#include <fftw3.h>

namespace cbl {
  
  /**
   *  @brief The namespace of the functions and classes used to
   *  construct <B> log-normal mocks </B> 
   *  
   *  The \e lognormal namespace contains all the functions and
   *  classes used to construct log-normal mocks
   */
  namespace lognormal {
    
    /**
     *  @class LogNormal LogNormal.h "Headers/LogNormal.h"
     *
     *  @brief The class LogNormal
     *
     *  This class is used to handle objects of type <EM> LogNormal
     *  </EM>
     */
    class LogNormal {
    
    protected:
 
      /// vector containing pointers to the log-normal realizations
      std::vector<std::shared_ptr<catalogue::Catalogue>> m_catalogue = {}; 
      
      /// the random catalogues used to construct the mask
      catalogue::Catalogue m_random;

      /// the assumed cosmological model
      cosmology::Cosmology m_cosmology;

      /// the mean total number of objects in the log-normal catalogues
      int m_nObjects; 

      /// the mean redshift the log-normal catalogues
      double m_redshift;  
    
      /// the bias of the log-normal density field catalogues
      double m_bias;

      /// the cell size in comoving scale
      double m_cell_size;
      
      /// true &rarr; real space; false &rarr; redshift space (only monopole distortions) 
      bool m_real;
    
      /// the method to compute the model power spectrum (i.e. the Boltzmann solver)
      std::string m_method_Pk;
    
      /// true &rarr; compute the non-linear power spectrum; false &rarr; compute the linear power spectrum
      bool m_NL;

    
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{
    
      /**
       *  @brief default constructor
       */
      LogNormal () = default;

      /**
       *  @brief default destructor
       */
      ~LogNormal () = default;

      /**
       *  @brief constructor 
       *
       *  @param random input random catalogue (should be much larger
       *  than the random catalogue used to measure the two-point
       *  correlation function)
       *
       *  @param cosmology the assumed cosmological model
       *
       *  @param nObjects mean total number of objects in the
       *  log-normal catalogues
       *
       *  @param redshift mean redshift the log-normal catalogues
       *
       *  @param bias bias of the log-normal density field catalogues
       *
       *  @param cell_size the cell size of the density field in
       *  comoving coordinates
       
       *  @param real true &rarr; real space; false &rarr; redshift
       *  space (only monopole distortions)
       *
       *  @param method_Pk the method to compute the model power
       *  spectrum (i.e. the Boltzmann solver)
       *
       *  @param NL true &rarr; compute the non-linear power spectrum;
       *  false &rarr; compute the linear power spectrum
       */
      LogNormal (const catalogue::Catalogue random, const cosmology::Cosmology cosmology, const int nObjects, const double redshift, const double bias, const double cell_size, const bool real=true, const std::string method_Pk="CAMB", const bool NL=false)
	: m_random(random), m_cosmology(cosmology), m_nObjects(nObjects), m_redshift(redshift), m_bias(bias), m_cell_size(cell_size), m_real(real), m_method_Pk(method_Pk), m_NL(NL) {}

      ///@}
      

      /**
       *  @name Functions to get the private members of the class
       */
      ///@{
      
      /**
       *  @brief get the private member LogNormal::m_LNCat[i]
       *
       *  @param i index of the log-normal realization
       *
       *  @return the i-th log-normal realization
       */
      std::shared_ptr<catalogue::Catalogue> catalogue (const size_t i);
      
      /**
       *  @brief get the private member LogNormal::m_nObjects
       *
       *  @return the mean total number of objects in the log-normal
       *  catalogues
       */
      int nObjects () const { return m_nObjects; }

      /**
       *  @brief get the private member LogNormal::m_redshift
       *
       *  @return the mean redshift the log-normal catalogues
       *  catalogues
       */
      int redshift () const { return m_redshift; }

      /**
       *  @brief get the private member LogNormal::m_bias
       *
       *  @return the bias of the log-normal density field catalogues
       */
      double bias () const { return m_bias; }

      /**
       *  @brief get the private member LogNormal::m_real
       *
       *  @return true &rarr; real space; false &rarr; redshift space
       *  (only monopole distortions)
       */
      bool real () const { return m_real; }

      /**
       *  @brief get the private member LogNormal::m_method_Pk
       *
       *  @return the method to compute the model power spectrum
       *  (i.e. the Boltzmann solver)
       */
      std::string method_Pk () const { return m_method_Pk; }

      /**
       *  @brief get the private member LogNormal::m_NL
       *
       *  @return true &rarr; compute the non-linear power spectrum;
       *  false &rarr; compute the linear power spectrum
       */
      bool NL () const { return m_NL; }

      ///@}
      

      /**
       *  @name Functions to set the private members of the class
       */
      ///@{
      
      /**
       *  @brief set the private member LogNormal::m_nObjects
       *
       *  @param nObjects the mean total number of objects in the
       *  log-normal catalogues
       */
      void set_nObjects (const int nObjects) { m_nObjects = nObjects; }

      /**
       *  @brief set the private member LogNormal::m_redshift
       *
       *  @param redshift the mean redshift the log-normal catalogues
       *  catalogues
       */
      void set_redshift (const double redshift) { m_redshift = redshift; }

      /**
       *  @brief set the private member LogNormal::m_bias
       *
       *  @param bias the bias of the log-normal density field
       *  catalogues
       */
      void set_bias (const double bias) { m_bias = bias; }

      /**
       *  @brief set the private member LogNormal::m_real
       *
       *  @param real true &rarr; real space; false &rarr; redshift
       *  space (only monopole distortions)
       */
      void set_real (const bool real) { m_real = real; }

      /**
       *  @brief set the private member LogNormal::m_method_Pk
       *
       *  @param method_Pk the method to compute the model power
       *  spectrum (i.e. the Boltzmann solver)
       */
      void set_method_Pk (const std::string method_Pk) { m_method_Pk = method_Pk; }

      /**
       *  @brief set the private member LogNormal::m_NL
       *
       *  @param NL true &rarr; compute the non-linear power spectrum;
       *  false &rarr; compute the linear power spectrum
       */
      void set_NL (const bool NL) { m_NL = NL; }

      ///@}
      

      /**
       *  @name Functions to generate the log-normal mock catalogues
       */
      ///@{
      
      /**
       *  @brief generate the log-normal mock catalogues
       *
       *  @param n_lognormal_mocks number of log-normal mock
       *  catalogues to be constructed
       *
       *  @param output_dir the output directory
       * 
       *  @param filename the prefix of the ouput file containing the
       *  LogNormal realizations
       *
       *  @param start the starting index of the mock to be created
       *
       *  @param seed the seed for random number generation
       */
      void generate (const int n_lognormal_mocks, const std::string output_dir, const std::string filename="lognormal", const int start=1, const int seed=3213);

      ///@}
      
    };

  }

}

#endif
