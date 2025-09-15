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
 *  @file Headers/MapCalculator.h
 *
 *  @brief The class MapCalculator
 *
 *  This file defines the interface of the class MapCalculator, used to 
 *  build the sky map for angular density fluctuations (ADF), angular 
 *  redshift fluctuations (ARF) and shot noise random realisations (SN)
 *
 *  @author Matteo Santini
 *
 *  @author matteo.santini7@studio.unibo.it
 */

#ifndef __MAP__
#define __MAP__


#include "Measure.h"
#include "Catalogue.h"
#include "Func.h"
#include <healpix_base.h>


// ===================================================================================================


namespace cbl {

  namespace measure {
   
    /**
     *  @brief Class for <B> HEALPix map generation </B>
     *  
     *  The \e measure::MapCalculator namespace contains all the functions
     *  and classes for building HEALPix maps for ADF and ARF
     */
    class MapCalculator : public Measure { 

    private:
        
      /**
      *  @param m_ns: the resolution nside of Healpix (must be a power of 2) 
      */
      int m_ns;

      /**
      *  @param m_scheme: scheme for the HEALPix map. RING (default) or NESTED 
      */
      Healpix_Ordering_Scheme m_scheme = RING; //RING is default

      /**
      *  @param m_nsg: number of sigmas for source selection (speeds up the computation)
      */
      int m_nsg;

      /**
      *  @param m_zobs: central redshifts of the redshift shells
      */
      std::vector<double> m_zobs;

      /**
      *  @param m_sgzar: width of the gaussian shells
      */
      std::vector<double> m_sgzar;

      /**
      *  @param m_seed: seed for random number generation
      */
      int m_seed = 123;
      	
    protected:
      /**
      *  @param m_data: points to the catalogue of data
      */
      std::shared_ptr<catalogue::Catalogue> m_data;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      MapCalculator () = default;

      /**
       *  @brief constructor with specified parameters.
       *
       *  @param ns_value resolution nside of Healpix (must be a power of 2).
       * 
       *  @param nsg_value number of sigmas for source selection.
       * 
       *  @param zobs central redshifts of the redshift shells.
       * 
       *  @param sgzar width of the Gaussian shells.
       * 
       *  @param data catalogue of data.
       */
      MapCalculator (int ns_value, int nsg_value, std::vector<double>& zobs, std::vector<double>& sgzar, const catalogue::Catalogue data) 
      : m_ns(ns_value), m_nsg(nsg_value), m_zobs(zobs), m_sgzar(sgzar) , m_data(std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(data)))) {};
        
      /**
       *  @brief constructor with specified parameters.
       *
       *  @param ns_value resolution nside of Healpix (must be a power of 2).
       * 
       *  @param scheme scheme for the HEALPix map (RING or NESTED).
       * 
       *  @param nsg_value number of sigmas for source selection.
       * 
       *  @param zobs central redshifts of the redshift shells.
       * 
       *  @param sgzar width of the Gaussian shells.
       * 
       *  @param data catalogue of data.
       */
      MapCalculator (int ns_value, Healpix_Ordering_Scheme scheme, int nsg_value, std::vector<double>& zobs, std::vector<double>& sgzar, const catalogue::Catalogue data) 
      : m_ns(ns_value), m_scheme(scheme), m_nsg(nsg_value), m_zobs(zobs), m_sgzar(sgzar) , m_data(std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(data)))) {};        
        
      /**
       *  @brief constructor with specified parameters.
       *
       *  @param ns_value resolution nside of Healpix (must be a power of 2).
       * 
       *  @param scheme scheme for the HEALPix map (RING or NESTED).
       * 
       *  @param nsg_value number of sigmas for source selection.
       * 
       *  @param zobs central redshifts of the redshift shells.
       * 
       *  @param sgzar width of the Gaussian shells.
       * 
       *  @param seed seed for random number generation.
       * 
       *  @param data catalogue of data.
       */
      MapCalculator (int ns_value, Healpix_Ordering_Scheme scheme, int nsg_value, std::vector<double>& zobs, std::vector<double>& sgzar, const int seed, const catalogue::Catalogue data) 
      : m_ns(ns_value), m_scheme(scheme), m_nsg(nsg_value), m_zobs(zobs), m_sgzar(sgzar) , m_seed(seed), m_data(std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(data)))) {};

      //@}

      /**
       *  @name Functions to get the private members of the class
       */
      ///@{

      /**
       *  @brief get the value of the private variable m_ns 
       *  @return the value of m_ns
       */
      int getNside () const {return m_ns;};

      /**
       *  @brief get the value of the private variable m_scheme 
       *  @return the value of m_scheme
       */
      Healpix_Ordering_Scheme getScheme () {return m_scheme;};

      /**
       *  @brief get the value of the private variable m_nsg 
       *  @return the value of m_nsg
       */
      int getNsigma () const {return m_nsg;};

      /**
       *  @brief get the values of the private variable m_zobs
       *  @return the value of m_zobs
       */
      std::vector<double> getZobs () const {return m_zobs;};

      /**
       *  @brief get the values of the private variable m_sgzar
       *  @return the value of m_sgzar
       */
      std::vector<double> getWidths () const {return m_sgzar;};

      /**
       *  @brief get the value of the private variable m_seed 
       *  @return the value of m_seed
       */
      int getSeed () const {return m_seed;};

      /**
       *  @brief get the values of the private variable m_data 
       *  @return pointer to the member m_data of the class Catalogue
       */
      std::shared_ptr<catalogue::Catalogue> getCatalogue () const {return m_data;};
        
      ///@}
      
      /**
       *  @name Functions to set the private members of the class
       */
      ///@{

      /**
       *  @brief set the value of the private variable m_ns 
       *  @param ns resolution nside of Healpix (must be a power of 2).
       */
      void setNside (const int ns) {m_ns = ns;};

      /**
       *  @brief set the value of the private variable m_scheme 
       *  @param scheme scheme for the HEALPix map (RING or NESTED).
       */
      void setScheme (Healpix_Ordering_Scheme scheme) {m_scheme = scheme;};

      /**
       *  @brief set the value of the private variable m_nsg 
       *  @param nsg number of sigmas for source selection.
       */
      void setNsigma (const int nsg) {m_nsg=nsg;};

      /**
       *  @brief set the values of the private variable m_zobs
       *  @param zobs central redshifts of the redshift shells.
       */
      void setZobs (const std::vector<double> zobs) {m_zobs=zobs;};

      /**
       *  @brief set the values of the private variable m_sgzar
       *  @param sgzar width of the Gaussian shells.
       */
      void setWidths (const std::vector<double> sgzar) {m_sgzar=sgzar;};

      /**
      *  @brief set the value of the private variable m_seed 
      *  @param seed seed for random number generation.
      */
      void setSeed (const int seed) {m_seed=seed;};
        
      /**
       *  @brief set the values of the private variable m_data 
       *  @param data catalogue of data.
       */
      void setCatalogue (const catalogue::Catalogue data)
       {m_data = std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(data)));};

      //@}

      /**
       *  @name Functions to retrieve HEALPix functions
       */
      ///@{

      /**
       *  @brief performs the conversion from angles to pixels
       * to build the HEALPix map. The angles are treated in radians units
       * and they refer to spherical polar coordinates theta phi with domain 
       * [0, pi[ and [0, 2pi[ 
       *  @return integer array of pixels
       */
      std::vector<int> HealpixMapping () const;

      //@}

      /**
       *  @name Functions to generate HEALPix maps
       */
      ///@{

      /**
       *  @brief builds the map for the angular density fluctuations, given the input mask
       *  @param mskin: Healpix sky mask
       *  @return 3D tensor containing the values associated to the map's pixels
       */
      std::vector<std::vector<std::vector<double>>> ADF (std::vector<double> const& mskin) const;

      /**
       *  @brief builds the map for the angular density fluctuations and the angular redshift fluctuations,
       *  given the input mask
       *  @param mskin: Healpix sky mask
       *  @return standard tuple containing the two maps
       */
      std::tuple<
      std::vector<std::vector<std::vector<double>>> ,
      std::vector<std::vector<std::vector<double>>>
      > ADF_ARF (std::vector<double> const& mskin) const;

      /**
       *  @brief builds the map for the redshifted angular density fluctuations, given the input mask
       *  @param mskin: Healpix sky mask
       *  @param zH: Hubble drift redshifts
       *  @return 3D tensor containing the values associated to the map's pixels
       */
      std::vector<std::vector<std::vector<double>>> ADFr (std::vector<double> const& mskin, std::vector<double> const& zH) const;

      /**
       *  @brief builds the map for the redshifted angular density fluctuations and the redshifted 
       * angular redshift fluctuations, given the input mask
       *  @param mskin: Healpix sky mask
       *  @param zH: Hubble drift redshifts
       *  @return standard tuple containing the two maps
       */
      std::tuple<
      std::vector<std::vector<std::vector<double>>> ,
      std::vector<std::vector<std::vector<double>>>
      > ADFr_ARFr (std::vector<double> const& mskin, std::vector<double> const& zH) const;

      /**
       *  @brief builds the map for the shot noise, given the input mask
       *  @param mskin: Healpix sky mask
       *  @param nran: number of mocks computing purely random ADF/ARF mocks
       *  @return 5D tensor containing 2 * nran maps
       */
      std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> SN (std::vector<double> const& mskin, const int nran) const;
       
      /**
      *  @brief prints a map on an output file
      *  @param map: your map as a 3D tensor
      *  @param filename: the name of the output file 
      */
      void print_map (const std::vector<std::vector<std::vector<double>>>& map, const std::string& filename) const;

      ///@}

    };

  }

}

#endif