/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli, Michele Moresco         *
 *  and Alfonso Veropalumbo                                         *
 *                                                                  *
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
 *  @file Headers/Lib/ThreePointCorrelation_comoving_connected.h
 *
 *  @brief The class ThreePointCorrelation_comoving_connected
 *
 *  This file defines the interface of the class
 *  ThreePointCorrelation_comoving_connected, used to measure the
 *  connected three-point correlation function in comoving coordinates
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#ifndef __THREEPOINTCOMCON__
#define __THREEPOINTCOMCON__ 


#include "ThreePointCorrelation.h"


// ===================================================================================================


namespace cosmobl {
  
  namespace threept {
    
    /**
     *  @class ThreePointCorrelation_comoving_connected ThreePointCorrelation_comoving_connected.h
     * "Headers/Lib/ThreePointCorrelation_comoving_connected.h"
     *
     *  @brief The class ThreePointCorrelation_comoving_connected
     *
     *  This is the base class used to measure the connected three-point
     *  correlation function in comoving coordinates
     */
    class ThreePointCorrelation_comoving_connected : public ThreePointCorrelation {

    protected :
    
      /**
       *  @name Three-point correlation function data
       */
      ///@{
    
      /// scale bins
      vector<double> m_scale;

      /// binned connected three-point correlation function
      vector<double> m_zeta;

      /// error on the binned connected three-point correlation function
      vector<double> m_error;
    
      ///@}

    
    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class
       *  ThreePointCorrelation_comoving_connected
       */
      ThreePointCorrelation_comoving_connected () {}

      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param tripletType the triplet type; it can be:
       *  _comoving_theta_, _comoving_side_
       *  @param side_s the size of r<SUB>12</SUB>
       *  @param side_u the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
       *  @param perc_increase the ratio
       *  &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
       *  @param nbins number of bins
       *  @return object of class
       *  ThreePointCorrelation_comoving_connected
       */
      ThreePointCorrelation_comoving_connected (const catalogue::Catalogue data, const catalogue::Catalogue random, const triplets::TripletType tripletType, const double side_s, const double side_u, const double perc_increase, const int nbins)
	: ThreePointCorrelation(data, random) { set_parameters(tripletType, side_s, side_u, perc_increase, nbins); }

      /**
       *  @brief default destructor
       *  @return none
       */
      ~ThreePointCorrelation_comoving_connected () {}

      ///@}

      
      /**
       *  @name Member functions to set the binning parameters
       */
      ///@{
      
      /**
       *  @brief set the binning parameters
       *  @param tripletType the triplet type; it can be:
       *  _comoving_theta_, _comoving_side_
       *  @param side_s the size of r<SUB>12</SUB>
       *  @param side_u the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
       *  @param perc_increase the ratio
       *  &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
       *  @param nbins number of bins
       *  @return none
       */
      void set_parameters (const triplets::TripletType tripletType, const double side_s, const double side_u, const double perc_increase, const int nbins);

      ///@}

      
      /**
       *  @name Member functions to get protected parameters
       */
      ///@{

      /**
       *  @brief get the protected member
       *  ThreePointCorrelation_comoving_connected::m_scale
       *  @return the scale bins
       */
      vector<double> scale () const override { return m_scale; }

      /**
       *  @brief get the protected member
       *  ThreePointCorrelation_comoving_connected::m_zeta
       *  @return the binned connected three-point correlation
       *  function
       */
      vector<double> zeta () const override { return m_zeta; }

      /**
       *  @brief get the protected member
       *  ThreePointCorrelation_comoving_connected::m_error
       *  @return the error on the connected three-point
       *  correlation function
       */
      vector<double> error () const override { return m_error; }

      ///@}


      /**
       *  @name Member functions to measure the three-point correlation function
       */
      ///@{

      /**
       * @brief method to measure the three-point correlation function
       *
       * @param dir_output_triplets name of the output directory used to
       * store the number of triplets
       * 
       * @param dir_input_triplets name of the input directories
       * containing the number of triplets
       *
       * @param count_ddd 1 &rarr; count the data-data-data
       * triplets; 0 &rarr; read the data-data-data triplets
       * from a file
       *
       * @param count_rrr 1 &rarr; count the random-random-random
       * triplets; 0 &rarr; read the random-random-random triplets
       * from a file
       *
       * @param count_ddr 1 &rarr; count the data-data-random
       * triplets; 0 &rarr; read the data-data-random triplets
       * from a file
       *
       * @param count_drr 1 &rarr; count the data-random-random
       * triplets; 0 &rarr; read the data-random-random triplets
       * from a file
       *
       * @param tcount 1 &rarr; activate the CPU time counter; 0
       * &rarr; no time counter
       *
       * @return none
       */
      void measure (const string dir_output_triplets, const vector<string> dir_input_triplets={}, const int count_ddd=1, const int count_rrr=1, const int count_ddr=1, const int count_drr=1, const bool tcount=0) override;
    
      ///@}

    
      /**
       *  @name Input/Output methods
       */
      ///@{
    
      /**
       *  @brief write the monopole of the two-point correlation
       *  function
       *  @param dir output directory
       *  @param file output file
       *  @return none
       */
      void write (const string dir, const string file) const override;
    
      ///@}
      
    };
  }
}

#endif
