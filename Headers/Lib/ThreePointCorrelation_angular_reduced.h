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
 *  @file Headers/Lib/ThreePointCorrelation_angular_reduced.h
 *
 *  @brief The class ThreePointCorrelation_angular_reduced
 *
 *  This file defines the interface of the class
 *  ThreePointCorrelation_angular_reduced, used to measure the
 *  connected three-point correlation function in angular coordinates
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#ifndef __THREEPOINTANGRED__
#define __THREEPOINTANGRED__ 


#include "ThreePointCorrelation_angular_connected.h"


// ===================================================================================================


namespace cosmobl {
  
  namespace threept {
    
    /**
     *  @class ThreePointCorrelation_angular_reduced
     *  ThreePointCorrelation_angular_reduced.h
     *  "Headers/Lib/ThreePointCorrelation_angular_reduced.h"
     *
     *  @brief The class ThreePointCorrelation_angular_reduced
     *
     *  This is the base class used to measure the connected three-point
     *  correlation function in angular coordinates
     *
     *  @warning This class has not been fully implemented yet
     */
    class ThreePointCorrelation_angular_reduced : public ThreePointCorrelation_angular_connected {

    protected :
    
      /**
       *  @name Three-point correlation function data
       */
      ///@{
    
      /// angular bins
      vector<double> m_theta;

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
       *  ThreePointCorrelation_angular_reduced
       *
       *  @warning This method has not been implemented yet
       */
      ThreePointCorrelation_angular_reduced () { ErrorMsg("Work in progress..."); }

      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param side_s the size of r<SUB>12</SUB>
       *  @param side_u the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
       *  @param perc_increase the ratio
       *  &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
       *  @param nbins number of bins
       *  @return object of class ThreePointCorrelation_angular_reduced
       *
       *  @warning This method has not been implemented yet
       */
      ThreePointCorrelation_angular_reduced (const catalogue::Catalogue data, const catalogue::Catalogue random, const double side_s, const double side_u, const double perc_increase, const int nbins)
	: ThreePointCorrelation_angular_connected(data, random, side_s, side_u, perc_increase, nbins) { ErrorMsg("Work in progress..."); }

      /**
       *  @brief default destructor
       *  @return none
       */
      ~ThreePointCorrelation_angular_reduced () {}

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
       * @param dir_output_2pt name of the output directory used to
       * store the two-point correlation functions
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
       *
       * @warning This method has not been implemented yet 
       */
      void measure (const string dir_output_triplets, const string dir_output_2pt, const vector<string> dir_input_triplets={}, const int count_ddd=1, const int count_rrr=1, const int count_ddr=1, const int count_drr=1, const bool tcount=0) override;
    
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
       *  @param connected 0 &rarr; write the reducted 3pt correlation
       *  function; 1 &rarr; write both the reduced and connected 3pt
       *  correlation function
       *  @return none
       *
       *  @warning This method has not been implemented yet
       */
      void write (const string dir, const string file, const bool connected) const override;
    
      ///@}
      
    };
  }
}

#endif
