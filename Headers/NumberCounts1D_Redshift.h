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
 *  @file Headers/NumberCounts1D_Redshift.h
 *
 *  @brief The class NumberCounts1D_Redshift
 *
 *  This file defines the interface of the class
 *  NumberCounts1D_Redshift, used to measure the redshift number
 *  counts
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __NCOUNTS1DR__
#define __NCOUNTS1DR__


#include "NumberCounts1D.h"


// ===================================================================================================


namespace cbl {

  namespace measure {

    /**
     *  @brief The namespace of the <B> number counts
     *  </B>
     *  
     *  The \e measure::numbercounts namespace contains all the functions and
     *  classes to measure the number counts
     */
    namespace numbercounts {

      /**
       *  @class NumberCounts1D_Redshift NumberCounts1D_Redshift.h
       *  "Headers/NumberCounts1D_Redshift.h"
       *
       *  @brief The class NumberCounts1D_Redshift
       *
       *  This is the base class used to measure the 
       *  redshift number counts
       *
       */
      class NumberCounts1D_Redshift : public NumberCounts1D {

	public:

	  /**
	   *  @name Constructors/destructors
	   */
	  ///@{

	  /**
	   *  @brief default constructor
	   *
	   *  @return object of class NumberCounts1D_Redshift
	   */
	  NumberCounts1D_Redshift () {}

	  /**
	   *  @brief default destructor
	   *  @return none
	   */
	  virtual ~NumberCounts1D_Redshift () = default;

	  /**
	   *  @brief constructor
	   *
	   *  @param data object of class Catalogue 
	   *
	   *  @param nbins the number of bins
	   *
	   *  @param minVar minimum range 
	   *  
	   *  @param maxVar maximmum range
	   *  
	   *  @param shift the shift of the bin
	   *  
	   *  @param hist_type the type of histogram
	   *
	   *  @param  fact factor used to normalized the distribution
	   *
	   *  @return object of class NumberCounts1D_Redshift
	   */
	  NumberCounts1D_Redshift (const catalogue::Catalogue data, const size_t nbins, const double minVar=par::defaultDouble, const double maxVar=par::defaultDouble, const double shift = 0.5, const glob::HistogramType hist_type=glob::HistogramType::_N_V_, const double fact = 1.) : NumberCounts1D(catalogue::Var::_Redshift_, BinType::_linear_, data, nbins, minVar, maxVar, shift, hist_type, fact) {}

	  ///@}

      }; 
    }
  }
}

#endif
