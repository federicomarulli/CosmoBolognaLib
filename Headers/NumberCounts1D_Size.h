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
 *  @file Headers/NumberCounts1D_Size.h
 *
 *  @brief The class NumberCounts1D_Size
 *
 *  This file defines the interface of the class NumberCounts1D_Size,
 *  used to measure the mass number counts
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __NCOUNTS1DS__
#define __NCOUNTS1DS__


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
       *  @class NumberCounts1D_Size NumberCounts1D_Size.h
       *  "Headers/NumberCounts1D_Size.h"
       *
       *  @brief The class NumberCounts1D_Size
       *
       *  This is the base class used to measure the 
       *  mass number counts
       *
       */
      class NumberCounts1D_Size : public NumberCounts1D {

	public:

	  /**
	   *  @name Constructors/destructors
	   */
	  ///@{

	  /**
	   *  @brief default constructor
	   *
	   *  1D_Size
	   */
	  NumberCounts1D_Size () {}

	  /**
	   *  @brief default destructor
	   *  
	   */
	  virtual ~NumberCounts1D_Size () = default;


	  /**
	   *  @brief constructor
	   *
	   *  @param data object of class Catalogue 
	   *
	   *  @param nbins the number of bins
	   *
	   *  @param minVar minimum range 
	   *  
	   *  @param maxVar maximum range
	   *  
	   *  @param shift the shift of the bin
	   *  
	   *  @param hist_type the type of histogram
	   *
	   *  @param fact factor used to normalized the distribution
	   *
	   *  1D_Size
	   */
	  NumberCounts1D_Size (const catalogue::Catalogue data, const size_t nbins, const double minVar=par::defaultDouble, const double maxVar=par::defaultDouble, const double shift = 0.5, const glob::HistogramType hist_type=glob::HistogramType::_dn_dlnV_, const double fact = 1.) : NumberCounts1D(catalogue::Var::_Radius_, BinType::_logarithmic_, data, nbins, minVar, maxVar, shift, hist_type, fact) {} 

	  ///@}

      }; 
    }
  }
}

#endif
