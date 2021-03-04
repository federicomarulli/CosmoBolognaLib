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
 *  @file Headers/NumberCounts1D_Mass.h
 *
 *  @brief The class NumberCounts1D_Mass
 *
 *  This file defines the interface of the class NumberCounts1D_Mass,
 *  used to measure the mass number counts
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __NCOUNTS1DM__
#define __NCOUNTS1DM__


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
       *  @class NumberCounts1D_Mass NumberCounts1D_Mass.h
       *  "Headers/NumberCounts1D_Mass.h"
       *
       *  @brief The class NumberCounts1D_Mass
       *
       *  This is the base class used to measure the 
       *  mass number counts
       *
       */
      class NumberCounts1D_Mass : public NumberCounts1D {

	public:

	  /**
	   *  @name Constructors/destructors
	   */
	  ///@{

	  /**
	   *  @brief default constructor
	   *
	   *  1D_Mass
	   */
	  NumberCounts1D_Mass () {}

	  /**
	   *  @brief default destructor
	   *  
	   */
	  virtual ~NumberCounts1D_Mass () = default;


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
	   *  @param bin_type bin type
	   *
	   *
	   *  1D_Mass
	   */
	  NumberCounts1D_Mass (const catalogue::Catalogue data, const size_t nbins, const double minVar=par::defaultDouble, const double maxVar=par::defaultDouble, const double shift = 0.5, const glob::HistogramType hist_type=glob::HistogramType::_N_V_, const double fact = 1., const BinType bin_type=BinType::_logarithmic_) : NumberCounts1D(catalogue::Var::_Mass_, bin_type, data, nbins, minVar, maxVar, shift, hist_type, fact) {}
	  
	  /**
	   *  @brief constructor
	   *
	   *  @param data object of class Catalogue 
	   *
	   *  @param vec_edges bin edges
	   *
	   *  @param hist_type the type of histogram
	   *
	   *  @param fact factor used to normalized the distribution
	   *
	   *  1D_Mass
	   */
	  NumberCounts1D_Mass (const catalogue::Catalogue data, const std::vector<double> vec_edges, const glob::HistogramType hist_type=glob::HistogramType::_N_V_, const double fact = 1.) : NumberCounts1D(catalogue::Var::_Mass_, vec_edges, data, hist_type, fact) {} 

	  ///@}

      }; 
    }
  }
}

#endif
