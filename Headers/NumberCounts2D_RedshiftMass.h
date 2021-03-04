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
 *  @file Headers/NumberCounts2D_RedshiftMass.h
 *
 *  @brief The class NumberCounts2D_RedshiftMass
 *
 *  This file defines the interface of the class NumberCounts2D_RedshiftMass,
 *  used to measure the mass-redshift number counts
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __NCOUNTS2DMR__
#define __NCOUNTS2DMR__


#include "NumberCounts2D.h"


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
       *  @class NumberCounts2D_RedshiftMass NumberCounts2D_RedshiftMass.h
       *  "Headers/NumberCounts2D_RedshiftMass.h"
       *
       *  @brief The class NumberCounts2D_Redshift
       *
       *  This is the base class used to measure the 
       *  mass-redshift number counts
       */
      class NumberCounts2D_RedshiftMass : public NumberCounts2D {

	public:

	  /**
	   *  @name Constructors/destructors
	   */
	  ///@{

	  /**
	   *  @brief default constructor
	   */
	  NumberCounts2D_RedshiftMass () {}

	  /**
	   *  @brief default destructor
	   */
	  virtual ~NumberCounts2D_RedshiftMass () = default;

	  /**
	   * @brief constructor
	   * 
	   * @param data object of class Catalogue 
	   *
	   * @param nbins1 the number of bins for the first variable
	   *
	   * @param nbins2 the number of bins for the second variable
	   *
	   * @param minVar1 minimum range  for the first variable
	   * 
	   * @param maxVar1 maximmum range for the first variable
	   * 
	   * @param minVar2 minimum range for the second variable 
	   * 
	   * @param maxVar2 maximmum range for the second variable
	   * 
	   * @param shift1 bin shift for the first variable
	   * 
	   * @param shift2 bin shift for the second variable
	   * 
	   * @param hist_type the type of histogram
	   *
	   * @param fact factor used to normalized the distribution
	   *
	   * @param bin_type1 the var1 bin type
	   *
	   * @param bin_type2 the var2 bin type
	   */
	  NumberCounts2D_RedshiftMass (const catalogue::Catalogue data, const size_t nbins1, const size_t nbins2, const double minVar1=par::defaultDouble, const double maxVar1=par::defaultDouble, const double minVar2=par::defaultDouble, const double maxVar2=par::defaultDouble, const double shift1=0.5, const double shift2=0.5, const glob::HistogramType hist_type=glob::HistogramType::_N_V_, const double fact = 1., const BinType bin_type1=BinType::_linear_, const BinType bin_type2=BinType::_logarithmic_) : NumberCounts2D (catalogue::Var::_Redshift_, bin_type1, catalogue::Var::_Mass_, bin_type2, data, nbins1, nbins2, minVar1, maxVar1, minVar2, maxVar2, shift1, shift2, hist_type, fact) {} 
	  
	  /**
	   * @brief constructor
	   * 
	   * @param data object of class Catalogue
	   *
	   * @param vec_edges1 the bin edges for var1
	   *
	   * @param vec_edges2 the bin edges for var2
	   * 
	   * @param hist_type the type of histogram
	   *
	   * @param fact factor used to normalized the distribution
	   *
	   */
	  NumberCounts2D_RedshiftMass (const catalogue::Catalogue data, const std::vector<double> vec_edges1, const std::vector<double> vec_edges2, const glob::HistogramType hist_type=glob::HistogramType::_N_V_, const double fact = 1.) : NumberCounts2D (catalogue::Var::_Redshift_, catalogue::Var::_Mass_, vec_edges1, vec_edges2, data, hist_type, fact) {} 

	  ///@}

      }; 
    }
  }
}

#endif
