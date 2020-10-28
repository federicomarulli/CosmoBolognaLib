/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/Modelling_NumberCounts2D.h
 *
 *  @brief The class Modelling_NumberCounts2D
 *
 *  This file defines the interface of the class Modelling_NumberCounts2D, 
 *  used to  model 2D number counts
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGNC2D__
#define __MODELLINGNC2D__

#include "Modelling2D.h"
#include "Modelling_NumberCounts.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    /**
     *  @brief The namespace of the <B> number counts
     *  modelling </B>
     *  
     *  The \e modelling::numbercounts namespace contains all the functions
     *  and classes to model number counts
     */
    namespace numbercounts {
    
      /**
       *  @class Modelling_NumberCounts2D
       *  Modelling_NumberCounts2D.h
       *  "Headers/Modelling_NumberCounts2D.h"
       *
       *  @brief The class Modelling_NumberCounts2D
       *
       *  This file defines the interface of the base class
       *  Modelling_NumberCounts2D, used for modelling 
       *  2D number counts measurements
       *
       */
      class Modelling_NumberCounts2D : public Modelling_NumberCounts, public Modelling2D
      {

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  _NumberCounts2D
	 */
	Modelling_NumberCounts2D () = default;
	
	/**
	 *  @brief constuctor
	 *  @param nc the number counts to model
	 *  _NumberCounts2D
	 */
	Modelling_NumberCounts2D (const std::shared_ptr<cbl::measure::numbercounts::NumberCounts> nc) 
	  : Modelling_NumberCounts (nc) {m_data = nc->dataset();}
	
	/**
	 *  @brief constuctor
	 *  @param dataset the number counts dataset
	 *  @param hist_type the histogram type
	 *  @param fact the normalization factor
	 *
	 *  _NumberCounts2D
	 */
	Modelling_NumberCounts2D (const std::shared_ptr<cbl::data::Data> dataset, glob::HistogramType hist_type, double fact)
	  : Modelling_NumberCounts (hist_type, fact) {m_data = dataset;}
	
	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~Modelling_NumberCounts2D () = default;

	///@}

      };
    }
  }
}

#endif
