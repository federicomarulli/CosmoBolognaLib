/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/Models/ModelCosmology.h
 *
 *  @brief The class ModelCosmology
 *
 *  This file defines the interface of the class ModelCosmology
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELCOSM__
#define __MODELCOSM__

#include "Model.h"

namespace cosmobl {
  
  namespace modelling {

    /**
     *  @class ModelCosmology ModelCosmology.h "Headers/Lib/ModelCosmology.h"
     *
     *  @brief The class ModelCosmology
     *
     *  This class is used to model cosmological parameters
     */
    class ModelCosmology : public Model
    {

    protected:
    
      /// number of model parameters
      vector<cosmobl::CosmoPar> m_cosmological_parameters;

    
    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{ 

      /**
       *  @brief default constructor
       *
       *  @return object of class Model
       */
      ModelCosmology () {}

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      virtual ~Model() {}

      ///@}
    };

  }
  
}

#endif
