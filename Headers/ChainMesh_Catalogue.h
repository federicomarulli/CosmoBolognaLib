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
 ********************************************************************/

/**
 *  @file Headers/ChainMesh_Catalogue.h
 *
 *  @brief Implementation of the chain-mesh data structure
 *
 *  This file contains the implementation of the chain-mesh 
 *  for Catalogue object
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "Catalogue.h"


#ifndef __CHAINMESHCAT__
#define __CHAINMESHCAT__


namespace cbl {

  namespace chainmesh {
    
    /**
     *  @class ChainMesh_Catalogue ChainMesh_Catalogue.h
     * "Headers/ChainMesh_Catalogue.h"
     *
     *  @brief The class ChainMesh_Catalogue
     *
     *  This class is used to handle objects of type <EM>
     *  ChainMesh_catalogue </EM>
     */
    class ChainMesh_Catalogue : public ChainMesh
    {
    
    private:
    
      /// pointer to catalogue used for the chain-mesh
      std::shared_ptr<catalogue::Catalogue> m_catalogue;

    public:
      /**
       *  @brief default constructor
       *  _Catalogue
       */
      ChainMesh_Catalogue () = default;

      /**
       *  @overload
       *  @brief constructor 
       *  @param cell_size the cell size
       *  @param cat pointer to an object of class Catalogue
       *  @param rmax the maximum separation
       *  @param rmin the minimum separation
       *  _Catalogue
       */
      ChainMesh_Catalogue (const double cell_size, std::shared_ptr<catalogue::Catalogue> cat, const double rmax, const double rmin=-1.) : ChainMesh(cell_size, 3)
	{
	  set_par(cell_size, cat, rmax, rmin); 
	}

      /**
       *  @brief default destructor
       */
      ~ChainMesh_Catalogue () = default;

      /**
       *  @brief function that set parameters for the chain-mesh 
       *  @param cell_size storing the cell size
       *  @param cat pointer to an object Catalogue
       *  @param rmax the maximum separation
       *  @param rmin the minimum separation
       */
      void set_par (const double cell_size, std::shared_ptr<catalogue::Catalogue> cat, const double rmax, const double rmin=-1.);

      /**
       *  @brief order the catalogue according to the input vector
       *  @param order vector used to order the catalogue
       */
      void get_order (std::vector<int> &order) const;

      /**
       *  @brief get list of objects close to the input
       *  @param object the center object
       *  @param ii -1 &rarr; takes all close objects; ii > -1 &rarr; takes
       *  only objects of index > ii
       *  @return vector of close objects
       */
      std::vector<std::shared_ptr<catalogue::Object>> object_list (std::shared_ptr<catalogue::Object> object, const int ii=-1);

      /**
       *  @brief get the internal variable m_catalogue
       *  @return the internal variable m_catalogue
       */
      std::shared_ptr<catalogue::Catalogue> catalogue () const { return m_catalogue; }
      
    };
  }
}

#endif
