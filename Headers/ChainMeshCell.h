/********************************************************************
 *  Copyright (C) 2022 by Federico Marulli and Simone Sartori       *
 *  federico.marulli3@unibo.it simone.sartori5@studio.unibo.it      *
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
 *  @file Headers/ChainMeshCell.h
 *
 *  @brief The class ChainMeshCell 
 *
 *  This file defines the interface of the class ChainMeshCell, used to handle
 *  objects of type <EM> chain mesh cells </EM>
 *
 *  @authors Federico Marulli, Simone Sartori
 *
 *  @authors federico.marulli3@unibo.it, simone.sartori5@studio.unibo.it
 */

#ifndef __CHAINMESHCELL__
#define __CHAINMESHCELL__ 


// ===================================================================================================


namespace cbl {

  namespace catalogue {

    /**
     *  @class ChainMeshCell ChainMeshCell.h "Headers/ChainMeshCell.h"
     *
     *  @brief The class ChainMeshCell
     *
     *  This class is used to handle objects of type <EM> ChainMeshCell
     *  </EM>
     */
    class ChainMeshCell : public Object {

    protected :

    /// the particles in the cell
    std::vector<unsigned int> m_part;   

    /// the cells near the object, sorted for distance
    std::vector<std::vector<unsigned int>> m_nearCells; 

    public:
    
    /**
     *  @name Constructors/destructors
     */
    ///@{   

    /**
     *  @brief default constructor
     *  
     */
    ChainMeshCell() : Object(), m_part({}), m_nearCells({}) {}

    /**
     *  @brief constructor that uses comoving coordinates
     *
     *  @param ID the object ID
     *
     *  @param part the particles in the cell
     *
     *  @param nearCells the cells near the object, sorted for distance
     *  
     */
    ChainMeshCell(const int ID, const std::vector<unsigned int> part={}, std::vector<std::vector<unsigned int>> nearCells={}) 
      : Object(ID), m_part(part), m_nearCells(nearCells) {}

      /**
       * @brief function that allows copying private variables of the class 
       * when an object of class Catalogue is copied
       * 
       * @return a shared pointer to the Object
       *
       */
      std::shared_ptr<Object> getShared() {
        return std::make_shared<ChainMeshCell>(*this);
      }
      
      /**
       *  @brief default destructor
       */
      ~ChainMeshCell() = default;   
      ///@}
    
  
      /**
       *  @name Member functions used to get the private members 
       */
      ///@{
      
      /**
       *  @brief get the protected member \e m_part
       *  @return the particles in the cell
       */
      std::vector<unsigned int> part() const override
      { return m_part; }      

      /**
       *  @brief get the protected member \e m_nearCells
       *  @return the cells near the object, sorted for distance
       */
      std::vector<std::vector<unsigned int>> nearCells() const override
      { return m_nearCells; }     
      ///@}


      /**
       *  @name Member functions used to set the private members 
       */
      ///@{
      
      /**
       *  @brief set the protected member \e m_part
       *  @param part the particles in the cell
       */
      void set_part(const std::vector<unsigned int> part={}) override
      { m_part = part; }

      /**
       *  @brief set the protected member \e m_nearCells
       *  @param nearCells the cells near the object, sorted for distance
       */
      void set_nearCells(const std::vector<std::vector<unsigned int>> nearCells = {}) override
      { m_nearCells = nearCells; }      

      ///@}     
      /**
       *  @name Member functions used to check if the private members are set 
       */
      ///@{
      
      /**
       *  @brief check if the protected member \e m_part is set
       *
       *  @return true if the particles in the cell are set; false otherwise
       */
      bool isSet_part() 
      { return (cbl::isSet(m_part)) ? true : false; }

      /**
       *  @brief check if the protected member \e m_nearCells is set
       *  
       *  @return true if the the cells near the object are set; false otherwise
       */
      bool isSet_nearCells() 
      { return (cbl::isSet(m_nearCells)) ? true : false; }      
      ///@}
    
    };
  }
}

#endif
