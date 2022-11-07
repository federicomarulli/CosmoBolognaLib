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
 *  @file Headers/CatalogueChainMesh.h
 *
 *  @brief The class catalogue::CatalogueChainMesh 
 *
 *  This file defines the interface of the class CatalogueChainMesh, used to handle
 *  objects of type <EM> CatalogueChainMesh </EM>
 *
 *  @author Simone Sartori, Federico Marulli
 *
 *  @author simone.sartori5@studio.unibo.it,  federico.marulli3@unibo.it
 */

#include "Catalogue.h"

#ifndef __CATALOGUECHAINMESH__
#define __CATALOGUECHAINMESH__

namespace cbl {

  namespace catalogue {
    
    /**
     *  @class Halo CatalogueChainMesh.h "Headers/CatalogueChainMesh.h"
     *
     *  @brief The class CatalogueChainMesh
     *
     *  This class is used to handle objects of type <EM> CatalogueChainMesh
     *  </EM>
     */
    class CatalogueChainMesh : public Catalogue {

    private :

      // size of cells in CatalogueChainMesh
      double m_cellsize;

      /// limits of the catalogue
      std::vector<std::vector<double>> m_lim;

      /// pointer to catalogue of particles in the chain-mesh
      std::shared_ptr<cbl::catalogue::Catalogue> m_part_catalogue;

      /// pointer to a second catalogue of particles, counterpart of the particle of the first catalogue;
      std::shared_ptr<cbl::catalogue::Catalogue> m_part_catalogue2={};   

      /// number of cells for side
      unsigned int m_dimension;

    public :
      /**
       *  @brief default constructor
       *  
       */
      CatalogueChainMesh () = default;

      /**
       *  @brief Copy constructor for the CatalogueChainMesh class. The copy constructor of 
       *  the parent class Catalogue is called
			 *  @param obj catalogue to be copied
       */
      CatalogueChainMesh (const CatalogueChainMesh& obj) : Catalogue (obj) {
        m_cellsize = obj.m_cellsize;
        m_lim = obj.m_lim;
        m_part_catalogue = obj.m_part_catalogue;
        m_part_catalogue2 = obj.m_part_catalogue2;
        m_dimension = obj.m_dimension;
      };

      /**
       *  @brief Move constructor for the CatalogueChainMesh class. The copy constructor of 
       *  the parent class Catalogue is called
			 *  @param obj catalogue to be moved
       */
      CatalogueChainMesh (CatalogueChainMesh&& obj) : Catalogue (obj) {
        m_cellsize = std::move(obj.m_cellsize);
        m_lim = std::move(obj.m_lim);
        m_part_catalogue = std::move(obj.m_part_catalogue);
        m_part_catalogue2 = std::move(obj.m_part_catalogue2);
        m_dimension = std::move(obj.m_dimension);
      };

      /**
       *  @brief Copy assignment for the CatalogueChainMesh class. The copy assignment of 
       *  the parent class Catalogue is called
			 *  @param obj catalogue to be copied
       *  @return object of type CatalogueChainMesh
       */
      CatalogueChainMesh& operator=(const CatalogueChainMesh& obj) 
      {
        Catalogue::operator=(obj);
        m_cellsize = obj.m_cellsize;
        m_lim = obj.m_lim;
        m_part_catalogue = obj.m_part_catalogue;
        m_part_catalogue2 = obj.m_part_catalogue2;
        m_dimension = obj.m_dimension;
        return *this;
      }

      /**
       *  @brief Move assignment for the CatalogueChainMesh class. The move assignment of 
       *  the parent class Catalogue is called
			 *  @param obj catalogue to be moved
       *  @return object of type CatalogueChainMesh
       */
      CatalogueChainMesh& operator=(CatalogueChainMesh&& obj) noexcept
      { 
        if (this == &obj) return *this;
        Catalogue::operator=(obj);
        m_cellsize = std::move(obj.m_cellsize);
        m_lim = std::move(obj.m_lim);
        m_part_catalogue = std::move(obj.m_part_catalogue);
        m_part_catalogue2 = std::move(obj.m_part_catalogue2);
        m_dimension = std::move(obj.m_dimension);
        return *this;
      } 

      /**
       *  @brief default destructor
       */
      ~CatalogueChainMesh () = default;

			/**
       *  @brief constructor that creates a catalogue of CatalogueChainMeshs in which the particles of 
       *  another catalogue are subdivided 
       *
       *  @param variables name of the variables that compose the space
       *
       *  @param cellsize the size of the cells
       *  
       *  @param cat shared_ptr of a catalogue of object to subdivide in the cells of the created catalogue
			 * 
			 *  @param cat2 shared_ptr of an optional second catalogue of object to provide further space limits
       *
       */
      CatalogueChainMesh (std::vector<catalogue::Var> variables, double cellsize, std::shared_ptr<cbl::catalogue::Catalogue> cat, std::shared_ptr<cbl::catalogue::Catalogue> cat2={});
    

			/**
       *  @name Member functions used to operate on the catalogue
       */
      ///@{


      /**
       *  @brief return the closer particle to a selected one (CatalogueChainMesh)
       *  @param pos position of the start
       *  @return a vector containing the index of the closer particle and the index of the cell in which it is located
       */    
      std::vector<unsigned int> Closer_object (const std::vector<double> pos); 

      /**
       *  @brief return the close particles to a selected one (CatalogueChainMesh)
       *  @param pos position of the start
       *  @param rMax maximum radius
       *  @param rMin minimum radius 
       *  @return a vector containing the index of the close particles 
       */ 
      std::vector<unsigned int> Close_objects (const std::vector<double> pos, const double rMax, const double rMin=0.); 

      /**
       *  @brief return the number of close particles to a selected one (CatalogueChainMesh)
       *  @param pos position of the start
       *  @param rMax maximum radius
       *  @param rMin minimum radius 
       *  @return the number of the particles in the distance range
       */ 
      unsigned int Count_objects (const std::vector<double> pos, const double rMax, const double rMin=0.); 

      /**
       *  @brief return an array with the N particles nearest the selected object, sorted by distance. (CatalogueChainMesh)
       *  @param pos position of the starting object
       *  @param N number of particles researched 
       */
      std::vector<unsigned int> N_nearest_objects (const std::vector<double> pos, const unsigned int N);

      /**
       *  @brief return an array with the N particles nearest the selected object, sorted by distance. (CatalogueChainMesh)
       *  @param obj index of the starting object
       *  @param N number of particles researched 
       */
      std::vector<unsigned int> N_nearest_objects (const unsigned int obj, const unsigned int N);

      /**
       *  @brief return a matrix with the N particles nearest each object of the catalogue, sorted by distance. (CatalogueChainMesh)
       *  @param N number of particles researched 
       */
      std::vector<std::vector<unsigned int>> N_nearest_objects_cat (const unsigned int N); 


      /**
       *  @brief find the cell of a particle stored in \e m_part and delete the particle i (CatalogueChainMesh)
       *  @param index index of the particle
       */
      void deletePart (const unsigned int index);

      /**
       *  @brief delete the particle part stored in \e m_part, variable of a given object of index i (CatalogueChainMesh)
       *  @param i index of the cell
       *  @param part index of the particle
       */
      void deletePart (const unsigned int i, const unsigned int p);

      ///@}

      /**
       *  @brief get the private member \e m_cellsize
       *  @return the cellsize of ChainMesh
       */
      double cellsize() const { return m_cellsize; }

      /**
       * @brief get the catalogue's limits 
       * @return private variable m_lim
       */
      std::vector<std::vector<double>> lim() const {return m_lim; }

      /**
       * @brief number of cells for side
       * @return private variable m_dimension
       */
      unsigned int dimension() const {return m_dimension; }

		};
	}
}

#endif
