/********************************************************************
 *  Copyright (C) 2024 by Federico Marulli & Simone Sartori         *
 *  federico.marulli3@unibo.it sartori@cppm.in2p3.fr                *
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
 *******************************************************************/

/**
 *  @file Headers/MeshGrid.h
 *
 *  @brief Implementation of the chain-mesh data structure
 *
 *  This file defines the interface of the class MeshGrid, used for
 *  the mesh-grid method
 *
 *  @authors Simone Sartori
 *
 *  @authors sartori@cppm.in2p3.fr
 */

#include "Kernel.h"

#ifndef __MESHGRID__
#define __MESHGRID__

namespace cbl {

  /**
   *  @brief The namespace of the functions and classes used for the
   *  <B> chain-mesh method </B>
   *
   *  The \e meshgrid namespace contains all the main functions and
   *  classes of the meshgrid method used for subdividing 
   *  the space in order to access in a faster way to the particles
   */
  namespace meshgrid {

    /**
     *  @class MeshGrid MeshGrid.h "Headers/MeshGrid.h"
     *
     *  @brief The class MeshGrid
     *
     *  This class is used to handle objects of type <EM> MeshGrid
     *  </EM>
     */
    class MeshGrid {

    protected:

      /// X coordinate vector
      std::vector<double> m_X;

      /// Y coordinate vector
      std::vector<double> m_Y;

      /// Z coordinate vector
      std::vector<double> m_Z;

      /// number of cells for each side (cubic box)
      std::vector<unsigned int> m_nCells;

      /// the size of the cell 
      double m_cellsize;

      /// limits of the box 
      std::vector<std::vector<double>> m_lims;

      /// cells mask for particle searching
      std::vector<std::vector<std::vector<int>>> m_mask;

      /// the mesh grid 
      std::vector<std::vector<std::vector<std::vector<unsigned int>>>> m_grid;

    public:

      /**
       *  @brief default constructor
       *
       */
      MeshGrid () = default;

      /**
       *  @brief default destructor
       */
      ~MeshGrid () = default;

      /**
       *  @brief constructor which initialize the mesh grid
       *  @param X X coordinates of the objects
       *  @param Y Y coordinates of the objects
       *  @param Z Z coordinates of the objects
       *  @param cellsize size of the cells
       *  @param lims cartesian limits of the catalogue (calculated if not provided)
       */
      MeshGrid (const std::vector<double> X, const std::vector<double> Y, const std::vector<double> Z, double cellsize, std::vector<std::vector<double>> lims={});

      /**
       *  @brief return the index of the closest object from a given point in the space (could refer to any point in the space inside the limits)
       *  @param X X coordinate 
       *  @param Y Y coordinate
       *  @param Z Z coordinate
       *  @return index of the closest object
       */
      unsigned int closerObject (const double X, const double Y, const double Z);

      /**
       *  @brief return the index of the closest object from a given particle (could refer to any point in the space inside the limits)
       *  @param index index of the particle
       *  @return index of the closest object
       */
      unsigned int closerObject (const unsigned int index);

      /**
       *  @brief return the indexes of the N closest objects from a given point in the space 
       *  (could refer to any point in the space inside the limits)
       *  @param N Number of close objects
       *  @param X X coordinate 
       *  @param Y Y coordinate
       *  @param Z Z coordinate
       *  @return vector of indexes of the N closest objects
       */
      std::vector<unsigned int> N_closerObjects (const unsigned int N, const double X, const double Y, const double Z);

      /**
       *  @brief return the indexes of the N closest objects from a given particle (could refer to any point in the space inside the limits)
       *  @param N Number of close objects
       *  @param index index of the particle
       *  @return vector of indexes of the N closest objects
       */
      std::vector<unsigned int> N_closerObjects (const unsigned int N, const unsigned int index);

      /**
       *  @brief return a vector of indexes of the close objects from a given point on the space (could refer to any point in the space inside the limits)
       *  @param X X coordinate 
       *  @param Y Y coordinate
       *  @param Z Z coordinate
       *  @param Rmax maximum radius from the point
       *  @param Rmin minimum radius from the point (0 by default)
       *  @return vector of indexes
       */
      std::vector<unsigned int> closeObjects (const double X, const double Y, const double Z, const double Rmax, const double Rmin=0.);

      /**
       *  @brief return a vector of indexes of the close objects from a given particle 
       *  @param index index of the particle
       *  @param Rmax maximum radius from the point
       *  @param Rmin minimum radius from the point (0 by default)
       *  @return vector of indexes
       */
      std::vector<unsigned int> closeObjects (const unsigned int index, const double Rmax, const double Rmin=0.);

      /**
       *  @brief return a vector of 3 indexes (of the grid) of the cell in which the coordinates fell 
       *  @param X X coordinate 
       *  @param Y Y coordinate
       *  @param Z Z coordinate
       *  @return vector of three indexes 
       */
      std::vector<unsigned int> get_Cell(const double X, const double Y, const double Z) 
      {
        if (X<m_lims[0][0] || X>m_lims[0][1]) ErrorCBL("X coordinate outside boundaries", "getCell", "MeshGrid.h");
        if (Y<m_lims[1][0] || Y>m_lims[1][1]) ErrorCBL("Y coordinate outside boundaries", "getCell", "MeshGrid.h");
        if (Z<m_lims[2][0] || Z>m_lims[2][1]) ErrorCBL("Z coordinate outside boundaries", "getCell", "MeshGrid.h");

        unsigned int Xind = floor((X-m_lims[0][0])/m_cellsize);
        unsigned int Yind = floor((Y-m_lims[1][0])/m_cellsize);
        unsigned int Zind = floor((Z-m_lims[2][0])/m_cellsize);

        return {Xind, Yind, Zind};
      }

      /**
       *  @brief return a vector of indexes of the objects of the cell in which the coordinates fell 
       *  @param X X coordinate 
       *  @param Y Y coordinate
       *  @param Z Z coordinate
       *  @return vector of indexes 
       */
      std::vector<unsigned int> get_ParticlesInCell(const double X, const double Y, const double Z) 
      {
        std::vector<unsigned int> cellCord = get_Cell(X, Y, Z);
        return m_grid[cellCord[0]][cellCord[1]][cellCord[2]];
      }

      /**
       *  @brief delete an object of a given index (the index MUST be the index related at the object in the moment 
       *  in which the function is called, NOT the original index if the function was already called)
       *  @param index index of the object
       */
      void deleteParticle(const unsigned int index)
      {
        if (index>=m_X.size()) ErrorCBL("index too high", "deleteParticle", "MeshGrid.h");
        std::vector<unsigned int> cellCord = get_Cell(m_X[index], m_Y[index], m_Z[index]);
        m_grid[cellCord[0]][cellCord[1]][cellCord[2]].erase(std::remove(m_grid[cellCord[0]][cellCord[1]][cellCord[2]].begin(), m_grid[cellCord[0]][cellCord[1]][cellCord[2]].end(), index), m_grid[cellCord[0]][cellCord[1]][cellCord[2]].end());
      }

      /**
       *  @brief add an object of given X,Y,Z coordinates
       *  @param X X coordinate 
       *  @param Y Y coordinate
       *  @param Z Z coordinate
       */
      void addParticle(const double X, const double Y, const double Z)
      {
        if (X<m_lims[0][0] || X>m_lims[0][1]) ErrorCBL("X coordinate outside boundaries", "addParticle", "MeshGrid.h");
        if (Y<m_lims[1][0] || Y>m_lims[1][1]) ErrorCBL("Y coordinate outside boundaries", "addParticle", "MeshGrid.h");
        if (Z<m_lims[2][0] || Z>m_lims[2][1]) ErrorCBL("Z coordinate outside boundaries", "addParticle", "MeshGrid.h");

        std::vector<unsigned int> cellCord = get_Cell(X, Y, Z);
        m_grid[cellCord[0]][cellCord[1]][cellCord[2]].push_back(m_Z.size());

        m_X.push_back(X);
        m_Y.push_back(Y);
        m_Z.push_back(Z);

      }      
    };
  }
}

#endif
