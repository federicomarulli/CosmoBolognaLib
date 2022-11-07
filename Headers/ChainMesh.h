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
 *******************************************************************/

/**
 *  @file Headers/ChainMesh.h
 *
 *  @brief Implementation of the chain-mesh data structure
 *
 *  This file defines the interface of the class ChainMesh, used for
 *  the chain-mesh method
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "Kernel.h"

#ifndef __CHAINMESH__
#define __CHAINMESH__


namespace cbl {

  /**
   *  @brief The namespace of the functions and classes used for the
   *  <B> chain-mesh method </B>
   *
   *  The \e chainmesh namespace contains all the main functions and
   *  classes of the chain-mesh method used for counting pairs and
   *  triplets
   */
  namespace chainmesh {

    /**
     *  @class ChainMesh ChainMesh.h "Headers/ChainMesh.h"
     *
     *  @brief The class ChainMesh
     *
     *  This class is used to handle objects of type <EM> ChainMesh
     *  </EM>
     */
    class ChainMesh {

    protected:

      /// the number of dimension
      int m_nDim;

      /// indexes in the i-th cell
      std::vector<long> m_cell_to_index;

      /// the size of the cell in unit of the interested quantity
      double m_cell_size;

      /// list of internal use
      std::vector<long> m_List;

      /// array containing the last particle of the chain-mesh in each cell
      std::vector<long> m_Label;

      /// index list of internal use
      std::vector<std::vector<long> > m_List_index;

      /// Min,Max limits of variable(s) used for the chain-mesh
      std::vector<std::vector<double> > m_Lim;

      /// Max-Min of variable(s) used for the chain-mesh
      std::vector<double> m_Delta;

      /// number of cell(s) for variable(s)
      std::vector<long> m_nCell;

      /// the list of cell around a generic center
      std::vector<long> m_search_region;

      /// the total number of cells
      long m_nCell_tot;

      /// the number of non-empty cells
      long m_nCell_NonEmpty;

      /// the total number of non-empty cells
      std::vector<long> m_NonEmpty_Cells;

      /// multiplicative factor to compute the projected number of cells 
      std::vector<long> m_multCell;

      /// the vector containing the sample points for the n-dim interpolation
      std::vector<std::vector<double>> m_points;

      /// the vector containing the values of the n-dim function on the sample points
      std::vector<double> m_values;

      /// the vector containing the extremals of the sample points coordinates
      std::vector<std::vector<double>> m_extremals;

      /// vetors of differences between the extremals of the sample points coordinates
      std::vector<double> m_delta;

      /// the maximum radius for the search of close points in the chain mesh
      double m_rMAX;

      /// the minimum radius for the search of close points in the chain mesh
      double m_rMIN;

    public:

      /**
       *  @brief default constructor
       *
       */
      ChainMesh () = default;

      /**
       *  @brief constructor
       *  @param cell_size double storing the cell size
       *  @param nDim the number of dimensions
       *
       */
      ChainMesh (const double cell_size, const long nDim);

      /**
       *  @brief default destructor
       */
      ~ChainMesh () = default;

      /**
       *  @brief function that set parameters for the chain-mesh
       *  @param cell_size double storing the cell size
       *  @param nDim the number of dimensions
       */
      void set_par (const double cell_size, const long nDim);

      /**
       *  @brief get the private member ChainMesh::m_nCell_tot
       *  @return total number of cells
       */
      long nCell() const { return m_nCell_tot; }

      /**
       *  @brief get the private member ChainMesh::m_nCell_NonEmpty
       *  @return number of non-empty cells
       */
      long nCell_NonEmpty() const { return m_nCell_NonEmpty; }

      /**
       *  @brief get the private member ChainMesh::m_NonEmpty_Cells
       *  @return indexes of non-empty cells
       */
      std::vector<long> NonEmpty_Cells() const { return m_NonEmpty_Cells; }

      /**
       * @brief get the index of the cell given the object coordinates
       * @param center the object coordinates
       * @return the cell index
       */
      long pos_to_index (const std::vector<double> center) const;

      /**
       * @brief get the unique index of the cell given the n indices
       * @param indx vector of the indices of the nD space cell
       * @return the cell unique index
       */
      long inds_to_index (const std::vector<long> indx) const;

      /**
       * @brief get the n indices given the unique index
       * @param index the unique index
       * @param nn the number of cells along the box axis
       * @param indx vector of the indices of the nD space cell
       */
      void index_to_inds (const long index, const std::vector<long> nn, std::vector<long> &indx) const;

      /**
       * @brief create the chain mesh
       * @param data the vector containing the coordinate of the object
       * @param rMax the maximum radius, to set the interal variable m_search_region
       * @param rMin the minimum radius, to set the interal variable m_search_region
       * @param nMAX maximum number of cells
       * @param nMIN minimum number of cells
       */
      void create_chain_mesh (const std::vector<std::vector<double> > data, const double rMax, const double rMin=-1., const long nMAX=300, const long nMIN=10);

      /**
       * @brief create the chain mesh
       * @param data the vector containing the coordinate of the object
       */
      void create_chain_mesh_m2 (const std::vector<std::vector<double> > data);

      /**
       * @brief set the internal variable m_search_region, the list
       * of cell around a generic center
       * @param r_max the maximum radius
       * @param r_min the minimum radius
       */
      void get_searching_region (const double r_max, const double r_min=-1.);

      /**
       * @brief get the indeces of the objects close to a cell
       * @param cell_index the cell index
       * @param ii the minimum index given in output
       * @return vector containing the index of the objects inside the cell
       */
      std::vector<long> close_objects_cell (const int cell_index, long ii=-1) const;

      /**
       * @brief get the indeces of the objects close to an object
       * @param center coordinates of an object 
       * @param ii the minimum index given in output
       * @return vector containing the index of the objects inside the cell
       */
      std::vector<long> close_objects (std::vector<double> center, long ii=-1) const;

      /**
       * @brief function to set a normalized (square/cubic) grid from
       * a sample of points, used for the N-dim interpolation
       *
       * @param points sample points 
       *
       * @param values the data on the sample points
       *
       * @param rMAX the maximum radius, to set the internal variable
       * m_search_region
       */
      void normalize (std::vector<std::vector<double>> points, std::vector<double> values, const double rMAX);

      /**
       * @brief N-dim interpolation of a set of N coordinates on a
       * normalised grid (see normalize)
       *
       * @param xi the input coordinates that has to be
       * interpolated
       *
       * @param distNum the number of closest points on which the
       * average will be done
       *
       * @return the interpolated value at the input coordinates
       */
      double interpolate (std::vector<double> xi, const int distNum);

      /**
       * @brief N-dim interpolation
       *
       * @param points sample points
       *
       * @param values sample values
       *
       * @param xi the coordinates to be interpolated
       *
       * @param distNum the number of closest points on which the
       * average will be done
       *
       * @param rMAX the maximum radius, to set the interal variable
       * m_search_region
       *
       * @return values interolated at the input coordinates
       */
      std::vector<double> interpolate (std::vector<std::vector<double>> points, std::vector<double> values, std::vector<std::vector<double>> xi, const int distNum, const double rMAX);

      /**
       * @brief get the index of the object inside a cell
       *
       * @param cell_index the cell index
       *
       * @return vector containing the index of the objects inside the
       * cell
       */
      std::vector<long> get_list (const long cell_index) const;

    };

    
    /**
     *  @class ChainMesh1D ChainMesh.h "Headers/ChainMesh.h"
     *
     *  @brief The class ChainMesh1D
     *
     *  This class is used to handle objects of type <EM> ChainMesh1D
     *  </EM>
     */
    class ChainMesh1D : public ChainMesh
    {
    public:
      /**
       *  @brief default constructor
       *
       */
      ChainMesh1D () = default;

      /**
       *  @brief default destructor
       */
      ~ChainMesh1D () = default;

      /**
       *  @brief function that set parameters for the chain-mesh
       *  @param cell_size double storing the cell size
       *  @param xx the vector with the variable used for the chain-mesh
       *  @param rMAX the maximum separation
       *  @param rMIN the minimum separation
       *  @param nMAX the allowed maximum number of chain-mesh cells in each dimension
       *  @param nMIN the allowed minimum number of chain-mesh cells in each dimension
       */
      void set_par (const double cell_size, const std::vector<double> xx, const double rMAX, const double rMIN=-1., const long nMAX=300, const long nMIN=0);

      /**
       *  @brief constructor
       *  @param cell_size double storing the cell size
       *  @param xx the vector with the variable used for the chain-mesh
       *  @param rMAX the maximum separation
       *  @param rMIN the minimum separation
       *  @param nMAX the allowed maximum number of chain-mesh cells in each dimension
       *  @param nMIN the allowed minimum number of chain-mesh cells in each dimension
       */
      ChainMesh1D (const double cell_size, const std::vector<double> xx, const double rMAX, const double rMIN=-1., const long nMAX=300, const long nMIN=0);
    };

    /**
     *  @class ChainMesh2D ChainMesh.h "Headers/ChainMesh.h"
     *
     *  @brief The class ChainMesh2D
     *
     *  This class is used to handle objects of type <EM> ChainMesh2D
     *  </EM>
     */
    class ChainMesh2D : public ChainMesh
    {
    public:
      /**
       *  @brief default constructor
       *  1D
       */
      ChainMesh2D () = default;

      /**
       *  @brief default destructor
       */
      ~ChainMesh2D () = default;

      /**
       *  @brief function that set parameters for the chain-mesh
       *  @param cell_size double storing the cell size
       *  @param xx the vector with the first variable used for the chain-mesh
       *  @param yy the vector with the second variable used for the chain-mesh
       *  @param rMAX the maximum separation
       *  @param rMIN the minimum separation
       *  @param nMAX the allowed maximum number of chain-mesh cells in each dimension
       *  @param nMIN the allowed minimum number of chain-mesh cells in each dimension
       */
      void set_par (const double cell_size, const std::vector<double> xx, const std::vector<double> yy, const double rMAX, const double rMIN=-1., const long nMAX=300, const long nMIN=0);

      /**
       *  @brief constructor
       *  @param cell_size double storing the cell size
       *  @param xx the vector with the first variable used for the chain-mesh
       *  @param yy the vector with the second variable used for the chain-mesh
       *  @param rMAX the maximum separation
       *  @param rMIN the minimum separation
       *  @param nMAX the allowed maximum number of chain-mesh cells in each dimension
       *  @param nMIN the allowed minimum number of chain-mesh cells in each dimension
       */
      ChainMesh2D (const double cell_size, const std::vector<double> xx, const std::vector<double> yy, const double rMAX, const double rMIN=-1., const long nMAX=300, const long nMIN=0);
    };

    /**
     *  @class ChainMesh3D ChainMesh.h "Headers/ChainMesh.h"
     *
     *  @brief The class ChainMesh3D
     *
     *  This class is used to handle objects of type <EM> ChainMesh3D
     *  </EM>
     */
    class ChainMesh3D : public ChainMesh
    {
    public:

      /**
       *  @brief default constructor
       *
       */
      ChainMesh3D () = default;

      /**
       *  @brief default destructor
       */
      ~ChainMesh3D () = default;

      /**
       *  @brief function that set parameters for the chain-mesh
       *  @param cell_size double storing the cell size
       *  @param xx the vector with the first variable used for the chain-mesh
       *  @param yy the vector with the second variable used for the chain-mesh
       *  @param zz the vector with the third variable used for the chain-mesh
       *  @param rMAX the maximum separation
       *  @param rMIN the minimum separation
       *  @param nMAX the allowed maximum number of chain-mesh cells
       *  in each dimension
       *  @param nMIN the allowed minimum number of chain-mesh cells
       *  in each dimension
       */
      void set_par (const double cell_size, const std::vector<double> xx, const std::vector<double> yy, const std::vector<double> zz, const double rMAX, const double rMIN=-1., const long nMAX=300, const long nMIN=0);

      /**
       *  @brief constructor
       *  @param cell_size double storing the cell size
       *  @param xx the vector with the first variable used for the chain-mesh
       *  @param yy the vector with the second variable used for the chain-mesh
       *  @param zz the vector with the third variable used for the chain-mesh
       *  @param rMAX the maximum separation
       *  @param rMIN the minimum separation
       *  @param nMAX the allowed maximum number of chain-mesh cells
       *  in each dimension
       *  @param nMIN the allowed minimum number of chain-mesh cells
       *  in each dimension
       */
      ChainMesh3D (const double cell_size, const std::vector<double> xx, const std::vector<double> yy, const std::vector<double> zz, const double rMAX, const double rMIN=-1., const long nMAX=300, const long nMIN=0);
    };

  }
}

#endif
