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
 *  @file Headers/Lib/ChainMesh.h
 *
 *  @brief Implementation of the chain-mesh data structure
 *
 *  This file defines the interface of the class ChainMesh, used for
 *  the chain-mesh method
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Cosmology.h"


#ifndef __CHAINMESH__
#define __CHAINMESH__


namespace cosmobl {

  /**
   *  @class ChainMesh ChainMesh.h "Headers/Lib/ChainMesh.h"
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
    vector<long> m_cell_to_index;

    /// the size of the cell in unit of the interested quantity
    double m_cell_size;

    /// list of internal use
    vector<long> m_List;
    
    /// array containing the last particle of the chain-mesh in each cell
    vector<long> m_Label;

    /// index list of internal use
    vector<vector<long> > m_List_index;

    /// Min,Max limits of variable(s) used for the chain-mesh
    vector<vector<double> > m_Lim;

    /// Max-Min of variable(s) used for the chain-mesh
    vector<double> m_Delta;

    /// number of cell(s) for variable(s) 
    vector<long> m_nCell;

    /// the list of cell around a generic centre
    vector<long> m_search_region;

    /// the total number of cells
    long m_nCell_tot;

  public:
    
    /**
     *  @brief default constructor
     *  @return object of class ChainMesh
     */
    ChainMesh () {}

    /**
     *  @brief default destructor
     *  @return none
     */
    ~ChainMesh () {}

    /**
     *  @brief function that set parameters for the chain-mesh 
     *  @param cell_size double storing the cell size
     *  @param nDim the number of dimensions
     *  @return none
     */
    void set_par (const double, const long);

    /**
     *  @brief constructor 
     *  @param cell_size double storing the cell size
     *  @param nDim the number of dimensions
     *  @return object of class ChainMesh
     */
    ChainMesh (double const, long const);
    
    /**
     *  @brief get the private member ChainMesh::m_nCell_tot
     *  @return total number of cells 
     */
    long nCell() const { return m_nCell_tot; }
    
    long pos_to_index (const vector<double>) const;
    
    long inds_to_index (const vector<long>) const;
    
    void index_to_inds (const long, const vector<long>, vector<long> &) const;

    void create_chain_mesh (const vector<vector<double> >, const double , const long nMIN=0, const long nMAX=300);
    
    void create_chain_mesh_m2 (const vector<vector<double> >);

    void get_searching_region (const double, const double r_min = -1);
    
    vector<long> close_objects (vector<double>, long ii=-1) const; 

    vector<long> get_list (const long) const;

  };

  /**
   *  @class ChainMesh1D ChainMesh.h "Headers/Lib/ChainMesh.h"
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
     *  @return an object of class ChainMesh1D
     */
    ChainMesh1D () {}

    /**
     *  @brief default destructor
     *  @return none
     */
    ~ChainMesh1D () {}

    /**
     *  @brief function that set parameters for the chain-mesh 
     *  @param cell_size double storing the cell size
     *  @param xx the vector with the variable used for the chain-mesh
     *  @param rMAX the maximum separation
     *  @param nMIN the allowed minimum number of chain-mesh cells in each dimension 
     *  @param nMAX the allowed maximum number of chain-mesh cells in each dimension 
     *  @return none
     */
    void set_par (const double, const vector<double>, const double, const long nMIN=0, const long nMAX=300);

    /**
     *  @brief constructor 
     *  @param cell_size double storing the cell size
     *  @param xx the vector with the variable used for the chain-mesh
     *  @param rMAX the maximum separation
     *  @param nMIN the allowed minimum number of chain-mesh cells in each dimension 
     *  @param nMAX the allowed maximum number of chain-mesh cells in each dimension 
     *  @return object of the class ChainMesh1D
     */
    ChainMesh1D (const double, const vector<double>, const double, const long nMIN=0, const long nMAX=300);
  };

  /**
   *  @class ChainMesh2D ChainMesh.h "Headers/Lib/ChainMesh.h"
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
     *  @return object of class ChainMesh1D
     */
    ChainMesh2D () {}

    /**
     *  @brief default destructor
     *  @return none
     */
    ~ChainMesh2D () {}

    /**
     *  @brief function that set parameters for the chain-mesh 
     *  @param cell_size double storing the cell size
     *  @param xx the vector with the first variable used for the chain-mesh
     *  @param yy the vector with the second variable used for the chain-mesh
     *  @param rMAX the maximum separation
     *  @param nMIN the allowed minimum number of chain-mesh cells in each dimension 
     *  @param nMAX the allowed maximum number of chain-mesh cells in each dimension 
     *  @return none
     */
    void set_par (const double, const vector<double>, const vector<double>, const double, const long nMIN=0, const long nMAX=300);

    /**
     *  @brief constructor 
     *  @param cell_size double storing the cell size
     *  @param xx the vector with the first variable used for the chain-mesh
     *  @param yy the vector with the second variable used for the chain-mesh
     *  @param rMAX the maximum separation
     *  @param nMIN the allowed minimum number of chain-mesh cells in each dimension 
     *  @param nMAX the allowed maximum number of chain-mesh cells in each dimension 
     *  @return object of the class ChainMesh2D
     */
    ChainMesh2D (const double, const vector<double>, const vector<double>, const double, const long nMIN=0, const long nMAX=300);
  };

  /**
   *  @class ChainMesh3D ChainMesh.h "Headers/Lib/ChainMesh.h"
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
     *  @return an object of class ChainMesh3D
     */
    ChainMesh3D () {}

    /**
     *  @brief default destructor
     *  @return none
     */
    ~ChainMesh3D () {}

    /**
     *  @brief function that set parameters for the chain-mesh 
     *  @param cell_size double storing the cell size
     *  @param xx the vector with the first variable used for the chain-mesh
     *  @param yy the vector with the second variable used for the chain-mesh
     *  @param zz the vector with the third variable used for the chain-mesh
     *  @param rMAX the maximum separation
     *  @param nMIN the allowed minimum number of chain-mesh cells in each dimension 
     *  @param nMAX the allowed maximum number of chain-mesh cells in each dimension 
     *  @return none
     */
    void set_par (const double, const vector<double>, const vector<double>, const vector<double>, const double, const long nMIN=0, const long nMAX=300);

    /**
     *  @brief constructor 
     *  @param cell_size double storing the cell size
     *  @param xx the vector with the first variable used for the chain-mesh
     *  @param yy the vector with the second variable used for the chain-mesh
     *  @param zz the vector with the third variable used for the chain-mesh
     *  @param rMAX the maximum separation
     *  @param nMIN the allowed minimum number of chain-mesh cells in each dimension 
     *  @param nMAX the allowed maximum number of chain-mesh cells in each dimension 
     *  @return object of the class ChainMesh3D
     */
    ChainMesh3D (const double, const vector<double>, const vector<double>, const vector<double>, const double, const long nMIN=0, const long nMAX=300);
  };

}

#endif
