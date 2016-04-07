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
 *  @file Catalogue/ChainMesh_Catalogue.cpp
 *
 *  @brief Methods of the class ChainMesh_Catalogue 
 *
 *  This file contains the implementation of the methods of the class
 *  ChainMesh_Catalogue, to fast search between objects in Catalogues
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "ChainMesh_Catalogue.h"

using namespace cosmobl;
using namespace catalogue;


// ============================================================================


void cosmobl::catalogue::ChainMesh_Catalogue::set_par (const double cell_size, shared_ptr<Catalogue> cat, const double rmax) 
{
  ChainMesh::set_par(cell_size, 3);
  
  m_catalogue = cat;

  vector<vector<double>> data = {m_catalogue->var(Var::_X_), m_catalogue->var(Var::_Y_), m_catalogue->var(Var::_Z_)}; 

  create_chain_mesh(data, rmax);

  vector<int> order;
  get_order(order);
  m_catalogue->Order(order);

  data.erase(data.begin(), data.end());
  data = {m_catalogue->var(Var::_X_), m_catalogue->var(Var::_Y_), m_catalogue->var(Var::_Z_)}; 

  create_chain_mesh(data, rmax);
}


// ============================================================================


cosmobl::catalogue::ChainMesh_Catalogue::ChainMesh_Catalogue(const double cell_size, shared_ptr<Catalogue> cat, const double rmax)
  : ChainMesh(cell_size, 3)
{
  set_par(cell_size, cat, rmax); 
}


// ============================================================================


void cosmobl::catalogue::ChainMesh_Catalogue::get_order (vector<int> &order) const
{
  order.erase(order.begin(),order.end());
  vector<long> Label_temp = m_Label;
  vector<long> List_temp = m_List;

  for (int i=0; i<m_nCell_tot; i++) {
    vector<int> vv;
    int j = Label_temp[i];
    while (j>-1) {
      vv.push_back(j);
      j = List_temp[j];
    }

    reverse(vv.begin(), vv.end());
    
    for (int j=0; j<int(vv.size()); j++)
      order.push_back(vv[j]);
  }
}


// ============================================================================


vector<shared_ptr<Object> > cosmobl::catalogue::ChainMesh_Catalogue::object_list (shared_ptr<Object> object, const int ii)
{
  vector<shared_ptr<Object> > obj_list;

  vector<double> center = object->coords();
  int center_indx = pos_to_index(center);

  for (unsigned int i=0; i<m_search_region.size(); i++) {
    int k = min(max(m_search_region[i]+center_indx, (long)0), m_nCell_tot-1);
    int j = m_Label[k];

    while (j>-1 && j>=ii) {
      obj_list.push_back(m_catalogue->catalogue_object(j));
      j = m_List[j];
    }
    
  }
  
  return obj_list;
}

