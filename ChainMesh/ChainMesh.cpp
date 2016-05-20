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
 *  @file ChainMesh/ChainMesh.cpp
 *
 *  @brief Methods of the class ChainMesh
 *
 *  This file contains the implementation of the chain-mesh method for
 *  n-dimensional data
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "ChainMesh.h"

using namespace cosmobl;


// ============================================================================


void cosmobl::ChainMesh::set_par (const double cell_size, const long nDim)  
{
  m_nDim = nDim;
  m_cell_size = cell_size;
  
  if (m_cell_size <= 0) {
    string Err = "Error in cosmobl::ChainMesh::set_par of ChainMesh.cpp: forbidden value for cell_size = "+conv(cell_size, par::fDP2);
    ErrorMsg(Err);
  }

  m_Lim.resize(m_nDim, vector<double> (2,0));
  m_Delta.resize(m_nDim);
  m_nCell.resize(m_nDim);
  m_cell_to_index.resize(m_nDim);
  m_nCell_tot = 1;
}


// ============================================================================


cosmobl::ChainMesh::ChainMesh (const double cell_size, const long nDim) : m_nDim(nDim), m_cell_size(cell_size) 
{
  set_par(cell_size, nDim);
}


// ============================================================================


long cosmobl::ChainMesh::pos_to_index (const vector<double> centre) const
{
  vector<long> indx(m_nDim);

  for (long j=0; j<m_nDim; j++)
    indx[j] = min(long((centre[j]-m_Lim[j][0])/m_cell_size),m_nCell[j]-1);
  
  long indx_tot = indx[m_nDim-1];

  for (long j=m_nDim-2; j>-1; j--) {
    long mult = 1;
    for (long k=j+1; k<m_nDim; k++) mult *= m_nCell[k];
    indx_tot += mult*indx[j];
  }

  return indx_tot;
}


// ============================================================================


long cosmobl::ChainMesh::inds_to_index (const vector<long> indx) const
{
  long indx_tot = indx[m_nDim-1];

  for (int j=m_nDim-2; j>-1; j--) {
    long mult = 1;
    for (int k=j+1; k<m_nDim; k++) mult *= m_nCell[k];
    indx_tot += mult*indx[j];
  }

  return indx_tot;
}


// ============================================================================


void cosmobl::ChainMesh::index_to_inds (const long index, const vector<long> nn, vector<long> &indx) const
{
  indx.resize(m_nDim, 0);
  long mult = 1;
  long sum = 0;
  for (int i=1; i<m_nDim; i++) mult *= nn[i];
  indx[0] = floor((index-sum)/mult); 
  sum = indx[0]*mult;

  for (int i=1; i<m_nDim; i++) {
    mult = 1;
    for (int k=i+1; k<m_nDim; k++) mult *= nn[k];
    indx[i] = floor((index-sum)/mult); 
    sum += indx[i]*mult;
  }
}


// ============================================================================


void cosmobl::ChainMesh::create_chain_mesh (const vector<vector<double> > data, const double rMAX, const long nMIN, const long nMAX) 
{ 
  // setting stuff, generalized for n(=nDim) dimensions 
  
  long nObj = data[0].size();
  m_List.erase(m_List.begin(), m_List.end()); m_List.resize(nObj, -1);
    
  m_nCell_tot = pow(nMAX,3)+1;
  
  double fact = 1.;

  check_memory(0.9, 1, "cosmobl::ChainMesh::create_chain_mesh of ChainMesh.cpp");
  
  while (m_nCell_tot>pow(nMAX,3) || m_nCell_tot<nMIN
	 || !check_memory(0.9, 0, "cosmobl::ChainMesh::create_chain_mesh of ChainMesh.cpp")) {
 
    m_nCell_tot = 1;
    m_cell_size *= fact;
    for (int i=0; i<m_nDim; i++) {
      m_Lim[i][0] = Min(data[i])-1.05*rMAX; 
      m_Lim[i][1] = Max(data[i])+1.05*rMAX;
      m_Delta[i] = m_Lim[i][1]-m_Lim[i][0];
      m_nCell[i] = nint(m_Delta[i]/m_cell_size);
      m_nCell_tot *= m_nCell[i];
      m_cell_to_index[i] = m_Delta[i]/m_nCell[i];
    }
    fact *= 1.1;

    m_Label.erase(m_Label.begin(), m_Label.end()); 

    if (m_nCell_tot<10) {
      string ERR = "Error in cosmobl::ChainMesh::create_chain_mesh! m_nCell_tot = "+conv(m_nCell_tot, par::fINT)+"!";
      ErrorMsg(ERR);
    }
    
  }
 
  if (m_nCell_tot>pow(nMAX,3) || m_nCell_tot<nMIN) {
    string ERR = "Error in cosmobl::ChainMesh::create_chain_mesh! m_nCell_tot = "+conv(m_nCell_tot, par::fINT)+", possible memory problems!";
    ErrorMsg(ERR);
  }

  m_Label.resize(m_nCell_tot, -1);
  
  for (long i=0; i<nObj; i++) {
    vector<double> centre(m_nDim);
    for (int j=0; j<m_nDim; j++)
      centre[j] = data[j][i];
    long indx_tot = pos_to_index(centre);
    m_List[i] = m_Label[indx_tot];
    m_Label[indx_tot] = i;
  }

  get_searching_region(rMAX);

}


// ============================================================================


void cosmobl::ChainMesh::create_chain_mesh_m2 (const vector<vector<double> > data) 
{
  // Setting stuff, generalized for n(=nDim) dimensions
  long nObj = data[0].size();

  for (int i=0; i<m_nDim; i++) {
    m_Lim[i][0] = Min(data[i]); m_Lim[i][1] = Max(data[i]);
    m_Delta[i] = m_Lim[i][1]-m_Lim[i][0];
    m_nCell[i] = nint(m_Delta[i]/m_cell_size);
    m_nCell_tot *= m_nCell[i];
  }

  m_List_index.erase(m_List_index.begin(),m_List_index.end()); m_List_index.resize(m_nCell_tot);

  for (long i=0; i<nObj; i++) {
    long indx_tot = pos_to_index(data[i]);
    m_List_index[indx_tot].push_back(i);
  }

}


// ============================================================================


void cosmobl::ChainMesh::get_searching_region (const double r_max, const double r_min) 
{
   int n_max = nint(r_max/m_cell_size);
   int n_min = nint(r_min/m_cell_size)-2; // exclusive
   vector<long> n_incl(m_nDim);

   for (int i=0; i<m_nDim; i++)
     n_incl[i] = 2*n_max+1;
   
   long sz_region = pow(2*n_max+1, m_nDim);

   m_search_region.erase(m_search_region.begin(),m_search_region.end());
   m_search_region.resize(sz_region,0);
   
   for (long i=0; i<sz_region; i++) {
     vector<long> indx(m_nDim);
     index_to_inds(i,n_incl,indx);
     for (int i=0; i<m_nDim; i++) indx[i]-=n_max;
     m_search_region[i] = inds_to_index(indx);
   }

   if (r_min >0 && n_min>0) {
     vector<long> n_excl(m_nDim);
     for (int i=0; i<m_nDim; i++)
       n_excl[i] = 2*n_min+1;

     long sz_excl_region = pow(2*n_min+1,m_nDim);

     for (long i=0; i<sz_excl_region; i++) {
       vector<long> indx(m_nDim);
       index_to_inds(i, n_excl, indx);
       for (int i=0; i<m_nDim; i++) indx[i] -= n_min;
       long veto = inds_to_index(indx);
       m_search_region.erase(remove(m_search_region.begin(), m_search_region.end(), veto), m_search_region.end());
     }
   }
}


// ============================================================================


vector<long> cosmobl::ChainMesh::close_objects (const vector<double> centre, const long ii) const
{
  // r2 != -1 ---> search in a nDim annulus from r1 to r2
  // r2 == -1 ---> search in a nDim sphere from centre to r1

  vector<long> list;
  
  if (long(centre.size()) != m_nDim) ErrorMsg("Error in ChainMesh::get_list : point must have same dimensions of the chain-mesh");

   long centre_indx = pos_to_index(centre);
   
   for (unsigned long i=0; i<m_search_region.size(); i++) {

     long k = min(max(m_search_region[i]+centre_indx, (long)0), m_nCell_tot-1);
     long j = m_Label[k];

     while (j>-1 && j>=ii) {
       list.push_back(j);
       j = m_List[j];
     }
   }

   return list;
}


// ============================================================================


vector<long> cosmobl::ChainMesh::get_list (const long cell_index) const
{
  vector<long> list;
  long j = m_Label[cell_index];

  while (j>-1) {
    list.push_back(j);
    j = m_List[j];
  }

  return list;
}


// ============================================================================


void cosmobl::ChainMesh1D::set_par (const double cell_size, const vector<double> xx, const double rMAX, const long nMIN, const long nMAX) 
{
  ChainMesh::set_par(cell_size, 1);
  
  vector<vector<double>> data; 
  data.push_back(xx);
  create_chain_mesh(data, rMAX, nMIN, nMAX);   
}


// ============================================================================


cosmobl::ChainMesh1D::ChainMesh1D (const double cell_size, const vector<double> xx, const double rMAX, const long nMIN, const long nMAX) : ChainMesh(cell_size,1)
{
  set_par(cell_size, xx, rMAX, nMIN, nMAX);
}


// ============================================================================


void cosmobl::ChainMesh2D::set_par (const double cell_size, const vector<double> xx, const vector<double> yy, const double rMAX, const long nMIN, const long nMAX) 
{
  ChainMesh::set_par(cell_size, 2);
  
  vector<vector<double>> data; 
  data.push_back(xx);
  data.push_back(yy);
  create_chain_mesh(data, rMAX, nMIN, nMAX);   
}


// ============================================================================


cosmobl::ChainMesh2D::ChainMesh2D (const double cell_size, const vector<double> xx, const vector<double> yy, const double rMAX, const long nMIN, const long nMAX) : ChainMesh(cell_size,2)
{
  set_par(cell_size, xx, yy, rMAX, nMIN, nMAX);
}


// ============================================================================


void cosmobl::ChainMesh3D::set_par (const double cell_size, const vector<double> xx, const vector<double> yy, const vector<double> zz, const double rMAX, const long nMIN, const long nMAX) 
{
  ChainMesh::set_par(cell_size, 3);
  vector<vector<double>> data;
  data.push_back(xx);
  data.push_back(yy);
  data.push_back(zz);
  create_chain_mesh(data, rMAX, nMIN, nMAX);
}


// ============================================================================


cosmobl::ChainMesh3D::ChainMesh3D (const double cell_size, const vector<double> xx, const vector<double> yy, const vector<double> zz, const double rMAX, const long nMIN, const long nMAX) : ChainMesh(cell_size, 3)
{
  set_par(cell_size, xx, yy, zz, rMAX, nMIN, nMAX);
}


