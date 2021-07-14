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
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "ChainMesh.h"


using namespace std;

using namespace cbl;


// ============================================================================


void cbl::chainmesh::ChainMesh::set_par (const double cell_size, const long nDim)
{
  m_nDim = nDim;
  m_cell_size = cell_size;

  if (m_cell_size<=0)
    ErrorCBL("forbidden value for cell_size = "+conv(cell_size, par::fDP2), "set_par", "ChainMesh.cpp");

  m_Lim.resize(m_nDim, vector<double> (2,0));
  m_Delta.resize(m_nDim);
  m_nCell.resize(m_nDim);
  m_cell_to_index.resize(m_nDim);
  m_nCell_tot = 1;
}


// ============================================================================


cbl::chainmesh::ChainMesh::ChainMesh (const double cell_size, const long nDim) : m_nDim(nDim), m_cell_size(cell_size)
{
  set_par(cell_size, nDim);
}


// ============================================================================


long cbl::chainmesh::ChainMesh::pos_to_index (const vector<double> center) const
{
  vector<long> indx(m_nDim);

  for (long j=0; j<m_nDim; j++)
    indx[j] = min(long((center[j]-m_Lim[j][0])/m_cell_size), m_nCell[j]-1);

  long indx_tot = indx[m_nDim-1];

  for (long j=m_nDim-2; j>-1; j--) {
    long mult = 1;
    for (long k=j+1; k<m_nDim; k++) mult *= m_nCell[k];
    indx_tot += mult*indx[j];
  }

  return indx_tot;
}


// ============================================================================


long cbl::chainmesh::ChainMesh::inds_to_index (const vector<long> indx) const
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


void cbl::chainmesh::ChainMesh::index_to_inds (const long index, const vector<long> nn, vector<long> &indx) const
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


void cbl::chainmesh::ChainMesh::create_chain_mesh (const vector<vector<double> > data, const double rMAX, const double rMIN, const long nMAX, const long nMIN)
{
  // setting stuff, generalized for n(=nDim) dimensions

  m_rMAX = rMAX;  
  m_rMIN = rMIN;

  long nObj = data[0].size();
  m_List.erase(m_List.begin(), m_List.end()); m_List.resize(nObj, -1);
  m_NonEmpty_Cells.erase(m_NonEmpty_Cells.begin(), m_NonEmpty_Cells.end());

  m_nCell_tot = pow(nMAX,3)+1;
  m_nCell_NonEmpty = 0;
  m_multCell.erase(m_multCell.begin(), m_multCell.end()); m_multCell.resize(m_nDim, 1);

  double fact = 1.;

  check_memory(3.0, true, "cbl::chainmesh::ChainMesh::create_chain_mesh of ChainMesh.cpp");

  while (m_nCell_tot>pow(nMAX, 3) || m_nCell_tot<nMIN
	 || !check_memory(3.0, false, "cbl::chainmesh::ChainMesh::create_chain_mesh of ChainMesh.cpp")) {

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
    for (int i=1; i<m_nDim; i++)
      m_multCell[m_nDim-1-i] = m_nCell[m_nDim-i]*m_multCell[i-1];
    
    fact *= 1.1;

    m_Label.erase(m_Label.begin(), m_Label.end());

  }

  m_Label.resize(m_nCell_tot, -1);

  for (long i=0; i<nObj; i++) {
    vector<double> center(m_nDim);
    for (int j=0; j<m_nDim; j++)
      center[j] = data[j][i];
    long indx_tot = pos_to_index(center);
    m_NonEmpty_Cells.push_back(indx_tot);
    m_List[i] = m_Label[indx_tot];
    m_Label[indx_tot] = i;
  }

  get_searching_region(rMAX, rMIN);
  m_NonEmpty_Cells = different_elements (m_NonEmpty_Cells);
  m_nCell_NonEmpty = (long)m_NonEmpty_Cells.size();

}


// ============================================================================


void cbl::chainmesh::ChainMesh::create_chain_mesh_m2 (const vector<vector<double> > data)
{
  // setting stuff, generalized for n(=nDim) dimensions
  long nObj = data[0].size();

  for (int i=0; i<m_nDim; i++) {
    m_Lim[i][0] = Min(data[i]); m_Lim[i][1] = Max(data[i]);
    m_Delta[i] = m_Lim[i][1]-m_Lim[i][0];
    m_nCell[i] = nint(m_Delta[i]/m_cell_size);
    m_nCell_tot *= m_nCell[i];
  }

  m_List_index.erase(m_List_index.begin(),m_List_index.end()); m_List_index.resize(m_nCell_tot);
  m_NonEmpty_Cells.erase(m_NonEmpty_Cells.begin(), m_NonEmpty_Cells.end());
  m_nCell_NonEmpty = 0;

  for (long i=0; i<nObj; i++) {
    long indx_tot = pos_to_index(data[i]);
    m_NonEmpty_Cells.push_back(indx_tot);
    m_List_index[indx_tot].push_back(i);
  }

  m_NonEmpty_Cells = different_elements (m_NonEmpty_Cells);
  m_nCell_NonEmpty = (long)m_NonEmpty_Cells.size();
}


// ============================================================================


void cbl::chainmesh::ChainMesh::get_searching_region (const double r_max, const double r_min)
{

  int n_max = nint(r_max/m_cell_size);
  n_max = (n_max*m_cell_size<r_max) ? n_max+1 : n_max;

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

  if (r_min>0 && n_min>0) {
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


vector<long> cbl::chainmesh::ChainMesh::close_objects_cell (const int cell_index, const long ii) const
{
  // r2 != -1 ---> search in a nDim annulus from r1 to r2
  // r2 == -1 ---> search in a nDim sphere from center to r1

  vector<long> list;

  for (size_t i=0; i<m_search_region.size(); i++) {

    long k = min(max(m_search_region[i]+cell_index, (long)0), m_nCell_tot-1);

    if (k!=ii) {
      long j = m_Label[k];
      while (j>-1) {
	list.push_back(j);
	j = m_List[j];
      }
    }
  }

  return list;
}


// ============================================================================


vector<long> cbl::chainmesh::ChainMesh::close_objects (const vector<double> center, const long ii) const
{
  vector<long> list;

  if (long(center.size())!=m_nDim) ErrorCBL("point must have same dimensions of the chain-mesh", "close_objects", "ChainMesh.cpp");

  long center_indx = pos_to_index(center);

  for (unsigned long i=0; i<m_search_region.size(); i++) {

    long k = min(max(m_search_region[i]+center_indx, (long)0), m_nCell_tot-1);
    long j = m_Label[k];
    
    // NEW PART  -> consider only the cells inside the searching radius (+ those intersected)
    if (m_cell_size>=m_rMAX/m_nDim) {
      vector<int> indx(m_nDim, 0);

      for (int ii=0; ii<m_nDim; ii++) {
	indx[ii] = floor(k/m_multCell[ii]);
	k -= m_multCell[ii]*indx[ii];
      }
      
      double dist = 0.;
      for (int ii=0; ii<m_nDim; ii++)
	dist += ((indx[ii]+0.5)*m_cell_size-center[ii]+m_Lim[ii][0])*((indx[ii]+0.5)*m_cell_size-center[ii]+m_Lim[ii][0]);
      
      dist = sqrt(dist);
      
      double half_diag = m_cell_size*0.5*sqrt(m_nDim);
      double max_radius = m_rMAX+half_diag;
      double min_radius = (m_rMIN>0.) ? m_rMIN-half_diag : 0.;
      
      if (dist<max_radius && dist>min_radius)
	while (j>-1 && j>ii) {
	  list.push_back(j);
	  j = m_List[j];
	}
    }
    else
      while (j>-1 && j>ii) {
	list.push_back(j);
	j = m_List[j];
      }
  }

  return list;
}

// ============================================================================


void cbl::chainmesh::ChainMesh::normalize (std::vector<std::vector<double>> points, std::vector<double> values, const double rMAX)
{
  if (points[0].size()!=values.size()) ErrorCBL("the size of points is different from the size of values", "norm_grid", "ChainMesh.cpp");

  const int nparams = points.size();

  m_values = values;
  m_rMAX = rMAX;

  std::vector<std::vector<double>> extremals(2, std::vector<double> (nparams, 0));
  std::vector<double> delta (nparams, 0.);

  for(int N=0; N<nparams; N++) {
    extremals[0][N] = cbl::Max(points[N]);
    extremals[1][N] = cbl::Min(points[N]);
  }

  for(int N=0; N<nparams; N++) {
    delta[N] = extremals[0][N]-extremals[1][N];
    if (delta[N]==0) ErrorCBL("the grid cannot be normalised", "norm_grid", "ChainMesh.cpp");
    for(size_t ii=0; ii<points[0].size(); ii++){
      points[N][ii] = 100. - (200./(delta[N]))*(extremals[0][N] - points[N][ii]);
    }
  }

  m_points = points;
  m_delta = delta;
  m_extremals = extremals;

  long nMAX = floor(200/m_cell_size)+10*rMAX;

  create_chain_mesh(m_points, m_rMAX, -1., nMAX, 0);
}


// ============================================================================


double cbl::chainmesh::ChainMesh::interpolate (std::vector<double> xi, const int distNum)
{
  if (m_points.size() != xi.size()) ErrorCBL("the dimension of points is different from the dimension of xi", "interp_coord", "ChainMesh.cpp");
  if (m_delta.size()==0.) ErrorCBL("the chain mesh grid is not normalised!", "interp_coord", "ChainMesh.cpp");

  const int nparams = m_points.size();

  // stretch the interpolation sample
  for (int N=0; N<nparams; N++)
    xi[N] = 100. - (200./(m_delta[N]))*(m_extremals[0][N] - xi[N]);
  
  double interp_value = cbl::Min(m_values);
  std::vector<long> close;
  std::vector<double> distances;
  std::vector<double> indices;
  double dist, sum = 0.;
  double sum_distances = 0.;
  int npoints = 0;

  close = close_objects(xi);

  if (close.size()>0) {
    for (auto&&k : close) {
      dist = 0.;
      for(int N=0; N<nparams; N++)
	dist += pow((xi[N]-m_points[N][k]), 2);
      dist = sqrt(dist);
      distances.push_back(dist);
      indices.push_back(k);
    }
    
    cbl::sort_2vectors(distances.begin(), indices.begin(), distances.size());
    
    for (auto&&k : indices) {
      if (distances[npoints]>0.) {
	sum += m_values[k]/distances[npoints];
	sum_distances += 1./distances[npoints];
	npoints++;
      }
      if (npoints>distNum) break;
    }
    
    interp_value = sum/sum_distances;
  }
  
  distances.clear();
  indices.clear();
  
  return interp_value;
}


// ============================================================================


vector<double> cbl::chainmesh::ChainMesh::interpolate (std::vector<std::vector<double>> points, std::vector<double> values, std::vector<std::vector<double>> xi, const int distNum, const double rMAX)
{
  if (points.size() != xi.size()) ErrorCBL("the dimension of points is different from the dimension of xi", "interpolate", "ChainMesh.cpp");
  if (points[0].size() != values.size()) ErrorCBL("the size of points is different from the size of values", "interpolate", "ChainMesh.cpp");

  const int nparams = points.size();
  const int points_values = points[0].size();
  const int xi_values = xi[0].size();
  std::vector<std::vector<double>> unique_vect(nparams, std::vector<double> (points_values+xi_values, 0));

  std::vector<std::vector<double>> extremals(2, std::vector<double> (nparams, 0));
  double delta = 0.;
  for (int N=0; N<nparams; N++) {
    extremals[0][N] = std::max(cbl::Max(points[N]), cbl::Max(xi[N]));
    extremals[1][N] = std::min(cbl::Min(points[N]), cbl::Min(xi[N]));
  }

  for (int N=0; N<nparams; N++){
    delta = extremals[0][N]-extremals[1][N];
    
    for (int ii=0; ii<xi_values; ii++) {
      xi[N][ii] = 100. - (200./(delta))*(extremals[0][N] - xi[N][ii]);
      unique_vect[N][ii] = xi[N][ii];
    }
    
    for (int ii=0; ii<points_values; ii++) {
      points[N][ii] = 100. - (200./(delta))*(extremals[0][N] - points[N][ii]);
      unique_vect[N][xi_values+ii] = points[N][ii];
    }
    
  }

  long nMAX = floor(200/m_cell_size)+10*rMAX;

  std::vector<double> interp_values(xi_values, cbl::Min(values));

  create_chain_mesh(unique_vect, rMAX, -1., nMAX, 0);

#pragma omp parallel num_threads(omp_get_max_threads())
  { 
    std::vector<double> center(nparams);
    std::vector<long> close;
    std::vector<double> distances;
    std::vector<double> indices;
    double dist, sum_posterior = 0.;
    int npoints = 0;
  
#pragma omp for schedule(dynamic)
    for (int ii=0; ii<xi_values; ii++) {
      sum_posterior = 0.;
      npoints = 0;
      for (int N=0; N<nparams; N++) center[N] = xi[N][ii];
      close = close_objects(center, xi_values);
      if (close.size()>0) {
	for (auto&&k : close) {
	  dist = 0.;
	  for (int N=0; N<nparams; N++) dist += pow((center[N]-unique_vect[N][k]), 2);
	  dist = sqrt(dist);
	  distances.push_back(dist);
	  indices.push_back(k-xi_values);
	}

	cbl::sort_2vectors(distances.begin(), indices.begin(), distances.size());

	for (auto&&k : indices) {
	  sum_posterior += values[k];
	  npoints++;
	  if (npoints>distNum) break;
	}
	interp_values[ii] = sum_posterior/npoints;
      }

      distances.clear();
      indices.clear();
    }
  }
  
  return interp_values;
}


// ============================================================================


vector<long> cbl::chainmesh::ChainMesh::get_list (const long cell_index) const
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


void cbl::chainmesh::ChainMesh1D::set_par (const double cell_size, const vector<double> xx, const double rMAX, const double rMIN, const long nMAX, const long nMIN)
{
  ChainMesh::set_par(cell_size, 1);

  vector<vector<double>> data;
  data.push_back(xx);
  create_chain_mesh(data, rMAX, rMIN, nMAX, nMIN);
}


// ============================================================================


cbl::chainmesh::ChainMesh1D::ChainMesh1D (const double cell_size, const std::vector<double> xx, const double rMAX, const double rMIN, const long nMAX, const long nMIN) : ChainMesh(cell_size,1)
{
  set_par(cell_size, xx, rMAX, rMIN, nMAX, nMIN);
}


// ============================================================================


void cbl::chainmesh::ChainMesh2D::set_par (const double cell_size, const std::vector<double> xx, const std::vector<double> yy, const double rMAX, const double rMIN, const long nMAX, const long nMIN)
{
  ChainMesh::set_par(cell_size, 2);

  vector<vector<double>> data;
  data.push_back(xx);
  data.push_back(yy);
  create_chain_mesh(data, rMAX, rMIN, nMAX, nMIN);
}


// ============================================================================


cbl::chainmesh::ChainMesh2D::ChainMesh2D (const double cell_size, const std::vector<double> xx, const std::vector<double> yy, const double rMAX, const double rMIN, const long nMAX, const long nMIN) : ChainMesh(cell_size,2)
{
  set_par(cell_size, xx, yy, rMAX, rMIN, nMAX, nMIN);
}


// ============================================================================


void cbl::chainmesh::ChainMesh3D::set_par (const double cell_size, const std::vector<double> xx, const std::vector<double> yy, const std::vector<double> zz, const double rMAX, const double rMIN, const long nMAX, const long nMIN)
{
  ChainMesh::set_par(cell_size, 3);
  vector<vector<double>> data;
  data.push_back(xx);
  data.push_back(yy);
  data.push_back(zz);
  create_chain_mesh(data, rMAX, rMIN, nMAX, nMIN);
}


// ============================================================================


cbl::chainmesh::ChainMesh3D::ChainMesh3D (const double cell_size, const std::vector<double> xx, const std::vector<double> yy, const std::vector<double> zz, const double rMAX, const double rMIN, const long nMAX, const long nMIN) : ChainMesh(cell_size, 3)
{
  set_par(cell_size, xx, yy, zz, rMAX, rMIN, nMAX, nMIN);
}
