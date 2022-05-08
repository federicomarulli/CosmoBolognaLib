/********************************************************************
 *  Copyright (C) 2022 by Federico Marulli, Simone Sartori          *
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
 *  @file Catalogue/CatalogueChainMesh.cpp
 *
 *  @brief Methods of the class Catalogue to construct chainmesh catalogues
 *
 *  This file contains the implementation of the methods of the class
 *  Catalogue used to create chainmesh catalogues
 *
 *  @author Federico Marulli, Simone Sartori
 *
 *  @author federico.marulli3@unibo.it, simone.sartori5@studio.unibo.it
 */

#include "CatalogueChainMesh.h"

using namespace std;
using namespace cbl;


// ============================================================================

cbl::catalogue::CatalogueChainMesh::CatalogueChainMesh (const std::function<double(double, double, double, double, double, double)> dist_function, double cellsize, std::shared_ptr<cbl::catalogue::Catalogue> cat)
{ 	
	m_part_catalogue = cat;
	vector<vector<double>> data = { m_part_catalogue->var(Var::_X_), m_part_catalogue->var(Var::_Y_), m_part_catalogue->var(Var::_Z_) }; 
	unsigned int nCells=0;
	vector<double> delta(3);
	double start_cellsize = cellsize;
	m_lim.resize(3);

	for (unsigned int i=0; i<3; i++) {
		m_lim[i].resize(2);
		m_lim[i][0] = *min_element(data[i].begin(), data[i].end()) - 0.05*start_cellsize;
		m_lim[i][1] = *max_element(data[i].begin(), data[i].end()) + 0.05*start_cellsize;
		delta[i] = m_lim[i][1]-m_lim[i][0];
	}
	
	double delta_max = *max_element(delta.begin(),delta.end());
	cellsize = delta_max/2;
	vector<vector<vector<vector<unsigned int>>>> cells(2, vector<vector<vector<unsigned int>>>(2, vector<vector<unsigned int>>(2))); 

	while(check_memory(2.0, false, "cbl::catalogue::Catalogue::Catalogue of ChainMeshCatalogue.cpp")==1 && cellsize > start_cellsize) {
		cellsize *= 0.9;
		nCells = delta_max/cellsize + 1;
		m_dimension = nCells;
		cells.resize(nCells);
		for (unsigned int i=0; i<nCells; i++) {
			cells[i].resize(nCells);
			for (unsigned int j=0; j<nCells; j++) { 
				cells[i][j].resize(nCells);
			}
		}
	}

	m_cellsize = cellsize;

	for (unsigned int i=0; i<data[0].size(); i++) {
		vector<unsigned int> data_inds(3);
		for (unsigned int j=0; j<3; j++) data_inds[j] = ((unsigned int)((data[j][i] - m_lim[j][0])/cellsize));
		cells[data_inds[0]][data_inds[1]][data_inds[2]].emplace_back(i);
	}

	unsigned int dim_nearCells = nCells*sqrt(3.)+1;

	// include the objects in the catalogue
	comovingCoordinates coord;

	for (unsigned int i=0; i<nCells; i++) {
		coord.xx = (double)i;
		for (unsigned int j=0; j<nCells; j++) {
				coord.yy = (double)j;
			for (unsigned int k=0; k<nCells; k++) {
					coord.zz = (double)k;
				unsigned int index = k+j*nCells+i*nCells*nCells;
				vector<vector<unsigned int>> nearCells(dim_nearCells);

				for (unsigned int ii=0; ii<nCells; ii++) {
					for (unsigned int jj=0; jj<nCells; jj++) {
						for (unsigned int kk=0; kk<nCells; kk++) {
							unsigned int dist_indx = (dist_function((double)i, ((double)i-(double)ii==0) ? (double) ii : ( ((double)i-(double)ii>=0) ? (double)ii+0.5 : (double)ii-0.5), 
													(double)j, ((double)j-(double)jj==0) ? (double) jj : ( ((double)j-(double)jj>=0) ? (double)jj+0.5 : (double)jj-0.5), 
													(double)k, ((double)k-(double)kk==0) ? (double) kk : ( ((double)k-(double)kk>=0) ? (double)kk+0.5 : (double)kk-0.5)));   
							unsigned int index2 = kk+jj*nCells+ii*nCells*nCells;
							nearCells[dist_indx].emplace_back(index2);
						}
					}
				}
				check_memory(4.0, true, "cbl::catalogue::Catalogue::Catalogue of ChainMeshCatalogue.cpp. Increase cellsize");
				add_object(move(Object::Create(coord, index, cells[i][j][k], nearCells)));
			}
		}
	}
}

	//VISUALIZER
	/*
	 for (unsigned int i=0; i<nCells[0]; i++) {
		coord.xx = (double)i;
		for (unsigned int j=0; j<nCells[1]; j++) {
				coord.yy = (double)j;
			for (unsigned int k=0; k<nCells[2]; k++) {
					coord.zz = (double)k;
				vector<vector<unsigned int>> nearCells(dim_nearCells);

				for (unsigned int ii=0; ii<nCells[0]; ii++) {
					for (unsigned int jj=0; jj<nCells[1]; jj++) {
						for (unsigned int kk=0; kk<nCells[2]; kk++) {
							unsigned int dist_indx = (dist_function((double)i, ((double)i-(double)ii==0) ? (double) ii : ( ((double)i-(double)ii>=0) ? (double)ii+0.5 : (double)ii-0.5), 
													(double)j, ((double)j-(double)jj==0) ? (double) jj : ( ((double)j-(double)jj>=0) ? (double)jj+0.5 : (double)jj-0.5), 
													(double)k, ((double)k-(double)kk==0) ? (double) kk : ( ((double)k-(double)kk>=0) ? (double)kk+0.5 : (double)kk-0.5)));              
							unsigned int index2 = kk+jj*nCells[2]+ii*nCells[2]*nCells[1];
							cout << setw(3) << index2 << "(" << dist_indx << ")  ";
						}
						cout << endl;
					}
					cout << endl;
				}
				cout << endl;
				cout << "----------------" << endl;
				cout << endl;
			}
		}
	}
	*/


// ============================================================================

cbl::catalogue::CatalogueChainMesh::CatalogueChainMesh (const std::function<double(double, double, double, double, double, double)> dist_function,  double cellsize, std::shared_ptr<cbl::catalogue::Catalogue> cat, std::shared_ptr<cbl::catalogue::Catalogue> cat2)
{ 
	m_part_catalogue = cat;
	m_part_catalogue2 = cat2;
	vector<vector<double>> data = { cat->var(Var::_X_), cat->var(Var::_Y_), cat->var(Var::_Z_) }; 
	vector<vector<double>> data2 = {};
	if (cat2->nObjects()>0) data2 = { cat2->var(Var::_X_), cat2->var(Var::_Y_), cat2->var(Var::_Z_) }; 
	unsigned int nCells=0;
	vector<unsigned int> num_cells(3);
	vector<double> delta(3);
	double start_cellsize = cellsize;
	m_lim.resize(3);

	for (unsigned int i=0; i<3; i++) {
		m_lim[i].resize(2);
		m_lim[i][0] = *min_element(data[i].begin(), data[i].end());
		m_lim[i][0] = (cat2->nObjects() > 0) ? min(m_lim[i][0], *min_element(data2[i].begin(), data2[i].end()))-0.05*cellsize :  m_lim[i][0]-0.05*cellsize;
		m_lim[i][1] = *max_element(data[i].begin(), data[i].end());
		m_lim[i][1] = (cat2->nObjects() > 0) ? max(m_lim[i][1], *max_element(data2[i].begin(), data2[i].end()))+0.05*cellsize :  m_lim[i][1]+0.05*cellsize;   
		delta[i] = m_lim[i][1]-m_lim[i][0];
	}

	m_part_catalogue = cat;
	m_part_catalogue2 = cat2;

	double delta_max = *max_element(delta.begin(),delta.end());
	cellsize = delta_max/2;
	vector<vector<vector<vector<unsigned int>>>> cells(2, vector<vector<vector<unsigned int>>>(2, vector<vector<unsigned int>>(2))); 

	while(check_memory(2.0, false, "cbl::catalogue::Catalogue::Catalogue of ChainMeshCatalogue.cpp")==1 && cellsize > start_cellsize) {
		cellsize *= 0.9;
		nCells = delta_max/cellsize + 1;
		m_dimension = nCells;
		cells.resize(nCells);
		for (unsigned int i=0; i<nCells; i++) {
			cells[i].resize(nCells);
			for (unsigned int j=0; j<nCells; j++) { 
				cells[i][j].resize(nCells);
			}
		}
	}

	m_cellsize = cellsize;

	for (unsigned int i=0; i<data[0].size(); i++) {
		vector<unsigned int> data_inds(3);
		for (unsigned int j=0; j<3; j++) data_inds[j] = ((unsigned int)((data[j][i] - m_lim[j][0])/cellsize));
		cells[data_inds[0]][data_inds[1]][data_inds[2]].push_back(i);
	}


	unsigned int dim_nearCells = nCells*sqrt(3.)+1;

	// include the objects in the catalogue

	comovingCoordinates coord;

	for (unsigned int i=0; i<nCells; i++) {
		coord.xx = (double)i;
		for (unsigned int j=0; j<nCells; j++) {
				coord.yy = (double)j;
			for (unsigned int k=0; k<nCells; k++) {
					coord.zz = (double)k;
				unsigned int index = k+j*nCells+i*nCells*nCells;
				vector<vector<unsigned int>> nearCells(dim_nearCells);

				for (unsigned int ii=0; ii<nCells; ii++) {
					for (unsigned int jj=0; jj<nCells; jj++) {
						for (unsigned int kk=0; kk<nCells; kk++) {
							unsigned int dist_indx = (dist_function((double)i, ((double)i-(double)ii==0) ? (double) ii : ( ((double)i-(double)ii>=0) ? (double)ii+0.5 : (double)ii-0.5), 
													(double)j, ((double)j-(double)jj==0) ? (double) jj : ( ((double)j-(double)jj>=0) ? (double)jj+0.5 : (double)jj-0.5), 
													(double)k, ((double)k-(double)kk==0) ? (double) kk : ( ((double)k-(double)kk>=0) ? (double)kk+0.5 : (double)kk-0.5)));   
							unsigned int index2 = kk+jj*nCells+ii*nCells*nCells;
							nearCells[dist_indx].push_back(index2);
						}
					}
				}
				check_memory(4.0, true, "cbl::catalogue::Catalogue::Catalogue of ChainMeshCatalogue.cpp. Increase cellsize");
				add_object(move(Object::Create(coord, index, cells[i][j][k], nearCells)));
			}
		}
	}
}

// =========================================================================

std::vector<unsigned int> cbl::catalogue::CatalogueChainMesh::Closer_object (const std::vector<double> pos) 
{
  vector<int> inds(3);
  for (int i=0; i<3; i++) inds[i] = (int)((pos[i] - m_lim[i][0])/m_cellsize);

	int center_index = inds[2] + inds[1]*m_dimension + inds[0]*m_dimension*m_dimension;

  unsigned int n_start_max = 1;
	unsigned int n_start_min = 0;
	unsigned int n_start = n_start_max;
  vector<unsigned int> cObject(2);
	double distance = (unsigned int)m_cellsize*10000;
  vector<vector<unsigned int>> nCells = nearCells(center_index);
	double radius = m_cellsize;

	int add = 3;

	while(distance == (unsigned int)m_cellsize*10000)
	{
		for(unsigned int i = n_start_min; i<n_start_max; i++)
		{
			if (i == nCells.size()) break;
			for(unsigned int j : nCells[i])
			{
				for (unsigned int k : part(j))
				{
					double new_distance = Euclidean_distance(pos[0], m_part_catalogue->xx(k), pos[1],  m_part_catalogue->yy(k), pos[2],  m_part_catalogue->zz(k));
					if(new_distance < distance && new_distance < radius && new_distance > n_start_min*m_cellsize) 
					{	
						distance = new_distance;
						cObject[0] = k;
						cObject[1] = j; 
					}
				}
			} 				
		}
				
		n_start_max += 1;
		n_start_min = (n_start_max <= n_start+add) ? 0 : n_start_min+1;
		radius +=  m_cellsize;
	}

  return cObject;
}

// =========================================================================

std::vector<unsigned int> cbl::catalogue::CatalogueChainMesh::Close_objects (const std::vector<double> pos, const double rMax, const double rMin) 
{
  vector<int> inds(3);
  for (int i=0; i<3; i++) inds[i] = (int)((pos[i] - m_lim[i][0])/m_cellsize);
	int center_index = inds[2] + inds[1]*m_dimension + inds[0]*m_dimension*m_dimension;

  unsigned int n_max = (rMax/m_cellsize)+1;
	n_max = (m_cellsize*n_max < rMax) ?  n_max+2 : n_max+1; 
  unsigned int n_min = max(nint(rMin/m_cellsize)-3, 0); 

  vector<vector<unsigned int>> nCells = nearCells(center_index);
	vector<unsigned int> objects;

	for(unsigned int i = n_min; i<n_max; i++)
	{
		if (i == nCells.size()) break;
		for(unsigned int j : nCells[i])
		{
			for (unsigned int k : part(j))
			{
				double new_distance = Euclidean_distance(pos[0], m_part_catalogue->xx(k), pos[1],  m_part_catalogue->yy(k), pos[2],  m_part_catalogue->zz(k));
				if(new_distance < rMax && new_distance >= rMin) 
				{	
					objects.emplace_back(k);
				}
			}
		} 				
	}

  return objects;
}

// =========================================================================

unsigned int cbl::catalogue::CatalogueChainMesh::Count_objects (const std::vector<double> pos, const double rMax, const double rMin) 
{
  vector<int> inds(3);
  for (int i=0; i<3; i++) inds[i] = (int)((pos[i] - m_lim[i][0])/m_cellsize);

	int center_index = inds[2] + inds[1]*m_dimension + inds[0]*m_dimension*m_dimension;

	int add = 3;

  unsigned int n_max = nint(rMax/m_cellsize);
	n_max = (m_cellsize*n_max < rMax) ?  n_max+2 : n_max+1; 

  unsigned int n_min = max(nint(rMin/m_cellsize)-add, 0); 

  vector<vector<unsigned int>> nCells = nearCells(center_index);
	unsigned int objects = 0;

	for(unsigned int i = n_min; i<n_max; i++)
	{
		if (i == nCells.size()) break;
		for(unsigned int j : nCells[i])
		{
			for (unsigned int k : part(j))
			{
				double new_distance = Euclidean_distance(pos[0], m_part_catalogue->xx(k), pos[1],  m_part_catalogue->yy(k), pos[2],  m_part_catalogue->zz(k));
				if(new_distance < rMax && new_distance >= rMin) 
				{	
					objects+=1;
				}
			}
		} 				
	}

  return objects;
}

// =========================================================================

std::vector<unsigned int> cbl::catalogue::CatalogueChainMesh::N_nearest_objects (const std::vector<double> pos, const unsigned int N) 
{
  vector<int> inds(3);
  for (int i=0; i<3; i++) inds[i] = (int)((pos[i] - m_lim[i][0])/m_cellsize);

	int center_index = inds[2] + inds[1]*m_dimension + inds[0]*m_dimension*m_dimension;

	vector<unsigned int> Object={}, temp_Object={};						
	vector<double> distance_obj={}, temp_distance_obj={};;
	vector<vector<unsigned int>> nCells = nearCells(center_index);
	double radius=m_cellsize, new_distance;
	size_t index=0;
	vector<bool> mask;
	while (Object.size() < N) {
		
		for (size_t i=0; i<temp_Object.size(); i++) if (temp_distance_obj[i]<radius*(index+1) && mask[i]==false) mask[i] = true; 

		for(unsigned int j : nCells[index]) {
			for (unsigned int k : part(j)) {
				new_distance = Euclidean_distance(pos[0], m_part_catalogue->xx(k), pos[1],  m_part_catalogue->yy(k), pos[2],  m_part_catalogue->zz(k));
				temp_Object.push_back(k);
				temp_distance_obj.push_back(new_distance);
				if(new_distance < radius*(index+1)) mask.emplace_back(true);
				else mask.emplace_back(false);
			}
		}	 			

		if (index==nCells.size()-1) {
			Object = temp_Object;
			distance_obj = temp_distance_obj;
			break;
		} 
		else {
			unsigned int Ntrue=count(mask.begin(), mask.end(), true);
			if (Ntrue>=N) {
				for (size_t i=0; i<temp_Object.size(); i++) {
					if (mask[i]==true) {
						Object.emplace_back(temp_Object[i]);
						distance_obj.emplace_back(temp_distance_obj[i]);
					}
				}
			} 	
		}

		index++;			
	}
	
	std::vector<double> doublecObj(Object.begin(), Object.end());
	sort_2vectors(distance_obj.begin(), doublecObj.begin(), Object.size());
	doublecObj.erase(doublecObj.begin()+N, doublecObj.end());
	std::vector<unsigned int> cObject(doublecObj.begin(), doublecObj.end());

  return cObject;
}

// =========================================================================

std::vector<unsigned int> cbl::catalogue::CatalogueChainMesh::N_nearest_objects (const unsigned int obj, const unsigned int N) 
{
  vector<int> inds(3);
	vector<double> pos = {m_part_catalogue->xx(obj), m_part_catalogue->yy(obj), m_part_catalogue->zz(obj)};
	for (int i=0; i<3; i++) inds[i] = (int)((pos[i] - m_lim[i][0])/m_cellsize);
	int center_index = inds[2] + inds[1]*m_dimension + inds[0]*m_dimension*m_dimension;

	vector<unsigned int> Object={}, temp_Object={};						
	vector<double> distance_obj={}, temp_distance_obj={};;
	vector<vector<unsigned int>> nCells = nearCells(center_index);
	double radius=m_cellsize, new_distance;
	size_t index=0;
	vector<bool> mask;
	while (Object.size() < N) {
		
		for (size_t i=0; i<temp_Object.size(); i++) if (temp_distance_obj[i]<radius*(index+1) && mask[i]==false) mask[i] = true; 

		for(unsigned int j : nCells[index]) {
			for (unsigned int k : part(j)) {
				new_distance = Euclidean_distance(pos[0], m_part_catalogue->xx(k), pos[1],  m_part_catalogue->yy(k), pos[2],  m_part_catalogue->zz(k));
				temp_Object.push_back(k);
				temp_distance_obj.push_back(new_distance);
				if(new_distance < radius*(index+1)) mask.emplace_back(true);
				else mask.emplace_back(false);
			}
		}	 			

		if (index==nCells.size()-1) {
			Object = temp_Object;
			distance_obj = temp_distance_obj;
			break;
		} 
		else {
			unsigned int Ntrue=count(mask.begin(), mask.end(), true);
			if (Ntrue>=N) {
				for (size_t i=0; i<temp_Object.size(); i++) {
					if (mask[i]==true) {
						Object.emplace_back(temp_Object[i]);
						distance_obj.emplace_back(temp_distance_obj[i]);
					}
				}
			} 	
		}

		index++;			
	}

	std::vector<double> doublecObj(Object.begin(), Object.end());
	sort_2vectors(distance_obj.begin(), doublecObj.begin(), Object.size());
	doublecObj.erase(doublecObj.begin()+N, doublecObj.end());
	std::vector<unsigned int> cObject(doublecObj.begin(), doublecObj.end());

  return cObject;
}

// =========================================================================

std::vector<std::vector<unsigned int>> cbl::catalogue::CatalogueChainMesh::N_nearest_objects_cat (const unsigned int N) 
{
	vector<vector<unsigned int>> near_part_cat(m_part_catalogue->nObjects());
#pragma omp parallel num_threads(omp_get_max_threads())
	{	    
#pragma omp for schedule(static)
		for (size_t i=0; i<m_part_catalogue->nObjects(); i++) {
			near_part_cat[i] = N_nearest_objects(i, N);
		}
	}

	return near_part_cat;
}

// =========================================================================

void cbl::catalogue::CatalogueChainMesh::deletePart (const unsigned int index) 
{ 
  vector<int> inds(3);
	vector<double> pos = {m_part_catalogue -> xx(index), m_part_catalogue -> yy(index), m_part_catalogue -> zz(index)};
	 
  for (int i=0; i<3; i++) inds[i] = (int)((pos[i] - m_lim[i][0])/m_cellsize);
	

	int center_index = inds[2] + inds[1]*m_dimension + inds[0]*m_dimension*m_dimension;
	vector<unsigned int> particles = part(center_index);
  particles.erase(remove(particles.begin(), particles.end(), index), particles.end());

  set_part(center_index, particles);
	
}  

// =========================================================================

void cbl::catalogue::CatalogueChainMesh::deletePart (const unsigned int i, const unsigned int p) 
{ 

	vector<unsigned int> particles = part(i);
  particles.erase(remove(particles.begin(), particles.end(), p), particles.end());
  set_part(i, particles);
	
}   
