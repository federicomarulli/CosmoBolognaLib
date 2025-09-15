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
 *  @file MeshGrid/MeshGrid.cpp
 *
 *  @brief Methods of the class MeshGrid
 *
 *  This file contains the implementation of the mesh-grid method for
 *  3D data
 *
 *  @authors Simone Sartori
 *
 *  @authors sartori@cppm.in2p3.fr
 */

#include "MeshGrid.h"

using namespace std;
using namespace cbl;


// ===========================================================================


cbl::meshgrid::MeshGrid::MeshGrid (const std::vector<double> X, const std::vector<double> Y, const std::vector<double> Z, double cellsize, std::vector<std::vector<double>> lims)
{
  m_X = X;
  m_Y = Y;
  m_Z = Z;

  m_cellsize = cellsize;

  // limits setting 
  vector<vector<double>> temp_lims(3, vector<double>(2));
  if (lims.empty()) {
    temp_lims[0][0] = *min_element(X.begin(), X.end())-0.001;
    temp_lims[1][0] = *min_element(Y.begin(), Y.end())-0.001;
    temp_lims[2][0] = *min_element(Z.begin(), Z.end())-0.001;
    temp_lims[0][1] = *max_element(X.begin(), X.end())+0.001;
    temp_lims[1][1] = *max_element(Y.begin(), Y.end())+0.001;
    temp_lims[2][1] = *max_element(Z.begin(), Z.end())+0.001;
  }
  else {
    if (lims.size()!=3 || lims[0].size()!=2) ErrorCBL("The limits matrix should be 3x2", "MeshGrid", "MeshGrid.cpp");
    temp_lims = lims;
  }

  m_nCells.resize(3);
  m_nCells[0] = ceil((temp_lims[0][1]-temp_lims[0][0])/m_cellsize);
  m_nCells[1] = ceil((temp_lims[1][1]-temp_lims[1][0])/m_cellsize);
  m_nCells[2] = ceil((temp_lims[2][1]-temp_lims[2][0])/m_cellsize);

  m_lims.resize(3, std::vector<double>(2));
  m_lims[0][0] = temp_lims[0][0] - (m_nCells[0]*cellsize-(temp_lims[0][1]-temp_lims[0][0]))/2;
  m_lims[1][0] = temp_lims[1][0] - (m_nCells[1]*cellsize-(temp_lims[1][1]-temp_lims[1][0]))/2;
  m_lims[2][0] = temp_lims[2][0] - (m_nCells[2]*cellsize-(temp_lims[2][1]-temp_lims[2][0]))/2;
  m_lims[0][1] = temp_lims[0][1] + (m_nCells[0]*cellsize-(temp_lims[0][1]-temp_lims[0][0]))/2;
  m_lims[1][1] = temp_lims[1][1] + (m_nCells[1]*cellsize-(temp_lims[1][1]-temp_lims[1][0]))/2;
  m_lims[2][1] = temp_lims[2][1] + (m_nCells[2]*cellsize-(temp_lims[2][1]-temp_lims[2][0]))/2;

  // mask creation 
  unsigned int dim_mask = ceil(sqrt(3)+sqrt(m_nCells[0]*m_nCells[0]+m_nCells[1]*m_nCells[1]+m_nCells[2]*m_nCells[2]));

  m_mask.resize(dim_mask, std::vector<std::vector<int>>());

  for (int i=-(int)m_nCells[0]; i<(int)m_nCells[0]; i++) {
    for (int j=-(int)m_nCells[1]; j<(int)m_nCells[1]; j++) {
      for (int k=-(int)m_nCells[2]; k<(int)m_nCells[2]; k++) {
        int dx = (i<0 ? i+1 : (i>0 ? i-1 : 0));
        int dy = (j<0 ? j+1 : (j>0 ? j-1 : 0));
        int dz = (k<0 ? k+1 : (k>0 ? k-1 : 0));
        unsigned int indexes = floor(sqrt(dx*dx+dy*dy+dz*dz));
        m_mask[indexes].push_back({i,j,k});
      }
    }
  }

  // grid creation 
  m_grid.resize(m_nCells[0], vector<vector<vector<unsigned int>>>(m_nCells[1], vector<vector<unsigned int>>(m_nCells[2], vector<unsigned int>())));

  for (size_t i=0; i<m_X.size(); i++) {
    unsigned int Xind = floor((m_X[i]-m_lims[0][0])/m_cellsize);
    unsigned int Yind = floor((m_Y[i]-m_lims[1][0])/m_cellsize);
    unsigned int Zind = floor((m_Z[i]-m_lims[2][0])/m_cellsize);
    m_grid[Xind][Yind][Zind].push_back(i);
  }
}


// ===========================================================================


unsigned int cbl::meshgrid::MeshGrid::closerObject (const double X, const double Y, const double Z)
{
  if (X<m_lims[0][0] || X>m_lims[0][1]) ErrorCBL("X coordinate outside boundaries", "closerObject", "MeshGrid.cpp");
  if (Y<m_lims[1][0] || Y>m_lims[1][1]) ErrorCBL("Y coordinate outside boundaries", "closerObject", "MeshGrid.cpp");
  if (Z<m_lims[2][0] || Z>m_lims[2][1]) ErrorCBL("Z coordinate outside boundaries", "closerObject", "MeshGrid.cpp");

  unsigned int Xind = floor((X-m_lims[0][0])/m_cellsize);
  unsigned int Yind = floor((Y-m_lims[1][0])/m_cellsize);
  unsigned int Zind = floor((Z-m_lims[2][0])/m_cellsize);

  double dist = 1e30;
  unsigned int maskInd = 0;
  unsigned int finalInd = 4294967295;
  while (dist == 1e30) {
    for (size_t i=0; i<m_mask[maskInd].size(); i++) {
      int Xind_temp = m_mask[maskInd][i][0]+Xind;
      int Yind_temp = m_mask[maskInd][i][1]+Yind;
      int Zind_temp = m_mask[maskInd][i][2]+Zind;
      if (Xind_temp >= 0 && Xind_temp<(int)m_nCells[0] && Yind_temp >= 0 && Yind_temp<(int)m_nCells[1] && Zind_temp >= 0 && Zind_temp<(int)m_nCells[2]) {
        for (auto& part : m_grid[Xind_temp][Yind_temp][Zind_temp]) {
          double dist_temp = sqrt((m_X[part]-X)*(m_X[part]-X)+(m_Y[part]-Y)*(m_Y[part]-Y)+(m_Z[part]-Z)*(m_Z[part]-Z));
          if (dist_temp < dist) {
            dist = dist_temp;
            finalInd = part;
          }
        }
      }
      else continue;
    }
  }

  return finalInd;

}


// ===========================================================================


unsigned int cbl::meshgrid::MeshGrid::closerObject (const unsigned int index)
{
  if (index >= m_X.size()) ErrorCBL("index too high", "closerObject", "MeshGrid.cpp");
  
  double X = m_X[index];
  double Y = m_Y[index];
  double Z = m_Z[index];

  unsigned int Xind = floor((X-m_lims[0][0])/m_cellsize);
  unsigned int Yind = floor((Y-m_lims[1][0])/m_cellsize);
  unsigned int Zind = floor((Z-m_lims[2][0])/m_cellsize);

  double dist = 1e30;
  unsigned int maskInd = 0;
  unsigned int finalInd = 4294967295;
  while (dist == 1e30) {
    for (size_t i=0; i<m_mask[maskInd].size(); i++) {
      int Xind_temp = m_mask[maskInd][i][0]+Xind;
      int Yind_temp = m_mask[maskInd][i][1]+Yind;
      int Zind_temp = m_mask[maskInd][i][2]+Zind;
      if (Xind_temp >= 0 && Xind_temp<(int)m_nCells[0] && Yind_temp >= 0 && Yind_temp<(int)m_nCells[1] && Zind_temp >= 0 && Zind_temp<(int)m_nCells[2]) {
        for (auto& part : m_grid[Xind_temp][Yind_temp][Zind_temp]) {
          if (part != index) {
            double dist_temp = sqrt((m_X[part]-X)*(m_X[part]-X)+(m_Y[part]-Y)*(m_Y[part]-Y)+(m_Z[part]-Z)*(m_Z[part]-Z));
            if (dist_temp < dist) {
              dist = dist_temp;
              finalInd = part;
            }
          }
        }
      }
      else continue;
    }
  }

  return finalInd;
}


// ===========================================================================


vector<unsigned int> cbl::meshgrid::MeshGrid::N_closerObjects (const unsigned int N, const double X, const double Y, const double Z)
{
  if (X<m_lims[0][0] || X>m_lims[0][1]) ErrorCBL("X coordinate outside boundaries", "N_closerObjects", "MeshGrid.cpp");
  if (Y<m_lims[1][0] || Y>m_lims[1][1]) ErrorCBL("Y coordinate outside boundaries", "N_closerObjects", "MeshGrid.cpp");
  if (Z<m_lims[2][0] || Z>m_lims[2][1]) ErrorCBL("Z coordinate outside boundaries", "N_closerObjects", "MeshGrid.cpp");
  if (N>m_X.size()) ErrorCBL("N bigger then the number of objects", "N_closerObjects", "MeshGrid.cpp");

  unsigned int Xind = floor((X-m_lims[0][0])/m_cellsize);
  unsigned int Yind = floor((Y-m_lims[1][0])/m_cellsize);
  unsigned int Zind = floor((Z-m_lims[2][0])/m_cellsize);
  vector<unsigned int> finalInd(N);
  unsigned int partConsidered = 0, allPart = 0;
  unsigned int maskIndexMax = 1;
  double distMax = m_cellsize;
  while (partConsidered < N) {
    partConsidered = 0;
    allPart = 0;
    for (unsigned int maskInd=0; maskInd<maskIndexMax; maskInd++) {
      for (size_t i=0; i<m_mask[maskInd].size(); i++) {
        int Xind_temp = m_mask[maskInd][i][0]+Xind;
        int Yind_temp = m_mask[maskInd][i][1]+Yind;
        int Zind_temp = m_mask[maskInd][i][2]+Zind;
        if (Xind_temp >= 0 && Xind_temp<(int)m_nCells[0] && Yind_temp >= 0 && Yind_temp<(int)m_nCells[1] && Zind_temp >= 0 && Zind_temp<(int)m_nCells[2]) {
          for (auto& part : m_grid[Xind_temp][Yind_temp][Zind_temp]) {
            double dist_temp = sqrt((m_X[part]-X)*(m_X[part]-X)+(m_Y[part]-Y)*(m_Y[part]-Y)+(m_Z[part]-Z)*(m_Z[part]-Z));
            allPart++;
            if (dist_temp < distMax) 
              partConsidered++;
          }
        }
        else continue;
      }
    }
    distMax += m_cellsize;
    maskIndexMax++;
  }

  vector<pair<double, unsigned int>> objects(allPart);
  unsigned int cont=0;
  for (unsigned int maskInd=0; maskInd<maskIndexMax-1; maskInd++) {
    for (size_t i=0; i<m_mask[maskInd].size(); i++) {
      int Xind_temp = m_mask[maskInd][i][0]+Xind;
      int Yind_temp = m_mask[maskInd][i][1]+Yind;
      int Zind_temp = m_mask[maskInd][i][2]+Zind;
      if (Xind_temp >= 0 && Xind_temp<(int)m_nCells[0] && Yind_temp >= 0 && Yind_temp<(int)m_nCells[1] && Zind_temp >= 0 && Zind_temp<(int)m_nCells[2]) {
        for (auto& part : m_grid[Xind_temp][Yind_temp][Zind_temp]) {
          objects[cont].first = sqrt((m_X[part]-X)*(m_X[part]-X)+(m_Y[part]-Y)*(m_Y[part]-Y)+(m_Z[part]-Z)*(m_Z[part]-Z));
          objects[cont].second = part;
          cont++;
        }
      }
      else continue;
    }
  }  

  sort(objects.begin(), objects.end(), [](const auto &a, const auto &b) 
  {return a.first < b.first; });

  for (unsigned int i=0; i<N; i++) 
    finalInd[i] = objects[i].second;

  return finalInd;

}


// ===========================================================================


vector<unsigned int> cbl::meshgrid::MeshGrid::N_closerObjects (const unsigned int N, const unsigned int index)
{
  if (index >= m_X.size()) ErrorCBL("index too high", "N_closerObjects", "MeshGrid.cpp");
  if (N>m_X.size()) ErrorCBL("N bigger then the number of objects", "N_closerObjects", "MeshGrid.cpp");

  double X = m_X[index];
  double Y = m_Y[index];
  double Z = m_Z[index];

  unsigned int Xind = floor((X-m_lims[0][0])/m_cellsize);
  unsigned int Yind = floor((Y-m_lims[1][0])/m_cellsize);
  unsigned int Zind = floor((Z-m_lims[2][0])/m_cellsize);

  vector<unsigned int> finalInd(N);
  unsigned int partConsidered = 0, allPart = 0;
  unsigned int maskIndexMax = 1;
  double distMax = m_cellsize;
  while (partConsidered < N) {
    partConsidered = 0;
    allPart = 0;
    for (unsigned int maskInd=0; maskInd<maskIndexMax; maskInd++) {
      for (size_t i=0; i<m_mask[maskInd].size(); i++) {
        int Xind_temp = m_mask[maskInd][i][0]+Xind;
        int Yind_temp = m_mask[maskInd][i][1]+Yind;
        int Zind_temp = m_mask[maskInd][i][2]+Zind;
        if (Xind_temp >= 0 && Xind_temp<(int)m_nCells[0] && Yind_temp >= 0 && Yind_temp<(int)m_nCells[1] && Zind_temp >= 0 && Zind_temp<(int)m_nCells[2]) {
          for (auto& part : m_grid[Xind_temp][Yind_temp][Zind_temp]) {
            if (part != index) {
              double dist_temp = sqrt((m_X[part]-X)*(m_X[part]-X)+(m_Y[part]-Y)*(m_Y[part]-Y)+(m_Z[part]-Z)*(m_Z[part]-Z));
              allPart++;
              if (dist_temp < distMax) 
                partConsidered++;
            }
          }
        }
        else continue;
      }
    }
    distMax += m_cellsize;
    maskIndexMax++;
  }

  vector<pair<double, unsigned int>> objects(allPart);
  unsigned int cont=0;
  for (unsigned int maskInd=0; maskInd<maskIndexMax-1; maskInd++) {
    for (size_t i=0; i<m_mask[maskInd].size(); i++) {
      int Xind_temp = m_mask[maskInd][i][0]+Xind;
      int Yind_temp = m_mask[maskInd][i][1]+Yind;
      int Zind_temp = m_mask[maskInd][i][2]+Zind;
      if (Xind_temp >= 0 && Xind_temp<(int)m_nCells[0] && Yind_temp >= 0 && Yind_temp<(int)m_nCells[1] && Zind_temp >= 0 && Zind_temp<(int)m_nCells[2]) {
        for (auto& part : m_grid[Xind_temp][Yind_temp][Zind_temp]) {
          if (part != index) {
            objects[cont].first = sqrt((m_X[part]-X)*(m_X[part]-X)+(m_Y[part]-Y)*(m_Y[part]-Y)+(m_Z[part]-Z)*(m_Z[part]-Z));
            objects[cont].second = part;
            cont++;
          }
        }
      }
      else continue;
    }
  }

  sort(objects.begin(), objects.end(), [](const auto &a, const auto &b) 
  {return a.first < b.first; });

  for (unsigned int i=0; i<N; i++) 
    finalInd[i] = objects[i].second;

  return finalInd;

}


// ===========================================================================


vector<unsigned int> cbl::meshgrid::MeshGrid::closeObjects (const double X, const double Y, const double Z, const double Rmax, const double Rmin)
{
  if (X<m_lims[0][0] || X>m_lims[0][1]) ErrorCBL("X coordinate outside boundaries", "closeObjects", "MeshGrid.cpp");
  if (Y<m_lims[1][0] || Y>m_lims[1][1]) ErrorCBL("Y coordinate outside boundaries", "closeObjects", "MeshGrid.cpp");
  if (Z<m_lims[2][0] || Z>m_lims[2][1]) ErrorCBL("Z coordinate outside boundaries", "closeObjects", "MeshGrid.cpp");
  if (Rmax < 0.) ErrorCBL("Rmax < 0.", "closeObjects", "MeshGrid.cpp");
  if (Rmin < 0.) ErrorCBL("Rmin < 0.", "closeObjects", "MeshGrid.cpp");
  if (Rmax < Rmin) ErrorCBL("Rmax < Rmin", "closeObjects", "MeshGrid.cpp");

  unsigned int Xind = floor((X-m_lims[0][0])/m_cellsize);
  unsigned int Yind = floor((Y-m_lims[1][0])/m_cellsize);
  unsigned int Zind = floor((Z-m_lims[2][0])/m_cellsize);

  unsigned int maskIndexMax = floor(Rmax/m_cellsize);
  unsigned int maskIndexMin = std::max((int)floor(Rmin/m_cellsize)-3, (int)0);

  vector<unsigned int> finalInd;
  for (unsigned int maskInd=maskIndexMin; maskInd<maskIndexMax+1; maskInd++) {
    for (size_t i=0; i<m_mask[maskInd].size(); i++) {
      int Xind_temp = m_mask[maskInd][i][0]+Xind;
      int Yind_temp = m_mask[maskInd][i][1]+Yind;
      int Zind_temp = m_mask[maskInd][i][2]+Zind;
      if (Xind_temp >= 0 && Xind_temp<(int)m_nCells[0] && Yind_temp >= 0 && Yind_temp<(int)m_nCells[1] && Zind_temp >= 0 && Zind_temp<(int)m_nCells[2]) {
        for (auto& part : m_grid[Xind_temp][Yind_temp][Zind_temp]) {
          double dist_temp = sqrt((m_X[part]-X)*(m_X[part]-X)+(m_Y[part]-Y)*(m_Y[part]-Y)+(m_Z[part]-Z)*(m_Z[part]-Z));
          if (dist_temp <= Rmax && dist_temp >= Rmin) 
            finalInd.push_back(part);
        }
      }
      else continue;
    }
  }

  return finalInd;
}


// ===========================================================================


vector<unsigned int> cbl::meshgrid::MeshGrid::closeObjects (const unsigned int index, const double Rmax, const double Rmin)
{
  if (index >= m_X.size()) ErrorCBL("index too high", "closerObject", "MeshGrid.cpp");
  
  double X = m_X[index];
  double Y = m_Y[index];
  double Z = m_Z[index];

  unsigned int Xind = floor((X-m_lims[0][0])/m_cellsize);
  unsigned int Yind = floor((Y-m_lims[1][0])/m_cellsize);
  unsigned int Zind = floor((Z-m_lims[2][0])/m_cellsize);

  unsigned int maskIndexMax = floor(Rmax/m_cellsize);
  unsigned int maskIndexMin = std::max((int)floor(Rmin/m_cellsize)-3, (int)0);

  vector<unsigned int> finalInd;
  for (unsigned int maskInd=maskIndexMin; maskInd<maskIndexMax+1; maskInd++) {
    for (size_t i=0; i<m_mask[maskInd].size(); i++) {
      int Xind_temp = m_mask[maskInd][i][0]+Xind;
      int Yind_temp = m_mask[maskInd][i][1]+Yind;
      int Zind_temp = m_mask[maskInd][i][2]+Zind;
      if (Xind_temp >= 0 && Xind_temp<(int)m_nCells[0] && Yind_temp >= 0 && Yind_temp<(int)m_nCells[1] && Zind_temp >= 0 && Zind_temp<(int)m_nCells[2]) {
        for (auto& part : m_grid[Xind_temp][Yind_temp][Zind_temp]) {
          double dist_temp = sqrt((m_X[part]-X)*(m_X[part]-X)+(m_Y[part]-Y)*(m_Y[part]-Y)+(m_Z[part]-Z)*(m_Z[part]-Z));
          if (dist_temp <= Rmax && dist_temp > Rmin) 
            finalInd.push_back(part);
        }
      }
      else continue;
    }
  }

  return finalInd;
}
