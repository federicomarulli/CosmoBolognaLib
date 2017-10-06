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
 *  @file GlobalFunc/SubSample.cpp
 *
 *  @brief Functions for dividing a catalogue in sub-samples
 *
 *  This file contains the implementation of a set of functions to
 *  divide a catalogue in sub-samples
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "GlobalFunc.h"

using namespace cosmobl;


// ============================================================================


void cosmobl::set_ObjectRegion_SubBoxes (catalogue::Catalogue &data, catalogue::Catalogue &random, const int nx, const int ny, const int nz)
{
  const double xMin = data.Min(catalogue::Var::_X_);
  const double yMin = data.Min(catalogue::Var::_Y_);
  const double zMin = data.Min(catalogue::Var::_Z_);
  
  const double Cell_X = (data.Max(catalogue::Var::_X_)-xMin)/nx;
  const double Cell_Y = (data.Max(catalogue::Var::_Y_)-yMin)/ny;
  const double Cell_Z = (data.Max(catalogue::Var::_Z_)-zMin)/nz;

#pragma omp parallel num_threads(omp_get_max_threads())
  {
    
#pragma omp for schedule(static, 2) 
    for (size_t i=0; i<data.nObjects(); i++) {
      const int i1 = min(int((data.xx(i)-xMin)/Cell_X), nx-1);
      const int j1 = min(int((data.yy(i)-yMin)/Cell_Y), ny-1);
      const int z1 = min(int((data.zz(i)-zMin)/Cell_Z), nz-1);
      const int index = z1+nz*(j1+ny*i1);
      data.catalogue_object(i)->set_region(index);
    }

#pragma omp for schedule(static, 2) 
    for (size_t i=0; i<random.nObjects(); i++) {
      const int i1 = min(int((random.xx(i)-xMin)/Cell_X), nx-1);
      const int j1 = min(int((random.yy(i)-yMin)/Cell_Y), ny-1);
      const int z1 = min(int((random.zz(i)-zMin)/Cell_Z), nz-1);
      const int index = z1+nz*(j1+ny*i1);
      random.catalogue_object(i)->set_region(index);
    }
  }

  cosmobl::check_regions(data, random);
}


// ============================================================================


void cosmobl::set_ObjectRegion_mangle (catalogue::Catalogue &data, catalogue::Catalogue &random, const int nSamples, const string polygonfile)
{
  string mangle_dir = fullpath(par::DirCosmo)+"/External/mangle/";

  string mangle_working_dir = mangle_dir+"output/";
  string mkdir = "mkdir -p "+mangle_working_dir;
  if (system(mkdir.c_str())) {} 

  string out_cat = mangle_working_dir+"data.dat";
  string out_ran = mangle_working_dir+"ran.dat";

  ofstream fout(out_cat.c_str()); checkIO(fout, out_cat); fout.precision(10);
  
  for (size_t i=0; i<data.nObjects(); i++)
    fout << data.ra(i) << " " << data.dec(i) << endl;
  fout.clear(); fout.close();
   
  fout.open(out_ran.c_str()); checkIO(fout, out_ran);
  for (size_t i=0; i<random.nObjects(); i++)
    fout << random.ra(i) << " " << random.dec(i) << endl;
  fout.clear(); fout.close();
  
  string cmd = mangle_dir+"bin/polyid -ur "+polygonfile+" "+out_cat+" "+out_cat+".id";
  if (system(cmd.c_str())) {}
  cmd = mangle_dir+"bin/polyid -ur "+polygonfile+" "+out_ran+" "+out_ran+".id";
  if (system(cmd.c_str())) {}
   
  vector<int> poly_data, poly_random, poly_list;

  string line;
  string in_cat = out_cat+".id";
  string in_ran = out_ran+".id";

  ifstream fin(in_cat.c_str()); checkIO(fin, in_cat); 
  getline(fin, line);
  while (getline(fin, line)) {
    stringstream ss(line); double NUM; int pp=-100;
    ss >> NUM; 
    ss >> NUM; 
    ss >> pp;
    if (pp==-100) ErrorCBL("Error in cosmobl::set_ObjectRegion_mangle!");
    poly_data.push_back(pp);
  }
  fin.clear(); fin.close();

  fin.open(in_ran.c_str()); checkIO(fin, in_ran); 
  getline(fin, line);
  while (getline(fin, line)) {
    stringstream ss(line); double NUM; int pp = -100;
    ss >> NUM; 
    ss >> NUM; 
    ss >> pp;
    if (pp==-100) ErrorCBL("Error in cosmobl::set_ObjectRegion_mangle!");
    poly_random.push_back(pp);
    poly_list.push_back(pp);
  }
  fin.clear(); fin.close();

  vector<int>::iterator it = poly_list.begin();
  sort(poly_list.begin(), poly_list.end());
  it = unique(poly_list.begin(), poly_list.end());
  poly_list.resize(distance(poly_list.begin(), it));

  int nPoly = poly_list.size();

  vector<int> boundaries(nSamples+1, 0);
  boundaries[0] = Min(poly_list); boundaries[nSamples] = Max(poly_list)+100;

  for (int i=1; i<nSamples; i++)
    boundaries[i] = poly_list[i*int(nPoly/(nSamples))];
 
  for (size_t i=1; i<boundaries.size(); i++) {
    for (size_t j=0; j<poly_data.size(); j++) 
      if (poly_data[j]>=boundaries[i-1] && poly_data[j] <boundaries[i])
	data.catalogue_object(j)->set_region(i-1);
    
    for (size_t j=0; j<poly_random.size(); j++) 
      if (poly_random[j]>=boundaries[i-1] && poly_random[j]<boundaries[i]) 
	random.catalogue_object(j)->set_region(i-1);
  }
  
  string RM = "rm -rf "+mangle_working_dir;
  if (system(RM.c_str())) {}

  cosmobl::check_regions(data, random);
}


// ============================================================================


void cosmobl::set_ObjectRegion_RaDec (catalogue::Catalogue &data, catalogue::Catalogue &random, const double Cell_size)
{
  vector<double> Lim;
  
  double Cell_sz = radians(Cell_size);

  vector<double> data_x = data.var(catalogue::Var::_RA_);
  vector<double> data_y = data.var(catalogue::Var::_Dec_);

  vector<double> random_x = random.var(catalogue::Var::_RA_);
  vector<double> random_y = random.var(catalogue::Var::_Dec_);

  for (size_t i=0; i<data.nObjects(); i++)
    data_x[i] *= cos(data_y[i]);

  for (size_t i=0; i<random.nObjects(); i++)
    random_x[i] *= cos(random_y[i]);

  Lim.push_back(Min(data_x));
  Lim.push_back(Max(data_x));

  Lim.push_back(Min(data_y));
  Lim.push_back(Max(data_y));

  const int nDec = (Lim[3]-Lim[2])/Cell_sz;

  vector<double> cell_size_x(nDec);
  vector<int> n_cells_x(nDec);
  vector<vector<int>> cells;
  int n = 0;

  for (int i=0; i<nDec; i++) {
    
    double cd = Lim[2]+(i+0.5)*Cell_sz;
    n_cells_x[i] = ceil((Lim[1]-Lim[0])*cos(cd)/Cell_sz);
    cell_size_x[i] = (Lim[1]-Lim[0])/n_cells_x[i];

    vector<int> vv(n_cells_x[i]);
    for (int j=0; j<n_cells_x[i]; j++) {
      vv[j] = n;
      n ++;
    }
    
    cells.push_back(vv);
  }

  
#pragma omp parallel num_threads(omp_get_max_threads())
  {
    
#pragma omp for schedule(static, 2) 
    for (size_t i=0; i<data.nObjects(); i++) {
      int j1 = min(int((data_y[i]-Lim[2])/Cell_sz), nDec-1);
      int i1 = min(int((data_x[i]-Lim[0])/cell_size_x[j1]), n_cells_x[j1]-1);
      data.catalogue_object(i)->set_region(cells[j1][i1]);
    }

#pragma omp for schedule(static, 2) 
    for (size_t i=0; i<random.nObjects(); i++) {
      int j1 = min(int((random_y[i]-Lim[2])/Cell_sz), nDec-1);
      int i1 = min(int((random_x[i]-Lim[0])/cell_size_x[j1]), n_cells_x[j1]-1);
      random.catalogue_object(i)->set_region(cells[j1][i1]);
    }
  }

  cosmobl::check_regions(data, random);
}


// ============================================================================


void cosmobl::check_regions (catalogue::Catalogue &data, catalogue::Catalogue &random)
{
  coutCBL << "Checking if the regions have been assigned correctly..." << endl;
  
  vector<long> data_regions = data.region_list();
  vector<long> random_regions = random.region_list();
  
  // check if data and random catalogues have the same number of regions
  if (data_regions.size() != random_regions.size()) 
    ErrorCBL("Error in check_regions of Subsample.cpp, data and random have different number of regions: data_regions.size() = "+conv(data_regions.size(), par::fINT)+", random_regions.size = "+conv(random_regions.size(), par::fINT));
  
  // check if data and random catalogues have the same regions
  size_t nRegions = data_regions.size();
  for (size_t i=0; i<nRegions; ++i)
    if (data_regions[i] != random_regions[i])
      ErrorCBL("Error in check_regions of Subsample.cpp, data and random have regions");
  
  // check if regions are consecutive and starting from 0
  bool cons = true;
  for (size_t i=0; i<nRegions; ++i) 
    if (data_regions[i] != (long)i)
      cons = false;
  
  // make regions 
  if (!cons) {
    map<long, long> regions;
    for (size_t i=0; i<nRegions; ++i)
      regions[data_regions[i]] = (long)i;
      
    for (size_t i=0; i<data.nObjects(); ++i) {
      long region = data.catalogue_object(i)->region();
      auto search = regions.find(region);
      data.catalogue_object(i)->set_region(search->second);
    }
      
    for (size_t i=0; i<random.nObjects(); ++i) {
      long region = random.catalogue_object(i)->region();
      auto search = regions.find(region);
      random.catalogue_object(i)->set_region(search->second);
    } 
  }

}
