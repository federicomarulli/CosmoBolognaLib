/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

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


void cosmobl::set_ObjectRegion_SubBoxes (shared_ptr<Catalogue> data, shared_ptr<Catalogue> random, const int nx, const int ny, const int nz)
{
  vector<double> Lim;
  data->MinMax_var(Var::_XX_, Lim, 0);
  data->MinMax_var(Var::_YY_, Lim, 0);
  data->MinMax_var(Var::_ZZ_, Lim, 0);

  double Cell_X = (Lim[1]-Lim[0])/nx;
  double Cell_Y = (Lim[3]-Lim[2])/ny;
  double Cell_Z = (Lim[5]-Lim[4])/nz;

#pragma omp parallel num_threads(omp_get_max_threads())
  {
    
#pragma omp for schedule(static, 2) 
    for (int i=0; i<data->nObjects(); i++) {
      int i1 = min(int((data->xx(i)-Lim[0])/Cell_X), nx-1);
      int j1 = min(int((data->yy(i)-Lim[2])/Cell_Y), ny-1);
      int z1 = min(int((data->zz(i)-Lim[4])/Cell_Z), nz-1);
      int index = z1+nz*(j1+ny*i1);
      data->object(i)->set_region(index);
    }

#pragma omp for schedule(static, 2) 
    for (int i=0; i<random->nObjects(); i++) {
      int i1 = min(int((random->xx(i)-Lim[0])/Cell_X), nx-1);
      int j1 = min(int((random->yy(i)-Lim[2])/Cell_Y), ny-1);
      int z1 = min(int((random->zz(i)-Lim[4])/Cell_Z), nz-1);
      int index = z1+nz*(j1+ny*i1);
      random->object(i)->set_region(index);
    }
  }
  
}


// ============================================================================


void cosmobl::set_ObjectRegion_mangle (shared_ptr<Catalogue> data, shared_ptr<Catalogue> random, const int nSamples, const string polygonfile, const string dir)
{
  string temp_dir = dir+"temp/";
  string mangle_dir = par::DirCosmo+"/CatalogueAnalysis/RandomCatalogue/mangle/";
  string cmd = "mkdir -p "+temp_dir; if(system(cmd.c_str())){}
  
  string out_cat = temp_dir+"data";
  string out_ran = temp_dir+"ran";

  ofstream fout(out_cat.c_str()); checkIO(out_cat, 0); fout.precision(10);
  
  for (int i=0; i<data->nObjects(); i++)
    fout << data->ra(i) << " " << data->dec(i) << endl;
  fout.clear(); fout.close();
   
  fout.open(out_ran.c_str()); checkIO(out_ran, 0);
  for (int i=0; i<random->nObjects(); i++)
    fout << random->ra(i) << " " << random->dec(i) << endl;
  fout.clear(); fout.close();
  
  cmd = mangle_dir+"bin/polyid -ur "+polygonfile+" "+out_cat+" "+out_cat+".id";
  if (system(cmd.c_str())) {}
  cmd = mangle_dir+"bin/polyid -ur "+polygonfile+" "+out_ran+" "+out_ran+".id";
  if (system(cmd.c_str())) {}
   
   
  vector<int> poly_data, poly_random, poly_list;

  string line;
  string in_cat = out_cat+".id";
  string in_ran = out_ran+".id";

  ifstream fin(in_cat.c_str()); checkIO(in_cat, 1); 
  getline(fin, line);
  while (getline(fin, line)) {
    stringstream ss(line); double NUM; int pp=-100;
    ss >> NUM; 
    ss >> NUM; 
    ss >> pp;
    if (pp==-100) ErrorMsg("Error in cosmobl::set_ObjectRegion_mangle!");
    poly_data.push_back(pp);
  }
  fin.clear(); fin.close();

  fin.open(in_ran.c_str()); checkIO(in_ran, 1); 
  getline(fin, line);
  while (getline(fin, line)) {
    stringstream ss(line); double NUM; int pp = -100;
    ss >> NUM; 
    ss >> NUM; 
    ss >> pp;
    if (pp==-100) ErrorMsg("Error in cosmobl::set_ObjectRegion_mangle!");
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
	data->object(j)->set_region(i-1);
    
    for (size_t j=0; j<poly_random.size(); j++) 
      if (poly_random[j]>=boundaries[i-1] && poly_random[j]<boundaries[i]) 
	random->object(j)->set_region(i-1);
  }
  
  string RM = "rm -rf "+dir+"temp/";
  if (system(RM.c_str())) {}
}
