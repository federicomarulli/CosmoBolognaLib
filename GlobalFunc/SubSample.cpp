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

using namespace std;

using namespace cbl;


// ============================================================================


void cbl::set_ObjectRegion_SubBoxes (catalogue::Catalogue &data, const int nx, const int ny, const int nz)
{
  const double xMin = data.Min(catalogue::Var::_X_);
  const double yMin = data.Min(catalogue::Var::_Y_);
  const double zMin = data.Min(catalogue::Var::_Z_);
  
  const double Cell_X = (data.Max(catalogue::Var::_X_)-xMin)/nx;
  const double Cell_Y = (data.Max(catalogue::Var::_Y_)-yMin)/ny;
  const double Cell_Z = (data.Max(catalogue::Var::_Z_)-zMin)/nz;

  vector<long> dataReg(data.nObjects());

#pragma omp parallel num_threads(omp_get_max_threads())
  {
    
#pragma omp for schedule(static, 2) 
    for (size_t i=0; i<data.nObjects(); i++) {
      const int i1 = min(int((data.xx(i)-xMin)/Cell_X), nx-1);
      const int j1 = min(int((data.yy(i)-yMin)/Cell_Y), ny-1);
      const int z1 = min(int((data.zz(i)-zMin)/Cell_Z), nz-1);
      const int index = z1+nz*(j1+ny*i1);
      dataReg[i] = index;
    }
  }

  data.set_region(dataReg, nx*ny*nz);
}


// ============================================================================


void cbl::set_ObjectRegion_mangle (catalogue::Catalogue &data, const int nSamples, const std::string polygonfile)
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
  
  string cmd = mangle_dir+"bin/polyid -ur "+polygonfile+" "+out_cat+" "+out_cat+".id";
  if (system(cmd.c_str())) {}
   
  vector<int> poly_data, poly_list;

  string line;
  string in_cat = out_cat+".id";

  ifstream fin(in_cat.c_str()); checkIO(fin, in_cat); 
  getline(fin, line);
  while (getline(fin, line)) {
    stringstream ss(line); double NUM; int pp=-100;
    ss >> NUM; 
    ss >> NUM; 
    ss >> pp;
    if (pp==-100) ErrorCBL("", "set_ObjectRegion_mangle", "GlobalFunc/SubSample.cpp");
    poly_data.push_back(pp);
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
 
  vector<long> dataReg(data.nObjects());

  for (size_t i=1; i<boundaries.size(); i++) {
    for (size_t j=0; j<poly_data.size(); j++) 
      if (poly_data[j]>=boundaries[i-1] && poly_data[j] <boundaries[i])
	dataReg[j] = i-1;
  }

  data.set_region(dataReg, nSamples);
  
  string RM = "rm -rf "+mangle_working_dir;
  if (system(RM.c_str())) {}

}


// ============================================================================


void cbl::set_ObjectRegion_RaDec (catalogue::Catalogue &data, const double Cell_size)
{
  vector<double> Lim;
  
  double Cell_sz = radians(Cell_size);

  vector<double> data_x = data.var(catalogue::Var::_RA_);
  vector<double> data_y = data.var(catalogue::Var::_Dec_);

  for (size_t i=0; i<data.nObjects(); i++)
    data_x[i] *= cos(data_y[i]);

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

  vector<long> dataReg(data.nObjects());
  
#pragma omp parallel num_threads(omp_get_max_threads())
  {

#pragma omp for schedule(static, 2) 
    for (size_t i=0; i<data.nObjects(); i++) {
      int j1 = min(int((data_y[i]-Lim[2])/Cell_sz), nDec-1);
      int i1 = min(int((data_x[i]-Lim[0])/cell_size_x[j1]), n_cells_x[j1]-1);
      dataReg[i] = cells[j1][i1];
    }

  }

  data.set_region(dataReg, n);
}



// ============================================================================


void cbl::set_ObjectRegion_SubBoxes (catalogue::Catalogue &data, catalogue::Catalogue &random, const int nx, const int ny, const int nz)
{
  coutCBL << "I'm putting data and random objects in box-sized regions."<<endl;

  const double xMin = data.Min(catalogue::Var::_X_);
  const double yMin = data.Min(catalogue::Var::_Y_);
  const double zMin = data.Min(catalogue::Var::_Z_);
  
  const double Cell_X = (data.Max(catalogue::Var::_X_)-xMin)/nx;
  const double Cell_Y = (data.Max(catalogue::Var::_Y_)-yMin)/ny;
  const double Cell_Z = (data.Max(catalogue::Var::_Z_)-zMin)/nz;

  vector<long> dataReg(data.nObjects());
  vector<long> randReg(random.nObjects());

#pragma omp parallel num_threads(omp_get_max_threads())
  {
    
#pragma omp for schedule(static, 2) 
    for (size_t i=0; i<data.nObjects(); i++) {
      const int i1 = min(int((data.xx(i)-xMin)/Cell_X), nx-1);
      const int j1 = min(int((data.yy(i)-yMin)/Cell_Y), ny-1);
      const int z1 = min(int((data.zz(i)-zMin)/Cell_Z), nz-1);
      const int index = z1+nz*(j1+ny*i1);
      dataReg[i] = index;
    }

#pragma omp for schedule(static, 2) 
    for (size_t i=0; i<random.nObjects(); i++) {
      const int i1 = min(int((random.xx(i)-xMin)/Cell_X), nx-1);
      const int j1 = min(int((random.yy(i)-yMin)/Cell_Y), ny-1);
      const int z1 = min(int((random.zz(i)-zMin)/Cell_Z), nz-1);
      const int index = z1+nz*(j1+ny*i1);
      randReg[i] = index;
    }
  }

  data.set_region(dataReg, nx*ny*nz);
  random.set_region(randReg, nx*ny*nz);

  cbl::check_regions(data, random);

  coutCBL << "Done!" << endl;
}


// ============================================================================


void cbl::set_ObjectRegion_mangle (catalogue::Catalogue &data, catalogue::Catalogue &random, const int nSamples, const std::string polygonfile)
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
    if (pp==-100) ErrorCBL("", "set_ObjectRegion_mangle", "GlobalFunc/SubSample.cpp");
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
    if (pp==-100) ErrorCBL("", "set_ObjectRegion_mangle", "GlobalFunc/SubSample.cpp");
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

  vector<long> dataReg(data.nObjects());
  vector<long> randReg(random.nObjects());
 
  for (size_t i=1; i<boundaries.size(); i++) {
    for (size_t j=0; j<poly_data.size(); j++) 
      if (poly_data[j]>=boundaries[i-1] && poly_data[j] <boundaries[i])
	dataReg[j] = i-1;
    
    for (size_t j=0; j<poly_random.size(); j++) 
      if (poly_random[j]>=boundaries[i-1] && poly_random[j]<boundaries[i]) 
	randReg[j] = i-1;
  }
  
  string RM = "rm -rf "+mangle_working_dir;
  if (system(RM.c_str())) {}

  data.set_region(dataReg, nSamples);
  random.set_region(randReg, nSamples);

  cbl::check_regions(data, random);
}


// ============================================================================


void cbl::set_ObjectRegion_RaDec (catalogue::Catalogue &data, catalogue::Catalogue &random, const double Cell_size)
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

  vector<long> dataReg(data.nObjects());
  vector<long> randReg(random.nObjects());

#pragma omp parallel num_threads(omp_get_max_threads())
  {
    
#pragma omp for schedule(static, 2) 
    for (size_t i=0; i<data.nObjects(); i++) {
      int j1 = min(int((data_y[i]-Lim[2])/Cell_sz), nDec-1);
      int i1 = min(int((data_x[i]-Lim[0])/cell_size_x[j1]), n_cells_x[j1]-1);
      dataReg[i] = cells[j1][i1];
    }

#pragma omp for schedule(static, 2) 
    for (size_t i=0; i<random.nObjects(); i++) {
      int j1 = min(int((random_y[i]-Lim[2])/Cell_sz), nDec-1);
      int i1 = min(int((random_x[i]-Lim[0])/cell_size_x[j1]), n_cells_x[j1]-1);
      randReg[i] = cells[j1][i1];
    }
  }

  data.set_region(dataReg, n);
  random.set_region(randReg, n);

  cbl::check_regions(data, random);
}


// ============================================================================


void cbl::check_regions (catalogue::Catalogue &data, catalogue::Catalogue &random)
{
  coutCBL << "Checking if the regions have been assigned correctly..." << endl;

  // check if data and random catalogues have the same number of regions
  // nRegions is a "proxy" for the region geometry
  if (data.nRegions()!=random.nRegions()) 
    ErrorCBL("data and random have different number of regions: data_regions.size() = "+conv(data.nRegions(), par::fINT)+", random_regions.size = "+conv(random.nRegions(), par::fINT)+": please set the number of regions through set_region_number() or check your inputs!", "check_regions", "GlobalFunc/SubSample.cpp"); 

  size_t nRegions = data.nRegions();

  // count how many objects fall in the regions
  vector<long> dataObj_region(nRegions, 0);
  vector<long> data_region = data.region();
  
  for (size_t i=0; i<data.nObjects(); ++i) {
    checkDim(dataObj_region, data_region[i], "dataObj_region", false); 
    dataObj_region[data_region[i]] ++;
  }
  
  vector<long> randObj_region(nRegions, 0);
  vector<long> rand_region = random.region();
  for (size_t i=0; i<random.nObjects(); ++i) {
    checkDim(randObj_region, rand_region[i], "randObj_region", false); 
    randObj_region[rand_region[i]] ++;
  }

  // empty random regions are not allowed!
  // rearrange regions
  size_t region_eff = 0;
  std::map<long, long> region_list_eff; 
  for (size_t i=0; i<nRegions; i++) {
    if (randObj_region[i]!=0) {
      region_list_eff.insert(pair<long, long>(i, region_eff));
      region_eff ++;
    }
    else {
      // throw error if random region is empty and data region is not
      if (dataObj_region[i]!=0)
	ErrorCBL("the "+conv(i, par::fINT)+" region is empty in the random sample, but contains "+conv(dataObj_region[i], par::fINT)+" data points; this is not allowed, Please check your inputs!", "check_regions", "GlobalFunc/SubSample.cpp");
    }
  }

  if (region_eff!=nRegions) {
    coutCBL << "Found "+conv(region_eff, par::fINT)+" non-empty regions" << endl;
    coutCBL << "Rearranging regions in data and random sample..." << endl;
    for (size_t i=0; i<data.nObjects(); ++i) {
      checkDim(data_region, i, "data_region", false);
      if (int(region_list_eff.size())<data_region[i]) ErrorCBL("the dimension of region_list_eff is: " + conv(region_list_eff.size(), par::fINT) + " ( < " + conv(data_region[i], par::fINT) + " )", "check_regions", "SubSample.cpp"); 
      data_region[i] = region_list_eff[data_region[i]];
    }
    
    for (size_t i=0; i<random.nObjects(); ++i) {
      checkDim(rand_region, i, "rand_region", false);
      if (int(region_list_eff.size())<rand_region[i]) ErrorCBL("the dimension of region_list_eff is: " + conv(region_list_eff.size(), par::fINT) + " ( < " + conv(rand_region[i], par::fINT) + " )", "check_regions", "SubSample.cpp"); 
      rand_region[i] = region_list_eff[rand_region[i]];
    }
    
    data.set_region(data_region, region_eff+1);
    random.set_region(rand_region, region_eff+1);
    coutCBL << "Done!" << endl;
  }

  coutCBL << "End check regions!" << endl;

}


// ============================================================================


void cbl::set_ObjectRegion_SDSS_stripes (catalogue::Catalogue &data, catalogue::Catalogue &random)
{
  vector<double> lambda, eta, random_lambda, random_eta;
  vector<int> stripe, random_stripe, str_u, random_str_u;

  eq2sdss(data.var(catalogue::Var::_RA_), data.var(catalogue::Var::_Dec_), lambda, eta);
  sdss_stripe (eta, lambda, stripe, str_u);
  
  eq2sdss(random.var(catalogue::Var::_RA_), random.var(catalogue::Var::_Dec_), random_lambda, random_eta);
  sdss_stripe (random_eta, random_lambda, random_stripe, random_str_u);

  if (!isDimEqual(str_u, random_str_u))
    ErrorCBL("data and random catalogues have different stripes!", "set_ObjectRegion_SDSS_stripes", "GlobalFunc/SubSample.cpp");

  for (size_t i=0; i<data.nObjects(); i++)
    data.set_var(i, catalogue::Var::_Region_, stripe[i]);

  for (size_t i=0; i<random.nObjects(); i++)
    random.set_var(i, catalogue::Var::_Region_, random_stripe[i]);

  data.set_region_number(str_u.size());
  random.set_region_number(random_str_u.size());
} 
