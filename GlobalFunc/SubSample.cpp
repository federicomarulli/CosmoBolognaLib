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
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "GlobalFunc.h"

using namespace std;

using namespace cbl;


// ============================================================================


void cbl::set_ObjectRegion_Tiles_Redshift (catalogue::Catalogue &data, catalogue::Catalogue &random, const int nz)
{  
  coutCBL << "I'm putting data and random objects in regions given by R.A.-Dec tiles." << endl;

  if (nz <= 0)
    ErrorCBL("nz must be >0.", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");

  // Check if the necessary quantities are properly set 
  for (size_t i=0; i<data.nObjects(); i++) {

    if (data.isSetVar(i, catalogue::Var::_Region_) == false)
      ErrorCBL("The tile number for the object "+cbl::conv(i,cbl::par::fINT)+" in the data catalogue is not set.", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");
    if (data.var(i, catalogue::Var::_Region_) < 0)
      ErrorCBL("The tile number for the object "+cbl::conv(i,cbl::par::fINT)+" in the data catalogue is <0. The tile numbers must be all the integers between 0 and N, where N is the highest tile number.", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");

  }
  
  for (size_t i=0; i<random.nObjects(); i++) {

    if (random.isSetVar(i, catalogue::Var::_Region_) == false)
      ErrorCBL("The tile number for the object "+cbl::conv(i,cbl::par::fINT)+" in the random catalogue is not set.", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");
    if (random.var(i, catalogue::Var::_Region_) < 0)
      ErrorCBL("The tile number for the object "+cbl::conv(i,cbl::par::fINT)+" in the random catalogue is <0. The tile numbers must be all the integers between 0 and N, where N is the highest tile number.", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");

  }

  // Re-set the region numbers, so that they become
  // all the integers between 0 and N, where N is the
  // highest number among the regions

  std::vector<long int> observedRegions (data.nObjects(), -1);
  std::vector<long int> randomRegions (random.nObjects(), -1);

  std::vector<bool> ok_data_tile (data.nObjects(), false);
  std::vector<bool> ok_random_tile (random.nObjects(), false);

  long int value = 0;
  
  for (size_t i=0; i<data.nObjects(); i++) {

    bool update_value = false;

    for (size_t j=0; j<data.nObjects(); j++)
      if (data.region(j) == data.region(i) && ok_data_tile[j] == false) {
	observedRegions[j] = value;
	ok_data_tile[j] = true;

	update_value = true;
      }
    
    for (size_t j=0; j<random.nObjects(); j++)
      if (random.region(j) == data.region(i) && ok_random_tile[j] == false) {
	randomRegions[j] = value;
	ok_random_tile[j] = true;
      }

    if (update_value)
      value ++;
    
  }

  for (size_t i=0; i<random.nObjects(); i++)
    
    if (ok_random_tile[i] == false) {
      
      const long int new_value = cbl::Max(randomRegions)+1;
      
      for (size_t j=0; j<random.nObjects(); j++)
	if (random.region(j) == random.region(i)) {
	  randomRegions[j] = new_value;
	  ok_random_tile[j] = true;
	}
    }

  // If the number of random regions is higher
  // than the number of data regions, remove the
  // random objects in excess
  
  if (cbl::Max(randomRegions) > cbl::Max(observedRegions)) {

    const int original_nObj = random.nObjects();
    
    for (size_t i=0; i<random.nObjects(); i++) {

      const int idx = original_nObj-1-i;
      
      if ((int)(randomRegions[idx]) > (int)(cbl::Max(observedRegions))) {
	random.remove_object(idx);
	randomRegions.erase(randomRegions.begin()+idx);
      }
      
    }
  }

  // Divide the samples in redshift sub-regions, by separating
  // the objects in the same tile but in different redshift regions.

  if (nz > 1) {
  
    const double zMin = data.Min(catalogue::Var::_Redshift_);
    const double Cell_z = (data.Max(catalogue::Var::_Redshift_)-zMin)/nz;

    int next_max_tile_number = cbl::Max(observedRegions)+1;
    std::vector<std::vector<long int>> tile_z_idx (cbl::Max(observedRegions)+1, std::vector<long int>(nz, -1));

    std::vector<long int> dummy_observedRegions = observedRegions;
  
    for (size_t i=0; i<observedRegions.size(); i++) {

      int z_idx = min(int((data.redshift(i)-zMin)/Cell_z), nz-1);
    
      if (z_idx != 0) {
      
	if (tile_z_idx[dummy_observedRegions[i]][z_idx] != -1)
	  observedRegions[i] = tile_z_idx[dummy_observedRegions[i]][z_idx];
      
	else {
	  observedRegions[i] = next_max_tile_number; // Given the same tile number, if the z cell is different also the region number is different.
	  tile_z_idx[dummy_observedRegions[i]][z_idx] = next_max_tile_number;
	  next_max_tile_number ++;	
	}
      
      }
    
    }

    // Reassign the regions to the random objects
    
    std::vector<long int> dummy_randomRegions = randomRegions;
  
    for (size_t i=0; i<randomRegions.size(); i++) {
      int z_idx = min(int((random.redshift(i)-zMin)/Cell_z), nz-1);
      randomRegions[i] = tile_z_idx[dummy_randomRegions[i]][z_idx];    
    }

    // It might happen that, given a tile, no objects lie in
    // the first redshift bin. In this way, indices are lost.
    // In the loop below we recover such indices, i.e. we
    // "fill the gaps" by changing the existing indices.

    std::vector<long int> sorted_regions = cbl::different_elements(observedRegions);
    std::sort(sorted_regions.begin(), sorted_regions.end());

    std::vector<long int> missed_index;
    int idx = -1;
    
    for (size_t i=0; i<sorted_regions.size(); i++) {
      idx ++;
      if ((size_t)(sorted_regions[idx]) != i) {
	missed_index.emplace_back(i);
	idx --;
      }
    }
    
    for (size_t i=0; i<missed_index.size(); i++) {

      const long int max_idx = sorted_regions[sorted_regions.size()-1];
    
      for (size_t j=0; j<observedRegions.size(); j++)
	if (observedRegions[j] == max_idx)
	  observedRegions[j] = missed_index[i];

      for (size_t j=0; j<randomRegions.size(); j++)
	if (randomRegions[j] == max_idx)
	  randomRegions[j] = missed_index[i];

      sorted_regions.erase(sorted_regions.end()-1);
      
    }
    
  }

  for (size_t i=0; i<random.nObjects(); i++)
    if (randomRegions[i] == -1)
      ErrorCBL("Some random objects fall in redshift bins that are not populated by any data object!", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");  

  // Set the regions
  const int nRegions_data = (int)((cbl::different_elements(observedRegions)).size());
  const int nRegions_random = (int)((cbl::different_elements(randomRegions)).size());
  
  data.set_region(observedRegions, nRegions_data);
  random.set_region(randomRegions, nRegions_random);

  cbl::check_regions(data, random);

  coutCBL << "Done!" << endl;
}


// ============================================================================


void cbl::set_ObjectRegion_Tiles_Redshift (catalogue::Catalogue &data, catalogue::Catalogue &random, const int nz, const double tile_width_RA, const double tile_width_Dec, const bool write_tiles, const std::string dir_tiles, const std::string file_tiles)
{
  ErrorCBL("Work in progress!", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");
  
  coutCBL << "I'm putting data and random objects in regions given by R.A.-Dec tiles." << endl;

  if (nz <= 0)
    ErrorCBL("nz must be >0.", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");

  const double RA_hw = 0.5 * tile_width_RA * (par::pi/180.); // half tile width in radians along R.A.
  const double Dec_hw = 0.5 * tile_width_Dec * (par::pi/180.);

  // Set the vector of tile numbers
  std::vector<long int> dummy_tiles = data.region();
  std::sort(dummy_tiles.begin(), dummy_tiles.end());
  std::vector<long int> unique_tile_numbers = cbl::different_elements (dummy_tiles);
  const int n_tiles = (int)(unique_tile_numbers.size());

  cout<<endl; coutCBL << "I'm setting the regions in " << n_tiles << " tiles and " << nz << " redshift sub-sample(s)." << endl;

  // Define the vectors of min/max R.A./Dec,
  // and check if the internal variables are properly set
  std::vector<double> RA_min (n_tiles);
  std::vector<double> RA_max (n_tiles);
  std::vector<double> Dec_min (n_tiles);
  std::vector<double> Dec_max (n_tiles);

  std::vector<bool> isSet_region (n_tiles);
  
  for (size_t i=0; i<data.nObjects(); i++) {

    if (data.isSetVar(i, catalogue::Var::_TileRA_) == false)
      ErrorCBL("The tile central R.A. for the object "+cbl::conv(i,cbl::par::fINT)+" in the original catalogue (data) is not set.", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");
    if (data.isSetVar(i, catalogue::Var::_TileDec_) == false)
      ErrorCBL("The tile central Dec for the object "+cbl::conv(i,cbl::par::fINT)+" in the original catalogue (data) is not set.", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");
    if (data.isSetVar(i, catalogue::Var::_Region_) == false)
      ErrorCBL("The tile number for the object "+cbl::conv(i,cbl::par::fINT)+" in the original catalogue (data) is not set.", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");
    if (data.var(i, catalogue::Var::_Region_) < 0)
      ErrorCBL("The tile number for the object "+cbl::conv(i,cbl::par::fINT)+" in the original catalogue (data) is <0. The tile numbers must be all the integers between 0 and N, where N is the highest tile number.", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");

    if (data.region(i) < n_tiles)
      isSet_region[data.region(i)] = true;
    else
      ErrorCBL("The tile number "+cbl::conv(i,cbl::par::fINT)+" in the original catalogue (data) is greater than the number of tiles. The tile numbers must be all the integers between 0 and N, where N is the highest tile number.", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp"); 

    // Set min/max R.A./Dec
    double candidate_RA_min = data.ra_tile(i) - RA_hw/cos(std::abs(data.dec_tile(i)));
    if (candidate_RA_min < 0)
      candidate_RA_min = 2*cbl::par::pi + candidate_RA_min;

    double candidate_RA_max = data.ra_tile(i) + RA_hw/cos(std::abs(data.dec_tile(i)));
    if (candidate_RA_max > 2*cbl::par::pi)
      candidate_RA_max = candidate_RA_max - 2*cbl::par::pi;
    
    double candidate_Dec_min = data.dec_tile(i) - Dec_hw;
    //if (candidate_Dec_min < - 0.5*cbl::par::pi)
    //  candidate_Dec_min = - 0.5*cbl::par::pi + (candidate_Dec_min - 0.5*cbl::par::pi);
    
    double candidate_Dec_max = data.dec_tile(i) + Dec_hw;
    //if (candidate_Dec_max > 0.5*cbl::par::pi)
    //  candidate_Dec_max = 0.5*cbl::par::pi - (candidate_Dec_max - 0.5*cbl::par::pi);
    
    RA_min[data.region(i)] = candidate_RA_min;
    RA_max[data.region(i)] = candidate_RA_max;
    Dec_min[data.region(i)] = candidate_Dec_min;
    Dec_max[data.region(i)] = candidate_Dec_max;

  }

  for (int i=0; i<n_tiles; i++)
    if (isSet_region[i] == false)
      ErrorCBL("The tile number "+cbl::conv(i,cbl::par::fINT)+" in the original catalogue (data) is missing. The tile numbers must be all the integers between 0 and N, where N is the highest tile number.", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp"); 

  std::vector<bool>().swap(isSet_region);

  // Write the tiles file
  if (write_tiles) {
    std::string mkdir = "mkdir -p "+dir_tiles; if (system(mkdir.c_str())) {}
    std::ofstream myfile; myfile.open(dir_tiles+file_tiles);
    myfile << "RA_min   RA_max   Dec_min   Dec_max   [all in degrees]" <<std::endl;
    for (int i=0; i<n_tiles; i++)
      myfile << RA_min[i]*(180./cbl::par::pi) << std::setw(20) << RA_max[i]*(180./cbl::par::pi) << std::setw(20) << Dec_min[i]*(180./cbl::par::pi) << std::setw(20) << Dec_max[i]*(180./cbl::par::pi) << std::endl;
    myfile.close();
    coutCBL<<"I wrote the file "+dir_tiles+file_tiles<<std::endl;
  }

  // Divide the observed sample in redshift sub-regions, by separating
  // the objects in the same tile but in different redshift regions.  
  const double zMin = data.Min(catalogue::Var::_Redshift_);
  const double Cell_z = (data.Max(catalogue::Var::_Redshift_)-zMin)/nz;
  
  std::vector<long int> observedRegions (data.nObjects());
  int next_max_tile_number = data.Max(catalogue::Var::_Region_) + 1;

  std::vector<std::vector<long int>> tile_z_idx (data.Max(catalogue::Var::_Region_)+1, std::vector<long int>(nz, -1));
    
  for (size_t i=0; i<data.nObjects(); i++) {   

    if (RA_max[data.region(i)] > RA_min[data.region(i)]) {
      
      if (RA_min[data.region(i)] < data.ra(i) && RA_max[data.region(i)] > data.ra(i) && Dec_min[data.region(i)] < data.dec(i) && Dec_max[data.region(i)] > data.dec(i))
	continue;
      else
	ErrorCBL("The data object "+cbl::conv(i,cbl::par::fINT)+", with R.A.="+cbl::conv(data.ra(i)*(180./cbl::par::pi),cbl::par::fDP5)+" deg and Dec="+cbl::conv(data.dec(i)*(180./cbl::par::pi),cbl::par::fDP5)+" deg, does not fall within its original tile with central R.A.="+cbl::conv(data.ra_tile(i)*(180./cbl::par::pi),cbl::par::fDP5)+" deg and central Dec="+cbl::conv(data.dec_tile(i)*(180./cbl::par::pi),cbl::par::fDP5)+" deg (the candidate tile has R.A.=["+cbl::conv(RA_min[data.region(i)]*(180./cbl::par::pi),cbl::par::fDP5)+","+cbl::conv(RA_max[data.region(i)]*(180./cbl::par::pi),cbl::par::fDP5)+"] and Dec=["+cbl::conv(Dec_min[data.region(i)]*(180./cbl::par::pi),cbl::par::fDP5)+","+cbl::conv(Dec_max[data.region(i)]*(180./cbl::par::pi),cbl::par::fDP5)+"]). Try to set another value for the tile width.", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");

    } else {
      
      if ( (RA_min[data.region(i)] < data.ra(i) || RA_max[data.region(i)] > data.ra(i)) && (Dec_min[data.region(i)] < data.dec(i) && Dec_max[data.region(i)] > data.dec(i)) )
	continue;
      else
        ErrorCBL("The data object "+cbl::conv(i,cbl::par::fINT)+", with R.A.="+cbl::conv(data.ra(i)*(180./cbl::par::pi),cbl::par::fDP5)+" deg and Dec="+cbl::conv(data.dec(i)*(180./cbl::par::pi),cbl::par::fDP5)+" deg, does not fall within its original tile with central R.A.="+cbl::conv(data.ra_tile(i)*(180./cbl::par::pi),cbl::par::fDP5)+" deg and central Dec="+cbl::conv(data.dec_tile(i)*(180./cbl::par::pi),cbl::par::fDP5)+" deg (the candidate tile has R.A.=["+cbl::conv(RA_min[data.region(i)]*(180./cbl::par::pi),cbl::par::fDP5)+","+cbl::conv(RA_max[data.region(i)]*(180./cbl::par::pi),cbl::par::fDP5)+"] and Dec=["+cbl::conv(Dec_min[data.region(i)]*(180./cbl::par::pi),cbl::par::fDP5)+","+cbl::conv(Dec_max[data.region(i)]*(180./cbl::par::pi),cbl::par::fDP5)+"]). Try to set another value for the tile width.", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");
      
    }

    int z_idx = min(int((data.redshift(i)-zMin)/Cell_z), nz-1);
    
    if (z_idx == 0)
      observedRegions[i] = data.region(i);
    
    else {
      
      if (tile_z_idx[data.region(i)][z_idx] != -1)
	observedRegions[i] = tile_z_idx[data.region(i)][z_idx];
      
      else {
	observedRegions[i] = next_max_tile_number; // Given the same tile number, if the z cell is different also the region number is different.
	tile_z_idx[data.region(i)][z_idx] = next_max_tile_number;
	next_max_tile_number ++;	
      }
      
    }
    
  }

  std::vector<std::vector<long int>>().swap(tile_z_idx);
  
  // Define the number of regions
  std::vector<long int> unique_region_numbers = cbl::different_elements (observedRegions);
  const int nRegions = unique_region_numbers.size();

  std::vector<long int>().swap(unique_region_numbers);

  // Define the random catalogue regions  
  std::vector<long int> randomRegions (random.nObjects());
  
  std::vector<bool> isSet_random_region (random.nObjects());

  next_max_tile_number = data.Max(catalogue::Var::_Region_) + 1; // Set it equal to the original value
  
  std::vector<std::vector<long int>> tile_z_idx_random (cbl::Max(unique_tile_numbers)+1, std::vector<long int>(nz, -1));
  
  for (size_t i=0; i<random.nObjects(); i++) {

    for (int j=0; j<n_tiles; j++) {

      if (Dec_min[j] < random.dec(i) && Dec_max[j] > random.dec(i)) {

	int z_idx = min(int((random.redshift(i)-zMin)/Cell_z), nz-1);

	if (RA_max[j] > RA_min[j]) { // To manage the cases in which a tile lies at the zero of R.A.

	  if (RA_min[j] < random.ra(i) && RA_max[j] > random.ra(i)) {

	    if (isSet_random_region[i])
	      WarningMsgCBL("The region of the "+cbl::conv((int)(i),cbl::par::fINT)+"-th random object, with R.A.="+cbl::conv(random.ra(i)*(180./cbl::par::pi),cbl::par::fDP5)+" deg and Dec="+cbl::conv(random.dec(i)*(180./cbl::par::pi),cbl::par::fDP5)+" deg, is already set!", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");

	    if (z_idx == 0)
	      randomRegions[i] = unique_tile_numbers[j];
    
	    else {
      
	      if (tile_z_idx_random[j][z_idx] != -1)
		randomRegions[i] = tile_z_idx_random[j][z_idx];
      
	      else {
		randomRegions[i] = next_max_tile_number;
		tile_z_idx_random[j][z_idx] = next_max_tile_number;
		next_max_tile_number ++;	
	      }      
	    }
	    isSet_random_region[i] = true;
	  }
	}
	else {

	  if (RA_min[j] < random.ra(i) || RA_max[j] > random.ra(i)) { // To manage the cases in which a tile lies at the zero of R.A.

	    if (isSet_random_region[i])
	      WarningMsgCBL("The region of the "+cbl::conv((int)(i),cbl::par::fINT)+"-th random object, with R.A.="+cbl::conv(random.ra(i)*(180./cbl::par::pi),cbl::par::fDP5)+" deg and Dec="+cbl::conv(random.dec(i)*(180./cbl::par::pi),cbl::par::fDP5)+" deg, is already set!", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");

	    if (z_idx == 0)
	      randomRegions[i] = unique_tile_numbers[j];
    
	    else {
      
	      if (tile_z_idx_random[j][z_idx] != -1)
		randomRegions[i] = tile_z_idx_random[j][z_idx];
      
	      else {
		randomRegions[i] = next_max_tile_number;
		tile_z_idx_random[j][z_idx] = next_max_tile_number;
		next_max_tile_number ++;	
	      }      
	    }
	    isSet_random_region[i] = true;
	  }
	}
      }
    }
  }

  for (size_t i=0; i<random.nObjects(); i++)
    if (isSet_random_region[i] == false)
    ErrorCBL("Region not assigned for the random object "+cbl::conv((int)(i),cbl::par::fINT)+". The random catalogue must account for the survey masks.", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");

  const int n1 = (int)(cbl::different_elements(observedRegions).size());
  const int n2 = (int)(cbl::different_elements(randomRegions).size());
  if (n1 < n2)
    ErrorCBL("", "set_ObjectRegion_Tiles_Redshift", "GlobalFunc/SubSample.cpp");

  // Set the regions
  data.set_region(observedRegions, nRegions);
  random.set_region(randomRegions, nRegions);

  cbl::check_regions(data, random);

  coutCBL << "Done!" << endl;
}


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


vector<double> cbl::colatitude (std::vector<double> latitude)
{
  vector<double> colatitude(latitude.size());
  
  for (size_t i=0; i<latitude.size(); i++)
    colatitude[i] = cbl::par::pi/2-latitude[i];

  return colatitude;
}


// ============================================================================


void cbl::set_ObjectRegion_RaDec (catalogue::Catalogue &data, const int nCells_Ra, const int nCells_Dec, const bool use_colatitude)
{
  vector<double> data_x = data.var(catalogue::Var::_RA_);
  vector<double> data_y = (use_colatitude) ? cbl::colatitude(data.var(catalogue::Var::_Dec_)) : data.var(catalogue::Var::_Dec_);
  vector<double> cos_data_y(data.nObjects(), 0);
  for (size_t i=0; i<data.nObjects(); i++)
    cos_data_y[i] = cos(data_y[i]);

  double min_ra = Min(data_x);
  double max_ra = Max(data_x);
  double deltaRa = (max_ra-min_ra)/nCells_Ra;

  double min_cdec = Min(cos_data_y);
  double max_cdec = Max(cos_data_y);
  double deltaCDec = (max_cdec-min_cdec)/nCells_Dec;

  int nCells = nCells_Ra*nCells_Dec;
  double Area = deltaRa*deltaCDec*nCells;

  coutCBL << "Survey area is: " << Area << endl;
  coutCBL << "Number of cells will be: " << nCells << endl;

  vector<long> dataReg(data.nObjects());

#pragma omp parallel num_threads(omp_get_max_threads())
  {
    
#pragma omp for schedule(static, 2) 
    for (size_t i=0; i<data.nObjects(); i++) {
      int j1 = min(int((cos_data_y[i]-min_cdec)/deltaCDec), nCells_Dec-1);
      int i1 = min(int((data_x[i]-min_ra)/deltaRa), nCells_Ra-1);
      dataReg[i] = i1*nCells_Dec+j1;
    }
  }

  data.set_region(dataReg, nCells);
}



// ============================================================================


void cbl::set_ObjectRegion_RaDec (catalogue::Catalogue &data, catalogue::Catalogue &random, const int nCells_Ra, const int nCells_Dec, const bool use_colatitude)
{
  vector<double> data_x = data.var(catalogue::Var::_RA_);
  vector<double> data_y = (use_colatitude) ? cbl::colatitude(data.var(catalogue::Var::_Dec_)) : data.var(catalogue::Var::_Dec_);
  vector<double> cos_data_y(data.nObjects(), 0);
  for (size_t i=0; i<data.nObjects(); i++) 
    cos_data_y[i] = cos(data_y[i]);
  

  vector<double> random_x = random.var(catalogue::Var::_RA_);
  vector<double> random_y = (use_colatitude) ? cbl::colatitude(random.var(catalogue::Var::_Dec_)) : random.var(catalogue::Var::_Dec_);
  vector<double> cos_random_y(random.nObjects(), 0);
  for (size_t i=0; i<random.nObjects(); i++)
    cos_random_y[i] = cos(random_y[i]);

  double min_ra = Min(random_x);
  double max_ra = Max(random_x);
  double deltaRa = (max_ra-min_ra)/nCells_Ra;

  double min_cdec = Min(cos_random_y);
  double max_cdec = Max(cos_random_y);
  double deltaCDec = (max_cdec-min_cdec)/nCells_Dec;

  int nCells = nCells_Ra*nCells_Dec;
  double Area = deltaRa*deltaCDec*nCells;

  coutCBL << "Survey area is: " << Area << endl;
  coutCBL << "Number of cells will be: " << nCells << endl;

  vector<long> dataReg(data.nObjects());
  vector<long> randReg(random.nObjects());

#pragma omp parallel num_threads(omp_get_max_threads())
  {
    
#pragma omp for schedule(static, 2) 
    for (size_t i=0; i<data.nObjects(); i++) {
      int j1 = min(int((cos_data_y[i]-min_cdec)/deltaCDec), nCells_Dec-1);
      int i1 = min(int((data_x[i]-min_ra)/deltaRa), nCells_Ra-1);
      dataReg[i] = i1*nCells_Dec+j1;
    }

#pragma omp for schedule(static, 2) 
    for (size_t i=0; i<random.nObjects(); i++) {
      int j1 = min(int((cos_random_y[i]-min_cdec)/deltaCDec), nCells_Dec-1);
      int i1 = min(int((random_x[i]-min_ra)/deltaRa), nCells_Ra-1);
      randReg[i] = i1*nCells_Dec+j1;
    }
  }

  data.set_region(dataReg, nCells);
  random.set_region(randReg, nCells);

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
