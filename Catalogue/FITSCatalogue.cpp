/********************************************************************
 *  Copyright (C) 2018 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Catalogue/FITSCatalogue.cpp
 *
 *  @brief Methods of the class Catalogue to construct catalogues from
 *  FITS files
 *
 *  This file contains the implementation of the constructors of the
 *  class Catalogue used to manage catalogue files in FITS format
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unibo.it
 */

#include "FITSwrapper.h"
#include "Catalogue.h"

using namespace std;

using namespace cbl;


// ============================================================================


cbl::catalogue::Catalogue::Catalogue (const ObjectType objectType, const CoordinateType coordinateType, const std::vector<std::string> file, const std::vector<std::string> column_names, const bool read_weights, const bool read_regions, const double nSub, const double fact, const cosmology::Cosmology &cosm, const CoordinateUnits inputUnits, const int seed)
{ 
  // parameters for random numbers used in case nSub!=1
  random::UniformRandomNumbers ran(0., 1., seed);

  
  // read the input catalogue files

  for (size_t dd=0; dd<file.size(); ++dd) {

    coutCBL << "Reading the catalogue: " << file[dd] << endl;

    // read the columns from the table searching by names
    vector<vector<double>> table = wrapper::ccfits::read_table_fits(file[dd], column_names, 1, 1.);
    
    if (((read_weights || read_regions) && table.size()<3) || ((read_weights && read_regions) && table.size()<4))
      ErrorCBL("the number of columns in the input FITS file is wrong!", "Catalogue", "FITSCatalogue.cpp"); 
    
    // include the objects in the catalogue
    
    for (size_t i=0; i<table[0].size(); ++i) {
 
      if (ran()<nSub) { // extract a subsample
	
	if (coordinateType==CoordinateType::_comoving_) { // comoving coordinates (x, y, z)
	  comovingCoordinates coord = {table[0][i]*fact, table[1][i]*fact, table[2][i]*fact};
	  m_object.push_back(move(Object::Create(objectType, coord, (read_weights) ? table[3][i] : 1., (read_regions) ? (long)table[(read_weights) ? 4 : 3][i] : 1)));
	}
	
	else if (coordinateType==CoordinateType::_observed_) { // observed coordinates (R.A., Dec, redshift)
	  if (table[2][i]*fact>0) {
	    observedCoordinates coord = {table[0][i]*fact, table[1][i]*fact, table[2][i]*fact};
	    m_object.push_back(move(Object::Create(objectType, coord, inputUnits, cosm, (read_weights) ? table[3][i] : 1., (read_regions) ? (long)table[(read_weights) ? 4 : 3][i] : 1)));
	  }
	  else WarningMsgCBL("the object "+conv(i, par::fINT)+" has z = "+conv(table[2][i]*fact, par::fDP2)+", and it will be not included in the catalogue!", "Catalogue", "FITSCatalogue.cpp");
	}

	else ErrorCBL("coordinateType is not valid!", "Catalogue", "FITSCatalogue.cpp");

      }
      
    }
  }
}

// ============================================================================


cbl::catalogue::Catalogue::Catalogue (const ObjectType objectType, const CoordinateType coordinateType, const std::vector<std::string> file, const std::vector<std::string> column_names, const std::vector<Var> attribute, const double nSub, const double fact, const cosmology::Cosmology &cosm, const CoordinateUnits inputUnits, const int seed)
{
  // preliminary check on vector sizes
  size_t nvar;
  if (attribute.size()==column_names.size()) nvar = attribute.size();
  else ErrorCBL("Column_names vector and attribute vector must have equal size!", "Catalogue", "FITSCatalogue.cpp");

  const int num_threads = (nvar>size_t(omp_get_max_threads())) ? omp_get_max_threads() : nvar;
  
  unordered_map<int, Var> varMap;
  for (size_t ii=0; ii<nvar; ii++)
    varMap.insert({ii, attribute[ii]});

  // parameters for random numbers used in case nSub!=1
  random::UniformRandomNumbers ran(0., 1., seed);
  
  // read the input catalogue files

  for (size_t dd=0; dd<file.size(); ++dd) {

    coutCBL << "Reading the catalogue: " << file[dd] << endl;

    // read the columns from the table searching by names
    vector<vector<double>> table = wrapper::ccfits::read_table_fits(file[dd], column_names, 1, 1.);
    
    // prepare default coordinates
    comovingCoordinates defaultComovingCoord = {par::defaultDouble, par::defaultDouble, par::defaultDouble};
    observedCoordinates defaultObservedCoord = {par::defaultDouble, -1., 0.1};
    
    // include the objects in the catalogue
    
    for (size_t i=0; i<table[0].size(); ++i) {

      size_t prev_nObj = nObjects();
 
      if (ran()<nSub) { // extract a subsample
	
	if (coordinateType==cbl::CoordinateType::_comoving_) 
	  m_object.push_back(move(Object::Create(objectType, defaultComovingCoord, 1.)));
	
	else if (coordinateType==cbl::CoordinateType::_observed_)
	  m_object.push_back(move(Object::Create(objectType, defaultObservedCoord, inputUnits, cosm, 1.)));

	
#pragma omp parallel num_threads(num_threads)
	{	    
#pragma omp for schedule(dynamic)
	  for (size_t i=0; i<nvar; ++i) {
	
	    for (size_t ii=prev_nObj; ii<nObjects(); ii++) {
	      double temp = ((varMap[i]==Var::_RA_) || (varMap[i]==Var::_Dec_)) ? radians(table[i][ii], inputUnits) : table[i][ii];
	      set_var(ii, varMap[i], ((varMap[i]==Var::_X_) || (varMap[i]==Var::_Y_) || (varMap[i]==Var::_Z_)) ? temp*fact : temp, cosm);
	    }
	    
	  }
	}
  
            
      }
    }
  }
}
