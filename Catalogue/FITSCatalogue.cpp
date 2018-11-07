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
 *  @author federico.marulli3@unbo.it
 */

#include "FITSwrapper.h"
#include "Catalogue.h"

using namespace std;

using namespace cbl;


// ============================================================================


cbl::catalogue::Catalogue::Catalogue (const ObjectType objectType, const CoordinateType coordinateType, const std::vector<std::string> file, const std::vector<std::string> column_names, const int indCoord1, const int indCoord2, const int indCoord3, const int indWeight, const int indRegion, const int next, const double fill_value, const double nSub, const double fact, const cosmology::Cosmology &cosm, const CoordinateUnits inputUnits, const int seed)
{
  // parameters for random numbers used in case nSub!=1
  random::UniformRandomNumbers ran(0., 1., seed);

  
  // read the input catalogue files

  for (size_t dd=0; dd<file.size(); ++dd) {
    
    // read the columns from the table searching by names
    vector<vector<double>> table = ccfitswrapper::read_table_fits(file[dd], column_names, next, fill_value);

    
    // include the objects in the catalogue
    
    for (size_t i=0; i<table[0].size(); ++i) {

      if (ran()<nSub) { // extract a subsample
	
	if (coordinateType==CoordinateType::_comoving_) { // comoving coordinates (x, y, z)
	  comovingCoordinates coord = {table[indCoord1][i]*fact, table[indCoord2][i]*fact, table[indCoord3][i]*fact};
	  m_object.push_back(move(Object::Create(objectType, coord, (indWeight>-1) ? table[indWeight][i] : 1., (indRegion>-1) ? (long)table[indRegion][i] : 1)));
	}
	else if (coordinateType==CoordinateType::_observed_) { // observed coordinates (R.A., Dec, redshift)
	  observedCoordinates coord = {table[indCoord1][i]*fact, table[indCoord2][i]*fact, table[indCoord3][i]*fact};
	  m_object.push_back(move(Object::Create(objectType, coord, inputUnits, cosm, (indWeight>-1) ? table[indWeight][i] : 1., (indRegion>-1) ? (long)table[indRegion][i] : 1)));
	}
	else ErrorCBL("Error in cbl::catalogue::Catalogue::Catalogue() in FITSCatalogue.cpp: coordinateType is not valid!");

      }
      
    }
  }
}
