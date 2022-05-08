/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Catalogue/RandomCatalogueVIPERS.cpp
 *
 *  @brief Methods of the class Catalogue to construct random
 *  catalogues for VIPERS
 *
 *  This file contains the implementation of the methods of the class
 *  Catalogue used to create random catalogues for VIPERS
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unibo.it
 */

#include "Catalogue.h"

using namespace std;

using namespace cbl;


// ============================================================================


/// @cond extrandom

cbl::catalogue::Catalogue::Catalogue (const RandomType type, const string WField, const bool isSpectroscopic, const Catalogue catalogue, const Catalogue catalogue_for_nz, const double N_R, const cosmology::Cosmology &cosm, const int step_redshift, const vector<double> lim, const double redshift_min, const double redshift_max, const bool do_convol, const double sigma, const bool use_venice, const bool do_zdistr_with_venice, const string file_random, const string mask, const string pointing_file, const string dir_venice, const int seed) 
{
  if (type!=RandomType::_createRandom_VIPERS_) ErrorCBL("the random catalogue has to be of type _VIPERS_ !", "Catalogue", "RandomCatalogueVIPERS.cpp");
  
  coutCBL << par::col_green << "I'm creating the random catalogue..." << par::col_default << endl;

  
  // ----- compute the redshift distribution of the input (spectroscopic) catalogue ----- 
  
  vector<double> xx, yy, err;

  string file_nz = par::DirLoc+"file_nz";
  
  if (isSpectroscopic) {
    vector<double> redshift = catalogue_for_nz.var(Var::_Redshift_); 
    vector<double> weight = catalogue_for_nz.var(Var::_Weight_); 
    const double weightedN = catalogue_for_nz.weightedN();  
    distribution(xx, yy, err, redshift, weight, step_redshift, true, file_nz, weightedN, cbl::Min(redshift), cbl::Max(redshift), "Linear", do_convol, sigma);
  }
  
  
  // ----- read the pointings ----- 
  
  vector<string> field, pointing;

  char WW = (WField=="W1") ? '1' : '4';
  
  if (isSpectroscopic) {
    
    string file_pointings = par::DirCosmo+"External/VIPERS/masks/pointings/"+pointing_file;
    ifstream fin_pointings(file_pointings); checkIO(fin_pointings, file_pointings);
    
    string line, Field, Quadrant;
    
    while (getline(fin_pointings, line)) {
      stringstream ss(line);
      ss >> Field; ss >> Quadrant;

      if (Field.at(1)==WW || WField=="W1W4") {
	field.emplace_back(Field);
	pointing.emplace_back(Field+Quadrant);
      }

    }
    
    fin_pointings.clear(); fin_pointings.close();
  }
  

  // ----- run venice ----- 
  
  vector<string> file_input;
  
  if (use_venice) {
    
    vector<string> where = {"inside"};
    if (!isSpectroscopic) where.push_back("outside");

    const int nRandom = catalogue.nObjects()*N_R*10.;
    const string Dir_venice = par::DirCosmo+"External/VIPERS/"+dir_venice+"/";

    // compile venice
    const string MAKE = /*"make clean -C "+Dir_venice+" && "+*/ "make CC=gcc -C "+Dir_venice+" ; ";
    coutCBL << endl << "--> " << MAKE << endl << endl;
    if (system(MAKE.c_str())) {}
    
    if (!isSpectroscopic) { 
      for (size_t i=0; i<where.size(); ++i) {
	file_input.emplace_back(par::DirLoc+"temp"+conv(i, par::fINT));
	string DO = Dir_venice+"venice -m "+mask+" -r -f "+where[i]+" -npart "+conv(nRandom, par::fINT)+" -o "+par::DirLoc+"temp"+conv(i, par::fINT);
	if (do_zdistr_with_venice && isSpectroscopic) DO += " -nz "+file_nz;
	coutCBL << endl << "--> " << MAKE << endl << endl;
	if (system(DO.c_str())) {}
      }
    }
    
    else {
      
      for (size_t i=0; i<where.size(); ++i) {
      
	file_input.emplace_back(par::DirLoc+"temp"+conv(i, par::fINT));


	// set how many random galaxies put in each pointing

	random::UniformRandomNumbers_Int random(0, pointing.size(), seed);

	vector<int> rnd_pnt(pointing.size(), 0);
	
	for (int i=0; i<nRandom; ++i)
	  rnd_pnt[random()] ++;
	
	
	// use each single pointing (separately)
      
	string filelist, file, Mask;
  
	for (size_t mm=0; mm<pointing.size(); ++mm) {
	
	  file = par::DirLoc+"temp"+conv(i, par::fINT)+"_"+conv(mm, par::fINT);

	  Mask = mask+"mask_"+pointing[mm]+".reg";
	
	  string DO = par::DirCosmo+"External/VIPERS/"+dir_venice+"/venice -m "+Mask+" -r -f "+where[i]+" -npart "+conv(rnd_pnt[mm], par::fINT)+" -o "+file;
	  if (do_zdistr_with_venice && isSpectroscopic) DO += " -nz "+file_nz;

	  //coutCBL << endl << "--> " << DO << endl << endl;
	  if (system(DO.c_str())) {}
	
	  filelist += file+" ";

	
	  // add the pointing names

	  ifstream fin(file); checkIO(fin, file);
	  string line; vector<string> Line;
	  while (getline(fin, line)) 
	    if (line.find("#") == string::npos) 
	      Line.emplace_back(line);
	  fin.clear(); fin.close();
	
	  ofstream fout(file); checkIO(fout, file);
	  for (auto &&LL : Line)
	    fout << LL << "   " << pointing[mm] << endl;
	  fout.clear(); fout.close();
	
	}

	string Merge = "cat "+filelist+" > "+par::DirLoc+"temp"+conv(i, par::fINT)+"; rm -f "+par::DirLoc+"temp*_*";
	if (system(Merge.c_str())) {}
      
      }
    }
  }
  
  else
    file_input.emplace_back(file_random);


  // ----- read the random catalogues created by venice ----- 
  
  for (size_t ff=0; ff<file_input.size(); ++ff) {
    
    coutCBL << "reading the field " << file_input[ff] << "..." << endl;
    
    ifstream fin(file_input[ff].c_str()); checkIO(fin, file_input[ff]);
    string line;
    getline(fin, line);
    
    double RA, DEC, REDSHIFT;
    string FIELD;
    
    if (do_zdistr_with_venice && isSpectroscopic) { // use venice to extract the redshifts of the random objects
      while (getline(fin, line)) 
	if (line.find("#") == string::npos) { // skip a comment
	  stringstream ss(line);
	  ss >> RA; ss >> DEC; ss >> REDSHIFT; ss >> FIELD;
	  if (redshift_min<REDSHIFT && REDSHIFT<redshift_max && lim[0]<RA && RA<lim[1] && lim[2]<DEC && DEC<lim[3]) {
	    observedCoordinates coord = {RA, DEC, REDSHIFT};
	    m_object.push_back(move(Object::Create(ObjectType::_Random_, coord, CoordinateUnits::_degrees_, cosm, 1., 0, -1, FIELD)));
	  }
	}
    }
    
    else { // extract the redshifts of the random objects directly
	  
      vector<double> random_ra, random_dec, random_redshift;
      vector<string> field;

      // read R.A., Dec, and the observed field (used for jackknife and bootstrap)
      while (getline(fin, line)) 
	if (line.find("#") == string::npos) { // skip a comment
	  stringstream ss(line);
	  ss >> RA; ss >> DEC; ss >> FIELD;
	  if (lim[0]<RA && RA<lim[1] && lim[2]<DEC && DEC<lim[3]) {
	    random_ra.push_back(RA);
	    random_dec.push_back(DEC);
	    field.emplace_back(FIELD);
	  }
	}

      if (isSpectroscopic) {
	// extract the random redshifts from the distribution
	vector<double> redshift = catalogue_for_nz.var(Var::_Redshift_); 
	double zmin = max(cbl::Min(redshift), redshift_min);
	double zmax = min(cbl::Max(redshift), redshift_max);
	random_redshift = vector_from_distribution(random_ra.size(), xx, yy, zmin, zmax, seed);
      }
      else
	for (size_t i=0; i<random_ra.size(); ++i)
	  random_redshift.emplace_back(1.);
      
      // construct the objects
      for (size_t i=0; i<random_ra.size(); ++i) {
	observedCoordinates coord = {random_ra[i], random_dec[i], random_redshift[i]};
	m_object.push_back(move(Object::Create(ObjectType::_Random_, coord, CoordinateUnits::_degrees_, cosm, 1., 0, -1, field[i])));
      }
      
    }
    
    fin.clear(); fin.close();
  }
  
 
  if (isSpectroscopic) 
    computeComovingCoordinates(cosm);
  
  
  // ----- resize the random catalogue to have exactly N_R*catalogue.nObjects() objects -----

  size_t nFieldsB = nFields();
  
  std::shuffle(begin(m_object), end(m_object), default_random_engine{});
  m_object.resize(N_R*catalogue.nObjects());
  

  // ----- check the number of fields -----
  
  if (nFields() != nFieldsB)
    ErrorCBL("nFields() = "+conv(nFields(), par::fINT)+" != nFieldsB = "+conv(nFieldsB, par::fINT)+" --> the random sample should be probably enlarged", "Catalogue", "RandomCatalogueVIPERS.cpp");

  if (catalogue.nFields()>nFields())
    ErrorCBL("catalogue.nFields() = "+conv(catalogue.nFields(), par::fINT)+" > nFields = "+conv(nFields(), par::fINT)+" --> the random sample should be probably enlarged", "Catalogue", "RandomCatalogueVIPERS.cpp");
  
  coutCBL <<"total number of fields in the catalogue = "<<catalogue.nFields() << endl;
  coutCBL <<"total number of fields in the random = "<<nFields()<<endl;


  // ----- remove temporary files -----
  
  string RM = "rm -rf "+par::DirLoc+"temp* "+file_nz;
  if (system(RM.c_str())) {}
  
}

/// @endcond
