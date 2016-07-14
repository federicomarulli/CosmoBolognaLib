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
 *  @file Catalogue/RandomCatalogue.cpp
 *
 *  @brief Methods of the class Catalogue to construct random catalogues
 *
 *  This file contains the implementation of the methods of the class
 *  Catalogue, used to create random catalogues
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Catalogue.h"

using namespace cosmobl;


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const RandomType type, const Catalogue catalogue, const double N_R)
{
  if (type!=_createRandom_box_) ErrorMsg("Error in cosmobl::catalogue::Catalogue::Catalogue : the random catalogue has to be cubic!");
  
  cout << "I'm creating a random catalogue with cubic geometry..." << endl;

  int nRandom = int(N_R*catalogue.nObjects());

  double Xmin = catalogue.Min(Var::_X_), Xmax = catalogue.Max(Var::_X_); 
  double Ymin = catalogue.Min(Var::_Y_), Ymax = catalogue.Max(Var::_Y_);
  double Zmin = catalogue.Min(Var::_Z_), Zmax = catalogue.Max(Var::_Z_);

  if (Xmin>Xmax || Ymin>Ymax || Zmin>Zmax) ErrorMsg("Error in Catalogue::Catalogue of Catalogue.cpp!");

  default_random_engine gen;
  uniform_real_distribution<float> ran(0., 1.);

  for (int i=0; i<nRandom; i++) {
    comovingCoordinates coord;
    coord.xx = ran(gen)*(Xmax-Xmin)+Xmin;
    coord.yy = ran(gen)*(Ymax-Ymin)+Ymin;
    coord.zz = ran(gen)*(Zmax-Zmin)+Zmin;
    m_sample.push_back(move(Object::Create(_RandomObject_, coord)));
  }

}


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const RandomType type, const cosmology::Cosmology &real_cosm, const cosmology::Cosmology &test_cosm, const string dir_in, const double Zguess_min, const double Zguess_max)
{
  if (type!=_createRandom_box_) ErrorMsg("Error in cosmobl::catalogue::Catalogue::Catalogue : the random catalogue has to be cubic!");

  cout << "I'm creating a random catalogue with warped cubic geometry, due to geometric distortions..." << endl;

  string file_in = dir_in+"random_catalogue";
  cout << "file_in = " << file_in << endl;
  ifstream fin (file_in.c_str()); checkIO (file_in,1);
  
  string line; double NUM;
  vector<double> random_X, random_Y, random_Z;
  
  while (getline(fin, line)) {
    stringstream ss(line);
    ss >> NUM; random_X.push_back(NUM);
    ss >> NUM; random_Y.push_back(NUM);
    ss >> NUM; random_Z.push_back(NUM);
  }
  
  fin.clear(); fin.close();


  // ------ warp the random box introducing geometrical distortios ------ 
    
  vector<double> random_ra(random_X.size()), random_dec(random_X.size()), random_dc(random_X.size()), random_dc_warped(random_X.size()), random_redshift(random_X.size()); 

  polar_coord(random_X, random_Y, random_Z, random_ra, random_dec, random_dc);

  for (size_t i=0; i<random_X.size(); i++) random_redshift[i] = real_cosm.Redshift(random_dc[i], Zguess_min, Zguess_max);
 
  for (size_t i=0; i<random_X.size(); i++) random_dc_warped[i] = test_cosm.D_C(random_redshift[i]);
    
  cartesian_coord(random_ra, random_dec, random_dc_warped, random_X, random_Y, random_Z);
    
  
  // ------ create the new random catalogue ------

  for (size_t i=0; i<random_X.size(); i++) {
    comovingCoordinates coord = {random_X[i], random_Y[i], random_Z[i]};
    m_sample.push_back(move(Object::Create(_RandomObject_, coord)));
  }
  
}


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const RandomType type, const Catalogue catalogue, const double N_R, const int nbin, const cosmology::Cosmology &cosm, const bool conv, const double sigma, const int seed)
{
  vector<double> ra = catalogue.var(_RA_), raB = ra;
  vector<double> dec = catalogue.var(_Dec_), decB = dec;

  if (type==_createRandom_shuffle_) {

    cout << "I'm creating a random catalogue with the 'shuffle' method..." << endl;
    
    // the R.A.-Dec coordinates of the random catalogue will be the same as the ones of the real catalogue
    for (int i=0; i<max(0, int(N_R)-1); i++) { // use N_R more random objects
      ra.insert(ra.end(), raB.begin(), raB.end());
      dec.insert(dec.end(), decB.begin(), decB.end());
    }

  }
  
  else if (type==_createRandom_square_) {
    cout << "I'm creating a random catalogue in a RA-Dec square..." << endl;

    int nRandom = int(catalogue.nObjects()*N_R);

    ra.erase(ra.begin(), ra.end());
    dec.erase(dec.begin(), dec.end());
        
    default_random_engine generator(seed);
    uniform_real_distribution<double> distribution;
    auto ran = bind(distribution, generator);

    double ra_min = catalogue.Min(Var::_RA_),
      ra_max = catalogue.Max(Var::_RA_),
      sin_dec_min = sin(catalogue.Min(Var::_Dec_)),
      sin_dec_max = sin(catalogue.Max(Var::_Dec_));

    for (int i=0; i<nRandom; i++) {
      ra.push_back((ra_max-ra_min)*ran()+ra_min);
      dec.push_back(asin((sin_dec_max-sin_dec_min)*ran()+sin_dec_min));
    }

  }

  else
    ErrorMsg("Error in cosmobl::catalogue::Catalogue::Catalogue : the random catalogue has to be cubic!");
  
  
  // extract the redshift distribution  
  
  vector<double> redshift = catalogue.var(_Redshift_);
  vector<double> xx, yy, err;
  distribution(xx, yy, err, redshift, {}, nbin, true, par::defaultString, 1., cosmobl::Min(redshift), cosmobl::Max(redshift), true, conv, sigma);

  
  // extract the random redshifts from the distribution
  
  vector<double> random_redshift;
  fill_distr(ra.size(), xx, yy, random_redshift, catalogue.Min(Var::_Redshift_), catalogue.Max(Var::_Redshift_), 13);

  
  // constructing the random sample
  
  for (size_t i=0; i<random_redshift.size(); i++) {
    observedCoordinates coord = {ra[i], dec[i], random_redshift[i]};
    m_sample.push_back(move(Object::Create(_RandomObject_, coord, cosm)));
  }
  
}


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const RandomType type, const Catalogue catalogue, const double N_R, const int nbin, const double Angle, const vector<double> redshift, const cosmology::Cosmology &cosm, const bool conv, const double sigma, const int seed)
{
  if (type!=_createRandom_cone_) ErrorMsg("Error in cosmobl::catalogue::Catalogue::Catalogue : the random catalogue has to be conic!");
  
  cout << "I'm creating a random catalogue in a cone..." << endl;

  int nRandom = int(catalogue.nObjects()*N_R);
  
  double redshift_min = catalogue.Min(Var::_Redshift_);
  double redshift_max = catalogue.Max(Var::_Redshift_);

  double D_min = cosm.D_C(redshift_min);
  double D_max = cosm.D_C(redshift_max);

  
  vector<double> DD;
 
  default_random_engine gen(seed);
  uniform_real_distribution<float> ran(0., 1.);

  
  // extract the redshift distribution  
  
  vector<double> xx, yy, err;
  distribution(xx, yy, err, redshift, {}, nbin, true, par::defaultString, 1., cosmobl::Min(redshift), cosmobl::Max(redshift), true, conv, sigma);

  
  // extract the random redshifts from the distribution
  
  vector<double> random_redshift;
  fill_distr(nRandom, xx, yy, random_redshift, catalogue.Min(Var::_Redshift_), catalogue.Max(Var::_Redshift_), 13);


  // constructing the random sample
  
  double rMax = cosm.D_C(redshift_max)*tan(Angle);
  
  double Dm = cosmobl::Min(DD)*0.999;
  double DM = cosmobl::Max(DD)*1.001;

  double fact = 1./DM;

  double r1, r2;
  
  for (size_t i=0; i<DD.size(); i++) {

    comovingCoordinates coord;
    
    bool OK = 0;
    while (!OK) {

      coord.xx = fact*DD[i]*rMax*2.*(ran(gen)-0.5);
      coord.yy = DD[i];
      coord.zz = fact*DD[i]*rMax*2.*(ran(gen)-0.5);

      r1 = sqrt(coord.xx*coord.xx+coord.zz*coord.zz);
      r2 = sqrt(coord.xx*coord.xx+coord.zz*coord.zz+DD[i]*DD[i]);
      
      if (r1<fact*DD[i]*rMax && Dm<r2 && r2<DM) OK = 1;
    }

    m_sample.push_back(move(Object::Create(_RandomObject_, coord)));
  }
  
}


// ============================================================================

/// @cond extrandom

cosmobl::catalogue::Catalogue::Catalogue (const RandomType type, const double N_R, const cosmology::Cosmology &cosm, const int step_redshift, const vector<double> lim, const vector<double> redshift, const vector<double> weight, const double redshift_min, const double redshift_max, const bool venice, const string where, string file_random, const string mask, const string dir_venice) 
{
  if (type!=_createRandom_VIPERS_) ErrorMsg("Error in cosmobl::catalogue::Catalogue::Catalogue : the random catalogue has to be of type _VIPERS_ !");
  
  cout << par::col_green << "I'm creating the random catalogue..." << par::col_default << endl;

  int nRandom = int(redshift.size()*N_R);
  
  vector<double> random_redshift, random_DD;


  // compute the redshift distribution and store

  vector<double> xx, fx, err;
  string file_nz = "file_nz";
  distribution(xx, fx, err, redshift, weight, step_redshift, 1, file_nz, 1., cosmobl::Min(redshift)*0.9, cosmobl::Max(redshift)*1.1, 1, 1, 0.05);
  
    
  // construct the random catalogue using venice
  
  if (venice) {
    file_random = "temp";
    string DO = par::DirCosmo+"External/VIPERS/"+dir_venice+"/venice -m "+mask+" -r -f "+where+" -npart "+conv(nRandom, par::fINT)+" -o "+file_random+" -nz "+file_nz;
    cout << endl << "--> " << DO << endl << endl;
    if (system(DO.c_str())) {}
  }
  
  ifstream fin(file_random.c_str()); checkIO(file_random, 1);
  string line;
  getline(fin, line);

  double RA, DEC, REDSHIFT;
  
  while (fin >> RA >> DEC >> REDSHIFT) 
    if (redshift_min<REDSHIFT && REDSHIFT<redshift_max && lim[0]<RA && RA<lim[1] && lim[2]<DEC && DEC<lim[3]) {
      observedCoordinates coord = {radians(RA), radians(DEC), REDSHIFT};
      m_sample.push_back(move(Object::Create(_RandomObject_, coord, cosm)));
    }
  
  computeComovingCoordinates(cosm);
  
  if (system ("rm -f file_nz temp")) {};
}

/// @endcond
