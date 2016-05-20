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

  double XX, YY, ZZ;

  for (int i=0; i<nRandom; i++) {
    XX = ran(gen)*(Xmax-Xmin)+Xmin;
    YY = ran(gen)*(Ymax-Ymin)+Ymin;
    ZZ = ran(gen)*(Zmax-Zmin)+Zmin;
    m_sample.push_back(move(Object::Create(_RandomObject_, XX, YY, ZZ)));
  }
  
}


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const RandomType type, const Cosmology &real_cosm, const Cosmology &test_cosm, const string dir_in, const string dir_out, const double Zguess_min, const double Zguess_max)
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

  for (size_t i=0; i<random_X.size(); i++) 
    m_sample.push_back(move(Object::Create(_RandomObject_, random_X[i], random_Y[i], random_Z[i])));
  
}


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const RandomType type, const Catalogue catalogue, const Cosmology &cosm, const int N_R, const int nbin, const bool conv, const double sigma, const int seed)
{
  vector<double> ra = catalogue.var(_RA_), raB = ra;
  vector<double> dec = catalogue.var(_Dec_), decB = dec;

  if (type==_createRandom_shuffle_) {

    cout << "I'm creating a random catalogue with the 'shuffle' method..." << endl;
    
    // the R.A.-Dec coordinates of the random catalogue will be the same as the ones of the real catalogue
    for (int i=0; i<max(0, N_R-1); i++) { // use N_R more random objects
      ra.insert(ra.end(), raB.begin(), raB.end());
      dec.insert(dec.end(), decB.begin(), decB.end());
    }

  }
  else if(type==_createRandom_box_RADec_){
    cout << "I'm creating a random catalogue in a RA-Dec box..." << endl;

    int nRandom = int(catalogue.nObjects()*N_R);

    ra.erase(ra.begin(),ra.end());
    dec.erase(dec.begin(), dec.end());
        
    default_random_engine generator(seed);
    uniform_real_distribution<double> distribution;
    auto ran = bind(distribution,generator);

    double ra_min = Min(Var::_RA_); double ra_max = Max(Var::_RA_);
    double sin_dec_min = sin(Min(Var::_Dec_)); double sin_dec_max = sin(Max(Var::_Dec_));

    for(int i=0;i<nRandom;i++){
      ra.push_back((ra_max-ra_min)*ran()+ra_min);
      dec.push_back(asin((sin_dec_max-sin_dec_min)*ran()+sin_dec_min));
    }

  }
  else
    ErrorMsg("Error in cosmobl::catalogue::Catalogue::Catalogue : the random catalogue has to be cubic!");

  
  // smoothing the redshift distribution by Gaussian kernel with 0 mean and sigma standard deviation 
  
  vector<double> redshift = catalogue.var(_Redshift_), redshiftB = redshift;
  vector<double> xx, yy, err;
  distribution(xx, yy, err, redshift, {}, nbin, 1, par::defaultString, 1, cosmobl::Min(redshift), cosmobl::Max(redshift), 1, conv, sigma);

  
  // extracting the redshifts from the smoothed distribution
  
  vector<double> random_redshift;
  random_numbers (ra.size(), 13, xx, yy, random_redshift);
  
  
  // constructing the random sample
  
  for (size_t i=0; i<random_redshift.size(); i++) 
    m_sample.push_back(move(Object::Create(_RandomObject_, ra[i], dec[i], random_redshift[i], cosm)));
}


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const RandomType type, const Catalogue catalogue, const int nRandom, const Cosmology &cosm, const double Angle, const int step_redshift, const vector<double> redshift, vector<double> &dc, vector<double> &convol, const int idum)
{
  if (type!=_createRandom_cone_) ErrorMsg("Error in cosmobl::catalogue::Catalogue::Catalogue : the random catalogue has to be conic!");
  
  cout << "I'm creating a random catalogue in a cone..." << endl;
  
  double redshift_min = catalogue.Min(Var::_Redshift_);
  double redshift_max = catalogue.Max(Var::_Redshift_);

  double D_min = cosm.D_C(redshift_min);
  double D_max = cosm.D_C(redshift_max);

  
  vector<double> DD;
 
  default_random_engine gen;
  uniform_real_distribution<float> ran(0., 1.);

  
  if (step_redshift>0) {
    
    // distribution along the line of sight
    
    redshift_min *= 0.9999;
    redshift_max *= 1.0001;
    double delta_redshift = (redshift_max-redshift_min)/step_redshift;
    double redshift1 = redshift_min;
    double redshift2 = redshift1+delta_redshift;
    vector<int> nObj(step_redshift,0);
    for (int i=0; i<step_redshift; i++) {
      for (size_t oo=0; oo<redshift.size(); oo++) 
	if (redshift1<redshift[oo] && redshift[oo]<redshift2) 
	  nObj[i] ++;
      redshift1 = redshift2;
      redshift2 += delta_redshift;
    }
  
    int N_R = int(double(nRandom)/double(redshift.size()));
    
    int nTOT = 0;
    for (int i=0; i<step_redshift; i++) {
      nObj[i] *= N_R;
      nTOT += nObj[i];
    }
    nObj[step_redshift-1] -= (nTOT-nRandom);

    redshift1 = redshift_min;
    redshift2 = redshift1+delta_redshift;
    
    for (int i=0; i<step_redshift; i++) {
      for (int j=0; j<nObj[i]; j++) {
	double REDSHIFT = redshift1+(redshift2-redshift1)*ran(gen);
	DD.push_back(cosm.D_C(REDSHIFT));
      }
      redshift1 = redshift2;
      redshift2 += delta_redshift;
    }
  
  }

  // from the convolvolution of N(D_C)

  else if (step_redshift==0) 
    random_numbers(nRandom, idum, dc, convol, DD, D_min, D_max);

  else ErrorMsg("Error in random_catalogue_cone of RandomCatalogue.cpp!");


  double XX, ZZ;

  double rMax = cosm.D_C(redshift_max)*tan(Angle);
  
  double Dm = cosmobl::Min(DD)*0.999;
  double DM = cosmobl::Max(DD)*1.001;

  double fact = 1./DM;

  double r1, r2;
  
  for (size_t i=0; i<DD.size(); i++) {

    bool OK = 0;
    while (!OK) {

      XX = fact*DD[i]*rMax*2.*(ran(gen)-0.5);
      ZZ = fact*DD[i]*rMax*2.*(ran(gen)-0.5);

      r1 = sqrt(XX*XX+ZZ*ZZ);
      r2 = sqrt(XX*XX+ZZ*ZZ+DD[i]*DD[i]);
      
      if (r1<fact*DD[i]*rMax && Dm<r2 && r2<DM) OK = 1;
    }

    m_sample.push_back(move(Object::Create(_RandomObject_, XX, DD[i], ZZ)));
  }
  
}


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const RandomType type, const Catalogue catalogue, const int nRandom, const Cosmology &cosm, const string dir, const int step_redshift, const vector<double> redshift, vector<double> &dc, vector<double> &convol, const int idum)
{
  if (type!=_createRandom_mock_) ErrorMsg("Error in cosmobl::catalogue::Catalogue::Catalogue : the random catalogue has to be of type _RandomMock_ !");
  
  cout << "I'm creating a random catalogue in a cone (for a simulated catalogue)..." << endl;
  
  double ra_min = catalogue.Min(Var::_RA_), ra_max = catalogue.Max(Var::_RA_);
  double dec_min = catalogue.Min(Var::_Dec_), dec_max = catalogue.Max(Var::_Dec_);
  double sin_dec_min = sin(dec_min), sin_dec_max=sin(dec_max);
  double redshift_min = catalogue.Min(Var::_Redshift_), redshift_max = catalogue.Max(Var::_Redshift_);

  cout << ra_min << " < ra < " << ra_max << "   " << dec_min << " < dec < " << dec_max << "   " << redshift_min << " < z < " << redshift_max << endl;

  vector<double> DD;
 
  default_random_engine gen;
  uniform_real_distribution<float> ran(0., 1.);

  if (step_redshift>0) {
    
    // distribution along the line of sight

    redshift_min *= 0.9999;
    redshift_max *= 1.0001;
    double delta_redshift = (redshift_max-redshift_min)/step_redshift;
    double redshift1 = redshift_min;
    double redshift2 = redshift1+delta_redshift;
    vector<int> nObj(step_redshift, 0);
    for (int i=0; i<step_redshift; i++) {
      for (size_t oo=0; oo<redshift.size(); oo++) 
	if (redshift1<redshift[oo] && redshift[oo]<redshift2) 
	  nObj[i] ++;
      redshift1 = redshift2;
      redshift2 += delta_redshift;
    }

    int N_R = int(double(nRandom)/double(redshift.size()));
    int nTOT = 0;
    for (int i=0; i<step_redshift; i++) {
      nObj[i] *= N_R;
      nTOT += nObj[i];
    }
    nObj[step_redshift-1] -= (nTOT-nRandom);

    redshift1 = redshift_min;
    redshift2 = redshift1+delta_redshift;
    for (int i=0; i<step_redshift; i++) {
      for (int j=0; j<nObj[i]; j++) {
	double REDSHIFT = redshift1+(redshift2-redshift1)*ran(gen);
	DD.push_back(cosm.D_C(REDSHIFT));
      }
      redshift1 = redshift2;
      redshift2 += delta_redshift;
    }

  }

  // from the convolvolution of N(D_C)

  else if (step_redshift==0) {
    double D_min = cosm.D_C(redshift_min);
    double D_max = cosm.D_C(redshift_max);
    random_numbers(nObjects(), idum, dc, convol, DD, D_min, D_max);
  }

  else ErrorMsg("Error in random_catalogue_mock of RandomCatalogue.cpp!");
  if (DD.size()==0) ErrorMsg("Error in random_catalogue_mock of RandomCatalogue.cpp: DD.size()=0!");

  double XX, YY, ZZ, RA, DEC;

  for (size_t i=0; i<DD.size(); i++) {
    RA = ra_min+(ra_max-ra_min)*ran(gen);
    DEC = asin(sin_dec_min+(sin_dec_max-sin_dec_min)*ran(gen));

    XX = DD[i]*cos(DEC)*sin(RA);
    YY = DD[i]*cos(DEC)*cos(RA);
    ZZ = DD[i]*sin(DEC);

    m_sample.push_back(move(Object::Create(_RandomObject_, XX, YY, ZZ)));
  }
  
}


// ============================================================================

/// @cond extrandom

cosmobl::catalogue::Catalogue::Catalogue (const RandomType type, const int nRandom, const Cosmology &cosm, const string dir_out, const int step_redshift, const vector<double> dc, const vector<double> convol, const vector<double> lim, const vector<double> redshift, const bool venice, string file_random, const string mask, const string dir_venice, const int idum) 
{
  if (type!=_createRandom_VIPERS_) ErrorMsg("Error in cosmobl::catalogue::Catalogue::Catalogue : the random catalogue has to be of type _VIPERS_ !");
  
  cout << par::col_green << "I'm creating the random catalogue..." << par::col_default << endl;

  vector<double> random_redshift, random_DD;


  // compute the redshift distribution and store

  vector<double> xx, fx, err;
  double redshift_min = cosmobl::Min(redshift)*0.9999;
  double redshift_max = cosmobl::Max(redshift)*1.0001;
  string file_nz = "file_nz";
  distribution(xx, fx, err, redshift, {}, step_redshift, 1, file_nz, 1., redshift_min, redshift_max, 1, 1, 0.05);
  
    
  // construct the random catalogue using venice
  
  if (venice) {
    file_random = "temp";
    string DO = par::DirCosmo+"External/VIPERS/"+dir_venice+"/venice -m "+mask+" -r -f inside -npart "+conv(nRandom, par::fINT)+" -o "+file_random+" -nz "+file_nz;
    cout << endl << "--> " << DO << endl << endl;
    if (system(DO.c_str())) {}
  }

  ifstream fin(file_random.c_str()); checkIO(file_random, 1);
  string line;
  getline(fin, line);

  double RA, DEC, REDSHIFT;
  
  while (fin >> RA >> DEC >> REDSHIFT) 
    if (redshift_min<REDSHIFT && REDSHIFT<redshift_max && lim[0]<RA && RA<lim[1] && lim[2]<DEC && DEC<lim[3]) 
      m_sample.push_back(move(Object::Create(_RandomObject_, radians(RA), radians(DEC), REDSHIFT, cosm)));

  computeComovingCoordinates(cosm);
  
  if (system ("rm -f file_nz temp")) {};
  
}

/// @endcond
