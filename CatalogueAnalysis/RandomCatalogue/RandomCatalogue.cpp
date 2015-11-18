/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli                         *
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
 *  @file RandomCatalogue/RandomCatalogue.cpp
 *
 *  @brief Functions for random catalogues
 *
 *  This file contains the implementation of a set of functions to
 *  create random catalogues
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "RandomCatalogue.h"
using namespace cosmobl;


// ============================================================================


shared_ptr<Catalogue> cosmobl::random_catalogue_fromFile (const vector<string> file, const double nSub) 
{
  cout <<"I'm reading the random catalogue..."<<endl;

  shared_ptr<Catalogue> random_catalogue (new Catalogue);
  
  vector<double> random_X, random_Y,  random_Z;

  default_random_engine gen;
  uniform_real_distribution<float> ran(0., 1.);

  for (size_t dd=0; dd<file.size(); dd++) {
    string file_in = file[dd];
    cout <<"file_in = "<<file_in<<endl;
    ifstream finr (file_in.c_str()); checkIO (file_in,1);
    
    string line; double NUM;
    while (getline(finr, line)) {
      if (ran(gen)<nSub) {
	stringstream ss(line);
	ss>>NUM; random_X.push_back(NUM);
	ss>>NUM; random_Y.push_back(NUM);
	ss>>NUM; random_Z.push_back(NUM);
      }
    }
    finr.clear(); finr.close(); 
  }


  vector<shared_ptr<Object>> data_random;

  for (size_t i=0; i<random_X.size(); i++) {
    shared_ptr<RandomObject> SMP(new RandomObject(random_X[i], random_Y[i], random_Z[i]));
    data_random.push_back(move(SMP));
  }
  
  random_catalogue->add_objects(data_random);

  return random_catalogue;
}


// ============================================================================


shared_ptr<Catalogue> cosmobl::random_catalogue_radecred_fromFile (const string file_in, const double z_min, const double z_max, const Cosmology &cosm, const double nSub, const double fact) 
{      
  cout <<"I'm reading the random catalogue: "<<file_in<<"..."<<endl;
  
  shared_ptr<Catalogue> random_catalogue (new Catalogue);
  
  vector<double> random_X, random_Y, random_Z, random_ra, random_dec, random_redshift, random_dc;

  ifstream finr (file_in.c_str());  checkIO (file_in,1);
 
  string line;
  double NUM, RA, DEC;

  default_random_engine gen;
  uniform_real_distribution<float> ran(0., 1.);

  while (getline(finr,line)) {
    stringstream ss(line);

    NUM = -1; ss>>NUM; 
    if (NUM>0) {
      RA = NUM*fact;
      ss>>NUM; DEC = NUM*fact;
      ss>>NUM; 
    
      if (z_min*0.99999<NUM && NUM<z_max*1.00001 && ran(gen)<nSub) {
	random_ra.push_back(RA);
	random_dec.push_back(DEC);
	random_redshift.push_back(NUM);
      }
    }
  }
  finr.clear(); finr.close(); 

  
  vector<shared_ptr<Object>> data_random;

  for (size_t i=0; i<random_ra.size(); i++) {
    shared_ptr<RandomObject> SMP(new RandomObject(random_ra[i], random_dec[i], random_redshift[i], cosm));
    data_random.push_back(move(SMP));
  }
  
  random_catalogue->add_objects(data_random);

  return random_catalogue;
}


// ============================================================================


shared_ptr<Catalogue> cosmobl::random_catalogue_box (const shared_ptr<Catalogue> catalogue, const int nRandom)
{
  // ------ create the random catalogue ------ 
  
  cout <<"I'm creating a random catalogue..."<<endl;

  shared_ptr<Catalogue> random_catalogue (new Catalogue);

  vector<double> random_X, random_Y, random_Z;

  vector<double> Lim; 
  catalogue->MinMax_var(Var::_XX_, Lim, 0);
  catalogue->MinMax_var(Var::_YY_, Lim, 0);
  catalogue->MinMax_var(Var::_ZZ_, Lim, 0);
  
  double Xmin = Lim[0]; double Xmax = Lim[1];
  double Ymin = Lim[2]; double Ymax = Lim[3];
  double Zmin = Lim[4]; double Zmax = Lim[5];

  if (Xmin>Xmax || Ymin>Ymax || Zmin>Zmax) ErrorMsg("Error in random_catalogue_box of RandomCatalogue.cpp!");

  default_random_engine gen;
  uniform_real_distribution<float> ran(0., 1.);

  for (int i=0; i<nRandom; i++) {
    random_X.push_back(ran(gen)*(Xmax-Xmin)+Xmin);
    random_Y.push_back(ran(gen)*(Ymax-Ymin)+Ymin);
    random_Z.push_back(ran(gen)*(Zmax-Zmin)+Zmin);
  }

  vector<shared_ptr<Object>> data_random;

  for (size_t i=0; i<random_X.size(); i++) {
    shared_ptr<RandomObject> SMP(new RandomObject(random_X[i], random_Y[i], random_Z[i]));
    data_random.push_back(move(SMP));
  }

  random_catalogue->add_objects(data_random);

  return random_catalogue;
}


// ============================================================================


shared_ptr<Catalogue> cosmobl::warped_random_catalogue (const Cosmology &real_cosm, const Cosmology &test_cosm, const string dir_in, const string dir_out, const double Zguess_min, const double Zguess_max)
{  
  cout <<"I'm creating a random catalogue..."<<endl;

  shared_ptr<Catalogue> random_catalogue (new Catalogue);
  
  vector<double> random_X, random_Y, random_Z;


  // ------ read the random data in the real (or fiducial) cosmology ------   
  
  string file_in = dir_in+"random_catalogue";
  cout <<"file_in = "<<file_in<<endl;
  ifstream fin (file_in.c_str()); checkIO (file_in,1);
  
  string line; double NUM;
  
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

  vector<shared_ptr<Object>> data_random;

  for (size_t i=0; i<random_X.size(); i++) {
    shared_ptr<RandomObject> SMP(new RandomObject(random_X[i], random_Y[i], random_Z[i]));
    data_random.push_back(move(SMP));
  }

  random_catalogue->add_objects(data_random);
  
  return random_catalogue;
}


// ============================================================================


shared_ptr<Catalogue> cosmobl::random_catalogue_cone (const shared_ptr<Catalogue> catalogue, const int nRandom, const Cosmology &cosm, const double Angle, const int step_redshift, const vector<double> redshift, vector<double> &dc, vector<double> &convol, const int idum)
{
  cout <<"I'm creating a random catalogue..."<<endl;
  
  shared_ptr<Catalogue> random_catalogue (new Catalogue);
  
  double redshift_min = catalogue->MinMax_var(Var::_REDSHIFT_)[0];
  double redshift_max = catalogue->MinMax_var(Var::_REDSHIFT_)[1];

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
  
  double Dm = Min(DD)*0.999;
  double DM = Max(DD)*1.001;

  double fact = 1./DM;

  double r1, r2;
  vector<double> random_X, random_Y, random_Z;
  
  for (size_t i=0; i<DD.size(); i++) {

    bool OK = 0;
    while (!OK) {

      XX = fact*DD[i]*rMax*2.*(ran(gen)-0.5);
      ZZ = fact*DD[i]*rMax*2.*(ran(gen)-0.5);

      r1 = sqrt(XX*XX+ZZ*ZZ);
      r2 = sqrt(XX*XX+ZZ*ZZ+DD[i]*DD[i]);

      if (r1<fact*DD[i]*rMax && Dm<r2 && r2<DM) OK = 1;
    }

    random_X.push_back(XX);
    random_Y.push_back(DD[i]);
    random_Z.push_back(ZZ);
  }


  vector<shared_ptr<Object>> data_random;

  for (size_t i=0; i<random_X.size(); i++) {
    shared_ptr<RandomObject> SMP(new RandomObject(random_X[i], random_Y[i], random_Z[i]));
    data_random.push_back(move(SMP));
  }

  random_catalogue->add_objects(data_random);

  return random_catalogue;
}


// ============================================================================


shared_ptr<Catalogue> cosmobl::random_catalogue_mock (shared_ptr<Catalogue> catalogue, int &nRandom, Cosmology &cosm, string &dir, int &step_redshift, vector<double> redshift, vector<double> &dc, vector<double> &convol, int idum)
{
  cout <<"I'm creating a random catalogue..."<<endl;
  
  shared_ptr<Catalogue> random_catalogue (new Catalogue);

  vector<double> random_X, random_Y, random_Z;

  vector<double> Lim; 
  catalogue->MinMax_var(Var::_RA_, Lim, 0);
  catalogue->MinMax_var(Var::_DEC_, Lim, 0);
  catalogue->MinMax_var(Var::_REDSHIFT_, Lim, 0);
  
  double ra_min = Lim[0]; double ra_max = Lim[1];
  double dec_min = Lim[2]; double dec_max = Lim[3];
  double redshift_min = Lim[4]; double redshift_max = Lim[5];

  cout <<ra_min<<" < ra < "<<ra_max<<"   "<<dec_min<<" < dec < "<<dec_max<<"   "<<redshift_min<<" < z < "<<redshift_max<<endl;

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
    random_numbers (random_catalogue->nObjects(), idum, dc, convol, DD, D_min, D_max);
  }

  else ErrorMsg("Error in random_catalogue_mock of RandomCatalogue.cpp!");
  if (DD.size()==0) ErrorMsg("Error in random_catalogue_mock of RandomCatalogue.cpp: DD.size()=0!");

  double XX, YY, ZZ, RA, DEC;
  vector<double> random_ra, random_dec, random_redshift;

  for (size_t i=0; i<DD.size(); i++) {
    RA = ra_min+(ra_max-ra_min)*ran(gen);
    DEC = dec_min+(dec_max-dec_min)*ran(gen);

    XX = DD[i]*cos(DEC)*sin(RA);
    YY = DD[i]*cos(DEC)*cos(RA);
    ZZ = DD[i]*sin(DEC);

    random_X.push_back(XX);
    random_Y.push_back(YY);
    random_Z.push_back(ZZ);
    random_ra.push_back(RA);
    random_dec.push_back(DEC);
    random_redshift.push_back(cosm.Redshift(DD[i], redshift_min, redshift_max));
  }

  vector<shared_ptr<Object>> data_random;

  for (size_t i=0; i<random_X.size(); i++) {
    shared_ptr<RandomObject> SMP(new RandomObject(random_X[i], random_Y[i], random_Z[i]));
    data_random.push_back(move(SMP));
  }

  random_catalogue->add_objects(data_random);
  
  
  string file_out = dir+"random_catalogue";
  ofstream fout (file_out.c_str()); checkIO (file_out,0);
  for (int i=0; i<int(random_X.size()); i++) 
    fout <<random_X[i]<<"   "<<random_Y[i]<<"   "<<random_Z[i]<<"   "<<random_ra[i]<<"   "<<random_dec[i]<<"   "<<random_redshift[i]<<endl;
  fout.clear(); fout.close(); cout <<"I wrote the file "<<file_out<<endl;
  cout <<"...done (with "<<random_X.size()<<" objects)"<<endl;

  return random_catalogue;
}


// ============================================================================

/// @cond extrandom

shared_ptr<Catalogue> cosmobl::random_catalogue_mock_cone (const shared_ptr<Catalogue> catalogue, const int nRandom, const Cosmology &cosm, const string dir, const int step_redshift, vector<double> redshift)
{
  cout <<"I'm creating a random catalogue..."<<endl;

  shared_ptr<Catalogue> random_catalogue (new Catalogue);
  
  vector<double> random_X, random_Y, random_Z;
  
  vector<double> Lim; 
  catalogue->MinMax_var(Var::_RA_, Lim, 0);
  catalogue->MinMax_var(Var::_DEC_, Lim, 0);
  catalogue->MinMax_var(Var::_REDSHIFT_, Lim, 0);

  double ra_min = Lim[0]; double ra_max = Lim[1];
  double dec_min = Lim[2]; double dec_max = Lim[3];
  double redshift_min = Lim[4]; double redshift_max = Lim[5];

  cout <<ra_min<<" < ra < "<<ra_max<<"   "<<dec_min<<" < dec < "<<dec_max<<"   "<<redshift_min<<" < z < "<<redshift_max<<endl;
   
  double XX, YY, ZZ, RA, DEC, DD, REDSHIFT;
  vector<double> random_ra, random_dec, random_redshift;
  
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
    
	REDSHIFT = redshift1+(redshift2-redshift1)*ran(gen);
	DD = cosm.D_C(REDSHIFT);

	bool DO = 1;
	while (DO) {

	  RA = ra_min+(ra_max-ra_min)*ran(gen);
	  DEC = dec_min+(dec_max-dec_min)*ran(gen);
	  //RAD = sqrt(pow(RA,2)+pow(DEC,2));
	
	  //if (RAD<ra_max && RAD<dec_max) {
	  XX = DD*cos(DEC)*sin(RA);
	  YY = DD*cos(DEC)*cos(RA);
	  ZZ = DD*sin(DEC);
	    
	  random_X.push_back(XX);
	  random_Y.push_back(YY);
	  random_Z.push_back(ZZ);
	  random_ra.push_back(RA);
	  random_dec.push_back(DEC);
	  random_redshift.push_back(REDSHIFT);
	  DO = 0;
	  //}
	}
      }
   
      redshift1 = redshift2;
      redshift2 += delta_redshift;
    } 
  }

  else { //shuffling the redshifts of the catalogue
  
    if (int(redshift.size())!=nRandom) ErrorMsg("Error in random_catalogue_mock_Euclid of RandomCatalogue.cpp!");

    random_shuffle (redshift.begin(), redshift.end());
    
    int index = 0;

    while (int(random_X.size())<nRandom) {

      REDSHIFT = redshift[index];
      DD = cosm.D_C(REDSHIFT);

      RA = ra_min+(ra_max-ra_min)*ran(gen);
      DEC = dec_min+(dec_max-dec_min)*ran(gen);
      //RAD = sqrt(pow(RA,2)+pow(DEC,2));
  
      //if (RAD<ra_max && RAD<dec_max) {
      XX = DD*cos(DEC)*sin(RA);
      YY = DD*cos(DEC)*cos(RA);
      ZZ = DD*sin(DEC);
	
      random_X.push_back(XX);
      random_Y.push_back(YY);
      random_Z.push_back(ZZ);
      random_ra.push_back(RA);
      random_dec.push_back(DEC);
      random_redshift.push_back(REDSHIFT);
      index ++;
      //}
    }
  }


  string file_out = dir+"random_catalogue";
  ofstream fout (file_out.c_str()); checkIO (file_out,0);
  for (int i=0; i<int(random_X.size()); i++) 
    fout <<random_X[i]<<"   "<<random_Y[i]<<"   "<<random_Z[i]<<"   "<<random_ra[i]<<"   "<<random_dec[i]<<"   "<<random_redshift[i]<<endl;
  fout.clear(); fout.close(); cout <<"I wrote the file "<<file_out<<endl;
  cout <<"...done (with "<<random_X.size()<<" objects)"<<endl;
  
  random_catalogue->set_var(Var::_XX_, random_X); random_catalogue->set_var(Var::_YY_, random_Y); random_catalogue->set_var(Var::_ZZ_, random_Z);

  return random_catalogue;
}


// ============================================================================


shared_ptr<Catalogue> cosmobl::random_catalogue_VIPERS (const int nRandom, const Cosmology &cosm, const string dir_out, const int step_redshift, const vector<double> dc, const vector<double> convol, const vector<double> lim, const vector<double> redshift, const bool venice, string file_random, const string mask, const string dir_venice, const int idum) 
{
  cout <<"I'm creating a random catalogue..."<<endl;

  shared_ptr<Catalogue> random_catalogue (new Catalogue);
  
  vector<double> random_X, random_Y, random_Z, random_ra, random_dec, random_redshift, random_DD;

  
  // construct the random catalogue using venice
  
  if (venice) {
    file_random = "temp";
    string DO = par::DirCosmo+"CatalogueAnalysis/RandomCatalogue/VIPERS/"+dir_venice+"/venice -m "+mask+" -r -f inside -npart "+conv(nRandom,par::fINT)+" -o "+file_random;
    if (system (DO.c_str())) {};
  }


  // count the number of random galaxies 
  
  cout <<"I'm reading the random catalogue: "<<file_random<<"..."<<endl;


  ifstream fin(file_random.c_str()); checkIO(file_random,1);
  string line;
  int ind = -1; // (the first line is not counted)
  while (getline(fin, line)) ind ++;
  fin.clear();
  int nRan = ind;

  
  default_random_engine gen;
  uniform_real_distribution<float> ran(0., 1.);

  // CHECK!!!
  double redshift_min = Min(redshift)*0.9999;
  double redshift_max = Max(redshift)*1.0001;


  if (step_redshift>0) {
    // distribution along the line of sight
    
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
    
    int N_R = int(double(nRan)/double(redshift.size()));
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
	random_redshift.push_back(REDSHIFT);
	random_DD.push_back(cosm.D_C(REDSHIFT));
      }
      redshift1 = redshift2;
      redshift2 += delta_redshift;
    }
  }


  // from the convolvolution of N(D_C) (with the Vmax method or with a Gaussian filter)

  else if (step_redshift==0 || step_redshift==-1) {
    double D_min = cosm.D_C(redshift_min);
    double D_max = cosm.D_C(redshift_max);
    random_numbers (nRan, idum, dc, convol, random_DD, D_min, D_max);
   
    double z_min = 0.3, z_max = 3.; // check!!!
    for (size_t i=0; i<random_DD.size(); i++)
      random_redshift.push_back(cosm.Redshift(random_DD[i], z_min, z_max));
  }

  else if (step_redshift==-2) {
    for (int i=0; i<nRan; i++) {
      int ind = ran(gen)*(redshift.size()-1);
      double REDSHIFT = redshift[ind];
      random_redshift.push_back(REDSHIFT);
      random_DD.push_back(cosm.D_C(REDSHIFT));
    }
  }

  else ErrorMsg("Error in cosmobl::random_catalogue_VIPERS of RandomCatalogue.cpp!");


  double RA, DEC, XX, YY, ZZ;

  fin.seekg(ios::beg);
  getline(fin, line);

  for (size_t i=0; i<random_DD.size(); i++) {
    fin >>RA>>DEC;
    if (lim[0]<RA && RA<lim[1] && lim[2]<DEC && DEC<lim[3]) {
      //RA *= par::pi/180.;
      //DEC *= par::pi/180.;

      XX = random_DD[i]*cos(DEC)*sin(RA);
      YY = random_DD[i]*cos(DEC)*cos(RA);
      ZZ = random_DD[i]*sin(DEC);

      //if (ran(gen)>0.8) {
      random_X.push_back(XX);
      random_Y.push_back(YY);
      random_Z.push_back(ZZ);
      random_ra.push_back(RA);
      random_dec.push_back(DEC);
      //}

    }
  }
  
  nRan = int(random_X.size());

  
  vector<shared_ptr<Object>> data_random;

  for (size_t i=0; i<random_X.size(); i++) {
    shared_ptr<RandomObject> SMP(new RandomObject(random_X[i], random_Y[i], random_Z[i]));
    data_random.push_back(move(SMP));
  }
  
  random_catalogue->add_objects(data_random);

  if (system ("rm -f temp")) {};
  
  return random_catalogue;
}


// ============================================================================


/* ======== Alfonso Veropalumbo ======== */

void cosmobl::random_redshift_distribution (const shared_ptr<Catalogue> random, const shared_ptr<Catalogue> data, const string dir_random, const int nbin, const bool convolution, const double sigma)
{
  if (random->nObjects()==0) ErrorMsg("Error in random_redshift_distribution of RandomCatalogue.cpp: random nObject=0!");

   int nRan = random->nObjects();

   // distribution vectors

   vector<double> var, distr;
   vector<double> random_redshift;
   vector<double> zminmax;

   data->MinMax_var(Var::_REDSHIFT_, zminmax);
   data->var_distr(Var::_REDSHIFT_, var, distr, nbin, 1, "NULL", -1.e30, -1.e30, convolution, sigma);

   fill_distr(nRan, var, distr, random_redshift, zminmax[0], zminmax[1], 112312);
   
   random->set_var(Var::_REDSHIFT_, random_redshift);

   vector<double> random_ra = random->var(Var::_RA_), random_dec = random->var(Var::_DEC_);
   
   string conv_string = (convolution) ? "_convolution_nBin="+conv(nbin,par::fINT)+"_sigma="+conv(sigma,par::fDP2) : "_distribution_nBin="+conv(nbin,par::fINT);
   string out_random = dir_random+"random_catalogue"+conv_string;
   ofstream fout (out_random.c_str()); checkIO (out_random,0);

   for (int i=0; i<nRan; i++)
     fout << random_ra[i] << " " << random_dec[i] << " " << random_redshift[i] << endl;
   
   fout.clear(); fout.close();
}


// ============================================================================


shared_ptr<Catalogue> cosmobl::random_sdss_angular_distribution (const int nRandom, const string out_dir, const string chunk_list, const bool veto)
{
  shared_ptr<Catalogue> random_catalogue (new Catalogue);

  double d2r = par::pi/180.;

  ifstream fin (chunk_list.c_str()); checkIO (chunk_list,1);

  vector<string> chunk;
  string data;

  while (fin >> data) 
    chunk.push_back(data);
   
  fin.clear(); fin.close();

  unique_unsorted(chunk);
   

  // MANGLE directories

  string mangle_dir = par::DirCosmo+"/CatalogueAnalysis/RandomCatalogue/mangle/";
  string mangle_bin = mangle_dir+"bin/";
  string mangle_plate = mangle_dir+"mask/util/plate.list";
  string mangle_reject = mangle_dir+"mask/reject/";

  string cmd = "mkdir -p "+out_dir+"/temp/"; if (system(cmd.c_str())) {};


  // read plates
   
  string out_plate_list = out_dir+"temp/plate.circles";
  ofstream fout (out_plate_list.c_str()); checkIO (out_plate_list,0);

  fin.open (mangle_plate.c_str()); checkIO (mangle_plate,1);
  getline (fin, data); // skip the header

  while (getline(fin, data)) {
    stringstream ss(data); string ch; double RA, DEC; 
    ss >> RA; ss >> DEC; ss >> ch;
     
    for (size_t i=0; i<chunk.size(); i++)
      if(ch==chunk[i])
	fout << RA << "  " << DEC << "  " << 1.49 << endl; // 1.49 degree is the radius of BOSS
  }

  fout.clear(); fout.close(); 
  fin.clear(); fin.close();


  // pixelize

  string out_pix = out_dir+"temp/pix.out";
  string pixelize = mangle_bin+"pixelize -icd -op "+out_plate_list+" "+out_pix;
  if (system(pixelize.c_str())) {};


  // snap

  string out_snap = out_dir+"temp/snap.out";
  string snap = mangle_bin+"snap "+out_pix+" "+out_snap;
  if (system(snap.c_str())) {};


  // balkanize

  string out_balk = out_dir+"survey_footprint.pix";
  string balkanize = mangle_bin+"balkanize "+out_snap+" "+out_balk;
  if (system(balkanize.c_str())) {};


  // ransack

  string random_file = out_dir+"random_catalogue->angular";
   
  string ransack = mangle_bin+"ransack -r"+conv(nRandom,par::fINT)+" "+out_balk+" "+random_file;
  if (system(ransack.c_str())) {};

  cmd = "sed -i 1d "+random_file;
  if (system(cmd.c_str())) {};

  vector<double> random_ra, random_dec;
   
  if (veto) {
    vector<string> veto_mask(4);
    veto_mask[0] = mangle_reject+"badfield_mask_unphot-ugriz_pix.ply";
    veto_mask[1] = mangle_reject+"centerpost_mask.ply";
    veto_mask[2] = mangle_reject+"bright_star_mask_pix.ply";
    veto_mask[3] = mangle_reject+"collision_priority_mask.ply";
     
    for (size_t i=0; i<veto_mask.size(); i++) {
      string random_temp = out_dir+"random_temp";
      string polyid = mangle_bin+"polyid "+veto_mask[i]+" "+random_file+" "+out_dir+"random_temp";
      if (system(polyid.c_str())) {};

      cmd ="sed -i 1d "+random_temp; // delete a useless header
      if (system(cmd.c_str())) {};

      fin.open (random_temp.c_str());
      fout.open (random_file.c_str());
       
      while (getline(fin,data)) {
	double RA,DEC,PP=-1; stringstream ss(data);
	ss >> RA; ss >> DEC; ss >> PP;
	if (PP==-1)
	  fout << RA << " " << DEC << endl;	 
      }
      fin.clear(); fin.close();
      fout.clear(); fout.close();
    }
  }

  fin.open (random_file.c_str()); checkIO (random_file,1);
  while (getline(fin,data)) {
    stringstream ss(data); double NUM;
    ss >> NUM; random_ra.push_back(NUM*d2r);
    ss >> NUM; random_dec.push_back(NUM*d2r);
  }
  fin.clear(); fin.close();

  random_catalogue->resize(random_ra.size());
 
  random_catalogue->set_var(Var::_RA_,random_ra);
  random_catalogue->set_var(Var::_DEC_,random_dec);
   
  string RM = "rm -rf "+out_dir+"temp";
  if (system(RM.c_str())) {};

  return random_catalogue;
}

/// @endcond
