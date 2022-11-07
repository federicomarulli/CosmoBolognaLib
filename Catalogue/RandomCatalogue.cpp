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
 *  Catalogue used to create random catalogues
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unibo.it
 */

#include "Catalogue.h"
#include "Histogram.h"

using namespace std;

using namespace cbl;


// ============================================================================


cbl::catalogue::Catalogue::Catalogue (const RandomType type, const cosmology::Cosmology &real_cosm, const cosmology::Cosmology &test_cosm, const std::string dir_in, const double Zguess_min, const double Zguess_max)
{
  if (type!=RandomType::_createRandom_box_) ErrorCBL("the random catalogue has to be cubic!", "Catalogue", "RandomCatalogue.cpp");

  coutCBL << "I'm creating a random catalogue with warped cubic geometry, due to geometric distortions: the undistorted random catalogue is read from a file..." << endl;

  string file_in = dir_in+"random_catalogue";
  coutCBL << "file_in = " << file_in << endl;
  ifstream fin(file_in.c_str()); checkIO(fin, file_in);
  
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
    m_object.push_back(move(Object::Create(ObjectType::_Random_, coord)));
  }

}


// ============================================================================

cbl::catalogue::Catalogue::Catalogue (const RandomType type, const Catalogue catalogue, const double N_R, const int nbin, const cosmology::Cosmology &cosm, const bool conv, const double sigma, const std::vector<double> redshift, const std::vector<double> RA, const std::vector<double> Dec, const int z_ndigits, const int seed)
{
  size_t nRandom = int(N_R*catalogue.nObjects());

  vector<double> ra, dec;
  if (RA.size()>0) ra = RA;
  if (Dec.size()>0) dec = Dec;
  if (ra.size()!=dec.size()) ErrorCBL("the dimension of observed random coordinates read in input is wrong!", "Catalogue", "RandomCatalogue.cpp");
  if (ra.size()>0 && ra.size()!=nRandom) nRandom = ra.size();
  
  if (type==RandomType::_createRandom_box_) {
    
    coutCBL << "* * * I'm creating a random catalogue with cubic geometry * * *" << endl;

    double Xmin = catalogue.Min(Var::_X_), Xmax = catalogue.Max(Var::_X_); 
    double Ymin = catalogue.Min(Var::_Y_), Ymax = catalogue.Max(Var::_Y_);
    double Zmin = catalogue.Min(Var::_Z_), Zmax = catalogue.Max(Var::_Z_);

    if (Xmin>Xmax || Ymin>Ymax || Zmin>Zmax) ErrorCBL("wrong values of the coordinates in the construction of the random catalogue; the following conditions have to be satisfied: Xmin<=Xmax, Ymin<=Ymax and Zmin<=Zmax!", "Catalogue", "RandomCatalogue.cpp");

    random::UniformRandomNumbers ran(0., 1., seed);
    
    for (size_t i=0; i<nRandom; ++i) {
      comovingCoordinates coord;
      coord.xx = ran()*(Xmax-Xmin)+Xmin;
      coord.yy = ran()*(Ymax-Ymin)+Ymin;
      coord.zz = ran()*(Zmax-Zmin)+Zmin;
      m_object.push_back(move(Object::Create(ObjectType::_Random_, coord)));
    }

  }

  else if (type==RandomType::_createRandom_shuffleTOT_) {

    if (ra.size()==0) {
      ra = catalogue.var(Var::_RA_);
      dec = catalogue.var(Var::_Dec_);
      vector<double> raB = ra, decB = dec;
      
      // the R.A.-Dec coordinates of the random catalogue will be the same as the ones of the real catalogue
      for (int i=0; i<max(0, int(N_R)-1); i++) { // use N_R times more random objects
        ra.insert(ra.end(), raB.begin(), raB.end());
        dec.insert(dec.end(), decB.begin(), decB.end());
      }
    }
    
    random::UniformRandomNumbers ran(0., catalogue.nObjects()-1, seed);
    
    // constructing the random sample
    for (size_t i=0; i<ra.size(); i++) {
      observedCoordinates coord = {ra[i], dec[i], round_to_precision(catalogue.redshift(ran()), z_ndigits)};
      m_object.push_back(move(Object::Create(ObjectType::_Random_, coord, cosm)));
    }
  }

  else {
    
    if (ra.size()==0) {

      ra = catalogue.var(Var::_RA_);
      dec = catalogue.var(Var::_Dec_);
      vector<double> raB = ra , decB = dec;
      
      if (type==RandomType::_createRandom_shuffle_) {
    
        coutCBL << "I'm creating a random catalogue with the 'shuffle' method..." << endl;
    
        // the R.A.-Dec coordinates of the random catalogue will be the same as the ones of the real catalogue
              
        for (int i=0; i<max(0, int(N_R)-1); i++) { // use N_R times more random objects
          ra.insert(ra.end(), raB.begin(), raB.end());
          dec.insert(dec.end(), decB.begin(), decB.end());
        }
      }
      
      else if (type==RandomType::_createRandom_square_) {
        coutCBL << "I'm creating a random catalogue in a RA-Dec square..." << endl;

        ra.erase(ra.begin(), ra.end());
        dec.erase(dec.begin(), dec.end());
              
        random::UniformRandomNumbers ran(0., 1., seed);

        const double ra_min = catalogue.Min(Var::_RA_),
        ra_max = catalogue.Max(Var::_RA_),
        sin_dec_min = sin(catalogue.Min(Var::_Dec_)),
        sin_dec_max = sin(catalogue.Max(Var::_Dec_));
    
        for (size_t i=0; i<nRandom; i++) {
          ra.push_back((ra_max-ra_min)*ran()+ra_min);
          dec.push_back(asin((sin_dec_max-sin_dec_min)*ran()+sin_dec_min));
    
          /*
            const double cos_ra_min = cos(catalogue.Min(Var::_RA_)),
            cos_ra_max = cos(catalogue.Max(Var::_RA_)), 
            dec_min = catalogue.Min(Var::_Dec_),
            dec_max = catalogue.Max(Var::_Dec_);
                    
            for (size_t i=0; i<nRandom; i++) {
            ra.push_back(acos(cos_ra_min+(cos_ra_max-cos_ra_min)*ran()));
            dec.push_back((dec_max-dec_min)*ran()+dec_min);
          */
    
         }
      }
      else
        ErrorCBL("the chosen random catalogue type is not allowed!", "Catalogue", "RandomCatalogue.cpp");
    }
    
  
    // compute the redshifts of the random objects
  
    vector<double> random_redshift = (redshift.size()>0) ? redshift : catalogue.var(Var::_Redshift_);
    
    if (cbl::Max(random_redshift)-cbl::Min(random_redshift)>1.e-30) {
    
      // extract the redshift distribution of the input catalogue  
      vector<double> xx, yy, err;
      distribution(xx, yy, err, random_redshift, {}, nbin, true, par::defaultString, 1., cbl::Min(random_redshift), cbl::Max(random_redshift), "Linear", conv, sigma);
            
      // extract the random redshifts from the distribution
      random_redshift.resize(ra.size());
      random_redshift = vector_from_distribution(ra.size(), xx, yy, catalogue.Min(Var::_Redshift_), catalogue.Max(Var::_Redshift_), 13);
      
    }
  
    // constructing the random sample
    for (size_t i=0; i<random_redshift.size(); i++) {
      observedCoordinates coord = {ra[i], dec[i], round_to_precision(random_redshift[i], z_ndigits)};
      m_object.push_back(move(Object::Create(ObjectType::_Random_, coord, cosm)));
    }  
  }
}


// ============================================================================

cbl::catalogue::Catalogue::Catalogue (const RandomType type, const int N_R, const double z_step, const Catalogue catalogue, const cosmology::Cosmology &cosm, const std::vector<double> RA, const std::vector<double> Dec, const double sigma_kernel, const int nbins, const int z_ndigits, const int times_default, double times_change, const double tolerance, const std::string out_path_nz, const std::string out_file_nz, const int seed)
{

  size_t nRandom = int(N_R*catalogue.nObjects());

  vector<double> ra, dec;
  if (RA.size()>0) ra = RA;
  if (Dec.size()>0) dec = Dec;
  if (ra.size()!=dec.size()) ErrorCBL("the dimension of observed random coordinates read in input is wrong!", "Catalogue", "RandomCatalogue.cpp");
  if (ra.size()>0 && ra.size()!=nRandom) nRandom = ra.size();

  if (ra.size()==0) {
    ra = catalogue.var(Var::_RA_);
    dec = catalogue.var(Var::_Dec_);
    vector<double> raB = ra, decB = dec;
      
    // the R.A.-Dec coordinates of the random catalogue will be the same as the ones of the real catalogue
    for (int i=0; i<max(0, int(N_R)-1); i++) { // use N_R times more random objects
      ra.insert(ra.end(), raB.begin(), raB.end());
      dec.insert(dec.end(), decB.begin(), decB.end());
    }
  }

  if (type==RandomType::_createRandom_shuffleTOT_) {

    random::UniformRandomNumbers ran(0., catalogue.nObjects()-1, seed);
    
    // construct the random sample
    for (size_t i=0; i<ra.size(); i++) {
      observedCoordinates coord = {ra[i], dec[i], round_to_precision(catalogue.redshift(ran()), z_ndigits)};
      m_object.push_back(move(Object::Create(ObjectType::_Random_, coord, cosm)));
    }
  }

  else if (type==RandomType::_createRandom_shuffle_) {

    // build the preliminary redshift bin edges
    const double zMin = catalogue.Min(cbl::catalogue::Var::_Redshift_);
    const double zMax = catalogue.Max(cbl::catalogue::Var::_Redshift_);
    std::vector<double> z_binEdges = {zMin};

    double z_position = zMin;
    while (z_position <= zMax) {
      z_position += z_step;
      z_binEdges.emplace_back(z_position);
    }

    // find the redshift bins containing at least one object in the data sample
    std::vector<bool> good_zbin(z_binEdges.size()-1,false);
    for (size_t j=0; j<z_binEdges.size()-1; j++)    
      for (size_t i=0; i<catalogue.nObjects(); i++)
        if (catalogue.redshift(i)>=z_binEdges[j] && catalogue.redshift(i)<z_binEdges[j+1]) {
          good_zbin[j] = true;
          break;    
        }

    std::vector<double> good_z_binEdges;
    for (size_t j=0; j<z_binEdges.size()-1; j++) {
      if (good_zbin[j]) {

        const bool found = std::find(good_z_binEdges.begin(), good_z_binEdges.end(), z_binEdges[j]) != good_z_binEdges.end();

        if (found==false)
          good_z_binEdges.emplace_back(z_binEdges[j]);
        good_z_binEdges.emplace_back(z_binEdges[j+1]);
      }    
    }
  
    // convolve the redshift distribution of the data
    std::vector<double> z_preliminary = catalogue.var(Var::_Redshift_);
    vector<double> xx, yy, err;
    distribution(xx, yy, err, z_preliminary, {}, nbins, true, par::defaultString, 1., zMin, zMax, "Linear", true, sigma_kernel);

    // loop until the random redshift distribution is ok
    std::vector<double> z;
    int times = times_default;
    bool ok_random = false, been_higher = false, been_lower = false;

    coutCBL << "I'm searching for the best number of extractions from the smoothed n(z)..." << endl;
    int i_counter = 0;

    while (ok_random==false) {
      z.resize(0);
      i_counter ++;
      if (i_counter == 100)
        ErrorCBL("Reached 100 iterations: I couldn't find a good way to build the random redshifts! Play with the following parameters: times_default, times_change, tolerance.", "Catalogue", "RandomCatalogue.cpp");
    
      // extract the random redshifts from the smoothed distribution
      z_preliminary.resize(ra.size()*times);
      z_preliminary = vector_from_distribution(ra.size()*times, xx, yy, zMin, zMax, seed);

      // select only the random redshifts in the good z bins
      for (size_t i=0; i<z_preliminary.size(); i++) {
        for (size_t j=0; j<good_z_binEdges.size()-1; j++) {
          if ( (good_z_binEdges[j+1]-good_z_binEdges[j]<1.5*z_step) && 
               (good_z_binEdges[j]<=z_preliminary[i]) && (z_preliminary[i]<good_z_binEdges[j+1]) ) {
            z.emplace_back(z_preliminary[i]);
          }      
        }
      }
    
      if (z.size() < ra.size()) {
        been_lower = true;
        times = (int)(times*(1+times_change));
      } 
      else if (z.size() > ra.size()*(1+tolerance)) {
        been_higher = true;
        times = (int)(times*(1-times_change));
      } else
      ok_random = true;
      if (been_higher && been_lower) {
        been_lower = false;
        been_higher = false;
        times_change -= 0.1;
        if (times_change<=0.015)
          ErrorCBL("I couldn't find a good way to build the random redshifts! Play with the following parameters: times_default, times_change, tolerance.", "Catalogue", "RandomCatalogue.cpp");
      }
    }
    coutCBL << "I've built the random redshifts vector." << endl;

    // uniformly erase the random redshifts in excess
    const int nExcess = (int)(z.size()-ra.size());

    for (int i=0; i<nExcess; i++) {
      random::UniformRandomNumbers_Int rand(0., z.size()-1, seed);
      z.erase(z.begin()+rand());
    }

    // If requested, write the n(z) file
    if (out_path_nz != par::defaultString && out_file_nz != par::defaultString) {

      // Normalise the original smoothed distribution (without holes)
      std::vector<double> yy_norm (yy.size());
      double yy_sum = 0, bin_width = xx[1]-xx[0];

      for (size_t j=0; j<xx.size(); j++)
        yy_sum += yy[j];

      for (size_t j=0; j<xx.size(); j++)
        yy_norm[j] = yy[j]/yy_sum/bin_width;

      // Build the normalised histogram of the extracted n(z) (without holes)      
      glob::Histogram1D hist_obj;
      hist_obj.set(nbins, xx[0], xx[xx.size()-1], 0.5, BinType::_linear_);

      std::vector<double> weight (z_preliminary.size(), 1.);
      hist_obj.put(z_preliminary, weight);

      std::vector<double> hist = hist_obj.operator()(glob::HistogramType::_N_V_);

      yy_sum = 0;
      for (size_t j=0; j<hist.size(); j++)
        yy_sum += hist[j];

      std::vector<double> yy2_norm (hist.size());
      for (size_t j=0; j<hist.size(); j++)
        yy2_norm[j] = hist[j]/yy_sum/bin_width;

      // Build the normalised histogram of the final n(z) (with holes)      
      hist_obj.set(nbins, xx[0], xx[xx.size()-1], 0.5, BinType::_linear_);
      weight.resize(z.size(), 1.);
      hist_obj.put(z, weight);

      hist = hist_obj.operator()(glob::HistogramType::_N_V_);

      yy_sum = 0;
      for (size_t j=0; j<hist.size(); j++)
        yy_sum += hist[j];

      std::vector<double> yy3_norm (hist.size());
      for (size_t j=0; j<hist.size(); j++)
        yy3_norm[j] = hist[j]/yy_sum/bin_width;
    
      // Write the file
      std::string mkdir = "mkdir -p "+out_path_nz;
      if (system(mkdir.c_str())) {}

      std::ofstream myfile;
      myfile.open(out_path_nz+out_file_nz);
    
      for (size_t j=0; j<xx.size(); j++)
        myfile << xx[j] << std::setw(20) << yy_norm[j] << std::setw(20) << yy2_norm[j] << std::setw(20) << yy3_norm[j] << std::endl;
      myfile.close();
      coutCBL<<"I wrote the file "+out_path_nz+out_file_nz<<std::endl;     
    }

    // finalise the random Catalogue instance
    for (size_t i=0; i<ra.size(); i++) {
      observedCoordinates coord = {ra[i], dec[i], round_to_precision(z[i], z_ndigits)};
      m_object.push_back(move(Object::Create(ObjectType::_Random_, coord, cosm)));
    }
  }
  else
    ErrorCBL("the chosen random catalogue type is not allowed!", "Catalogue", "RandomCatalogue.cpp"); 
}


// ============================================================================


cbl::catalogue::Catalogue::Catalogue (const RandomType type, const Catalogue catalogue, const double N_R, const int nbin, const double Angle, const std::vector<double> redshift, const cosmology::Cosmology &cosm, const bool conv, const double sigma, const int seed)
{
  if (type!=RandomType::_createRandom_cone_) ErrorCBL("the random catalogue has to be conic!", "Catalogue", "RandomCatalogue.cpp");
  
  coutCBL << "I'm creating a random catalogue in a cone..." << endl;
  const int nRandom = int(catalogue.nObjects()*N_R);
  random::UniformRandomNumbers ran(0., 1., seed);
  
  // extract the redshift distribution  
  
  vector<double> xx, yy, err;
  distribution(xx, yy, err, redshift, {}, nbin, true, par::defaultString, 1., cbl::Min(redshift), cbl::Max(redshift), "Linear", conv, sigma);
  
  // extract the random redshifts from the distribution
  
  vector<double> random_redshift = vector_from_distribution(nRandom, xx, yy, catalogue.Min(Var::_Redshift_), catalogue.Max(Var::_Redshift_), 13);

  // set the comoving distances
  
  vector<double> DD;
  for (auto &&red : random_redshift) 
    DD.emplace_back(cosm.D_C(red));

  // constructing the random sample
  
  const double rMax = cosm.D_C(catalogue.Max(Var::_Redshift_))*tan(Angle);
  const double Dm = cbl::Min(DD)*0.999;
  const double DM = cbl::Max(DD)*1.001;
  const double fact = 1./DM;
  
  for (size_t i=0; i<DD.size(); ++i) {
    comovingCoordinates coord;
    bool OK = false;
    while (!OK) {

      coord.xx = fact*DD[i]*rMax*2.*(ran()-0.5);
      coord.yy = DD[i];
      coord.zz = fact*DD[i]*rMax*2.*(ran()-0.5);

      const double r1 = sqrt(coord.xx*coord.xx+coord.zz*coord.zz);
      const double r2 = sqrt(coord.xx*coord.xx+coord.zz*coord.zz+DD[i]*DD[i]);
      
      if (r1<fact*DD[i]*rMax && Dm<r2 && r2<DM) OK = true;
    }

    m_object.push_back(move(Object::Create(ObjectType::_Random_, coord)));
  }

  computePolarCoordinates(); 
  
}


// ============================================================================


cbl::catalogue::Catalogue::Catalogue (const RandomType type, const std::vector<std::string> mangle_mask, const Catalogue catalogue, const double N_R, const int nbin, const cosmology::Cosmology cosm, const bool conv, const double sigma, const int seed)
{
  if (type!=RandomType::_createRandom_MANGLE_) ErrorCBL("the random catalogue has to be of type _MANGLE_!", "Catalogue", "RandomCatalogue.cpp");

  cbl::Path path;
  string mangle_dir = path.DirCosmo()+"/External/mangle/";

  string mangle_working_dir = mangle_dir+"output/";
  string mkdir = "mkdir -p "+mangle_working_dir;
  if (system(mkdir.c_str())) {}

  string mangle_mask_list;
  for (size_t i=0; i<mangle_mask.size(); i++)
    mangle_mask_list += mangle_mask[i]+" ";

  string random_output = mangle_working_dir+"temporary_random_output.dat";


  // generate the angular random distribution with mangle

  double nrandom = N_R*catalogue.nObjects();
  string nrandom_str = cbl::conv(nrandom, par::fDP1);

  string ransack = mangle_dir+"/bin/ransack -r"+nrandom_str+" "+mangle_mask_list+" "+random_output;
  if (system(ransack.c_str())) {}

  vector<double> random_ra, random_dec, random_redshift;
  
  ifstream fin(random_output.c_str()); checkIO(fin, random_output);
  string line;

  getline(fin, line);

  while (getline(fin, line)) {
    stringstream ss(line);
    double ra, dec;
    ss >> ra; ss >> dec;
    random_ra.push_back(ra);
    random_dec.push_back(dec);
  }
  fin.clear(); fin.close();

  
  // generate the redshift distribution from the input catalogue
    
  vector<double> redshift = catalogue.var(Var::_Redshift_);
  vector<double> xx, yy, err;
  distribution(xx, yy, err, redshift, {}, nbin, true, par::defaultString, 1., cbl::Min(redshift), cbl::Max(redshift), "Linear", conv, sigma);

  
  // extract the random redshifts from the distribution
  
  random_redshift = vector_from_distribution(random_ra.size(), xx, yy, catalogue.Min(Var::_Redshift_), catalogue.Max(Var::_Redshift_), seed);

  for (size_t i =0; i<random_ra.size(); i++) {
    observedCoordinates coord = {radians(random_ra[i]), radians(random_dec[i]), random_redshift[i]};
    m_object.push_back(move(Object::Create(ObjectType::_Random_, coord, cosm)));
  }
  
}


// ============================================================================


cbl::catalogue::Catalogue::Catalogue (const RandomType type, const Catalogue catalogue, const double N_R, const bool dndz_per_stripe, const int nbin, const cosmology::Cosmology cosm, const bool conv, const double sigma, const int seed)
{
  if (type!=RandomType::_createRandom_SDSS_stripes_) ErrorCBL("the random catalogue has to be of type _SDSS_stripe_!", "Catalogue", "RandomCatalogue.cpp");
  
  vector<double> lambda, eta;
  vector<int> stripe, stripe_list;

  vector<double> ra, dec;

  for (size_t i=0; i<catalogue.nObjects(); i++) {
    ra.push_back(catalogue.ra(i));
    dec.push_back(catalogue.dec(i));
  }

  eq2sdss(ra, dec, lambda, eta);

  sdss_stripe(eta, lambda, stripe, stripe_list);

  std::map<int, vector<double>> lambda_stripe, eta_stripe, redshift_stripe;
  std::map<int, int> nObj_stripe;

  for (size_t i=0; i<catalogue.nObjects(); i++) {
    nObj_stripe[stripe[i]] ++;
    lambda_stripe[stripe[i]].push_back(lambda[i]);
    eta_stripe[stripe[i]].push_back(eta[i]);
    if (dndz_per_stripe)
      redshift_stripe[stripe[i]].push_back(catalogue.redshift(i));
  }
    
  
  // extract random points
  
  random::UniformRandomNumbers random(0., 1., seed);

  vector<double> random_lambda, random_eta, random_ra, random_dec, random_redshift;

  for (auto &&ss : stripe_list) {
    const double MinL = sin(cbl::Min(lambda_stripe[ss])*par::pi/180);
    const double MaxL = sin(cbl::Max(lambda_stripe[ss])*par::pi/180);
    const double deltaL = MaxL-MinL;
    const double MinE = cbl::Min(eta_stripe[ss]);
    const double MaxE = cbl::Max(eta_stripe[ss]);
    const double deltaE = MaxE-MinE;

    const int nObj = lambda_stripe[ss].size();
    const int nRan = int(N_R*nObj);

    for (int j=0; j<nRan; j++) {
      random_lambda.push_back(asin(MinL+deltaL*random())*180/par::pi);
      random_eta.push_back(MinE+deltaE*random());
      if (dndz_per_stripe)
        random_redshift.push_back(redshift_stripe[ss][int(random()*(nObj-1))]);
    }
  }

  sdss2eq(random_lambda, random_eta, random_ra, random_dec);

  if (!dndz_per_stripe) {
    // generate the redshift distribution from the input catalogue

    vector<double> redshift = catalogue.var(Var::_Redshift_);
    vector<double> xx, yy, err;
    distribution(xx, yy, err, redshift, {}, nbin, true, par::defaultString, 1., cbl::Min(redshift), cbl::Max(redshift), "Linear", conv, sigma);

    // extract the random redshifts from the distribution

    random_redshift = vector_from_distribution(random_ra.size(), xx, yy, catalogue.Min(Var::_Redshift_), catalogue.Max(Var::_Redshift_), seed);
  }

  for (size_t i =0; i<random_ra.size(); i++) {
    observedCoordinates coord = {radians(random_ra[i]), radians(random_dec[i]), random_redshift[i]};
    m_object.push_back(move(Object::Create(ObjectType::_Random_, coord, cosm)));
  }
}

//=========================================================================

cbl::catalogue::Catalogue::Catalogue(const RandomType type, Catalogue catalogue, const double N_R, const cosmology::Cosmology cosm, const std::vector<double> RA_range, const std::vector<double> DEC_range, const unsigned int nbin, const int seed)
{
  if (type != RandomType::_createRandom_homogeneous_LC_)
    ErrorCBL("the random catalogue has to be of type _homogeneous_LC_!", "Catalogue", "RandomCatalogue.cpp");
  
  auto prop = catalogue.compute_catalogueProperties_lightCone(cosm, RA_range, DEC_range, nbin);
  vector<double> nObj(nbin+1);
  nObj[0] = 1;
  vector<double> zbins = linear_bin_vector(nbin+1, catalogue.Min(Var::_Redshift_), catalogue.Max(Var::_Redshift_));
  for (size_t i=1; i<nbin+1; i++) 
    nObj[i]=nObj[i-1]+prop[1][i-1];
  
  std::default_random_engine generator;
  generator.seed(seed);
  std::uniform_int_distribution<int> distribution(1,catalogue.nObjects());
  random::UniformRandomNumbers randomDec(sin(DEC_range[0]), sin(DEC_range[1]), seed+1);
  random::UniformRandomNumbers randomRA(RA_range[0], RA_range[1], seed+2);

  for (int i=0; i<int(catalogue.nObjects()*N_R); i++) {
    observedCoordinates coord = {randomRA(), asin(randomDec()), interpolated(distribution(generator), nObj, zbins, "Linear")};
    m_object.push_back(move(Object::Create(ObjectType::_Random_, coord, cosm)));
  }

}

//=========================================================================

void cbl::catalogue::Catalogue::equalize_random_lightCone (cbl::catalogue::Catalogue tracer_catalogue, cbl::cosmology::Cosmology cosm, const std::vector<double> RA_range, const std::vector<double> DEC_range, const int seed) 
{
  if (nObjects()>tracer_catalogue.nObjects()) {
    std::default_random_engine generator;
    generator.seed(seed);
    std::uniform_int_distribution<int> distribution(0,nObjects()-1);

    vector<bool> mask(nObjects(), false);

    size_t erased=0;
    while (erased<nObjects()-tracer_catalogue.nObjects()) {
      int rand = distribution(generator);
      if (mask[rand] == false) {
        mask[rand] = true;
        erased++;
      }
    }

    remove_objects(mask);
  }
  else {

    auto prop = tracer_catalogue.compute_catalogueProperties_lightCone(cosm, RA_range, DEC_range, 100);
    vector<double> nObj = {1.};
    vector<double> zbins = linear_bin_vector(101, cbl::Min(var(Var::_Redshift_)), cbl::Max(var(Var::_Redshift_)));
    for (size_t i=0; i<100; i++) 
      nObj.emplace_back(prop[1][i]);

    std::default_random_engine generator;
    generator.seed(seed);
    std::uniform_int_distribution<int> distribution(1,tracer_catalogue.nObjects());
    random::UniformRandomNumbers randomDec(sin(DEC_range[0]), sin(DEC_range[1]), seed+1);
    random::UniformRandomNumbers randomRA(RA_range[0], RA_range[1], seed+2);

    int nObjec = nObjects();
    for (int i=0; i<int(tracer_catalogue.nObjects()-nObjec); i++) {
      observedCoordinates coord = {randomRA(), asin(randomDec()), interpolated(distribution(generator), nObj, zbins, "Linear")};
      add_object(move(Object::Create(ObjectType::_Random_, coord, cosm)));
    }
  }
}

//=========================================================================

void cbl::catalogue::Catalogue::equalize_random_box (cbl::catalogue::Catalogue tracer_catalogue, const int seed) 
{
  if (nObjects()>tracer_catalogue.nObjects()) {
    std::default_random_engine generator;
    generator.seed(seed);
    std::uniform_int_distribution<int> distribution(0,nObjects()-1);

    vector<bool> mask(nObjects(), false);

    size_t erased=0;
    while (erased<nObjects()-tracer_catalogue.nObjects()) {
      int rand = distribution(generator);
      if (mask[rand] == false) {
        mask[rand] = true;
        erased++;
      }
    }

    remove_objects(mask);
  }
  else {
    random::UniformRandomNumbers randomX(tracer_catalogue.Min(Var::_X_), tracer_catalogue.Max(Var::_X_), seed);
    random::UniformRandomNumbers randomY(tracer_catalogue.Min(Var::_Y_), tracer_catalogue.Max(Var::_Y_), seed+1);
    random::UniformRandomNumbers randomZ(tracer_catalogue.Min(Var::_Z_), tracer_catalogue.Max(Var::_Z_), seed+2);

    int nObjec = nObjects();
    for (int i=0; i<int(tracer_catalogue.nObjects()-nObjec); i++) {
      comovingCoordinates coord = {randomX(), randomY(), randomZ()};
      add_object(move(Object::Create(ObjectType::_Random_, coord)));
    }
  }
}
