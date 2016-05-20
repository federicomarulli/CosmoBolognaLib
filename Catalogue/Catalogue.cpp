/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Catalogue/Catalogue.cpp
 *
 *  @brief Methods of the class Catalogue 
 *
 *  This file contains the implementation of the methods of the class
 *  Catalogue, used to handle catalogues of astronomical sources
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Field3D.h"
#include "Catalogue.h"

using namespace cosmobl;
using namespace catalogue;


// ============================================================================

/// @cond template

template cosmobl::catalogue::Catalogue::Catalogue (vector<cosmobl::catalogue::RandomObject>);
template cosmobl::catalogue::Catalogue::Catalogue (vector<cosmobl::catalogue::Mock>);
template cosmobl::catalogue::Catalogue::Catalogue (vector<cosmobl::catalogue::Halo>);
template cosmobl::catalogue::Catalogue::Catalogue (vector<cosmobl::catalogue::Galaxy>);
template cosmobl::catalogue::Catalogue::Catalogue (vector<cosmobl::catalogue::Cluster>);
template cosmobl::catalogue::Catalogue::Catalogue (vector<cosmobl::catalogue::Void>);

template void cosmobl::catalogue::Catalogue::add_object (cosmobl::catalogue::RandomObject);
template void cosmobl::catalogue::Catalogue::add_object (cosmobl::catalogue::Mock);
template void cosmobl::catalogue::Catalogue::add_object (cosmobl::catalogue::Halo);
template void cosmobl::catalogue::Catalogue::add_object (cosmobl::catalogue::Galaxy);
template void cosmobl::catalogue::Catalogue::add_object (cosmobl::catalogue::Cluster);
template void cosmobl::catalogue::Catalogue::add_object (cosmobl::catalogue::Void);

template void cosmobl::catalogue::Catalogue::add_objects (vector<cosmobl::catalogue::RandomObject>);
template void cosmobl::catalogue::Catalogue::add_objects (vector<cosmobl::catalogue::Mock>);
template void cosmobl::catalogue::Catalogue::add_objects (vector<cosmobl::catalogue::Halo>);
template void cosmobl::catalogue::Catalogue::add_objects (vector<cosmobl::catalogue::Galaxy>);
template void cosmobl::catalogue::Catalogue::add_objects (vector<cosmobl::catalogue::Cluster>);
template void cosmobl::catalogue::Catalogue::add_objects (vector<cosmobl::catalogue::Void>);

template void cosmobl::catalogue::Catalogue::remove_objects (vector<cosmobl::catalogue::RandomObject>);
template void cosmobl::catalogue::Catalogue::remove_objects (vector<cosmobl::catalogue::Mock>);
template void cosmobl::catalogue::Catalogue::remove_objects (vector<cosmobl::catalogue::Halo>);
template void cosmobl::catalogue::Catalogue::remove_objects (vector<cosmobl::catalogue::Galaxy>);
template void cosmobl::catalogue::Catalogue::remove_objects (vector<cosmobl::catalogue::Cluster>);
template void cosmobl::catalogue::Catalogue::remove_objects (vector<cosmobl::catalogue::Void>);

/// @endcond

// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const ObjType type, const vector<double> xx, const vector<double> yy, const vector<double> zz, vector<double> weight)
{
  if (weight.size() == 0)
    weight.resize(xx.size(), 1);

  if (!(xx.size()==yy.size() && yy.size()==zz.size() && zz.size()==weight.size()))
    ErrorMsg("Error in Catalogue::Catalogue() of Catalogue.cpp: coordinates with different dimension!"); 
  
  for (size_t i=0; i<xx.size(); i++)
    m_sample.push_back(move(Object::Create(type, xx[i], yy[i], zz[i], weight[i])));
}


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const Catalogue& cat)
{
  m_sample = cat.m_sample;
  m_index = cat.m_index;
}


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const ObjType type, const vector<double> ra, const vector<double> dec, const vector<double> redshift, const Cosmology &cosm, const CoordUnits inputUnits, vector<double> weight)
{
  // conversion of coordinate units into radiants
  vector<double> RA(ra.size()), DEC(dec.size());
  for (size_t i=0; i<ra.size(); i++) 
    RA[i] = (inputUnits==_radians_) ? ra[i] : radians(ra[i], inputUnits);
  for (size_t i=0; i<dec.size(); i++) 
    DEC[i] = (inputUnits==_radians_) ? dec[i] : radians(dec[i], inputUnits);
 
  // (weight=1 if weights are not used)
  if (weight.size()==0)
    weight.resize(RA.size(), 1);

  // check the vector dimensions
  if (!(RA.size()==DEC.size() && DEC.size()==redshift.size() && redshift.size()==weight.size()))
    ErrorMsg("Error in Catalogue::Catalogue() of Catalogue.cpp: coordinates with different dimension!"); 

  // include the objects in the catalogue
  for (size_t i=0; i<RA.size(); i++)
    m_sample.push_back(move(Object::Create(type, RA[i], DEC[i], redshift[i], cosm, weight[i])));

}
  


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const ObjType type, const vector<string> file, const int col_X, const int col_Y, const int col_Z, const int col_Weight, const double nSub) 
{
  default_random_engine gen;
  uniform_real_distribution<float> ran(0., 1.);
    
  string line; double XX, YY, ZZ, WW, VV;
    
  for (size_t dd=0; dd<file.size(); dd++) {

    string file_in = file[dd];
    cout << "I'm reading the catalogue: " << file_in << endl;
    ifstream finr (file_in.c_str()); checkIO(file_in, 1);
    
    while (getline(finr, line)) {
      if (ran(gen)<nSub) {
	stringstream ss(line);
	vector<double> val; while (ss>>VV) val.push_back(VV);
	checkDim(val, col_X, "val", 0); checkDim(val, col_Y, "val", 0); checkDim(val, col_Z, "val", 0);
	XX = val[col_X]; YY = val[col_Y]; ZZ = val[col_Z];
	WW = (col_Weight!=-1) ? val[col_Weight] : 1.;
        m_sample.push_back(move(Object::Create(type, XX, YY, ZZ, WW)));
      }
    }
    
    finr.clear(); finr.close(); 
  }
}


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const ObjType type, const vector<string> file, const Cosmology &cosm, const int col_RA, const int col_Dec, const int col_redshift, const int col_Weight, const CoordUnits inputUnits, const double nSub, const double fact) 
{ 
  default_random_engine gen;
  uniform_real_distribution<float> ran(0., 1.);
 
  string line; double RA, DEC, ZZ, WW, VV;
  
  for (size_t dd=0; dd<file.size(); dd++) {

    string file_in = file[dd];
    cout << "I'm reading the catalogue: " << file_in << endl;
    ifstream finr (file_in.c_str()); checkIO(file_in, 1);

    while (getline(finr, line)) {
      if (ran(gen)<nSub) {
	stringstream ss(line);
	vector<double> val; while (ss>>VV) val.push_back(VV);
	checkDim(val, col_RA, "val", 0); checkDim(val, col_Dec, "val", 0); checkDim(val, col_redshift, "val", 0);
	RA = (inputUnits==_radians_) ? val[col_RA]*fact : radians(val[col_RA], inputUnits)*fact;
	DEC = (inputUnits==_radians_) ? val[col_Dec]*fact : radians(val[col_Dec], inputUnits)*fact;
	ZZ = val[col_redshift];
	WW = (col_Weight!=-1) ? val[col_Weight] : 1.;
	m_sample.push_back(move(Object::Create(type, RA, DEC, ZZ, cosm, WW)));
      }
    }
    
    finr.clear(); finr.close(); 
  }
}


// ============================================================================


int cosmobl::catalogue::Catalogue::Nregion () const
{ 
  vector<int> regions; 
  for (int i=0; i<nObjects(); i++) regions.push_back(m_sample[i]->region()); 
  sort(regions.begin(), regions.end());
  vector<int>::iterator it = unique(regions.begin(), regions.end());
  regions.resize(std::distance(regions.begin(), it));
  return regions.size();
}


// ============================================================================


vector<long> cosmobl::catalogue::Catalogue::get_region_list () const
{ 
  vector<long> regions; 
  for (int i=0; i<nObjects(); i++) regions.push_back(m_sample[i]->region()); 
  sort(regions.begin(), regions.end());
  vector<long>::iterator it = unique(regions.begin(), regions.end());
  regions.resize(std::distance(regions.begin(), it));
  return regions;
}


// ============================================================================


vector<double> cosmobl::catalogue::Catalogue::var (Var var_name) const
{ 
  vector<double> vv(m_sample.size(), 0.);
  
  switch (var_name) {

  case Var::_X_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->xx();
    break;

  case Var::_Y_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->yy();
    break;

  case Var::_Z_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->zz();
    break;

  case Var::_RA_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->ra();
    break;

  case Var::_Dec_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->dec();
    break;

  case Var::_Redshift_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->redshift();
    break;

  case Var::_Dc_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->dc();
    break;

  case Var::_Weight_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->weight();
    break;

  case Var::_Mass_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->mass();
    break;

  case Var::_Magnitude_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->magnitude();
    break;

  case Var::_Richness_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->richness();
    break;

  case Var::_Vx_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->vx();
    break;
  
  case Var::_Vy_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->vy();
    break;
  
  case Var::_Vz_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->vz();
    break;

  case Var::_Region_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->region();
    break; 
  
  case Var::_Generic_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->generic();
    break;

  case Var::_Radius_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->radius();
    break;

  default:
    ErrorMsg("Error in cosmobl::catalogue::Catalogue::var of Catalogue.cpp: no such a variable in the list!");
  }
  
  return vv;
}


// ============================================================================


void cosmobl::catalogue::Catalogue::set_var (const Var var_name, const vector<double> _var)
{
  if (m_sample.size()!=_var.size()) ErrorMsg ("Error in cosmobl::catalogue::Catalogue::set_var of Catalogue.cpp: different sizes!");
  
  switch (var_name) {

  case Var::_X_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_xx(_var[i]);
    break;

  case Var::_Y_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_yy(_var[i]);
    break;

  case Var::_Z_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_zz(_var[i]);
    break;

  case Var::_RA_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_ra(_var[i]);
    break;

  case Var::_Dec_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_dec(_var[i]);
    break;

  case Var::_Redshift_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_redshift(_var[i]);
    break;

  case Var::_Dc_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_dc(_var[i]);
    break;

  case Var::_Weight_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_weight(_var[i]);
    break;

  case Var::_Mass_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_mass(_var[i]);
    break;

  case Var::_Magnitude_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_magnitude(_var[i]);
    break;

  case Var::_Richness_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_richness(_var[i]);
    break;

  case Var::_Vx_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_vx(_var[i]);
    break;
  
  case Var::_Vy_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_vy(_var[i]);
    break;
  
  case Var::_Vz_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_vz(_var[i]);
    break;

  case Var::_Region_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_region(_var[i]);
    break;

  case Var::_Generic_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_generic(_var[i]);
    break;

  case Var::_Radius_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_radius(_var[i]);
    break;

  default:
    ErrorMsg("Error in cosmobl::catalogue::Catalogue::set_var of Catalogue.cpp: no such a variable in the list!");
  }

}


// ============================================================================


void cosmobl::catalogue::Catalogue::stats_var (const Var var_name, vector<double> &stats) const
{
  stats.erase(stats.begin(), stats.end());
  stats.resize(4);
  
  stats[0] = Average(var(var_name)); 
  stats[2] = Sigma(var(var_name));
  
  stats[1] = Quartile(var(var_name))[1];
  stats[3] = Quartile(var(var_name))[2]-Quartile(var(var_name))[0];
}


// ============================================================================


void cosmobl::catalogue::Catalogue::stats_var (const vector<Var> var_name, vector<vector<double> > &stats) const
{
  stats.erase(stats.begin(),stats.end());

  for (unsigned int i=0; i<var_name.size(); i++) {

    vector<double> stats_temp;
    stats_var(var_name[i],stats_temp);

    stats.push_back(stats_temp);
  }
}


// ============================================================================


void cosmobl::catalogue::Catalogue::var_distr (const Var var_name, vector<double> &_var, vector<double> &dist, vector<double> &err, const int nbin, const bool linear, const string file_out, const double Volume, const bool norm, const double V1, const double V2, const bool bin_type, const bool convolution, const double sigma) const
{
  distribution(_var, dist, err, var(var_name), var(Var::_Weight_), nbin, linear, file_out, (norm) ? Volume*weightedN() : Volume, V1, V2, bin_type, convolution, sigma);
}


// ============================================================================


void cosmobl::catalogue::Catalogue::computeComovingCoordinates (const Cosmology &cosm, const CoordUnits inputUnits)
{
  double red, xx, yy, zz;

  
  // ----- unit conversion -----
  
  vector<double> RA(nObjects()), DEC(nObjects());

  for (int i=0; i<nObjects(); i++) {
    RA[i] = (inputUnits==_radians_) ? ra(i) : radians(ra(i), inputUnits);
    DEC[i] = (inputUnits==_radians_) ? dec(i) : radians(dec(i), inputUnits);
  }

  
  // ----- compute comoving coordinates -----
  
  for (int i=0; i<nObjects(); i++) {

    red = redshift(i);
    m_sample[i]->set_dc(cosm.D_C(red));
    
    cartesian_coord(RA[i], DEC[i], dc(i), xx, yy, zz);
    
    m_sample[i]->set_xx(xx); 
    m_sample[i]->set_yy(yy); 
    m_sample[i]->set_zz(zz);
  }
}


// ============================================================================


void cosmobl::catalogue::Catalogue::computePolarCoordinates (const CoordUnits outputUnits)
{
  double ra, dec, dc;

  
  // ----- compute polar coordinates in radians -----
  
  for (int i=0; i<nObjects(); i++) {
    polar_coord(xx(i), yy(i), zz(i), ra, dec, dc);
    m_sample[i]->set_ra(ra); 
    m_sample[i]->set_dec(dec); 
    m_sample[i]->set_dc(dc);
  }

  
  // ----- unit conversion -----
  
  if (outputUnits!=_radians_) {

    if (outputUnits==_degrees_)
      for (int i=0; i<nObjects(); i++) {
	m_sample[i]->set_ra(degrees(ra));
	m_sample[i]->set_dec(degrees(dec)); 		    
      }

    else if (outputUnits==_arcseconds_)
      for (int i=0; i<nObjects(); i++) {
	m_sample[i]->set_ra(arcseconds(ra));
	m_sample[i]->set_dec(arcseconds(dec)); 		    
      }

    else if (outputUnits==_arcminutes_)
      for (int i=0; i<nObjects(); i++) {
	m_sample[i]->set_ra(arcminutes(ra));
	m_sample[i]->set_dec(arcminutes(dec)); 		    
      }

    else ErrorMsg("Error in cosmobl::catalogue::Catalogue::computePolarCoordinates: outputUnits type not allowed!");
  }
  
}

// ============================================================================


void cosmobl::catalogue::Catalogue::computePolarCoordinates (const Cosmology &cosm, const double z1, const double z2, const CoordUnits outputUnits)
{
  double ra, dec, dc;

  
  // ----- compute polar coordinates in radians -----

  for (int i=0; i<nObjects(); i++) {
    polar_coord(xx(i), yy(i), zz(i), ra, dec, dc);
    m_sample[i]->set_ra(ra); 
    m_sample[i]->set_dec(dec); 
    m_sample[i]->set_dc(dc);
    m_sample[i]->set_redshift(cosm.Redshift(dc, z1, z2));
  }

  
  // ----- unit conversion -----

  if (outputUnits!=_radians_) {

    if (outputUnits==_degrees_)
      for (int i=0; i<nObjects(); i++) {
	m_sample[i]->set_ra(degrees(ra));
	m_sample[i]->set_dec(degrees(dec)); 		    
      }

    else if (outputUnits==_arcseconds_)
      for (int i=0; i<nObjects(); i++) {
	m_sample[i]->set_ra(arcseconds(ra));
	m_sample[i]->set_dec(arcseconds(dec)); 		    
      }

    else if (outputUnits==_arcminutes_)
      for (int i=0; i<nObjects(); i++) {
	m_sample[i]->set_ra(arcminutes(ra));
	m_sample[i]->set_dec(arcminutes(dec)); 		    
      }

    else ErrorMsg("Error in cosmobl::catalogue::Catalogue::computePolarCoordinates: outputUnits type not allowed!");
  }
  
}

// ============================================================================


void cosmobl::catalogue::Catalogue::normalizeComovingCoordinates () 
{
  for (int i=0; i<nObjects(); i++) { 
    m_sample[i]->set_xx(xx(i)/dc(i)); 
    m_sample[i]->set_yy(yy(i)/dc(i)); 
    m_sample[i]->set_zz(zz(i)/dc(i));
  }
}


// ============================================================================


void cosmobl::catalogue::Catalogue::restoreComovingCoordinates ()
{
  for (int i=0; i<nObjects(); i++) {
    m_sample[i]->set_xx(xx(i)*dc(i)); 
    m_sample[i]->set_yy(yy(i)*dc(i)); 
    m_sample[i]->set_zz(zz(i)*dc(i));
  }
}


// ============================================================================


void cosmobl::catalogue::Catalogue::Order (const vector<int> vv) 
{
  int nObj = m_sample.size();

  if (int(vv.size())!=nObj) ErrorMsg("Error in Catalogue::Order()!");
 
  vector<shared_ptr<Object>> obj(nObj);
  
  m_index.resize(nObj);
  
  for (size_t i=0; i<vv.size(); i++) {
    m_index[i] = vv[i];
    obj[i] = m_sample[vv[i]];
  }

  m_sample = obj;
}


// ============================================================================


void cosmobl::catalogue::Catalogue::Order () 
{ 
  size_t nObj = m_sample.size();
  
  vector<shared_ptr<Object>> obj(nObj);
  
  if (m_index.size() != nObj) 
    ErrorMsg("Error in Catalogue::Order() of Catalogue.cpp, order not found!");
  
  obj = m_sample;
  
  for (size_t i=0; i<nObj; i++) 
    m_sample[i] = obj[m_index[i]];
}


// ============================================================================


double cosmobl::catalogue::Catalogue::weightedN () const
{
  double nn = 0.;
  for (unsigned int i=0; i<m_sample.size(); i++)
    nn += m_sample[i]->weight();
  return nn;
}


// ============================================================================


void cosmobl::catalogue::Catalogue::write_comoving_coordinates (const string file_output) const
{
  if (m_sample.size()==0) ErrorMsg("Error in cosmobl::catalogue::Catalogue::write_comoving_coordinates: m_sample.size()=0!");

  ofstream fout(file_output.c_str()); checkIO(file_output, 0);
  
  for (int i=0; i<nObjects(); i++) 
    fout << xx(i) << "   " << yy(i) << "   " << zz(i) << endl;

  cout << "I wrote the file: " << file_output << endl;
  fout.clear(); fout.close();
}


// ============================================================================


void cosmobl::catalogue::Catalogue::write_obs_coordinates (const string file_output) const 
{
  if (m_sample.size()==0) ErrorMsg("Error in cosmobl::catalogue::Catalogue::write_obs_coordinates: m_sample.size()=0!");
  
  ofstream fout(file_output.c_str()); checkIO(file_output, 0);
  
  if (!isSet(ra(0)) || !isSet(dec(0)) || !isSet(redshift(0)))
    ErrorMsg("Error in cosmobl::catalogue::Catalogue::write_obs_coords of Catalogue.cpp! Polar coordinates are not set!");

  for (int i=0; i<nObjects(); i++) 
    fout << ra(i) << "   " << dec(i) << "   " << redshift(i) << endl;
  
  cout << "I wrote the file: " << file_output << endl;
  fout.clear(); fout.close();
}


// ============================================================================


void cosmobl::catalogue::Catalogue::write_coordinates (const string file_output) const 
{
  ofstream fout(file_output.c_str()); checkIO(file_output, 0);

  if (!isSet(ra(0)) || !isSet(dec(0)) || !isSet(redshift(0)) || !isSet(dc(0)))
    ErrorMsg("Error in cosmobl::catalogue::Catalogue::write_coordinates of Catalogue.cpp! Polar coordinates are not set!");
  
  for (int i=0; i<nObjects(); i++) 
    fout << xx(i) << "   " << yy(i) << "   " << zz(i) << "   " << ra(i) << "   " << dec(i) << "   " << redshift(i) << "   " << dc(i) << endl;
  
  cout << "I wrote the file: " << file_output << endl;
  fout.clear(); fout.close();
}


// ============================================================================


shared_ptr<Catalogue> cosmobl::catalogue::Catalogue::cut (const Var var_name, const double down, const double up, const bool excl)
{
  vector<shared_ptr<Object> > objects;
  vector<double> vvar = var(var_name);
  vector<int> w(vvar.size());

  for (size_t i=0; i<m_sample.size(); i++){
      w[i] = (excl) ? 1 : 0;
    if (vvar[i] >= down && vvar[i] < up)
      w[i] = (excl) ? 0 : 1;
  }

  for (size_t i=0; i<m_sample.size(); i++)
    if (w[i]==1)
      objects.push_back(m_sample[i]);
 
  shared_ptr<Catalogue> cat(new Catalogue{objects});
  return cat;
}


// ============================================================================


double cosmobl::catalogue::Catalogue::distance (const int i, shared_ptr<Object> obj) const
{
  return sqrt((m_sample[i]->xx()-obj->xx())*(m_sample[i]->xx()-obj->xx())+
	      (m_sample[i]->yy()-obj->yy())*(m_sample[i]->yy()-obj->yy())+
	      (m_sample[i]->zz()-obj->zz())*(m_sample[i]->zz()-obj->zz()));
}


// ============================================================================


double cosmobl::catalogue::Catalogue::angsep_xyz (const int i, shared_ptr<Object> obj) const
{ 
  return 2.*asin(0.5*sqrt((m_sample[i]->xx()-obj->xx())*(m_sample[i]->xx()-obj->xx())+
			  (m_sample[i]->yy()-obj->yy())*(m_sample[i]->yy()-obj->yy())+
			  (m_sample[i]->zz()-obj->zz())*(m_sample[i]->zz()-obj->zz())));
}
    

// ============================================================================


shared_ptr<Catalogue> cosmobl::catalogue::Catalogue::smooth (const double gridsize, const vector<Var> vars, const int SUB)
{
  shared_ptr<Catalogue> cat {new Catalogue(*this)};
  
  if (gridsize<1.e-30) return cat;
  
  double rMAX = 0.;
  
  vector<shared_ptr<Object>> sample;

  
  // ----------------------------------------------------------------------------
  // ----- subdivide the catalogue in SUB^3 sub-catalogues, to avoid ------------
  // ----- memory problems possibly caused by too large chain-mesh vectors ------
  // ----------------------------------------------------------------------------

  cout <<"Please wait, I'm subdividing the catalogue in "<<pow(SUB, 3)<<" sub-catalogues..."<<endl;
 
  int nx = SUB, ny = SUB, nz = SUB;
  
  double Cell_X = (Max(Var::_X_)-Min(Var::_X_))/nx;
  double Cell_Y = (Max(Var::_Y_)-Min(Var::_Y_))/ny;
  double Cell_Z = (Max(Var::_Z_)-Min(Var::_Z_))/nz;

  for (int i=0; i<cat->nObjects(); i++) {
    int i1 = min(int((cat->xx(i)-Min(Var::_X_))/Cell_X), nx-1);
    int j1 = min(int((cat->yy(i)-Min(Var::_Y_))/Cell_Y), ny-1);
    int z1 = min(int((cat->zz(i)-Min(Var::_Z_))/Cell_Z), nz-1);
    int index = z1+nz*(j1+ny*i1);
    cat->catalogue_object(i)->set_region(index);
  }
  
  vector<long> region_list = cat->get_region_list();
  int nRegions = region_list.size();
  
  vector<shared_ptr<Catalogue>> subSamples(nRegions);
  
  for (int i=0; i<nRegions; i++) {
    double start = region_list[i];
    double stop = start+1;
    subSamples[i] = cut(Var::_Region_, start, stop);
  }

  
  // -----------------------------------------
  // ----- smooth all the sub-catalogues -----
  // -----------------------------------------

  for (int rr=0; rr<nRegions; rr++) {
    
    vector<double> _xx = subSamples[rr]->var(Var::_X_), _yy = subSamples[rr]->var(Var::_Y_), _zz = subSamples[rr]->var(Var::_Z_);
    
    ChainMesh3D ll(gridsize, _xx, _yy, _zz, rMAX, (long)-1.e5, (long)1.e5);
   
    
    for (int i=0; i<ll.nCell(); i++) {
      vector<long> list = ll.get_list(i);
     
      int nObj = list.size();
      if (nObj>0) {
	
	double XX = 0., YY = 0., ZZ = 0., RA = 0., DEC = 0., REDSHIFT = 0., WEIGHT = 0.;

	for (size_t j=0; j<list.size(); j++) {
	  XX += subSamples[rr]->xx(list[j]);
	  YY += subSamples[rr]->yy(list[j]);
	  ZZ += subSamples[rr]->zz(list[j]);
	  RA += subSamples[rr]->ra(list[j]);
	  DEC += subSamples[rr]->dec(list[j]);
	  REDSHIFT += subSamples[rr]->redshift(list[j]);
	  WEIGHT += subSamples[rr]->weight(list[j]);
	}
	
	shared_ptr<Object> obj{new Object()};
	XX /= nObj; obj->set_xx(XX);
	YY /= nObj; obj->set_yy(YY);
	ZZ /= nObj; obj->set_zz(ZZ);
	RA /= nObj; obj->set_ra(RA);
	DEC /= nObj; obj->set_dec(DEC);
	REDSHIFT /= nObj; obj->set_redshift(REDSHIFT);
	obj->set_weight(WEIGHT);
	sample.push_back(obj);

      }
    }
    
  }
  
  shared_ptr<Catalogue> cat_new(new Catalogue{sample});
  return cat_new;
}


// ============================================================================


int cosmobl::catalogue::Catalogue::nObjects_condition (const Var var_name, const double down, const double up, const bool excl)
{
  int nObjw = 0;
  vector<double> vvar = var(var_name);

  for (size_t i=0; i<m_sample.size(); i++)

    if (vvar[i] >= down && vvar[i] < up) 
      nObjw ++;
    
  nObjw = (excl) ? weightedN()-nObjw : nObjw;

  return nObjw;
}


// ============================================================================


double cosmobl::catalogue::Catalogue::weightedN_condition (const Var var_name, const double down, const double up, const bool excl)
{
  double nObjw = 0;
  vector<double> vvar = var(var_name);

  for (size_t i=0; i<m_sample.size(); i++)
    if (vvar[i] >= down && vvar[i] < up)
      nObjw += weight(i);
    
  nObjw = (excl) ? weightedN()-nObjw : nObjw;

  return nObjw;
}


// ============================================================================


ScalarField3D cosmobl::catalogue::Catalogue::density_field (const double cell_size, const int interpolation_type, const double kernel_radius, const bool useMass) const
{
  ScalarField3D density(cell_size, Min(Var::_X_), Max(Var::_X_), Min(Var::_Y_), Max(Var::_Y_), Min(Var::_Z_), Max(Var::_Z_));

  double deltaX = density.deltaX();
  double deltaY = density.deltaY();
  double deltaZ = density.deltaZ();
  int nx = density.nx();
  int ny = density.ny();
  int nz = density.nz();

  long int nCells = density.nCells();
  double VolumeCell_inv = pow(deltaX*deltaY*deltaZ,-1);

  for(int i=0;i<nObjects();i++){
    int i1 = min(int((xx(i)-Min(Var::_X_))/deltaX),nx-1);
    int j1 = min(int((yy(i)-Min(Var::_Y_))/deltaY),ny-1);
    int k1 = min(int((zz(i)-Min(Var::_Z_))/deltaZ),nz-1);

    double w = (useMass) ? mass(i)*VolumeCell_inv : VolumeCell_inv;

    if (interpolation_type==0){
      density.set_ScalarField(w, i1, j1, k1,1);
    }
    else if(interpolation_type==1){
      double dx_samecell, dy_samecell, dz_samecell;
      int dx_index, dy_index, dz_index;

      dx_samecell = (xx(i)-(i1*deltaX+Min(Var::_X_)))/deltaX;
      dy_samecell = (yy(i)-(j1*deltaY+Min(Var::_Y_)))/deltaY;
      dz_samecell = (zz(i)-(k1*deltaZ+Min(Var::_Z_)))/deltaZ;

      dx_index = (dx_samecell <0.5) ? -1 : 1 ;
      dy_index = (dy_samecell <0.5) ? -1 : 1 ;
      dz_index = (dz_samecell <0.5) ? -1 : 1 ;
      
      dx_samecell = (dx_samecell <0.5) ? (xx(i)+0.5*deltaX-(i1*deltaX+Min(Var::_X_)))/deltaX: 1-(xx(i)-0.5*deltaX-(i1*deltaX+Min(Var::_X_)))/deltaX;
      dy_samecell = (dy_samecell <0.5) ? (yy(i)+0.5*deltaY-(j1*deltaY+Min(Var::_Y_)))/deltaY: 1-(yy(i)-0.5*deltaY-(j1*deltaY+Min(Var::_Y_)))/deltaY;
      dz_samecell = (dz_samecell <0.5) ? (zz(i)+0.5*deltaZ-(k1*deltaZ+Min(Var::_Z_)))/deltaZ: 1-(zz(i)-0.5*deltaZ-(k1*deltaZ+Min(Var::_Z_)))/deltaZ;
      
      if(dx_samecell<0 || dy_samecell<0 || dz_samecell<0){
	cout << dx_index << " " <<  dx_samecell << " " <<dy_index << " " << dy_samecell << " " <<dz_index << " " << dz_samecell << endl;
      }

      double ww = 0.;
      int i2 = ((i1+dx_index>-1) & (i1+dx_index<nx)) ? i1+dx_index : i1;
      int j2 = ((j1+dy_index>-1) & (j1+dy_index<ny)) ? j1+dy_index : j1;
      int k2 = ((k1+dz_index>-1) & (k1+dz_index<nz)) ? k1+dz_index : k1;
      
      ww += w*(dx_samecell*dy_samecell*dz_samecell);
      density.set_ScalarField(w*(dx_samecell*dy_samecell*dz_samecell), i1, j1, k1, 1);
      
      ww += w*(dx_samecell*dy_samecell*(1.-dz_samecell));
      density.set_ScalarField(w*(dx_samecell*dy_samecell*(1.-dz_samecell)), i1, j1, k2, 1);
      
      ww += w*(dx_samecell*(1.-dy_samecell)*dz_samecell);
      density.set_ScalarField(w*(dx_samecell*(1.-dy_samecell)*dz_samecell), i1, j2, k1, 1);
      
      ww += w*(dx_samecell*(1.-dy_samecell)*(1.-dz_samecell));
      density.set_ScalarField(w*(dx_samecell*(1.-dy_samecell)*(1.-dz_samecell)), i1, j2, k2, 1);
      
      ww += w*((1.-dx_samecell)*dy_samecell*dz_samecell);
      density.set_ScalarField(w*((1.-dx_samecell)*dy_samecell*dz_samecell), i2, j1, k1, 1);
      
      ww += w*((1.-dx_samecell)*dy_samecell*(1.-dz_samecell));
      density.set_ScalarField(w*((1.-dx_samecell)*dy_samecell*(1.-dz_samecell)), i2, j1, k2, 1);
      
      ww += w*((1.-dx_samecell)*(1.-dy_samecell)*dz_samecell);
      density.set_ScalarField(w*((1.-dx_samecell)*(1.-dy_samecell)*dz_samecell), i2, j2, k1, 1);

      ww += w*((1.-dx_samecell)*(1.-dy_samecell)*(1.-dz_samecell));
      density.set_ScalarField(w*((1.-dx_samecell)*(1.-dy_samecell)*(1.-dz_samecell)), i2, j2, k2, 1);

    }
  }

  if (kernel_radius>0)
    density.GaussianConvolutionField(kernel_radius);
  
  double norm = 0;
  for (int i=0; i<nx; i++) 
    for (int j=0; j<ny; j++) 
      for (int k=0; k<nz; k++) 
	norm += density.ScalarField(i, j, k);

  norm = norm/nCells;

  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      for (int k=0; k<nz; k++) {
	double val = density.ScalarField(i, j, k)/norm-1;
	density.set_ScalarField(val, i, j, k);
      }
    }
  }

  return density;
}


// ============================================================================


ScalarField3D cosmobl::catalogue::Catalogue::density_field (const double cell_size, const Catalogue mask_catalogue, const int interpolation_type, const double kernel_radius, const bool useMass) const
{
  ScalarField3D mask_density = mask_catalogue.density_field(cell_size, interpolation_type, 0, 0);

  ScalarField3D density = density_field(cell_size, interpolation_type, 10, useMass);

  if (density.nx() != mask_density.nx() || density.ny() != mask_density.ny() || density.nz() != mask_density.nz())
  { ErrorMsg("Error in density_field, mask_catalogue is not correct"); } 
  
  for (int i=0; i<density.nx(); i++) 
    for (int j=0; j<density.ny(); j++) 
      for (int k=0; k<density.nz(); k++) {
	double val = density.ScalarField(i, j, k)+1;
	double mask_val = mask_density.ScalarField(i, j, k)+1;
	val = val*mask_val-1;
	density.set_ScalarField(val, i, j, k);
      }
    
  return density;
}
