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
using namespace chainmesh;


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

template void cosmobl::catalogue::Catalogue::replace_objects (vector<cosmobl::catalogue::RandomObject>);
template void cosmobl::catalogue::Catalogue::replace_objects (vector<cosmobl::catalogue::Mock>);
template void cosmobl::catalogue::Catalogue::replace_objects (vector<cosmobl::catalogue::Halo>);
template void cosmobl::catalogue::Catalogue::replace_objects (vector<cosmobl::catalogue::Galaxy>);
template void cosmobl::catalogue::Catalogue::replace_objects (vector<cosmobl::catalogue::Cluster>);
template void cosmobl::catalogue::Catalogue::replace_objects (vector<cosmobl::catalogue::Void>);

/// @endcond

// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const ObjType objType, const CoordType coordType, const vector<double> coord1, const vector<double> coord2, const vector<double> coord3, const vector<double> weight, const cosmology::Cosmology &cosm, const CoordUnits inputUnits)
{ 
  // check the vector dimensions
  if (!(coord1.size()==coord2.size() && coord2.size()==coord3.size()))
    ErrorCBL("Error in Catalogue::Catalogue() of Catalogue.cpp: coordinates with different dimensions!"); 
  
  // weight[]=1, if are not used
  vector<double> _weight = weight;
  if (_weight.size()==0) _weight.resize(coord1.size(), 1);
  
  
  // include the objects in the catalogue
  
  for (size_t i=0; i<coord1.size(); ++i) {

    if (coordType==_comovingCoordinates_) { // comoving coordinates (x, y, z)
      comovingCoordinates coord = {coord1[i], coord2[i], coord3[i]};
      m_object.push_back(move(Object::Create(objType, coord, _weight[i])));
    }
    else if (coordType==_observedCoordinates_) { // observed coordinates (R.A., Dec, redshift)
      observedCoordinates coord = {coord1[i], coord2[i], coord3[i]};
      m_object.push_back(move(Object::Create(objType, coord, inputUnits, cosm, _weight[i])));
    }
    else ErrorCBL("Error in Catalogue::Catalogue() of Catalogue.cpp: coordType is not valid!");

  }
  
}


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const ObjType objType, const CoordType coordType, const vector<string> file, const int col1, const int col2, const int col3, const int colWeight, const int colRegion, const double nSub, const double fact, const cosmology::Cosmology &cosm, const CoordUnits inputUnits) 
{ 
  // parameters for random numbers used in case nSub!=1
  
  default_random_engine gen;
  uniform_real_distribution<float> ran(0., 1.);
  
  
  // read the input catalogue files
  
  for (size_t dd=0; dd<file.size(); ++dd) {

    string line, file_in = file[dd];
    coutCBL << "I'm reading the catalogue: " << file_in << endl;
    ifstream finr(file_in.c_str()); checkIO(finr, file_in);

    double Weight, Value; long Region;

    while (getline(finr, line)) { // read the lines
      
      if (ran(gen)<nSub) { // extract a subsample
	
	stringstream ss(line);
	vector<double> value; while (ss>>Value) value.emplace_back(Value);
	checkDim(value, col1, "value", false); checkDim(value, col2, "value", false); 

	if (coordType==_comovingCoordinates_) { // comoving coordinates (x, y, z)
	  checkDim(value, col3, "value", false);
	  comovingCoordinates coord;
	  coord.xx = value[col1-1]*fact;
	  coord.yy = value[col2-1]*fact;
	  coord.zz = value[col3-1]*fact;
	  Weight = (colWeight!=-1) ? value[colWeight-1] : 1.;
	  Region = (long)value[colRegion-1];
	  m_object.push_back(move(Object::Create(objType, coord, Weight, Region)));
	}

	else if (coordType==_observedCoordinates_) { // observed coordinates (R.A., Dec (redshift))
	  observedCoordinates coord;
	  coord.ra = value[col1-1]*fact;
	  coord.dec = value[col2-1]*fact;
	  coord.redshift = ((int)value.size()>=col3) ? value[col3-1] : 1.;
	  Weight = (colWeight!=-1) ? value[colWeight-1] : 1.;
	  Region = (long)value[colRegion-1];
	  m_object.push_back(move(Object::Create(objType, coord, inputUnits, cosm, Weight, Region)));
	}
	
	else ErrorCBL("Error in Catalogue::Catalogue() of Catalogue.cpp: coordType is not valid!");	
      }
    }
    
    finr.clear(); finr.close();
  }
}


// ============================================================================


vector<long> cosmobl::catalogue::Catalogue::region () const
{ 
  vector<long> vv(m_object.size());
  
  for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->region();

  return vv;
}


// ============================================================================


vector<string> cosmobl::catalogue::Catalogue::field () const
{ 
  vector<string> vv(m_object.size());
  
  for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->field();

  return vv;
}

// ============================================================================

    
vector<double> cosmobl::catalogue::Catalogue::var (Var var_name) const
{ 
  vector<double> vv(m_object.size(), 0.);
  
  switch (var_name) {

  case Var::_X_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->xx();
    break;

  case Var::_Y_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->yy();
    break;

  case Var::_Z_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->zz();
    break;

  case Var::_RA_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->ra();
    break;

  case Var::_Dec_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->dec();
    break;

  case Var::_Redshift_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->redshift();
    break;

  case Var::_Dc_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->dc();
    break;

  case Var::_Weight_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->weight();
    break;

  case Var::_Mass_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->mass();
    break;

  case Var::_Magnitude_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->magnitude();
    break;

  case Var::_Richness_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->richness();
    break;

  case Var::_Vx_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->vx();
    break;
  
  case Var::_Vy_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->vy();
    break;
  
  case Var::_Vz_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->vz();
    break;

  case Var::_Region_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->region();
    break; 
  
  case Var::_Generic_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->generic();
    break;

  case Var::_Radius_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->radius();
    break;

  case Var::_X_displacement_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->radius();
    break;
    
  case Var::_Y_displacement_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->radius();
    break;
    
  case Var::_Z_displacement_:
    for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->radius();
    break;

  default:
    ErrorCBL("Error in cosmobl::catalogue::Catalogue::var of Catalogue.cpp: no such a variable in the list!");
  }
  
  return vv;
}


// ============================================================================


void cosmobl::catalogue::Catalogue::set_region (const vector<long> region)
{
  for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_region(region[i]);
}


// ============================================================================


void cosmobl::catalogue::Catalogue::set_field (const vector<string> field)
{
  for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_field(field[i]);
}


// ============================================================================


void cosmobl::catalogue::Catalogue::set_var (const Var var_name, const vector<double> var)
{
  if (m_object.size()!=var.size()) ErrorCBL("Error in cosmobl::catalogue::Catalogue::set_var of Catalogue.cpp: different sizes!");
  
  switch (var_name) {

  case Var::_X_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_xx(var[i]);
    break;

  case Var::_Y_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_yy(var[i]);
    break;

  case Var::_Z_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_zz(var[i]);
    break;

  case Var::_RA_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_ra(var[i]);
    break;

  case Var::_Dec_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_dec(var[i]);
    break;

  case Var::_Redshift_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_redshift(var[i]);
    break;

  case Var::_Dc_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_dc(var[i]);
    break;

  case Var::_Weight_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_weight(var[i]);
    break;

  case Var::_Mass_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_mass(var[i]);
    break;

  case Var::_Magnitude_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_magnitude(var[i]);
    break;

  case Var::_Richness_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_richness(var[i]);
    break;

  case Var::_Vx_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_vx(var[i]);
    break;
  
  case Var::_Vy_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_vy(var[i]);
    break;
  
  case Var::_Vz_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_vz(var[i]);
    break;

  case Var::_Region_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_region(var[i]);
    break;
    
  case Var::_Generic_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_generic(var[i]);
    break;

  case Var::_Radius_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_radius(var[i]);
    break;

  case Var::_X_displacement_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_x_displacement();
    break;
    
  case Var::_Y_displacement_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_y_displacement();
    break;
    
  case Var::_Z_displacement_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_z_displacement();
    break;

  default:
    ErrorCBL("Error in cosmobl::catalogue::Catalogue::set_var of Catalogue.cpp: no such a variable in the list!");
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


void cosmobl::catalogue::Catalogue::stats_var (const vector<Var> var_name, vector<vector<double>> &stats) const
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


void cosmobl::catalogue::Catalogue::computeComovingCoordinates (const cosmology::Cosmology &cosm, const CoordUnits inputUnits)
{
  double red, xx, yy, zz;

  
  // ----- unit conversion -----
  
  vector<double> RA(nObjects()), DEC(nObjects());

  for (size_t i=0; i<nObjects(); ++i) {
    RA[i] = (inputUnits==_radians_) ? ra(i) : radians(ra(i), inputUnits);
    DEC[i] = (inputUnits==_radians_) ? dec(i) : radians(dec(i), inputUnits);
  }

  
  // ----- compute comoving coordinates -----
  
  for (size_t i=0; i<nObjects(); ++i) {

    red = redshift(i);
    m_object[i]->set_dc(cosm.D_C(red));
    
    cartesian_coord(RA[i], DEC[i], dc(i), xx, yy, zz);
    
    m_object[i]->set_xx(xx); 
    m_object[i]->set_yy(yy); 
    m_object[i]->set_zz(zz);
  }
}


// ============================================================================


void cosmobl::catalogue::Catalogue::computePolarCoordinates (const CoordUnits outputUnits)
{
  double ra, dec, dc;

  
  // ----- compute polar coordinates in radians -----
  
  for (size_t i=0; i<nObjects(); ++i) {
    polar_coord(xx(i), yy(i), zz(i), ra, dec, dc);
    m_object[i]->set_ra(ra); 
    m_object[i]->set_dec(dec); 
    m_object[i]->set_dc(dc);
  }

  
  // ----- unit conversion -----
  
  if (outputUnits!=_radians_) {

    if (outputUnits==_degrees_)
      for (size_t i=0; i<nObjects(); ++i) {
	m_object[i]->set_ra(degrees(ra));
	m_object[i]->set_dec(degrees(dec));
      }

    else if (outputUnits==_arcseconds_)
      for (size_t i=0; i<nObjects(); ++i) {
	m_object[i]->set_ra(arcseconds(ra));
	m_object[i]->set_dec(arcseconds(dec)); 		    
      }

    else if (outputUnits==_arcminutes_)
      for (size_t i=0; i<nObjects(); ++i) {
	m_object[i]->set_ra(arcminutes(ra));
	m_object[i]->set_dec(arcminutes(dec)); 		    
      }

    else ErrorCBL("Error in cosmobl::catalogue::Catalogue::computePolarCoordinates: outputUnits type not allowed!");
  }
  
}

// ============================================================================


void cosmobl::catalogue::Catalogue::computePolarCoordinates (const cosmology::Cosmology &cosm, const double z1, const double z2, const CoordUnits outputUnits)
{
  double ra, dec, dc;

  
  // ----- compute polar coordinates in radians -----

  for (size_t i=0; i<nObjects(); ++i) {
    polar_coord(xx(i), yy(i), zz(i), ra, dec, dc);
    m_object[i]->set_ra(ra); 
    m_object[i]->set_dec(dec); 
    m_object[i]->set_dc(dc);
    m_object[i]->set_redshift(cosm.Redshift(dc, z1, z2));
  }

  
  // ----- unit conversion -----

  if (outputUnits!=_radians_) {

    if (outputUnits==_degrees_)
      for (size_t i=0; i<nObjects(); ++i) {
	m_object[i]->set_ra(degrees(ra));
	m_object[i]->set_dec(degrees(dec)); 		    
      }

    else if (outputUnits==_arcseconds_)
      for (size_t i=0; i<nObjects(); ++i) {
	m_object[i]->set_ra(arcseconds(ra));
	m_object[i]->set_dec(arcseconds(dec)); 		    
      }

    else if (outputUnits==_arcminutes_)
      for (size_t i=0; i<nObjects(); ++i) {
	m_object[i]->set_ra(arcminutes(ra));
	m_object[i]->set_dec(arcminutes(dec)); 		    
      }

    else ErrorCBL("Error in cosmobl::catalogue::Catalogue::computePolarCoordinates: outputUnits type not allowed!");
  }
  
}

// ============================================================================


void cosmobl::catalogue::Catalogue::normalizeComovingCoordinates () 
{
  for (size_t i=0; i<nObjects(); ++i) { 
    m_object[i]->set_xx(xx(i)/dc(i)); 
    m_object[i]->set_yy(yy(i)/dc(i)); 
    m_object[i]->set_zz(zz(i)/dc(i));
  }
}


// ============================================================================


void cosmobl::catalogue::Catalogue::restoreComovingCoordinates ()
{
  for (size_t i=0; i<nObjects(); ++i) {
    m_object[i]->set_xx(xx(i)*dc(i)); 
    m_object[i]->set_yy(yy(i)*dc(i)); 
    m_object[i]->set_zz(zz(i)*dc(i));
  }
}


// ============================================================================


void cosmobl::catalogue::Catalogue::Order (const vector<int> vv) 
{
  int nObj = m_object.size();

  if (int(vv.size())!=nObj) ErrorCBL("Error in Catalogue::Order()!");
 
  vector<shared_ptr<Object>> obj(nObj);
  
  m_index.resize(nObj);
  
  for (size_t i=0; i<vv.size(); i++) {
    m_index[i] = vv[i];
    obj[i] = m_object[vv[i]];
  }

  m_object = obj;
}


// ============================================================================


void cosmobl::catalogue::Catalogue::Order () 
{ 
  size_t nObj = m_object.size();
  
  vector<shared_ptr<Object>> obj(nObj);
  
  if (m_index.size() != nObj) 
    ErrorCBL("Error in Catalogue::Order() of Catalogue.cpp, order not found!");
  
  obj = m_object;
  
  for (size_t i=0; i<nObj; i++) 
    m_object[i] = obj[m_index[i]];
}


// ============================================================================


double cosmobl::catalogue::Catalogue::weightedN () const
{
  double nn = 0.;
  for (size_t i=0; i<m_object.size(); i++)
    nn += m_object[i]->weight();
  return nn;
}


// ============================================================================


void cosmobl::catalogue::Catalogue::write_comoving_coordinates (const string outputFile) const
{
  if (m_object.size()==0) ErrorCBL("Error in cosmobl::catalogue::Catalogue::write_comoving_coordinates: m_object.size()=0!");

  ofstream fout(outputFile.c_str()); checkIO(fout, outputFile);
  
  for (size_t i=0; i<nObjects(); ++i) 
    fout << xx(i) << "   " << yy(i) << "   " << zz(i) << endl;

  coutCBL << "I wrote the file: " << outputFile << endl;
  fout.clear(); fout.close();
}


// ============================================================================


void cosmobl::catalogue::Catalogue::write_obs_coordinates (const string outputFile) const 
{
  if (m_object.size()==0) ErrorCBL("Error in cosmobl::catalogue::Catalogue::write_obs_coordinates: m_object.size()=0!");
  
  ofstream fout(outputFile.c_str()); checkIO(fout, outputFile);
  
  if (!isSet(ra(0)) || !isSet(dec(0)) || !isSet(redshift(0)))
    ErrorCBL("Error in cosmobl::catalogue::Catalogue::write_obs_coords of Catalogue.cpp! Polar coordinates are not set!");

  for (size_t i=0; i<nObjects(); ++i) 
    fout << ra(i) << "   " << dec(i) << "   " << redshift(i) << endl;
  
  coutCBL << "I wrote the file: " << outputFile << endl;
  fout.clear(); fout.close();
}


// ============================================================================


void cosmobl::catalogue::Catalogue::write_data (const string outputFile) const 
{
  ofstream fout(outputFile.c_str()); checkIO(fout, outputFile);

  if (!isSet(ra(0)) || !isSet(dec(0)) || !isSet(redshift(0)) || !isSet(dc(0)))
    ErrorCBL("Error in cosmobl::catalogue::Catalogue::write_coordinates of Catalogue.cpp! Polar coordinates are not set!");

  if (!isSet(region(0)))
    for (size_t i=0; i<nObjects(); ++i) 
      fout << xx(i) << "   " << yy(i) << "   " << zz(i) << "   " << ra(i) << "   " << dec(i) << "   " << redshift(i) << "   " << dc(i) << endl;
  else
    for (size_t i=0; i<nObjects(); ++i) 
      fout << xx(i) << "   " << yy(i) << "   " << zz(i) << "   " << ra(i) << "   " << dec(i) << "   " << redshift(i) << "   " << dc(i) << "   " << region(i) <<endl;
      
  coutCBL << "I wrote the file: " << outputFile << endl;
  fout.clear(); fout.close();
}


// ============================================================================


Catalogue cosmobl::catalogue::Catalogue::cut (const Var var_name, const double down, const double up, const bool excl)
{
  vector<shared_ptr<Object> > objects;
  vector<double> vvar = var(var_name);
  vector<int> w(vvar.size());

  for (size_t i=0; i<m_object.size(); i++) {
    w[i] = (excl) ? 1 : 0;
    if (vvar[i] >= down && vvar[i] < up)
      w[i] = (excl) ? 0 : 1;
  }
  
  for (size_t i=0; i<m_object.size(); i++)
    if (w[i]==1)
      objects.push_back(m_object[i]);
 
  return Catalogue {objects};
}


// ============================================================================


double cosmobl::catalogue::Catalogue::distance (const int i, shared_ptr<Object> obj) const
{
  return sqrt((m_object[i]->xx()-obj->xx())*(m_object[i]->xx()-obj->xx())+
	      (m_object[i]->yy()-obj->yy())*(m_object[i]->yy()-obj->yy())+
	      (m_object[i]->zz()-obj->zz())*(m_object[i]->zz()-obj->zz()));
}


// ============================================================================


double cosmobl::catalogue::Catalogue::angsep_xyz (const int i, shared_ptr<Object> obj) const
{ 
  return 2.*asin(0.5*sqrt((m_object[i]->xx()-obj->xx())*(m_object[i]->xx()-obj->xx())+
			  (m_object[i]->yy()-obj->yy())*(m_object[i]->yy()-obj->yy())+
			  (m_object[i]->zz()-obj->zz())*(m_object[i]->zz()-obj->zz())));
}
    

// ============================================================================


shared_ptr<Catalogue> cosmobl::catalogue::Catalogue::smooth (const double gridsize, const vector<Var> vars, const int SUB)
{
  (void)vars;
  
  shared_ptr<Catalogue> cat {new Catalogue(*this)};
  
  if (gridsize<1.e-30) return cat;
  
  double rMAX = 0.;
  
  vector<shared_ptr<Object>> sample;

  
  // ----------------------------------------------------------------------------
  // ----- subdivide the catalogue in SUB^3 sub-catalogues, to avoid ------------
  // ----- memory problems possibly caused by too large chain-mesh vectors ------
  // ----------------------------------------------------------------------------

  coutCBL <<"Please wait, I'm subdividing the catalogue in "<<pow(SUB, 3)<<" sub-catalogues..."<<endl;
 
  int nx = SUB, ny = SUB, nz = SUB;
  
  double Cell_X = (Max(Var::_X_)-Min(Var::_X_))/nx;
  double Cell_Y = (Max(Var::_Y_)-Min(Var::_Y_))/ny;
  double Cell_Z = (Max(Var::_Z_)-Min(Var::_Z_))/nz;

  for (size_t i=0; i<cat->nObjects(); i++) {
    int i1 = min(int((cat->xx(i)-Min(Var::_X_))/Cell_X), nx-1);
    int j1 = min(int((cat->yy(i)-Min(Var::_Y_))/Cell_Y), ny-1);
    int z1 = min(int((cat->zz(i)-Min(Var::_Z_))/Cell_Z), nz-1);
    int index = z1+nz*(j1+ny*i1);
    cat->catalogue_object(i)->set_region(index);
  }
  
  size_t nRegions = cat->nRegions();
  
  vector<Catalogue> subSamples(nRegions);
  
  for (size_t i=0; i<nRegions; ++i) {
    double start = (double)cat->region_list()[i];
    double stop = start+1;
    subSamples[i] = cut(Var::_Region_, start, stop);
  }

  
  // -----------------------------------------
  // ----- smooth all the sub-catalogues -----
  // -----------------------------------------

  for (size_t rr=0; rr<nRegions; ++rr) {
    
    vector<double> _xx = subSamples[rr].var(Var::_X_), _yy = subSamples[rr].var(Var::_Y_), _zz = subSamples[rr].var(Var::_Z_);
    
    ChainMesh3D ll(gridsize, _xx, _yy, _zz, rMAX, (long)-1.e5, (long)1.e5);
   
    
    for (int i=0; i<ll.nCell(); i++) {
      vector<long> list = ll.get_list(i);
     
      int nObj = list.size();
      if (nObj>0) {
	
	double XX = 0., YY = 0., ZZ = 0., RA = 0., DEC = 0., REDSHIFT = 0., WEIGHT = 0.;

	for (size_t j=0; j<list.size(); j++) {
	  XX += subSamples[rr].xx(list[j]);
	  YY += subSamples[rr].yy(list[j]);
	  ZZ += subSamples[rr].zz(list[j]);
	  RA += subSamples[rr].ra(list[j]);
	  DEC += subSamples[rr].dec(list[j]);
	  REDSHIFT += subSamples[rr].redshift(list[j]);
	  WEIGHT += subSamples[rr].weight(list[j]);
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

  for (size_t i=0; i<m_object.size(); i++)

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

  for (size_t i=0; i<m_object.size(); i++)
    if (vvar[i] >= down && vvar[i] < up)
      nObjw += weight(i);
    
  nObjw = (excl) ? weightedN()-nObjw : nObjw;

  return nObjw;
}


// ============================================================================


data::ScalarField3D cosmobl::catalogue::Catalogue::density_field (const double cell_size, const int interpolation_type, const double kernel_radius, const bool useMass) const
{
  data::ScalarField3D density(cell_size, Min(Var::_X_), Max(Var::_X_), Min(Var::_Y_), Max(Var::_Y_), Min(Var::_Z_), Max(Var::_Z_));

  double deltaX = density.deltaX();
  double deltaY = density.deltaY();
  double deltaZ = density.deltaZ();
  int nx = density.nx();
  int ny = density.ny();
  int nz = density.nz();

  long int nCells = density.nCells();
  double nCells_inv = 1./nCells;

  for (size_t i=0; i<nObjects(); ++i) {
    int i1 = min(int((xx(i)-density.MinX())/deltaX),nx-1);
    int j1 = min(int((yy(i)-density.MinY())/deltaY),ny-1);
    int k1 = min(int((zz(i)-density.MinZ())/deltaZ),nz-1);

    double w = (useMass) ? mass(i) : 1.;

    if (interpolation_type==0) {
      density.set_ScalarField(w, i1, j1, k1,1);
    }
    else if (interpolation_type==1) {
      double dx_samecell, dy_samecell, dz_samecell;
      int dx_index, dy_index, dz_index;

      dx_samecell = (xx(i)-(i1*deltaX+density.MinX()))/deltaX;
      dy_samecell = (yy(i)-(j1*deltaY+density.MinY()))/deltaY;
      dz_samecell = (zz(i)-(k1*deltaZ+density.MinZ()))/deltaZ;

      dx_index = (dx_samecell <0.5) ? -1 : 1 ;
      dy_index = (dy_samecell <0.5) ? -1 : 1 ;
      dz_index = (dz_samecell <0.5) ? -1 : 1 ;
      
      dx_samecell = (dx_samecell<0.5) ? (xx(i)+0.5*deltaX-(i1*deltaX+density.MinX()))/deltaX : 1-(xx(i)-0.5*deltaX-(i1*deltaX+density.MinX()))/deltaX;
      dy_samecell = (dy_samecell<0.5) ? (yy(i)+0.5*deltaY-(j1*deltaY+density.MinY()))/deltaY : 1-(yy(i)-0.5*deltaY-(j1*deltaY+density.MinY()))/deltaY;
      dz_samecell = (dz_samecell<0.5) ? (zz(i)+0.5*deltaZ-(k1*deltaZ+density.MinZ()))/deltaZ : 1-(zz(i)-0.5*deltaZ-(k1*deltaZ+density.MinZ()))/deltaZ;
      
      if (dx_samecell<0 || dy_samecell<0 || dz_samecell<0)
	coutCBL << dx_index << " " << dx_samecell << " " << dy_index << " " << dy_samecell << " " << dz_index << " " << dz_samecell << endl;

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

  double wMean = 0;
  for (int i=0; i<nx; i++) 
    for (int j=0; j<ny; j++) 
      for (int k=0; k<nz; k++) 
	wMean += density.ScalarField(i, j, k);

  wMean *= nCells_inv;

  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      for (int k=0; k<nz; k++) {
	double val = density.ScalarField(i, j, k)/wMean-1;
	density.set_ScalarField(val, i, j, k);
      }
    }
  }

  return density;
}


// ============================================================================


data::ScalarField3D cosmobl::catalogue::Catalogue::density_field (const double cell_size, const Catalogue mask_catalogue, const int interpolation_type, const double kernel_radius, const bool useMass) const
{
  (void)kernel_radius;
  
  data::ScalarField3D mask_density = mask_catalogue.density_field(cell_size, interpolation_type, 0, 0);

  data::ScalarField3D density = density_field(cell_size, interpolation_type, 10, useMass);

  if (density.nx() != mask_density.nx() || density.ny() != mask_density.ny() || density.nz() != mask_density.nz())
    ErrorCBL("Error in density_field, mask_catalogue is not correct"); 
  
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


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const Catalogue input_catalogue, const Catalogue target_catalogue, const Var var_name1, const int nbin1, const int seed) 
{
  vector<double> fvar_input(nbin1, 0);
  vector<double> fvar_target(nbin1, 0);

  vector<double> input_var = input_catalogue.var(var_name1);
  vector<double> target_var = target_catalogue.var(var_name1);

  double Vmin = target_catalogue.Min(var_name1);
  double Vmax = target_catalogue.Max(var_name1);

  double binSize_inv = pow((Vmax-Vmin)/nbin1, -1);

  for (size_t i=0; i<input_var.size(); i++)
    if (input_var[i] < Vmax && Vmin < input_var[i]) {
      int occ = max(0, min(int((input_var[i]-Vmin)*binSize_inv), nbin1));
      fvar_input[occ] ++;
    }

  for (size_t i=0; i<target_var.size(); i++)
    if (target_var[i] < Vmax && Vmin < target_var[i]) {
      int occ = max(0, min(int((target_var[i]-Vmin)*binSize_inv), nbin1));
      fvar_target[occ] ++;
    }

  random::UniformRandomNumbers ran(0., 1., seed);

  for (size_t i=0; i<input_var.size(); i++) 
    if (input_var[i] < Vmax && Vmin < input_var[i]) {
      int occ = max(0, min(int((input_var[i]-Vmin)*binSize_inv), nbin1));
      if (ran() < fvar_target[occ]/fvar_input[occ])
	m_object.push_back(shared_ptr<Object>(input_catalogue.catalogue_object(i)));
    }
}


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const Catalogue input_catalogue, const Catalogue target_catalogue, const Var var_name1, const int nbin1, const Var var_name2, const int nbin2, const int seed) 
{
  vector< vector<double> > fvars_input(nbin1, vector<double>(nbin2, 0));
  vector< vector<double> > fvars_target(nbin1, vector<double>(nbin2, 0));

  vector<double> input_var1 = input_catalogue.var(var_name1);
  vector<double> target_var1 = target_catalogue.var(var_name1);
  vector<double> input_var2 = input_catalogue.var(var_name2);
  vector<double> target_var2 = target_catalogue.var(var_name2);

  double V1min = target_catalogue.Min(var_name1);
  double V1max = target_catalogue.Max(var_name1);

  double V2min = target_catalogue.Min(var_name2);
  double V2max = target_catalogue.Max(var_name2);

  double binSize1_inv = pow((V1max-V1min)/nbin1,-1);
  double binSize2_inv = pow((V2max-V2min)/nbin2,-1);

  for (size_t i=0; i<input_var1.size(); i++)
    if ( (input_var1[i] < V1max && V1min < input_var1[i]) && (input_var2[i] < V2max && V2min < input_var2[i])) {
      int occ1 = max(0, min(int((input_var1[i]-V1min)*binSize1_inv), nbin1));
      int occ2 = max(0, min(int((input_var2[i]-V2min)*binSize2_inv), nbin2));
      fvars_input[occ1][occ2] ++;
    }

  for (size_t i=0; i<target_var1.size(); i++)
    if ( (target_var1[i] < V1max && V1min < target_var2[i]) && (target_var2[i] < V2max && V2min < target_var2[i])) {
      int occ1 = max(0, min(int((target_var1[i]-V1min)*binSize1_inv), nbin1));
      int occ2 = max(0, min(int((target_var2[i]-V2min)*binSize2_inv), nbin2));
      fvars_target[occ1][occ2] ++;
    }


  random::UniformRandomNumbers ran(0., 1., seed);

  for (size_t i=0; i<input_var1.size(); i++)
    if ((input_var1[i] < V1max && V1min < input_var1[i]) && (input_var2[i] < V2max && V2min < input_var2[i])) {
      int occ1 = max(0, min(int((input_var1[i]-V1min)*binSize1_inv), nbin1));
      int occ2 = max(0, min(int((input_var2[i]-V2min)*binSize2_inv), nbin2));
      if (ran() < fvars_target[occ1][occ2]/fvars_input[occ1][occ2])
	m_object.push_back(shared_ptr<Object>(input_catalogue.catalogue_object(i)));
    }
  
} 
