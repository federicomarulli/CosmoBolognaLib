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
 *  @author federico.marulli3@unibo.it
 */

#include "Field3D.h"
#include "Catalogue.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace chainmesh;


// ============================================================================

/// @cond template

template cbl::catalogue::Catalogue::Catalogue (vector<cbl::catalogue::RandomObject>);
template cbl::catalogue::Catalogue::Catalogue (vector<cbl::catalogue::Mock>);
template cbl::catalogue::Catalogue::Catalogue (vector<cbl::catalogue::Halo>);
template cbl::catalogue::Catalogue::Catalogue (vector<cbl::catalogue::Galaxy>);
template cbl::catalogue::Catalogue::Catalogue (vector<cbl::catalogue::Cluster>);
template cbl::catalogue::Catalogue::Catalogue (vector<cbl::catalogue::Void>);
template cbl::catalogue::Catalogue::Catalogue (vector<cbl::catalogue::HostHalo>);

template void cbl::catalogue::Catalogue::add_object (cbl::catalogue::RandomObject);
template void cbl::catalogue::Catalogue::add_object (cbl::catalogue::Mock);
template void cbl::catalogue::Catalogue::add_object (cbl::catalogue::Halo);
template void cbl::catalogue::Catalogue::add_object (cbl::catalogue::Galaxy);
template void cbl::catalogue::Catalogue::add_object (cbl::catalogue::Cluster);
template void cbl::catalogue::Catalogue::add_object (cbl::catalogue::Void);
template void cbl::catalogue::Catalogue::add_object (cbl::catalogue::HostHalo);

template void cbl::catalogue::Catalogue::add_objects (vector<cbl::catalogue::RandomObject>);
template void cbl::catalogue::Catalogue::add_objects (vector<cbl::catalogue::Mock>);
template void cbl::catalogue::Catalogue::add_objects (vector<cbl::catalogue::Halo>);
template void cbl::catalogue::Catalogue::add_objects (vector<cbl::catalogue::Galaxy>);
template void cbl::catalogue::Catalogue::add_objects (vector<cbl::catalogue::Cluster>);
template void cbl::catalogue::Catalogue::add_objects (vector<cbl::catalogue::Void>);
template void cbl::catalogue::Catalogue::add_objects (vector<cbl::catalogue::HostHalo>);

template void cbl::catalogue::Catalogue::replace_objects (vector<cbl::catalogue::RandomObject>);
template void cbl::catalogue::Catalogue::replace_objects (vector<cbl::catalogue::Mock>);
template void cbl::catalogue::Catalogue::replace_objects (vector<cbl::catalogue::Halo>);
template void cbl::catalogue::Catalogue::replace_objects (vector<cbl::catalogue::Galaxy>);
template void cbl::catalogue::Catalogue::replace_objects (vector<cbl::catalogue::Cluster>);
template void cbl::catalogue::Catalogue::replace_objects (vector<cbl::catalogue::Void>);
template void cbl::catalogue::Catalogue::replace_objects (vector<cbl::catalogue::HostHalo>);

/// @endcond


// ============================================================================


cbl::catalogue::Catalogue::Catalogue (const ObjectType objectType, const CoordinateType coordinateType, const std::vector<double> coord1, const std::vector<double> coord2, const std::vector<double> coord3, const std::vector<double> weight, const cosmology::Cosmology &cosm, const CoordinateUnits inputUnits)
{ 
  // check the vector dimensions
  if (!(coord1.size()==coord2.size() && coord2.size()==coord3.size()))
    ErrorCBL("coordinates with different dimensions!", "Catalogue", "Catalogue.cpp"); 
  
  // weight[]=1, if are not used
  vector<double> _weight = weight;
  if (_weight.size()==0) _weight.resize(coord1.size(), 1);
  
  
  // include the objects in the catalogue
  
  for (size_t i=0; i<coord1.size(); ++i) {

    if (coordinateType==cbl::CoordinateType::_comoving_) { // comoving coordinates (x, y, z)
      comovingCoordinates coord = {coord1[i], coord2[i], coord3[i]};
      m_object.push_back(move(Object::Create(objectType, coord, _weight[i])));
    }
    else if (coordinateType==cbl::CoordinateType::_observed_) { // observed coordinates (R.A., Dec, redshift)
      observedCoordinates coord = {coord1[i], coord2[i], coord3[i]};
      m_object.push_back(move(Object::Create(objectType, coord, inputUnits, cosm, _weight[i])));
    }
    else ErrorCBL("CoordinateType is not valid!", "Catalogue", "Catalogue.cpp");

  }
  
}


// ============================================================================


cbl::catalogue::Catalogue::Catalogue (const ObjectType objectType, const CoordinateType coordinateType, const std::vector<std::string> file, const int col1, const int col2, const int col3, const int colWeight, const int colRegion, const double nSub, const double fact, const cosmology::Cosmology &cosm, const CoordinateUnits inputUnits, const CharEncode charEncode, const std::string comment, const int seed) 
{ 
  // parameters for random numbers used in case nSub!=1
  random::UniformRandomNumbers ran(0., 1., seed);
  
  // read the input catalogue files
  
  for (size_t dd=0; dd<file.size(); ++dd) {

    string line, file_in = file[dd];
         
    if (charEncode==CharEncode::_ascii_) {
	  
      coutCBL << "I'm reading the catalogue: " << file_in << endl;
      ifstream finr(file_in.c_str()); checkIO(finr, file_in);
      
      double Weight, Value; long Region;

      long maxRegion = -100000000; // check!
      
      while (getline(finr, line)) { // read the lines

	if (line.find(comment)==string::npos) { // skip a comment

	  if (ran()<nSub) { // extract a subsample
	  
	    stringstream ss(line);
	    vector<double> value; while (ss>>Value) value.emplace_back(Value);
	    checkDim(value, col1, "value", false); checkDim(value, col2, "value", false); 

	    if (coordinateType==cbl::CoordinateType::_comoving_) { // comoving coordinates (x, y, z)
	      checkDim(value, col3, "value", false);
	      comovingCoordinates coord;
	      coord.xx = value[col1-1]*fact;
	      coord.yy = value[col2-1]*fact;
	      coord.zz = value[col3-1]*fact;
	      Weight = (colWeight!=-1 && colWeight-1<(int)value.size()) ? value[colWeight-1] : 1.;
	      Region = (colRegion!=-1 && colRegion-1<(int)value.size()) ? (long)value[colRegion-1] : 0;
	      m_object.push_back(move(Object::Create(objectType, coord, Weight, Region)));
	    }

	    else if (coordinateType==cbl::CoordinateType::_observed_) { // observed coordinates (R.A., Dec (redshift))
	      observedCoordinates coord;
	      coord.ra = value[col1-1]*fact;
	      coord.dec = value[col2-1]*fact;
	      coord.redshift = ((int)value.size()>=col3) ? value[col3-1] : 1.;
	      Weight = (colWeight!=-1 && colWeight-1<(int)value.size()) ? value[colWeight-1] : 1.;
	      Region = (colRegion!=-1 && colRegion-1<(int)value.size()) ? (long)value[colRegion-1] : 0;
	      m_object.push_back(move(Object::Create(objectType, coord, inputUnits, cosm, Weight, Region)));
	    }
	  
	    else ErrorCBL("CoordinateType is not valid!", "Catalogue", "Catalogue.cpp");	

	    maxRegion = max(maxRegion, Region);
	  }
	}
      }

      if (colRegion>0) {
	m_nRegions = maxRegion+1;
	WarningMsgCBL("The total number of region is "+conv(m_nRegions, par::fINT)+", deducted from the maximum value of the input regions; if needed, this number can be changed with the function set_region_number()", "Catalogue", "Catalogue.cpp");
      }
      
      finr.clear(); finr.close();
    }
          
    else if (charEncode==CharEncode::_binary_) {
            
      // read the input catalogue files
	
      coutCBL << "I'm reading the catalogue: " << file_in << endl;

      short num_bin;
      float val;
      
      ifstream finr(file_in.c_str(), ios::in|ios::binary|ios::ate); checkIO(finr, file_in);	

      comovingCoordinates coord;

      if (finr.is_open())
	finr.seekg(0, ios::beg);
      
      finr.read((char*)(&num_bin), 2);

      int n_blocks = num_bin;
      
      for (int i=1; i<=n_blocks; ++i) {
        finr.read((char*)(&num_bin), 4);
        int n_objs = num_bin;

        for (int j=1; j<=n_objs; ++j) {
          finr.read((char*)(&val), 4);
          coord.xx = (val)*fact;
          finr.read((char*)(&val), 4);
          coord.yy = (val)*fact;
          finr.read((char*)(&val), 4);
          coord.zz = (val)*fact;
	  // Weight = (colWeight!=-1 && colWeight-1<value.size()) ? value[colWeight-1] : 1.;
	  // Region = (colRegion!=-1 && colRegion-1<value.size()) ? (long)value[colRegion-1] : 0;
          if (ran()<nSub) m_object.push_back(move(Object::Create(objectType, coord)));
	  // if (ran()<nSub) m_object.push_back(move(Object::Create(objectType, coord, Weight, Region)));
        }
	
        finr.read((char*)(&num_bin), 4);

	int n_objs2 = num_bin;
	
        if (n_objs2!=n_objs) ErrorCBL("wrong reading of input binary file", "Catalogue", "Catalogue.cpp");
      }
      
      finr.clear(); finr.close();
    }
        
    else ErrorCBL("charEncode is not valid!", "Catalogue", "Catalogue.cpp");
  }
}


// ============================================================================


cbl::catalogue::Catalogue::Catalogue (const ObjectType objectType, const CoordinateType coordinateType, const std::vector<Var> attribute, const std::vector<int> column, const std::vector<std::string> file, const int comments, const double nSub, const double fact, const cosmology::Cosmology &cosm, const CoordinateUnits inputUnits, const int seed) 
{
  // preliminary check on vector sizes
  size_t nvar;
  if (attribute.size()==column.size()) nvar = attribute.size();
  else ErrorCBL("Column vector and attribute vector must have equal size!", "Catalogue", "Catalogue.cpp");

  if (!(std::is_sorted(column.begin(), column.end()))) ErrorCBL("Column vector must be sort in ascending order!", "Catalogue", "Catalogue.cpp");
  
  // std::unordered_map to map columns and attributes
  unordered_map<int, Var> varMap;
  for (size_t ii=0; ii<nvar; ii++)
    varMap.insert({column[ii], attribute[ii]});
  
  // parameters for random numbers used in case nSub!=1
  random::UniformRandomNumbers ran(0., 1., seed);
  
  // read the input catalogue files
 
  for (size_t dd=0; dd<file.size(); ++dd) {

    string line, file_in = file[dd];
	  
    coutCBL << "I'm reading the catalogue: " << file_in << endl;
    ifstream finr(file_in.c_str()); checkIO(finr, file_in);

    // prepare default coordinates
    comovingCoordinates defaultComovingCoord = { par::defaultDouble, par::defaultDouble, par::defaultDouble};
    observedCoordinates defaultObservedCoord = { par::defaultDouble, -1., 0.1};
  
    for (int cc=0; cc<comments; cc++) getline(finr, line); // ignore commented lines at the beginning of file
      
    while (getline(finr, line)) { // read the lines

      if (ran()<nSub) { // extract a subsample
	  
	if (coordinateType==cbl::CoordinateType::_comoving_) 
	  m_object.push_back(move(Object::Create(objectType, defaultComovingCoord, 1.)));
	  
	else if (coordinateType==cbl::CoordinateType::_observed_)
	  m_object.push_back(move(Object::Create(objectType, defaultObservedCoord, inputUnits, cosm, 1.)));

	else ErrorCBL("CoordinateType is not valid!", "Catalogue", "Catalogue.cpp");
	
	stringstream ss(line);
	
	double Value_d;
	int Value_i;
	size_t ii = nObjects()-1;
	int index = 0;

	for (int jj=1; jj<=cbl::Max(column); jj++) {
	  if (varMap[column[index]]==Var::_ID_) {
	    ss>>Value_i;
	    if (std::find(column.begin(), column.end(), jj)!=column.end()) {
	      set_var(ii, varMap[column[index]], Value_i);
	      index++;
	    }
	  }
	  else {
	    ss>>Value_d;
	    if (std::find(column.begin(), column.end(), jj)!=column.end()) {
	      if ((varMap[column[index]]==Var::_RA_) || (varMap[column[index]]==Var::_Dec_)) Value_d = radians(Value_d, inputUnits);
	      set_var(ii,
		      varMap[column[index]],
		      ((varMap[column[index]]==Var::_X_) || (varMap[column[index]]==Var::_Y_) || (varMap[column[index]]==Var::_Z_)) ?
		      Value_d*fact : Value_d,
		      cosm);
	      index++;
	    }
	  }
	  	  	  
	}	
      }
    }
    finr.clear(); finr.close();
  }
}
    
// ============================================================================


cbl::catalogue::Catalogue::Catalogue (const ObjectType objectType, const std::vector<Var> attribute, const std::vector<int> column, const std::vector<std::string> file, const int comments, const double nSub, const int seed) 
{ 
  // parameters for random numbers used in case nSub!=1
  random::UniformRandomNumbers ran(0., 1., seed);

  // read the input catalogue files
 
  vector<vector<double>> vars;

  for (size_t dd=0; dd<file.size(); ++dd) {

    string line, file_in = file[dd];

    coutCBL << "I'm reading the catalogue: " << file_in << endl;
    ifstream finr(file_in.c_str()); checkIO(finr, file_in);

    // start reading catalogue
    for (int cc=0; cc<comments; cc++) getline(finr, line); // ignore commented lines at the beginning of file
    while (getline(finr, line)) { // read the lines

      if (ran()<nSub) { // extract a subsample

	stringstream ss(line);
	vector<double> num; double NN = -1.e30;
	while (ss>>NN) num.push_back(NN);

	m_object.push_back(move(Object::Create(objectType)));

	vector<double> vv;
	for (size_t i=0; i<column.size(); i++)
	  vv.push_back(num[column[i]-1]);
	
	vars.push_back(vv);
      }
    }
      
    finr.clear(); finr.close();
  }

  vars = transpose(vars);

  for (size_t i=0; i<vars.size(); i++)
    set_var(attribute[i], vars[i]);
}


// ============================================================================


size_t cbl::catalogue::Catalogue::nRegions () 
{
  if (m_nRegions==0) {
    m_nRegions = Max(Var::_Region_)+1;

    WarningMsgCBL("The total number of region is "+conv(m_nRegions, par::fINT)+", deducted from the maximum value of the input regions; if needed, this number can be changed with the function set_region_number()", "Catalogue", "Catalogue.cpp");
  }

  return m_nRegions;
}


// ============================================================================


std::vector<long> cbl::catalogue::Catalogue::region_list () const
{
  vector<long> vv(m_nRegions);

  for (size_t i=0; i<m_nRegions; ++i) vv[i] = static_cast<long>(i);

  return vv;
}

// ============================================================================


std::vector<long> cbl::catalogue::Catalogue::region () const
{ 
  vector<long> vv(m_object.size());
  
  for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->region();

  return vv;
}


// ============================================================================


std::vector<std::string> cbl::catalogue::Catalogue::field () const
{ 
  vector<string> vv(m_object.size());
  
  for (size_t i=0; i<nObjects(); ++i) vv[i] = m_object[i]->field();

  return vv;
}


// ============================================================================


double cbl::catalogue::Catalogue::var (int index, Var var_name) const
{
  double vv = 0;
  
  switch (var_name) {

  case Var::_X_:
    vv = m_object[index]->xx();
    break;

  case Var::_Y_:
    vv = m_object[index]->yy();
    break;

  case Var::_Z_:
    vv = m_object[index]->zz();
    break;

  case Var::_RA_:
    vv = m_object[index]->ra();
    break;

  case Var::_Dec_:
    vv = m_object[index]->dec();
    break;
    
  case Var::_SN_:
    vv = m_object[index]->sn();
    break;

  case Var::_Redshift_:
    vv = m_object[index]->redshift();
    break;
    
  case Var::_RedshiftMin_:
    vv = m_object[index]->redshiftMin();
    break;
    
  case Var::_RedshiftMax_:
    vv = m_object[index]->redshiftMax();
    break;
    
  case Var::_Shear1_:
    vv = m_object[index]->shear1();
    break;
    
  case Var::_Shear2_:
    vv = m_object[index]->shear2();
    break;
    
  case Var::_ODDS_:
    vv = m_object[index]->odds();
    break;
    
  case Var::_LensingWeight_:
    vv = m_object[index]->lensingWeight();
    break;
    
  case Var::_LensingCalib_:
    vv = m_object[index]->lensingCalib();
    break;

  case Var::_Dc_:
    vv = m_object[index]->dc();
    break;

  case Var::_Weight_:
    vv = m_object[index]->weight();
    break;

  case Var::_Mass_:
    vv = m_object[index]->mass();
    break;

  case Var::_Magnitude_:
    vv = m_object[index]->magnitude();
    break;
    
  case Var::_MagnitudeU_:
    vv = m_object[index]->magnitudeU();
    break;
    
  case Var::_MagnitudeG_:
    vv = m_object[index]->magnitudeG();
    break;
    
  case Var::_MagnitudeR_:
    vv = m_object[index]->magnitudeR();
    break;
    
  case Var::_MagnitudeI_:
    vv = m_object[index]->magnitudeI();
    break;

  case Var::_SFR_:
    vv = m_object[index]->SFR();
    break;

  case Var::_sSFR_:
    vv = m_object[index]->sSFR();
    break;

  case Var::_MassProxy_:
    vv = m_object[index]->mass_proxy();
    break;
    
  case Var::_MassProxyError_:
    vv = m_object[index]->mass_proxy_error();
    break;

  case Var::_Vx_:
    vv = m_object[index]->vx();
    break;
  
  case Var::_Vy_:
    vv = m_object[index]->vy();
    break;
  
  case Var::_Vz_:
    vv = m_object[index]->vz();
    break;

  case Var::_Region_:
    vv = m_object[index]->region();
    break; 
  
  case Var::_Generic_:
    vv = m_object[index]->generic();
    break;

  case Var::_Radius_:
    vv = m_object[index]->radius();
    break;
    
  case Var::_DensityContrast_:
    vv = m_object[index]->densityContrast();
    break;

  case Var::_CentralDensity_:
    vv = m_object[index]->centralDensity();
    break;

  case Var::_X_displacement_:
    vv = m_object[index]->x_displacement();
    break;
    
  case Var::_Y_displacement_:
    vv = m_object[index]->y_displacement();
    break;
    
  case Var::_Z_displacement_:
    vv = m_object[index]->z_displacement();
    break;

  case Var::_MassEstimate_:
    vv = m_object[index]->mass_estimate();
    break;
    
  case Var::_RadiusEstimate_:
    vv = m_object[index]->radius_estimate();
    break;
    
  case Var::_VeldispEstimate_:
    vv = m_object[index]->veldisp_estimate();
    break;

  case Var::_XCM_:
    vv = m_object[index]->xcm();
    break;
    
  case Var::_YCM_:
    vv = m_object[index]->ycm();
    break;
    
  case Var::_ZCM_:
    vv = m_object[index]->zcm();
    break;

  case Var::_XSpin_:
    vv = m_object[index]->spin_x();
    break;
    
  case Var::_YSpin_:
    vv = m_object[index]->spin_y();
    break;
    
  case Var::_ZSpin_:
    vv = m_object[index]->spin_z();
    break;

  case Var::_VelDisp_:
    vv = m_object[index]->veldisp();
    break;
    
  case Var::_Vmax_:
    vv = m_object[index]->vmax();
    break;
    
  case Var::_VmaxRad_:
    vv = m_object[index]->vmax_rad();
    break;
    
  case Var::_TotMass_:
    vv = m_object[index]->tot_mass();
    break;

  case Var::_ID_:
    vv = m_object[index]->ID();
    break;
    
  case Var::_Nsub_:
    vv = m_object[index]->nsub();
    break;
    
  case Var::_Parent_:
    vv = m_object[index]->parent();
    break;

  default:
    ErrorCBL("no such a variable in the list!", "var", "Catalogue.cpp");
  }
  
  return vv;
}


// ============================================================================


std::vector<double> cbl::catalogue::Catalogue::var (Var var_name) const
{ 
  vector<double> vv(m_object.size(), 0.);
  
  for (size_t i=0; i<nObjects(); ++i) vv[i] = var(i, var_name);

  return vv;
}


// ============================================================================


bool cbl::catalogue::Catalogue::isSetVar (int index, Var var_name) const
{
  if (var_name==Var::_X_)
    return m_object[index]->isSet_xx();

  else if (var_name==Var::_Y_)
    return m_object[index]->isSet_yy();

  else if (var_name==Var::_Z_)
    return m_object[index]->isSet_zz();

  else if (var_name==Var::_RA_)
    return m_object[index]->isSet_ra();

  else if (var_name==Var::_Dec_)
    return m_object[index]->isSet_dec();
    
  else if (var_name==Var::_SN_)
    return m_object[index]->isSet_sn();

  else if (var_name==Var::_Redshift_)
    return m_object[index]->isSet_redshift();
    
  else if (var_name==Var::_RedshiftMin_)
    return m_object[index]->isSet_redshiftMin();
    
  else if (var_name==Var::_RedshiftMax_)
    return m_object[index]->isSet_redshiftMax();
    
  else if (var_name==Var::_Shear1_)
    return m_object[index]->isSet_shear1();
    
  else if (var_name==Var::_Shear2_)
    return m_object[index]->isSet_shear2();
    
  else if (var_name==Var::_ODDS_)
    return m_object[index]->isSet_odds();
    
  else if (var_name==Var::_LensingWeight_)
    return m_object[index]->isSet_lensingWeight();
    
  else if (var_name==Var::_LensingCalib_)
    return m_object[index]->isSet_lensingCalib();

  else if (var_name==Var::_Dc_)
    return m_object[index]->isSet_dc();

  else if (var_name==Var::_Weight_)
    return m_object[index]->isSet_weight();

  else if (var_name==Var::_Mass_)
    return m_object[index]->isSet_mass();

  else if (var_name==Var::_Magnitude_)
    return m_object[index]->isSet_magnitude();
    
  else if (var_name==Var::_MagnitudeU_)
    return m_object[index]->isSet_magnitudeU();
    
  else if (var_name==Var::_MagnitudeG_)
    return m_object[index]->isSet_magnitudeG();
    
  else if (var_name==Var::_MagnitudeR_)
    return m_object[index]->isSet_magnitudeR();
    
  else if (var_name==Var::_MagnitudeI_)
    return m_object[index]->isSet_magnitudeI();

  else if (var_name==Var::_SFR_)
    return m_object[index]->isSet_SFR();

  else if (var_name==Var::_sSFR_)
    return m_object[index]->isSet_sSFR();
    
  else if (var_name==Var::_MassProxy_)
    return m_object[index]->isSet_mass_proxy();

  else if (var_name==Var::_MassProxyError_)
    return m_object[index]->isSet_mass_proxy_error();

  else if (var_name==Var::_Vx_)
    return m_object[index]->isSet_vx();

  else if (var_name==Var::_Vy_)
    return m_object[index]->isSet_vy();

  else if (var_name==Var::_Vz_)
    return m_object[index]->isSet_vz();

  else if (var_name==Var::_Region_)
    return m_object[index]->isSet_region();

  else if (var_name==Var::_Generic_)
    return m_object[index]->isSet_generic();

  else if (var_name==Var::_Radius_)
    return m_object[index]->isSet_radius();

  else if (var_name==Var::_DensityContrast_)
    return m_object[index]->isSet_densityContrast();

  else if (var_name==Var::_CentralDensity_)
    return m_object[index]->isSet_centralDensity();

  else if (var_name==Var::_X_displacement_)
    return m_object[index]->isSet_x_displacement();

  else if (var_name==Var::_Y_displacement_)
    return m_object[index]->isSet_y_displacement();
  
  else if (var_name==Var::_Z_displacement_)
    return m_object[index]->isSet_z_displacement();

  else if (var_name==Var::_MassEstimate_)
    return m_object[index]->isSet_mass_estimate();

  else if (var_name==Var::_RadiusEstimate_)
    return m_object[index]->isSet_radius_estimate();

  else if (var_name==Var::_VeldispEstimate_)
    return m_object[index]->isSet_veldisp_estimate();

  else if (var_name==Var::_XCM_)
    return m_object[index]->isSet_xcm();

  else if (var_name==Var::_YCM_)
    return m_object[index]->isSet_ycm();

  else if (var_name==Var::_ZCM_)
    return m_object[index]->isSet_zcm();

  else if (var_name==Var::_XSpin_)
    return m_object[index]->isSet_spin_x();

  else if (var_name==Var::_YSpin_)
    return m_object[index]->isSet_spin_y();

  else if (var_name==Var::_ZSpin_)
    return m_object[index]->isSet_spin_z();

  else if (var_name==Var::_VelDisp_)
    return m_object[index]->isSet_veldisp();

  else if (var_name==Var::_Vmax_)
    return m_object[index]->isSet_vmax();

  else if (var_name==Var::_VmaxRad_)
    return m_object[index]->isSet_vmax_rad();

  else if (var_name==Var::_TotMass_)
    return m_object[index]->isSet_tot_mass();

  else if (var_name==Var::_ID_)
    return m_object[index]->isSet_ID();

  else
    return ErrorCBL("no such a variable in the list!", "isSetVar", "Catalogue.cpp");
}


// ============================================================================


bool cbl::catalogue::Catalogue::isSetVar (Var var_name) const
{
  bool ret = true;

  size_t i = 0;
  while (ret==true && i<nObjects()) 
    ret = isSetVar(i++, var_name);

  return ret;
}


// ============================================================================


void cbl::catalogue::Catalogue::set_region (const std::vector<long> region, const int nRegions)
{
  for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_region(region[i]);

  set_region_number( static_cast<size_t>( (nRegions<0) ? cbl::Max<long>(region)+1 : nRegions) );
}


// ============================================================================


void cbl::catalogue::Catalogue::set_region_number (const size_t nRegions)
{
  for (size_t i=0; i<m_object.size(); i++)
    if (m_object[i]->region()>=static_cast<int>(nRegions))
      ErrorCBL("region index for object "+conv(i, par::fINT)+" is larger than input number of regions! "+conv(m_object[i]->region(), par::fINT)+" >= "+conv(nRegions, par::fINT), "set_region_number", "Catalogue.cpp");

  m_nRegions = nRegions;
}


// ============================================================================


void cbl::catalogue::Catalogue::set_field (const std::vector<std::string> field)
{
  for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_field(field[i]);
}


// ============================================================================


void cbl::catalogue::Catalogue::set_var (const int index, const Var var_name, const double value, const cosmology::Cosmology cosmology, const bool update_coordinates)
{
  switch (var_name) {

  case Var::_X_:
    m_object[index]->set_xx(value);
    break;

  case Var::_Y_:
    m_object[index]->set_yy(value);
    break;

  case Var::_Z_:
    m_object[index]->set_zz(value);
    break;

  case Var::_RA_:
    m_object[index]->set_ra(value);
    break;

  case Var::_Dec_:
    m_object[index]->set_dec(value);
    break;
    
  case Var::_SN_:
    m_object[index]->set_sn(value);
    break;

  case Var::_Redshift_:
    m_object[index]->set_redshift(value, cosmology, update_coordinates);
    break;
    
  case Var::_RedshiftMin_:
    m_object[index]->set_redshiftMin(value);
    break;
    
  case Var::_RedshiftMax_:
    m_object[index]->set_redshiftMax(value);
    break;
    
  case Var::_Shear1_:
    m_object[index]->set_shear1(value);
    break;
    
  case Var::_Shear2_:
    m_object[index]->set_shear2(value);
    break;
    
  case Var::_ODDS_:
    m_object[index]->set_odds(value);
    break;
    
  case Var::_LensingWeight_:
    m_object[index]->set_lensingWeight(value);
    break;
    
  case Var::_LensingCalib_:
    m_object[index]->set_lensingCalib(value);
    break;    

  case Var::_Dc_:
    m_object[index]->set_dc(value);
    break;

  case Var::_Weight_:
    m_object[index]->set_weight(value);
    break;

  case Var::_Mass_:
    m_object[index]->set_mass(value);
    break;

  case Var::_Magnitude_:
    m_object[index]->set_magnitude(value);
    break;
    
  case Var::_MagnitudeU_:
    m_object[index]->set_magnitudeU(value);
    break;
    
  case Var::_MagnitudeG_:
    m_object[index]->set_magnitudeG(value);
    break;
    
  case Var::_MagnitudeR_:
    m_object[index]->set_magnitudeR(value);
    break;
    
  case Var::_MagnitudeI_:
    m_object[index]->set_magnitudeI(value);
    break;

  case Var::_SFR_:
    m_object[index]->set_SFR(value);
    break;

  case Var::_sSFR_:
    m_object[index]->set_sSFR(value);
    break;
    
  case Var::_MassProxy_:
    m_object[index]->set_mass_proxy(value);
    break;
    
  case Var::_MassProxyError_:
    m_object[index]->set_mass_proxy_error(value);
    break;

  case Var::_Vx_:
    m_object[index]->set_vx(value);
    break;
  
  case Var::_Vy_:
    m_object[index]->set_vy(value);
    break;
  
  case Var::_Vz_:
    m_object[index]->set_vz(value);
    break;

  case Var::_Region_:
    m_object[index]->set_region(value);
    break;
    
  case Var::_Generic_:
    m_object[index]->set_generic(value);
    break;

  case Var::_Radius_:
    m_object[index]->set_radius(value);
    break;

  case Var::_CentralDensity_:
    m_object[index]->set_centralDensity(value);
    break;

  case Var::_DensityContrast_:
    m_object[index]->set_densityContrast(value);
    break;

  case Var::_X_displacement_:
    m_object[index]->set_x_displacement(value);
    break;
    
  case Var::_Y_displacement_:
    m_object[index]->set_y_displacement(value);
    break;
    
  case Var::_Z_displacement_:
    m_object[index]->set_z_displacement(value);
    break;
    
  case Var::_MassEstimate_:
    m_object[index]->set_mass_estimate(value);
    break;
    
  case Var::_RadiusEstimate_:
    m_object[index]->set_radius_estimate(value);
    break;
    
  case Var::_VeldispEstimate_:
    m_object[index]->set_veldisp_estimate(value);
    break;
    
  case Var::_XCM_:
    m_object[index]->set_xcm(value);
    break;
    
  case Var::_YCM_:
    m_object[index]->set_ycm(value);
    break;
    
  case Var::_ZCM_:
    m_object[index]->set_zcm(value);
    break;
    
  case Var::_XSpin_:
    m_object[index]->set_spin_x(value);
    break;
    
  case Var::_YSpin_:
    m_object[index]->set_spin_y(value);
    break;
    
  case Var::_ZSpin_:
    m_object[index]->set_spin_z(value);
    break;
    
  case Var::_VelDisp_:
    m_object[index]->set_veldisp(value);
    break;
    
  case Var::_Vmax_:
    m_object[index]->set_vmax(value);
    break;
    
  case Var::_VmaxRad_:
    m_object[index]->set_vmax_rad(value);
    break;
    
  case Var::_TotMass_:
    m_object[index]->set_tot_mass(value);
    break;
    
  default:
    ErrorCBL("no such a variable in the list!", "set_var", "Catalogue.cpp");
  }

}


// ============================================================================


void cbl::catalogue::Catalogue::set_var (const int index, const Var var_name, const int value, const cosmology::Cosmology cosmology)
{
  switch (var_name) {

  case Var::_ID_:
    m_object[index]->set_ID(value);
    (void)cosmology;
    break;
    
  default:
    ErrorCBL("no such a variable in the list!", "set_var", "Catalogue.cpp");
  }

}
    

// ============================================================================


void cbl::catalogue::Catalogue::set_var (const Var var_name, const std::vector<double> var, const cosmology::Cosmology cosmology, const bool update_coordinates)
{
  if (m_object.size()!=var.size()) ErrorCBL("m_object.size()!=var.size()!", "set_var", "Catalogue.cpp");
  
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
    
  case Var::_SN_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_sn(var[i]);
    break;

  case Var::_Redshift_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_redshift(var[i], cosmology, update_coordinates);
    break;
    
  case Var::_RedshiftMin_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_redshiftMin(var[i]);
    break;
    
  case Var::_RedshiftMax_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_redshiftMax(var[i]);
    break;
    
  case Var::_Shear1_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_shear1(var[i]);
    break;
    
  case Var::_Shear2_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_shear2(var[i]);
    break;
    
  case Var::_ODDS_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_odds(var[i]);
    break;
    
  case Var::_LensingWeight_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_lensingWeight(var[i]);
    break;
    
  case Var::_LensingCalib_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_lensingCalib(var[i]);
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
    
  case Var::_MagnitudeU_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_magnitudeU(var[i]);
    break;
    
  case Var::_MagnitudeG_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_magnitudeG(var[i]);
    break;
    
  case Var::_MagnitudeR_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_magnitudeR(var[i]);
    break;
    
  case Var::_MagnitudeI_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_magnitudeI(var[i]);
    break;

  case Var::_SFR_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_SFR(var[i]);
    break;

  case Var::_sSFR_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_sSFR(var[i]);
    break;
    
  case Var::_MassProxy_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_mass_proxy(var[i]);
    break;
    
  case Var::_MassProxyError_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_mass_proxy_error(var[i]);
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
    {
      vector<long> regList(nObjects());
      for (size_t i=0; i<nObjects(); ++i) regList[i] = static_cast<long>(var[i]);
      set_region(regList);
      break;
    }
    
  case Var::_Generic_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_generic(var[i]);
    break;

  case Var::_Radius_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_radius(var[i]);
    break;

  case Var::_CentralDensity_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_centralDensity(var[i]);
    break;

  case Var::_DensityContrast_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_densityContrast(var[i]);
    break;

  case Var::_X_displacement_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_x_displacement(var[i]);
    break;
    
  case Var::_Y_displacement_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_y_displacement(var[i]);
    break;
    
  case Var::_Z_displacement_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_z_displacement(var[i]);
    break;
    
  case Var::_MassEstimate_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_mass_estimate(var[i]);
    break;
    
  case Var::_RadiusEstimate_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_radius_estimate(var[i]);
    break;
    
  case Var::_VeldispEstimate_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_veldisp_estimate(var[i]);
    break;
    
  case Var::_XCM_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_xcm(var[i]);
    break;
    
  case Var::_YCM_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_ycm(var[i]);
    break;
    
  case Var::_ZCM_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_zcm(var[i]);
    break;
    
  case Var::_XSpin_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_spin_x(var[i]);
    break;
    
  case Var::_YSpin_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_spin_y(var[i]);
    break;
    
  case Var::_ZSpin_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_spin_z(var[i]);
    break;
    
  case Var::_VelDisp_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_veldisp(var[i]);
    break;
    
  case Var::_Vmax_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_vmax(var[i]);
    break;
    
  case Var::_VmaxRad_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_vmax_rad(var[i]);
    break;
    
  case Var::_TotMass_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_tot_mass(var[i]);
    break;

  default:
    ErrorCBL("no such a variable in the list!", "set_var", "Catalogue.cpp");
  }

}

// ============================================================================


void cbl::catalogue::Catalogue::set_var (const Var var_name, const std::vector<int> var, const cosmology::Cosmology cosmology)
{
  if (m_object.size()!=var.size()) ErrorCBL("m_object.size()!=var.size()!", "set_var", "Catalogue.cpp");
  
  switch (var_name) {
    
  case Var::_ID_:
    for (size_t i=0; i<nObjects(); ++i) m_object[i]->set_ID(var[i]);
    (void)cosmology;
    break;

  default:
    ErrorCBL("no such a variable in the list!", "set_var", "Catalogue.cpp");
  }

}

// ============================================================================


void cbl::catalogue::Catalogue::stats_var (const Var var_name, std::vector<double> &stats) const
{
  stats.erase(stats.begin(), stats.end());
  stats.resize(4);
  
  stats[0] = Average(var(var_name)); 
  stats[2] = Sigma(var(var_name));
  
  stats[1] = Quartile(var(var_name))[1];
  stats[3] = Quartile(var(var_name))[2]-Quartile(var(var_name))[0];
}


// ============================================================================


void cbl::catalogue::Catalogue::stats_var (const std::vector<Var> var_name, std::vector<std::vector<double>> &stats) const
{
  stats.erase(stats.begin(),stats.end());

  for (unsigned int i=0; i<var_name.size(); i++) {

    vector<double> stats_temp;
    stats_var(var_name[i],stats_temp);

    stats.push_back(stats_temp);
  }
}


// ============================================================================


void cbl::catalogue::Catalogue::var_distr (const Var var_name, std::vector<double> &_var, std::vector<double> &dist, std::vector<double> &err, const int nbin, const bool linear, const std::string file_out, const double Volume, const bool norm, const double V1, const double V2, const std::string bin_type, const bool convolution, const double sigma) const
{
  distribution(_var, dist, err, var(var_name), var(Var::_Weight_), nbin, linear, file_out, (norm) ? Volume*weightedN() : Volume, V1, V2, bin_type, convolution, sigma);
}


// ============================================================================


void cbl::catalogue::Catalogue::computeComovingCoordinates (const cosmology::Cosmology &cosm, const CoordinateUnits inputUnits)
{
  double red, xx, yy, zz;

  
  // ----- unit conversion -----
  
  vector<double> RA(nObjects()), DEC(nObjects());

  for (size_t i=0; i<nObjects(); ++i) {
    RA[i] = (inputUnits==CoordinateUnits::_radians_) ? ra(i) : radians(ra(i), inputUnits);
    DEC[i] = (inputUnits==CoordinateUnits::_radians_) ? dec(i) : radians(dec(i), inputUnits);
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


void cbl::catalogue::Catalogue::computePolarCoordinates (const CoordinateUnits outputUnits)
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
  
  if (outputUnits!=CoordinateUnits::_radians_) {    
    
    if (outputUnits==CoordinateUnits::_degrees_)
      for (size_t i=0; i<nObjects(); ++i) {
        ra = m_object[i]->ra();
        dec = m_object[i]->dec();
        m_object[i]->set_ra(degrees(ra));
        m_object[i]->set_dec(degrees(dec));
      }

    else if (outputUnits==CoordinateUnits::_arcseconds_)
      for (size_t i=0; i<nObjects(); ++i) {
        ra = m_object[i]->ra();
        dec = m_object[i]->dec();
        m_object[i]->set_ra(arcseconds(ra));
        m_object[i]->set_dec(arcseconds(dec)); 		    
      }

    else if (outputUnits==CoordinateUnits::_arcminutes_)
      for (size_t i=0; i<nObjects(); ++i) {
        ra = m_object[i]->ra();
        dec = m_object[i]->dec();
        m_object[i]->set_ra(arcminutes(ra));
        m_object[i]->set_dec(arcminutes(dec)); 		    
      }

    else ErrorCBL("outputUnits type not allowed!", "computePolarCoordinates", "Catalogue.cpp");
  }
  
}

// ============================================================================


void cbl::catalogue::Catalogue::computePolarCoordinates (const cosmology::Cosmology &cosmology, const double z1, const double z2, const CoordinateUnits outputUnits)
{
  double ra, dec, dc;

  
  // ----- compute polar coordinates in radians -----

  for (size_t i=0; i<nObjects(); ++i) {
    polar_coord(xx(i), yy(i), zz(i), ra, dec, dc);
    m_object[i]->set_ra(ra); 
    m_object[i]->set_dec(dec); 
    m_object[i]->set_dc(dc);
    m_object[i]->set_redshift(cosmology.Redshift(dc, z1, z2), cosmology);
  }

  
  // ----- unit conversion -----

  if (outputUnits!=CoordinateUnits::_radians_) {

    if (outputUnits==CoordinateUnits::_degrees_)
      for (size_t i=0; i<nObjects(); ++i) {
	m_object[i]->set_ra(degrees(ra));
	m_object[i]->set_dec(degrees(dec)); 		    
      }

    else if (outputUnits==CoordinateUnits::_arcseconds_)
      for (size_t i=0; i<nObjects(); ++i) {
	m_object[i]->set_ra(arcseconds(ra));
	m_object[i]->set_dec(arcseconds(dec)); 		    
      }

    else if (outputUnits==CoordinateUnits::_arcminutes_)
      for (size_t i=0; i<nObjects(); ++i) {
	m_object[i]->set_ra(arcminutes(ra));
	m_object[i]->set_dec(arcminutes(dec)); 		    
      }

    else ErrorCBL("outputUnits type not allowed!", "computePolarCoordinates", "Catalogue.cpp");
  }
  
}

// ============================================================================


void cbl::catalogue::Catalogue::normalizeComovingCoordinates () 
{
  for (size_t i=0; i<nObjects(); ++i) { 
    m_object[i]->set_xx(xx(i)/dc(i)); 
    m_object[i]->set_yy(yy(i)/dc(i)); 
    m_object[i]->set_zz(zz(i)/dc(i));
  }
}


// ============================================================================


void cbl::catalogue::Catalogue::restoreComovingCoordinates ()
{
  for (size_t i=0; i<nObjects(); ++i) {
    m_object[i]->set_xx(xx(i)*dc(i)); 
    m_object[i]->set_yy(yy(i)*dc(i)); 
    m_object[i]->set_zz(zz(i)*dc(i));
  }
}


// ============================================================================


void cbl::catalogue::Catalogue::Order (const std::vector<int> vv) 
{
  const int nObj = m_object.size();
  
  if (int(vv.size())!=nObj)
    ErrorCBL("vv.size()="+conv(vv.size(), par::fINT)+" and nObj="+conv(nObj, par::fINT)+" must be equal!", "Order", "Catalogue.cpp");
 
  vector<shared_ptr<Object>> obj(nObj);
  
  m_index.resize(nObj);
  
  for (size_t i=0; i<vv.size(); i++) {
    m_index[i] = vv[i];
    obj[i] = m_object[vv[i]];
  }

  m_object = obj;
}


// ============================================================================


void cbl::catalogue::Catalogue::Order () 
{
  if (m_index.size()>0) {
  
    const size_t nObj = m_object.size();

    if (m_index.size()!=nObj) 
      ErrorCBL("m_index.size()="+conv(m_index.size(), par::fINT)+" and nObj="+conv(nObj, par::fINT)+" must be equal!", "Order", "Catalogue.cpp");
  
    vector<shared_ptr<Object>> obj(nObj);
  
    obj = m_object;
  
    for (size_t i=0; i<nObj; i++) 
      m_object[i] = obj[m_index[i]];

  }
}


// ============================================================================


double cbl::catalogue::Catalogue::weightedN () const
{
  double nn = 0.;
  for (size_t i=0; i<m_object.size(); i++)
    nn += m_object[i]->weight();
  return nn;
}


// ============================================================================


void cbl::catalogue::Catalogue::write_comoving_coordinates (const std::string outputFile) const
{
  if (m_object.size()==0) ErrorCBL("m_object.size()=0!", "write_comoving_coordinates", "Catalogue.cpp");

  coutCBL << "I'm writing the file: " << outputFile << "..." << endl;
  
  ofstream fout(outputFile.c_str()); checkIO(fout, outputFile);
  
  for (size_t i=0; i<nObjects(); ++i) 
    fout << xx(i) << "   " << yy(i) << "   " << zz(i) << endl;

  coutCBL << "I wrote the file: " << outputFile << endl;
  fout.clear(); fout.close();
}


// ============================================================================


void cbl::catalogue::Catalogue::write_obs_coordinates (const std::string outputFile) const 
{
  if (m_object.size()==0) ErrorCBL("m_object.size()=0!", "write_obs_coordinates", "Catalogue.cpp");

  coutCBL << "I'm writing the file: " << outputFile << "..." << endl;
  
  ofstream fout(outputFile.c_str()); checkIO(fout, outputFile);
  
  if (!isSet(ra(0)) || !isSet(dec(0)) || !isSet(redshift(0)))
    ErrorCBL("polar coordinates are not set!", "write_obs_coordinates", "Catalogue.cpp");

  for (size_t i=0; i<nObjects(); ++i) 
    fout << ra(i) << "   " << dec(i) << "   " << redshift(i) << endl;
  
  coutCBL << "I wrote the file: " << outputFile << endl;
  fout.clear(); fout.close();
}


// ============================================================================


void cbl::catalogue::Catalogue::write_data (const std::string outputFile, const std::vector<Var> var_name) const 
{
  coutCBL << "I'm writing the file: " << outputFile << "..." << endl;
  
  ofstream fout(outputFile.c_str()); checkIO(fout, outputFile);

  if (var_name.size()==0) {
    
    if (!isSet(ra(0)) || !isSet(dec(0)) || !isSet(redshift(0)) || !isSet(dc(0)))
      ErrorCBL("polar coordinates are not set!", "write_data", "Catalogue.cpp");
    
    if (!isSet(region(0)))
      for (size_t i=0; i<nObjects(); ++i) 
	fout << xx(i) << "   " << yy(i) << "   " << zz(i) << "   " << ra(i) << "   " << dec(i) << "   " << redshift(i) << "   " << dc(i) << endl;
    else
      for (size_t i=0; i<nObjects(); ++i) 
	fout << xx(i) << "   " << yy(i) << "   " << zz(i) << "   " << ra(i) << "   " << dec(i) << "   " << redshift(i) << "   " << dc(i) << "   " << region(i) <<endl;

  }

  else {
    vector<vector<double>> data;
    for (size_t j=0; j<var_name.size(); j++)
      data.push_back(var(var_name[j]));

    for (size_t i=0; i<nObjects(); ++i) {
      for (size_t j=0; j<data.size(); j++)
	fout << data[j][i] << "   ";
      fout << endl;
    }
    
  }

  coutCBL << "I wrote the file: " << outputFile << endl;
  fout.clear(); fout.close();
}


// ============================================================================


Catalogue cbl::catalogue::Catalogue::sub_catalogue (const mask_function mask, const bool excl) const
{
  vector<shared_ptr<Object>> objects;
  vector<bool> w(m_object.size());
  
  bool fact = (excl) ? false : true;

  for (size_t i=0; i<m_object.size(); i++) 
    w[i] = mask(m_object[i])&&fact;

  for (size_t i=0; i<m_object.size(); i++)
    if (w[i])
      objects.push_back(m_object[i]);
 
  return Catalogue{objects};
}


// ============================================================================


Catalogue cbl::catalogue::Catalogue::sub_catalogue (const cbl::catalogue::MaskObject &mask, const bool excl) const
{
  vector<shared_ptr<Object>> objects;
  vector<bool> w(m_object.size());
  
  bool fact = (excl) ? false : true;

  for (size_t i=0; i<m_object.size(); i++) 
    w[i] = mask(m_object[i])&&fact;

  for (size_t i=0; i<m_object.size(); i++)
    if (w[i])
      objects.push_back(m_object[i]);
 
  return Catalogue{objects};
}


// ============================================================================


Catalogue cbl::catalogue::Catalogue::sub_catalogue (const Var var_name, const double down, const double up, const bool excl) const
{
  vector<shared_ptr<Object>> objects;
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
 
  return Catalogue{objects};
}


// ============================================================================


Catalogue cbl::catalogue::Catalogue::mangle_cut (const std::string mangle_mask, const bool excl) const
{

  vector<shared_ptr<Object>> objects;

  string mangle_dir = par::DirCosmo+"/External/mangle/";

  string mangle_working_dir = mangle_dir+"output/";
  string mkdir = "mkdir -p "+mangle_working_dir;
  if (system(mkdir.c_str())) {}

  string input = mangle_working_dir+"temporary.dat";
  string output = mangle_working_dir+"temporary_output.dat";

  write_obs_coordinates(input);

  string polyid = mangle_dir+"/bin/polyid -ur "+mangle_mask+" "+input+" "+output;
  if (system(polyid.c_str())) {}

  ifstream fin(output);
  string line;

  int nn=0;
  getline(fin, line);

  while(getline(fin, line))
    {
      stringstream ss(line);
      vector<double> num; double NN = par::defaultDouble;
      while (ss>>NN) num.push_back(NN);
    
      if( !excl && num.size()>2)
	objects.push_back(m_object[nn]);
      if( excl && num.size()<3)
	objects.push_back(m_object[nn]);
      nn++;
    }
  fin.clear(); fin.close();
  
  string rm = "rm "+input+" "+output;
  if (system(rm.c_str())) {}

  return Catalogue{objects};
}

// ============================================================================


Catalogue cbl::catalogue::Catalogue::diluted_catalogue (const double nSub, const int seed) const
{
  if (nSub<=0 || nSub>1 || !isfinite(nSub)) ErrorCBL("nSub must be in the range (0,1] !", "diluted_catalogue", "Catalogue.cpp");
  
  // copy the catalogue into a new one 
  auto diluted_catalogue = *this;
  
  // set the index vector that will be used to remove the objects
  vector<bool> index(nObjects(), false);

  // nObjects()*(1-nSub) will be removed
  for (size_t i=0; i<nObjects()*(1-nSub); ++i) index[i] = true;
  
  // shuffle the indexes of the objects that will be removed
  default_random_engine engine(seed);
  std::shuffle(begin(index), end(index), engine);
  
  // dilute the new catalogue
  diluted_catalogue.remove_objects(index);
  
  return diluted_catalogue;
}

  
// ============================================================================

  
double cbl::catalogue::Catalogue::distance (const int i, shared_ptr<Object> obj) const
{
  return sqrt((m_object[i]->xx()-obj->xx())*(m_object[i]->xx()-obj->xx())+
	      (m_object[i]->yy()-obj->yy())*(m_object[i]->yy()-obj->yy())+
	      (m_object[i]->zz()-obj->zz())*(m_object[i]->zz()-obj->zz()));
}


// ============================================================================


double cbl::catalogue::Catalogue::angsep_xyz (const int i, shared_ptr<Object> obj) const
{ 
  return 2.*asin(0.5*sqrt((m_object[i]->xx()-obj->xx())*(m_object[i]->xx()-obj->xx())+
			  (m_object[i]->yy()-obj->yy())*(m_object[i]->yy()-obj->yy())+
			  (m_object[i]->zz()-obj->zz())*(m_object[i]->zz()-obj->zz())));
}
    

// ============================================================================


shared_ptr<Catalogue> cbl::catalogue::Catalogue::smooth (const double gridsize, const cosmology::Cosmology cosmology, const std::vector<Var> vars, const int SUB)
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
  double MinX = Min(Var::_X_);
  double MinY = Min(Var::_Y_);
  double MinZ = Min(Var::_Z_);

  vector<long> regions(cat->nObjects());

  for (size_t i=0; i<cat->nObjects(); i++) {
    int i1 = min(int((cat->xx(i)-MinX)/Cell_X), nx-1);
    int j1 = min(int((cat->yy(i)-MinY)/Cell_Y), ny-1);
    int z1 = min(int((cat->zz(i)-MinZ)/Cell_Z), nz-1);
    long index = z1+nz*(j1+ny*i1);
    regions[i] = index;
  }

  cat->set_region(regions);

  size_t nRegions = cat->nRegions();

  vector<Catalogue> subSamples(nRegions);
  
  for (size_t i=0; i<nRegions; ++i) {
    double start = (double)cat->region_list()[i];
    double stop = start+1;
    subSamples[i] = sub_catalogue(Var::_Region_, start, stop);
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
	REDSHIFT /= nObj; obj->set_redshift(REDSHIFT, cosmology);
	obj->set_weight(WEIGHT);
	sample.push_back(obj);

      }
    }
    
  }
  
  shared_ptr<Catalogue> cat_new(new Catalogue{sample});
  return cat_new;
}


// ============================================================================


int cbl::catalogue::Catalogue::nObjects_condition (const Var var_name, const double down, const double up, const bool excl)
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


double cbl::catalogue::Catalogue::weightedN_condition (const Var var_name, const double down, const double up, const bool excl)
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


data::ScalarField3D cbl::catalogue::Catalogue::counts_in_cell (const double cell_size, const int interpolation_type, const bool useMass, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ) const
{
  double _minX = (minX>par::defaultDouble) ? minX : Min(Var::_X_)-3*cell_size;
  double _maxX = (maxX>par::defaultDouble) ? maxX : Max(Var::_X_)+3*cell_size;
  double _minY = (minY>par::defaultDouble) ? minY : Min(Var::_Y_)-3*cell_size;
  double _maxY = (maxY>par::defaultDouble) ? maxY : Max(Var::_Y_)+3*cell_size;
  double _minZ = (minZ>par::defaultDouble) ? minZ : Min(Var::_Z_)-3*cell_size;
  double _maxZ = (maxZ>par::defaultDouble) ? maxZ : Max(Var::_Z_)+3*cell_size;

  data::ScalarField3D density(cell_size, _minX, _maxX, _minY, _maxY, _minZ, _maxZ);

  double deltaX = density.deltaX();
  double deltaY = density.deltaY();
  double deltaZ = density.deltaZ();
  int nx = density.nx();
  int ny = density.ny();
  int nz = density.nz();

  for (size_t i=0; i<nObjects(); ++i) {
    int i1 = min(int((xx(i)-density.MinX())/deltaX),nx-1);
    int j1 = min(int((yy(i)-density.MinY())/deltaY),ny-1);
    int k1 = min(int((zz(i)-density.MinZ())/deltaZ),nz-1);

    double w = (useMass) ? mass(i)*weight(i) : weight(i);

    if (interpolation_type==0) {
      density.set_ScalarField(w, i1, j1, k1, 1);
    }
    else if (interpolation_type==1) {

      double dx = (xx(i)-density.XX(i1))/deltaX;
      double dy = (yy(i)-density.YY(j1))/deltaY;
      double dz = (zz(i)-density.ZZ(k1))/deltaZ;

      int i2 = (dx<0) ? i1-1 : i1+1;
      int j2 = (dy<0) ? j1-1 : j1+1;
      int k2 = (dz<0) ? k1-1 : k1+1;

      dx = fabs(dx);
      dy = fabs(dy);
      dz = fabs(dz);

      double tx = 1-dx;
      double ty = 1-dy;
      double tz = 1-dz;

      double ww = 0.;
      ww += w*tx*ty*tz;
      ww += w*dx*ty*tz;
      ww += w*tx*dy*tz;
      ww += w*tx*ty*dz;
      ww += w*dx*dy*tz;
      ww += w*dx*ty*dz;
      ww += w*tx*dy*dz;
      ww += w*dx*dy*dz;

      density.set_ScalarField(w*tx*ty*tz, i1, j1, k1, 1);
      density.set_ScalarField(w*dx*ty*tz, i2, j1, k1, 1);
      density.set_ScalarField(w*tx*dy*tz, i1, j2, k1, 1);
      density.set_ScalarField(w*tx*ty*dz, i1, j1, k2, 1);
      density.set_ScalarField(w*dx*dy*tz, i2, j2, k1, 1);
      density.set_ScalarField(w*dx*ty*dz, i2, j1, k2, 1);
      density.set_ScalarField(w*tx*dy*dz, i1, j2, k2, 1);
      density.set_ScalarField(w*dx*dy*dz, i2, j2, k2, 1);

    }
  }

  return density;
}


// ============================================================================


data::ScalarField3D cbl::catalogue::Catalogue::density_field (const double cell_size, const Catalogue mask_catalogue, const int interpolation_type, const double kernel_radius, const bool useMass) const
{

  data::ScalarField3D data_cic = counts_in_cell(cell_size, interpolation_type, useMass);

  data::ScalarField3D mask_cic = mask_catalogue.counts_in_cell(cell_size, 0 /*interpolation_type*/, useMass, data_cic.MinX(), data_cic.MaxX(), data_cic.MinY(), data_cic.MaxY(), data_cic.MinZ(), data_cic.MaxZ());


  data::ScalarField3D density(cell_size, data_cic.MinX(), data_cic.MaxX(), data_cic.MinY(), data_cic.MaxY(), data_cic.MinZ(), data_cic.MaxZ());

  double data_tot=0, random_tot=0;
  int nrandom = 0;

  for (int i=0; i<density.nx(); i++) 
    for (int j=0; j<density.ny(); j++) 
      for (int k=0; k<density.nz(); k++) {
	if (mask_cic.ScalarField(i, j, k)>0) {
	  random_tot += mask_cic.ScalarField(i, j, k);	
	  nrandom ++;
	}
      }

  double mean_random = random_tot/nrandom;
  coutCBL << "Mean random objects " << mean_random << " in " << nrandom << " cells " << endl;

  random_tot=0;
  int masked_cells=0;

  for (int i=0; i<density.nx(); i++) 
    for (int j=0; j<density.ny(); j++) 
      for (int k=0; k<density.nz(); k++) {
	if (mask_cic.ScalarField(i, j, k)>0.) {
	  data_tot += data_cic.ScalarField(i, j, k);
	  random_tot += mask_cic.ScalarField(i, j, k);
	}
	else if (mask_cic.ScalarField(i, j, k)>0 && mask_cic.ScalarField(i, j, k)<0.1*mean_random)
	  {
	    masked_cells ++;
	    data_cic.set_ScalarField(0, i, j, k, 0);
	    mask_cic.set_ScalarField(0, i, j, k, 0);
	  }
      }

  coutCBL << "Masked " << masked_cells << "/" << nrandom << " for bad random coverage " << endl;

  double norm = int(random_tot)/data_tot;
  for (int i=0; i<density.nx(); i++) 
    for (int j=0; j<density.ny(); j++) 
      for (int k=0; k<density.nz(); k++) {
        double val=0;
	val = (mask_cic.ScalarField(i, j, k)>0) ? data_cic.ScalarField(i,j,k)/mask_cic.ScalarField(i,j,k)*norm-1 : 0;
        density.set_ScalarField(val, i, j, k);
      }
  
  if (kernel_radius>0)
    density.GaussianConvolutionField(kernel_radius);

  return density;
}


// ============================================================================


cbl::catalogue::Catalogue::Catalogue (const Catalogue input_catalogue, const Catalogue target_catalogue, const Var var_name, const int nbin, const int seed) 
{
  vector<double> fvar_input(nbin, 0);
  vector<double> fvar_target(nbin, 0);

  vector<double> input_var = input_catalogue.var(var_name);
  vector<double> target_var = target_catalogue.var(var_name);

  double Vmin = target_catalogue.Min(var_name);
  double Vmax = target_catalogue.Max(var_name);

  double binSize_inv = pow((Vmax-Vmin)/nbin, -1);

  for (size_t i=0; i<input_var.size(); i++)
    if (input_var[i] < Vmax && Vmin < input_var[i]) {
      int occ = max(0, min(int((input_var[i]-Vmin)*binSize_inv), nbin));
      fvar_input[occ] ++;
    }

  for (size_t i=0; i<target_var.size(); i++)
    if (target_var[i] < Vmax && Vmin < target_var[i]) {
      int occ = max(0, min(int((target_var[i]-Vmin)*binSize_inv), nbin));
      fvar_target[occ] ++;
    }

  random::UniformRandomNumbers ran(0., 1., seed);

  for (size_t i=0; i<input_var.size(); i++) 
    if (input_var[i] < Vmax && Vmin < input_var[i]) {
      int occ = max(0, min(int((input_var[i]-Vmin)*binSize_inv), nbin));
      if (ran() < fvar_target[occ]/fvar_input[occ])
	m_object.push_back(shared_ptr<Object>(input_catalogue.catalogue_object(i)));
    }
}


// ============================================================================


cbl::catalogue::Catalogue::Catalogue (const Catalogue input_catalogue, const Catalogue target_catalogue, const Var var_name1, const int nbin1, const Var var_name2, const int nbin2, const int seed) 
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
      const int occ1 = max(0, min(int((input_var1[i]-V1min)*binSize1_inv), nbin1));
      const int occ2 = max(0, min(int((input_var2[i]-V2min)*binSize2_inv), nbin2));
      if (ran() < fvars_target[occ1][occ2]/fvars_input[occ1][occ2])
	m_object.push_back(shared_ptr<Object>(input_catalogue.catalogue_object(i)));
    }
  
} 


// ============================================================================


void cbl::catalogue::Catalogue::remove_objects (const std::vector<bool> index)
{
  if (index.size() != m_object.size()) ErrorCBL ("argument size not valid!", "remove_objects", "Catalogue.cpp");

  decltype(m_object) object_temp;
  
  for (size_t ii = 0; ii<index.size(); ii++) 
    if (!index[ii]) object_temp.emplace_back(m_object[ii]);
  
  m_object.swap(object_temp);
}


// ============================================================================


void cbl::catalogue::Catalogue::swap_objects (const int ind1, const int ind2)
{
  shared_ptr<Object> temp = m_object[ind1];
  m_object[ind1] = m_object[ind2];
  m_object[ind2] = temp;
}


// ============================================================================


void cbl::catalogue::Catalogue::sort (const Var var_name, const bool increasing)
{
  coutCBL << "I'm sorting the catalogue..." << endl;
  
  vector<double> variable = var(var_name);
  bool swap = true;
  while (swap) {
    swap = false;
    for (size_t i = 0; i<nObjects()-1; ++i) {
      if (increasing) {
	if (variable[i] > variable[i+1]) {
	  double temp = variable[i];
	  variable[i] = variable[i+1];
	  variable[i+1] = temp;
	  swap_objects(i, i+1);
	  swap = true;
	}
      }
      else {
	if (variable[i] < variable[i+1]) {
	  double temp = variable[i];
	  variable[i] = variable[i+1];
	  variable[i+1] = temp;
	  swap_objects(i, i+1);
	  swap = true;
	}
      }
    }
  }
}


// ============================================================================


void cbl::catalogue::Catalogue::shuffle (const int seed)
{
  coutCBL << "I'm shuffling objects in the catalogue..." << endl;

  /// pseudo-random numbers generator
  std::mt19937_64 generator;
  generator.seed(seed);

  std::shuffle(m_object.begin(), m_object.end(), generator);
}



// ============================================================================


void cbl::catalogue::Catalogue::compute_catalogueProperties (const double boxside)
{
  m_volume = (boxside > 0.) ? pow(boxside, 3.) :
    (cbl::Max(var(Var::_X_)) - cbl::Min(var(Var::_X_)))*
    (cbl::Max(var(Var::_Y_)) - cbl::Min(var(Var::_Y_)))*
    (cbl::Max(var(Var::_Z_)) - cbl::Min(var(Var::_Z_)));
  coutCBL << "Sample volume = " << m_volume << " (Mpc/h)^3" << endl;
  
  m_numdensity = m_object.size()/m_volume;
  coutCBL << "Sample density = " << m_numdensity << " (Mpc/h)^-3" << endl;
  
  m_mps = pow(m_numdensity, -1./3.);
  coutCBL << "Sample mps = " << m_mps << " Mpc/h" << endl; 
}

