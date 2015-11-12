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

#include "Catalogue.h"
using namespace cosmobl;


// ============================================================================


int cosmobl::Catalogue::Nregion ()
{ 
  vector<int> regions; 
  for (int i=0; i<nObjects(); i++) regions.push_back(m_sample[i]->region()); 
  sort(regions.begin(), regions.end());
  vector<int>::iterator it = unique(regions.begin(), regions.end());
  regions.resize(std::distance(regions.begin(), it));
  return regions.size();
}


// ============================================================================


vector<long> cosmobl::Catalogue::get_region_list ()
{ 
  vector<long> regions; 
  for (int i=0; i<nObjects(); i++) regions.push_back(m_sample[i]->region()); 
  sort(regions.begin(), regions.end());
  vector<long>::iterator it = unique(regions.begin(), regions.end());
  regions.resize(std::distance(regions.begin(), it));
  return regions;
}


// ============================================================================


vector<double> cosmobl::Catalogue::var (Var var_name) 
{
  vector<double> vv(m_sample.size(), 0.);
  
  switch (var_name) {

  case Var::_XX_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->xx();
    break;

  case Var::_YY_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->yy();
    break;

  case Var::_ZZ_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->zz();
    break;

  case Var::_RA_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->ra();
    break;

  case Var::_DEC_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->dec();
    break;

  case Var::_REDSHIFT_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->redshift();
    break;

  case Var::_DC_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->dc();
    break;

  case Var::_WEIGHT_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->weight();
    break;

  case Var::_MASS_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->mass();
    break;

  case Var::_MAGNITUDE_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->magnitude();
    break;

  case Var::_RICHNESS_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->richness();
    break;

  case Var::_VX_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->vx();
    break;
  
  case Var::_VY_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->vy();
    break;
  
  case Var::_VZ_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->vz();
    break;

  case Var::_REGION_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->region();
    break; 
  
  case Var::_GENERIC_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->generic();
    break;

  case Var::_RADIUS_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->radius();
    break;

  default:
    ErrorMsg("Error in cosmobl::Catalogue::var of Catalogue.cpp: no such a variable in the list!");
  }
  
  return vv;
}


// ============================================================================


void cosmobl::Catalogue::set_var (Var var_name, vector<double> _var)
{
  if (m_sample.size()!=_var.size()) ErrorMsg ("Error in cosmobl::Catalogue::set_var of Catalogue.cpp: different sizes!");
  
  switch (var_name) {

  case Var::_XX_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_xx(_var[i]);
    break;

  case Var::_YY_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_yy(_var[i]);
    break;

  case Var::_ZZ_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_zz(_var[i]);
    break;

  case Var::_RA_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_ra(_var[i]);
    break;

  case Var::_DEC_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_dec(_var[i]);
    break;

  case Var::_REDSHIFT_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_redshift(_var[i]);
    break;

  case Var::_DC_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_dc(_var[i]);
    break;

  case Var::_WEIGHT_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_weight(_var[i]);
    break;

  case Var::_MASS_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_mass(_var[i]);
    break;

  case Var::_MAGNITUDE_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_magnitude(_var[i]);
    break;

  case Var::_RICHNESS_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_richness(_var[i]);
    break;

  case Var::_VX_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_vx(_var[i]);
    break;
  
  case Var::_VY_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_vy(_var[i]);
    break;
  
  case Var::_VZ_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_vz(_var[i]);
    break;

  case Var::_REGION_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_region(_var[i]);
    break;

  case Var::_GENERIC_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_generic(_var[i]);
    break;

  case Var::_RADIUS_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_radius(_var[i]);
    break;

  default:
    ErrorMsg("Error in cosmobl::Catalogue::set_var of Catalogue.cpp: no such a variable in the list!");
  }

}


// ============================================================================


void cosmobl::Catalogue::MinMax_var (Var var_name, vector<double> &Lim, bool er)
{
  if (er) Lim.erase(Lim.begin(), Lim.end());

  Lim.push_back(Min(var(var_name))); 
  Lim.push_back(Max(var(var_name)));
}


// ============================================================================


void cosmobl::Catalogue::MinMax_var (vector<Var> var_name, vector<vector<double> > &Lim, bool er)
{
  if (er) Lim.erase(Lim.begin(),Lim.end());

  for (unsigned int i=0; i<var_name.size(); i++) 
    Lim.push_back(MinMax_var(var_name[i]));
}


// ============================================================================


vector<double> cosmobl::Catalogue::MinMax_var (Var var_name)
{
  vector<double> Lim; 

  Lim.push_back(Min(var(var_name))); Lim.push_back(Max(var(var_name)));

  return Lim;
}


// ============================================================================


void cosmobl::Catalogue::stats_var (Var var_name, vector<double> &stats)
{
  stats.erase(stats.begin(), stats.end());
  stats.resize(4);
  
  stats[0] = Average(var(var_name)); 
  stats[2] = Sigma(var(var_name));
  
  double first, third;
  Quartile(var(var_name), stats[1], first, third);
  
  stats[3] = third-first;
}


// ============================================================================


void cosmobl::Catalogue::stats_var (vector<Var> var_name, vector< vector<double> > &stats)
{
  stats.erase(stats.begin(),stats.end());

  for (unsigned int i=0; i<var_name.size(); i++) {

    vector<double> stats_temp;
    stats_var(var_name[i],stats_temp);

    stats.push_back(stats_temp);
  }
}


// ============================================================================


void cosmobl::Catalogue::var_distr (Var var_name, vector<double> &_var, vector<double> &dist, int nbin, bool linear, string file_out, double Volume, bool norm, double V1, double V2, bool bin_type, bool convolution, double sigma)
{ 
  distribution(_var, dist, var(var_name), var(Var::_WEIGHT_), nbin, linear, file_out, (norm) ? Volume*weightedN() : Volume, V1, V2, bin_type, convolution, sigma);
}


// ============================================================================


void cosmobl::Catalogue::computeComovingCoordinates (Cosmology &cosm)
{
  double red, xx, yy, zz;

  for (int i=0; i<nObjects(); i++) {

    red = redshift(i);
    m_sample[i]->set_dc(cosm.D_C(red));
    
    cartesian_coord(ra(i), dec(i), dc(i), &xx, &yy, &zz);
    
    m_sample[i]->set_xx(xx); 
    m_sample[i]->set_yy(yy); 
    m_sample[i]->set_zz(zz);
  }
}


// ============================================================================


void cosmobl::Catalogue::computePolarCoordinates ()
{
  double ra, dec, dc;

  for (int i=0; i<nObjects(); i++) {
    polar_coord(xx(i), yy(i), zz(i), &ra, &dec, &dc);
    m_sample[i]->set_ra(ra); 
    m_sample[i]->set_dec(dec); 
    m_sample[i]->set_dc(dc);
  }
}

// ============================================================================


void cosmobl::Catalogue::computePolarCoordinates (Cosmology &cosm, double z1, double z2)
{
  double ra, dec, dc;

  for (int i=0; i<nObjects(); i++) {
    polar_coord(xx(i), yy(i), zz(i), &ra, &dec, &dc);
    m_sample[i]->set_ra(ra); 
    m_sample[i]->set_dec(dec); 
    m_sample[i]->set_dc(dc);
    m_sample[i]->set_redshift(cosm.Redshift(dc, z1, z2));
  }
}

// ============================================================================


void cosmobl::Catalogue::normalizeComovingCoordinates () 
{
  for (int i=0; i<nObjects(); i++) { 
    m_sample[i]->set_xx(xx(i)/dc(i)); 
    m_sample[i]->set_yy(yy(i)/dc(i)); 
    m_sample[i]->set_zz(zz(i)/dc(i));
  }
}


// ============================================================================


void cosmobl::Catalogue::restoreComovingCoordinates ()
{
  for (int i=0; i<nObjects(); i++) {
    m_sample[i]->set_xx(xx(i)*dc(i)); 
    m_sample[i]->set_yy(yy(i)*dc(i)); 
    m_sample[i]->set_zz(zz(i)*dc(i));
  }
}


// ============================================================================


void cosmobl::Catalogue::Order (vector<int> vv) 
{
  int nObj = m_sample.size();

  if (int(vv.size())!=nObj) ErrorMsg("Error in cosmobl::Catalogue::Order!");
 
  vector<shared_ptr<Object>> obj(nObj);

  m_index.resize(nObj);
  
  for (unsigned int i=0; i<vv.size(); i++) {
    m_index[i] = vv[i];
    obj[i] = m_sample[vv[i]];
  }

  m_sample = obj;
}


// ============================================================================


void cosmobl::Catalogue::Order () 
{ 
  int nObj = m_sample.size();
  
  vector<shared_ptr<Object>> obj(nObj);
  
  m_index.resize(nObj);
  
  obj = m_sample;
  
  for (int i=0; i<nObj; i++) 
    m_sample[i] = obj[m_index[i]];
}


// ============================================================================


double cosmobl::Catalogue::weightedN () 
{
  double nn = 0.;
  for (unsigned int i=0; i<m_sample.size(); i++)
    nn += m_sample[i]->weight();
  return nn;
}


// ============================================================================


void cosmobl::Catalogue::write_coords (string &cat) 
{
  ofstream fout(cat.c_str()); checkIO(cat, 0);
 
  for (int i=0; i<nObjects(); i++) 
    fout << xx(i) << " " << yy(i) << " " << zz(i) << endl;

  cout <<"I wrote the file: "<<cat<<endl;
  fout.clear(); fout.close();
}


// ============================================================================


void cosmobl::Catalogue::write_obs_coords (string &cat) 
{
  ofstream fout(cat.c_str()); checkIO(cat, 0);

  for (int i=0; i<nObjects(); i++)
    fout << ra(i) << " " << dec(i) << " " << redshift(i) << " " << dc(i) << endl;
  
  cout <<"I wrote the file: "<<cat<<endl;
  fout.clear(); fout.close();
}


// ============================================================================


shared_ptr<Catalogue> cosmobl::Catalogue::cut (Var var_name, double &down, double &up, bool excl)
{
  vector<shared_ptr<Object> > objects;
  vector<double> vvar = var(var_name);

  for (size_t i=0; i<m_sample.size(); i++)   
    
    if (!excl) { 
      if (down<=vvar[i] && vvar[i]<up)
	objects.push_back(m_sample[i]);
    }
    else {
      if (down>vvar[i] || vvar[i]>up)
	objects.push_back(m_sample[i]);
    }
  
  shared_ptr<Catalogue> cat(new Catalogue{objects});
  return cat;
}


// ============================================================================


double cosmobl::Catalogue::distance (int i, shared_ptr<Object> obj)
{
  return sqrt((m_sample[i]->xx()-obj->xx())*(m_sample[i]->xx()-obj->xx())+
	      (m_sample[i]->yy()-obj->yy())*(m_sample[i]->yy()-obj->yy())+
	      (m_sample[i]->zz()-obj->zz())*(m_sample[i]->zz()-obj->zz()));
}


// ============================================================================


double cosmobl::Catalogue::angsep_xyz (int i, shared_ptr<Object> obj)
{ 
  return 2.*asin(0.5*sqrt((m_sample[i]->xx()-obj->xx())*(m_sample[i]->xx()-obj->xx())+
			  (m_sample[i]->yy()-obj->yy())*(m_sample[i]->yy()-obj->yy())+
			  (m_sample[i]->zz()-obj->zz())*(m_sample[i]->zz()-obj->zz())));
}
    

// ============================================================================


shared_ptr<Catalogue> cosmobl::Catalogue::smooth (double gridsize, vector<Var> vars, int SUB)
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
  
  vector<double> Lim;
  cat->MinMax_var(Var::_XX_, Lim, 0);
  cat->MinMax_var(Var::_YY_, Lim, 0);
  cat->MinMax_var(Var::_ZZ_, Lim, 0);

  int nx = SUB, ny = SUB, nz = SUB;
  
  double Cell_X = (Lim[1]-Lim[0])/nx;
  double Cell_Y = (Lim[3]-Lim[2])/ny;
  double Cell_Z = (Lim[5]-Lim[4])/nz;

  for (int i=0; i<cat->nObjects(); i++) {
    int i1 = min(int((cat->xx(i)-Lim[0])/Cell_X), nx-1);
    int j1 = min(int((cat->yy(i)-Lim[2])/Cell_Y), ny-1);
    int z1 = min(int((cat->zz(i)-Lim[4])/Cell_Z), nz-1);
    int index = z1+nz*(j1+ny*i1);
    cat->object(i)->set_region(index);
  }
  
  vector<long> region_list = cat->get_region_list();
  int nRegions = region_list.size();
  
  vector<shared_ptr<Catalogue>> subSamples(nRegions);
  
  for (int i=0; i<nRegions; i++) {
    double start = region_list[i];
    double stop = start+1;
    subSamples[i] = cut(Var::_REGION_, start, stop);
  }

  
  // -----------------------------------------
  // ----- smooth all the sub-catalogues -----
  // -----------------------------------------

  for (int rr=0; rr<nRegions; rr++) {
    
    vector<double> _xx = subSamples[rr]->var(Var::_XX_), _yy = subSamples[rr]->var(Var::_YY_), _zz = subSamples[rr]->var(Var::_ZZ_);
    
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
	
	shared_ptr<Object> obj{new GenericObject()};
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
