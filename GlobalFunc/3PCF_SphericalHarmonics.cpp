/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file GlobalFunc/3PCF_SphericalHarmonics.cpp
 *
 *  @brief Temporary function to compute the 3pcf following
 *  Slepian, Eisenstein 2015
 *
 *  This file contains temporary function to compute the 3pcf following
 *  Slepian, Eisenstein 2015
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "GlobalFunc.h"
#include <chrono>

/// time_point_type
using time_point_type=std::chrono::time_point < std::chrono::system_clock, std::chrono::milliseconds > ;

using namespace cosmobl;


// ============================================================================


void cosmobl::spherical_harmonics_coeff::initialize (const int _norder, const int _nbins)
{

  m_nbins = _nbins;
  m_norder = _norder;
  m_lmax = m_norder-1;
  m_n_sph = gsl_sf_legendre_array_n(m_lmax);

  m_n_sph_l.resize(m_norder, 0);

  int n=0;
  for (int i=0; i<m_norder; i++) {
    for (int j=0; j<i+1; j++)
      n++;
    m_n_sph_l[i] = n;
  }

  vector<vector<complex<double>>> _alm(m_nbins, vector<complex<double>>(m_n_sph, 0));

  m_alm = _alm;
}


// ============================================================================


void cosmobl::spherical_harmonics_coeff::reset ()
{
  for (int b=0; b<m_nbins; b++)
    for (int k1=0; k1<m_n_sph; k1++)
	m_alm[b][k1] = 0;
}


// ============================================================================


vector<complex<double>> cosmobl::spherical_harmonics_coeff::alm (const double xx, const double yy, const double zz)
{
  return spherical_harmonics_array (m_lmax, xx, yy, zz);
}


// ============================================================================


void cosmobl::spherical_harmonics_coeff::add (const vector<complex<double>> alm, const double ww, const int bin)
{
  for(int n=0; n<m_n_sph; n++)
    m_alm[bin][n] += ww*alm[n];
}


// ============================================================================


void cosmobl::spherical_harmonics_coeff::add (const double xx, const double yy, const double zz, const double ww, const int bin)
{
  vector<complex<double>> alm = spherical_harmonics_array (m_lmax, xx, yy, zz);
  for(int n=0; n<m_n_sph; n++)
    m_alm[bin][n] += ww*alm[n];
}


// ============================================================================


double cosmobl::spherical_harmonics_coeff::power (const int l, const int bin1, const int bin2)
{
  const int min_n = m_n_sph_l[l-1];
  double power = (real(min_n, bin1)*real(min_n, bin2)+imag(min_n, bin1)*imag(min_n, bin2));
  for (int m=1; m<l+1; m++) {
    const int pos = min_n+m;
    power += 2.*(real(pos, bin1)*real(pos, bin2)+imag(pos, bin1)*imag(pos, bin2));
  }

  return 4.*par::pi/(2*l+1)*power;
}


// ============================================================================


vector<double> cosmobl::count_triplets_SphericalHarmonics (const double r12_min, const double r12_max, const double r13_min, const double r13_max, const int norders, const cosmobl::catalogue::Catalogue catalogue)
{
  vector<double> zeta_l(norders, 0);

  shared_ptr<catalogue::Catalogue> cc(new catalogue::Catalogue(catalogue));

  chainmesh::ChainMesh_Catalogue cm;

  cm.set_par(r13_max*0.5, cc, r13_max);

  auto cat = cm.catalogue();

  const int nObjects = catalogue.nObjects();

  // start the multithreading parallelization
#pragma omp parallel num_threads(omp_get_max_threads())
  {
    spherical_harmonics_coeff alm(norders, 2);
    vector<double> _zeta_l(norders, 0);

    // parallelized loop
#pragma omp for schedule(static, 2)
    for (int i=0; i<nObjects; i++)
    {
      alm.reset();

      double ixx = cat->xx(i);
      double iyy = cat->yy(i);
      double izz = cat->zz(i);
      double iww = cat->weight(i);

      vector<vector<double>> r1, r2;

      vector<long> close_objects = cm.close_objects({ixx, iyy, izz}, -1);

      for (size_t j=0; j<close_objects.size(); j++) {
	double jxx = cat->xx(close_objects[j]);
	double jyy = cat->yy(close_objects[j]);
	double jzz = cat->zz(close_objects[j]);

	double xx = jxx-ixx;
	double yy = jyy-iyy;
	double zz = jzz-izz;

	double rr = sqrt(xx*xx+yy*yy+zz*zz);

	if (rr>=r12_min && rr<=r12_max)
	  alm.add (xx/rr, yy/rr, zz/rr, cat->weight(close_objects[j]), 0);

	if (rr>=r13_min && rr<=r13_max)
	  alm.add (xx/rr, yy/rr, zz/rr, cat->weight(close_objects[j]), 1);
      }

      for (int l=0; l<norders; l++)
	  _zeta_l[l] += iww*alm.power(l, 0, 1);

    }

#pragma omp critical
    {
      for (int l=0; l<norders; l++)
	zeta_l[l] += _zeta_l[l];
    }

  }

  return zeta_l;
}


// ============================================================================


vector<vector<vector<double>>> cosmobl::count_triplets_SphericalHarmonics (const double rmin, const double rmax, const int nbins, const int norders, const cosmobl::catalogue::Catalogue catalogue)
{
  vector<vector<vector<double>>> zeta_l(nbins, vector<vector<double>>(nbins, vector<double>(norders, 0)));

  shared_ptr<catalogue::Catalogue> cc(new catalogue::Catalogue(catalogue));

  chainmesh::ChainMesh_Catalogue cm;

  double deltaBin = (rmax-rmin)/double(nbins);
  double binSize_inv = 1./deltaBin;

  cm.set_par(rmax*0.5, cc, rmax*1.1);

  auto cat = cm.catalogue();

  const int nObjects = catalogue.nObjects();

  // start the multithreading parallelization
#pragma omp parallel num_threads(omp_get_max_threads())
  {
    spherical_harmonics_coeff alm(norders, nbins+1);

    vector<vector<vector<double>>> _zeta_l(nbins, vector<vector<double>>(nbins, vector<double>(norders, 0)));

    // parallelized loop
#pragma omp for schedule(static, 2)
    for (int i=0; i<nObjects; i++)
    {
      alm.reset();

      double ixx = cat->xx(i);
      double iyy = cat->yy(i);
      double izz = cat->zz(i);
      double iww = cat->weight(i);

      vector<long> close_objects = cm.close_objects({ixx, iyy, izz}, -1);

      for (size_t j=0; j<close_objects.size(); j++) {

	double jxx = cat->xx(close_objects[j]);
	double jyy = cat->yy(close_objects[j]);
	double jzz = cat->zz(close_objects[j]);
	double jww = cat->weight(close_objects[j]);

	double xx = jxx-ixx;
	double yy = jyy-iyy;
	double zz = jzz-izz;

	double rr = sqrt(xx*xx+yy*yy+zz*zz);

	if (rr>=rmin && rr<=rmax && i!=close_objects[j]) {

	  int jbin = int((rr-rmin)*binSize_inv);

	  alm.add (xx/rr, yy/rr, zz/rr, jww, jbin);
	}
      }

      for (int b1=0; b1<nbins; b1++)
	for (int b2=b1; b2>-1; b2--)
	  for (int l=0; l<norders; l++)
	    _zeta_l[b1][b2][l] += iww*alm.power(l, b1, b2);

    }
#pragma omp critical
    {
      for (int b1=0; b1<nbins; b1++)
	for (int b2=b1; b2>-1; b2--)
	  for (int l=0; l<norders; l++)
	    zeta_l[b1][b2][l] += _zeta_l[b1][b2][l];
    }


  }

  return zeta_l;
}


// ============================================================================


void cosmobl::zeta_SphericalHarmonics (const double rmin, const double rmax, const int nbins , const int norders, const cosmobl::catalogue::Catalogue catalogue, const cosmobl::catalogue::Catalogue random_catalogue)
{
  // ---------------------------------------------------------------
  // ---------------- construct the mixed catalogue ----------------
  // ---------------------------------------------------------------

  double N_R = double(random_catalogue.nObjects())/catalogue.nObjects();

  catalogue::Catalogue mixed_catalogue (catalogue);

  for(size_t i=0; i<random_catalogue.nObjects(); i++){
    double xx = random_catalogue.xx(i);
    double yy = random_catalogue.yy(i);
    double zz = random_catalogue.zz(i);
    double ww = -random_catalogue.weight(i)/N_R;

    auto obj = make_shared<catalogue::Object>( catalogue::Object());
    obj->set_xx(xx);
    obj->set_yy(yy);
    obj->set_zz(zz);
    obj->set_weight(ww);
    mixed_catalogue.add_object(obj);
  }

  // ------------------------------------------------------------------------------------
  // ---------------- measure triplet multipoles expansion for R and D-R ----------------
  // ------------------------------------------------------------------------------------

  coutCBL << endl;
  coutCBL << "Counting triplet multipoles DDD" << endl;
  vector<vector<vector<double>>> _DDD = count_triplets_SphericalHarmonics(rmin, rmax, nbins, norders, catalogue);
  coutCBL << "Done!" << endl;
  coutCBL << endl;
  coutCBL << "Counting triplet multipoles RRR" << endl;
  vector<vector<vector<double>>> _RRR = count_triplets_SphericalHarmonics(rmin, rmax, nbins, norders, random_catalogue);
  coutCBL << "Done!" << endl;
  coutCBL << endl;
  coutCBL << "Counting triplet multipoles NNN" << endl;
  vector<vector<vector<double>>> _NNN = count_triplets_SphericalHarmonics(rmin, rmax, nbins, norders, random_catalogue);
  coutCBL << "Done!" << endl;
  coutCBL << endl;

  
  for(int b1=0; b1<nbins; b1++){
    for(int b2=b1-1; b2>-1; b2--){
      vector<double> RRR = _RRR[b1][b2];
      vector<double> NNN = _NNN[b1][b2];

      for(int l=1; l<norders; l++){
	RRR[l] *= RRR[0];
	NNN[l] *= NNN[0];
      }

      for(int l=0; l<norders; l++){
	RRR[l] *= (2.*l+1)/2;
	NNN[l] *= (2.*l+1)/2;
      }

      // ------------------------------------------------------
      // ---------------- compute the matrix A ----------------
      // ------------------------------------------------------

      vector<double> fl = RRR;
      for(int i=0; i<norders; i++)
	fl[i] = fl[i]/RRR[0];

      vector<vector<double>> A(norders, vector<double>(norders, 0)), A_inverse;

      for(int k=0; k<norders; k++)
	for(int l=0; l<norders; l++)
	  for(int lp=1; lp<norders; lp++)
	    A[k][l] += (2*k+1)*pow(gsl_sf_coupling_3j(2*l, 2*lp, 2*k, 0, 0, 0),2)*fl[lp];

      for(int k=0; k<norders; k++)
	A[k][k] += 1.;

      invert_matrix(A, A_inverse);

      // ----------------------------------------------------
      // ---------------- compute the zeta_l ----------------
      // ----------------------------------------------------

      vector<double> zeta_l(norders, 0);

      for(int l=0; l<norders; l++)
	for(int k=0; k<norders; k++)
	  zeta_l[l] += NNN[k]*A_inverse[l][k]/RRR[0];

      // --------------------------------------------------
      // ---------------- compute the zeta ----------------
      // --------------------------------------------------

      vector<double> zeta(nbins, 0);
      cout << b1 << " " << b2 << " ";

      for (int i=0; i<nbins; i++){
	double angle = (i+0.5)/nbins*par::pi;

	double NNNt =0, RRRt=0;

	for (int l=0; l<norders; l++){
	  zeta[i] += zeta_l[l]*legendre_polynomial(cos(angle), l);
	  NNNt += NNN[l]*legendre_polynomial(cos(angle), l);
	  RRRt += RRR[l]*legendre_polynomial(cos(angle), l);
	  cout << zeta_l[l] << " ";
	}
      }
      cout << endl;
    }
  }

}


// ============================================================================


vector<double> cosmobl::zeta_SphericalHarmonics (const int nbins, const double side_s, const double side_u, const double perc_increase, const int norders, const cosmobl::catalogue::Catalogue catalogue, const cosmobl::catalogue::Catalogue random_catalogue, const string output_dir, const string output_file)
{
  double r12_min = side_s*(1-perc_increase); 
  double r12_max = side_s*(1+perc_increase); 
  double r13_min = side_s*side_u*(1-perc_increase); 
  double r13_max = side_s*side_u*(1+perc_increase);

  return zeta_SphericalHarmonics (nbins, r12_min, r12_max, r13_min, r13_max, norders, catalogue, random_catalogue, output_dir, output_file);
}


// ============================================================================


vector<double> cosmobl::zeta_SphericalHarmonics (const int nbins, const double r12_min, const double r12_max, const double r13_min, const double r13_max, const int norders, const cosmobl::catalogue::Catalogue catalogue, const cosmobl::catalogue::Catalogue random_catalogue, const string output_dir, const string output_file)
{

  // ---------------------------------------------------------------
  // ---------------- construct the mixed catalogue ----------------
  // ---------------------------------------------------------------

  double N_R = double(random_catalogue.nObjects())/catalogue.nObjects();

  catalogue::Catalogue ran_catalogue (random_catalogue);
  catalogue::Catalogue mixed_catalogue (catalogue);

  for (size_t i=0; i<random_catalogue.nObjects(); i++) {
    double xx = random_catalogue.xx(i);
    double yy = random_catalogue.yy(i);
    double zz = random_catalogue.zz(i);
    double ww = -random_catalogue.weight(i)/N_R;
    ran_catalogue.set_var(i, catalogue::Var::_Weight_, -ww);

    auto obj = make_shared<catalogue::Object>( catalogue::Object());
    obj->set_xx(xx);
    obj->set_yy(yy);
    obj->set_zz(zz);
    obj->set_weight(ww);
    mixed_catalogue.add_object(obj);
  }

  // ------------------------------------------------------------------------------------
  // ---------------- measure triplet multipoles expansion for R and D-R ----------------
  // ------------------------------------------------------------------------------------

  coutCBL << endl;
  coutCBL << "Counting triplet multipoles DDD" << endl;
  vector<double> DDD = count_triplets_SphericalHarmonics(r12_min, r12_max, r13_min, r13_max, norders, catalogue);
  coutCBL << "Done!" << endl;
  coutCBL << endl;
  coutCBL << "Counting triplet multipoles RRR" << endl;
  vector<double> RRR = count_triplets_SphericalHarmonics(r12_min, r12_max, r13_min, r13_max, norders, random_catalogue);
  coutCBL << "Done!" << endl;
  coutCBL << endl;
  coutCBL << "Counting triplet multipoles NNN" << endl;
  vector<double> NNN = count_triplets_SphericalHarmonics(r12_min, r12_max, r13_min, r13_max, norders, mixed_catalogue);
  coutCBL << "Done!" << endl;
  coutCBL << endl;

  for (int l=0; l<norders; l++) {
    RRR[l] *= (2.*l+1)/2;
    NNN[l] *= (2.*l+1)/2;
  }

  
  // ------------------------------------------------------
  // ---------------- compute the matrix A ----------------
  // ------------------------------------------------------

  vector<double> fl = RRR;
  for (int i=0; i<norders; i++)
    fl[i] = fl[i]/RRR[0];

  vector<vector<double>> A(norders, vector<double>(norders, 0)), A_inverse;

  for (int k=0; k<norders; k++)
    for (int l=0; l<norders; l++)
      for (int lp=1; lp<norders; lp++)
	A[k][l] += (2*k+1)*pow(gsl_sf_coupling_3j(2*l, 2*lp, 2*k, 0, 0, 0),2)*fl[lp];

  for (int k=0; k<norders; k++)
    A[k][k] += 1.;

  invert_matrix(A, A_inverse);

  
  // ----------------------------------------------------
  // ---------------- compute the zeta_l ----------------
  // ----------------------------------------------------

  vector<double> zeta_l(norders, 0);

  for (int l=0; l<norders; l++)
    for (int k=0; k<norders; k++)
      zeta_l[l] += NNN[k]*A_inverse[l][k]/RRR[0];

  
  // --------------------------------------------------
  // ---------------- compute the zeta ----------------
  // --------------------------------------------------

  string mkdir = "mkdir -p "+output_dir;
  if(system(mkdir.c_str())) {}
  
  const string file = output_dir+output_file;
  ofstream fout(file.c_str());
  vector<double> zeta(nbins, 0);

  for (int i=0; i<nbins; i++) {
    
    double angle = (i+0.5)/nbins*par::pi;
    fout << setprecision(10) << angle/par::pi;
    for (int l=0; l<norders; l++){
      zeta[i] += zeta_l[l]*legendre_polynomial(cos(angle), l);
      fout << setprecision(10) << " " << zeta[i];
    }

    fout << endl;
  }
  fout.clear(); fout.close();

  coutCBL << "I wrote the file: "<< file << endl;

  return zeta;
}


// ============================================================================


void cosmobl::count_triplets_SphericalHarmonics (vector<double> &NNN, vector<double> &RRR, const double r12_min, const double r12_max, const double r13_min, const double r13_max, const int norders, const cosmobl::catalogue::Catalogue catalogue)
{
  NNN.erase(NNN.begin(), NNN.end());
  NNN.resize(norders, 0);

  RRR.erase(RRR.begin(), RRR.end());
  RRR.resize(norders, 0);

  shared_ptr<catalogue::Catalogue> cc(new catalogue::Catalogue(catalogue));

  chainmesh::ChainMesh_Catalogue cm;

  cm.set_par(r13_max*0.5, cc, r13_max);

  auto cat = cm.catalogue();

  const int nObjects = catalogue.nObjects();

  // start the multithreading parallelization
#pragma omp parallel num_threads(omp_get_max_threads())
  {
    spherical_harmonics_coeff sph(norders, 2);
    spherical_harmonics_coeff sph_random(norders, 2);
    vector<double> _zeta(norders, 0);
    vector<double> _zeta_random(norders, 0);

    // parallelized loop
#pragma omp for schedule(static, 2)
    for (int i=0; i<nObjects; i++)
    {
      sph.reset();

      double ixx = cat->xx(i);
      double iyy = cat->yy(i);
      double izz = cat->zz(i);
      double iww = cat->weight(i);

      if(iww<0)
	sph_random.reset();

      vector<long> close_objects = cm.close_objects({ixx, iyy, izz}, -1);

      for (size_t j=0; j<close_objects.size(); j++) {
	double jxx = cat->xx(close_objects[j]);
	double jyy = cat->yy(close_objects[j]);
	double jzz = cat->zz(close_objects[j]);
	double jww = cat->weight(close_objects[j]);

	double xx = jxx-ixx;
	double yy = jyy-iyy;
	double zz = jzz-izz;

	double rr = sqrt(xx*xx+yy*yy+zz*zz);

	if (rr>=r12_min && rr<=r12_max){
	  vector<complex<double>> _alm = sph.alm(xx/rr, yy/rr, zz/rr);
	  sph.add (_alm, jww, 0);
	  if(iww<0 && jww <0)
	    sph_random.add(_alm, -jww, 0);
	}

	if (rr>=r13_min && rr<=r13_max){
	  vector<complex<double>> _alm = sph.alm(xx/rr, yy/rr, zz/rr);
	  sph.add (_alm, jww, 1);
	  if(iww<0 && jww <0)
	    sph_random.add(_alm, -jww, 1);
	}
      }

      for (int l=0; l<norders; l++){
	  _zeta[l] += iww*sph.power(l, 0, 1);
	  if(iww<0)
	    _zeta_random[l] -= iww*sph_random.power(l, 0, 1);
      }

    }

#pragma omp critical
    {
      for (int l=0; l<norders; l++){
	NNN[l] += _zeta[l];
	RRR[l] += _zeta_random[l];
      }
    }

  }
}


// ============================================================================


vector<double> cosmobl::zeta_SphericalHarmonics_AllInOne (const int nbins, const double side_s, const double side_u, const double perc_increase, const int norders, const cosmobl::catalogue::Catalogue catalogue, const cosmobl::catalogue::Catalogue random_catalogue, const string output_dir, const string output_file)
{
  double r12_min = side_s*(1-perc_increase); 
  double r12_max = side_s*(1+perc_increase); 
  double r13_min = side_s*side_u*(1-perc_increase); 
  double r13_max = side_s*side_u*(1+perc_increase);

  return zeta_SphericalHarmonics_AllInOne (nbins, r12_min, r12_max, r13_min, r13_max, norders, catalogue, random_catalogue, output_dir, output_file);
}


// ============================================================================


vector<double> cosmobl::zeta_SphericalHarmonics_AllInOne (const int nbins, const double r12_min, const double r12_max, const double r13_min, const double r13_max, const int norders, const cosmobl::catalogue::Catalogue catalogue, const cosmobl::catalogue::Catalogue random_catalogue, const string output_dir, const string output_file)
{

  // ---------------------------------------------------------------
  // ---------------- construct the mixed catalogue ----------------
  // ---------------------------------------------------------------

  double N_R = double(random_catalogue.nObjects())/catalogue.nObjects();

  catalogue::Catalogue ran_catalogue (random_catalogue);
  catalogue::Catalogue mixed_catalogue (catalogue);

  for (size_t i=0; i<random_catalogue.nObjects(); i++) {
    double xx = random_catalogue.xx(i);
    double yy = random_catalogue.yy(i);
    double zz = random_catalogue.zz(i);
    double ww = -random_catalogue.weight(i)/N_R;
    ran_catalogue.set_var(i, catalogue::Var::_Weight_, -ww);

    auto obj = make_shared<catalogue::Object>( catalogue::Object());
    obj->set_xx(xx);
    obj->set_yy(yy);
    obj->set_zz(zz);
    obj->set_weight(ww);
    mixed_catalogue.add_object(obj);
  }

  // ------------------------------------------------------------------------------------
  // ---------------- measure triplet multipoles expansion for R and D-R ----------------
  // ------------------------------------------------------------------------------------

  coutCBL << endl;
  coutCBL << "Counting triplets multipoles NNN, RRR" << endl;
  vector<double> NNN, RRR;
  count_triplets_SphericalHarmonics(NNN, RRR, r12_min, r12_max, r13_min, r13_max, norders, mixed_catalogue);
  coutCBL << "Done!" << endl;
  coutCBL << endl;

  for (int l=0; l<norders; l++) {
    RRR[l] *= (2.*l+1)/2;
    NNN[l] *= (2.*l+1)/2;
  }

  
  // ------------------------------------------------------
  // ---------------- compute the matrix A ----------------
  // ------------------------------------------------------

  vector<double> fl = RRR;
  for (int i=0; i<norders; i++)
    fl[i] = fl[i]/RRR[0];

  vector<vector<double>> A(norders, vector<double>(norders, 0)), A_inverse;

  for (int k=0; k<norders; k++)
    for (int l=0; l<norders; l++)
      for (int lp=1; lp<norders; lp++)
	A[k][l] += (2*k+1)*pow(gsl_sf_coupling_3j(2*l, 2*lp, 2*k, 0, 0, 0),2)*fl[lp];

  for (int k=0; k<norders; k++)
    A[k][k] += 1.;

  invert_matrix(A, A_inverse);

  
  // ----------------------------------------------------
  // ---------------- compute the zeta_l ----------------
  // ----------------------------------------------------

  vector<double> zeta_l(norders, 0);

  for (int l=0; l<norders; l++)
    for (int k=0; k<norders; k++)
      zeta_l[l] += NNN[k]*A_inverse[l][k]/RRR[0];

  
  // --------------------------------------------------
  // ---------------- compute the zeta ----------------
  // --------------------------------------------------

  string mkdir = "mkdir -p "+output_dir;
  if(system(mkdir.c_str())) {}
  
  const string file = output_dir+output_file;
  ofstream fout(file.c_str());
  vector<double> zeta(nbins, 0);

  for (int i=0; i<nbins; i++) {
    
    double angle = (i+0.5)/nbins*par::pi;
    fout << setprecision(10) << angle/par::pi;
    for (int l=0; l<norders; l++){
      zeta[i] += zeta_l[l]*legendre_polynomial(cos(angle), l);
      fout << setprecision(10) << " " << zeta[i];
    }

    fout << endl;
  }
  fout.clear(); fout.close();

  coutCBL << "I wrote the file: "<< file << endl;

  return zeta;
}
