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

using namespace std;

/// time_point_type
using time_point_type=std::chrono::time_point < std::chrono::system_clock, std::chrono::milliseconds > ;

using namespace cbl;


// ============================================================================


void cbl::glob::spherical_harmonics_coeff::initialize (const int norder, const int nbins)
{

  m_nbins = nbins;
  m_norder = norder;
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


void cbl::glob::spherical_harmonics_coeff::reset ()
{
  for (int b=0; b<m_nbins; b++)
    for (int k1=0; k1<m_n_sph; k1++)
	m_alm[b][k1] = 0;
}


// ============================================================================


std::vector<std::complex<double>> cbl::glob::spherical_harmonics_coeff::alm (const double xx, const double yy, const double zz)
{
  return spherical_harmonics_array (m_lmax, xx, yy, zz);
}


// ============================================================================


void cbl::glob::spherical_harmonics_coeff::add (const std::vector<std::complex<double>> alm, const double ww, const int bin)
{
  for(int n=0; n<m_n_sph; n++)
    m_alm[bin][n] += ww*alm[n];
}


// ============================================================================


void cbl::glob::spherical_harmonics_coeff::add (const double xx, const double yy, const double zz, const double ww, const int bin)
{
  vector<complex<double>> alm = spherical_harmonics_array (m_lmax, xx, yy, zz);
  for(int n=0; n<m_n_sph; n++)
    m_alm[bin][n] += ww*alm[n];
}


// ============================================================================


double cbl::glob::spherical_harmonics_coeff::power (const int l, const int bin1, const int bin2)
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


vector<double> cbl::glob::count_triplets_SphericalHarmonics (const double r12_min, const double r12_max, const double r13_min, const double r13_max, const int norders, const cbl::catalogue::Catalogue catalogue)
{
  vector<double> zeta_l(norders, 0);

  std::shared_ptr<catalogue::Catalogue> cc(new catalogue::Catalogue(catalogue));

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

	if(rr>0) {

	  if (rr>=r12_min && rr<=r12_max) 
	    alm.add (xx/rr, yy/rr, zz/rr, cat->weight(close_objects[j]), 0);

	  if (rr>=r13_min && rr<=r13_max) 
	    alm.add (xx/rr, yy/rr, zz/rr, cat->weight(close_objects[j]), 1);
	  
	}
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


void cbl::glob::count_triplets_SphericalHarmonics (std::vector<double> &pairs, std::vector<std::vector<std::vector<double>>> &triplets, const double rmin, const double rmax, const int nbins, const int norders, const cbl::catalogue::Catalogue catalogue)
{
  pairs.erase(pairs.begin(), pairs.end());
  triplets.erase(triplets.begin(), triplets.end());

  pairs.resize(nbins, 0);
  triplets.resize(nbins, vector<vector<double>>(nbins, vector<double>(norders, 0)));

  std::shared_ptr<catalogue::Catalogue> cc(new catalogue::Catalogue(catalogue));

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

    vector<double> _pairs(nbins, 0);
    vector<vector<vector<double>>> _triplets(nbins, vector<vector<double>>(nbins, vector<double>(norders, 0)));

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

	  int jbin = max(0, min(int((rr-rmin)*binSize_inv), nbins));

	  _pairs[jbin] += iww*jww;

	  alm.add (xx/rr, yy/rr, zz/rr, jww, jbin);
	}
      }

      for (int b1=0; b1<nbins; b1++)
	for (int b2=0; b2<nbins; b2++)
	  for (int l=0; l<norders; l++)
	    _triplets[b1][b2][l] += iww*alm.power(l, b1, b2);
    }
#pragma omp critical
    {
      for (int b1=0; b1<nbins; b1++){
	pairs[b1] += _pairs[b1];
	for (int b2=0; b2<nbins; b2++)
	  for (int l=0; l<norders; l++)
	    triplets[b1][b2][l] += _triplets[b1][b2][l];
      }
    }

  }
}


// ============================================================================


void cbl::glob::count_triplets_SphericalHarmonics (const double rmin, const double rmax, const int nbins, const int norders, const cbl::catalogue::Catalogue catalogue, const std::string output_dir, const std::string output_file_pairs, const std::string output_file_triplets)
{
  // ----------------------------------------------------------------------
  // ---------------- measure triplet multipoles expansion ----------------
  // ----------------------------------------------------------------------

  coutCBL << endl;
  coutCBL << "Counting triplets" << endl;
  vector<double> pairs;
  vector<vector<vector<double>>> triplets;
  count_triplets_SphericalHarmonics(pairs, triplets, rmin, rmax, nbins, norders, catalogue);
  coutCBL << "Done!" << endl;
  coutCBL << endl;

  const string mkdir = "mkdir -p "+output_dir;
  if(system(mkdir.c_str())) {}

  string file_out = output_dir+output_file_pairs;
  ofstream fout(file_out.c_str());

  for (int b1=0; b1<nbins; b1++)
    fout << b1 << " " << pairs[b1] << endl;
  fout.clear(); fout.close();

  file_out = output_dir+output_file_triplets;
  fout.open(file_out.c_str());

  for (int b1=0; b1<nbins; b1++)
    for (int b2=0; b2<nbins; b2++)
      for (int l=0; l<norders; l++) 
	fout << b1 << " " << b2 << " " << l << " " << setprecision(10) << triplets[b1][b2][l] <<endl;
  fout.clear(); fout.close();
}


// ============================================================================


std::vector<double> cbl::glob::zeta_SphericalHarmonics (const int nbins, const double side_s, const double side_u, const double perc_increase, const int norders, const cbl::catalogue::Catalogue catalogue, const cbl::catalogue::Catalogue random_catalogue, const std::string output_dir, const std::string output_file)
{
  double r12_min = side_s*(1-perc_increase); 
  double r12_max = side_s*(1+perc_increase); 
  double r13_min = side_s*side_u*(1-perc_increase); 
  double r13_max = side_s*side_u*(1+perc_increase);

  return zeta_SphericalHarmonics (nbins, r12_min, r12_max, r13_min, r13_max, norders, catalogue, random_catalogue, output_dir, output_file);
}


// ============================================================================


std::vector<double> cbl::glob::zeta_SphericalHarmonics (const int nbins, const double r12_min, const double r12_max, const double r13_min, const double r13_max, const int norders, const cbl::catalogue::Catalogue catalogue, const cbl::catalogue::Catalogue random_catalogue, const std::string output_dir, const std::string output_file)
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


void cbl::glob::count_triplets_SphericalHarmonics (std::vector<double> &NNN, std::vector<double> &RRR, const double r12_min, const double r12_max, const double r13_min, const double r13_max, const int norders, const cbl::catalogue::Catalogue catalogue)
{
  NNN.erase(NNN.begin(), NNN.end());
  NNN.resize(norders, 0);

  RRR.erase(RRR.begin(), RRR.end());
  RRR.resize(norders, 0);

  std::shared_ptr<catalogue::Catalogue> cc(new catalogue::Catalogue(catalogue));

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

	if (rr>r12_min && rr<r12_max){
	  vector<complex<double>> _alm = sph.alm(xx/rr, yy/rr, zz/rr);
	  sph.add (_alm, jww, 0);
	  if(iww<0 && jww <0)
	    sph_random.add(_alm, -jww, 0);
	}

	if (rr>r13_min && rr<r13_max){
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


std::vector<double> cbl::glob::zeta_SphericalHarmonics_AllInOne (const int nbins, const double side_s, const double side_u, const double perc_increase, const int norders, const cbl::catalogue::Catalogue catalogue, const cbl::catalogue::Catalogue random_catalogue, const std::string output_dir, const std::string output_file, const bool count_triplets, const std::string dir_triplets)
{
  double r12_min = side_s*(1-perc_increase); 
  double r12_max = side_s*(1+perc_increase); 
  double r13_min = side_s*side_u*(1-perc_increase); 
  double r13_max = side_s*side_u*(1+perc_increase);

  return zeta_SphericalHarmonics_AllInOne (nbins, r12_min, r12_max, r13_min, r13_max, norders, catalogue, random_catalogue, output_dir, output_file, count_triplets, dir_triplets);
}


// ============================================================================


std::vector<double> cbl::glob::zeta_SphericalHarmonics_AllInOne (const int nbins, const double r12_min, const double r12_max, const double r13_min, const double r13_max, const int norders, const cbl::catalogue::Catalogue catalogue, const cbl::catalogue::Catalogue random_catalogue, const std::string output_dir, const std::string output_file, const bool count_triplets, const std::string dir_triplets)
{

  // ---------------------------------------------------------------
  // ---------------- construct the mixed catalogue ----------------
  // ---------------------------------------------------------------

  double N_R = double(random_catalogue.weightedN())/catalogue.weightedN();

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

  ofstream fouttemp("mixed_catalogue.dat");

  for(size_t i=0; i<mixed_catalogue.nObjects(); i++)
    fouttemp << mixed_catalogue.xx(i) << " " << mixed_catalogue.yy(i) << " " << mixed_catalogue.zz(i) << " " << mixed_catalogue.weight(i) << endl;
  fouttemp.close();

  // ------------------------------------------------------------------------------------
  // ---------------- measure triplet multipoles expansion for R and D-R ----------------
  // ------------------------------------------------------------------------------------

  vector<double> NNN, RRR;

  if (count_triplets) {
    coutCBL << endl;
    coutCBL << "Counting triplets multipoles NNN, RRR" << endl;
    count_triplets_SphericalHarmonics(NNN, RRR, r12_min, r12_max, r13_min, r13_max, norders, mixed_catalogue);
    coutCBL << "Done!" << endl;
    coutCBL << endl;

    if(dir_triplets!=par::defaultString) {
      string mkdir = "mkdir -p "+dir_triplets;
      if (system(mkdir.c_str())) {}

      string triplet_file = dir_triplets+"triplets.dat";
      ofstream fout(triplet_file.c_str());

      for (int i=0; i<norders; i++) 
	fout << i << " " << NNN[i] << " " << RRR[i] << endl;

      fout.close();

      coutCBL << "I wrote the file " << triplet_file << endl;
      coutCBL << endl;
    }
  }
  else { 
    string triplet_file = dir_triplets+"triplets.dat";
    coutCBL << "I'm reading the file " << triplet_file << endl;

    ifstream fin (triplet_file.c_str());

    for (int i=0; i<norders; i++) {
      int l;
      double _NNN, _RRR;
      fin >> l >> _NNN >> _RRR;
      NNN.push_back(_NNN);
      RRR.push_back(_RRR);
    }

    fin.close();

    coutCBL << "Done" << endl;

    coutCBL << endl;
  }


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
  const string file2 = output_dir+"reconstructed_"+output_file;
  const string file3 = output_dir+"NNN_"+output_file;
  const string file4 = output_dir+"RRR_"+output_file;

  ofstream fout2(file2.c_str());
  ofstream fout3(file3.c_str());
  ofstream fout4(file4.c_str());

  vector<double> zeta(nbins, 0);

  for (int i=0; i<nbins; i++) {

    double mu_min = cos(double(i)/nbins*par::pi) ;
    double mu_max = cos(double(i+1)/nbins*par::pi);

    double angle = (i+0.5)/nbins*par::pi;
    fout2 << setprecision(10) << angle/par::pi;
    fout3 << setprecision(10) << angle/par::pi;
    fout4 << setprecision(10) << angle/par::pi;

    double nnn_theta = 0;
    double rrr_theta = 0;
    for (int l=0; l<norders; l++){
      double fact = cbl::Legendre_polynomial_mu_average (mu_min, mu_max, l);
      nnn_theta += NNN[l]*fact;
      rrr_theta += RRR[l]*fact;
      fout2 << setprecision(10) << " " << nnn_theta/rrr_theta;
      fout3 << setprecision(10) << " " << nnn_theta;
      fout4 << setprecision(10) << " " << rrr_theta;
    }

    fout2 << endl;
    fout3 << endl;
    fout4 << endl;
  }
  fout2.clear(); fout2.close();
  fout3.clear(); fout3.close();
  fout4.clear(); fout4.close();


  ofstream fout(file.c_str());
  for (int l=0; l<norders; l++)
    fout << l << " " << zeta_l[l] << endl;

  coutCBL << "I wrote the file: "<< file << endl;
  fout.clear(); fout.close();

  return zeta;
}


// ============================================================================


std::vector<double> cbl::glob::count_triplets_classic (const double r12_min, const double r12_max, const double r13_min, const double r13_max, const int nbins, const cbl::catalogue::Catalogue catalogue, const cbl::triplets::TripletType tripletType)
{
  double r12 = 0.5*(r12_max-r12_min); 
  double r13 = 0.5*(r13_max+r13_min); 
  const double r12_binSize = r12_max-r12_min;
  const double r13_binSize = r13_max-r13_min;

  shared_ptr<triplets::Triplet> tt = move(triplets::Triplet::Create(tripletType, r12, r12_binSize, r13, r13_binSize, nbins));
  int nn = 0;
  int np =0 ;

  std::shared_ptr<catalogue::Catalogue> cc(new catalogue::Catalogue(catalogue));

  chainmesh::ChainMesh_Catalogue cm;

  cm.set_par(r13_max*0.5, cc, r13_max);

  auto cat = cm.catalogue();

  const int nObjects = catalogue.nObjects();

  // start the multithreading parallelization
#pragma omp parallel num_threads(omp_get_max_threads())
  {
    // internal object used by each thread to handle triplets
    shared_ptr<triplets::Triplet> tt_thread = move(triplets::Triplet::Create(tripletType, r12, r12_binSize, r13, r13_binSize, nbins));
    int _nn =0;
    int _np =0;

    // parallelized loop
#pragma omp for schedule(static, 2)
    for (int i=0; i<nObjects; i++)
    {

      double ixx = cat->xx(i);
      double iyy = cat->yy(i);
      double izz = cat->zz(i);

      vector<int> r12_object, r13_object;

      vector<long> close_objects = cm.close_objects({ixx, iyy, izz}, -1);

      for (size_t j=0; j<close_objects.size(); j++) {
	double jxx = cat->xx(close_objects[j]);
	double jyy = cat->yy(close_objects[j]);
	double jzz = cat->zz(close_objects[j]);

	double xx = jxx-ixx;
	double yy = jyy-iyy;
	double zz = jzz-izz;

	double rr = sqrt(xx*xx+yy*yy+zz*zz);

	if(rr>0) {
	  if (rr>=r12_min && rr<=r12_max) 
	    r12_object.push_back(close_objects[j]);
          
	  if (rr>=r13_min && rr<=r13_max)
	    r13_object.push_back(close_objects[j]);
	}
      }
      _np += (int)r12_object.size();

      for (size_t j=0; j<r12_object.size(); j++)
	for (size_t k=0; k<r13_object.size(); k++) {
	  tt_thread->put(cat->operator[](i), cat->operator[](r12_object[j]), cat->operator[](r13_object[k]));
	  if ((r13_object[k] == r12_object[j]))
	    _nn+=cat->weight(i)*cat->weight(r12_object[j])*cat->weight(r13_object[k]);
	}
    }

    #pragma omp critical
    {
      // sum all the object triplets computed by each thread
      tt->Sum(tt_thread);
      nn += _nn;
      np += _np;
    }


  }

  coutCBL << "Number of pairs: " << np << endl;
  coutCBL << "Number of degenerate triangles: " << nn << endl;

  return tt->TT1D();
}
