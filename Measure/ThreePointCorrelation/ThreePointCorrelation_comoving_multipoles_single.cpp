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
 *  @file CosmoBolognaLib/Measure/ThreePointCorrelation/ThreePointCorrelation_comoving_multipoles_single.cpp
 *
 *  @brief Methods of the class
 *  ThreePointCorrelation_comoving_multipoles_single used to compute the 
 *  the multipoles of the three-point correlation function for a single configuration
 *
 *  This file contains the implementation of the methods of the class
 *  ThreePointCorrelation_comoving_multipoles_single used to compute the 
 *  the multipoles of the three-point correlation function for a single configuration
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "ChainMesh_Catalogue.h"
#include "SphericalHarmonics_Coefficients.h"
#include "ThreePointCorrelation_comoving_multipoles_single.h"

using namespace std;

using namespace cbl;

using namespace catalogue;
using namespace measure;
using namespace threept;
using namespace glob;



// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_single::m_count_triplets (std::vector<double> &NNN, std::vector<double> &RRR, const double r12_min, const double r12_max, const double r13_min, const double r13_max, const int norders, const catalogue::Catalogue& catalogue) const
{ 
  // timer 
  time_t start; time(&start);
  int dp = cout.precision();
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);

  const int nObjects = catalogue.nObjects();

  /// Erase pair and triplet vectors
  double NN, RR;

  NNN.erase(NNN.begin(), NNN.end());
  NNN.resize(norders, 0);

  RRR.erase(RRR.begin(), RRR.end());
  RRR.resize(norders, 0);

  auto cc = make_shared<Catalogue>(catalogue);
  chainmesh::ChainMesh_Catalogue cm;

  cm.set_par(r13_max*0.5, cc, r13_max);

  auto cat = cm.catalogue();

  // thread number
  int tid = 0;

  // start the multithreading parallelization
#pragma omp parallel num_threads(omp_get_max_threads()) private(tid)
  {
    tid = omp_get_thread_num();

    SphericalHarmonics_Coefficients sph(norders, 2);
    SphericalHarmonics_Coefficients sph_random(norders, 2);
    vector<double> _zeta(norders, 0);
    vector<double> _zeta_random(norders, 0);
    double _nn=0;
    double _rr=0;

    // parallelized loop
#pragma omp for schedule(static, 2)
    for (int i=0; i<nObjects; i++)
      {
	sph.reset();

	double ixx = cat->xx(i);
	double iyy = cat->yy(i);
	double izz = cat->zz(i);
	double iww = cat->weight(i);

	if(iww<0) // iww<0 RANDOM OBJECT!
	  sph_random.reset();

	vector<long> close_objects = cm.close_objects({ixx, iyy, izz}, -1);

	for (size_t j=0; j<close_objects.size(); j++) {

	  int inCommon  = 0;
	  int j_idx = close_objects[j];

	  double jxx = cat->xx(j_idx);
	  double jyy = cat->yy(j_idx);
	  double jzz = cat->zz(j_idx);
	  double jww = cat->weight(j_idx);

	  double xx = jxx-ixx;
	  double yy = jyy-iyy;
	  double zz = jzz-izz;

	  double rr = sqrt(xx*xx+yy*yy+zz*zz);

	  if (rr>r12_min && rr<r12_max){
	    vector<complex<double>> _alm = sph.alm(xx/rr, yy/rr, zz/rr);
	    sph.add (_alm, jww, 0);
	    if(iww<0 && jww <0)
	      sph_random.add(_alm, jww, 0);
	    inCommon += 1;
	  }

	  if (rr>r13_min && rr<r13_max){
	    vector<complex<double>> _alm = sph.alm(xx/rr, yy/rr, zz/rr);
	    sph.add (_alm, jww, 1);
	    if(iww<0 && jww <0)
	      sph_random.add(_alm, jww, 1);
	    inCommon += 1;
	  }

	  /// Check if a pair is in both bins
	  if (inCommon>1) {
	    _nn += iww*jww*jww;
	    if(iww<0 && jww <0)
	      _rr -= iww*jww*jww;
	  }
	}

	for (int l=0; l<norders; l++){
	  _zeta[l] += iww*sph.power(l, 0, 1);
	  if(iww<0)
	    _zeta_random[l] -= iww*sph_random.power(l, 0, 1);
	}

	// estimate the computational time and update the time count
	if (i==int(nObjects*0.25)) coutCBL << ".............25% completed" << endl;
	if (i==int(nObjects*0.5)) coutCBL << ".............50% completed" << endl;
	if (i==int(nObjects*0.75)) coutCBL << ".............75% completed"<< endl;   
      }

#pragma omp critical
    {
      NN += _nn;
      RR += _rr;
      for (int l=0; l<norders; l++){
	NNN[l] += _zeta[l]; 
	RRR[l] += _zeta_random[l];
      }
    }

  }

  for (int l=0; l<norders; l++){
    NNN[l] -= NN; 
    RRR[l] -= RR;
  }

  // show the time spent by the method
  time_t end; time (&end);
  double diff = difftime(end, start);
  if (tid==0) {
    if (diff<60) coutCBL << "   time spent to count the triplets: " << diff << " seconds" << endl ;
    else if (diff<3600) coutCBL << "   time spent to count the triplets: " << diff/60 << " minutes" << endl ;
    else coutCBL << "   time spent to count the triplets: " << diff/3600 << " hours" << endl ;
  }
  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(dp);

}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_single::m_write_triplets (const std::vector<double> TL, const std::string dir, const std::string file) const 
{
  string mkdir = "mkdir -p "+dir;
  if(system(mkdir.c_str())) {}

  string outFile = dir+file;
  ofstream fout (outFile.c_str());

  fout << "# r12   r13 degree triplets_l" << endl;

  for (size_t i=0; i<m_nOrders; i++)
    fout << setprecision(10) << (m_r12Max+m_r12Min)*0.5 << " " << (m_r13Max+m_r13Min)*0.5  << " "<<  i << " " << TL[i] << endl;

  fout.clear(); fout.close(); 

  coutCBL << "I wrote the file " << outFile << endl;
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_single::m_read_triplets (std::vector<double> &TL, const std::vector<std::string> dir, const std::string file) 
{
  TL.erase(TL.begin(), TL.end());

  string outFile = dir[0]+file;
  ifstream fin (outFile.c_str());

  string line;
  // skip header
  getline(fin, line);

  int ell;
  double TT, r12, r13;

  vector<double> triplets_ell;

  while(getline(fin, line))
    {
      stringstream ss(line);
      ss >> r12 >> r13 >> ell >> TT;

      if (r12!=(m_r12Max+m_r12Min)*0.5 || r13!=(m_r13Max+m_r13Min)*0.5)
	ErrorCBL("reading incompatible triplet file: check your inputs!", "m_read_triplets", "ThreePointCorrelation_comoving_multipoles_single.cpp");

      triplets_ell.push_back(TT);
    } 

  fin.clear(); fin.close(); 

  if (triplets_ell.size() < m_nOrders) {
    ErrorCBL("the size of input triplets is smaller than m_nOrders: check your inputs!", "m_read_triplets", "ThreePointCorrelation_comoving_multipoles_single.cpp");
  }
  else if (triplets_ell.size() < m_nOrders) {
    WarningMsgCBL("the size of input triplets is larger than m_nOrders; you can ignore this warning if it's what you want", "m_read_triplets", "ThreePointCorrelation_comoving_multipoles_single.cpp");
    for (size_t i=0; i<m_nOrders; i++)
      TL[i] = triplets_ell[i];
  }
  else {
    TL = triplets_ell;
  }

}


// ============================================================================


cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_single::ThreePointCorrelation_comoving_multipoles_single (const cbl::catalogue::Catalogue catalogue, const cbl::catalogue::Catalogue random_catalogue, const double r12Min, const double r12Max, const double r13Min, const double r13Max, const size_t nOrders, const double split, const size_t seed)
{
  set_parameters(r12Min, r12Max, r13Min, r13Max, nOrders);
  set_catalogues(catalogue, random_catalogue, split, seed);
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_single::set_parameters (const double r12Min, const double r12Max, const double r13Min, const double r13Max, const size_t nOrders)
{
  m_r12Min = r12Min;
  m_r12Max = r12Max;
  m_r13Min = r13Min;
  m_r13Max = r13Max;
  m_nOrders = nOrders;
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_single::measure (const ErrorType errorType, const std::string dir_output_triplets, const std::vector<std::string> dir_input_triplets, const int nResamplings, const bool count_triplets, const bool tcount, const int seed)
{
  coutCBL << "I'm computing the three-point correlation multipoles..." << endl;

  (void)errorType; 
  (void)nResamplings;
  (void)tcount;
  (void)seed;

  if (count_triplets) { // count triplets

    m_nnn.resize(m_nOrders, 0);
    m_rrr.resize(m_nOrders, 0);

    size_t start, end;
    for (size_t i=0; i<m_nSplit; i++) {

      if (m_nSplit>1)
	coutCBL << "Computing the three-point correlation multipoles, " << i+1 << "-th splitting" << endl; 
      
      start = (m_splitFactor<0) ? 0 :  m_data->nObjects()*(i*m_splitFactor);
      end = (m_splitFactor<0) ? m_random->nObjects() : m_data->nObjects()*((i+1)*m_splitFactor);

      auto joined = m_join_catalogues(*m_data, *m_random, start, end);
      
      vector<double> nnn, rrr;
      m_count_triplets(nnn, rrr, m_r12Min, m_r12Max, m_r13Min, m_r13Max, m_nOrders, *joined);

      for (size_t j=0; j<nnn.size(); j++) {
	m_nnn[j] += nnn[j];
	m_rrr[j] += rrr[j];
      }
    }

    for (size_t i=0; i<m_nnn.size(); i++) {
      m_nnn[i] /= m_nSplit;
      m_rrr[i] /= m_nSplit;
    }

    m_write_triplets(m_nnn, dir_output_triplets, "nnn_legendre.dat");
    m_write_triplets(m_rrr, dir_output_triplets, "rrr_legendre.dat");
  }
  else { //Read triplets from file
    m_read_triplets(m_nnn, dir_input_triplets, "nnn_legendre.dat");
    m_read_triplets(m_rrr, dir_input_triplets, "rrr_legendre.dat");
  }

  /// compute zeta

  m_zeta = m_SzapudiSzalay_multipoles (m_nnn, m_rrr);

  coutCBL << "Done!" << endl;
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_single::write (const std::string dir, const std::string file) const 
{
  string mkdir = "mkdir -p "+dir;
  if(system(mkdir.c_str())) {}

  string outFile = dir+file;
  ofstream fout (outFile.c_str());

  fout << "# r12  r13  degree  zeta_l" << endl;

  for (size_t i=0; i<m_nOrders; i++)
    fout << setprecision(10) <<  (m_r12Max+m_r12Min)*0.5 << " " << (m_r13Max+m_r13Min)*0.5 << " " << i << " " << m_zeta[i] <<  endl;

  fout.clear(); fout.close(); 

  coutCBL << "I wrote the file " << outFile << endl;
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_single::resum (const std::string dir, const std::string file, const cbl::triplets::TripletType tripletType, const int nBins, const bool bin) const
{
  std::function<double(const double, const double, const int)> legFunction;

  double min, max;
  string name;

  switch (tripletType) {

  case triplets::TripletType::_comoving_theta_: 
    {
      min = 0;
      max = par::pi;
      name = "theta";

      if (bin)
	legFunction = [&] (const double min, const double max, const int ell) {return cbl::Legendre_polynomial_theta_average(min, max, ell); };
      else
	legFunction = [&] (const double min, const double max, const int ell) {return cbl::legendre_polynomial(cos(0.5*(min+max)), ell); };

      break;
    }

  case triplets::TripletType::_comoving_costheta_:
    {
      min = -1;
      max = 1.;
      name = "costheta";

      if (bin)
	legFunction = [&] (const double min, const double max, const int ell) {return cbl::Legendre_polynomial_mu_average(min, max, ell); };
      else
	legFunction = [&] (const double min, const double max, const int ell) {return cbl::legendre_polynomial(0.5*(min+max), ell); };

      break;
    }

  case triplets::TripletType::_comoving_side_:
    {
      if (!bin) {
	WarningMsgCBL("Cannot use bin=false: Legendre polynomial averaging will be forced!", "resum", "ThreePointCorrelation_comoving_multipoles_single.cpp");
      }

      name = "r23";
      min = (m_r12Min != m_r13Min) ? m_r13Min-m_r12Max : 0.;
      max = m_r13Max+m_r12Max;
      legFunction = [&] (const double min, const double max, const int ell) { return 
									      Legendre_polynomial_triangles_average (m_r12Min, m_r12Max, m_r13Min, m_r13Max, min, max, ell, 1.e-3, 0, 1000); };
      break;
    }

  default:
    {
      ErrorCBL("wrong tripletType in input!", "resum", "ThreePointCorrelation_comoving_multipoles_single.cpp");
      break;
    }
  }

  double delta = (max-min)/nBins;

  string outFile = dir+file;

  ofstream fout(outFile.c_str());

  fout << "# r12  r13  " << name << "  zeta1  zeta2  nnn  rrr" << endl;

  for (int i=0; i<nBins; i++) {
    double _min = min+i*delta;
    double _max = min+(i+1)*delta;

    double nnn=0., rrr=0., zeta=0.;

    for (size_t ell=0; ell<m_nOrders; ell++) {
      double legVal = legFunction(_min, _max, ell);
      double norm = legFunction(_min, _max, 0);

      nnn += legVal*m_nnn[ell]*0.5*(2.*ell+1)/norm;
      rrr += legVal*m_rrr[ell]*0.5*(2.*ell+1)/norm;
      zeta += legVal*m_zeta[ell]/norm;
    }

    fout << 0.5*(m_r12Max+m_r12Min) << " " << 0.5*(m_r13Max+m_r13Min) << " " << (_max+_min)*0.5 << " " << zeta << " " << nnn/rrr << " " << nnn << " " << rrr << endl;
  }

  fout.clear(); fout.close();

  coutCBL << "I wrote the file " << outFile << endl;
}
