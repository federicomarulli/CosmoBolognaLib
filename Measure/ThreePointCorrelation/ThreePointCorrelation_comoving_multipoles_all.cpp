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
 *  @file CosmoBolognaLib/Measure/ThreePointCorrelation/ThreePointCorrelation_comoving_multipoles_all.cpp
 *
 *  @brief Methods of the class
 *  ThreePointCorrelation_comoving_multipoles_all used to compute the 
 *  the multipoles of the three-point correlation function for all configurations
 *
 *  This file contains the implementation of the methods of the class
 *  ThreePointCorrelation_comoving_multipoles_all used to compute the 
 *  the multipoles of the three-point correlation function for all configurations
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */


#include "ChainMesh_Catalogue.h"
#include "SphericalHarmonics_Coefficients.h"
#include "ThreePointCorrelation_comoving_multipoles_all.h"

using namespace std;

using namespace cbl;

using namespace catalogue;
using namespace measure;
using namespace threept;
using namespace glob;



// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_all::m_count_triplets (std::vector<double> &NNN, std::vector<double> &RRR, const double rmin, const double rmax, const int nbins, const int norders, const catalogue::Catalogue& catalogue) const
{
  // timer 
  time_t start; time(&start);
  int dp = cout.precision();
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);

  /// Erase pair and triplet vectors
  vector<double> NN(nbins, 0);
  NNN.erase(NNN.begin(), NNN.end());
  NNN.resize(nbins*nbins*norders);

  vector<double> RR(nbins, 0);
  RRR.erase(RRR.begin(), RRR.end());
  RRR.resize(nbins*nbins*norders);

  auto cc = make_shared<Catalogue>(catalogue);

  chainmesh::ChainMesh_Catalogue cm;

  double deltaBin = (rmax-rmin)/double(nbins);
  double binSize_inv = 1./deltaBin;

  cm.set_par(rmax*0.5, cc, rmax*1.1);

  auto cat = cm.catalogue();

  const int nObjects = catalogue.nObjects();

  // thread number
  int tid = 0;

  // start the multithreading parallelization
#pragma omp parallel num_threads(omp_get_max_threads())  private(tid)
  {
    tid = omp_get_thread_num();

    SphericalHarmonics_Coefficients alm_n(norders, nbins+1);

    vector<double> _NN(nbins+1, 0);
    vector<double> _NNN(nbins*nbins*norders);

    SphericalHarmonics_Coefficients alm_r(norders, nbins+1);
    vector<double> _RR(nbins+1, 0);
    vector<double> _RRR(nbins*nbins*norders);

    int index;

    // parallelized loop
#pragma omp for schedule(static, 2)
    for (int i=0; i<nObjects; i++)
      {

	alm_n.reset();

	double ixx = cat->xx(i);
	double iyy = cat->yy(i);
	double izz = cat->zz(i);
	double iww = cat->weight(i);

	if(iww<0)
	  alm_r.reset();

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

	    vector<complex<double>> _alm = alm_n.alm(xx/rr, yy/rr, zz/rr);

	    int jbin = max(0, min(int((rr-rmin)*binSize_inv), nbins));

	    _NN[jbin] += iww*jww*jww;

	    alm_n.add (_alm, jww, jbin);
	    if (iww<0 && jww <0) {
	      alm_r.add(_alm, jww, jbin);
	      _RR[jbin] -= iww*jww*jww;
	    }
	  }
	}

	for (int b1=0; b1<nbins; b1++)
	  for (int b2=0; b2<nbins; b2++)
	    for (int l=0; l<norders; l++) {
	      index = l+b2*norders+b1*norders*nbins;
	      _NNN[index] += iww*alm_n.power(l, b1, b2);
	      if (iww<0)
		_RRR[index] -= iww*alm_r.power(l, b1, b2);
	    }

	// estimate the computational time and update the time count
	if (i==int(nObjects*0.25)) coutCBL << ".............25% completed" << endl;
	if (i==int(nObjects*0.5)) coutCBL << ".............50% completed" << endl;
	if (i==int(nObjects*0.75)) coutCBL << ".............75% completed"<< endl;   
      }
#pragma omp critical
    {
      for (int b1=0; b1<nbins; b1++){
	NN[b1] += _NN[b1];
	RR[b1] += _RR[b1];
	for (int b2=0; b2<nbins; b2++)
	  for (int l=0; l<norders; l++) {
	    int index = l+b2*norders+b1*norders*nbins;
	    NNN[index] += _NNN[index];
	    RRR[index] += _RRR[index];
	  }
      }
    }
  }

  // isosceles correction
  
  for (int b1=0; b1<nbins; b1++){
    for (int l=0; l<norders; l++) {
      int index = l+b1*norders+b1*norders*nbins;
      NNN[index] -= NN[b1];
      RRR[index] -= RR[b1];
    }
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


void cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_all::m_write_triplets (const std::vector<double> TL, const std::string dir, const std::string file) const 
{
  string mkdir = "mkdir -p "+dir;
  if(system(mkdir.c_str())) {}
  string outFile = dir+file;
  ofstream fout (outFile.c_str());

  fout << "#r12 r13 degree triplets_l index" << endl;

  for (size_t i=0; i<m_nBins; i++) 
    for (size_t j=0; j<m_nBins; j++) 
      for (size_t ell=0; ell<m_nOrders; ell++) {
	int index = ell+j*m_nOrders+i*m_nOrders*m_nBins;
	fout << setprecision(10) <<  m_rMin+(i+0.5)*m_binSize << " " << m_rMin+(j+0.5)*m_binSize << " " << ell << " " << TL[index] << " " << index << endl;
      }

  fout.clear(); fout.close(); 

  coutCBL << "I wrote the file " << outFile << endl;
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_all::m_read_triplets (std::vector<double> &TL, const std::vector<std::string> dir, const std::string file) 
{
  TL.erase(TL.begin(), TL.end());

  string outFile = dir[0]+file;
  ifstream fin (outFile.c_str());

  string line;
  // skip header
  getline(fin, line);

  int ell, index;
  double TT, r12, r13;

  vector<double> triplets_ell;

  while(getline(fin, line))
    {
      stringstream ss(line);
      ss >> r12 >> r13 >> ell >> TT >> index;

      // Add a check here... to think about
      TL.push_back(TT);
    } 

  if (TL.size() != m_nBins*m_nBins*m_nOrders)
    ErrorCBL("wrong total number of triplets, check your input!", "m_read_tripelts ", "ThreePointCorrelation_comoving_multipoles_all.cpp");
}



// ============================================================================


cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_all::ThreePointCorrelation_comoving_multipoles_all (const cbl::catalogue::Catalogue catalogue, const cbl::catalogue::Catalogue random_catalogue, const double rMin, const double rMax, const double binSize, const size_t nOrders, const double split, const size_t seed)
{
  coutCBL << "Setting parameters..." << endl;
  set_parameters(rMin, rMax, binSize, nOrders);
  set_catalogues(catalogue, random_catalogue, split, seed);
  coutCBL << "Done!" << endl;
}


// ============================================================================


cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_all::ThreePointCorrelation_comoving_multipoles_all (const ThreePointCorrelation_comoving_multipoles_all &threept, const double newBinSize)
{
  double rMin = threept.m_rMin;
  double rMax = threept.m_rMax;
  double binSize = threept.m_binSize;
  size_t nOrders = threept.m_nOrders;
  size_t nBins1 = threept.m_nBins;
  (void)nOrders;
  (void)nBins1;

  // First check if newBinSize is a multiple of original binSize. Throw an error if it's not. 

  double integerPart, fracPart;
  double ratio = static_cast<double>(newBinSize/binSize);
  std::modf(ratio, &integerPart);
  fracPart = ratio-integerPart;

  if (fracPart!=0)
    ErrorCBL("New bin size must be a multiple of the input one", "ThreePointCorrelation_comoving_multipoles_all", "ThreePointCorrelation_comoving_multipoles_all.cpp");
  
  // First check if (rMax-rMin)/newBinSize is an integer. 

  double nBins2;
  int min_bin = 0;

  fracPart = 0.1;

  while (fracPart!=0) {
    nBins2 = static_cast<double> ((rMax-(rMin+min_bin*binSize))/newBinSize);

    std::modf(nBins2, &integerPart);
    fracPart = nBins2-integerPart;

    nBins2 = static_cast<size_t> (integerPart); 
    min_bin += 1;
  }
  min_bin -= 1;
  coutCBL << min_bin << " "  << nBins2 << " " << integerPart << " " << fracPart << " " << rMin << " " << rMax <<  endl;

  if (min_bin!=0) {
    WarningMsgCBL("Minimum separation has been set to "+conv(rMin+min_bin*binSize, par::fDP2)+", while the input one is "+conv(rMin, par::fDP2)+".", "ThreePointCorrelation_comoving_multipoles_all", "ThreePointCorrelation_comoving_multipoles_all.cpp");
  }

  // set the parameters
  set_parameters(rMin+min_bin*binSize, rMax, newBinSize, nOrders);

  /// copy catalogue-related variables
  m_data = threept.m_data;
  m_random = threept.m_random;
  m_joined = threept.m_joined;
  m_splitFactor = threept.m_splitFactor;
  m_nSplit = threept.m_nSplit;
    
  // rebin triplets

  m_nnn.erase(m_nnn.begin(), m_nnn.end());
  m_nnn.resize(m_nBins*m_nBins*m_nOrders);

  m_rrr.erase(m_rrr.begin(), m_rrr.end());
  m_rrr.resize(m_nBins*m_nBins*m_nOrders);

  for (size_t i=min_bin; i<nBins1; i++) {
    for (size_t j=min_bin; j<nBins1; j++) {
      int b1 = static_cast<int>((rMin+i*binSize-m_rMin)/m_binSize);
      int b2 = static_cast<int>((rMin+j*binSize-m_rMin)/m_binSize);

      for (size_t ell=0; ell<m_nOrders; ell++) {
	int old_pos = ell+j*m_nOrders+i*m_nOrders*nBins1;
	int new_pos = ell+b2*m_nOrders+b1*m_nOrders*m_nBins;
	m_nnn[new_pos] += threept.m_nnn[old_pos];
	m_rrr[new_pos] += threept.m_rrr[old_pos];
      }

    }
  }

  // compute zeta

  m_zeta.erase(m_zeta.begin(), m_zeta.end());

  for (size_t i=0; i<m_nBins; i++) 
    for (size_t j=0; j<m_nBins; j++)  {

      vector<double> nnn(m_nOrders), rrr(m_nOrders);
      for (size_t ell=0; ell<m_nOrders; ell++) {
	int index = ell+j*m_nOrders+i*m_nOrders*m_nBins;
	nnn[ell] = m_nnn[index];
	rrr[ell] = m_rrr[index];
      }

      vector<double> zeta = m_SzapudiSzalay_multipoles (nnn, rrr);

      for (size_t ell=0; ell<m_nOrders; ell++) {
	m_zeta.push_back(zeta[ell]);
      }
    }
  
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_all::set_parameters (const double rMin, const double rMax, const double binSize, const size_t nOrders)
{
  m_rMin = rMin;
  m_rMax = rMax;
  m_binSize = binSize;
  m_nBins = static_cast<size_t>((m_rMax-m_rMin)/m_binSize);
  m_nOrders = nOrders;
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_all::measure (const ErrorType errorType, const std::string dir_output_triplets, const std::vector<std::string> dir_input_triplets, const int nResamplings, const bool count_triplets, const bool tcount, const int seed)
{
  coutCBL << "I'm computing the three-point correlation multipoles..." << endl;

  (void)errorType; 
  (void)nResamplings;
  (void)tcount;
  (void)seed;

  if (count_triplets) { // count triplets

    m_nnn.resize(m_nOrders*m_nBins*m_nBins, 0);
    m_rrr.resize(m_nOrders*m_nBins*m_nBins, 0);

    size_t start, end;
    for (size_t i=0; i<m_nSplit; i++) {

      if (m_nSplit>1)
	coutCBL << "Computing the three-point correlation multipoles, " << i+1 << "-th splitting" << endl; 
      
      start = (m_splitFactor<0) ? 0 :  m_data->nObjects()*(i*m_splitFactor);
      end = (m_splitFactor<0) ? m_random->nObjects() : m_data->nObjects()*((i+1)*m_splitFactor);

      auto joined = m_join_catalogues(*m_data, *m_random, start, end);
      
      vector<double> nnn, rrr;
      m_count_triplets(nnn, rrr, m_rMin, m_rMax, m_nBins, m_nOrders, *joined);

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
  else {
    m_read_triplets(m_nnn, dir_input_triplets, "nnn_legendre.dat");
    m_read_triplets(m_rrr, dir_input_triplets, "rrr_legendre.dat");
  }

  /// compute zeta
  
  m_zeta.erase(m_zeta.begin(), m_zeta.end());

  for (size_t i=0; i<m_nBins; i++) 
    for (size_t j=0; j<m_nBins; j++)  {

      vector<double> nnn(m_nOrders), rrr(m_nOrders);
      for (size_t ell=0; ell<m_nOrders; ell++) {
	int index = ell+j*m_nOrders+i*m_nOrders*m_nBins;
	nnn[ell] = m_nnn[index];
	rrr[ell] = m_rrr[index];
      }

      vector<double> zeta = m_SzapudiSzalay_multipoles (nnn, rrr);

      for (size_t ell=0; ell<m_nOrders; ell++) {
	m_zeta.push_back(zeta[ell]);
      }
    }

  coutCBL << "Done!" << endl;
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_all::write (const std::string dir, const std::string file) const 
{
  string mkdir = "mkdir -p "+dir;
  if(system(mkdir.c_str())) {}

  string outFile = dir+file;
  ofstream fout (outFile.c_str());

  fout << "#r12 r13 degree zeta_l index" << endl;

  for (size_t i=0; i<m_nBins; i++) 
    for (size_t j=0; j<m_nBins; j++) 
      for (size_t ell=0; ell<m_nOrders; ell++) {
	int index = ell+j*m_nOrders+i*m_nOrders*m_nBins;
	fout << setprecision(10) <<  m_rMin+(i+0.5)*m_binSize << " " << m_rMin+(j+0.5)*m_binSize << " " << ell << " " << m_zeta[index] << " " << index << endl;
      }

  fout.clear(); fout.close(); 

  coutCBL << "I wrote the file " << outFile << endl;
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_multipoles_all::resum (const std::string dir, const std::string file, const cbl::triplets::TripletType tripletType, const int nBins, const bool bin) const
{
  (void)dir;
  (void)file;
  (void)tripletType;
  (void)nBins;
  (void)bin;
  ErrorCBL("the function has not been implemented yet!", "resum", "ThreePointCorrelation_comoving_multipoles_all.cpp", glob::ExitCode::_workInProgress_);
}
