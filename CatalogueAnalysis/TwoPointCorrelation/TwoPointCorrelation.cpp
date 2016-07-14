/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli                          *
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
 *  @file TwoPointCorrelation/TwoPointCorrelation.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation 
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation, used to handle catalogues of astronomical sources
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "TwoPointCorrelation1D_monopole.h"
#include "TwoPointCorrelation1D_angular.h"
#include "TwoPointCorrelation_deprojected.h"
#include "TwoPointCorrelation_multipoles.h"
#include "TwoPointCorrelation_wedges.h"

using namespace cosmobl;
using namespace catalogue;
using namespace chainmesh;
using namespace data; 
using namespace pairs;
using namespace twopt;


// ============================================================================


shared_ptr<TwoPointCorrelation> cosmobl::twopt::TwoPointCorrelation::Create (const TwoPType type, const Catalogue data, const Catalogue random, const binType binType, const double Min, const double Max, const int nbins, const double shift, const CoordUnits angularUnits, function<double(double)> angularWeight)
{
  if (type==_1D_monopole_) return move(unique_ptr<TwoPointCorrelation1D_monopole>(new TwoPointCorrelation1D_monopole(data, random, binType, Min, Max, nbins, shift, angularUnits, angularWeight)));
  else if (type==_1D_angular_) return move(unique_ptr<TwoPointCorrelation1D_angular>(new TwoPointCorrelation1D_angular(data, random, binType, Min, Max, nbins, shift, angularUnits, angularWeight)));
  
  else ErrorMsg("Error in cosmobl::twopt::TwoPointCorrelation::Create of TwoPointCorrelation.cpp: no such type of object, or error in the input parameters!");
  
  return NULL;
}


// ============================================================================


shared_ptr<TwoPointCorrelation> cosmobl::twopt::TwoPointCorrelation::Create (const TwoPType type, const Catalogue data, const Catalogue random, const binType binType, const double Min, const double Max, const double binSize, const double shift, const CoordUnits angularUnits, function<double(double)> angularWeight)
{
  if (type==_1D_monopole_) return move(unique_ptr<TwoPointCorrelation1D_monopole>(new TwoPointCorrelation1D_monopole(data, random, binType, Min, Max, binSize, shift, angularUnits, angularWeight)));
  else if (type==_1D_angular_) return move(unique_ptr<TwoPointCorrelation1D_angular>(new TwoPointCorrelation1D_angular(data, random, binType, Min, Max, binSize, shift, angularUnits, angularWeight)));
 
  else ErrorMsg("Error in cosmobl::twopt::TwoPointCorrelation::Create of TwoPointCorrelation.cpp: no such type of object, or error in the input parameters!!");
  
  return NULL;
}


// ============================================================================


shared_ptr<TwoPointCorrelation> cosmobl::twopt::TwoPointCorrelation::Create (const TwoPType type, const Catalogue data, const Catalogue random, const binType binType_D1, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const binType binType_D2, const double Min_D2, const double Max_D2, const int nbins_D2, const double shift_D2, const CoordUnits angularUnits, function<double(double)> angularWeight)
{
  if (type==_2D_Cartesian_) return move(unique_ptr<TwoPointCorrelation2D_cartesian>(new TwoPointCorrelation2D_cartesian(data, random, binType_D1, Min_D1, Max_D1, nbins_D1, shift_D1, binType_D2, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)));

  else if (type==_2D_polar_) return move(unique_ptr<TwoPointCorrelation2D_polar>(new TwoPointCorrelation2D_polar(data, random, binType_D1, Min_D1, Max_D1, nbins_D1, shift_D1, binType_D2, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)));

  else ErrorMsg("Error in cosmobl::twopt::TwoPointCorrelation::Create of TwoPointCorrelation.cpp: no such type of object, or error in the input parameters!!");
  
  return NULL;
}


// ============================================================================


shared_ptr<TwoPointCorrelation> cosmobl::twopt::TwoPointCorrelation::Create (const TwoPType type, const Catalogue data, const Catalogue random, const binType binType_D1, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const binType binType_D2, const double Min_D2, const double Max_D2, const double binSize_D2, const double shift_D2, const CoordUnits angularUnits, function<double(double)> angularWeight)
{
  if (type==_2D_Cartesian_) return move(unique_ptr<TwoPointCorrelation2D_cartesian>(new TwoPointCorrelation2D_cartesian(data, random, binType_D1, Min_D1, Max_D1, binSize_D1, shift_D1, binType_D2, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)));

  else if (type==_2D_polar_) return move(unique_ptr<TwoPointCorrelation2D_polar>(new TwoPointCorrelation2D_polar(data, random, binType_D1, Min_D1, Max_D1, binSize_D1, shift_D1, binType_D2, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)));

  else ErrorMsg("Error in cosmobl::twopt::TwoPointCorrelation::Create of TwoPointCorrelation.cpp: no such type of object, or error in the input parameters!!");
  
  return NULL;
}


// ============================================================================


shared_ptr<TwoPointCorrelation> cosmobl::twopt::TwoPointCorrelation::Create (const TwoPType type, const Catalogue data, const Catalogue random, const binType binType_D1, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const double Min_D2, const double Max_D2, const int nbins_D2, const double shift_D2, const double piMax_integral, const CoordUnits angularUnits, function<double(double)> angularWeight)
{
  if (type==_1D_projected_) return move(unique_ptr<TwoPointCorrelation_projected>(new TwoPointCorrelation_projected(data, random, binType_D1, Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, piMax_integral, angularUnits, angularWeight)));

  else if (type==_1D_deprojected_) return move(unique_ptr<TwoPointCorrelation_deprojected>(new TwoPointCorrelation_deprojected(data, random, Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, piMax_integral, angularUnits, angularWeight)));

  else if (type==_1D_multipoles_) return move(unique_ptr<TwoPointCorrelation_multipoles>(new TwoPointCorrelation_multipoles(data, random, binType_D1, Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)));

  else if (type==_1D_wedges_) return move(unique_ptr<TwoPointCorrelation_wedges>(new TwoPointCorrelation_wedges(data, random, binType_D1, Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)));

  else ErrorMsg("Error in cosmobl::twopt::TwoPointCorrelation::Create of TwoPointCorrelation.cpp: no such type of object, or error in the input parameters!!");
  
  return NULL;
}


// ============================================================================


shared_ptr<TwoPointCorrelation> cosmobl::twopt::TwoPointCorrelation::Create (const TwoPType type, const Catalogue data, const Catalogue random, const binType binType_D1, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const double Min_D2, const double Max_D2, const double binSize_D2, const double shift_D2, const double piMax_integral, const CoordUnits angularUnits, function<double(double)> angularWeight)
{
  if (type==_1D_projected_) return move(unique_ptr<TwoPointCorrelation_projected>(new TwoPointCorrelation_projected(data, random, binType_D1, Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, piMax_integral, angularUnits, angularWeight)));

  else if (type==_1D_deprojected_) return move(unique_ptr<TwoPointCorrelation_deprojected>(new TwoPointCorrelation_deprojected(data, random, Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, piMax_integral, angularUnits, angularWeight)));

  else if (type==_1D_multipoles_) return move(unique_ptr<TwoPointCorrelation_multipoles>(new TwoPointCorrelation_multipoles(data, random, binType_D1, Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)));

  else if (type==_1D_wedges_) return move(unique_ptr<TwoPointCorrelation_wedges>(new TwoPointCorrelation_wedges(data, random, binType_D1, Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)));

  else ErrorMsg("Error in cosmobl::twopt::TwoPointCorrelation::Create of TwoPointCorrelation.cpp: no such type of object, or error in the input parameters!!");
  
  return NULL;
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation::count_pairs (const shared_ptr<Catalogue> cat1, const ChainMesh_Catalogue &ChM, shared_ptr<Pair> pp, const bool cross, const bool tcount)  
{
  // timer 
  time_t start; time (&start);
  int dp = cout.precision();
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);
  
  // number of objects in the first catalogue
  int nObj = cat1->nObjects();

  // factor used by the timer
  float fact_count = 100./nObj;
  
  // pointer to the second catalogue (get from the chain mesh)
  shared_ptr<Catalogue> cat2 = ChM.catalogue();

  // thread number
  int tid = 0;
  
  // start the multithreading parallelization
#pragma omp parallel num_threads(omp_get_max_threads()) private(tid)
  
  {
    tid = omp_get_thread_num();

    // if (tid == 0) cout << "Number of threads = " << omp_get_num_threads() << endl;
    
    // internal object used by each thread to handle pairs
    shared_ptr<Pair> pp_thread = (pp->pairDim()==_1D_) ? move(Pair::Create(pp->pairType(), pp->sMin(), pp->sMax(), pp->nbins(), pp->shift(), pp->angularUnits(), pp->angularWeight()))
      : move(Pair::Create(pp->pairType(), pp->sMin_D1(), pp->sMax_D1(), pp->nbins_D1(), pp->shift_D1(), pp->sMin_D2(), pp->sMax_D2(), pp->nbins_D2(), pp->shift_D2(), pp->angularUnits(), pp->angularWeight()));
    
    
    // parallelized loop
#pragma omp for schedule(static, 2)

    // loop on the objects of the first catalogue    
    for (int i=0; i<nObj; i++) { 
      
      // get the indexes of objects in the second catalogue that are close to the object i of the first catalogue
      // (i.e. objects inside the close cells of the chain-mesh)
      vector<long> close_objects = ChM.close_objects(cat1->coordinates(i), (cross) ? -1 : (long)i);
      
      for (auto &&j : close_objects) // loop on the nearby objects 
	pp_thread->put(cat1->catalogue_object(i), cat2->catalogue_object(j)); // estimate the distance between the two objects and update the pair count
      
      // estimate the computational time and update the time count
      time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
      if (tcount && tid==0) { cout << "\r..." << float(i)*fact_count << "% completed (" << diff_temp << " seconds)\r"; cout.flush(); }    
      if (i==int(nObj*0.25)) cout << "......25% completed" << endl;
      if (i==int(nObj*0.5)) cout << "......50% completed" << endl;
      if (i==int(nObj*0.75)) cout << "......75% completed"<< endl;   
    }

#pragma omp critical
    {
      // sum all the object pairs computed by each thread
      pp->Sum(pp_thread);
    }
    
  }
  

  // show the time spent by the method
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (tid==0) {
    if (diff<60) cout << "   time spent to count the pairs: " << diff << " seconds" << endl ;
    else if (diff<3600) cout << "   time spent to count the pairs: " << diff/60 << " minutes" << endl ;
    else cout << "   time spent to count the pairs: " << diff/3600 << " hours" << endl ;
  }
  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(dp);
  
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation::count_allPairs (const TwoPType type, const string dir_output_pairs, const vector<string> dir_input_pairs, const int count_dd, const int count_rr, const int count_dr, const bool tcount)  
{
  // ----------- compute polar coordinates, if necessary ----------- 

  if (!isSet(m_data->var(Var::_RA_)) || !isSet(m_data->var(Var::_Dec_)) || !isSet(m_data->var(Var::_Dc_))) 
    m_data->computePolarCoordinates();

  if (!isSet(m_random->var(Var::_RA_)) || !isSet(m_random->var(Var::_Dec_)) || !isSet(m_random->var(Var::_Dc_))) 
    m_random->computePolarCoordinates();

  if (type == _1D_angular_){
    m_data->normalizeComovingCoordinates();
    m_random->normalizeComovingCoordinates();
  }

  
  // ----------- create the chain-mesh ----------- 

  double rMAX;

  if (type==_1D_monopole_)
    rMAX = m_dd->sMax();

  else if (type==_1D_angular_) {
    double xx, yy, zz;
    cartesian_coord(radians(m_dd->sMax(), m_dd->angularUnits()), radians(m_dd->sMax(), m_dd->angularUnits()), 1., xx, yy, zz);
    rMAX = max(xx, zz);
  }

  else rMAX = m_dd->sMax_D1();

  
  double cell_size = rMAX*0.1; // to be optimized!!!
  
  ChainMesh_Catalogue ChM_data, ChM_random;

  if (count_dd==1)
    ChM_data.set_par(cell_size, m_data, rMAX);
  
  if (count_rr==1 || count_dr==1) 
    ChM_random.set_par(cell_size, m_random, rMAX);
  
  
  // ----------- count the number of pairs or read them from file -----------

  string file;
  
  if (count_dd>0) cout << endl << par::col_green << "data-data" << par::col_default << endl;

  file = "dd.dat";
 
  if (count_dd==1) {
    count_pairs(m_data, ChM_data, m_dd, 0, tcount);
    if (dir_output_pairs!=par::defaultString) write_pairs(m_dd, dir_output_pairs, file);
  }
  else if (count_dd==0) read_pairs(m_dd, dir_input_pairs, file);
  
  
  if (count_rr>0) cout << endl << par::col_green << "random-random" << par::col_default << endl;

  file = "rr.dat";
 
  if (count_rr==1) {
    count_pairs(m_random, ChM_random, m_rr, 0, tcount);
    if (dir_output_pairs!=par::defaultString) write_pairs(m_rr, dir_output_pairs, file);
  }
  else if (count_rr==0) read_pairs(m_rr, dir_input_pairs, file);

  
  if (count_dr>0) cout << endl << par::col_green << "data-random" << par::col_default << endl;

  file = "dr.dat";
 
  if (count_dr==1) {
    count_pairs(m_data, ChM_random, m_dr, 1, tcount);
    if (dir_output_pairs!=par::defaultString) write_pairs(m_dr, dir_output_pairs, file);
  }
  else if (count_dr==0) read_pairs(m_dr, dir_input_pairs, file);

  if (count_dd==1)
    m_data->Order();
  
  if (count_rr==1 || count_dr==1)
    m_random->Order();

  if (type==_1D_angular_) {
    m_data->restoreComovingCoordinates();
    m_random->restoreComovingCoordinates();
  }
  
}


// ============================================================================


double TwoPointCorrelation::PoissonError (const double dd, const double rr, const double dr, const int nData, const int nRandom) const
{
  double norm1 = double(nRandom)*double(nRandom-1)/(double(nData)*double(nData-1));
  double norm2 = double(nRandom-1)/double(nData);
  
  return sqrt(pow(norm1*sqrt(dd)/rr,2)+pow(norm2*sqrt(dr)/rr,2)+pow((norm1*dd-norm2*dr)*pow(rr,-1.5),2));
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation::count_pairs_region (const shared_ptr<Catalogue> cat1, const ChainMesh_Catalogue &ChM, shared_ptr<Pair> pp, vector<shared_ptr<Pair>> pp_regions, const bool cross, const bool tcount)  
{
  // timer 
  time_t start; time (&start);
  int dp = cout.precision();
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);

  // number of objects in the first catalogue
  int nObj = cat1->nObjects();

  // factor used by the timer
  float fact_count = 100./nObj;

  // pointer to the second catalogue (get from the chain mesh)
  shared_ptr<Catalogue> cat2 = ChM.catalogue();

  // thread number
  int tid = 0;
   
#pragma omp parallel num_threads(omp_get_max_threads()) private(tid)

  {
    tid = omp_get_thread_num();

    vector<shared_ptr<Pair> > pp_thread(pp_regions.size());

    for (size_t i=0; i<pp_regions.size(); i++) 
      pp_thread[i] = (pp->pairDim()==_1D_) ? move(Pair::Create(pp->pairType(), pp->sMin(), pp->sMax(), pp->nbins(), pp->shift(), pp->angularUnits(), pp->angularWeight()))
	: move(Pair::Create(pp->pairType(), pp->sMin_D1(), pp->sMax_D1(), pp->nbins_D1(), pp->shift_D1(), pp->sMin_D2(), pp->sMax_D2(), pp->nbins_D2(), pp->shift_D2(), pp->angularUnits(), pp->angularWeight()));
    
    int nRegions = cat1->Nregion();
    
    // parallelized loop
#pragma omp for schedule(static, 2)
    for (int i=0; i<nObj; i++) {

      vector<long int> close_objects = ChM.close_objects(cat1->coordinates(i), (cross) ? -1 : i);

      for (auto &&j : close_objects) {      
	int reg1 = min(cat1->region(i), cat2->region(j));
	int reg2 = max(cat1->region(i), cat2->region(j));
	int index = reg1*nRegions+reg2-(reg1-1)*reg1/2-reg1; 
	pp_thread[index]->put(cat1->catalogue_object(i), cat2->catalogue_object(j));
      }
      
      // estimate the computational time and update the time count
      time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
      if (tcount && tid==0) { cout << "\r..." << float(i)*fact_count << "% completed (" << diff_temp << " seconds)\r"; cout.flush(); }    
      if (i==int(nObj*0.25)) cout << "......25% completed" << endl;
      if (i==int(nObj*0.5)) cout << "......50% completed" << endl;
      if (i==int(nObj*0.75)) cout << "......75% completed"<< endl;   
    }
    
#pragma omp critical
    {
      // sum all the object pairs computed by each thread
      for(size_t i =0;i<pp_regions.size();i++)
	pp_regions[i]->Sum(pp_thread[i]);
    }

  }
  

  // show the time spent by the method

  time_t end; time (&end);
  double diff = difftime(end,start);
  if (tid==0) {
    if (diff<60) cout << "   time spent to count the pairs: " << diff << " seconds" << endl ;
    else if (diff<3600) cout << "   time spent to count the pairs: " << diff/60 << " minutes" << endl ;
    else cout << "   time spent to count the pairs: " << diff/3600 << " hours" << endl ;
  }
  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(dp);

  
  // sum pairs for the entire sample

  switch (pp->pairDim()) {

    case _1D_:
      for (int i=0; i<pp->nbins(); i++)
	for (size_t r=0; r<pp_regions.size(); r++)
	  pp->add_PP1D(i,pp_regions[r]->PP1D(i));
      break;

    case _2D_:
      for (int i=0; i<pp->nbins_D1(); i++)
	for (int j=0; j<pp->nbins_D2(); j++)
	  for (size_t r=0; r<pp_regions.size(); r++)
	    pp->add_PP2D(i,j,pp_regions[r]->PP2D(i,j));
      break;
      
    default:
      ErrorMsg("Error in count_pairs_region of TwoPointCorrelation.cpp, no such type of Dimension for pair");
      break;
  }
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation::count_allPairs_region (vector<shared_ptr<Pair> > &dd_regions, vector<shared_ptr<Pair> > &rr_regions, vector<shared_ptr<Pair> > &dr_regions, const TwoPType type, const string dir_output_pairs, const vector<string> dir_input_pairs, const int count_dd, const int count_rr, const int count_dr, const bool tcount)  
{
  // ----------- compute polar coordinates, if necessary ----------- 

  if (!isSet(m_data->var(Var::_RA_)) || !isSet(m_data->var(Var::_Dec_)) || !isSet(m_data->var(Var::_Dc_))) 
    m_data->computePolarCoordinates();

  if (!isSet(m_random->var(Var::_RA_)) || !isSet(m_random->var(Var::_Dec_)) || !isSet(m_random->var(Var::_Dc_))) 
    m_random->computePolarCoordinates();

  if (type == _1D_angular_){
    m_data->normalizeComovingCoordinates();
    m_random->normalizeComovingCoordinates();
  }

  // ----------- create the chain-mesh ----------- 

  double rMAX = (type==_1D_monopole_ || type==_1D_angular_) ? m_dd->sMax() : m_dd->sMax_D1(); // check!!!

  double cell_size = rMAX*0.1; // to be optimized!!!

  ChainMesh_Catalogue ChM_data, ChM_random;

  if (count_dd==1)
    ChM_data.set_par(cell_size, m_data, rMAX);

  if (count_rr==1 || count_dr==1) 
    ChM_random.set_par(cell_size, m_random, rMAX);

  // ----------- Initialize pair vector for resampling ----------- 

  vector<long> region_list = m_data->get_region_list();
  int nRegions = region_list.size();
  int nP = nRegions*(nRegions+1)/2;

  dd_regions.erase(dd_regions.begin(), dd_regions.end());
  rr_regions.erase(rr_regions.begin(), rr_regions.end());
  dr_regions.erase(dr_regions.begin(), dr_regions.end());

  for (int i=0; i<nP; i++) {

    dd_regions.push_back((m_dd->pairDim()==_1D_) ? move(Pair::Create(m_dd->pairType(), m_dd->sMin(), m_dd->sMax(), m_dd->nbins(), m_dd->shift(), m_dd->angularUnits(), m_dd->angularWeight())) : move(Pair::Create(m_dd->pairType(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2(), m_dd->angularUnits(), m_dd->angularWeight())));
    
    rr_regions.push_back((m_rr->pairDim()==_1D_) ? move(Pair::Create(m_rr->pairType(), m_rr->sMin(), m_rr->sMax(), m_rr->nbins(), m_rr->shift(), m_rr->angularUnits(), m_rr->angularWeight())) : move(Pair::Create(m_rr->pairType(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2(), m_rr->angularUnits(), m_rr->angularWeight())));

    dr_regions.push_back((m_dr->pairDim()==_1D_) ? move(Pair::Create(m_dr->pairType(), m_dr->sMin(), m_dr->sMax(), m_dr->nbins(), m_dr->shift(), m_dr->angularUnits(), m_dr->angularWeight())) : move(Pair::Create(m_dr->pairType(), m_dr->sMin_D1(), m_dr->sMax_D1(), m_dr->nbins_D1(), m_dr->shift_D1(), m_dr->sMin_D2(), m_dr->sMax_D2(), m_dr->nbins_D2(), m_dr->shift_D2(), m_dr->angularUnits(), m_dr->angularWeight())));

  }


  // ----------- count the number of pairs or read them from file -----------

  string file, file_regions;

  cout << endl << par::col_green << "data-data" << par::col_default << endl;

  file = "dd.dat";
  file_regions = "dd_regions.dat";

  if (count_dd==1) {
    count_pairs_region(m_data, ChM_data, m_dd, dd_regions, 0, tcount);
    if (dir_output_pairs!=par::defaultString) {
      write_pairs(m_dd, dir_output_pairs, file);
      write_pairs(dd_regions,dir_output_pairs,file_regions);
    }
  }
  else if (count_dd==0) {
    read_pairs(m_dd, dir_input_pairs, file);
    read_pairs(dd_regions, dir_input_pairs, file_regions);
  }

  cout << endl << par::col_green << "random-random" << par::col_default << endl;

  file = "rr.dat";
  file_regions = "rr_regions.dat";

  if (count_rr==1) {
    count_pairs_region(m_random, ChM_random, m_rr, rr_regions, 0, tcount);
    if (dir_output_pairs!=par::defaultString) {
      write_pairs(m_rr, dir_output_pairs, file);
      write_pairs(rr_regions,dir_output_pairs,file_regions);
    }
  }
  else if (count_rr==0) {
    read_pairs(m_rr, dir_input_pairs, file);
    read_pairs(rr_regions, dir_input_pairs, file_regions);
  }

  cout << endl << par::col_green << "data-random" << par::col_default << endl;

  file = "dr.dat";
  file_regions = "dr_regions.dat";

  if (count_dr==1) {
    count_pairs_region(m_data, ChM_random, m_dr, dr_regions, 1, tcount);
    if (dir_output_pairs!=par::defaultString) {
      write_pairs(m_dr, dir_output_pairs, file);
      write_pairs(dr_regions, dir_output_pairs, file_regions);
    }
  }
  else if (count_dr==0) {
    read_pairs(m_dr, dir_input_pairs, file);
    read_pairs(dr_regions, dir_input_pairs, file_regions);
  }
  
  if (count_dd==1)
    m_data->Order();
  
  if (count_rr==1 || count_dr==1) 
    m_random->Order();
 
  if (type == _1D_angular_) {
    m_data->restoreComovingCoordinates();
    m_random->restoreComovingCoordinates();
  }

}




