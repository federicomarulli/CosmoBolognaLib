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
 *  TwoPointCorrelation, used to handle catalogues of astronomical
 *  sources
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unibo.it
 */

#include "TwoPointCorrelation.h"
#include "TwoPointCorrelation1D_angular.h"
#include "TwoPointCorrelation1D_monopole.h"
#include "TwoPointCorrelation1D_filtered.h"
#include "TwoPointCorrelation_deprojected.h"
#include "TwoPointCorrelation_multipoles_direct.h"
#include "TwoPointCorrelation_multipoles_integrated.h"
#include "TwoPointCorrelation_wedges.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace chainmesh;
using namespace data; 
using namespace pairs;
using namespace measure::twopt;


// ============================================================================


shared_ptr<TwoPointCorrelation> cbl::measure::twopt::TwoPointCorrelation::Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType, const double Min, const double Max, const int nbins, const double shift, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight, const bool compute_extra_info, const double random_dilution_fraction)
{
  if (type==TwoPType::_angular_) return move(unique_ptr<TwoPointCorrelation1D_angular>(new TwoPointCorrelation1D_angular(data, random, binType, Min, Max, nbins, shift, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));
  
  else if (type==TwoPType::_monopole_) return move(unique_ptr<TwoPointCorrelation1D_monopole>(new TwoPointCorrelation1D_monopole(data, random, binType, Min, Max, nbins, shift, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));

  else if (type==TwoPType::_multipoles_direct_) return move(unique_ptr<TwoPointCorrelation_multipoles_direct>(new TwoPointCorrelation_multipoles_direct(data, random, binType, Min, Max, nbins, shift, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));
  
  else ErrorCBL("no such type of object, or error in the input parameters!", "Create", "TwoPointCorrelation.cpp");
  
  return NULL;
}


// ============================================================================


shared_ptr<TwoPointCorrelation> cbl::measure::twopt::TwoPointCorrelation::Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType, const double Min, const double Max, const double binSize, const double shift, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight, const bool compute_extra_info, const double random_dilution_fraction)
{
  if (type==TwoPType::_angular_) return move(unique_ptr<TwoPointCorrelation1D_angular>(new TwoPointCorrelation1D_angular(data, random, binType, Min, Max, binSize, shift, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));
  
  else if (type==TwoPType::_monopole_) return move(unique_ptr<TwoPointCorrelation1D_monopole>(new TwoPointCorrelation1D_monopole(data, random, binType, Min, Max, binSize, shift, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));

  else if (type==TwoPType::_multipoles_direct_) return move(unique_ptr<TwoPointCorrelation_multipoles_direct>(new TwoPointCorrelation_multipoles_direct(data, random, binType, Min, Max, binSize, shift, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));
  
  else ErrorCBL("no such type of object, or error in the input parameters!", "Create", "TwoPointCorrelation.cpp");
  
  return NULL;
}


// ============================================================================


shared_ptr<TwoPointCorrelation> cbl::measure::twopt::TwoPointCorrelation::Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const double Min_D2, const double Max_D2, const int nbins_D2, const double shift_D2, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight, const bool compute_extra_info, const double random_dilution_fraction)
{
  if (type==TwoPType::_multipoles_integrated_) return move(unique_ptr<TwoPointCorrelation_multipoles_integrated>(new TwoPointCorrelation_multipoles_integrated(data, random, binType_D1, Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));
  
  else if (type==TwoPType::_filtered_) return move(unique_ptr<TwoPointCorrelation1D_filtered>(new TwoPointCorrelation1D_filtered(data, random, binType_D1, Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));
  
  else ErrorCBL("no such type of object, or error in the input parameters!", "Create", "TwoPointCorrelation.cpp");
  
  return NULL;
}


// ============================================================================


shared_ptr<TwoPointCorrelation> cbl::measure::twopt::TwoPointCorrelation::Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const double Min_D2, const double Max_D2, const double binSize_D2, const double shift_D2, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight, const bool compute_extra_info, const double random_dilution_fraction)
{
  if (type==TwoPType::_multipoles_integrated_) return move(unique_ptr<TwoPointCorrelation_multipoles_integrated>(new TwoPointCorrelation_multipoles_integrated(data, random, binType_D1, Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));

  else if (type==TwoPType::_filtered_) return move(unique_ptr<TwoPointCorrelation1D_filtered>(new TwoPointCorrelation1D_filtered(data, random, binType_D1, Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));
  
  else ErrorCBL("no such type of object, or error in the input parameters!", "Create", "TwoPointCorrelation.cpp");
  
  return NULL;
}


// ============================================================================


shared_ptr<TwoPointCorrelation> cbl::measure::twopt::TwoPointCorrelation::Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const double Min_D2, const double Max_D2, const int nbins_D2, const double shift_D2, const double piMax_integral, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight, const bool compute_extra_info, const double random_dilution_fraction)
{
  if (type==TwoPType::_projected_) return move(unique_ptr<TwoPointCorrelation_projected>(new TwoPointCorrelation_projected(data, random, binType_D1, Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, piMax_integral, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));

  else if (type==TwoPType::_deprojected_) return move(unique_ptr<TwoPointCorrelation_deprojected>(new TwoPointCorrelation_deprojected(data, random, Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, piMax_integral, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));
  
  else ErrorCBL("no such type of object, or error in the input parameters!", "Create", "TwoPointCorrelation.cpp");
  
  return NULL;
}


// ============================================================================


shared_ptr<TwoPointCorrelation> cbl::measure::twopt::TwoPointCorrelation::Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const double Min_D2, const double Max_D2, const double binSize_D2, const double shift_D2, const double piMax_integral, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight, const bool compute_extra_info, const double random_dilution_fraction)
{
  if (type==TwoPType::_projected_) return move(unique_ptr<TwoPointCorrelation_projected>(new TwoPointCorrelation_projected(data, random, binType_D1, Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, piMax_integral, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));

  else if (type==TwoPType::_deprojected_) return move(unique_ptr<TwoPointCorrelation_deprojected>(new TwoPointCorrelation_deprojected(data, random, Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, piMax_integral, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));
  
  else ErrorCBL("no such type of object, or error in the input parameters!", "Create", "TwoPointCorrelation.cpp");
  
  return NULL;
}


// ============================================================================


shared_ptr<TwoPointCorrelation> cbl::measure::twopt::TwoPointCorrelation::Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const int nWedges, const int nbins_D2, const double shift_D2, const std::vector<std::vector<double>> mu_integral_limits, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight, const bool compute_extra_info, const double random_dilution_fraction)
{
  if (type==TwoPType::_wedges_) return move(unique_ptr<TwoPointCorrelation_wedges>(new TwoPointCorrelation_wedges(data, random, binType_D1, Min_D1, Max_D1, nbins_D1, shift_D1, nWedges, nbins_D2, shift_D2, mu_integral_limits, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));
  
  else ErrorCBL("no such type of object, or error in the input parameters!", "Create", "TwoPointCorrelation.cpp");
  
  return NULL;
}


// ============================================================================


shared_ptr<TwoPointCorrelation> cbl::measure::twopt::TwoPointCorrelation::Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const int nWedges, const double binSize_D2, const double shift_D2, const std::vector<std::vector<double>> mu_integral_limits, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight, const bool compute_extra_info, const double random_dilution_fraction)
{
  if (type==TwoPType::_wedges_) return move(unique_ptr<TwoPointCorrelation_wedges>(new TwoPointCorrelation_wedges(data, random, binType_D1, Min_D1, Max_D1, binSize_D1, shift_D1, nWedges, binSize_D2, shift_D2, mu_integral_limits, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));
  
  else ErrorCBL("no such type of object, or error in the input parameters!", "Create", "TwoPointCorrelation.cpp");
  
  return NULL;
}


// ============================================================================


shared_ptr<TwoPointCorrelation> cbl::measure::twopt::TwoPointCorrelation::Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const BinType binType_D2, const double Min_D2, const double Max_D2, const int nbins_D2, const double shift_D2, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight, const bool compute_extra_info, const double random_dilution_fraction)
{
  if (type==TwoPType::_2D_Cartesian_) return move(unique_ptr<TwoPointCorrelation2D_cartesian>(new TwoPointCorrelation2D_cartesian(data, random, binType_D1, Min_D1, Max_D1, nbins_D1, shift_D1, binType_D2, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));

  else if (type==TwoPType::_2D_polar_) return move(unique_ptr<TwoPointCorrelation2D_polar>(new TwoPointCorrelation2D_polar(data, random, binType_D1, Min_D1, Max_D1, nbins_D1, shift_D1, binType_D2, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));

  else ErrorCBL("no such type of object, or error in the input parameters!", "Create", "TwoPointCorrelation.cpp");
  
  return NULL;
}


// ============================================================================

void cbl::measure::twopt::TwoPointCorrelation::resets ()
{
  // reset the data-data pairs
  m_dd->reset();
  // reset the random-random pairs
  m_rr->reset();
  // reset the data-random pairs
  m_dr->reset();
}

// ============================================================================


std::shared_ptr<TwoPointCorrelation> cbl::measure::twopt::TwoPointCorrelation::Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const BinType binType_D2, const double Min_D2, const double Max_D2, const double binSize_D2, const double shift_D2, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight, const bool compute_extra_info, const double random_dilution_fraction)
{
  if (type==TwoPType::_2D_Cartesian_) return move(unique_ptr<TwoPointCorrelation2D_cartesian>(new TwoPointCorrelation2D_cartesian(data, random, binType_D1, Min_D1, Max_D1, binSize_D1, shift_D1, binType_D2, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));

  else if (type==TwoPType::_2D_polar_) return move(unique_ptr<TwoPointCorrelation2D_polar>(new TwoPointCorrelation2D_polar(data, random, binType_D1, Min_D1, Max_D1, binSize_D1, shift_D1, binType_D2, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)));

  else ErrorCBL("no such type of object, or error in the input parameters!", "Create", "TwoPointCorrelation.cpp");
  
  return NULL;
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation::count_pairs (const std::shared_ptr<Catalogue> cat1, const ChainMesh_Catalogue &ChM, std::shared_ptr<Pair> pp, const bool cross, const bool tcount)
{ 
  // timer 
  time_t start; time(&start);
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

    // if (tid == 0) coutCBL << "Number of threads = " << omp_get_num_threads() << endl;

    
    // internal object used by each thread to handle pairs
    
    shared_ptr<Pair> pp_thread = (pp->pairDim()==Dim::_1D_) ? move(Pair::Create(pp->pairType(), pp->pairInfo(), pp->sMin(), pp->sMax(), pp->nbins(), pp->shift(), pp->angularUnits(), pp->angularWeight()))
      : move(Pair::Create(pp->pairType(), pp->pairInfo(), pp->sMin_D1(), pp->sMax_D1(), pp->nbins_D1(), pp->shift_D1(), pp->sMin_D2(), pp->sMax_D2(), pp->nbins_D2(), pp->shift_D2(), pp->angularUnits(), pp->angularWeight()));
    
    
    // parallelized loop
#pragma omp for schedule(static, 2)

    // loop on the objects of the first catalogue    
    for (int i=0; i<nObj; ++i) {

      // get the indexes of objects in the second catalogue that are close to the object i of the first catalogue
      // (i.e. objects inside the close cells of the chain-mesh)
      vector<long> close_objects = ChM.close_objects(cat1->coordinate(i), (cross) ? -1 : (long)i);
      
      // loop on the nearby objects
      for (auto &&j : close_objects)
	// estimate the distance between the two objects and update the pair count
	pp_thread->put(cat1->catalogue_object(i), cat2->catalogue_object(j));
    
      // estimate the computational time and update the time count
      time_t end_temp; time(&end_temp); double diff_temp = difftime(end_temp, start);
      if (tcount && tid==0) { coutCBL << "\r" << float(i)*fact_count << "% completed (" << diff_temp << " seconds)\r"; cout.flush(); }    
      if (i==int(nObj*0.25)) coutCBL << ".............25% completed   " << endl;
      if (i==int(nObj*0.5)) coutCBL << ".............50% completed   " << endl;
      if (i==int(nObj*0.75)) coutCBL << ".............75% completed   "<< endl;   
   
    }

#pragma omp critical
    {
      // sum all the object pairs computed by each thread
      pp->Sum(pp_thread);
    }
    
  }

  
  // show the time spent by the method
  time_t end; time (&end);
  double diff = difftime(end, start);
  if (tid==0) {
    if (diff<60) coutCBL << "   time spent to count the pairs: " << diff << " seconds" << endl ;
    else if (diff<3600) coutCBL << "   time spent to count the pairs: " << diff/60 << " minutes" << endl ;
    else coutCBL << "   time spent to count the pairs: " << diff/3600 << " hours" << endl ;
  }
  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(dp);

}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation::count_allPairs (const TwoPType type, const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const double fact)  
{
  // ----------- compute polar coordinates, if necessary ----------- 

  if (!m_data->isSetVar(Var::_RA_) || !m_data->isSetVar(Var::_Dec_) || !m_data->isSetVar(Var::_Dc_))
    m_data->computePolarCoordinates();
  
  if (!m_random->isSetVar(Var::_RA_) || !m_random->isSetVar(Var::_Dec_) || !m_random->isSetVar(Var::_Dc_))
    m_random->computePolarCoordinates();

  if (type==TwoPType::_angular_) {
    m_data->normalizeComovingCoordinates();
    m_random->normalizeComovingCoordinates();
  }

  
  // ----------- dilute the random catalogue used to compute the RR pairs (to improve the performance) ----------- 

  if (estimator==Estimator::_natural_ && m_random_dilution_fraction!=1.) {
    m_random_dilution_fraction = 1.;
    WarningMsgCBL("m_random_dilution_fraction = 1, since the random catalogue is not diluted when using the natural estimator!", "count_allPairs", "TwoPointCorrelation.cpp");
  }

  auto random_dil = make_shared<catalogue::Catalogue>(catalogue::Catalogue(move(m_random->diluted_catalogue(m_random_dilution_fraction))));
  
  
  // ----------- create the chain-mesh ----------- 

  double rMAX;
  double rMIN;

  if (type==TwoPType::_monopole_ || type==TwoPType::_multipoles_direct_ || type==TwoPType::_filtered_) {
    rMAX = m_dd->sMax();
    rMIN = m_dd->sMin();
  }

  else if (type==TwoPType::_angular_) {
    double xx, yy, zz;
    cartesian_coord(radians(m_dd->sMax(), m_dd->angularUnits()), radians(m_dd->sMax(), m_dd->angularUnits()), 1., xx, yy, zz);
    rMAX = sqrt(yy*yy+zz*zz);
    cartesian_coord(radians(m_dd->sMin(), m_dd->angularUnits()), radians(m_dd->sMin(), m_dd->angularUnits()), 1., xx, yy, zz);
    rMIN = sqrt(yy*yy+zz*zz);
  }

  else if (type==TwoPType::_2D_polar_ || type==TwoPType::_multipoles_integrated_ || type ==TwoPType::_wedges_) {
    rMAX = m_dd->sMax_D1();
    rMIN = m_dd->sMin_D1();
  }
  
  else if (type==TwoPType::_2D_Cartesian_ || type==TwoPType::_projected_ || type==TwoPType::_deprojected_) {
    rMAX = max(m_dd->sMax_D1(), m_dd->sMax_D2())*sqrt(2.);
    rMIN = max(m_dd->sMin_D1(), m_dd->sMin_D2())*sqrt(2.);
  }

  else
    ErrorCBL("the chosen two-point correlation function type is uknown!", "count_allPairs", "TwoPointCorrelation.cpp");


  double cell_size = fact*rMAX; // to be optimized!!!
  
  ChainMesh_Catalogue ChM_data, ChM_random, ChM_random_dil;
  
  if (count_dd)
    ChM_data.set_par(cell_size, m_data, rMAX, rMIN);
  
  if (count_rr)
    ChM_random_dil.set_par(cell_size, random_dil, rMAX, rMIN);

  if (count_dr)
    ChM_random.set_par(cell_size, m_random, rMAX, rMIN);

  // ----------- reset the pair counts --------------------------------------
  resets();
  
  // ----------- count the number of pairs or read them from file -----------

  string file;
  
  cout << endl; coutCBL << par::col_green << "data-data" << par::col_default << endl;
  file = "dd.dat";
 
  if (count_dd) {
    count_pairs(m_data, ChM_data, m_dd, false, tcount);
    if (dir_output_pairs!=par::defaultString) write_pairs(m_dd, dir_output_pairs, file);
  }
  else read_pairs(m_dd, dir_input_pairs, file);
  

  cout << endl; coutCBL << par::col_green << "random-random" << par::col_default << endl;
  file = "rr.dat";
 
  if (count_rr) {
    count_pairs(random_dil, ChM_random_dil, m_rr, false, tcount);
    if (dir_output_pairs!=par::defaultString) write_pairs(m_rr, dir_output_pairs, file);
  }
  else read_pairs(m_rr, dir_input_pairs, file);

  
  if (estimator==Estimator::_LandySzalay_) {
    
    cout << endl; coutCBL << par::col_green << "data-random" << par::col_default << endl; 
    file = "dr.dat";
    
    if (count_dr) {
      count_pairs(m_data, ChM_random, m_dr, true, tcount);
      if (dir_output_pairs!=par::defaultString) write_pairs(m_dr, dir_output_pairs, file);
    }
    else read_pairs(m_dr, dir_input_pairs, file);

  }
  
  if (count_dd)
    m_data->Order();
  
  if (count_rr || count_dr)
    m_random->Order();

  if (type==TwoPType::_angular_) {
    m_data->restoreComovingCoordinates();
    m_random->restoreComovingCoordinates();
  }

}


// ============================================================================


double cbl::measure::twopt::TwoPointCorrelation::PoissonError (const Estimator estimator, const double dd, const double rr, const double dr, const int nData, const int nRandom) const
{
  if (estimator!=Estimator::_natural_ && estimator!=Estimator::_LandySzalay_)
    ErrorCBL("The implementation of Poisson errors for the chosen estimator is not available yet!", "PoissonError", "TwoPointCorrelation.cpp", glob::ExitCode::_workInProgress_);

  double fR = m_random_dilution_fraction;
  if (estimator==Estimator::_natural_ && fR!=1) {
    fR = 1.;
    WarningMsgCBL("fR = 1, since the random catalogue is not diluted when using the natural estimator!", "PoissonError", "TwoPointCorrelation.cpp");
  }
  
  const double norm1 = double(nRandom)*double(nRandom-1)/(double(nData)*double(nData-1));
  const double norm2 = double(nRandom-1)/double(nData);

  const double T1 = pow(norm1*sqrt(dd)/rr, 2);
  const double T2 = pow(norm2*sqrt(dr)/rr, 2);
  const double T3 = pow((norm1*dd-norm2*dr)*pow(rr, -1.5), 2);

  if (min(T2, T3)>T1)
    WarningMsgCBL("enlarge the random sample, that dominates the Poisson errors!", "PoissonError", "TwoPointCorrelation.cpp");
  
  return pow(fR, 2)*sqrt(T1+T2+T3);
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation::count_pairs_region (const std::shared_ptr<Catalogue> cat1, const ChainMesh_Catalogue &ChM, std::shared_ptr<Pair> pp, std::vector<std::shared_ptr<Pair>> pp_regions, const bool cross, const bool tcount)  
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

  // set the number of regions (to be checked!)
  const int nRegions = cat2->nRegions();
    
#pragma omp parallel num_threads(omp_get_max_threads()) private(tid)
  {
    tid = omp_get_thread_num();

    vector<shared_ptr<Pair> > pp_thread(pp_regions.size());
    
    for (size_t i=0; i<pp_regions.size(); ++i) 
      pp_thread[i] = (pp->pairDim()==Dim::_1D_) ? move(Pair::Create(pp->pairType(), pp->pairInfo(), pp->sMin(), pp->sMax(), pp->nbins(), pp->shift(), pp->angularUnits(), pp->angularWeight()))
	: move(Pair::Create(pp->pairType(), pp->pairInfo(), pp->sMin_D1(), pp->sMax_D1(), pp->nbins_D1(), pp->shift_D1(), pp->sMin_D2(), pp->sMax_D2(), pp->nbins_D2(), pp->shift_D2(), pp->angularUnits(), pp->angularWeight()));
    
    
    // parallelized loop
#pragma omp for schedule(static, 2)
    for (int i=0; i<nObj; ++i) {  
      
      vector<long int> close_objects = ChM.close_objects(cat1->coordinate(i), (cross) ? -1 : i);

      for (auto &&j : close_objects) {      
	const int reg1 = (cross) ? cat1->region(i) : min(cat1->region(i), cat2->region(j));
	const int reg2 = (cross) ? cat2->region(j) : max(cat1->region(i), cat2->region(j));

	// giving an element from lower/upper triangular matrix
	const size_t index = (cross) ? reg1*nRegions+reg2 : reg1*nRegions+reg2-(reg1-1)*reg1/2-reg1;

	pp_thread[index]->put(cat1->catalogue_object(i), cat2->catalogue_object(j));
      }
      
      // estimate the computational time and update the time count
      time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
      if (tcount && tid==0) { coutCBL << "\r" << float(i)*fact_count << "% completed (" << diff_temp << " seconds)\r"; cout.flush(); }    
      if (i==int(nObj*0.25)) coutCBL << ".............25% completed   " << endl;
      if (i==int(nObj*0.5)) coutCBL << ".............50% completed   " << endl;
      if (i==int(nObj*0.75)) coutCBL << ".............75% completed   "<< endl;   
    }
    
#pragma omp critical
    {
      // sum all the object pairs computed by each thread
      for (size_t i=0; i<pp_regions.size(); ++i)
	pp_regions[i]->Sum(pp_thread[i]);
    }

  }

  
  // show the time spent by the method
  time_t end; time (&end);
  double diff = difftime(end, start);
  if (tid==0) {
    if (diff<60) coutCBL << "   time spent to count the pairs: " << diff << " seconds" << endl ;
    else if (diff<3600) coutCBL << "   time spent to count the pairs: " << diff/60 << " minutes" << endl ;
    else coutCBL << "   time spent to count the pairs: " << diff/3600 << " hours" << endl ;
  }
  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(dp);
  

  // sum the pairs of the entire sample

  switch (pp->pairDim()) {

    case Dim::_1D_:
      for (size_t i=0; i<pp->PP1D().size(); ++i)
        for (size_t r=0; r<pp_regions.size(); ++r)
          pp->add_data1D(i, pp_regions[r]);
      break;

    case Dim::_2D_:
      for (int i=0; i<pp->nbins_D1(); ++i)
        for (int j=0; j<pp->nbins_D2(); ++j)
          for (size_t r=0; r<pp_regions.size(); ++r)
            pp->add_data2D(i, j, pp_regions[r]);
      break;

    default:
      ErrorCBL("no such type of pair dimension!", "count_pairs_region", "TwoPointCorrelation.cpp");
      break;
  }
  
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation::count_pairs_region_test (const std::shared_ptr<Catalogue> cat1, const ChainMesh_Catalogue &ChM, std::shared_ptr<Pair> pp, std::vector<std::shared_ptr<Pair>> pp_res, const std::vector<double> weight, const bool cross, const bool tcount)  
{
  if(pp->pairDim()==Dim::_1D_)
    count_pairs_region_test_1D(cat1, ChM, pp, pp_res, weight, cross, tcount);
  else if(pp->pairDim()==Dim::_2D_)
    count_pairs_region_test_2D(cat1, ChM, pp, pp_res, weight, cross, tcount);
  else
    ErrorCBL("wrong pair type!", "count_pairs_region_test", "TwoPointCorrelation.cpp");
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation::count_pairs_region_test_1D (const std::shared_ptr<Catalogue> cat1, const ChainMesh_Catalogue &ChM, std::shared_ptr<Pair> pp, std::vector<std::shared_ptr<Pair>> pp_res, const std::vector<double> weight, const bool cross, const bool tcount)  
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

    shared_ptr<Pair> pp_thread;
    vector<shared_ptr<Pair> > pp_res_thread(pp_res.size());
    
    pp_thread = move(Pair::Create(pp->pairType(), pp->pairInfo(), pp->sMin(), pp->sMax(), pp->nbins(), pp->shift(), pp->angularUnits(), pp->angularWeight()));

    for (size_t i=0; i<pp_res.size(); ++i) 
      pp_res_thread[i] = move(Pair::Create(pp->pairType(), PairInfo::_standard_, pp->sMin(), pp->sMax(), pp->nbins(), pp->shift(), pp->angularUnits(), pp->angularWeight()));
    
    // parallelized loop
#pragma omp for schedule(static, 2)
    for (int i=0; i<nObj; ++i) {

      vector<long int> close_objects = ChM.close_objects(cat1->coordinate(i), (cross) ? -1 : i);

      for (auto &&j : close_objects) {      
	int kk;
	double wkk;

	pp_thread->get(cat1->catalogue_object(i), cat2->catalogue_object(j), kk, wkk);
	pp_thread->set(kk, wkk);

	vector<double> ww = weight;

	ww[cat1->region(i)] = 0;
	ww[cat2->region(j)] = 0;

	for(size_t k=0; k< weight.size(); k++)
	  pp_res_thread[k]->set(kk, wkk, ww[k]);
      }
      
      // estimate the computational time and update the time count
      time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
      if (tcount && tid==0) { coutCBL << "\r" << float(i)*fact_count << "% completed (" << diff_temp << " seconds)\r"; cout.flush(); }    
      if (i==int(nObj*0.25)) coutCBL << ".............25% completed   " << endl;
      if (i==int(nObj*0.5)) coutCBL << ".............50% completed   " << endl;
      if (i==int(nObj*0.75)) coutCBL << ".............75% completed   "<< endl;   
    }
    
#pragma omp critical
    {
      // sum all the object pairs computed by each thread
      pp->Sum(pp_thread);

      for (size_t i=0; i<pp_res.size(); ++i)
	pp_res[i]->Sum(pp_res_thread[i]);
    }

  }

  
  // show the time spent by the method
  time_t end; time (&end);
  double diff = difftime(end, start);
  if (tid==0) {
    if (diff<60) coutCBL << "   time spent to count the pairs: " << diff << " seconds" << endl ;
    else if (diff<3600) coutCBL << "   time spent to count the pairs: " << diff/60 << " minutes" << endl ;
    else coutCBL << "   time spent to count the pairs: " << diff/3600 << " hours" << endl ;
  }
  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(dp);

}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation::count_pairs_region_test_2D (const std::shared_ptr<Catalogue> cat1, const ChainMesh_Catalogue &ChM, std::shared_ptr<Pair> pp, std::vector<std::shared_ptr<Pair>> pp_res, const std::vector<double> weight, const bool cross, const bool tcount)  
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

    shared_ptr<Pair> pp_thread;
    vector<shared_ptr<Pair> > pp_res_thread(pp_res.size());
    
    pp_thread =  move(Pair::Create(pp->pairType(), pp->pairInfo(), pp->sMin_D1(), pp->sMax_D1(), pp->nbins_D1(), pp->shift_D1(), pp->sMin_D2(), pp->sMax_D2(), pp->nbins_D2(), pp->shift_D2(), pp->angularUnits(), pp->angularWeight()));

    for (size_t i=0; i<pp_res.size(); ++i) 
      pp_res_thread[i] =  move(Pair::Create(pp->pairType(), PairInfo::_standard_, pp->sMin_D1(), pp->sMax_D1(), pp->nbins_D1(), pp->shift_D1(), pp->sMin_D2(), pp->sMax_D2(), pp->nbins_D2(), pp->shift_D2(), pp->angularUnits(), pp->angularWeight()));
    
    // parallelized loop
#pragma omp for schedule(static, 2)
    for (int i=0; i<nObj; ++i) {

      vector<long int> close_objects = ChM.close_objects(cat1->coordinate(i), (cross) ? -1 : i);

      for (auto &&j : close_objects) {      
	int ir, jr;
	double wkk;

	pp_thread->get(cat1->catalogue_object(i), cat2->catalogue_object(j), ir, jr, wkk);
	pp_thread->set(ir, jr, wkk);

	int reg1 = (cross) ? cat1->region(i) : min(cat1->region(i), cat2->region(j));
	int reg2 = (cross) ? cat2->region(j) : max(cat1->region(i), cat2->region(j));

	vector<double> ww = weight;

	ww[reg1] = 0;
	ww[reg2] = 0;

	for(size_t k=0; k< weight.size(); k++)
	  pp_res_thread[k]->set(ir, jr, wkk, ww[k]);
      }
      
      // estimate the computational time and update the time count
      time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
      if (tcount && tid==0) { coutCBL << "\r" << float(i)*fact_count << "% completed (" << diff_temp << " seconds)\r"; cout.flush(); }    
      if (i==int(nObj*0.25)) coutCBL << ".............25% completed   " << endl;
      if (i==int(nObj*0.5)) coutCBL << ".............50% completed   " << endl;
      if (i==int(nObj*0.75)) coutCBL << ".............75% completed   "<< endl;   
    }
    
#pragma omp critical
    {
      // sum all the object pairs computed by each thread
      pp->Sum(pp_thread);

      for (size_t i=0; i<pp_res.size(); ++i)
	pp_res[i]->Sum(pp_res_thread[i]);
    }

  }

  
  // show the time spent by the method
  time_t end; time (&end);
  double diff = difftime(end, start);
  if (tid==0) {
    if (diff<60) coutCBL << "   time spent to count the pairs: " << diff << " seconds" << endl ;
    else if (diff<3600) coutCBL << "   time spent to count the pairs: " << diff/60 << " minutes" << endl ;
    else coutCBL << "   time spent to count the pairs: " << diff/3600 << " hours" << endl ;
  }
  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(dp);

}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation::count_allPairs_region (std::vector<std::shared_ptr<Pair> > &dd_regions, std::vector<std::shared_ptr<Pair> > &rr_regions, std::vector<std::shared_ptr<Pair> > &dr_regions, const TwoPType type, const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const double fact)  
{
  // ----------- compute polar coordinates, if necessary ----------- 

  if (!m_data->isSetVar(Var::_RA_) || !m_data->isSetVar(Var::_Dec_) || !m_data->isSetVar(Var::_Dc_))
    m_data->computePolarCoordinates();
  
  if (!m_random->isSetVar(Var::_RA_) || !m_random->isSetVar(Var::_Dec_) || !m_random->isSetVar(Var::_Dc_))
    m_random->computePolarCoordinates();

  if (type==TwoPType::_angular_) {
    m_data->normalizeComovingCoordinates();
    m_random->normalizeComovingCoordinates();
  }
  

  // ----------- dilute the random catalogue used to compute the RR pairs (to improve the performance) ----------- 

  if (estimator==Estimator::_natural_ && m_random_dilution_fraction!=1.) {
    m_random_dilution_fraction = 1.;
    WarningMsgCBL("m_random_dilution_fraction = 1, since the random catalogue is not diluted when using the natural estimator!", "count_allPairs_region", "TwoPointCorrelation.cpp");
  }

  auto random_dil = make_shared<catalogue::Catalogue>(catalogue::Catalogue(move(m_random->diluted_catalogue(m_random_dilution_fraction))));
  
  
  // ----------- create the chain-mesh ----------- 

  double rMAX;
  double rMIN;

  if (type==TwoPType::_monopole_ || type==TwoPType::_multipoles_direct_ || type==TwoPType::_filtered_) {
    rMAX = m_dd->sMax();
    rMIN = m_dd->sMin();
  }

  else if (type==TwoPType::_angular_) {
    double xx, yy, zz;
    cartesian_coord(radians(m_dd->sMax(), m_dd->angularUnits()), radians(m_dd->sMax(), m_dd->angularUnits()), 1., xx, yy, zz);
    rMAX = sqrt(yy*yy+zz*zz);
    cartesian_coord(radians(m_dd->sMin(), m_dd->angularUnits()), radians(m_dd->sMin(), m_dd->angularUnits()), 1., xx, yy, zz);
    rMIN = sqrt(yy*yy+zz*zz);
  }

  else if (type==TwoPType::_2D_polar_ || type==TwoPType::_multipoles_integrated_ || type ==TwoPType::_wedges_) {
    rMAX = m_dd->sMax_D1();
    rMIN = m_dd->sMin_D1();
  }
  
  else if (type==TwoPType::_2D_Cartesian_ || type==TwoPType::_projected_ || type==TwoPType::_deprojected_) {
    rMAX = max(m_dd->sMax_D1(), m_dd->sMax_D2())*sqrt(2.);
    rMIN = max(m_dd->sMin_D1(), m_dd->sMin_D2())*sqrt(2.);
  }

  else
    ErrorCBL("the chosen two-point correlation function type is uknown!", "count_allPairs_regions", "TwoPointCorrelation.cpp");
  
  
  double cell_size = fact*rMAX; // to be optimized!!!

  ChainMesh_Catalogue ChM_data, ChM_random, ChM_random_dil;

  if (count_dd)
    ChM_data.set_par(cell_size, m_data, rMAX, rMIN);

  if (count_rr || count_dr) 
    ChM_random.set_par(cell_size, m_random, rMAX, rMIN);

  if (count_dr)
    ChM_random_dil.set_par(cell_size, random_dil, rMAX, rMIN);

  
  // ----------- initialize the pair vectors used for resampling ----------- 
  
  const int nRegions = m_random->nRegions();

  if (nRegions<1)
    ErrorCBL("set regions in data and random catalogue to compute 2pcf Jackknife/Bootstrap error!", "count_allPairs_region", "TwoPointCorrelation.cpp");

  const int nP_auto = nRegions*(nRegions+1)*0.5;
  const int nP_cross = nRegions*nRegions;
  
  dd_regions.erase(dd_regions.begin(), dd_regions.end());
  rr_regions.erase(rr_regions.begin(), rr_regions.end());
  dr_regions.erase(dr_regions.begin(), dr_regions.end());

  for (int i=0; i<nP_auto; ++i) {
    
    dd_regions.push_back((m_dd->pairDim()==Dim::_1D_) ? move(Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin(), m_dd->sMax(), m_dd->nbins(), m_dd->shift(), m_dd->angularUnits(), m_dd->angularWeight())) : move(Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2(), m_dd->angularUnits(), m_dd->angularWeight())));
    
    rr_regions.push_back((m_rr->pairDim()==Dim::_1D_) ? move(Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin(), m_rr->sMax(), m_rr->nbins(), m_rr->shift(), m_rr->angularUnits(), m_rr->angularWeight())) : move(Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2(), m_rr->angularUnits(), m_rr->angularWeight())));
  }

  for (int i=0; i<nP_cross; ++i) {

    dr_regions.push_back((m_dr->pairDim()==Dim::_1D_) ? move(Pair::Create(m_dr->pairType(), m_dr->pairInfo(), m_dr->sMin(), m_dr->sMax(), m_dr->nbins(), m_dr->shift(), m_dr->angularUnits(), m_dr->angularWeight())) : move(Pair::Create(m_dr->pairType(), m_dr->pairInfo(), m_dr->sMin_D1(), m_dr->sMax_D1(), m_dr->nbins_D1(), m_dr->shift_D1(), m_dr->sMin_D2(), m_dr->sMax_D2(), m_dr->nbins_D2(), m_dr->shift_D2(), m_dr->angularUnits(), m_dr->angularWeight())));

  }

  // ----------- reset the pair counts --------------------------------------
  resets();

  // ----------- count the number of pairs or read them from file -----------

  string file, file_regions;

  cout << endl; coutCBL << par::col_green << "data-data" << par::col_default << endl;
 
  file = "dd.dat";
  file_regions = "dd_regions.dat";

  if (count_dd) {
    count_pairs_region(m_data, ChM_data, m_dd, dd_regions, false, tcount);
    if (dir_output_pairs!=par::defaultString) {
      write_pairs(m_dd, dir_output_pairs, file);
      write_pairs(dd_regions, dir_output_pairs, file_regions);
    }
  }
  else {
    read_pairs(m_dd, dir_input_pairs, file);
    read_pairs(dd_regions, dir_input_pairs, file_regions);
  }

  cout << endl; coutCBL << par::col_green << "random-random" << par::col_default << endl;

  
  file = "rr.dat";
  file_regions = "rr_regions.dat";

  if (count_rr) {
    count_pairs_region(random_dil, ChM_random_dil, m_rr, rr_regions, false, tcount);
    if (dir_output_pairs!=par::defaultString) {
      write_pairs(m_rr, dir_output_pairs, file);
      write_pairs(rr_regions, dir_output_pairs, file_regions);
    }
  }
  else {
    read_pairs(m_rr, dir_input_pairs, file);
    read_pairs(rr_regions, dir_input_pairs, file_regions);
  }

  if (estimator==Estimator::_LandySzalay_) {

    cout << endl; coutCBL << par::col_green << "data-random" << par::col_default << endl; 

    file = "dr.dat";
    file_regions = "dr_regions.dat";
    
    if (count_dr) {
      count_pairs_region(m_data, ChM_random, m_dr, dr_regions, true, tcount);
      if (dir_output_pairs!=par::defaultString) {
	write_pairs(m_dr, dir_output_pairs, file);
	write_pairs(dr_regions, dir_output_pairs, file_regions);
      }
    }
    else {
      read_pairs(m_dr, dir_input_pairs, file);
      read_pairs(dr_regions, dir_input_pairs, file_regions);
    }

  }
  
  if (count_dd)
    m_data->Order();
  
  if (count_rr || count_dr) 
    m_random->Order();
 
  if (type==TwoPType::_angular_) {
    m_data->restoreComovingCoordinates();
    m_random->restoreComovingCoordinates();
  }

}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation::count_allPairs_region_test (const TwoPType type, const std::vector<double> weight, const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const double fact)  
{
  // ----------- compute polar coordinates, if necessary ----------- 

  if (!m_data->isSetVar(Var::_RA_) || !m_data->isSetVar(Var::_Dec_) || !m_data->isSetVar(Var::_Dc_))
    m_data->computePolarCoordinates();
  
  if (!m_random->isSetVar(Var::_RA_) || !m_random->isSetVar(Var::_Dec_) || !m_random->isSetVar(Var::_Dc_))
    m_random->computePolarCoordinates();

  if (type==TwoPType::_angular_) {
    m_data->normalizeComovingCoordinates();
    m_random->normalizeComovingCoordinates();
  }
  

  // ----------- dilute the random catalogue used to compute the RR pairs (to improve the performance) ----------- 

  if (estimator==Estimator::_natural_ && m_random_dilution_fraction!=1.) {
    m_random_dilution_fraction = 1.;
    WarningMsgCBL("m_random_dilution_fraction = 1, since the random catalogue is not diluted when using the natural estimator!", "count_allPairs_region_test", "TwoPointCorrelation.cpp");
  }

  auto random_dil = make_shared<catalogue::Catalogue>(catalogue::Catalogue(move(m_random->diluted_catalogue(m_random_dilution_fraction))));
  
  
  // ----------- create the chain-mesh ----------- 

  double rMAX;
  double rMIN;

  if (type==TwoPType::_monopole_ || type==TwoPType::_filtered_) {
    rMAX = m_dd->sMax();
    rMIN = m_dd->sMin();
  }

  else if (type==TwoPType::_angular_) {
    double xx, yy, zz;
    cartesian_coord(radians(m_dd->sMax(), m_dd->angularUnits()), radians(m_dd->sMax(), m_dd->angularUnits()), 1., xx, yy, zz);
    rMAX = sqrt(yy*yy+zz*zz);
    cartesian_coord(radians(m_dd->sMin(), m_dd->angularUnits()), radians(m_dd->sMin(), m_dd->angularUnits()), 1., xx, yy, zz);
    rMIN = sqrt(yy*yy+zz*zz);
  }

  else if (type==TwoPType::_2D_polar_ || type==TwoPType::_multipoles_integrated_ || type ==TwoPType::_wedges_) {
    rMAX = m_dd->sMax_D1();
    rMIN = m_dd->sMin_D1();
  }
  
  else if (type==TwoPType::_2D_Cartesian_ || type==TwoPType::_projected_ || type==TwoPType::_deprojected_) {
    rMAX = max(m_dd->sMax_D1(), m_dd->sMax_D2())*sqrt(2.);
    rMIN = max(m_dd->sMin_D1(), m_dd->sMin_D2())*sqrt(2.);
  }
  
  else
    ErrorCBL("the chosen two-point correlation function type is uknown!", "count_allPairs_regions", "TwoPointCorrelation.cpp");
  
  
  double cell_size = fact*rMAX; // to be optimized!!!

  ChainMesh_Catalogue ChM_data, ChM_random, ChM_random_dil;

  if (count_dd)
    ChM_data.set_par(cell_size, m_data, rMAX, rMIN);

  if (count_rr || count_dr) 
    ChM_random.set_par(cell_size, m_random, rMAX, rMIN);

  if (count_dr)
    ChM_random_dil.set_par(cell_size, random_dil, rMAX, rMIN);

  
  // ----------- initialize the pair vectors used for resampling ----------- 

  int nResamplings = weight.size();

  for (int i=0; i<nResamplings; ++i) {

    m_dd_res.push_back((m_dd->pairDim()==Dim::_1D_) ? move(Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin(), m_dd->sMax(), m_dd->nbins(), m_dd->shift(), m_dd->angularUnits(), m_dd->angularWeight())) : move(Pair::Create(m_dd->pairType(), m_dd->pairInfo(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2(), m_dd->angularUnits(), m_dd->angularWeight())));

    m_rr_res.push_back((m_rr->pairDim()==Dim::_1D_) ? move(Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin(), m_rr->sMax(), m_rr->nbins(), m_rr->shift(), m_rr->angularUnits(), m_rr->angularWeight())) : move(Pair::Create(m_rr->pairType(), m_rr->pairInfo(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2(), m_rr->angularUnits(), m_rr->angularWeight())));

    m_dr_res.push_back((m_dr->pairDim()==Dim::_1D_) ? move(Pair::Create(m_dr->pairType(), m_dr->pairInfo(), m_dr->sMin(), m_dr->sMax(), m_dr->nbins(), m_dr->shift(), m_dr->angularUnits(), m_dr->angularWeight())) : move(Pair::Create(m_dr->pairType(), m_dr->pairInfo(), m_dr->sMin_D1(), m_dr->sMax_D1(), m_dr->nbins_D1(), m_dr->shift_D1(), m_dr->sMin_D2(), m_dr->sMax_D2(), m_dr->nbins_D2(), m_dr->shift_D2(), m_dr->angularUnits(), m_dr->angularWeight())));

  }


  // ----------- count the number of pairs or read them from file -----------

  string file, file_regions;

  cout << endl; coutCBL << par::col_green << "data-data" << par::col_default << endl;

  file = "dd.dat";
  file_regions = "dd_regions.dat";

  if (count_dd) {
    count_pairs_region_test(m_data, ChM_data, m_dd, m_dd_res, weight, false, tcount);
    if (dir_output_pairs!=par::defaultString) 
      write_pairs(m_dd, dir_output_pairs, file);
  }
  else 
    read_pairs(m_dd, dir_input_pairs, file);
  
  cout << endl; coutCBL << par::col_green << "random-random" << par::col_default << endl;

  file = "rr.dat";
  file_regions = "rr_regions.dat";

  if (count_rr) {
    count_pairs_region_test(random_dil, ChM_random_dil, m_rr, m_rr_res, weight, false, tcount);
    if (dir_output_pairs!=par::defaultString) 
      write_pairs(m_rr, dir_output_pairs, file);
    
  }
  else 
    read_pairs(m_rr, dir_input_pairs, file);
  

  if (estimator==Estimator::_LandySzalay_) {

    cout << endl; coutCBL << par::col_green << "data-random" << par::col_default << endl; 

    file = "dr.dat";
    file_regions = "dr_regions.dat";
    
    if (count_dr) {
      count_pairs_region_test(m_data, ChM_random, m_dr, m_dr_res, weight, true, tcount);
      if (dir_output_pairs!=par::defaultString) 
	write_pairs(m_dr, dir_output_pairs, file);
    }
    else 
      read_pairs(m_dr, dir_input_pairs, file);

  }
  
  if (count_dd)
    m_data->Order();
  
  if (count_rr || count_dr) 
    m_random->Order();
 
  if (type==TwoPType::_angular_) {
    m_data->restoreComovingCoordinates();
    m_random->restoreComovingCoordinates();
  }
  
}
