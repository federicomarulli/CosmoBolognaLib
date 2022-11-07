/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli, Michele Moresco         *
 *  and Alfonso Veropalumbo                                         *
 *                                                                  *
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
 *  @file
 *  CatalogueAnalysis/ThreePointCorrelation/ThreePointCorrelation.cpp
 *
 *  @brief Methods of the class ThreePointCorrelation used to measure
 *  the three-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  ThreePointCorrelation used to measure the three-point correlation
 *  function
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#include "ThreePointCorrelation.h"
#include "ThreePointCorrelation_angular_connected.h"
#include "ThreePointCorrelation_angular_reduced.h"
#include "ThreePointCorrelation_comoving_connected.h"
#include "ThreePointCorrelation_comoving_reduced.h"
#include "ThreePointCorrelation_comoving_multipoles_single.h"
#include "ThreePointCorrelation_comoving_multipoles_all.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace chainmesh;
using namespace triplets;
using namespace measure::threept;


// ============================================================================


std::shared_ptr<ThreePointCorrelation> cbl::measure::threept::ThreePointCorrelation::Create (const ThreePType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const triplets::TripletType tripletType, const double side_s, const double side_u, const double perc_increase, const int nbins)
{
  if (type==ThreePType::_angular_connected_) return move(unique_ptr<ThreePointCorrelation_angular_connected>(new ThreePointCorrelation_angular_connected(data, random, side_s, side_u, perc_increase, nbins)));
 
  if (type==ThreePType::_angular_reduced_) return move(unique_ptr<ThreePointCorrelation_angular_reduced>(new ThreePointCorrelation_angular_reduced(data, random, side_s, side_u, perc_increase, nbins)));
 
  if (type==ThreePType::_comoving_connected_) return move(unique_ptr<ThreePointCorrelation_comoving_connected>(new ThreePointCorrelation_comoving_connected(data, random, tripletType, side_s, side_u, perc_increase, nbins)));
 
  if (type==ThreePType::_comoving_reduced_) return move(unique_ptr<ThreePointCorrelation_comoving_reduced>(new ThreePointCorrelation_comoving_reduced(data, random, tripletType, side_s, side_u, perc_increase, nbins)));
 
  else ErrorCBL("no such type of object!", "Create", "ThreePointCorrelation.cpp");
  
  return NULL;
}


// ============================================================================


std::shared_ptr<ThreePointCorrelation> cbl::measure::threept::ThreePointCorrelation::Create (const ThreePType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const triplets::TripletType tripletType, const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int nbins)
{
  if (type==ThreePType::_angular_connected_) return move(unique_ptr<ThreePointCorrelation_angular_connected>(new ThreePointCorrelation_angular_connected(data, random, r12, r12_binSize, r13, r13_binSize, nbins)));
 
  if (type==ThreePType::_angular_reduced_) return move(unique_ptr<ThreePointCorrelation_angular_reduced>(new ThreePointCorrelation_angular_reduced(data, random, r12, r12_binSize, r13, r13_binSize, nbins)));
 
  if (type==ThreePType::_comoving_connected_) return move(unique_ptr<ThreePointCorrelation_comoving_connected>(new ThreePointCorrelation_comoving_connected(data, random, tripletType, r12, r12_binSize, r13, r13_binSize, nbins)));
 
  if (type==ThreePType::_comoving_reduced_) return move(unique_ptr<ThreePointCorrelation_comoving_reduced>(new ThreePointCorrelation_comoving_reduced(data, random, tripletType, r12, r12_binSize, r13, r13_binSize, nbins)));
 
  else ErrorCBL("no such type of object!", "Create", "ThreePointCorrelation.cpp");
  
  return NULL;
}


// ============================================================================


std::shared_ptr<ThreePointCorrelation> cbl::measure::threept::ThreePointCorrelation::Create (const catalogue::Catalogue data, const catalogue::Catalogue random, const double r12Min, const double r12Max, const double r13Min, const double r13Max, const int nOrders, const double split, const int seed)
{
  return move(unique_ptr<ThreePointCorrelation_comoving_multipoles_single>(new ThreePointCorrelation_comoving_multipoles_single(data, random, r12Min, r12Max, r13Min, r13Max, nOrders, split, seed)));
}


// ============================================================================


std::shared_ptr<ThreePointCorrelation> cbl::measure::threept::ThreePointCorrelation::Create (const catalogue::Catalogue data, const catalogue::Catalogue random, const double rMin, const double rMax, const double binSize, const int nOrders, const double split, const int seed)
{
  return move(unique_ptr<ThreePointCorrelation_comoving_multipoles_all>(new ThreePointCorrelation_comoving_multipoles_all(data, random, rMin, rMax, binSize, nOrders, split, seed)));
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation::count_triplets (const std::shared_ptr<Catalogue> cat1, const ChainMesh_Catalogue &ChainMesh_rMAX1, const ChainMesh_Catalogue &ChainMesh_rMAX2, std::shared_ptr<Triplet> tt, const bool tcount) 
{
  time_t start; time (&start);
  
  int nObj = cat1->nObjects();

  float fact_count = 100./nObj;
  
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);

  shared_ptr<Catalogue> cat2 = ChainMesh_rMAX1.catalogue();
  shared_ptr<Catalogue> cat3 = ChainMesh_rMAX2.catalogue();

  double r12_min = tt->r12()-0.5*tt->r12_binSize();
  double r12_max = tt->r12()+0.5*tt->r12_binSize();
  double r13_min = tt->r13()-0.5*tt->r13_binSize();
  double r13_max = tt->r13()+0.5*tt->r13_binSize();

  int tid = 0;
#pragma omp parallel private(tid)
  {
    tid = omp_get_thread_num();
    
    // if (tid == 0) coutCBL << "Number of threads = " << omp_get_num_threads() << endl;

    // internal object used by each thread to handle triplets
    shared_ptr<Triplet> tt_thread = move(Triplet::Create(tt->tripletType(), tt->r12(), tt->r12_binSize(), tt->r13(), tt->r13_binSize(), tt->nbins()));
    
#pragma omp for schedule(static, 2)  
    for (int i=0; i<nObj; i++) { // loop on the objects of the catalogue
    
      // get the indexes of the objects at r12
      vector<long> close_objects12 = ChainMesh_rMAX1.close_objects(cat1->coordinate(i), -1);

      // get the indexes of objects at r13
      vector<long> close_objects13 = ChainMesh_rMAX2.close_objects(cat1->coordinate(i), -1);

      vector<double> x2, y2, z2, r12, w2, x3, y3, z3, r13, w3;

      for (auto &&j : close_objects12) { // loop on the objects at r12 

	double rr = cat1->distance(i, cat2->catalogue_object(j));
	if (r12_min<rr && rr<r12_max){
	  x2.push_back(cat2->xx(j));
	  y2.push_back(cat2->yy(j));
	  z2.push_back(cat2->zz(j));
	  r12.push_back(rr);
	  w2.push_back(cat2->weight(j));
	}
      }

      for (auto &&k : close_objects13) { // loop on the objects at r13
	double rr = cat1->distance(i, cat3->catalogue_object(k));
	if (r13_min<rr && rr<r13_max){ 
	  x3.push_back(cat3->xx(k));
	  y3.push_back(cat3->yy(k));
	  z3.push_back(cat3->zz(k));
	  r13.push_back(rr);
	  w3.push_back(cat3->weight(k));
	}
      }

      for (size_t j=0; j<r12.size(); j++)
	for (size_t k=0; k<r13.size(); k++){
	  double r23 = sqrt( pow(x2[j]-x3[k],2)+pow(y2[j]-y3[k],2)+pow(z2[j]-z3[k],2));
	  double ww = cat1->weight(i)*w2[j]*w3[k];
	  tt_thread->put(r12[j], r13[k], r23, ww);
	}

      time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
      if (tcount && tid==0) { coutCBL <<"\r..."<<float(i)*fact_count<<"% completed  ("<<diff_temp/60<<" minutes)\r"; cout.flush(); }
      if (i==int(nObj*0.25)) coutCBL <<".....25% completed   "<<endl;
      if (i==int(nObj*0.5)) coutCBL <<".....50% completed   "<<endl;
      if (i==int(nObj*0.75)) coutCBL <<".....75% completed   "<<endl;
    }
    
#pragma omp critical
    {
      // sum all the object triplets computed by each thread
      tt->Sum(tt_thread);
    }
    
  }

  
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<3600) coutCBL <<"   time spent to count the triplets: "<<diff/60<<" minutes"<<endl<<endl;
  else coutCBL <<"   time spent to count the triplets: "<<diff/3600<<" hours"<<endl<<endl;

  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(6); 
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation::count_allTriplets (const std::string dir_output_triplets, const std::vector<std::string> dir_input_triplets, const bool count_ddd, const bool count_rrr, const bool count_ddr, const bool count_drr, const bool tcount, const double fact)
{
  // ----- double chain-mesh -----
  
  double rMAX1 = m_ddd->r12()+m_ddd->r12_binSize();
  double rMAX2 = m_ddd->r13()+m_ddd->r13_binSize();

  double cell_size1 = max(5., rMAX1*fact);
  double cell_size2 = max(5., rMAX2*fact);

  ChainMesh_Catalogue ChainMesh_data_rMAX1, ChainMesh_data_rMAX2, ChainMesh_random_rMAX1, ChainMesh_random_rMAX2;

  if (count_ddd || count_ddr || count_drr) {
    auto data_rMAX1 = make_shared<catalogue::Catalogue>(*m_data);
    auto data_rMAX2 = make_shared<catalogue::Catalogue>(*m_data);
    ChainMesh_data_rMAX1.set_par(cell_size1, data_rMAX1, rMAX1);
    ChainMesh_data_rMAX2.set_par(cell_size2, data_rMAX2, rMAX2);
  }
  if (count_rrr || count_ddr || count_drr) {
    auto random_rMAX1 = make_shared<catalogue::Catalogue>(*m_random);
    auto random_rMAX2 = make_shared<catalogue::Catalogue>(*m_random);
    ChainMesh_random_rMAX1.set_par(cell_size1, random_rMAX1, rMAX1);
    ChainMesh_random_rMAX2.set_par(cell_size2, random_rMAX2, rMAX2);
  }
  
  // ----------- count the number of triplets ----------- 

  string file;
  
  coutCBL << par::col_green << "data-data-data" << par::col_default << endl;
  
  file = "ddd.dat";
  
  if (count_ddd) { 
    count_triplets(m_data, ChainMesh_data_rMAX1, ChainMesh_data_rMAX2, m_ddd, tcount);
    if (dir_output_triplets!=par::defaultString) write_triplets(m_ddd, dir_output_triplets, file);
  } 
  else read_triplets (m_ddd, dir_input_triplets, file);


  coutCBL << par::col_green << "random-random-random" << par::col_default << endl;
  
  file = "rrr.dat";
  
  if (count_rrr==1) {
    count_triplets(m_random, ChainMesh_random_rMAX1, ChainMesh_random_rMAX2, m_rrr, tcount);
    if (dir_output_triplets!=par::defaultString) write_triplets(m_rrr, dir_output_triplets, file);
  } 
  else if (count_rrr==0) read_triplets (m_rrr, dir_input_triplets, file);
 

  coutCBL << par::col_green << "data-data-random" << par::col_default << endl;
  
  file = "ddr.dat";
  
  if (count_ddr) {

    shared_ptr<Triplet> ddr1 = move(Triplet::Create(m_ddr->tripletType(), m_ddr->r12(), m_ddr->r12_binSize(), m_ddr->r13(), m_ddr->r13_binSize(), m_ddr->nbins()));
    shared_ptr<Triplet> ddr2 = move(Triplet::Create(m_ddr->tripletType(), m_ddr->r12(), m_ddr->r12_binSize(), m_ddr->r13(), m_ddr->r13_binSize(), m_ddr->nbins()));
    shared_ptr<Triplet> ddr3 = move(Triplet::Create(m_ddr->tripletType(), m_ddr->r12(), m_ddr->r12_binSize(), m_ddr->r13(), m_ddr->r13_binSize(), m_ddr->nbins()));

    count_triplets(m_data, ChainMesh_data_rMAX1, ChainMesh_random_rMAX2, ddr1, tcount);
    count_triplets(m_data, ChainMesh_random_rMAX1, ChainMesh_data_rMAX2, ddr2, tcount);
    count_triplets(m_random, ChainMesh_data_rMAX1, ChainMesh_data_rMAX2, ddr3, tcount);

    m_ddr->Sum(ddr1); m_ddr->Sum(ddr2); m_ddr->Sum(ddr3); 
   
    if (dir_output_triplets!=par::defaultString) write_triplets (m_ddr, dir_output_triplets, file);
  } 

  else read_triplets(m_ddr, dir_input_triplets, file);

  
  coutCBL << par::col_green << "data-random-random" << par::col_default << endl;
  
  file = "drr.dat";
  
  if (count_drr) {

    shared_ptr<Triplet> drr1 = move(Triplet::Create(m_drr->tripletType(), m_drr->r12(), m_drr->r12_binSize(), m_drr->r13(), m_drr->r13_binSize(), m_drr->nbins()));
    shared_ptr<Triplet> drr2 = move(Triplet::Create(m_drr->tripletType(), m_drr->r12(), m_drr->r12_binSize(), m_drr->r13(), m_drr->r13_binSize(), m_drr->nbins()));
    shared_ptr<Triplet> drr3 = move(Triplet::Create(m_drr->tripletType(), m_drr->r12(), m_drr->r12_binSize(), m_drr->r13(), m_drr->r13_binSize(), m_drr->nbins()));
    
    count_triplets(m_random, ChainMesh_random_rMAX1, ChainMesh_data_rMAX2, drr1, tcount);
    count_triplets(m_random, ChainMesh_data_rMAX1, ChainMesh_random_rMAX2, drr2, tcount);
    count_triplets(m_data, ChainMesh_random_rMAX1, ChainMesh_random_rMAX2, drr3, tcount);

    m_drr->Sum(drr1); m_drr->Sum(drr2); m_drr->Sum(drr3);
    
    if (dir_output_triplets!=par::defaultString) write_triplets(m_drr, dir_output_triplets, file);
  } 

  else 
    read_triplets (m_drr, dir_input_triplets, file);
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation::count_triplets_region (const std::shared_ptr<Catalogue> cat1, const ChainMesh_Catalogue &ChainMesh_rMAX1, const ChainMesh_Catalogue &ChainMesh_rMAX2, std::shared_ptr<Triplet> tt, std::vector<std::shared_ptr<Triplet>> tt_region, const std::vector<std::vector<double>> weight, const bool tcount) 
{
  time_t start; time (&start);
  
  int nObj = cat1->nObjects();

  float fact_count = 100./nObj;
  
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);

  shared_ptr<Catalogue> cat2 = ChainMesh_rMAX1.catalogue();
  shared_ptr<Catalogue> cat3 = ChainMesh_rMAX2.catalogue();

  double r12_min = tt->r12()-0.5*tt->r12_binSize();
  double r12_max = tt->r12()+0.5*tt->r12_binSize();
  double r13_min = tt->r13()-0.5*tt->r13_binSize();
  double r13_max = tt->r13()+0.5*tt->r13_binSize();

  int tid = 0;
#pragma omp parallel private(tid)
  {
    tid = omp_get_thread_num();
    
    // if (tid == 0) coutCBL << "Number of threads = " << omp_get_num_threads() << endl;

    // internal object used by each thread to handle triplets
    shared_ptr<Triplet> tt_thread = move(Triplet::Create(tt->tripletType(), tt->r12(), tt->r12_binSize(), tt->r13(), tt->r13_binSize(), tt->nbins()));
    vector<shared_ptr<Triplet> > tt_region_thread(tt_region.size());

    for (size_t i=0; i<tt_region.size(); ++i) 
      tt_region_thread[i] = move(Triplet::Create(tt->tripletType(), tt->r12(), tt->r12_binSize(), tt->r13(), tt->r13_binSize() , tt->nbins()));

#pragma omp for schedule(static, 2)  
    for (int i=0; i<nObj; i++) { // loop on the objects of the catalogue

      int reg1 = cat1->region(i);
    
      // get the indexes of the objects at r12
      vector<long> close_objects12 = ChainMesh_rMAX1.close_objects(cat1->coordinate(i), -1);

      // get the indexes of objects at r13
      vector<long> close_objects13 = ChainMesh_rMAX2.close_objects(cat1->coordinate(i), -1);

      vector<double> x2, y2, z2, r12, w2, x3, y3, z3, r13, w3;
      vector<long> reg2, reg3;

      for (auto &&j : close_objects12) { // loop on the objects at r12 
	double rr = cat1->distance(i, cat2->catalogue_object(j));

	if (r12_min<rr && rr<r12_max) {
	  x2.push_back(cat2->xx(j));
	  y2.push_back(cat2->yy(j));
	  z2.push_back(cat2->zz(j));
	  r12.push_back(rr);
	  w2.push_back(cat2->weight(j));
	  reg2.push_back(cat2->region(j));
	}
      }

      for (auto &&k : close_objects13) { // loop on the objects at r12 
	double rr = cat1->distance(i, cat3->catalogue_object(k));

	if (r13_min<rr && rr<r13_max) {
	  x3.push_back(cat3->xx(k));
	  y3.push_back(cat3->yy(k));
	  z3.push_back(cat3->zz(k));
	  r13.push_back(rr);
	  w3.push_back(cat3->weight(k));
	  reg3.push_back(cat3->region(k));
	}
      }

      for (size_t j=0; j<r12.size(); j++)
	for (size_t k=0; k<r13.size(); k++){
	  int klin;
	  double r23 = sqrt( pow(x2[j]-x3[k],2)+pow(y2[j]-y3[k],2)+pow(z2[j]-z3[k],2));
	  double ww = cat1->weight(i)*w2[j]*w3[k];
	  tt_thread->get_triplet(r12[j], r13[k], r23, klin);
	  tt_thread->set_triplet(klin, ww);
	  for (size_t r=0; r<weight.size(); r++) {
	    double region_weight = weight[r][reg1]*weight[r][reg2[j]]*weight[r][reg3[k]];
	    tt_region_thread[r]->set_triplet(klin, ww*region_weight);
	  }
	}
      time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
      if (tcount && tid==0) { coutCBL <<"\r..."<<float(i)*fact_count<<"% completed  ("<<diff_temp/60<<" minutes)\r"; cout.flush(); }
      if (i==int(nObj*0.25)) coutCBL <<".....25% completed   "<<endl;
      if (i==int(nObj*0.5)) coutCBL <<".....50% completed   "<<endl;
      if (i==int(nObj*0.75)) coutCBL <<".....75% completed   "<<endl;
    }
    
#pragma omp critical
    {
      // sum all the object triplets computed by each thread
      tt->Sum(tt_thread);

      for (size_t k=0; k< weight.size(); k++)
	tt_region[k]->Sum(tt_region_thread[k]);

    }
    
  }

  
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<3600) coutCBL <<"   time spent to count the triplets: "<<diff/60<<" minutes"<<endl<<endl;
  else coutCBL <<"   time spent to count the triplets: "<<diff/3600<<" hours"<<endl<<endl;

  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(6); 
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation::count_allTriplets_region (const std::vector<std::vector<double>> weight, const std::string dir_output_triplets, const std::vector<std::string> dir_input_triplets, const bool count_ddd, const bool count_rrr, const bool count_ddr, const bool count_drr, const bool tcount, const double fact)
{
  // ----- double chain-mesh -----

  double rMAX1 = m_ddd->r12()+m_ddd->r12_binSize();
  double rMAX2 = m_ddd->r13()+m_ddd->r13_binSize();

  double cell_size1 = max(5., rMAX1*fact);
  double cell_size2 = max(5., rMAX2*fact);

  ChainMesh_Catalogue ChainMesh_data_rMAX1, ChainMesh_data_rMAX2, ChainMesh_random_rMAX1, ChainMesh_random_rMAX2;

  if (count_ddd || count_ddr || count_drr) {
    auto data_rMAX1 = make_shared<catalogue::Catalogue>(*m_data);
    auto data_rMAX2 = make_shared<catalogue::Catalogue>(*m_data);
    ChainMesh_data_rMAX1.set_par(cell_size1, data_rMAX1, rMAX1);
    ChainMesh_data_rMAX2.set_par(cell_size2, data_rMAX2, rMAX2);
  }
  if (count_rrr || count_ddr || count_drr) {
    auto random_rMAX1 = make_shared<catalogue::Catalogue>(*m_random);
    auto random_rMAX2 = make_shared<catalogue::Catalogue>(*m_random);
    ChainMesh_random_rMAX1.set_par(cell_size1, random_rMAX1, rMAX1);
    ChainMesh_random_rMAX2.set_par(cell_size2, random_rMAX2, rMAX2);
  }

  // ----------- initialize the triplet vectors used for resampling ----------- 

  m_ddd_regions.erase(m_ddd_regions.begin(), m_ddd_regions.end());
  m_rrr_regions.erase(m_rrr_regions.begin(), m_rrr_regions.end());
  m_ddr_regions.erase(m_ddr_regions.begin(), m_ddr_regions.end());
  m_drr_regions.erase(m_drr_regions.begin(), m_drr_regions.end());

  int nResamplings = weight.size();

  for (int i=0; i<nResamplings; ++i) {
    m_ddd_regions.push_back(move(Triplet::Create(m_ddd->tripletType(), m_ddd->r12(), m_ddd->r12_binSize(),  m_ddd->r13(), m_ddd->r13_binSize(), m_ddd->nbins())));
    m_rrr_regions.push_back(move(Triplet::Create(m_rrr->tripletType(), m_rrr->r12(), m_rrr->r12_binSize(),  m_rrr->r13(), m_rrr->r13_binSize(), m_rrr->nbins())));
    m_ddr_regions.push_back(move(Triplet::Create(m_ddr->tripletType(), m_ddr->r12(), m_ddr->r12_binSize(),  m_ddr->r13(), m_ddr->r13_binSize(), m_ddr->nbins())));
    m_drr_regions.push_back(move(Triplet::Create(m_drr->tripletType(), m_drr->r12(), m_drr->r12_binSize(),  m_drr->r13(), m_drr->r13_binSize(), m_drr->nbins())));
  }
  
  // ----------- count the number of triplets ----------- 

  string file;
  
  coutCBL << par::col_green << "data-data-data" << par::col_default << endl;
  
  file = "ddd.dat";
  
  if (count_ddd) { 
    count_triplets_region(m_data, ChainMesh_data_rMAX1, ChainMesh_data_rMAX2, m_ddd, m_ddd_regions, weight, tcount);
    if (dir_output_triplets!=par::defaultString) write_triplets(m_ddd, dir_output_triplets, file);
  } 
  else read_triplets (m_ddd, dir_input_triplets, file);


  coutCBL << par::col_green << "random-random-random" << par::col_default << endl;
  
  file = "rrr.dat";
  
  if (count_rrr==1) {
    count_triplets_region(m_random, ChainMesh_random_rMAX1, ChainMesh_random_rMAX2, m_rrr, m_rrr_regions, weight, tcount);
    if (dir_output_triplets!=par::defaultString) write_triplets(m_rrr, dir_output_triplets, file);
  } 
  else if (count_rrr==0) read_triplets (m_rrr, dir_input_triplets, file);
 

  coutCBL << par::col_green << "data-data-random" << par::col_default << endl;
  
  file = "ddr.dat";
  
  if (count_ddr) {

    shared_ptr<Triplet> ddr1 = move(Triplet::Create(m_ddr->tripletType(), m_ddr->r12(), m_ddr->r12_binSize(),  m_ddr->r13(), m_ddr->r13_binSize(), m_ddr->nbins()));
    shared_ptr<Triplet> ddr2 = move(Triplet::Create(m_ddr->tripletType(), m_ddr->r12(), m_ddr->r12_binSize(),  m_ddr->r13(), m_ddr->r13_binSize(), m_ddr->nbins()));
    shared_ptr<Triplet> ddr3 = move(Triplet::Create(m_ddr->tripletType(), m_ddr->r12(), m_ddr->r12_binSize(),  m_ddr->r13(), m_ddr->r13_binSize(), m_ddr->nbins()));

    vector<shared_ptr<Triplet>> ddr1_regions;
    vector<shared_ptr<Triplet>> ddr2_regions;
    vector<shared_ptr<Triplet>> ddr3_regions;

    for (int i=0; i<nResamplings; i++) {
      ddr1_regions.push_back(move(Triplet::Create(m_ddr->tripletType(), m_ddr->r12(), m_ddr->r12_binSize(),  m_ddr->r13(), m_ddr->r13_binSize(), m_ddr->nbins())));
      ddr2_regions.push_back(move(Triplet::Create(m_ddr->tripletType(), m_ddr->r12(), m_ddr->r12_binSize(),  m_ddr->r13(), m_ddr->r13_binSize(), m_ddr->nbins())));
      ddr3_regions.push_back(move(Triplet::Create(m_ddr->tripletType(), m_ddr->r12(), m_ddr->r12_binSize(),  m_ddr->r13(), m_ddr->r13_binSize(), m_ddr->nbins())));
    }
    
    count_triplets_region(m_data, ChainMesh_data_rMAX1, ChainMesh_random_rMAX2, ddr1, ddr1_regions, weight, tcount);
    count_triplets_region(m_data, ChainMesh_random_rMAX1, ChainMesh_data_rMAX2, ddr2, ddr2_regions, weight,tcount);
    count_triplets_region(m_random, ChainMesh_data_rMAX1, ChainMesh_data_rMAX2, ddr3, ddr3_regions, weight, tcount);

    m_ddr->Sum(ddr1); m_ddr->Sum(ddr2); m_ddr->Sum(ddr3); 

    for (int i=0; i<nResamplings; i++) {
      m_ddr_regions[i]->Sum(ddr1_regions[i]);
      m_ddr_regions[i]->Sum(ddr2_regions[i]);
      m_ddr_regions[i]->Sum(ddr3_regions[i]);
    }
   
    if (dir_output_triplets!=par::defaultString) write_triplets (m_ddr, dir_output_triplets, file);
  } 

  else read_triplets(m_ddr, dir_input_triplets, file);

  
  coutCBL << par::col_green << "data-random-random" << par::col_default << endl;
  
  file = "drr.dat";
  
  if (count_drr) {

    shared_ptr<Triplet> drr1 = move(Triplet::Create(m_drr->tripletType(), m_drr->r12(), m_drr->r12_binSize(),  m_drr->r13(), m_drr->r13_binSize(), m_drr->nbins()));
    shared_ptr<Triplet> drr2 = move(Triplet::Create(m_drr->tripletType(), m_drr->r12(), m_drr->r12_binSize(),  m_drr->r13(), m_drr->r13_binSize(), m_drr->nbins()));
    shared_ptr<Triplet> drr3 = move(Triplet::Create(m_drr->tripletType(), m_drr->r12(), m_drr->r12_binSize(),  m_drr->r13(), m_drr->r13_binSize(), m_drr->nbins()));
    
    vector<shared_ptr<Triplet>> drr1_regions;
    vector<shared_ptr<Triplet>> drr2_regions;
    vector<shared_ptr<Triplet>> drr3_regions;

    for (int i=0; i<nResamplings; i++) {
      drr1_regions.push_back(move(Triplet::Create(m_drr->tripletType(), m_drr->r12(), m_drr->r12_binSize(),  m_drr->r13(), m_drr->r13_binSize(), m_drr->nbins())));
      drr2_regions.push_back(move(Triplet::Create(m_drr->tripletType(), m_drr->r12(), m_drr->r12_binSize(),  m_drr->r13(), m_drr->r13_binSize(), m_drr->nbins())));
      drr3_regions.push_back(move(Triplet::Create(m_drr->tripletType(), m_drr->r12(), m_drr->r12_binSize(),  m_drr->r13(), m_drr->r13_binSize(), m_drr->nbins())));
    }
    
    count_triplets_region(m_random, ChainMesh_random_rMAX1, ChainMesh_data_rMAX2, drr1, drr1_regions, weight, tcount);
    count_triplets_region(m_random, ChainMesh_data_rMAX1, ChainMesh_random_rMAX2, drr2, drr1_regions, weight, tcount);
    count_triplets_region(m_data, ChainMesh_random_rMAX1, ChainMesh_random_rMAX2, drr3, drr1_regions, weight, tcount);

    m_drr->Sum(drr1); m_drr->Sum(drr2); m_drr->Sum(drr3); 

    for (int i=0; i<nResamplings; i++) {
      m_drr_regions[i]->Sum(drr1_regions[i]);
      m_drr_regions[i]->Sum(drr2_regions[i]);
      m_drr_regions[i]->Sum(drr3_regions[i]);
    }
    
    if (dir_output_triplets!=par::defaultString) write_triplets(m_drr, dir_output_triplets, file);
  } 

  else 
    read_triplets (m_drr, dir_input_triplets, file);
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation::write_triplets (std::shared_ptr<triplets::Triplet> TT, const std::string dir, const std::string file) const
{  
  string MK = "mkdir -p "+dir;
  if (system (MK.c_str())) {}
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  for (int i=0; i<TT->nbins(); i++) 
    fout << TT->TT1D(i) << endl;
  
  fout.clear(); fout.close(); coutCBL << "I wrote the file " << file << endl << endl;
}



// ============================================================================


void cbl::measure::threept::ThreePointCorrelation::read_triplets (std::shared_ptr<triplets::Triplet> TT, const std::vector<std::string> dir, const std::string file) 
{
  if (dir.size()==0)
    ErrorCBL("dir.size()=0!", "read_triplets", "TwoPointCorrelation.cpp");

  for (size_t dd=0; dd<dir.size(); dd++) {
        
    string file_in = dir[dd]+file; 
    coutCBL << "I'm reading the triplet file: " << file_in << endl;
    
    ifstream fin(file_in.c_str()); checkIO(fin, file_in);
   
    double pp;
    for (int i=0; i<TT->nbins(); i++) {
      fin >>pp;
      TT->add_TT1D(i, pp);
    }
    
    fin.clear(); fin.close(); coutCBL << "I read the file " << file_in << endl;
  }
}
 
