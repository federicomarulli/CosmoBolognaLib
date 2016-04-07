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
 *  @file CatalogueAnalysis/ThreePointCorrelation/ThreePointCorrelation.cpp
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
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#include "ThreePointCorrelation_angular_connected.h"
#include "ThreePointCorrelation_angular_reduced.h"
#include "ThreePointCorrelation_comoving_connected.h"
#include "ThreePointCorrelation_comoving_reduced.h"

using namespace cosmobl;
using namespace catalogue;
using namespace triplets;
using namespace threept;


// ============================================================================


shared_ptr<ThreePointCorrelation> cosmobl::threept::ThreePointCorrelation::Create (const ThreePType type, const Catalogue data, const Catalogue random, const TripletType tripletType, const double side_s, const double side_u, const double perc_increase, const int nbins)
{
  if (type==_angular_connected_) return move(unique_ptr<ThreePointCorrelation_angular_connected>(new ThreePointCorrelation_angular_connected(data, random, side_s, side_u, perc_increase, nbins)));
 
  if (type==_angular_reduced_) return move(unique_ptr<ThreePointCorrelation_angular_reduced>(new ThreePointCorrelation_angular_reduced(data, random, side_s, side_u, perc_increase, nbins)));
 
  if (type==_comoving_connected_) return move(unique_ptr<ThreePointCorrelation_comoving_connected>(new ThreePointCorrelation_comoving_connected(data, random, tripletType, side_s, side_u, perc_increase, nbins)));
 
  if (type==_comoving_reduced_) return move(unique_ptr<ThreePointCorrelation_comoving_reduced>(new ThreePointCorrelation_comoving_reduced(data, random, tripletType, side_s, side_u, perc_increase, nbins)));
 
  else ErrorMsg("Error in cosmobl::threept::ThreePointCorrelation::Create of ThreePointCorrelation.cpp: no such type of object!");
  
  return NULL;
}


// ============================================================================


void cosmobl::threept::ThreePointCorrelation::count_triplets (const shared_ptr<Catalogue> cat1, const ChainMesh_Catalogue &ChainMesh_rMAX1, const ChainMesh_Catalogue &ChainMesh_rMAX2, shared_ptr<Triplet> tt, const bool tcount) 
{
  time_t start; time (&start);
  
  int nObj = cat1->nObjects();

  float fact_count = 100./nObj;
  
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);

  shared_ptr<Catalogue> cat2 = ChainMesh_rMAX1.catalogue();
  shared_ptr<Catalogue> cat3 = ChainMesh_rMAX2.catalogue();

  int tid = 0;
#pragma omp parallel private(tid)
  {
    tid = omp_get_thread_num();
    
    // if (tid == 0) cout << "Number of threads = " << omp_get_num_threads() << endl;

    // internal object used by each thread to handle triplets
    shared_ptr<Triplet> tt_thread = move(Triplet::Create(tt->tripletType(), tt->side_s(), tt->side_u(), tt->perc_increase(), tt->nbins()));
    
#pragma omp for schedule(static, 2)  
    for (int i=0; i<nObj; i++) { // loop on the objects of the catalogue
    
      // get the indexes of the objects at r12
      vector<long> close_objects12 = ChainMesh_rMAX1.close_objects(cat1->coordinates(i), -1);
    
      for (auto &&j : close_objects12) { // loop on the objects at r12 

	double r12 = cat1->distance(i, cat2->catalogue_object(j));
      
	if (tt_thread->side_s()*(1-tt_thread->perc_increase())<r12 && r12<tt_thread->side_s()*(1+tt_thread->perc_increase())) {

	  // get the indexes of objects at r13
	  vector<long> close_objects13 = ChainMesh_rMAX2.close_objects(cat1->coordinates(i), -1);
	
	  for (auto &&k : close_objects13) { // loop on the objects at r13

	    double r13 = cat1->distance(i, cat3->catalogue_object(k));
	  
	    if (tt_thread->side_u()*tt_thread->side_s()*(1-tt_thread->perc_increase())<r13 && r13<tt_thread->side_u()*tt_thread->side_s()*(1+tt_thread->perc_increase())) {

	      double r23 = cat2->distance(j, cat3->catalogue_object(k));
	    
	      double ww = cat1->weight(i)*cat2->weight(j)*cat3->weight(k);
	    
	      tt->put(r12, r13, r23, ww);
	    }
	  }
	}
      }

      time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
      if (tcount && tid==0) { cout <<"\r..."<<float(i)*fact_count<<"% completed  ("<<diff_temp/60<<" minutes)\r"; cout.flush(); }
      if (i==int(nObj*0.25)) cout <<".....25% completed"<<endl;
      if (i==int(nObj*0.5)) cout <<".....50% completed"<<endl;
      if (i==int(nObj*0.75)) cout <<".....75% completed"<<endl;
    }
    
#pragma omp critical
    {
      // sum all the object triplets computed by each thread
      tt->Sum(tt_thread);
    }
    
  }

  
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<3600) cout <<"   time spent to count the triplets: "<<diff/60<<" minutes"<<endl<<endl;
  else cout <<"   time spent to count the triplets: "<<diff/3600<<" hours"<<endl<<endl;

  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(6); 
}


// ============================================================================


void cosmobl::threept::ThreePointCorrelation::count_allTriplets (const string dir_output_triplets, const vector<string> dir_input_triplets, const int count_ddd, const int count_rrr, const int count_ddr, const int count_drr, const bool tcount)
{
  // ----- double chain-mesh -----
  
  double rMAX1 = m_ddd->side_s()*(1.+2.*m_ddd->perc_increase());
  double rMAX2 = m_ddd->side_u()*m_ddd->side_s()*(1.+2.*m_ddd->perc_increase());
 
  double cell_size = rMAX2*0.1;
  
  ChainMesh_Catalogue ChainMesh_data_rMAX1, ChainMesh_data_rMAX2, ChainMesh_random_rMAX1, ChainMesh_random_rMAX2;
  
  if (count_ddd==1 || count_ddr==1) {
    ChainMesh_data_rMAX1.set_par(cell_size, m_data, rMAX2);
    ChainMesh_data_rMAX1.get_searching_region(rMAX1);
    ChainMesh_data_rMAX2.set_par(cell_size, m_data, rMAX2);
  }
  if (count_rrr==1 || count_drr==1) {
    ChainMesh_random_rMAX1.set_par(cell_size, m_random, rMAX2);
    ChainMesh_random_rMAX1.get_searching_region(rMAX1);
    ChainMesh_random_rMAX2.set_par(cell_size, m_random, rMAX2);
  }

  
  // ----------- count the number of triplets ----------- 

  string file;
  
  cout << endl << par::col_green << "data-data-data" << par::col_default << endl;
  
  file = "ddd.dat";
  
  if (count_ddd==1) { 
    count_triplets(m_data, ChainMesh_data_rMAX1, ChainMesh_data_rMAX2, m_ddd, tcount);
    if (dir_output_triplets!=par::defaultString) write_triplets(m_ddd, dir_output_triplets, file);
  } 
  else if (count_ddd==0) read_triplets (m_ddd, dir_input_triplets, file);


  cout << endl << par::col_green << "random-random-random" << par::col_default << endl;
  
  file = "rrr.dat";
  
  if (count_rrr==1) {
    count_triplets(m_random, ChainMesh_random_rMAX1, ChainMesh_random_rMAX2, m_rrr, tcount);
    if (dir_output_triplets!=par::defaultString) write_triplets(m_rrr, dir_output_triplets, file);
  } 
  else if (count_rrr==0) read_triplets (m_rrr, dir_input_triplets, file);
 

  cout << endl << par::col_green << "data-data-random" << par::col_default << endl;
  
  file = "ddr.dat";
  
  if (count_ddr==1) {

    shared_ptr<Triplet> ddr1 = move(Triplet::Create(m_ddr->tripletType(), m_ddr->side_s(), m_ddr->side_u(), m_ddr->perc_increase(), m_ddr->nbins()));
    shared_ptr<Triplet> ddr2 = move(Triplet::Create(m_ddr->tripletType(), m_ddr->side_s(), m_ddr->side_u(), m_ddr->perc_increase(), m_ddr->nbins()));
    shared_ptr<Triplet> ddr3 = move(Triplet::Create(m_ddr->tripletType(), m_ddr->side_s(), m_ddr->side_u(), m_ddr->perc_increase(), m_ddr->nbins()));

    count_triplets(m_data, ChainMesh_data_rMAX1, ChainMesh_random_rMAX2, ddr1, tcount);
    count_triplets(m_data, ChainMesh_random_rMAX1, ChainMesh_data_rMAX2, ddr2, tcount);
    count_triplets(m_random, ChainMesh_data_rMAX1, ChainMesh_data_rMAX2, ddr3, tcount);

    m_ddr->Sum(ddr1); m_ddr->Sum(ddr2); m_ddr->Sum(ddr3); 
   
    if (dir_output_triplets!=par::defaultString) write_triplets (m_ddr, dir_output_triplets, file);
  } 

  else if (count_ddr==0) read_triplets (m_ddr, dir_input_triplets, file);

  
  cout << endl << par::col_green << "data-random-random" << par::col_default << endl;
  
  file = "drr.dat";
  
  if (count_drr==1) {

    shared_ptr<Triplet> drr1 = move(Triplet::Create(m_drr->tripletType(), m_drr->side_s(), m_drr->side_u(), m_drr->perc_increase(), m_drr->nbins()));
    shared_ptr<Triplet> drr2 = move(Triplet::Create(m_drr->tripletType(), m_drr->side_s(), m_drr->side_u(), m_drr->perc_increase(), m_drr->nbins()));
    shared_ptr<Triplet> drr3 = move(Triplet::Create(m_drr->tripletType(), m_drr->side_s(), m_drr->side_u(), m_drr->perc_increase(), m_drr->nbins()));
    
    count_triplets(m_random, ChainMesh_random_rMAX1, ChainMesh_data_rMAX2, drr1, tcount);
    count_triplets(m_random, ChainMesh_data_rMAX1, ChainMesh_random_rMAX2, drr2, tcount);
    count_triplets(m_data, ChainMesh_random_rMAX1, ChainMesh_random_rMAX2, drr3, tcount);

    m_drr->Sum(drr1); m_drr->Sum(drr2); m_drr->Sum(drr3);
    
    if (dir_output_triplets!=par::defaultString) write_triplets (m_drr, dir_output_triplets, file);
  } 

  else if (count_drr==0) 
    read_triplets (m_drr, dir_input_triplets, file);
}


// ============================================================================


void cosmobl::threept::ThreePointCorrelation::write_triplets (shared_ptr<triplets::Triplet> TT, const string dir, const string file) const
{  
  string MK = "mkdir -p "+dir;
  if (system (MK.c_str())) {};
  
  string file_out = dir+file;
  ofstream fout (file_out.c_str()); checkIO(file_out, 0);

  for (int i=0; i<TT->nbins(); i++) 
    fout << TT->TT1D(i) << endl;
  
  fout.clear(); fout.close(); cout << "I wrote the file " << file << endl << endl;
}



// ============================================================================


void cosmobl::threept::ThreePointCorrelation::read_triplets (shared_ptr<triplets::Triplet> TT, const vector<string> dir, const string file) 
{
 if (dir.size()==0)
    ErrorMsg ("Error in cosmobl::twopt::TwoPointCorrelation1D::read_triplets of TwoPointCorrelation1D.cpp! dir.size()=0!");
      
  for (size_t dd=0; dd<dir.size(); dd++) {
        
    string file_in = dir[dd]+file; 
    cout << "I'm reading the triplet file: " << file_in << endl;
    
    ifstream fin(file_in.c_str()); checkIO(file_in, 1);
   
    double pp;
    for (int i=0; i<TT->nbins(); i++) {
      fin >>pp;
      TT->add_TT1D(i, pp);
    }
    
    fin.clear(); fin.close(); cout << "I read the file " << file_in << endl;
  }
}
 
