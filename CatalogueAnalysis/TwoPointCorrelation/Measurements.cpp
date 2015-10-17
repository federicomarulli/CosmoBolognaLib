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
 *  @file CatalogueAnalysis/TwoPointCorrelation/Measurements.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation used to measure
 *  the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation used to measure the two-point correlation
 *  function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "TwoPointCorrelation.h"
using namespace cosmobl;


// ============================================================================


void cosmobl::TwoPointCorrelation::count_pairs (shared_ptr<Catalogue> cat1, ChainMesh_Catalogue &ChM, Pairs &pp, bool cross, bool tcount)  
{
  time_t start; time (&start);

  int nObj = cat1->nObjects();

  float fact_count = 100./nObj;
  int dp = cout.precision();
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);

  
  void (Pairs::*Put)(shared_ptr<Object>, shared_ptr<Object>) = (m_do_all) ? &Pairs::put : &Pairs::put_log;
  
  shared_ptr<Catalogue> cat2 = ChM.catalogue();


  int tid = 0;
  
#pragma omp parallel num_threads(omp_get_max_threads()) private(tid) 
  {
    tid = omp_get_thread_num();

    if (tid == 0) {
      int nthreads = omp_get_num_threads();
      cout <<"Number of threads = "<<nthreads<<endl;
    }
    
#pragma omp for schedule(static, 2) 
    for (int i=0; i<nObj; i++) { // loop on the objects of the catalogue    
      
      // get the indexes of close objects (i.e. objects inside the close cells of the chain-mesh)
      vector<long> close_objects = ChM.close_objects(cat1->coordinates(i), (cross) ? -1 : (long)i);
      
      for (auto &&j : close_objects) // loop on the closed objects 
        (pp.*Put)(cat1->object(i), cat2->object(j)); // estimate the distance between the two objects and update the pair count
	
      // estimate the computational time and update the time count
      time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
      if (tcount && tid==0) { cout <<"\r..."<<float(i)*fact_count<<"% completed ("<<diff_temp/60<<" minutes)\r"; cout.flush(); }
      if (i==int(nObj*0.25)) cout <<"......25% completed"<<endl;
      if (i==int(nObj*0.5)) cout <<"......50% completed"<<endl;
      if (i==int(nObj*0.75)) cout <<"......75% completed"<<endl;   
    }
    
  }
  
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (tid==0) {
    if (diff<3600) cout <<"   time spent to count the pairs: "<<diff/60<<" minutes"<<endl<<endl;
    else cout <<"   time spent to count the pairs: "<<diff/3600<<" hours"<<endl<<endl;
  }
  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(dp);
}


// ============================================================================


void cosmobl::TwoPointCorrelation::count_pairs_regions (shared_ptr<Catalogue> cat1, ChainMesh_Catalogue &ChM, vector<shared_ptr<Pairs>> pp, bool cross, bool tcount)
{
  time_t start; time (&start);
  int dp = cout.precision();

  int nObj = cat1->nObjects();
  int nRegions = cat1->Nregion();

  float fact_count = 100./nObj;
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);

  shared_ptr<Catalogue> cat2 = ChM.catalogue();

  int tid = 0;
  
#pragma omp parallel num_threads(omp_get_max_threads()) private(tid)
  {
    tid = omp_get_thread_num();

    if (tid == 0) {
      int nthreads = omp_get_num_threads();
      cout <<"Number of threads = "<<nthreads<<endl;
    }

#pragma omp for schedule(static, 2) 
    for (int i=0; i<nObj; i++) {
    
      vector<long> close_objects = ChM.close_objects(cat1->coordinates(i), (cross) ? -1 : (long)i);
    
      for (auto &&j : close_objects) {      
       
	int reg1 = min(cat1->region(i), cat2->region(j));
	int reg2 = max(cat1->region(i), cat2->region(j));
      
	int index = reg1*nRegions+reg2-(reg1-1)*reg1/2-reg1; //reg1*nRegions+reg2;
      
	pp[index]->put(cat1->object(i), cat2->object(j));

      }
    
      // estimate the computational time and update the time count
      time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
      if (tcount && tid==0) { cout <<"\r..."<<float(i)*fact_count<<"% completed ("<<diff_temp/60<<" minutes)\r"; cout.flush(); }
      if (i==int(nObj*0.25)) cout <<"......25% completed"<<endl;
      if (i==int(nObj*0.5)) cout <<"......50% completed"<<endl;
      if (i==int(nObj*0.75)) cout <<"......75% completed"<<endl; 
    }

  }
  
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<3600) cout <<"   time spent to count the pairs: "<<diff/60<<" minutes"<<endl<<endl;
  else cout <<"   time spent to count the pairs: "<<diff/3600<<" hours"<<endl<<endl;
  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(dp);

} 


// ============================================================================

// measure the two-point correlation function and estimate analytic errors 

void cosmobl::TwoPointCorrelation::measure_xi (string dir_output_pairs, vector<string> dir_input_pairs, int count_gg, int count_rr, int count_gr, bool doGR, bool tcount)
{
  if (m_nGal==0 || m_nRan==0)  
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::measure_xi of Measurements.cpp! cat1->nObject()=0!");
 
  if (count_gg>-1 || count_rr>-1 || count_gr>-1) allocate_vectors_xi(doGR);

  // ----------- compute polar coordinates, if necessary ----------- 

  if (!m_read_only) {
    
    if (!isSet(m_data->var(Var::_RA_)) || !isSet(m_data->var(Var::_DEC_)) || !isSet(m_data->var(Var::_DC_))) 
      m_data->computePolarCoordinates();
    
    if (!isSet(m_random->var(Var::_RA_)) || !isSet(m_random->var(Var::_DEC_)) || !isSet(m_random->var(Var::_DC_))) 
      m_random->computePolarCoordinates();

  }
  
  
  // ----------- create the chain-mesh ----------- 
  
  double cell_size = m_rMAX_eff*0.1;

  ChainMesh_Catalogue ChM_data, ChM_random;

  if (!m_read_only) {

    if (count_gg==1)
      ChM_data.set_par(cell_size, m_data, m_rMAX_eff);

    if (count_rr==1 || count_gr==1) 
      ChM_random.set_par(cell_size, m_random, m_rMAX_eff);

  }
  
  
  // ----------- count the number of pairs ----------- 

  string file;

  
  // ===== Object-Object =====

  file = "gg";
  cout <<"Object-Object"<<endl;

  if (count_gg==1 && !m_read_only) {

    Pairs3D gg (m_nlogbins, m_nlinbins, m_ncosbins, m_rMIN, m_rMAX_eff, m_logbinsz, m_linbinsz, m_cosbinsz);

    count_pairs(m_data, ChM_data, gg, 0, tcount);       
    
    if (!isDimEqual(m_gg_log, gg.PPlog()) || !isDimEqual(m_gg_lin, gg.PPlin()) || 
	!isDimEqual(m_gg_2d, gg.PP2d()) || !isDimEqual(m_gg_slog, gg.PPslog()) || 
	!isDimEqual(m_gg_coslog, gg.PPcoslog()) || !isDimEqual(m_gg_coslin, gg.PPcoslin()))
      ErrorMsg("Error in of cosmobl::TwoPointCorrelation::measure_xi of Measurements.cpp!");

    m_gg_log = gg.PPlog(); m_gg_lin = gg.PPlin();
    m_gg_2d = gg.PP2d(); m_gg_slog = gg.PPslog(); m_gg_coslog = gg.PPcoslog(); m_gg_coslin = gg.PPcoslin();

    if (dir_output_pairs!="NULL") write_pairs (m_gg_log, m_gg_lin, m_gg_2d, m_gg_slog, m_gg_coslog, m_gg_coslin, dir_output_pairs, file);
    
    if (!m_read_only) m_data->Order();
  }

  else if (count_gg==0)
    read_pairs (m_gg_log, m_gg_lin, m_gg_2d, m_gg_slog, m_gg_coslog, m_gg_coslin, dir_input_pairs, file);
  
  
  // ===== Random-Random =====
  
  file = "rr";
  cout <<"Random-Random"<<endl; 
  
  if (count_rr==1 && !m_read_only) { 

    Pairs3D rr (m_nlogbins, m_nlinbins, m_ncosbins, m_rMIN, m_rMAX_eff, m_logbinsz, m_linbinsz, m_cosbinsz);
   
    count_pairs(m_random, ChM_random, rr, 0, tcount);
    
    if (!isDimEqual(m_rr_log, rr.PPlog()) || !isDimEqual(m_rr_lin, rr.PPlin()) || 
	!isDimEqual(m_rr_2d, rr.PP2d()) || !isDimEqual(m_rr_slog, rr.PPslog()) || 
	!isDimEqual(m_rr_coslog, rr.PPcoslog()) || !isDimEqual(m_rr_coslin, rr.PPcoslin()))
      ErrorMsg("Error in of cosmobl::TwoPointCorrelation::measure_xi of Measurements.cpp!");
  
    m_rr_log = rr.PPlog(); m_rr_lin = rr.PPlin();
    m_rr_2d = rr.PP2d(); m_rr_slog = rr.PPslog(); m_rr_coslog = rr.PPcoslog(); m_rr_coslin = rr.PPcoslin();

    if (dir_output_pairs!="NULL") write_pairs (m_rr_log, m_rr_lin, m_rr_2d, m_rr_slog, m_rr_coslog, m_rr_coslin, dir_output_pairs, file);

    if ((!doGR || count_gr!=1) && !m_read_only) m_random->Order();
  }

  else if (count_rr==0)
    read_pairs (m_rr_log, m_rr_lin, m_rr_2d, m_rr_slog, m_rr_coslog, m_rr_coslin, dir_input_pairs, file);


  // ===== Object-Random =====

  if (doGR) {

    file = "gr";
    cout <<"Object-Random"<<endl;

    if (count_gr==1 && !m_read_only) {

      Pairs3D gr (m_nlogbins, m_nlinbins, m_ncosbins, m_rMIN, m_rMAX_eff, m_logbinsz, m_linbinsz, m_cosbinsz);

      count_pairs(m_data, ChM_random, gr, 1, tcount); 

      if (!isDimEqual(m_gr_log, gr.PPlog()) || !isDimEqual(m_gr_lin, gr.PPlin()) || 
	!isDimEqual(m_gr_2d, gr.PP2d()) || !isDimEqual(m_gr_slog, gr.PPslog()) || 
	!isDimEqual(m_gr_coslog, gr.PPcoslog()) || !isDimEqual(m_gr_coslin, gr.PPcoslin()))
      ErrorMsg("Error in of cosmobl::TwoPointCorrelation::measure_xi of Measurements.cpp!");
  
      m_gr_log = gr.PPlog(); m_gr_lin = gr.PPlin();
      m_gr_2d = gr.PP2d(); m_gr_slog = gr.PPslog(); m_gr_coslog = gr.PPcoslog(); m_gr_coslin = gr.PPcoslin();
      
      if (dir_output_pairs!="NULL") write_pairs (m_gr_log, m_gr_lin, m_gr_2d, m_gr_slog, m_gr_coslog, m_gr_coslin, dir_output_pairs, file);
   
    }

    else if (count_gr==0) 
      read_pairs (m_gr_log, m_gr_lin, m_gr_2d, m_gr_slog, m_gr_coslog, m_gr_coslin, dir_input_pairs, file);
    
    if (!m_read_only) m_random->Order();
  }

  
  // ----------- compute the correlation functions ----------- 

  m_N_R = double(m_nRan)/m_nGal; 

  double norm = double(m_nRan)*double(m_nRan-1)/(m_nGal*double(m_nGal-1));
  double norm1 = double(m_nRan-1)/m_nGal;

  cout <<"nGal = "<<m_nGal<<endl;
  cout <<"nRan = "<<m_nRan<<endl;
  cout <<"N_R = "<<m_N_R<<endl;
  cout <<"norm = "<<norm<<endl;
  cout <<"norm1 = "<<norm1<<endl;
 

  // \xi(r) in logarithmic bins
  
  for (int i=0; i<m_nlogbins; i++) {
    m_rad_log[i] = pow(10.,i*m_logbinsz+m_shift_log+log10(m_rMIN));
    if (m_gg_log[i]>0 && m_rr_log[i]>0) {
      if (doGR) {
	m_xi_log[i] = max(-1.,norm*m_gg_log[i]/m_rr_log[i]-norm1*m_gr_log[i]/m_rr_log[i]+1.);
	m_error_xi_log[i] = Error(m_gg_log[i], m_rr_log[i], m_gr_log[i]);
      } else {
	m_xi_log[i] = max(-1.,norm*m_gg_log[i]/m_rr_log[i]-1.);
	m_error_xi_log[i] = Error(m_gg_log[i], m_rr_log[i]);
      }
    }
  }
 

  // \xi(r) in linear bins
 
  for (int i=0; i<m_nlinbins; i++) {
    m_rad_lin[i] = i*m_linbinsz+m_shift_lin+m_rMIN;
    if (m_gg_lin[i]>0 && m_rr_lin[i]>0) {
      if (doGR) {
	m_xi_lin[i] = max(-1.,norm*m_gg_lin[i]/m_rr_lin[i]-norm1*m_gr_lin[i]/m_rr_lin[i]+1.);
	m_error_xi_lin[i] = Error(m_gg_lin[i],m_rr_lin[i],m_gr_lin[i]);
      } else {
	m_xi_lin[i] = max(-1.,norm*m_gg_lin[i]/m_rr_lin[i]-1.);
	m_error_xi_lin[i] = Error(m_gg_lin[i],m_rr_lin[i]);
      }
    }
  }


  // \xi(r_p,\pi) in linear bins 

  for (int i=0; i<m_nlinbins; i++) 
    for (int j=0; j<m_nlinbins; j++) {
      if (m_gg_2d[i][j]>0 && m_rr_2d[i][j]>0) {
	if (doGR) {
	  m_xi_2d_lin[i][j] = max(-1.,norm*m_gg_2d[i][j]/m_rr_2d[i][j]-norm1*m_gr_2d[i][j]/m_rr_2d[i][j]+1.);
	  m_error_xi_2d_lin[i][j] = Error(m_gg_2d[i][j],m_rr_2d[i][j],m_gr_2d[i][j]);
	} else {
	  m_xi_2d_lin[i][j] = max(-1.,norm*m_gg_2d[i][j]/m_rr_2d[i][j]-1.);
	  m_error_xi_2d_lin[i][j] = Error(m_gg_2d[i][j],m_rr_2d[i][j]);
	}
      }
    }
  

  // \xi(r_p,\pi) in logarithmic-linear bins

  for (int i=0; i<m_nlogbins; i++) 
    for (int j=0; j<m_nlinbins; j++) 
      if (m_gg_slog[i][j]>0 && m_rr_slog[i][j]>0) {
	if (doGR) {
	  m_xi_2d_loglin[i][j] = max(-1.,norm*m_gg_slog[i][j]/m_rr_slog[i][j]-norm1*m_gr_slog[i][j]/m_rr_slog[i][j]+1.);
	  m_error_xi_2d_loglin[i][j] = Error(m_gg_slog[i][j],m_rr_slog[i][j],m_gr_slog[i][j]);
	} else {
	  m_xi_2d_loglin[i][j] = max(-1.,norm*m_gg_slog[i][j]/m_rr_slog[i][j]-1.);
	  m_error_xi_2d_loglin[i][j] = Error(m_gg_slog[i][j],m_rr_slog[i][j]);
	}
      }  


  // \xi(rr,cos(theta)) with rr in logarithmic bins 

  for (int i=0; i<m_nlogbins; i++) 
    for (int j=0; j<m_ncosbins; j++) {
      m_cos_lin[j] = j*m_cosbinsz+m_shift_cos;
      if (m_gg_coslog[i][j]>0 && m_rr_coslog[i][j]>0) {
	if (doGR) {
	  m_xi_coslog[i][j] = max(-1.,norm*m_gg_coslog[i][j]/m_rr_coslog[i][j]-norm1*m_gr_coslog[i][j]/m_rr_coslog[i][j]+1.);
	  m_error_xi_coslog[i][j] = Error(m_gg_coslog[i][j],m_rr_coslog[i][j],m_gr_coslog[i][j]);
	} else {
	  m_xi_coslog[i][j] = max(-1.,norm*m_gg_coslog[i][j]/m_rr_coslog[i][j]-1.);
	  m_error_xi_coslog[i][j] = Error(m_gg_coslog[i][j],m_rr_coslog[i][j]);
	}
      }
    }
  

  // \xi(rr,cos(theta)) with rr in linear bins 

  for (int i=0; i<m_nlinbins; i++) 
    for (int j=0; j<m_ncosbins; j++) {
      if (m_gg_coslin[i][j]>0 && m_rr_coslin[i][j]>0) {
	if (doGR) {
	  m_xi_coslin[i][j] = max(-1.,norm*m_gg_coslin[i][j]/m_rr_coslin[i][j]-norm1*m_gr_coslin[i][j]/m_rr_coslin[i][j]+1.);
	  m_error_xi_coslin[i][j] = Error(m_gg_coslin[i][j],m_rr_coslin[i][j],m_gr_coslin[i][j]);
	} else {
	  m_xi_coslin[i][j] = max(-1.,norm*m_gg_coslin[i][j]/m_rr_coslin[i][j]-1.);
	  m_error_xi_coslin[i][j] = Error(m_gg_coslin[i][j],m_rr_coslin[i][j]);
	}
      }
    }    

}


// ============================================================================

// measure the two-point correlation function and estimate the internal errors (i.e. jackknife or bootstrap)

void cosmobl::TwoPointCorrelation::measure_xi (string dir_output_pairs, string dir_output_pairs_subSamples, string dir_output_covariance, bool doJK, int nSamples, vector<string> dir_input_pairs, int count_gg, int count_rr, int count_gr, bool doGR, string suffix, bool tcount)
{
  if (m_data->region(0)<0 || m_random->region(0)<0)
      ErrorMsg("Error in cosmobl::TwoPointCorrelation::measure_xi of Measurements.cpp! The function set_ObjectRegion_SubBoxes has not been called!");
  
  if (count_gg>-1 || count_rr>-1 || count_gr>-1) allocate_vectors_xi(doGR);

  string MKDIR = "mkdir -p "+dir_output_pairs;
  if (system(MKDIR.c_str())) {};


  // ----------- compute polar coordinates, if necessary ----------- 

  if (!isSet(m_data->var(Var::_RA_)) || !isSet(m_data->var(Var::_DEC_)) || !isSet(m_data->var(Var::_DC_))) 
    m_data->computePolarCoordinates();

  if (!isSet(m_random->var(Var::_RA_)) || !isSet(m_random->var(Var::_DEC_)) || !isSet(m_random->var(Var::_DC_))) 
    m_random->computePolarCoordinates();

  if (!isSet(m_data->var(Var::_REGION_))) 
    ErrorMsg("Error in print of cosmobl::TwoPointCorrelation::measure_xi of Errors.cpp, no subSampling in data catalogue!");
  if (!isSet(m_random->var(Var::_REGION_))) 
    ErrorMsg("Error in print of cosmobl::TwoPointCorrelation::measure_xi of Errors.cpp, no subSampling in random catalogue!");

  vector<long> region_list = m_data->get_region_list();
  int nRegions = region_list.size();

  
  // ----------- create the chain-mesh ----------- 

  cout << "I'm creating the chain-mesh vectors..." << endl;
  
  double cell_size = m_rMAX_eff*0.1;

  ChainMesh_Catalogue ChM_data, ChM_random;
  
  if (count_gg==1)
    ChM_data.set_par(cell_size, m_data, m_rMAX_eff);
  
  if (count_rr==1 || count_gr==1) 
    ChM_random.set_par(cell_size, m_random, m_rMAX_eff);

  int nP = nRegions*(nRegions+1)/2;
  vector<shared_ptr<Pairs>> gg(nP), rr(nP), gr(nP);
  
  for (int i=0; i<nP; i++) {
    shared_ptr<Pairs3D> ggt {new Pairs3D(m_nlogbins, m_nlinbins, m_ncosbins, m_rMIN, m_rMAX_eff, m_logbinsz, m_linbinsz, m_cosbinsz)};
    shared_ptr<Pairs3D> rrt {new Pairs3D(m_nlogbins, m_nlinbins, m_ncosbins, m_rMIN, m_rMAX_eff, m_logbinsz, m_linbinsz, m_cosbinsz)};
    shared_ptr<Pairs3D> grt {new Pairs3D(m_nlogbins, m_nlinbins, m_ncosbins, m_rMIN, m_rMAX_eff, m_logbinsz, m_linbinsz, m_cosbinsz)};
    gg[i] = ggt;
    rr[i] = rrt;
    gr[i] = grt;
  }
  
  vector<shared_ptr<Catalogue>> data_subSamples(nRegions), random_subSamples(nRegions);
  
  for (int i=0; i<nRegions; i++) {
    double start = region_list[i];
    double stop = start+1;
    data_subSamples[i] = m_data->cut(Var::_REGION_, start, stop);
    random_subSamples[i] = m_random->cut(Var::_REGION_, start, stop);
  }
  
  
  // ===== Object-Object =====
  
  string file = "gg";
  cout <<"Object-Object"<<endl; 

  if (count_gg==1) {

    shared_ptr<Pairs3D> gg_tot = shared_ptr<Pairs3D> {new Pairs3D(m_nlogbins, m_nlinbins, m_ncosbins, m_rMIN, m_rMAX_eff, m_logbinsz, m_linbinsz, m_cosbinsz)};
  
    count_pairs_regions(m_data, ChM_data, gg, 0, tcount);       
    
    for (int i=0; i<nP; i++) gg_tot->sum(gg[i]);

    m_gg_log = gg_tot->PPlog(); m_gg_lin = gg_tot->PPlin();
    m_gg_2d = gg_tot->PP2d(); m_gg_slog = gg_tot->PPslog(); m_gg_coslog = gg_tot->PPcoslog(); m_gg_coslin = gg_tot->PPcoslin();

    if (dir_output_pairs!="NULL") 
      write_pairs(m_gg_log, m_gg_lin, m_gg_2d, m_gg_slog, m_gg_coslog, m_gg_coslin, dir_output_pairs, file);
    
    if (dir_output_pairs_subSamples!="NULL") 
      write_pairs_subSamples(gg, nRegions, dir_output_pairs_subSamples, file);
  }

  else if (count_gg==0) {
    read_pairs(m_gg_log, m_gg_lin, m_gg_2d, m_gg_slog, m_gg_coslog, m_gg_coslin, dir_input_pairs, file);
    read_pairs_subSamples(gg, nRegions, dir_input_pairs, file);
  }

  
  // ===== Random-Random =====

  file = "rr";
  cout <<"Random-Random"<<endl; 

  if (count_rr==1) {

    shared_ptr<Pairs3D> rr_tot = shared_ptr<Pairs3D> {new Pairs3D(m_nlogbins, m_nlinbins, m_ncosbins, m_rMIN, m_rMAX_eff, m_logbinsz, m_linbinsz, m_cosbinsz)};
    count_pairs_regions(m_random, ChM_random, rr, 0, tcount);       

    for (int i=0; i<nP; i++) rr_tot->sum(rr[i]);

    m_rr_log = rr_tot->PPlog(); m_rr_lin = rr_tot->PPlin();
    m_rr_2d = rr_tot->PP2d(); m_rr_slog = rr_tot->PPslog(); m_rr_coslog = rr_tot->PPcoslog(); m_rr_coslin = rr_tot->PPcoslin();
    
    if (dir_output_pairs!="NULL") 
      write_pairs(m_rr_log, m_rr_lin, m_rr_2d, m_rr_slog, m_rr_coslog, m_rr_coslin, dir_output_pairs, file);

    if (dir_output_pairs_subSamples!="NULL") 
      write_pairs_subSamples(rr, nRegions, dir_output_pairs_subSamples, file);
  }

  else if (count_rr==0) {
    read_pairs(m_rr_log, m_rr_lin, m_rr_2d, m_rr_slog, m_rr_coslog, m_rr_coslin, dir_input_pairs, file);
    read_pairs_subSamples(rr, nRegions, dir_input_pairs, file);
  }

  
  // ===== Object-Random =====
  
  if (doGR) {

    file = "gr";
    cout <<"Object-Random"<<endl;

    if (count_gr==1) {

      shared_ptr<Pairs3D> gr_tot = shared_ptr<Pairs3D> {new Pairs3D(m_nlogbins, m_nlinbins, m_ncosbins, m_rMIN, m_rMAX_eff, m_logbinsz, m_linbinsz, m_cosbinsz)};
      count_pairs_regions(m_data, ChM_random, gr, 1, tcount);
      
      for (int i=0; i<nP; i++) gr_tot->sum(gr[i]);

      m_gr_log = gr_tot->PPlog(); m_gr_lin = gr_tot->PPlin();
      m_gr_2d = gr_tot->PP2d(); m_gr_slog = gr_tot->PPslog(); m_gr_coslog = gr_tot->PPcoslog(); m_gr_coslin = gr_tot->PPcoslin();

      if (dir_output_pairs!="NULL") 
	write_pairs(m_gr_log, m_gr_lin, m_gr_2d, m_gr_slog, m_gr_coslog, m_gr_coslin, dir_output_pairs, file);
      
      if (dir_output_pairs_subSamples!="NULL") 
	write_pairs_subSamples(gr, nRegions, dir_output_pairs_subSamples, file);
    }
    
    else if (count_gr==0) {
      read_pairs(m_gr_log, m_gr_lin, m_gr_2d, m_gr_slog, m_gr_coslog, m_gr_coslin, dir_input_pairs, file);
      read_pairs_subSamples(gr, nRegions, dir_input_pairs, file);
    }

  }


  // ================= measure the xi(r) total ================= 
  
  measure_xi(dir_input_pairs, -1, -1, -1, doGR);

  nSamples = (doJK) ? nRegions : nSamples;
  
  m_twop_mock.resize(nSamples);

  Ran w(time(NULL));
  
  for (int i=0; i<nSamples; i++) {

    cout << endl << endl << "sample: " << i << " (of " << nSamples << ")" << endl << endl;
    
    vector<int> ww(nRegions, 0); 
    
    if (doJK) 
      for (int j=0; j<nRegions; j++)
	ww[j] = (i!=j) ? 1 : 0;
    
    else
      for (int j=0; j<nRegions; j++)
	ww[w.int64()%nRegions] ++; 
    
    shared_ptr<Pairs3D> gg_Sample = shared_ptr<Pairs3D> {new Pairs3D(m_nlogbins, m_nlinbins, m_ncosbins, m_rMIN, m_rMAX_eff, m_logbinsz, m_linbinsz, m_cosbinsz)};
    shared_ptr<Pairs3D> rr_Sample = shared_ptr<Pairs3D> {new Pairs3D(m_nlogbins, m_nlinbins, m_ncosbins, m_rMIN, m_rMAX_eff, m_logbinsz, m_linbinsz, m_cosbinsz)};
    shared_ptr<Pairs3D> gr_Sample = shared_ptr<Pairs3D> {new Pairs3D(m_nlogbins, m_nlinbins, m_ncosbins, m_rMIN, m_rMAX_eff, m_logbinsz, m_linbinsz, m_cosbinsz)};
    
    vector<shared_ptr<Object> > data_SS, random_SS;
    for (int cc=0; cc<nRegions; cc++) 
      for (int tt=0; tt<ww[cc]; tt++) {
	
	for (int kk=0; kk<data_subSamples[cc]->nObjects(); kk++)
	  data_SS.push_back(data_subSamples[cc]->object(kk));      
	
	for(int kk=0; kk<random_subSamples[cc]->nObjects(); kk++)
	  random_SS.push_back(random_subSamples[cc]->object(kk));	
      }
    
    shared_ptr<Catalogue> p_data(new Catalogue{data_SS}), p_random(new Catalogue{random_SS});

    
    for (int c1=0; c1<nRegions; c1++) 
      for (int c2=0; c2<nRegions; c2++) {
	
	int j = c1*nRegions+c2-(c1-1)*c1/2-c1; //int j = c1*nRegions+c2;
	double www = (c1==c2) ? double(ww[c1]) : double(ww[c1]*ww[c2]);
	
	gg_Sample->sum(gg[j], www);
	rr_Sample->sum(rr[j], www);
	gr_Sample->sum(gr[j], www);
	    
      }
    
    shared_ptr<TwoPointCorrelation> twop(new TwoPointCorrelation{p_data, p_random});
   
    twop->setParameters(m_rMIN, m_rMAX, m_logbinsz, m_linbinsz, m_cosbinsz);
    twop->allocate_vectors_xi(doGR);

    twop->set_gg_pairs(gg_Sample->PPlog(), gg_Sample->PPlin(), gg_Sample->PP2d(), gg_Sample->PPslog(), gg_Sample->PPcoslog(), gg_Sample->PPcoslin());
    twop->set_rr_pairs(rr_Sample->PPlog(), rr_Sample->PPlin(), rr_Sample->PP2d(), rr_Sample->PPslog(), rr_Sample->PPcoslog(), rr_Sample->PPcoslin());
    twop->set_gr_pairs(gr_Sample->PPlog(), gr_Sample->PPlin(), gr_Sample->PP2d(), gr_Sample->PPslog(), gr_Sample->PPcoslog(), gr_Sample->PPcoslin());
    twop->measure_xi(dir_input_pairs, -1, -1, -1, doGR);

    m_twop_mock[i] = twop;
  }


  get_covariance(dir_output_covariance, doJK, suffix);
   
}


// ============================================================================


void cosmobl::TwoPointCorrelation::measure_wtheta (string dir_output_pairs, vector<string> dir_input_pairs, int count_gg, int count_rr, int count_gr, bool doGR, bool tcount)
{
  allocate_vectors_ACF(doGR);

  shared_ptr<Catalogue> catdata = move(m_data);
  shared_ptr<Catalogue> catrandom = move(m_random);

  catdata->normalizeComovingCoordinates();
  catrandom->normalizeComovingCoordinates();

  
  // ----------- create the chain-mesh ----------- 
  
  double cell_size = m_thetaMAX_eff*0.1;

  ChainMesh_Catalogue ChM_data, ChM_random;
    
  if (count_gg==1)
    ChM_data.set_par(cell_size, m_data, m_thetaMAX_eff);
 
  if (count_rr==1 || count_gr==1)
    ChM_random.set_par(cell_size, m_random, m_thetaMAX_eff);
 

  // ----------- compute the correlation functions ----------- 
  
  string file;

  // --- Object-Object ---

  file = "gg_theta";
  cout <<"Object-Object"<<endl;
  
  if (count_gg==1) {

    Pairs2D gg (m_nlogthetabins, m_nlinthetabins, m_thetaMIN, m_thetaMAX_eff, m_logthetabinsz, m_linthetabinsz);

    count_pairs(catdata, ChM_data, gg, 0, tcount);      
 
    m_gg_theta_log = gg.PPlog(); m_gg_theta_lin = gg.PPlin();

    write_pairs (m_gg_theta_log, m_gg_theta_lin, dir_output_pairs, file);

  }

  else if (count_gg==0) 
    read_pairs (m_gg_theta_log, m_gg_theta_lin, dir_input_pairs, file);

   
  // --- Random-Random ---

  file = "rr_theta";
  cout <<"Random-Random"<<endl; 

  if (count_rr==1) { 

    Pairs2D rr (m_nlogthetabins, m_nlinthetabins, m_thetaMIN, m_thetaMAX_eff, m_logthetabinsz, m_linthetabinsz);

    count_pairs(catrandom, ChM_random, rr, 0, tcount);  	

    m_rr_theta_log = rr.PPlog(); m_rr_theta_lin = rr.PPlin();

    write_pairs (m_rr_theta_log, m_rr_theta_lin, dir_output_pairs, file);

  }

  else if (count_rr==0)
    read_pairs (m_rr_theta_log, m_rr_theta_lin, dir_input_pairs, file);
  

  // --- Object-Random ---

  if (doGR) {

    file = "gr_theta";
    cout <<"Object-Random"<<endl;

    if (count_gr==1) {

      Pairs2D gr (m_nlogthetabins, m_nlinthetabins, m_thetaMIN, m_thetaMAX_eff, m_logthetabinsz, m_linthetabinsz);

      count_pairs(catdata, ChM_random, gr, 1, tcount); 

      m_gr_theta_log = gr.PPlog(); m_gr_theta_lin = gr.PPlin();

      write_pairs (m_gr_theta_log, m_gr_theta_lin, dir_output_pairs, file);

    }

    else if (count_gr==0)
      read_pairs (m_gr_theta_log, m_gr_theta_lin, dir_input_pairs, file);

  }

  catdata->restoreComovingCoordinates();
  catrandom->restoreComovingCoordinates();

  m_N_R = double(m_nRan)/m_nGal; 

  double norm = double(m_nRan)*double(m_nRan-1)/(m_nGal*double(m_nGal-1));
  double norm1 = double(m_nRan-1)/m_nGal;

  cout <<"m_nGal = "<<m_nGal<<", m_data.size() = "<<m_data->nObjects()<<endl;
  cout <<"m_nRan = "<<m_nRan<<", N_R = "<<m_N_R<<endl;
  cout <<"norm = "<<norm<<endl;
  cout <<"norm1 = "<<norm1<<endl;
  

  // \xi(r) in logarithmic bins 

  for (int i=0; i<m_nlogthetabins; i++) {
    m_theta_log[i] = pow(10.,i*m_logthetabinsz+m_shift_theta_log+log10(m_thetaMIN));
    if (m_gg_theta_log[i]>0 && m_rr_theta_log[i]>0) {
      if (doGR) {
	m_wtheta_log[i] = max(-1.,norm*m_gg_theta_log[i]/m_rr_theta_log[i]-norm1*m_gr_theta_log[i]/m_rr_theta_log[i]+1.);
	m_error_wtheta_log[i] = 1./*Error(m_gg_theta_log[i],m_rr_theta_log[i],m_gr_theta_log[i])*/;
	WarningMsg("Attention: the error is not computed and set equal to 1");
      } else {
	m_wtheta_log[i] = max(-1.,norm*m_gg_theta_log[i]/m_rr_theta_log[i]-1.);
	m_error_wtheta_log[i] = 1./*Error(m_gg_theta_log[i],m_rr_theta_log[i])*/;
	WarningMsg("Attention: the error is not computed and set equal to 1");
      }
    }
  }


  // \xi(r) in linear bins
  for (int i=0; i<m_nlinthetabins; i++) {
    m_theta_lin[i] = i*m_linthetabinsz+m_shift_theta_lin+m_thetaMIN;
    if (m_gg_theta_lin[i]>0 && m_rr_theta_lin[i]>0) {
      if (doGR) {
	m_wtheta_lin[i] = max(-1.,norm*m_gg_theta_lin[i]/m_rr_theta_lin[i]-norm1*m_gr_theta_lin[i]/m_rr_theta_lin[i]+1.);
	m_error_wtheta_lin[i] = 1./*Error(m_gg_theta_lin[i],m_rr_theta_lin[i],m_gr_theta_lin[i])*/;
	WarningMsg("Attention: the error is not computed and set equal to 1");
      } else {
	m_wtheta_lin[i] = max(-1.,norm*m_gg_theta_lin[i]/m_rr_theta_lin[i]-1.);
	m_error_wtheta_lin[i] = 1./*Error(m_gg_theta_lin[i],m_rr_theta_lin[i])*/;
	WarningMsg("Attention: the error is not computed and set equal to 1");
      }
    }
  }
}


// ============================================================================


void cosmobl::TwoPointCorrelation::measure_projected_xi (double &pimax, string &dir)
{
  if (m_xi_2d_loglin.size()==0) 
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::measure_projected_xi of Measurements.cpp: xi(rp,pi) has not been measured yet!"); 
  if (pimax==0) 
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::measure_projected_xi of Measurements.cpp: pimax=0!");


  cout <<"I'm computing the projected xi..."<<endl;

  int pim = nint(pimax/m_linbinsz); // to convert from Mpc/h into the vector index
 
  if (pim>int(m_xi_2d_loglin[0].size())) 
    WarningMsg("Attention, I have to extrapolate, as pim = "+conv(pimax,par::fDP1)+" Mpc/h!");


  m_rp_proj.erase (m_rp_proj.begin(), m_rp_proj.end());
  m_xi_proj.erase (m_xi_proj.begin(), m_xi_proj.end());
  m_error_xi_proj.erase (m_error_xi_proj.begin(), m_error_xi_proj.end());
  for (int i=0; i<m_nlogbins; i++) {m_xi_proj.push_back(0.); m_error_xi_proj.push_back(0.);}
  
  for (int i=0; i<m_nlogbins; i++) {
    m_rp_proj.push_back(pow(10.,(i*m_logbinsz+m_shift_log+log10(m_rMIN))));
    
    vector<double> rrf, xif, errf;
    
    for (int j=0; j<pim; j++) {

      if (j<int(m_xi_2d_loglin[i].size())-2) {

	if (m_xi_2d_loglin[i][j]>-1.e20) { // check!!!
	  m_xi_proj[i] += 2.*m_linbinsz*m_xi_2d_loglin[i][j];
	  m_error_xi_proj[i] += pow(2.*m_linbinsz*m_error_xi_2d_loglin[i][j],2); // check!!!!

	  rrf.push_back(m_rad_lin[j]);
	  xif.push_back(m_xi_2d_loglin[i][j]);
	  errf.push_back(m_error_xi_2d_loglin[i][j]);
	}
      }
      
      else {

	ErrorMsg("work in progress... (see measure_projected_xi of Measurements.cpp)");
	
	/*
	double AA = -1.e30, BB = -1.e30, CC = -1.e30;
	if (AA<-1.e29) {
	  AA = 1., BB = 1., CC = 1.;
	  quad_fit (rrf, xif, errf, &AA, &BB, &CC);
	}
	cout <<"--> "<<AA<<"   "<<BB<<"   "<<CC<<endl;exit(1);

	double R1 = 0.001;
	double R2 = 100.;
	int step = 100;
	double DeltaR = (R2-R1)/step;
	double RR = R1;
	for (int kk=0; kk<step; kk++) {
	  cout <<RR<<"   "<<Pol2(RR,AA,BB,CC)<<endl;
	  RR += DeltaR;
	}

	//double RR = rad_lin[j];
	//double XI = Pol2(RR,AA,BB,CC);

	//xi_proj[i] += 2.*linbinsz*XI;
	//error_xi_proj[i] += 1.e20; // check!!!
	*/
      }
    }
  }
  
  for (unsigned int i=0; i<m_error_xi_proj.size(); i++) m_error_xi_proj[i] = sqrt(m_error_xi_proj[i]);
  
  if (dir!="NULL") {
    string file_proj = dir+"xi_projected";
    ofstream fout (file_proj.c_str()); checkIO(file_proj,0);
    
    for (unsigned int i=0; i<m_rp_proj.size(); i++) 
      fout <<m_rp_proj[i]<<"   "<<m_xi_proj[i]<<"   "<<m_error_xi_proj[i]<<endl;
    
    fout.clear(); fout.close();
    cout <<"I wrote the file: "<<file_proj<<endl;
  }

}


// ============================================================================


void cosmobl::TwoPointCorrelation::measure_projected_xi_Mocks (double pimax, string dir_samples)
{
  for (size_t i=0; i<m_twop_mock.size(); i++) {
    string dir = dir_samples+"sample_"+conv(i+1,par::fINT)+"/";
    string cmd ="mkdir -p "+dir; if (system(cmd.c_str())) {}
    m_twop_mock[i]->measure_projected_xi(pimax,dir);
  }
}
