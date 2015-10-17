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
 *  @file CatalogueAnalysis/ThreePointCorrelation/Measurements.cpp
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

#include "ThreePointCorrelation.h"
using namespace cosmobl;


// ============================================================================


void cosmobl::ThreePointCorrelation::measure_Q (string dir_output_triplets, string dir_output_2pt, vector<string> dir_input_triplets, int count_ggg, int count_rrr, int count_ggr, int count_grr, bool tcount) 
{   
  allocate_vectors_zeta();

  cout <<endl<<"***************************************************"<<endl;
  cout <<"I am now computing the 3-point correlation function"<<endl;
  cout <<"***************************************************"<<endl<<endl;

  
  // ----- double chain-mesh -----
  
  double rMAX1 = m_side_s*(1.+2.*m_perc_increase);
  double rMAX2 = m_side_u*m_side_s*(1.+2.*m_perc_increase);
 
  double cell_size = rMAX2*0.1;
  
  ChainMesh_Catalogue lkList_data_rMAX1, lkList_data_rMAX2, lkList_random_rMAX1, lkList_random_rMAX2;
  
  if (count_ggg==1 || count_ggr==1) {
    lkList_data_rMAX1.set_par(cell_size, m_data, rMAX2);
    lkList_data_rMAX1.get_searching_region(rMAX1);
    lkList_data_rMAX2.set_par(cell_size, m_data, rMAX2);
  }
  if (count_rrr==1 || count_grr==1) {
    lkList_random_rMAX1.set_par(cell_size, m_random, rMAX2);
    lkList_random_rMAX1.get_searching_region(rMAX1);
    lkList_random_rMAX2.set_par(cell_size, m_random, rMAX2);
  }

  // ----------- count the number of triplets ----------- 

  string file;

  // --- Object-Object-Object ---
  
  file = "ggg";
  cout <<"Object-Object-Object"<<endl;

  if (count_ggg==1) { 

    Triplets3D ggg {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    
    count_triplets(m_data, lkList_data_rMAX1, lkList_data_rMAX2, ggg, 1, tcount);

    m_ggg = ggg.TT();

    write_triplets(m_ggg, dir_output_triplets, file);
   
  } 

  else if (count_ggg==0) 
    read_triplets (m_ggg, dir_input_triplets, file);


  // --- Random-Random-Random ---

  file = "rrr";
  cout <<"Random-Random-Random"<<endl;

  if (count_rrr==1) {

    Triplets3D rrr {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};

    count_triplets(m_random, lkList_random_rMAX1, lkList_random_rMAX2, rrr, 1, tcount);

    m_rrr = rrr.TT();

    write_triplets(m_rrr, dir_output_triplets, file);

  } 

  else if (count_rrr==0) 
    read_triplets (m_rrr, dir_input_triplets, file);
 

  // --- Object-Object-Random ---

  file = "ggr";
  cout <<"Object-Object-Random"<<endl;

  if (count_ggr==1) {

    Triplets3D ggr1 {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    Triplets3D ggr2 {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    Triplets3D ggr3 {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    
    count_triplets(m_data, lkList_data_rMAX1, lkList_random_rMAX2, ggr1, 1, tcount);
    count_triplets(m_data, lkList_random_rMAX1, lkList_data_rMAX2, ggr2, 1, tcount);
    count_triplets(m_random, lkList_data_rMAX1, lkList_data_rMAX2, ggr3, 1, tcount);

    for (size_t i=0; i<m_ggr.size(); i++) m_ggr[i] = ggr1.TT()[i]+ggr2.TT()[i]+ggr3.TT()[i];

    write_triplets (m_ggr, dir_output_triplets, file);
    
  } 

  else if (count_ggr==0) 
    read_triplets (m_ggr, dir_input_triplets, file);
 

  // --- Random-Random-Random ---

  file = "grr";
  cout <<"Object-Random-Random"<<endl;

  if (count_grr==1) {

    Triplets3D grr1 {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    Triplets3D grr2 {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    Triplets3D grr3 {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    
    count_triplets(m_random, lkList_random_rMAX1, lkList_data_rMAX2, grr1, 1, tcount);
    count_triplets(m_random, lkList_data_rMAX1, lkList_random_rMAX2, grr2, 1, tcount);
    count_triplets(m_data, lkList_random_rMAX1, lkList_random_rMAX2, grr3, 1, tcount);

    for (size_t i=0; i<m_grr.size(); i++) m_grr[i] = grr1.TT()[i]+grr2.TT()[i]+grr3.TT()[i];

    write_triplets (m_grr, dir_output_triplets, file);

  } 

  else if (count_grr==0) 
    read_triplets (m_grr, dir_input_triplets, file);
 


  // ----------- compute the 3-point correlation function ----------- 

  double norm1 = (double(m_nGal)*double(m_nGal-1)*double(m_nGal-2))/6.;
  double norm2 = (double(m_nGal)*double(m_nGal-1)*double(m_nRan))*0.5;
  double norm3 = (double(m_nGal)*double(m_nRan)*double(m_nRan-1))*0.5;
  double norm4 = (double(m_nRan)*double(m_nRan-1)*double(m_nRan-2))/6.;

  cout <<endl<<"---------------------------------------------------"<<endl<<endl;
  cout <<"nGal = "<<m_nGal<<", data.size() = "<<m_data->nObjects()<<endl;
  cout <<"nRan = "<<m_nRan<<endl;
  cout <<"norm(GGG) = "<<norm1<<endl;
  cout <<"norm(GGR) = "<<norm2<<endl;
  cout <<"norm(GRR) = "<<norm3<<endl;
  cout <<"norm(RRR) = "<<norm4<<endl;


  // \zeta(r), connected three-point correlation function
  
  for (int i=0; i<m_nbins; i++) 
    if (m_ggg[i]>0 && m_rrr[i]>0) 
      m_zeta[i] = ((m_ggg[i]/norm1)/(m_rrr[i]/norm4))-3.*((m_ggr[i]/norm2)/(m_rrr[i]/norm4))+3.*(((m_grr[i]/norm3)/(m_rrr[i]/norm4)))-1.;
      

  // -------------------------------------------------

  cout <<endl<<"***************************************************"<<endl;
  cout <<"I am now computing the 2-point correlation function"<<endl;
  cout <<"***************************************************"<<endl<<endl;
 
  double _rMIN = m_side_s;
  double _rMAX = (1+m_side_u)*m_side_s*(1+2*m_perc_increase);
  double _logbinSize = 0.05;
  double _binSize = 1.;
  double _cosSize = 0.02;
  bool doGR = 1;
 
  setParameters(_rMIN, _rMAX, _logbinSize, _binSize, _cosSize);
 
  allocate_vectors_xi(doGR);

  measure_xi(dir_output_triplets, tcount);

  vector<double> log_r(m_nlogbins), log_xi_log(m_nlogbins);
  for (int i=0; i<m_nlogbins; i++) {
    log_r[i] = log10(m_rad_log[i]);
    log_xi_log[i] = log10(1.+m_xi_log[i]);
  }
    
  vector<double> values_interp(m_nbins+2, 0.), xi_real_lin(m_nbins+2, 0.);
  double err = -1, log_xi_real_lin;

  values_interp[0] = log10(m_side_s);
  values_interp[1] = log10(m_side_u*m_side_s);

  for (int i=0; i<m_nbins; i++) {
    double tmp_value = -1.;
    if (m_type_binning=="ang") {
      double theta=((i+0.5)*m_binsize);
      tmp_value = m_side_s*sqrt(1+m_side_u*m_side_u-2*m_side_u*cos(theta));
    }
    else if (m_type_binning=="lin")
      tmp_value = (m_side_s+((i+0.5)*m_binsize));
    values_interp[i+2] = log10(tmp_value);
  }

  string file_2pt = dir_output_2pt+"2ptCorrelation_3pt.dat";
  ofstream fout (file_2pt.c_str()); checkIO(file_2pt, 0);

  for (size_t i=0; i<values_interp.size(); i++) {
    interpolation_extrapolation(values_interp[i], log_r, log_xi_log, "Linear", -1, &log_xi_real_lin, &err);
    xi_real_lin[i] = pow(10.,log_xi_real_lin)-1.;
    //cout << pow(10.,values_interp[i]) << " --- " << xi_real_lin[i] << endl;
    fout << pow(10.,values_interp[i]) << "     " << xi_real_lin[i] << endl;
  }

  fout.clear(); fout.close(); cout <<"I wrote the file "<<file_2pt<<endl;

  
  // ----------- compute the reduced 3-point correlation function -----------

  cout <<endl<<"***********************************************************"<<endl;
  cout <<"I am now computing the reduced 3-point correlation function"<<endl;
  cout <<"***********************************************************"<<endl<<endl;

  for (int i=0; i<m_nbins; i++) 
    m_zeta_red[i] = m_zeta[i]/((xi_real_lin[0]*xi_real_lin[1])+(xi_real_lin[0]*xi_real_lin[i+2])+(xi_real_lin[1]*xi_real_lin[i+2]));
  
}


// ============================================================================


void cosmobl::ThreePointCorrelation::measure_Q_TEST (string dir_output_triplets, string dir_output_2pt, vector<string> dir_input_triplets, bool count, bool tcount) 
{

  allocate_vectors_zeta();

  cout <<endl<<"***************************************************"<<endl;
  cout <<"I am now computing the 3-point correlation function"<<endl;
  cout <<"***************************************************"<<endl<<endl;

  
  // ----- double chain-mesh -----
  
  double rMAX1 = m_side_s*(1.+2.*m_perc_increase);
  double rMAX2 = m_side_u*m_side_s*(1.+2.*m_perc_increase);
 
  double cell_size = rMAX2*0.1;
  
  ChainMesh_Catalogue lkList_data_rMAX1, lkList_data_rMAX2, lkList_random_rMAX1, lkList_random_rMAX2;

  if (count) {
    lkList_data_rMAX1.set_par(cell_size, m_data, rMAX2);
    lkList_data_rMAX1.get_searching_region(rMAX1);
    lkList_data_rMAX2.set_par(cell_size, m_data, rMAX2);
    lkList_random_rMAX1.set_par(cell_size, m_random, rMAX2);
    lkList_random_rMAX1.get_searching_region(rMAX1);
    lkList_random_rMAX2.set_par(cell_size, m_random, rMAX2);
  }

  
  // ----------- count the number of triplets ----------- 


  // --- Object-Object-Object and Object-Object-Random ---
  
  cout <<"Object-Object-Object and Object-Object-Random"<<endl;

  if (count) { 

    Triplets3D ggg {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    Triplets3D ggr {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};

    count_triplets_gg(ggg, ggr, lkList_data_rMAX1, lkList_data_rMAX2, lkList_random_rMAX2, 1, tcount);

    m_ggg = ggg.TT();
    
 
    // --- Random-Random-Random and Random-Random-Object ---

    cout <<"Random-Random-Random and Random-Random-Object"<<endl;

 
    Triplets3D rrg {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    Triplets3D rrr {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};

    count_triplets_rr(rrg, rrr, lkList_random_rMAX1, lkList_data_rMAX2, lkList_random_rMAX2, 1, tcount);
    
    m_rrr = rrr.TT();
    
  
    // --- Object-Random-Object and Object-Random-Random ---

    cout <<"Object-Object-Random and Object-Random-Random"<<endl;

 
    Triplets3D grg {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    Triplets3D grr {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};

    count_triplets_gr(grg, grr, lkList_random_rMAX1, lkList_data_rMAX2, lkList_random_rMAX2, 1, tcount);


    // --- Random-Object-Object and Random-Object-Random ---

    cout <<"Random-Object-Object and Random-Object-Random"<<endl;

    Triplets3D rgg {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    Triplets3D rgr {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};

    count_triplets_rg(rgg, rgr, lkList_data_rMAX1, lkList_data_rMAX2, lkList_random_rMAX2, 1, tcount);

    for (size_t i=0; i<m_ggr.size(); i++) m_ggr[i] = ggr.TT()[i] + grg.TT()[i] + rgg.TT()[i];
    for (size_t i=0; i<m_grr.size(); i++) m_grr[i] = grr.TT()[i] + rrg.TT()[i] + rgr.TT()[i];

    
    write_triplets(m_ggg, dir_output_triplets, "ggg");
    write_triplets(m_ggr, dir_output_triplets, "ggr");
    write_triplets(m_rrr, dir_output_triplets, "rrr");
    write_triplets(m_grr, dir_output_triplets, "grr");
  }

  else {
    read_triplets(m_ggg, dir_input_triplets, "ggg");
    read_triplets(m_ggr, dir_input_triplets, "ggr");
    read_triplets(m_rrr, dir_input_triplets, "rrr");
    read_triplets(m_grr, dir_input_triplets, "grr");
  }

  
  // ----------- compute the 3-point correlation function ----------- 

  double norm1 = (double(m_nGal)*double(m_nGal-1)*double(m_nGal-2))/6.;
  double norm2 = (double(m_nGal)*double(m_nGal-1)*double(m_nRan))*0.5;
  double norm3 = (double(m_nGal)*double(m_nRan)*double(m_nRan-1))*0.5;
  double norm4 = (double(m_nRan)*double(m_nRan-1)*double(m_nRan-2))/6.;

  cout <<endl<<"---------------------------------------------------"<<endl<<endl;
  cout <<"nGal = "<<m_nGal<<", data.size() = "<<m_data->nObjects()<<endl;
  cout <<"nRan = "<<m_nRan<<endl;
  cout <<"norm(GGG) = "<<norm1<<endl;
  cout <<"norm(GGR) = "<<norm2<<endl;
  cout <<"norm(GRR) = "<<norm3<<endl;
  cout <<"norm(RRR) = "<<norm4<<endl;


  // \zeta(r) in linear bins
  
  for (int i=0; i<m_nbins; i++) 
    if (m_ggg[i]>0 && m_rrr[i]>0) 
      m_zeta[i] = ((m_ggg[i]/norm1)/(m_rrr[i]/norm4))-3.*((m_ggr[i]/norm2)/(m_rrr[i]/norm4))+3.*(((m_grr[i]/norm3)/(m_rrr[i]/norm4)))-1.;
      

  // -------------------------------------------------

  cout <<endl<<"***************************************************"<<endl;
  cout <<"I am now computing the 2-point correlation function"<<endl;
  cout <<"***************************************************"<<endl<<endl;
 
  double _rMIN = m_side_s;
  double _rMAX = (1+m_side_u)*m_side_s*(1+2*m_perc_increase);
  double _logbinSize = 0.05;
  double _binSize = 1.;
  double _cosSize = 0.02;
  bool doGR = 1;
 
  setParameters(_rMIN, _rMAX, _logbinSize, _binSize, _cosSize);
 
  allocate_vectors_xi(doGR);

  measure_xi(dir_output_triplets, tcount);

  vector<double> log_r(m_nlogbins), log_xi_log(m_nlogbins);
  for (int i=0; i<m_nlogbins; i++) {
    log_r[i] = log10(m_rad_log[i]);
    log_xi_log[i] = log10(1.+m_xi_log[i]);
  }
    
  vector<double> values_interp(m_nbins+2, 0.), xi_real_lin(m_nbins+2, 0.);
  double err = -1, log_xi_real_lin;

  values_interp[0] = log10(m_side_s);
  values_interp[1] = log10(m_side_u*m_side_s);

  for (int i=0; i<m_nbins; i++) {
    double tmp_value = -1.;
    if (m_type_binning=="ang") {
      double theta=((i+0.5)*m_binsize);
      tmp_value = m_side_s*sqrt(1+m_side_u*m_side_u-2*m_side_u*cos(theta));
    }
    else if (m_type_binning=="lin")
      tmp_value = (m_side_s+((i+0.5)*m_binsize));
    values_interp[i+2] = log10(tmp_value);
  }

  string file_2pt = dir_output_2pt+"2ptCorrelation_3pt.dat";
  ofstream fout (file_2pt.c_str()); checkIO(file_2pt, 0);

  for (size_t i=0; i<values_interp.size(); i++) {
    interpolation_extrapolation(values_interp[i], log_r, log_xi_log, "Linear", -1, &log_xi_real_lin, &err);
    xi_real_lin[i] = pow(10.,log_xi_real_lin)-1.;
    //cout << pow(10.,values_interp[i]) << " --- " << xi_real_lin[i] << endl;
    fout << pow(10.,values_interp[i]) << "     " << xi_real_lin[i] << endl;
  }

  fout.clear(); fout.close(); cout <<"I wrote the file "<<file_2pt<<endl;

  
  // ----------- compute the reduced 3-point correlation function -----------

  cout <<endl<<"***********************************************************"<<endl;
  cout <<"I am now computing the reduced 3-point correlation function"<<endl;
  cout <<"***********************************************************"<<endl<<endl;

  for (int i=0; i<m_nbins; i++) 
    m_zeta_red[i] = m_zeta[i]/((xi_real_lin[0]*xi_real_lin[1])+(xi_real_lin[0]*xi_real_lin[i+2])+(xi_real_lin[1]*xi_real_lin[i+2]));
  
}


// ============================================================================


void cosmobl::ThreePointCorrelation::measure_Q_ang (string dir_output_triplets, string dir_output_2pt, vector<string> dir_input_triplets, int count_ggg, int count_rrr, int count_ggr, int count_grr, bool tcount) 
{

  allocate_vectors_zeta();

  cout <<endl<<"***************************************************"<<endl;
  cout <<"I am now computing the 3-point correlation function"<<endl;
  cout <<"***************************************************"<<endl<<endl;

  
  // ----- double chain-mesh -----
  
  double rMAX1 = m_side_s*(1.+2.*m_perc_increase);
  double rMAX2 = m_side_u*m_side_s*(1.+2.*m_perc_increase);

  double cell_size1 = rMAX1*0.1;
  double cell_size2 = rMAX2*0.1;
  
  ChainMesh_Catalogue lkList_data_rMAX1, lkList_data_rMAX2, lkList_random_rMAX1, lkList_random_rMAX2;

  if (count_ggg==1 || count_ggr==1) {
    lkList_data_rMAX1.set_par(cell_size1, m_data, rMAX1);
    lkList_data_rMAX2.set_par(cell_size2, m_data, rMAX2);
  }
  if (count_rrr==1 || count_grr==1) {
    lkList_random_rMAX1.set_par(cell_size1, m_random, rMAX1);
    lkList_random_rMAX2.set_par(cell_size2, m_random, rMAX2);
  }

  
  // ----------- count the number of triplets ----------- 

  string file;

  // --- Object-Object-Object ---
  
  file = "ggg";
  cout <<"Object-Object-Object"<<endl;

  if (count_ggg==1) { 

    Triplets3D ggg {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    
    count_triplets(m_data, lkList_data_rMAX1, lkList_data_rMAX2, ggg, 0, tcount);

    m_ggg = ggg.TT();

    write_triplets(m_ggg, dir_output_triplets, file);

  } 

  else if (count_ggg==0) 
    read_triplets (m_ggg, dir_input_triplets, file);


  // --- Random-Random-Random ---

  file = "rrr";
  cout <<"Random-Random-Random"<<endl;

  if (count_rrr==1) {

    Triplets3D rrr {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};

    count_triplets(m_random, lkList_random_rMAX1, lkList_random_rMAX2, rrr, 0, tcount);

    m_rrr = rrr.TT();

    write_triplets(m_rrr, dir_output_triplets, file);

  } 

  else if (count_rrr==0) 
    read_triplets (m_rrr, dir_input_triplets, file);
 

  // --- Object-Object-Random ---

  file = "ggr";
  cout <<"Object-Object-Random"<<endl;

  if (count_ggr==1) {

    Triplets3D ggr1 {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    Triplets3D ggr2 {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    Triplets3D ggr3 {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    
    count_triplets(m_data, lkList_data_rMAX1, lkList_random_rMAX2, ggr1, 0, tcount);
    count_triplets(m_data, lkList_random_rMAX1, lkList_data_rMAX2, ggr2, 0, tcount);
    count_triplets(m_random, lkList_data_rMAX1, lkList_data_rMAX2, ggr3, 0, tcount);

    for (size_t i=0; i<m_ggr.size(); i++) m_ggr[i] = ggr1.TT()[i]+ggr2.TT()[i]+ggr3.TT()[i];

    write_triplets (m_ggr, dir_output_triplets, file);
    
  } 

  else if (count_ggr==0) 
    read_triplets (m_ggr, dir_input_triplets, file);
 

  // --- Random-Random-Random ---

  file = "grr";
  cout <<"Object-Random-Random"<<endl;

  if (count_grr==1) {

    Triplets3D grr1 {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    Triplets3D grr2 {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    Triplets3D grr3 {m_type_binning, m_binsize, m_side_s, m_side_u, m_perc_increase};
    
    count_triplets(m_random, lkList_random_rMAX1, lkList_data_rMAX2, grr1, 0, tcount);
    count_triplets(m_random, lkList_data_rMAX1, lkList_random_rMAX2, grr2, 0, tcount);
    count_triplets(m_data, lkList_random_rMAX1, lkList_random_rMAX2, grr3, 0, tcount);

    for (size_t i=0; i<m_grr.size(); i++) m_grr[i] = grr1.TT()[i]+grr2.TT()[i]+grr3.TT()[i];

    write_triplets (m_grr, dir_output_triplets, file);

  } 

  else if (count_grr==0) 
    read_triplets (m_grr, dir_input_triplets, file);
 


  // ----------- compute the 3-point correlation function ----------- 

  double norm1 = (double(m_nGal)*double(m_nGal-1)*double(m_nGal-2))/6.;
  double norm2 = (double(m_nGal)*double(m_nGal-1)*double(m_nRan))*0.5;
  double norm3 = (double(m_nGal)*double(m_nRan)*double(m_nRan-1))*0.5;
  double norm4 = (double(m_nRan)*double(m_nRan-1)*double(m_nRan-2))/6.;

  cout <<endl<<"---------------------------------------------------"<<endl<<endl;
  cout <<"nGal = "<<m_nGal<<", data.size() = "<<m_data->nObjects()<<endl;
  cout <<"nRan = "<<m_nRan<<endl;
  cout <<"norm(GGG) = "<<norm1<<endl;
  cout <<"norm(GGR) = "<<norm2<<endl;
  cout <<"norm(GRR) = "<<norm3<<endl;
  cout <<"norm(RRR) = "<<norm4<<endl;


  // \zeta(r) in linear bins
  
  for (int i=0; i<m_nbins; i++) 
    if (m_ggg[i]>0 && m_rrr[i]>0) 
      m_zeta[i] = ((m_ggg[i]/norm1)/(m_rrr[i]/norm4))-3.*((m_ggr[i]/norm2)/(m_rrr[i]/norm4))+3.*(((m_grr[i]/norm3)/(m_rrr[i]/norm4)))-1.;
      

  // -------------------------------------------------

  ErrorMsg("Work in progres...");
  
  
  cout <<endl<<"***************************************************"<<endl;
  cout <<"I am now computing the 2-point correlation function"<<endl;
  cout <<"***************************************************"<<endl<<endl;
 
  double _rMIN = m_side_s;
  double _rMAX = (1+m_side_u)*m_side_s*(1+2*m_perc_increase);
  double _logbinSize = 0.05;
  double _binSize = 1.;
  double _cosSize = 0.02;
  bool doGR = 1;
 
  setParameters(_rMIN, _rMAX, _logbinSize, _binSize, _cosSize);
 
  allocate_vectors_xi(doGR);

  measure_xi(dir_output_triplets, tcount);

  vector<double> log_r(m_nlogbins), log_xi_log(m_nlogbins);
  for (int i=0; i<m_nlogbins; i++) {
    log_r[i] = log10(m_rad_log[i]);
    log_xi_log[i] = log10(1.+m_xi_log[i]);
  }

  vector<double> values_interp(m_nbins+2, 0.), xi_real_lin(m_nbins+2, 0.);
  double err = -1, log_xi_real_lin;

  values_interp[0] = log10(m_side_s);
  values_interp[1] = log10(m_side_u*m_side_s);

  for (int i=0; i<m_nbins; i++) {
    double tmp_value = -1.;
    if (m_type_binning=="ang") {
      double theta=((i+0.5)*m_binsize);
      tmp_value = m_side_s*sqrt(1+m_side_u*m_side_u-2*m_side_u*cos(theta));
    }
    else if (m_type_binning=="lin")
      tmp_value = (m_side_s+((i+0.5)*m_binsize));
    values_interp[i+2] = log10(tmp_value);
  }

  string file_2pt = dir_output_2pt+"2ptCorrelation_3pt.dat";
  ofstream fout (file_2pt.c_str()); checkIO(file_2pt, 0);

  for (size_t i=0; i<values_interp.size(); i++) {
    interpolation_extrapolation(values_interp[i], log_r, log_xi_log, "Linear", -1, &log_xi_real_lin, &err);
    xi_real_lin[i] = pow(10.,log_xi_real_lin)-1.;
    //cout << pow(10.,values_interp[i]) << " --- " << xi_real_lin[i] << endl;
    fout << pow(10.,values_interp[i]) << "     " << xi_real_lin[i] << endl;
  }

  fout.clear(); fout.close(); cout <<"I wrote the file "<<file_2pt<<endl;

  
  // ----------- compute the reduced 3-point correlation function -----------

  cout <<endl<<"***********************************************************"<<endl;
  cout <<"I am now computing the reduced 3-point correlation function"<<endl;
  cout <<"***********************************************************"<<endl<<endl;

  for (int i=0; i<m_nbins; i++) 
    m_zeta_red[i] = m_zeta[i]/((xi_real_lin[0]*xi_real_lin[1])+(xi_real_lin[0]*xi_real_lin[i+2])+(xi_real_lin[1]*xi_real_lin[i+2]));
  
}


// ============================================================================


void cosmobl::ThreePointCorrelation::count_triplets (shared_ptr<Catalogue> cat1, ChainMesh_Catalogue &lkList_rMAX1, ChainMesh_Catalogue &lkList_rMAX2, Triplets &tt, bool do3D, bool tcount) 
{
  time_t start; time (&start);
  
  int nObj = cat1->nObjects();

  float fact_count = 100./nObj;
  
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);
 
  if (!do3D) ErrorMsg("Work in progres...");
  //double (Catalogue::*Dist)(int, shared_ptr<Object>) = (do3D) ? &Catalogue::distance : &Catalogue::angsep_xyz;

  shared_ptr<Catalogue> cat2 = lkList_rMAX1.catalogue();
  shared_ptr<Catalogue> cat3 = lkList_rMAX2.catalogue();

  int tid = 0;
#pragma omp parallel private(tid)
  {
    tid = omp_get_thread_num();

    if (tid == 0) {
      int nthreads = omp_get_num_threads();
      cout <<"Number of threads = "<<nthreads<<endl;
    }
    
#pragma omp for schedule(static, 2)  
    for (int i=0; i<nObj; i++) { // loop on the objects of the catalogue
    
      // get the indexes of the objects at r12
      vector<long> close_objects12 = lkList_rMAX1.close_objects(cat1->coordinates(i), -1);
    
      for (auto &&j : close_objects12) { // loop on the objects at r12 

	double r12 = cat1->distance(i, cat2->object(j));
      
	if (m_side_s*(1-m_perc_increase)<r12 && r12<m_side_s*(1+m_perc_increase)) {

	  // get the indexes of objects at r13
	  vector<long> close_objects13 = lkList_rMAX2.close_objects(cat1->coordinates(i), -1);
	
	  for (auto &&k : close_objects13) { // loop on the objects at r13

	    double r13 = cat1->distance(i, cat3->object(k));
	  
	    if (m_side_u*m_side_s*(1-m_perc_increase)<r13 && r13<m_side_u*m_side_s*(1+m_perc_increase)) {

	      double r23 = cat2->distance(j, cat3->object(k));
	    
	      double ww = cat1->weight(i)*cat2->weight(j)*cat3->weight(k);
	    
	      tt.put(r12, r13, r23, ww);
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

  }

  
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<3600) cout <<"   time spent to count the triplets: "<<diff/60<<" minutes"<<endl<<endl;
  else cout <<"   time spent to count the triplets: "<<diff/3600<<" hours"<<endl<<endl;

  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(6); 
}


// ============================================================================


void cosmobl::ThreePointCorrelation::count_triplets_gg (Triplets &GGG, Triplets &GGR, ChainMesh_Catalogue &lkList_rMAX1,  ChainMesh_Catalogue &lkList1_rMAX2, ChainMesh_Catalogue &lkList2_rMAX2, bool do3D, bool tcount) 
{
  time_t start; time (&start);
  
  int nObj = m_data->nObjects();

  float fact_count = 100./nObj;
  
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);
 
  if (!do3D) ErrorMsg("Work in progres...");
  //double (Catalogue::*Dist)(int, shared_ptr<Object>) = (do3D) ? &Catalogue::distance : &Catalogue::angsep_xyz;

  
  for (int i=0; i<nObj; i++) { // loop on the objects of the catalogue    

    vector<long> close_objects12 = lkList_rMAX1.close_objects(m_data->coordinates(i), -1);

    for (auto &&j : close_objects12) { // loop on the objects at r12 

      double r12 = m_data->distance(i, m_data->object(j));
      
      if (m_side_s*(1-m_perc_increase)<r12 && r12<m_side_s*(1+m_perc_increase)) {

	
	// ----- GGG -----
	
	vector<long> close_objects13_GGG = lkList1_rMAX2.close_objects(m_data->coordinates(i), -1);      
	
	for (auto &&k : close_objects13_GGG) { 

	  double r13 = m_data->distance(i, m_data->object(k));
	  
	  if (m_side_u*m_side_s*(1-m_perc_increase)<r13 && r13<m_side_u*m_side_s*(1+m_perc_increase)) {

	    double r23 = m_data->distance(j, m_data->object(k));

	    double ww = m_data->weight(i)*m_data->weight(j)*m_data->weight(k);
	    
	    GGG.put(r12, r13, r23, ww);
	  }
	}


	// ----- GGR -----
	
	vector<long> close_objects13_GGR = lkList2_rMAX2.close_objects(m_data->coordinates(i), -1);      
	
	for (auto &&k : close_objects13_GGR) {

	  double r13 = m_data->distance(i, m_random->object(k));
	  
	  if (m_side_u*m_side_s*(1-m_perc_increase)<r13 && r13<m_side_u*m_side_s*(1+m_perc_increase)) {
	    
	    double r23 = m_data->distance(j, m_random->object(k));

	    double ww = m_data->weight(i)*m_data->weight(j)*m_random->weight(k);
	    
	    GGR.put(r12, r13, r23, ww);

	  }
	}
      }
    }
    

    time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
    if (tcount) { cout <<"\r..."<<float(i)*fact_count<<"% completed  ("<<diff_temp/60<<" minutes)\r"; cout.flush(); }

    if (i==int(nObj*0.25)) cout <<".....25% completed"<<endl;
    if (i==int(nObj*0.5)) cout <<".....50% completed"<<endl;
    if (i==int(nObj*0.75)) cout <<".....75% completed"<<endl;
  }
  
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<3600) cout <<"   time spent to count the triplets: "<<diff/60<<" minutes"<<endl<<endl;
  else cout <<"   time spent to count the triplets: "<<diff/3600<<" hours"<<endl<<endl;

  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(6);
}


// ============================================================================


void cosmobl::ThreePointCorrelation::count_triplets_rr (Triplets &RRG, Triplets &RRR, ChainMesh_Catalogue &lkList_rMAX1,  ChainMesh_Catalogue &lkList1_rMAX2, ChainMesh_Catalogue &lkList2_rMAX2, bool do3D, bool tcount) 
{
  time_t start; time (&start);
  
  int nObj = m_random->nObjects();

  float fact_count = 100./nObj;
  
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);
 
  if (!do3D) ErrorMsg("Work in progres...");
  //double (Catalogue::*Dist)(int, shared_ptr<Object>) = (do3D) ? &Catalogue::distance : &Catalogue::angsep_xyz;

  
  for (int i=0; i<nObj; i++) { // loop on the objects of the catalogue    

    vector<long> close_objects12 = lkList_rMAX1.close_objects(m_random->coordinates(i), -1);

    for (auto &&j : close_objects12) {

      double r12 = m_random->distance(i, m_random->object(j));
      
      if (m_side_s*(1-m_perc_increase)<r12 && r12<m_side_s*(1+m_perc_increase)) {


	// ----- RRG -----
	
	vector<long> close_objects13_RRG = lkList1_rMAX2.close_objects(m_random->coordinates(i), -1);      
	
	for (auto &&k : close_objects13_RRG) {

	  double r13 = m_random->distance(i, m_data->object(k));
	  
	  if (m_side_u*m_side_s*(1-m_perc_increase)<r13 && r13<m_side_u*m_side_s*(1+m_perc_increase)) {

	    double r23 = m_random->distance(j, m_data->object(k));

	    double ww = m_random->weight(i)*m_random->weight(j)*m_data->weight(k);
	    
	    RRG.put(r12, r13, r23, ww);
	  }
	}


	// ----- RRR -----
	
	vector<long> close_objects13_RRR = lkList2_rMAX2.close_objects(m_random->coordinates(i), -1);      
	
	for (auto &&k : close_objects13_RRR) {

	  double r13 = m_random->distance(i, m_random->object(k));
	  
	  if (m_side_u*m_side_s*(1-m_perc_increase)<r13 && r13<m_side_u*m_side_s*(1+m_perc_increase)) {

	    double r23 = m_random->distance(j, m_random->object(k));

	    double ww = m_random->weight(i)*m_random->weight(j)*m_random->weight(k);
	    
	    RRR.put(r12, r13, r23, ww);

	  }
	}
      }
    }
    

    time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
    if (tcount) { cout <<"\r..."<<float(i)*fact_count<<"% completed  ("<<diff_temp/60<<" minutes)\r"; cout.flush(); }

    if (i==int(nObj*0.25)) cout <<".....25% completed"<<endl;
    if (i==int(nObj*0.5)) cout <<".....50% completed"<<endl;
    if (i==int(nObj*0.75)) cout <<".....75% completed"<<endl;
  }
  
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<3600) cout <<"   time spent to count the triplets: "<<diff/60<<" minutes"<<endl<<endl;
  else cout <<"   time spent to count the triplets: "<<diff/3600<<" hours"<<endl<<endl;

  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(6);
}


// ============================================================================


void cosmobl::ThreePointCorrelation::count_triplets_gr (Triplets &GRG, Triplets &GRR, ChainMesh_Catalogue &lkList_rMAX1,  ChainMesh_Catalogue &lkList1_rMAX2, ChainMesh_Catalogue &lkList2_rMAX2, bool do3D, bool tcount) 
{
  time_t start; time (&start);
  
  int nObj = m_data->nObjects();

  float fact_count = 100./nObj;
  
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);
 
  if (!do3D) ErrorMsg("Work in progres...");
  //double (Catalogue::*Dist)(int, shared_ptr<Object>) = (do3D) ? &Catalogue::distance : &Catalogue::angsep_xyz;

  
  for (int i=0; i<nObj; i++) { // loop on the objects of the catalogue    

    vector<long> close_objects12 = lkList_rMAX1.close_objects(m_data->coordinates(i), -1);

    for (auto &&j : close_objects12) {

      double r12 = m_data->distance(i, m_random->object(j));
      
      if (m_side_s*(1-m_perc_increase)<r12 && r12<m_side_s*(1+m_perc_increase)) {


	// ----- GRG -----
	
	vector<long> close_objects13_GRG = lkList1_rMAX2.close_objects(m_data->coordinates(i), -1);      
	
	for (auto &&k : close_objects13_GRG) {

	  double r13 = m_data->distance(i, m_data->object(k));
	  
	  if (m_side_u*m_side_s*(1-m_perc_increase)<r13 && r13<m_side_u*m_side_s*(1+m_perc_increase)) {
	    
	    double r23 = m_random->distance(j, m_data->object(k));

	    double ww = m_data->weight(i)*m_random->weight(j)*m_data->weight(k);
	    
	    GRG.put(r12, r13, r23, ww);
	  }
	}


	// ----- GRR -----
	
	vector<long> close_objects13_GRR = lkList2_rMAX2.close_objects(m_data->coordinates(i), -1);      
	
	for (auto &&k : close_objects13_GRR) {

	  double r13 = m_data->distance(i, m_random->object(k));
	  
	  if (m_side_u*m_side_s*(1-m_perc_increase)<r13 && r13<m_side_u*m_side_s*(1+m_perc_increase)) {

	    double r23 = m_random->distance(j, m_random->object(k));

	    double ww = m_data->weight(i)*m_random->weight(j)*m_random->weight(k);
	    
	    GRR.put(r12, r13, r23, ww);

	  }
	}
      }
    }
    

    time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
    if (tcount) { cout <<"\r..."<<float(i)*fact_count<<"% completed  ("<<diff_temp/60<<" minutes)\r"; cout.flush(); }

    if (i==int(nObj*0.25)) cout <<".....25% completed"<<endl;
    if (i==int(nObj*0.5)) cout <<".....50% completed"<<endl;
    if (i==int(nObj*0.75)) cout <<".....75% completed"<<endl;
  }
  
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<3600) cout <<"   time spent to count the triplets: "<<diff/60<<" minutes"<<endl<<endl;
  else cout <<"   time spent to count the triplets: "<<diff/3600<<" hours"<<endl<<endl;

  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(6);
}


// ============================================================================


void cosmobl::ThreePointCorrelation::count_triplets_rg (Triplets &RGG, Triplets &RGR, ChainMesh_Catalogue &lkList_rMAX1,  ChainMesh_Catalogue &lkList1_rMAX2, ChainMesh_Catalogue &lkList2_rMAX2, bool do3D, bool tcount) 
{
  time_t start; time (&start);
  
  int nObj = m_data->nObjects();

  float fact_count = 100./nObj;
  
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);
 
  if (!do3D) ErrorMsg("Work in progres...");
  //double (Catalogue::*Dist)(int, shared_ptr<Object>) = (do3D) ? &Catalogue::distance : &Catalogue::angsep_xyz;

  
  for (int i=0; i<nObj; i++) { // loop on the objects of the catalogue    

    vector<long> close_objects12 = lkList_rMAX1.close_objects(m_random->coordinates(i), -1);

    for (auto &&j : close_objects12) {
      
      double r12 = m_random->distance(i, m_data->object(j));
      
      if (m_side_s*(1-m_perc_increase)<r12 && r12<m_side_s*(1+m_perc_increase)) {


	// ----- RGG -----
	
	vector<long> close_objects13_RGG = lkList1_rMAX2.close_objects(m_random->coordinates(i), -1);      
	
	for (auto &&k : close_objects13_RGG) {

	  double r13 = m_random->distance(i, m_data->object(k));
	  
	  if (m_side_u*m_side_s*(1-m_perc_increase)<r13 && r13<m_side_u*m_side_s*(1+m_perc_increase)) {

	    double r23 = m_data->distance(j, m_data->object(k));
	    
	    double ww = m_random->weight(i)*m_data->weight(j)*m_data->weight(k);
	    
	    RGG.put(r12, r13, r23, ww);
	  }
	}


	// ----- RGR -----
	
	vector<long> close_objects13_RGR = lkList2_rMAX2.close_objects(m_random->coordinates(i), -1);      
	
	for (auto &&k : close_objects13_RGR) {

	  double r13 = m_random->distance(i, m_random->object(k));
	  
	  if (m_side_u*m_side_s*(1-m_perc_increase)<r13 && r13<m_side_u*m_side_s*(1+m_perc_increase)) {

	    double r23 = m_data->distance(j, m_random->object(k));

	    double ww = m_random->weight(i)*m_data->weight(j)*m_random->weight(k);
	    
	    RGR.put(r12, r13, r23, ww);

	  }
	}
      }
    }
    

    time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
    if (tcount) { cout <<"\r..."<<float(i)*fact_count<<"% completed  ("<<diff_temp/60<<" minutes)\r"; cout.flush(); }

    if (i==int(nObj*0.25)) cout <<".....25% completed"<<endl;
    if (i==int(nObj*0.5)) cout <<".....50% completed"<<endl;
    if (i==int(nObj*0.75)) cout <<".....75% completed"<<endl;
  }
  
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<3600) cout <<"   time spent to count the triplets: "<<diff/60<<" minutes"<<endl<<endl;
  else cout <<"   time spent to count the triplets: "<<diff/3600<<" hours"<<endl<<endl;

  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(6);
}

