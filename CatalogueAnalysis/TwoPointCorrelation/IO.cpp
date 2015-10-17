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
 *  @file CatalogueAnalysis/TwoPointCorrelation/IO.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation used for
 *  Input/Output
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation used for Input/Output
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "TwoPointCorrelation.h"
using namespace cosmobl;


// ============================================================================


void cosmobl::TwoPointCorrelation::write_pairs (vector<double> &PP, vector<double> &PP_lin, vector< vector<double> > &PP_2d, vector< vector<double> > &PP_slog, vector< vector<double> > &PP_coslog, vector< vector<double> > &PP_coslin, string &dir, string &file) 
{  
  string MK = "mkdir -p "+dir;
  if (system (MK.c_str())) {};
  string file_out;
  ofstream fout;
  fout.setf(ios::fixed);

  file_out = dir+file; 
  fout.open (file_out.c_str()); checkIO(file_out,0);

  for (int i=0; i<m_nlogbins; i++) 
    fout <</*setprecision(0)<<*/PP[i]<<endl;
  fout.clear(); fout.close(); cout <<"I wrote the file "<<file_out<<endl;

  file_out = dir+file+"_lin"; 
  fout.open (file_out.c_str()); checkIO(file_out,0);
  
  for (int i=0; i<m_nlinbins; i++) 
    fout <</*setprecision(0)<<*/PP_lin[i]<<endl;
  fout.clear(); fout.close(); cout <<"I wrote the file "<<file_out<<endl;

  file_out = dir+file+"_2d"; 
  fout.open (file_out.c_str()); checkIO(file_out,0);

  for (int i=0; i<m_nlinbins; i++) 
    for (int j=0; j<m_nlinbins; j++) 
      fout <</*setprecision(0)<<*/PP_2d[i][j]<<endl;
  fout.clear(); fout.close(); cout <<"I wrote the file "<<file_out<<endl;

  file_out = dir+file+"_slog"; 
  fout.open (file_out.c_str()); checkIO(file_out,0);
  
  for (int i=0; i<m_nlogbins; i++) 
    for (int j=0; j<m_nlinbins; j++) 
      fout <</*setprecision(0)<<*/PP_slog[i][j]<<endl;
  fout.clear(); fout.close(); cout <<"I wrote the file "<<file_out<<endl;

  file_out = dir+file+"_coslog"; 
  fout.open (file_out.c_str()); checkIO(file_out,0);
  
  for (int i=0; i<m_nlogbins; i++) 
    for (int j=0; j<m_ncosbins; j++) 
      fout <</*setprecision(0)<<*/PP_coslog[i][j]<<endl;
  fout.clear(); fout.close(); cout <<"I wrote the file "<<file_out<<endl;

  file_out = dir+file+"_coslin"; 
  fout.open (file_out.c_str()); checkIO(file_out,0);
 
  for (int i=0; i<m_nlinbins; i++) 
    for (int j=0; j<m_ncosbins; j++) 
      fout <</*setprecision(0)<<*/PP_coslin[i][j]<<endl;
  fout.clear(); fout.close(); cout <<"I wrote the file "<<file_out<<endl;
}


// ============================================================================


void cosmobl::TwoPointCorrelation::read_pairs (vector<double> &PP, vector<double> &PP_lin, vector< vector<double> > &PP_2d, vector< vector<double> > &PP_slog, vector< vector<double> > &PP_coslog, vector< vector<double> > &PP_coslin, vector<string> &dir, string &file) 
{
  if (dir.size()==0)
    ErrorMsg ("Error in cosmobl::TwoPointCorrelation::read_pairs of IO.cpp! dir.size()=0!");
      
  for (unsigned int dd=0; dd<dir.size(); dd++) {
    
    string file_in;
    ifstream fin;
    
    file_in = dir[dd]+file; 
    cout <<"I'm reading the pair file: "<<file_in<<endl;
    fin.open (file_in.c_str()); checkIO(file_in,1);
   
    double pp;
    for (int i=0; i<m_nlogbins; i++) {
      fin >>pp;
      PP[i] += pp;
    }
    fin.clear(); fin.close(); //cout <<"I read the file "<<file_in<<endl;
    
    file_in = dir[dd]+file+"_lin"; 
    fin.open (file_in.c_str()); checkIO(file_in,1);
    
    for (int i=0; i<m_nlinbins; i++) {
      fin >>pp;
      PP_lin[i] += pp;
    }
    fin.clear(); fin.close(); //cout <<"I read the file "<<file_in<<endl;
    
    file_in = dir[dd]+file+"_2d"; 
    fin.open (file_in.c_str()); checkIO(file_in,1);
    
    for (int i=0; i<m_nlinbins; i++) 
      for (int j=0; j<m_nlinbins; j++) {
	fin >>pp;
	PP_2d[i][j] += pp;
      }
    fin.clear(); fin.close(); //cout <<"I read the file "<<file_in<<endl;
    
    file_in = dir[dd]+file+"_slog"; 
    fin.open (file_in.c_str()); checkIO(file_in,1);
    
    for (int i=0; i<m_nlogbins; i++) 
      for (int j=0; j<m_nlinbins; j++) {
	fin >>pp;
	PP_slog[i][j] += pp;
      }
    fin.clear(); fin.close(); //cout <<"I read the file "<<file_in<<endl;

    file_in = dir[dd]+file+"_coslog"; 
    fin.open (file_in.c_str()); checkIO(file_in,1);
 
    for (int i=0; i<m_nlogbins; i++) 
      for (int j=0; j<m_ncosbins; j++) {
	fin >>pp;
	PP_coslog[i][j] += pp;
      }
    fin.clear(); fin.close(); //cout <<"I read the file "<<file_in<<endl;

    file_in = dir[dd]+file+"_coslin"; 
    fin.open (file_in.c_str()); checkIO(file_in,1);
 
    for (int i=0; i<m_nlinbins; i++) 
      for (int j=0; j<m_ncosbins; j++) {
	fin >>pp;
	PP_coslin[i][j] += pp;
      }
    fin.clear(); fin.close(); //cout <<"I read the file "<<file_in<<endl;
  }
 
}


// ============================================================================


void cosmobl::TwoPointCorrelation::write_xi (string &dir_out, int rank) 
{    
  // ---------------------------------------------------
  // ---- check if the vector dimensions is correct ----
  // ---------------------------------------------------

  checkDim(m_rad_log, m_nlogbins, "rad_log");
  checkDim(m_rad_lin, m_nlinbins, "rad_lin");
  checkDim(m_xi_2d_lin, m_nlinbins, m_nlinbins, "xi_2d_lin");
  checkDim(m_error_xi_2d_lin, m_nlinbins, m_nlinbins, "error_xi_2d_lin");
  checkDim(m_xi_2d_loglin, m_nlogbins, m_nlinbins, "xi_2d_loglin");
  checkDim(m_error_xi_2d_loglin, m_nlogbins, m_nlinbins, "error_xi_2d_loglin");
  checkDim(m_xi_coslog, m_nlogbins, m_ncosbins, "xi_coslog");
  checkDim(m_error_xi_coslog, m_nlogbins, m_ncosbins, "error_xi_coslog");
  checkDim(m_xi_coslin, m_nlinbins, m_ncosbins, "xi_coslin");
  checkDim(m_error_xi_coslin, m_nlinbins, m_ncosbins, "error_xi_coslin");
  
  if (!m_read_only) {
    
    // -----------------------------------
    // ---- get the number of objects ----
    // -----------------------------------
    
    string file0 = dir_out+"nObjects";
    ofstream fout0 (file0.c_str()); checkIO(file0, 0);
    
    fout0 <<m_data->nObjects()<<"   "<<m_random->nObjects()<<endl;
    fout0.clear(); fout0.close();
  }

  
  string file1 = dir_out+"csi";
  string file2 = dir_out+"csi_lin";
  string file3 = dir_out+"csi_2d";
  string file4 = dir_out+"xi";
  string file44 = file4+"_SM", file444 = par::DirCosmo+"Func/file1"+conv(rank,par::fINT);
  string file5 = dir_out+"csi_slog";
  string file6 = dir_out+"csi_coslog";
  string file7 = dir_out+"csi_coslin";
  
  ofstream fout1 (file1.c_str()); checkIO(file1,0);
  ofstream fout2 (file2.c_str()); checkIO(file2,0);
  ofstream fout3 (file3.c_str()); checkIO(file3,0);
  ofstream fout4 (file4.c_str()); checkIO(file4,0);
  ofstream fout44 (file444.c_str()); checkIO(file444,0);
  ofstream fout5 (file5.c_str()); checkIO(file5,0);
  ofstream fout6 (file6.c_str()); checkIO(file6,0);
  ofstream fout7 (file7.c_str()); checkIO(file7,0);


  // -------------------------------------------------------
  // ---- 1D correlation function with logarithmic bins ----
  // -------------------------------------------------------

  for (int i=0; i<m_nlogbins; i++) 
    if (m_rad_log[i]>m_rMIN)
      fout1 <<m_rad_log[i]<<"   "<<m_xi_log[i]<<"   "<<m_error_xi_log[i]<<endl;

  fout1.close(); cout <<"I wrote the file: "<<file1<<endl;
  

  // --------------------------------------------------
  // ---- 1D correlation function with linear bins ----
  // --------------------------------------------------
  
  for (int i=0; i<m_nlinbins; i++) 
    if (m_rad_lin[i]>m_rMIN)
      fout2 <<m_rad_lin[i]<<"   "<<m_xi_lin[i]<<"   "<<m_error_xi_lin[i]<<endl;
  
  fout2.close(); cout <<"I wrote the file: "<<file2<<endl;
  

  // -----------------------------------------------------------------------------------
  // ---- 2D correlation function, xi(rp,pi), with linear bins for both coordinates ----
  // -----------------------------------------------------------------------------------

  for (int i=0; i<m_nlinbins; i++) {
    for (int j=0; j<m_nlinbins; j++) {
      double err_rel = (m_xi_2d_lin[i][j]>0) ? m_error_xi_2d_lin[i][j]/m_xi_2d_lin[i][j] : 0.;
      fout3 <<m_rad_lin[i]<<"   "<<m_rad_lin[j]<<"   "<<m_xi_2d_lin[i][j]<<"   "<<m_error_xi_2d_lin[i][j]<<"   "<<err_rel<<endl;
    }
    fout3 <<endl;
  }
  fout3.close(); cout <<"I wrote the file: "<<file3<<endl;
  

  // -------------------------------------------------------------------------------------------------------
  // ---- 2D correlation function, xi(rp,pi), in the 4 quadrants, with linear bins for both coordinates ----
  // -------------------------------------------------------------------------------------------------------

  int ind_i = 0, ind_j = 0;

  for (int i=0; i<m_nlinbins; i++) {
    for (int j=0; j<m_nlinbins; j++) {
      if ((m_nlinbins-(i+1))*m_linbinsz+m_shift_lin>0 && (m_nlinbins-(j+1))*m_linbinsz+m_shift_lin>0) {
	if (m_xi_2d_lin[m_nlinbins-i-1][m_nlinbins-j-1]>0) fout4 <<-m_rad_lin[m_nlinbins-(i+1)]<<"   "<<-m_rad_lin[m_nlinbins-(j+1)]<<"   "<<m_xi_2d_lin[m_nlinbins-i-1][m_nlinbins-j-1]<<endl;
	else fout4 <<-m_rad_lin[m_nlinbins-(i+1)]<<"   "<<-m_rad_lin[m_nlinbins-(j+1)]<<"    0."<<endl;
	fout44 <<ind_i<<"   "<<ind_j++<<"   "<<m_xi_2d_lin[m_nlinbins-i-1][m_nlinbins-j-1]<<endl;
      }
    }
    for (int j=0; j<m_nlinbins; j++) {
      if ((m_nlinbins-(i+1))*m_linbinsz+m_shift_lin>0) {
	if (m_xi_2d_lin[m_nlinbins-i-1][j]>0) fout4 <<-m_rad_lin[m_nlinbins-(i+1)]<<"   "<<m_rad_lin[j]<<"   "<<m_xi_2d_lin[m_nlinbins-i-1][j]<<endl;
	else fout4 <<-m_rad_lin[m_nlinbins-(i+1)]<<"   "<<m_rad_lin[j]<<"    0."<<endl;
	fout44 <<ind_i<<"   "<<ind_j++<<"   "<<m_xi_2d_lin[m_nlinbins-i-1][j]<<endl;
      }      
    }
    fout4 <<endl;
    ind_i++;
  }

  for (int i=0; i<m_nlinbins; i++) {
    for (int j=0; j<m_nlinbins; j++) {
      if ((m_nlinbins-(j+1))*m_linbinsz+m_shift_lin>0) {
	if (m_xi_2d_lin[i][m_nlinbins-j-1]>0) fout4 <<m_rad_lin[i]<<"   "<<-m_rad_lin[m_nlinbins-(j+1)]<<"   "<<m_xi_2d_lin[i][m_nlinbins-j-1]<<endl;
	else fout4 <<m_rad_lin[i]<<"   "<<-m_rad_lin[m_nlinbins-(j+1)]<<"    0."<<endl;
	fout44 <<ind_i<<"   "<<ind_j++<<"   "<<m_xi_2d_lin[i][m_nlinbins-j-1]<<endl;
      }
    }    
    for (int j=0; j<m_nlinbins; j++) {
      if (m_xi_2d_lin[i][j]>0) fout4 <<m_rad_lin[i]<<"   "<<m_rad_lin[j]<<"   "<<m_xi_2d_lin[i][j]<<endl;
      else fout4 <<m_rad_lin[i]<<"   "<<m_rad_lin[j]<<"    0."<<endl;
      fout44 <<ind_i<<"   "<<ind_j++<<"   "<<m_xi_2d_lin[i][j]<<endl;
    }
    fout4 <<endl;
    ind_i++;
  }
  fout4.close(); cout <<"I wrote the file: "<<file4<<endl;
  fout44.close(); 


  // convert the file into the SM format

  string dir_CAT = par::DirCosmo+"Func/";

  int nn = m_nlinbins*2;

  string CONV = dir_CAT+"conv "+conv(nn,par::fINT)+" "+conv(nn,par::fINT)+" "+dir_CAT+" "+conv(rank,par::fINT);
  if (system (CONV.c_str())) {};
  string MV = "cp "+dir_CAT+"file2* "+file44+"; rm "+dir_CAT+"file* -f";
  if (system (MV.c_str())) {};
  cout <<"I wrote the file: "<<file44<<endl;
  

  // -----------------------------------------------------------------
  // ---- 2D correlation function with linear-logarithmic binning ----
  // -----------------------------------------------------------------
      
  for (int i=0; i<m_nlogbins; i++) {
    for (int j=0; j<m_nlinbins; j++) 
      fout5 <<m_rad_log[i]<<"   "<<m_rad_lin[j]<<"   "<<m_xi_2d_loglin[i][j]<<"   "<<m_error_xi_2d_loglin[i][j]<<endl;
    fout5 <<endl;
  }
  fout5.close(); cout <<"I wrote the file: "<<file5<<endl;


  // -----------------------------------------------------------------------------------------------------
  // ---- 2D correlation function in polar coordinates, xi(r,cos(theta)), with logarithmic bins for r ----
  // -----------------------------------------------------------------------------------------------------
      
  for (int i=0; i<m_nlogbins; i++) {
    for (int j=0; j<m_ncosbins; j++) 
      fout6 <<m_rad_log[i]<<"   "<<m_cos_lin[j]<<"   "<<m_xi_coslog[i][j]<<"   "<<m_error_xi_coslog[i][j]<<endl;
    fout6 <<endl;
  }
  fout6.close(); cout <<"I wrote the file: "<<file6<<endl;

  
  // ------------------------------------------------------------------------------------------------
  // ---- 2D correlation function in polar coordinates, xi(r,cos(theta)), with linear bins for r ----
  // ------------------------------------------------------------------------------------------------
      
  for (int i=0; i<m_nlinbins; i++) {
    for (int j=0; j<m_ncosbins; j++) 
      fout7 <<m_rad_lin[i]<<"   "<<m_cos_lin[j]<<"   "<<m_xi_coslin[i][j]<<"   "<<m_error_xi_coslin[i][j]<<endl;
    fout7 <<endl;
  }
  fout7.close(); cout <<"I wrote the file: "<<file7<<endl;

}


// ============================================================================


void cosmobl::TwoPointCorrelation::write_xi (vector< vector<double> > &Xi, string &file, int rank) 
{
  string dir_CAT = par::DirCosmo+"Catalogue/";
  string file1 = dir_CAT+"file1"+conv(rank,par::fINT);
  string file2 = dir_CAT+"file2"+conv(rank,par::fINT);


  // file for SM
  
  ofstream fout_sm (file1.c_str()); checkIO(file1,0);
  
  int ind_i = 0, ind_j = 0;
  for (int i=0; i<m_nlinbins; i++) {
    for (int j=0; j<m_nlinbins; j++) 
      fout_sm <<ind_i<<"   "<<ind_j++<<"   "<<Xi[m_nlinbins-i-1][m_nlinbins-j-1]<<endl;
    for (int j=0; j<m_nlinbins; j++) 
      fout_sm <<ind_i<<"   "<<ind_j++<<"   "<<Xi[m_nlinbins-i-1][j]<<endl;
    ind_i ++;
  }

  for (int i=0; i<m_nlinbins; i++) {
    for (int j=0; j<m_nlinbins; j++)
      fout_sm <<ind_i<<"   "<<ind_j++<<"   "<<Xi[i][m_nlinbins-j-1]<<endl;
    for (int j=0; j<m_nlinbins; j++)
      fout_sm <<ind_i<<"   "<<ind_j++<<"   "<<Xi[i][j]<<endl;
    ind_i ++;
  }
  
  string file_sm = file+"_xsm.dat";

  int nn = m_nlinbins*2; 
  
  string CONV = dir_CAT+"conv "+conv(nn,par::fINT)+" "+conv(nn,par::fINT)+" "+dir_CAT;
  if (system (CONV.c_str())) {};
  string MV = "cp "+file2+" "+file_sm+"; rm -f "+file2+" "+file1;
  if (system (MV.c_str())) {};
  cout <<"I wrote the file: "<<file_sm<<endl;


  // file for gnuplot

  string file_gnu = file+"_xgnu.dat";
  ofstream fout_gnu (file_gnu.c_str()); checkIO(file_gnu,0);
  
  double rp, pi;
  for (int i=0; i<m_nlinbins; i++) {
    rp = m_rad_lin[m_nlinbins-(i+1)];
    for (int j=0; j<m_nlinbins; j++) {
      pi = m_rad_lin[m_nlinbins-(j+1)];
      fout_gnu <<-rp<<"   "<<-pi<<"   "<<Xi[m_nlinbins-i-1][m_nlinbins-j-1]<<endl;
    }
    for (int j=0; j<m_nlinbins; j++) {
      pi = m_rad_lin[j];
      fout_gnu <<-rp<<"   "<<pi<<"   "<<Xi[m_nlinbins-i-1][j]<<endl;
    }
    fout_gnu <<endl;
  }

  for (int i=0; i<m_nlinbins; i++) {
    rp = m_rad_lin[i];
    for (int j=0; j<m_nlinbins; j++) {
      pi = m_rad_lin[m_nlinbins-(j+1)];
      fout_gnu <<rp<<"   "<<-pi<<"   "<<Xi[i][m_nlinbins-j-1]<<endl;
    }    
    for (int j=0; j<m_nlinbins; j++) {
      pi = m_rad_lin[j];
      fout_gnu <<rp<<"   "<<pi<<"   "<<Xi[i][j]<<endl;
    }
    fout_gnu <<endl;
  }
  fout_gnu.clear(); fout_gnu.close(); 
  cout <<"I wrote the file: "<<file_gnu<<endl;
  
}


// ============================================================================


void cosmobl::TwoPointCorrelation::write_multipoles (string &file_log, string &file_lin)
{
  // ---------------------------------------------
  // ---- multipoles with logarithmic binning ----
  // ---------------------------------------------

  ofstream fout_log (file_log.c_str()); checkIO(file_log,0);
  
  for (unsigned int i=0; i<m_xi0_log.size(); i++) 
    fout_log <<m_rad_log[i]<<"   "<<m_xi0_log[i]<<"   "<<m_error_xi0_log[i]<<"   "<<m_xi2_log[i]<<"   "<<m_error_xi2_log[i]<<"   "<<m_xi4_log[i]<<"   "<<m_error_xi4_log[i]<<endl;
  
  fout_log.clear(); fout_log.close(); cout <<"I wrote the file: "<<file_log<<endl;


  // ---------------------------------------------
  // ---- multipoles with linear binning ----
  // ---------------------------------------------

  ofstream fout_lin (file_lin.c_str()); checkIO(file_lin,0);
  
  for (unsigned int i=0; i<m_xi0_lin.size(); i++) 
    fout_lin <<m_rad_lin[i]<<"   "<<m_xi0_lin[i]<<"   "<<m_error_xi0_lin[i]<<"   "<<m_xi2_lin[i]<<"   "<<m_error_xi2_lin[i]<<"   "<<m_xi4_lin[i]<<"   "<<m_error_xi4_lin[i]<<endl;
  
  fout_lin.clear(); fout_lin.close(); cout <<"I wrote the file: "<<file_lin<<endl;
}


// ============================================================================


void cosmobl::TwoPointCorrelation::read_multipoles (string &file_log, string &file_lin)
{
  erase_multipoles ();

  double RR, Xi0, EXi0, Xi2, EXi2, Xi4, EXi4;


  // ---------------------------------------------
  // ---- multipoles with logarithmic binning ----
  // ---------------------------------------------

  ifstream fin_log (file_log.c_str()); checkIO(file_log,1);
  
  while (fin_log >>RR>>Xi0>>EXi0>>Xi2>>EXi2>>Xi4>>EXi4) {
    m_xi0_log.push_back(Xi0);
    m_error_xi0_log.push_back(EXi0);
    m_xi2_log.push_back(Xi2);
    m_error_xi2_log.push_back(EXi2);
    m_xi4_log.push_back(Xi4);
    m_error_xi4_log.push_back(EXi4);
  }
  
  fin_log.clear(); fin_log.close(); cout <<"I read the file: "<<file_log<<endl;


  // ---------------------------------------------
  // ---- multipoles with linear binning ----
  // ---------------------------------------------

  ifstream fin_lin (file_lin.c_str()); checkIO(file_lin,1);
  
  while (fin_lin >>RR>>Xi0>>EXi0>>Xi2>>EXi2>>Xi4>>EXi4) {
    m_xi0_lin.push_back(Xi0);
    m_error_xi0_lin.push_back(EXi0);
    m_xi2_lin.push_back(Xi2);
    m_error_xi2_lin.push_back(EXi2);
    m_xi4_lin.push_back(Xi4);
    m_error_xi4_lin.push_back(EXi4);
  }
  
  fin_lin.clear(); fin_lin.close(); cout <<"I read the file: "<<file_lin<<endl;
}


// ============================================================================


void cosmobl::TwoPointCorrelation::write_wtheta (string &dir_out) 
{    

  // number of objects
  
  string file0 = dir_out+"nObjects";
  ofstream fout0 (file0.c_str()); checkIO(file0,0);

  fout0 <<m_data->nObjects()<<"   "<<m_random->nObjects()<<endl;
  fout0.clear(); fout0.close();

  string file1 = dir_out+"ACF";
  string file2 = dir_out+"ACF_lin";
  
  ofstream fout1 (file1.c_str()); checkIO(file1,0);
  ofstream fout2 (file2.c_str()); checkIO(file2,0);


  // -------------------------------------------------------
  // ---- 1D correlation function with logarithmic bins ----
  // -------------------------------------------------------

  for (int i=0; i<m_nlogthetabins; i++) 
    if (m_theta_log[i]>m_thetaMIN)
      fout1 << m_theta_log[i]*180/par::pi <<"   " << m_wtheta_log[i] << "   " << m_error_wtheta_log[i] << " " << m_gg_theta_log[i] << " " << m_rr_theta_log[i] << " " << m_gr_theta_log[i] << endl;

  fout1.close(); cout <<"I wrote the file: "<<file1<<endl;
  

  // --------------------------------------------------
  // ---- 1D correlation function with linear bins ----
  // --------------------------------------------------
  
  for (int i=0; i<m_nlinthetabins; i++) 
    fout2 << m_theta_lin[i]*180/par::pi << "   " << m_wtheta_lin[i] << "   " << m_error_wtheta_lin[i] << " " << m_gg_theta_lin[i] << " " << m_rr_theta_lin[i] << " " << m_gr_theta_lin[i] << endl;
  
  fout2.close(); cout <<"I wrote the file: "<<file2<<endl;

}


// ============================================================================


void cosmobl::TwoPointCorrelation::write_pairs (vector<double> &PP, vector<double> &PP_lin, string &dir, string &file) 
{  
  string MK = "mkdir -p "+dir;
  if (system (MK.c_str())) {};
  string file_out;
  ofstream fout;
  fout.setf(ios::fixed);

  file_out = dir+file; 
  fout.open (file_out.c_str()); checkIO(file_out,0);

  for (int i=0; i<m_nlogthetabins; i++) 
    fout <<PP[i]<<endl;
  fout.clear(); fout.close(); cout <<"I wrote the file "<<file_out<<endl;

  file_out = dir+file+"_lin"; 
  fout.open (file_out.c_str()); checkIO(file_out,0);
  
  for (int i=0; i<m_nlinthetabins; i++) 
    fout <<PP_lin[i]<<endl;
  fout.clear(); fout.close(); cout <<"I wrote the file "<<file_out<<endl;

}


// ============================================================================


void cosmobl::TwoPointCorrelation::read_pairs (vector<double> &PP, vector<double> &PP_lin, vector<string> &dir, string &file) 
{
  for (unsigned int dd=0; dd<dir.size(); dd++) {
  
    string file_in;
    ifstream fin;
    
    file_in = dir[dd]+file; 
    cout <<"I'm reading the pair file: "<<file_in<<endl;
    fin.open (file_in.c_str()); checkIO(file_in,1);
   
    double pp;
    for (int i=0; i<m_nlogthetabins; i++) {
      fin >>pp;
      PP[i] += pp;
    }
    fin.clear(); fin.close(); //cout <<"I read the file "<<file_in<<endl;
    
    file_in = dir[dd]+file+"_lin"; 

    cout <<"I'm reading the pair file: "<<file_in<<endl;
    fin.open (file_in.c_str()); checkIO(file_in,1);
    
    for (int i=0; i<m_nlinthetabins; i++) {
      fin >>pp;
      PP_lin[i] += pp;
    }
    fin.clear(); fin.close(); //cout <<"I read the file "<<file_in<<endl;
    
  }

}


// ============================================================================


void cosmobl::TwoPointCorrelation::read_bias (string &file_bias, bool proj)
{
  ifstream fin (file_bias.c_str()); checkIO(file_bias,1);
  
  double RR, BL, EBL, BNL, EBNL, BLs, EBLs, BNLs, EBNLs;

  while (fin >>RR>>BL>>EBL>>BNL>>EBNL>>BLs>>EBLs>>BNLs>>EBNLs) {
    if (RR>0) {
      if (proj) {
	m_rr_bias_lin_wp.push_back(RR);
	m_rr_bias_nl_wp.push_back(RR);
	m_bias_lin_wp.push_back(BL);
	m_error_bias_lin_wp.push_back(EBL);
	m_bias_nl_wp.push_back(BNL);
	m_error_bias_nl_wp.push_back(EBNL);
      }
      else {
	m_rr_bias_lin_xi.push_back(RR);
	m_rr_bias_nl_xi.push_back(RR);
	m_bias_lin_xi.push_back(BL);
	m_error_bias_lin_xi.push_back(EBL);
	m_bias_nl_xi.push_back(BNL);
	m_error_bias_nl_xi.push_back(EBNL);
      }
    }
  }
  
  fin.clear(); fin.close();
}


// ============================================================================


void cosmobl::TwoPointCorrelation::write_pairs_subSamples (vector<shared_ptr<Pairs>> pp, int nRegions, string &dir, string &file)
{
  cout << endl << "I'm writing the pairs of each sub-sample (" << dir+file << "_*pairs)..." << endl << endl;
  
  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {};
  
  string ff = dir+file+"_log_pairs";  
  ofstream fout1 (ff.c_str()); checkIO(ff, 0);
  ff = dir+file+"_lin_pairs"; 
  ofstream fout2 (ff.c_str()); checkIO(ff, 0);
  ff = dir+file+"_2d_pairs"; 
  ofstream fout3 (ff.c_str()); checkIO(ff, 0);
  ff = dir+file+"_loglin_pairs"; 
  ofstream fout4 (ff.c_str()); checkIO(ff, 0);
  ff = dir+file+"_coslog_pairs"; 
  ofstream fout5 (ff.c_str()); checkIO(ff, 0);
  ff = dir+file+"_coslin_pairs"; 
  ofstream fout6 (ff.c_str()); checkIO(ff, 0);


  for (int i=0; i<nRegions; i++) {
   
    for (int j=i; j<nRegions; j++) {

      int index = i*nRegions-(i-1)*i/2+j-i;
      
      for (int r1=0; r1<m_nlogbins; r1++)
	if (pp[index]->PPlog(r1)>0)
	  fout1 << i << " " << j << " " << r1 << " " << pp[index]->PPlog(r1) << endl;

      for (int r1=0; r1<m_nlinbins; r1++)
	if (pp[index]->PPlin(r1)>0)
	  fout2 << i << " " << j << " " << r1 << " " <<  pp[index]->PPlin(r1) << endl;

      for (int r1=0; r1<m_nlinbins; r1++) 
	for (int r2=0; r2<m_nlinbins; r2++) 
	  if (pp[index]->PP2d(r1,r2)>0)
	    fout3 << i << " " << j << " " << r1 << " " << r2 << " " << pp[index]->PP2d(r1,r2) << endl;

      for (int r1=0; r1<m_nlogbins; r1++) 
	for (int r2=0; r2<m_nlinbins; r2++) 
	  if (pp[index]->PPslog(r1,r2)>0)
	    fout4 << i << " " << j << " " << r1 << " " << r2 << " " << pp[index]->PPslog(r1,r2) << endl;

      for (int r1=0; r1<m_nlogbins; r1++) 
	for (int r2=0; r2<m_ncosbins; r2++) 
	  if (pp[index]->PPcoslog(r1,r2)>0)
	    fout5 << i << " " << j << " " << r1 << " " << r2 << " " << pp[index]->PPcoslog(r1,r2) << endl;

      for (int r1=0; r1<m_nlinbins; r1++) 
	for (int r2=0; r2<m_ncosbins; r2++) 
	  if (pp[index]->PPcoslin(r1,r2)>0)
	    fout6 << i << " " << j << " " << r1 << " " << r2 << " " << pp[index]->PPcoslin(r1,r2) << endl;
      
    }
  }

  fout1.clear(); fout1.close();
  fout2.clear(); fout2.close();
  fout3.clear(); fout3.close();
  fout4.clear(); fout4.close();
  fout5.clear(); fout5.close();
  fout6.clear(); fout6.close();
}


// ============================================================================


void cosmobl::TwoPointCorrelation::read_pairs_subSamples (vector<shared_ptr<Pairs>> pp, int nRegions, vector<string> &dir, string &file)
{

  for (unsigned int dd=0; dd<dir.size(); dd++) {

    string ff = dir[dd]+file+"_log_pairs"; 
    ifstream fin1(ff.c_str()); checkIO(ff, 1);
    ff = dir[dd]+file+"_lin_pairs"; 
    ifstream fin2(ff.c_str()); checkIO(ff, 1);
    ff = dir[dd]+file+"_2d_pairs";
    ifstream fin3(ff.c_str()); checkIO(ff, 1);
    ff = dir[dd]+file+"_loglin_pairs";
    ifstream fin4(ff.c_str()); checkIO(ff, 1);
    ff = dir[dd]+file+"_coslog_pairs"; 
    ifstream fin5(ff.c_str()); checkIO(ff, 1);
    ff = dir[dd]+file+"_coslin_pairs"; 
    ifstream fin6(ff.c_str()); checkIO(ff, 1);

    double pairs = 0.;
    int I, J, bin1, bin2, index;
    
    while (fin1 >> I >> J >> bin1 >> pairs) {
      index = I*nRegions-(I-1)*I/2+J-I;
      pairs += pp[index]->PPlog(bin1);
      pp[index]->set_PPlog(bin1,pairs);
    }
    fin1.clear(); fin1.close();
    
    pairs = 0.;
    while (fin2 >> I >> J >> bin1 >> pairs) {
      index = I*nRegions-(I-1)*I/2+J-I;
      pairs += pp[index]->PPlin(bin1);
      pp[index]->set_PPlin(bin1,pairs);
    }
    fin2.clear(); fin2.close();

    pairs = 0.;
    while (fin3 >> I >> J >> bin1 >> bin2 >> pairs) {
      index = I*nRegions-(I-1)*I/2+J-I;
      pairs += pp[index]->PP2d(bin1,bin2);
      pp[index]->set_PP2d(bin1,bin2,pairs);
    }
    fin3.clear(); fin3.close();
    
    pairs = 0.;
    while (fin4 >> I >> J >> bin1 >> bin2 >> pairs) {
      index = I*nRegions-(I-1)*I/2+J-I;
      pairs += pp[index]->PPslog(bin1,bin2);
      pp[index]->set_PPslog(bin1,bin2,pairs);
    }
    fin4.clear(); fin4.close();
    
    pairs = 0.;
    while (fin5 >> I >> J >> bin1 >> bin2 >> pairs) {
      index = I*nRegions-(I-1)*I/2+J-I;
      pairs += pp[index]->PPcoslog(bin1,bin2);
      pp[index]->set_PPcoslog(bin1,bin2,pairs);
    }
    fin5.clear(); fin5.close();
    
    pairs = 0.;
    while (fin6 >> I >> J >> bin1 >> bin2 >> pairs) {
      index = I*nRegions-(I-1)*I/2+J-I;
      pairs += pp[index]->PPcoslin(bin1,bin2);
      pp[index]->set_PPcoslin(bin1,bin2,pairs);
    }
    fin6.clear(); fin6.close();
  }

}


// ============================================================================


void cosmobl::TwoPointCorrelation::write_xi_Mocks (string dir_samples)
{
  for (size_t i=0; i<m_twop_mock.size(); i++) {
    string dir = dir_samples+"sample_"+conv(i+1, par::fINT)+"/";
    string cmd = "mkdir -p "+dir; if (system(cmd.c_str())) {}
    m_twop_mock[i]->write_xi(dir);
  } 
}
