/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli                         *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file CatalogueAnalysis/ModelTwoPointCorrelation/IO.cpp
 *
 *  @brief Methods of the class ModelTwoPointCorrelation used for
 *  Input/Output
 *
 *  This file contains the implementation of the methods of the class
 *  ModelTwoPointCorrelation used for Input/Output
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "ModelTwoPointCorrelation.h"
using namespace cosmobl;


// ============================================================================


void cosmobl::ModelTwoPointCorrelation::write_xi2D_DispersionModelXiMeasured (string &file, double &beta, double &sigma12, Cosmology &cosm, double &redshift, bool &NL, int &XiReal, double &pimax, string &dir_realxi, double &rApp, int &FV, bool bias_nl, double bA, double v_min, double v_max, int step_v, int rank)
{ 
  // ----- compute the real space xi(r) ----- 

  if (m_TwoP->sizeof_xi_real_lin()==0) {
    m_TwoP->derive_real_xi (XiReal, pimax, dir_realxi, rApp);

    vector<double> rad_real_lin, xi_real_lin_extr;
    for (int i=0; i<m_TwoP->sizeof_xi_real_lin_extr(); i++) 
      if (m_TwoP->rad_real_lin(i)>0) {
	rad_real_lin.push_back(m_TwoP->rad_real_lin(i));
	xi_real_lin_extr.push_back(m_TwoP->xi_real_lin_extr(i));
      }
  }

  vector< vector<double> > Xi2D (m_TwoP->sizeof_xi_real_lin(), vector<double>(m_TwoP->sizeof_xi_real_lin(),0.));

  m_lim_index_fit.erase(m_lim_index_fit.begin(), m_lim_index_fit.end());

  compute_xi2D_modelXiMeasured (Xi2D, beta, sigma12, cosm, redshift, NL, XiReal, pimax, dir_realxi, rApp, FV, bias_nl, bA, v_min, v_max, step_v); 

  write_map (Xi2D, file, rank);

}


// ============================================================================


void cosmobl::ModelTwoPointCorrelation::write_multipoles_DispersionModelXiMeasured (string &file_lin, string &file_log, double &beta, double &sigma12, Cosmology &cosm, double &redshift, bool &NL, int &XiReal, double &pimax, string &dir_realxi, double &rApp, int &FV, bool effective, double v_min, double v_max, int step_v, int step_cos, __attribute__((unused)) int rank)
{ 
  vector<double> xi0_lin, xi2_lin, xi4_lin, xi0_log, xi2_log, xi4_log;

  if (effective) compute_effective_multipoles_modelXiMeasured (xi0_lin, xi2_lin, xi4_lin, xi0_log, xi2_log, xi4_log, beta, sigma12, cosm, redshift, NL, XiReal, pimax, dir_realxi, rApp, FV, v_min, v_max, step_v); 
  
  else compute_multipoles_modelXiMeasured (xi0_lin, xi2_lin, xi4_lin, xi0_log, xi2_log, xi4_log, beta, sigma12, cosm, redshift, NL, XiReal, pimax, dir_realxi, rApp, FV, v_min, v_max, step_v, step_cos); 

  
  ofstream fout_log (file_log.c_str()); checkIO(file_log,0);
  ofstream fout_lin (file_lin.c_str()); checkIO(file_lin,0);
  
  for (unsigned int i=0; i<xi0_log.size(); i++) 
    fout_log <<m_TwoP->rad_log(i)<<"   "<<xi0_log[i]<<"   "<<xi2_log[i]<<"   "<<xi4_log[i]<<endl;

  for (unsigned int i=0; i<xi0_lin.size(); i++)
    fout_lin <<m_TwoP->rad_lin(i)<<"   "<<xi0_lin[i]<<"   "<<xi2_lin[i]<<"   "<<xi4_lin[i]<<endl;
  
  fout_log.clear(); fout_log.close(); cout <<"I wrote the file: "<<file_log<<endl;
  fout_lin.clear(); fout_lin.close(); cout <<"I wrote the file: "<<file_lin<<endl;
}
 

// ============================================================================


void cosmobl::ModelTwoPointCorrelation::write_xi2D_DispersionModel (string &file, double &f_sigma8, double &bias_sigma8, double &sigma12, Cosmology &cosm, double &redshift, bool &NL, string &author, string &Model, int &FV, bool bias_nl, double bA, bool xiType, bool xiNL, double v_min, double v_max, int step_v, int norm, double r_min, double r_max, double k_min, double k_max, double aa, bool GSL, double prec, string file_par, int rank)
{
  cout <<endl<<"I'm computing xi(rp,pi) using the dispersion model..."<<endl;
  
  int dim = m_TwoP->sizeof_xi_lin()-1;
  vector< vector<double> > Xi2D (dim, vector<double>(dim,0.));

  m_lim_index_fit.erase(m_lim_index_fit.begin(), m_lim_index_fit.end());

  compute_xi2D_model (Xi2D, f_sigma8, bias_sigma8, sigma12, cosm, redshift, NL, author, Model, FV, bias_nl, bA, xiType, xiNL, v_min, v_max, step_v, norm, r_min, r_max, k_min, k_max, aa, GSL, prec,  file_par); 
  
  write_map (Xi2D, file, rank);
 
}


// ============================================================================


void cosmobl::ModelTwoPointCorrelation::write_xi2D_CWModel (string &file, double &beta, double &bias_lin, double &bA, double &sigmav0, double &cmu, double &cs1, double &cs2, Cosmology &cosm, double &redshift, string &Model, bool BAO, bool xiType, bool xiNL, double r_min, double r_max, double v_min, double v_max, int step_v, double k_min, double k_max, double x_min, double x_max, int step_x, double aa, bool GSL, double prec, string file_par, int rank)
{
  cout <<endl<<"I'm computing xi(rp,pi) using the Chuang&Wang model..."<<endl;

  string author = "CAMB";
  double k_star = cosm.k_star (author, redshift, Model, k_max, file_par); 

  string author1 = "EisensteinHu"; 
  string author2 = "CAMB";
  vector<double> rr1, Xi1, rr2, Xi2, Xi1_, Xi1__, Xi2_, Xi2__;
  cosm.get_xi (rr1, Xi1, author1, redshift, Model, xiType, k_star, xiNL, 1, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
  cosm.get_xi (rr2, Xi2, author2, redshift, Model, xiType, k_star, xiNL, 1, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
  cosm.get_barred_xi (rr1, Xi1, Xi1_, Xi1__, author1, redshift, xiType, k_star, xiNL, 1, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
  cosm.get_barred_xi (rr2, Xi2, Xi2_, Xi2__, author2, redshift, xiType, k_star, xiNL, 1, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);

  vector< vector<double> > Xi2D (m_TwoP->sizeof_xi_real_lin(), vector<double>(m_TwoP->sizeof_xi_real_lin(),0.));

  for (int i=m_lim_index_fit[0]; i<m_lim_index_fit[1]; i++) {
    double rp = m_TwoP->rad_lin(i); 
    for (int j=m_lim_index_fit[2]; j<m_lim_index_fit[3]; j++) {
      double pi = m_TwoP->rad_lin(j);

      Xi2D[i][j] = cosm.xi2D_CW (rp, pi, beta, bias_lin, bA, sigmav0, cmu, cs1, cs2, redshift, rr1, Xi1, rr2, Xi2, Xi1_, Xi1__, Xi2_, Xi2__, Model, BAO, xiType, k_star, xiNL, r_min, r_max, v_min, v_max, step_v, k_min, k_max, x_min, x_max, step_x, aa, GSL, prec, file_par); 
      
    }
  }

  write_map (Xi2D, file, rank);
}


// ============================================================================


void cosmobl::ModelTwoPointCorrelation::write_map (vector< vector<double> > Xi, string &file, int rank)
{
  // ----- file for SM ----- 

  string dir_CAT = par::DirCosmo+"Func/";

  string file1 = dir_CAT+"file1"+conv(rank,par::fINT);
  string file2 = dir_CAT+"file2"+conv(rank,par::fINT);

  ofstream fout_sm (file1.c_str()); checkIO (file1,0);
  
  int ind_i = 0, ind_j = 0;
  for (unsigned int i=0; i<Xi.size(); i++) {
    for (unsigned int j=0; j<Xi.size(); j++) 
      fout_sm <<ind_i<<"   "<<ind_j++<<"   "<<Xi[Xi.size()-i-1][Xi.size()-j-1]<<endl;
    for (unsigned int j=0; j<Xi.size(); j++) 
      fout_sm <<ind_i<<"   "<<ind_j++<<"   "<<Xi[Xi.size()-i-1][j]<<endl;
    ind_i ++;
  }

  for (unsigned int i=0; i<Xi.size(); i++) {
    for (unsigned int j=0; j<Xi.size(); j++)
      fout_sm <<ind_i<<"   "<<ind_j++<<"   "<<Xi[i][Xi.size()-j-1]<<endl;
    for (unsigned int j=0; j<Xi.size(); j++)
      fout_sm <<ind_i<<"   "<<ind_j++<<"   "<<Xi[i][j]<<endl;
    ind_i ++;
  }

  string file_sm = file+"_xsm.dat";

  int nn = Xi.size()*2; 
  
  string CONV = dir_CAT+"conv "+conv(nn,par::fINT)+" "+conv(nn,par::fINT)+" "+dir_CAT+" "+conv(rank,par::fINT);
  if (system (CONV.c_str())) {};
  string MV = "cp "+file2+" "+file_sm+"; rm -f "+file2+" "+file1;
  if (system (MV.c_str())) {};
  cout <<"I wrote the file: "<<file_sm<<endl;


  // ----- file for gnuplot ----- 

  string file_gnu = file+"_xgnu.dat";
  ofstream fout_gnu (file_gnu.c_str()); checkIO (file_gnu,0);
  
  double rp, pi;
  for (unsigned int i=0; i<Xi.size(); i++) {
    rp = (Xi.size()-(i+1))*m_TwoP->linbinsz()+m_TwoP->shift_lin()+m_TwoP->rMIN();
    for (unsigned int j=0; j<Xi.size(); j++) {
      pi = (Xi.size()-(j+1))*m_TwoP->linbinsz()+m_TwoP->shift_lin()+m_TwoP->rMIN();
      fout_gnu <<-rp<<"   "<<-pi<<"   "<<Xi[Xi.size()-i-1][Xi.size()-j-1]<<endl;
    }
    for (unsigned int j=0; j<Xi.size(); j++) {
      pi = j*m_TwoP->linbinsz()+m_TwoP->shift_lin()+m_TwoP->rMIN();
      fout_gnu <<-rp<<"   "<<pi<<"   "<<Xi[Xi.size()-i-1][j]<<endl;
    }
    fout_gnu <<endl;
  }

  for (unsigned int i=0; i<Xi.size(); i++) {
    rp = i*m_TwoP->linbinsz()+m_TwoP->shift_lin()+m_TwoP->rMIN();
    for (unsigned int j=0; j<Xi.size(); j++) {
      pi = (Xi.size()-(j+1))*m_TwoP->linbinsz()+m_TwoP->shift_lin()+m_TwoP->rMIN();
      fout_gnu <<rp<<"   "<<-pi<<"   "<<Xi[i][Xi.size()-j-1]<<endl;
    }    
    for (unsigned int j=0; j<Xi.size(); j++) {
      pi = j*m_TwoP->linbinsz()+m_TwoP->shift_lin()+m_TwoP->rMIN();
      fout_gnu <<rp<<"   "<<pi<<"   "<<Xi[i][j]<<endl;
    }
    fout_gnu <<endl;
  }
  fout_gnu.clear(); fout_gnu.close(); 
  cout <<"I wrote the file: "<<file_gnu<<endl;
}

