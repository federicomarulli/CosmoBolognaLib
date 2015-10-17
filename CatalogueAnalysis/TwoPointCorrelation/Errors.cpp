/*******************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo *
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
 *  @file CatalogueAnalysis/TwoPointCorrelation/Errors.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation used to compute
 *  the errors and the covariance matrix
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation used to computate the errors and the
 *  covariance matrix
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "TwoPointCorrelation.h"

using namespace cosmobl;


// ============================================================================


double TwoPointCorrelation::Error (double &GG, double &RR) 
{
  double normGG = 2./(m_nGal*(m_nGal-1.));
  double normRR = 2./(m_nRan*(m_nRan-1.));
  
  double GGn = GG*normGG;
  double RRn = RR*normRR;
  double xi = GGn/RRn;

  double fact = normRR/normGG*RR*(1.+xi)+4./m_nGal*pow(normRR*RR/normGG*(1.+xi),2.);

  double ERR = (RR>0 && fact>0) ? normGG/(normRR*RR)*sqrt(fact) : 1.e30;

  ERR *= sqrt(3); // check!!!
 
  return ERR;
}


// ============================================================================


double TwoPointCorrelation::Error (double &GG, double &RR, double &GR) 
{
  double normGG = 2./(m_nGal*(m_nGal-1.));
  double normRR = 2./(m_nRan*(m_nRan-1.));
  double normGR = 1./(m_nGal*m_nRan);

  double GGn = GG*normGG;
  double RRn = RR*normRR;
  double GRn = GR*normGR;
  double xi = max(-0.999,GGn/RRn-2.*GRn/RRn+1.); // check!!!

  double fact = normRR/normGG*RR*(1.+xi)+4./m_nGal*pow(normRR*RR/normGG*(1.+xi),2.);
  
  double ERR = (RR>0 && fact>0) ? normGG/(normRR*RR)*sqrt(fact) : 1.e30;

  ERR *= sqrt(3); // check!!!

  return ERR;
}


// ============================================================================


void TwoPointCorrelation::measure_covariance_1D (vector<string> &file_in, vector<double> &mean, vector< vector<double> > &covariance, string file_out)
{
  vector<double> rad;
  covariance_matrix(file_in, rad,mean, covariance);
 
  if (file_out!="NULL") {
    ofstream fout(file_out.c_str()); checkIO(file_out, 0);
    
    for (unsigned int i=0; i<covariance.size(); i++) {
      for (unsigned int j=0; j<covariance[i].size(); j++) 
	fout <<rad[i]<<"   "<<rad[j]<<"   "<<covariance[i][j]<<endl;
      fout <<endl;
    }
    fout.clear(); fout.close(); cout <<"I wrote the file: "<<file_out<<endl;
  }

}


// ============================================================================


void cosmobl::TwoPointCorrelation::get_covariance (string &dir_out_covariance, bool &doJK, string suffix)
{
  string MKDIR = "mkdir -p "+dir_out_covariance; if (system(MKDIR.c_str())) {}
  cout << "---> " << MKDIR << endl;
  
  int nXi = m_twop_mock.size();
  vector<vector<double>> xi_log_tot(nXi,vector<double>(m_nlogbins, 0)),
    xi_lin_tot(nXi,vector<double>(m_nlinbins, 0)),
    xi_proj_tot(nXi,vector<double>(m_nlogbins, 0));

  vector<vector<vector<double>>> xi_2d_tot(nXi,vector<vector<double>> (m_nlinbins, vector<double> (m_nlinbins, 0)));

  for (int i=0; i<nXi; i++) {

    for (int j=0; j<m_nlogbins; j++) 
      xi_log_tot[i][j] = m_twop_mock[i]->xi_log(j);
    
    for (int j=0; j<m_nlinbins; j++) {
      xi_lin_tot[i][j] = m_twop_mock[i]->xi_lin(j);

      for (int k=0; k<m_nlinbins; k++)
	xi_2d_tot[i][j][k] = m_twop_mock[i]->xi_2d_lin(j,k);
    }
  }

  ofstream fout;
  string fileout;
  
  
  // covariance log

  fileout = (suffix!="NULL") ? dir_out_covariance+"covariance_log_"+suffix : dir_out_covariance+"covariance_log";
  fout.open(fileout.c_str()); checkIO(fileout, 0);
  vector<vector<double>> cov_log;
  covariance_matrix(xi_log_tot, cov_log, doJK);

  for (int i=0; i<m_nlogbins; i++) {
    m_error_xi_log[i] = pow(cov_log[i][i],0.5);
    for (int j=0; j<m_nlogbins; j++) 
      fout << rad_log(i) << " " << rad_log(j) << " " << cov_log[i][j] << endl;
    fout << endl;
  }
  fout.clear(); fout.close(); cout << "I wrote the file: " << fileout << endl;


  //covariance lin
  
  fileout = (suffix!="NULL") ? dir_out_covariance+"covariance_lin_"+suffix : dir_out_covariance+"covariance_lin";
  fout.open(fileout.c_str()); checkIO(fileout, 0);
  vector<vector<double> > cov_lin;
  covariance_matrix(xi_lin_tot, cov_lin, doJK);

  for (int i=0; i<m_nlinbins; i++) {
    m_error_xi_lin[i] = sqrt(cov_lin[i][i]);
    for (int j=0; j<m_nlinbins; j++)
      fout << rad_lin(i) << " " << rad_lin(j) << " " << cov_lin[i][j] << endl;
  }
  fout.clear(); fout.close(); cout << "I wrote the file: " << fileout << endl;


  // error xi_2d

  for (int i=0; i<m_nlinbins; i++) 
    for (int j=0; j<m_nlinbins; j++) {
      
      vector<double> vect;
      for (int k=0; k<nXi; k++)
	vect.push_back(xi_2d_tot[k][i][j]);      

      m_error_xi_2d_lin[i][j] = Sigma(vect);
    }
}

