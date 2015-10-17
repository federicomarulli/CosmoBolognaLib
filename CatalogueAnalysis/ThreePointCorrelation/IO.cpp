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
 *  @file CatalogueAnalysis/ThreePointCorrelation/IO.cpp
 *
 *  @brief Methods of the class ThreePointCorrelation used for
 *  Input/Output
 *
 *  This file contains the implementation of the methods of the class
 *  ThreePointCorrelation used for Input/Output
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#include "ThreePointCorrelation.h"
using namespace cosmobl;


// ============================================================================


void cosmobl::ThreePointCorrelation::write_triplets (vector<double> &TT, string &dir, string file) 
{  
  string MK = "mkdir -p "+dir;
  if (system (MK.c_str())) {};
  
  string file_out;
  if (m_type_binning=="ang") {
    file_out = dir+file+"_ang"; 
  }
  if (m_type_binning=="lin") {
    file_out = dir+file+"_lin";
  }

  ofstream fout (file_out.c_str()); checkIO (file_out,0);
  for (int i=0; i<m_nbins; i++) 
    fout <<TT[i]<<endl;
  fout.clear(); fout.close(); cout <<"I wrote the file "<<file_out<<endl<<endl;
}


// ============================================================================


void cosmobl::ThreePointCorrelation::read_triplets (vector<double> &TT, vector<string> &dir, string file) 
{
  for (unsigned int dd=0; dd<dir.size(); dd++) {

    string file_in;
    if (m_type_binning=="ang") {
      file_in = dir[dd]+file+"_ang";
    }
    if (m_type_binning=="lin") {
      file_in = dir[dd]+file+"_lin"; 
    }
    ifstream fin (file_in.c_str()); checkIO (file_in,1);

    double tt;
    for (int i=0; i<m_nbins; i++) {
      fin >>tt;
      TT[i] += tt;
    }
    
    fin.clear(); fin.close(); //cout <<"I read the file "<<file_in<<endl<<endl;

  }
}
 

// ============================================================================
 
  
void cosmobl::ThreePointCorrelation::write_Q (string &dir) 
{    
  
  // number of objects
  
  string file0 = dir+"nObjects";
  ofstream fout0 (file0.c_str()); checkIO (file0, 0);
  fout0 <<m_data->nObjects()<<"   "<<m_random->nObjects()<<endl;
  fout0.clear(); fout0.close();


  /*
    output files that will be generated:
    zeta  : Angular 3-points correlation function with linear bins
    zeta_red  : Angular reduced 3-points correlation function with linear bins
  */

  string file1, file2;

  if (m_type_binning=="ang") {
    file1 = dir+"zeta_ang";
    file2 = dir+"zeta_red_ang";
  }
  else if (m_type_binning=="lin") {
    file1 = dir+"zeta_lin";
    file2 = dir+"zeta_red_lin";
  }

  
  ofstream fout1 (file1.c_str()); checkIO (file1,0);
  ofstream fout2 (file2.c_str()); checkIO (file2,0);

  if (m_type_binning=="ang") {  
    double theta;

    // file: zeta (angular bin): 
    // c1: theta (in pi units)
    // c2: zeta(theta)

    for (int i=0; i<m_nbins; i++) {
      theta = ((i+0.5)*m_binsize/par::pi);
      fout1 <<theta<<"   "<<m_zeta[i]<<"   "<<m_error_zeta[i]<<endl;
    }
    fout1.close(); cout <<endl<<"--------------------------------------"<<endl<<"I wrote the file: "<<file1<<endl;

    // file: zeta_red (reduced 3-pt, angular bin):
    // c1: theta (in pi units)
    // c2: zeta_red(theta)
 
    for (int i=0; i<m_nbins; i++) {
      theta = ((i+0.5)*m_binsize/par::pi);
      fout2 <<theta<<"   "<<m_zeta_red[i]<<"   "<<m_error_zeta[i]<<endl;
    }
    fout2.close(); cout <<endl<<"--------------------------------------"<<endl<<"I wrote the file: "<<file2<<endl;
  }

  else if (m_type_binning=="lin") {
    double side;
  
    // file: zeta (linear bin):
    // c1: side (in Mpc/h units)
    // c2: zeta(r)
  
    for (int i=0; i<m_nbins; i++) {
      side = (m_side_s+((i+0.5)*m_binsize));
      fout1 <<side<<"   "<<m_zeta[i]<<"   "<<m_error_zeta[i]<<endl;
    } 
    fout1.close(); cout <<endl<<"--------------------------------------"<<endl<<"I wrote the file: "<<file1<<endl;
    
    // file: zeta_red (reduced 3-pt, linear bin):
    // c1: side (in Mpc/h units)
    // c2: zeta_red(r)

    for (int i=0; i<m_nbins; i++) {
      side = (m_side_s+((i+0.5)*m_binsize));
      fout2 <<side<<"   "<<m_zeta_red[i]<<"   "<<m_error_zeta[i]<<endl;
    }
    fout2.close(); cout <<endl<<"--------------------------------------"<<endl<<"I wrote the file: "<<file2<<endl;
  }

}  
