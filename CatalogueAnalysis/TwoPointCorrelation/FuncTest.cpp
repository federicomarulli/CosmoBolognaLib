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
 *  @file CatalogueAnalysis/TwoPointCorrelation/FuncTest.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation used for testing
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation used for testing
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "TwoPointCorrelation.h"
using namespace cosmobl;


// ============================================================================


void cosmobl::TwoPointCorrelation::count_pairs_direct (Catalogue &cat1, Catalogue &cat2) 
{
  time_t start; time (&start);
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(1);
  
  vector<double> x1 = cat1.var(Var::_XX_), y1 = cat1.var(Var::_YY_), z1 = cat1.var(Var::_ZZ_),
    x2 = cat2.var(Var::_XX_), y2 = cat2.var(Var::_YY_), z2 = cat2.var(Var::_ZZ_); 
  
  float fact_count = 100./float(x1.size());
  

  vector<double> RR; double rr;
  for (int i=0; i<int(x1.size()); i++) { 
    for (int j=i; j<int(x2.size()); j++) { 
      rr = sqrt((x1[i]-x2[j])*(x1[i]-x2[j])+(y1[i]-y2[j])*(y1[i]-y2[j])+(z1[i]-z2[j])*(z1[i]-z2[j]));
      if (log10(m_rMIN)<log10(rr) && rr<m_rMAX_eff) RR.push_back(rr);
    }
    time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
    cout <<"\r..."<<float(i)*fact_count<<"% completed  ("<<diff_temp/60<<" minutes)\r"; cout.flush();
  }
  
  cout.precision(6); 
  double lgr1, lgr2, logbinsz_inv = 1./m_logbinsz;
  int num;
   
  for (int i=0; i<m_nlogbins; i++) {
    lgr1 = (i*m_logbinsz+m_shift_log+log10(m_rMIN))-m_logbinsz*0.5;
    lgr2 = (i*m_logbinsz+m_shift_log+log10(m_rMIN))+m_logbinsz*0.5;
    num = 0;    

    if (lgr2>log10(m_rMAX_eff)) {
      string Err = "Error in cosmobl::TwoPointCorrelation::count_pairs_direct of FuncTest.cpp: r2 = " + conv(pow(10.,lgr2),par::fDP3) + " >rMAX_eff = " + conv(m_rMAX_eff,par::fDP3);
      ErrorMsg(Err);
    }

    for (unsigned int k=0; k<RR.size(); k++) {
      if (lgr1<=log10(RR[k]) && log10(RR[k])<lgr2) num ++;
      if (int((log10(RR[k])-log10(m_rMIN))*logbinsz_inv)==i) 
	if (pow(10.,lgr1)>RR[k] || RR[k]>pow(10.,lgr2)) ErrorMsg("Error in cosmobl::TwoPointCorrelation::count_pairs_direct of FuncTest.cpp!");
    }

    cout.setf(ios::fixed); cout <<setprecision(9)<<" "<<pow(10,i*m_logbinsz+m_shift_log+log10(m_rMIN))<<"   "<<setprecision(0)<<num*2<<endl;
  }

  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<3600) cout <<"   time spent to count the pairs: "<<diff/60<<" minutes"<<endl<<endl;
  else cout <<"   time spent to count the pairs: "<<diff/3600<<" hours"<<endl<<endl;

  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(6);  
  
}


