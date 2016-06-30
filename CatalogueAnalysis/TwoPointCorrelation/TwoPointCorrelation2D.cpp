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
 *  @file CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation2D.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation2D used to
 *  measure the monopole of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation2D used to measure the monopole of the
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "TwoPointCorrelation2D.h"
#include "Data2D.h"

using namespace cosmobl;
using namespace catalogue;
using namespace chainmesh;
using namespace data;
using namespace pairs;
using namespace twopt;


// ===================================================================================================


void cosmobl::twopt::TwoPointCorrelation2D::write_pairs (const shared_ptr<pairs::Pair> PP, const string dir, const string file) const 
{  
  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {};
  
  string file_out = dir+file;
  ofstream fout (file_out.c_str()); checkIO(file_out, 0);
  
  fout.setf(ios::fixed);

  for (int i=0; i<PP->nbins_D1(); i++)
    for (int j=0; j<PP->nbins_D2(); j++) 
      fout << setprecision(5) << PP->PP2D(i, j) << endl;
  
  fout.clear(); fout.close(); cout << "I wrote the file " << file_out << endl;
  
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation2D::read_pairs (shared_ptr<pairs::Pair> PP, const vector<string> dir, const string file) const
{
  if (dir.size()==0)
    ErrorMsg ("Error in cosmobl::twopt::TwoPointCorrelation2D::read_pairs of TwoPointCorrelation2D.cpp! dir.size()=0!");
      
  for (size_t dd=0; dd<dir.size(); dd++) {
        
    string file_in = dir[dd]+file; 
    cout << "I'm reading the pair file: " << file_in << endl;
    
    ifstream fin(file_in.c_str()); checkIO(file_in, 1);
   
    double pp;
    for (int i=0; i<PP->nbins_D1(); i++) 
      for (int j=0; j<PP->nbins_D2(); j++) {
	fin >>pp;
	PP->add_PP2D(i, j, pp);
      }
    
    fin.clear(); fin.close(); cout << "I read the file " << file_in << endl;
  }
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation2D::write_pairs (const vector<shared_ptr<pairs::Pair> > PP, const string dir, const string file) const 
{  

  int nRegions = m_data->get_region_list().size();

  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {};
  
  string file_out = dir+file;  
  ofstream fout (file_out.c_str()); checkIO(file_out, 0);
  fout.setf(ios::fixed);

  for (int i=0; i<nRegions; i++) {
   
    for (int j=i; j<nRegions; j++) {

      int index = i*nRegions-(i-1)*i/2+j-i;
      
      for (int r1=0; r1<PP[index]->nbins_D1(); r1++)
	for (int r2=0; r2<PP[index]->nbins_D2(); r2++)
	  if(PP[index]->PP2D(r1,r2)>0)
	    fout << i << " " << j << " " << r1 << " " << r2 << " " << setprecision(5) << PP[index]->PP2D(r1,r2) << endl;

    }
  }

  fout.clear(); fout.close();
  
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation2D::read_pairs (vector<shared_ptr<pairs::Pair> > PP, const vector<string> dir, const string file) const
{

  int nRegions = m_data->get_region_list().size();
  for (unsigned int dd=0; dd<dir.size(); dd++) {

    string ff = dir[dd]+file; 
    cout << "I'm reading the pair file: " << ff << endl;
    ifstream fin(ff.c_str()); checkIO(ff, 1);

    int I, J, bin1, bin2, index;
    double pairs;

    while (fin >> I >> J >> bin1 >> bin2 >> pairs) {
      index = I*nRegions-(I-1)*I/2+J-I;
      PP[index]->add_PP2D(bin1,bin2,pairs);
    }
    fin.clear(); fin.close();

    fin.clear(); fin.close(); cout << "I read the file " << ff << endl;

  }
}


// ============================================================================


shared_ptr<Data> cosmobl::twopt::TwoPointCorrelation2D::NaturalEstimatorTwoP (const shared_ptr<pairs::Pair> dd, const shared_ptr<pairs::Pair> rr, const int nData, const int nRandom) {

  vector<double> scale_D1,scale_D2;
  vector<vector<double> > xi,error;
  scale_D1.resize(m_dd->nbins_D1()); 
  scale_D2.resize(m_dd->nbins_D2()); 

  xi.resize(m_dd->nbins_D1(),vector<double>(m_dd->nbins_D2(),0));
  error.resize(m_dd->nbins_D1(),vector<double>(m_dd->nbins_D2(),0));

  double norm = double(nRandom)*double(nRandom-1)/(double(nData)*double(nData-1));

  for (int i=0; i<dd->nbins_D1(); i++) {
    scale_D1[i] = dd->scale_D1(i);
    for (int j=0; j<dd->nbins_D2(); j++) {
      scale_D2[j] = dd->scale_D2(j);

      xi[i][j] = -1.;
      error[i][j] = 1000.;

      if (dd->PP2D(i,j)>0 && rr->PP2D(i,j)>0) {

	xi[i][j] = max(-1., norm*dd->PP2D(i,j)/rr->PP2D(i,j)-1.);

	error[i][j]= PoissonError(dd->PP2D(i,j), rr->PP2D(i,j), 0, nData, nRandom); 
      }
    }
  }

  return move(unique_ptr<Data2D>(new Data2D(scale_D1,scale_D2,xi,error)));
}


// ============================================================================


shared_ptr<Data> cosmobl::twopt::TwoPointCorrelation2D::LandySzalayEstimatorTwoP (const shared_ptr<pairs::Pair> dd, const shared_ptr<pairs::Pair> rr, shared_ptr<pairs::Pair> dr, int nData, int nRandom) {

  vector<double> scale_D1,scale_D2;
  vector<vector<double> > xi,error;
  scale_D1.resize(m_dd->nbins_D1()); 
  scale_D2.resize(m_dd->nbins_D2()); 

  xi.resize(m_dd->nbins_D1(),vector<double>(m_dd->nbins_D2(),0));
  error.resize(m_dd->nbins_D1(),vector<double>(m_dd->nbins_D2(),0));
  
  double norm1 = double(nRandom)*double(nRandom-1)/(double(nData)*double(nData-1));
  double norm2 = double(nRandom-1)/double(nData);

  for (int i=0; i<dd->nbins_D1(); i++) {
    scale_D1[i] = dd->scale_D1(i);
    for (int j=0; j<dd->nbins_D2(); j++) {
      scale_D2[j] = dd->scale_D2(j);

      xi[i][j] = -1.;
      error[i][j] = 1000.;

      if (dd->PP2D(i,j)>0 && rr->PP2D(i,j)>0) {

	xi[i][j] = max(-1., norm1*dd->PP2D(i,j)/rr->PP2D(i,j)-norm2*dr->PP2D(i,j)/rr->PP2D(i,j)+1.);

	error[i][j]= PoissonError(dd->PP2D(i,j), rr->PP2D(i,j), dr->PP2D(i,j),nData,nRandom); 
      }
    }
  }

  return move(unique_ptr<Data2D>(new Data2D(scale_D1,scale_D2,xi,error)));
}


// ============================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation2D::XiJackknife(const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr)
{
  vector<long> region_list = m_data->get_region_list();
  int nRegions = region_list.size();

  vector<shared_ptr<Data> > data;

  for(int i=0;i<nRegions;i++){
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2(), m_rr->angularUnits(), m_rr->angularWeight());

    double nData_SS = m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    double nRandom_SS = m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);

    vector<int> w(nRegions, 1);
    w[i] = 0;

    for (int j=0; j<nRegions; j++) {
      for (int k=j; k<nRegions; k++) {
	int index = j*nRegions-(j-1)*j/2+k-j;
	double ww = w[j]*w[k];
	if (ww>0) {	
	  for (int bin1=0;bin1<dd_SS->nbins_D1();bin1++) {
	    for (int bin2=0;bin2<dd_SS->nbins_D2();bin2++) {
	      dd_SS->add_PP2D(bin1 ,bin2, dd[index]->PP2D(bin1,bin2));
	      rr_SS->add_PP2D(bin1,bin2 , rr[index]->PP2D(bin1,bin2));
	    }
	  }
	}
      }
    }

    data.push_back(move(NaturalEstimatorTwoP(dd_SS, rr_SS, nData_SS, nRandom_SS)));
  }
  return data;

}


// ============================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation2D::XiJackknife(const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr)
{
  vector<long> region_list = m_data->get_region_list();
  int nRegions = region_list.size();

  vector<shared_ptr<Data> > data;

  for(int i=0;i<nRegions;i++){

    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2(), m_rr->angularUnits(), m_rr->angularWeight());
    auto dr_SS = Pair::Create(m_dr->pairType(), m_dr->sMin_D1(), m_dr->sMax_D1(), m_dr->nbins_D1(), m_dr->shift_D1(), m_dr->sMin_D2(), m_dr->sMax_D2(), m_dr->nbins_D2(), m_dr->shift_D2(), m_dr->angularUnits(), m_dr->angularWeight());

    double nData_SS = m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    double nRandom_SS = m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);

    vector<int> w(nRegions, 1);
    w[i] = 0;
    
    for (int j=0; j<nRegions; j++) {
      for (int k=j; k<nRegions; k++) {
	int index = j*nRegions-(j-1)*j/2+k-j;
	double ww = w[j]*w[k];
	if (ww>0) {	
	  for (int bin1=0;bin1<dd_SS->nbins_D1();bin1++) {
	    for (int bin2=0;bin2<dd_SS->nbins_D2();bin2++) {
	      dd_SS->add_PP2D(bin1, bin2, dd[index]->PP2D(bin1,bin2));
	      rr_SS->add_PP2D(bin1, bin2, rr[index]->PP2D(bin1,bin2));
	      dr_SS->add_PP2D(bin1, bin2, dr[index]->PP2D(bin1,bin2));
	    }
	  }
	}
      }
    }
   
    data.push_back(move(LandySzalayEstimatorTwoP(dd_SS, rr_SS, dr_SS, nData_SS, nRandom_SS)));
    /*
    for(int aa=0;aa<dd_SS->nbins_D1();aa++)
      for(int bb=0;bb<dd_SS->nbins_D2();bb++)
	cout << dd_SS->PP2D(aa,bb) << " " << rr_SS->PP2D(aa,bb) << " " << dr_SS->PP2D(aa,bb) << endl;
    exit(0);
*/
  }
  return data;

}


// ============================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation2D::XiBootstrap(const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr)
{
  vector<long> region_list = m_data->get_region_list();
  int nRegions = region_list.size();

  vector<shared_ptr<Data> > data;
  vector<double> nData_reg, nRandom_reg;

  for(int i=0;i<nMocks;i++){
    nData_reg.push_back(m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg.push_back(m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
  }

  uniform_int_distribution<int> uni(0, nRegions-1);
  default_random_engine rng;
  int val=3; //See Norberg et al. 2009

  for(int i=0;i<nMocks;i++){
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2(), m_rr->angularUnits(), m_rr->angularWeight());
    
    double nData_SS = 0., nRandom_SS = 0.;

    vector<int> w(nRegions, 0);
    for (int n=0; n<val*nRegions; n++)
      w[uni(rng)] +=1;

    for (int j=0; j<nRegions; j++) {
      nData_SS += w[j]*nData_reg[j];
      nRandom_SS += w[j]*nRandom_reg[j];

      for (int k=j; k<nRegions; k++) {
	int index = j*nRegions-(j-1)*j/2+k-j;
	double ww = (k==j) ? w[k] : w[j]*w[k];
	if (ww>0) {
	  for (int bin1=0;bin1<dd_SS->nbins_D1();bin1++) {
	    for (int bin2=0;bin2<dd_SS->nbins_D2();bin2++) {
	      dd_SS->add_PP2D(bin1 ,bin2, dd[index]->PP2D(bin1,bin2));
	      rr_SS->add_PP2D(bin1,bin2 , rr[index]->PP2D(bin1,bin2));
	    }
	  }
	}
      }
    }

    data.push_back(move(NaturalEstimatorTwoP(dd_SS, rr_SS, nData_SS, nRandom_SS)));
  }
  return data;

}


// ============================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation2D::XiBootstrap(const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr)
{
  vector<long> region_list = m_data->get_region_list();
  int nRegions = region_list.size();

  vector<shared_ptr<Data> > data;
  vector<double> nData_reg, nRandom_reg;

  for(int i=0;i<nMocks;i++){
    nData_reg.push_back(m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg.push_back(m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
  }

  uniform_int_distribution<int> uni(0, nRegions-1);
  default_random_engine rng;
  int val=3; //See Norberg et al. 2009

  for(int i=0;i<nMocks;i++){
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2(), m_rr->angularUnits(), m_rr->angularWeight());
    auto dr_SS = Pair::Create(m_dr->pairType(), m_dr->sMin_D1(), m_dr->sMax_D1(), m_dr->nbins_D1(), m_dr->shift_D1(), m_dr->sMin_D2(), m_dr->sMax_D2(), m_dr->nbins_D2(), m_dr->shift_D2(), m_dr->angularUnits(), m_dr->angularWeight());

    double nData_SS=0, nRandom_SS=0;

    vector<int> w(nRegions, 0);
    for (int n=0; n<val*nRegions; n++)
      w[uni(rng)] +=1;

    for (int j=0; j<nRegions; j++) {
      nData_SS += w[j]*nData_reg[j];
      nRandom_SS += w[j]*nRandom_reg[j];

      for (int k=j; k<nRegions; k++) {
	int index = j*nRegions-(j-1)*j/2+k-j;
	double ww = (k==j) ? w[k] : w[j]*w[k];
	if (ww>0) {
	  for (int bin1=0;bin1<dd_SS->nbins_D1();bin1++) {
	    for (int bin2=0;bin2<dd_SS->nbins_D2();bin2++) {
	      dd_SS->add_PP2D(bin1, bin2, ww*dd[index]->PP2D(bin1,bin2));
	      rr_SS->add_PP2D(bin1, bin2, ww*rr[index]->PP2D(bin1,bin2));
	      dr_SS->add_PP2D(bin1, bin2, ww*dr[index]->PP2D(bin1,bin2));
	    }
	  }
	}
      }
    }

    data.push_back(move(LandySzalayEstimatorTwoP(dd_SS, rr_SS, dr_SS, nData_SS, nRandom_SS)));
  }
  return data;

}
