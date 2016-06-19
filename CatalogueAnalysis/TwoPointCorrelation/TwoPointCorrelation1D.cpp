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
 *  @file CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation1D.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation1D used to
 *  measure the monopole of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation1D used to measure the monopole of the
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation1D.h"
#include "Data1D.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;
using namespace twopt;


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::write_pairs (const shared_ptr<pairs::Pair> PP, const string dir, const string file) const 
{  
  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {};
  
  string file_out = dir+file;
  ofstream fout (file_out.c_str()); checkIO(file_out, 0);
  
  fout.setf(ios::fixed);

  for (int i=0; i<PP->nbins(); i++) 
    fout << setprecision(5) << PP->PP1D(i) << endl;
  
  fout.clear(); fout.close(); cout << "I wrote the file " << file_out << endl;
  
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::read_pairs (shared_ptr<pairs::Pair> PP, const vector<string> dir, const string file) const
{
  if (dir.size()==0)
    ErrorMsg ("Error in cosmobl::twopt::TwoPointCorrelation1D::read_pairs of TwoPointCorrelation1D.cpp! dir.size()=0!");
      
  for (size_t dd=0; dd<dir.size(); dd++) {
        
    string file_in = dir[dd]+file; 
    cout << "I'm reading the pair file: " << file_in << endl;
    
    ifstream fin(file_in.c_str()); checkIO(file_in, 1);
   
    double pp;
    for (int i=0; i<PP->nbins(); i++) {
      fin >>pp;
      PP->add_PP1D(i, pp);
    }
    
    fin.clear(); fin.close(); cout << "I read the file " << file_in << endl;
  }
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::write_pairs (const vector<shared_ptr<pairs::Pair> > PP, const string dir, const string file) const 
{  
  int nRegions = m_data->get_region_list().size();
  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {};
  
  string file_out = dir+file;  
  ofstream fout (file_out.c_str()); checkIO(file_out, 0);
  fout.setf(ios::fixed);

  for (int i=0; i<nRegions; i++) {
   
    for (int j=i; j<nRegions; j++) {

      int index = i*nRegions+j-(i-1)*i/2-i;

      for (int r1=0; r1<PP[index]->nbins(); r1++)
	if(PP[index]->PP1D(r1)>0)
	  fout << i << " " << j << " " << r1 << " " << setprecision(5) << PP[index]->PP1D(r1) << endl;
    }
  }

  fout.clear(); fout.close();
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::read_pairs (vector<shared_ptr<pairs::Pair> > PP, const vector<string> dir, const string file) const
{
  int nRegions = m_data->get_region_list().size();

  for (size_t dd=0; dd<dir.size(); dd++) {

    string ff = dir[dd]+file; 
    cout << "I'm reading the pair file: " << ff << endl;
    ifstream fin(ff.c_str()); checkIO(ff, 1);

    int I, J, bin1, index;
    double pairs;

    while (fin >> I >> J >> bin1 >> pairs) {
      index = I*nRegions-(I-1)*I/2+J-I;
      PP[index]->add_PP1D(bin1,pairs);
    }
    
    fin.clear(); fin.close(); cout << "I read the file " << ff << endl;
  }
}


// ============================================================================


shared_ptr<Data> cosmobl::twopt::TwoPointCorrelation1D::NaturalEstimatorTwoP (const shared_ptr<pairs::Pair> dd, const shared_ptr<pairs::Pair> rr, const int nData, const int nRandom)
{
  vector<double> rad(m_dd->nbins()), xi(m_dd->nbins(), -1.), error(m_dd->nbins(), 1000.);
  
  double norm = double(nRandom)*double(nRandom-1)/(double(nData)*double(nData-1));

  for (int i=0; i<dd->nbins(); i++) {

    rad[i] = dd->scale(i);
    
    if (dd->PP1D(i)>0 && rr->PP1D(i)>0) {
      xi[i] = max(-1., norm*dd->PP1D(i)/rr->PP1D(i)-1.);
      error[i] = PoissonError(dd->PP1D(i), rr->PP1D(i), 0, nData, nRandom); 
    }
    
  }

  return move(unique_ptr<Data1D>(new Data1D(rad, xi, error)));
}


// ============================================================================


shared_ptr<Data> cosmobl::twopt::TwoPointCorrelation1D::LandySzalayEstimatorTwoP (const shared_ptr<pairs::Pair> dd, const shared_ptr<pairs::Pair> rr, const shared_ptr<pairs::Pair> dr, const int nData, const int nRandom)
{
  vector<double> rad(m_dd->nbins()), xi(m_dd->nbins(), -1.), error(m_dd->nbins(), 1000.);
  
  double norm1 = double(nRandom)*double(nRandom-1)/(double(nData)*double(nData-1));
  double norm2 = double(nRandom-1)/double(nData);
  
  for (int i=0; i<dd->nbins(); i++) {
  
    rad[i] = dd->scale(i);
    
    if (dd->PP1D(i)>0 && rr->PP1D(i)>0) {
      xi[i] = max(-1., norm1*dd->PP1D(i)/rr->PP1D(i)-norm2*dr->PP1D(i)/rr->PP1D(i)+1.);
      error[i] = PoissonError(dd->PP1D(i), rr->PP1D(i), dr->PP1D(i), nData, nRandom); 
    }
  }

  return move(unique_ptr<Data1D>(new Data1D(rad, xi, error)));
}


// ============================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation1D::XiJackknife (const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr)
{
  vector<long> region_list = m_data->get_region_list();
  int nRegions =  region_list.size();

  vector<shared_ptr<Data> > data;

  for (int i=0; i<nRegions; i++) {
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->sMin(), m_dd->sMax(), m_dd->nbins(), m_dd->shift(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->sMin(), m_rr->sMax(), m_rr->nbins(), m_rr->shift(), m_rr->angularUnits(), m_rr->angularWeight());

    double nData_SS = m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    double nRandom_SS = m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);

    vector<int> w(nRegions, 1);
    w[i] = 0;

    for (int j=0; j<nRegions; j++) {
      for (int k=j; k<nRegions; k++) {
	int index = j*nRegions-(j-1)*j/2+k-j;
	double ww = w[j]*w[k];
	if (ww>0) {	
	  for (int bin=0;bin<dd_SS->nbins();bin++) {
	    dd_SS->add_PP1D(bin, dd[index]->PP1D(bin));
	    rr_SS->add_PP1D(bin, rr[index]->PP1D(bin));
	  }
	}
      }
    }

    data.push_back(move(NaturalEstimatorTwoP(dd_SS, rr_SS, nData_SS, nRandom_SS)));
  }
  
  return data;
}


// ============================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation1D::XiJackknife (const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr)
{
  vector<long> region_list = m_data->get_region_list();
  int nRegions =  region_list.size();

  vector<shared_ptr<Data> > data;

  for (int i=0; i<nRegions; i++) {
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->sMin(), m_dd->sMax(), m_dd->nbins(), m_dd->shift(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->sMin(), m_rr->sMax(), m_rr->nbins(), m_rr->shift(), m_rr->angularUnits(), m_rr->angularWeight());
    auto dr_SS = Pair::Create(m_dr->pairType(), m_dr->sMin(), m_dr->sMax(), m_dr->nbins(), m_dr->shift(), m_dr->angularUnits(), m_dr->angularWeight());

    double nData_SS = m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);
    double nRandom_SS = m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 1);

    vector<int> w(nRegions, 1);
    w[i] = 0;

    for (int j=0; j<nRegions; j++) {
      for (int k=j; k<nRegions; k++) {
	int index = j*nRegions+k-(j-1)*j/2-j;
	double ww = w[j]*w[k];
	if (ww>0) {	
	  for (int bin=0;bin<dd_SS->nbins();bin++) {
	    dd_SS->add_PP1D(bin, dd[index]->PP1D(bin));
	    rr_SS->add_PP1D(bin, rr[index]->PP1D(bin));
	    dr_SS->add_PP1D(bin, dr[index]->PP1D(bin));
	  }
	}
      }
    }

    data.push_back(move(LandySzalayEstimatorTwoP(dd_SS, rr_SS, dr_SS, nData_SS, nRandom_SS)));
  }
  
  return data;
}


// ============================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation1D::XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr)
{
  vector<long> region_list = m_data->get_region_list();
  int nRegions = region_list.size();

  vector<shared_ptr<Data> > data;
  vector<double> nData_reg, nRandom_reg;

  for (int i=0; i<nMocks; i++) {
    nData_reg.push_back(m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg.push_back(m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
  }

  uniform_int_distribution<int> uni(0, nRegions-1);
  default_random_engine rng;
  int val=1; //See Norberg et al. 2009

  for (int i=0; i<nMocks; i++) {
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->sMin(), m_dd->sMax(), m_dd->nbins(), m_dd->shift(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->sMin(), m_rr->sMax(), m_rr->nbins(), m_rr->shift(), m_rr->angularUnits(), m_rr->angularWeight());

    double nData_SS=0, nRandom_SS=0;

    vector<int> w(nRegions, 0);
    for (int n=0; n<val*nRegions; n++)
      w[uni(rng)] ++;

    for (int j=0; j<nRegions; j++) {
      nData_SS += w[j]*nData_reg[j];
      nRandom_SS += w[j]*nRandom_reg[j];

      for (int k=j; k<nRegions; k++) {
	int index = j*nRegions+k-(j-1)*j/2-j;
	double ww = (k==j) ? w[k] : w[j]*w[k];
	if (ww>0) {
	  for (int bin=0; bin<dd_SS->nbins(); bin++) {
	    dd_SS->add_PP1D(bin, ww*dd[index]->PP1D(bin));
	    rr_SS->add_PP1D(bin, ww*rr[index]->PP1D(bin));
	  }
	}
      }
    }

    data.push_back(move(NaturalEstimatorTwoP(dd_SS, rr_SS, nData_SS, nRandom_SS)));
  }
  return data;

}


// ============================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation1D::XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr)
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
  int val = 1; // see Norberg et al. 2009

  for(int i=0;i<nMocks;i++){
    auto dd_SS = Pair::Create(m_dd->pairType(), m_dd->sMin(), m_dd->sMax(), m_dd->nbins(), m_dd->shift(), m_dd->angularUnits(), m_dd->angularWeight());
    auto rr_SS = Pair::Create(m_rr->pairType(), m_rr->sMin(), m_rr->sMax(), m_rr->nbins(), m_rr->shift(), m_rr->angularUnits(), m_rr->angularWeight());
    auto dr_SS = Pair::Create(m_dr->pairType(), m_dr->sMin(), m_dr->sMax(), m_dr->nbins(), m_dr->shift(), m_dr->angularUnits(), m_dr->angularWeight());

    double nData_SS=0, nRandom_SS=0;

    vector<int> w(nRegions, 0);
    for (int n=0; n<val*nRegions; n++)
      w[uni(rng)] +=1;

    for (int j=0; j<nRegions; j++) {
      nData_SS += w[j]*nData_reg[j];
      nRandom_SS += w[j]*nRandom_reg[j];

      for (int k=j; k<nRegions; k++) {
	int index = j*nRegions+k-(j-1)*j/2-j;
	double ww = (k==j) ? w[k] : w[j]*w[k];
	if (ww>0) {
	  for (int bin=0; bin<dd_SS->nbins(); bin++) {
	    dd_SS->add_PP1D(bin, ww*dd[index]->PP1D(bin));
	    rr_SS->add_PP1D(bin, ww*rr[index]->PP1D(bin));
	    dr_SS->add_PP1D(bin, ww*dr[index]->PP1D(bin));
	  }
	}
      }
    }

    data.push_back(move(LandySzalayEstimatorTwoP(dd_SS, rr_SS, dr_SS, nData_SS, nRandom_SS)));
  }
  
  return data;
}

// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::read_covariance_matrix (const string dir, const string file)
{
  m_dataset->set_covariance(dir+file);
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::write_covariance_matrix (const string dir, const string file) const
{
  cout << "HH" << endl;
  m_dataset->write_covariance(dir, file, "r");
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::compute_covariance_matrix (const vector<shared_ptr<Data>> xi_collection, const bool doJK)
{
  vector<vector<double>> xi;

  for(size_t i=0;i<xi_collection.size();i++)
    xi.push_back(xi_collection[i]->fx());

  vector<vector<double> > cov_mat;
  cosmobl::covariance_matrix(xi, cov_mat, doJK);
  
  m_dataset->set_covariance(cov_mat);
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::compute_covariance_matrix (const vector<string> file_xi, const bool doJK)
{
  vector<double> rad, mean;
  vector<vector<double> > cov_mat;

  cosmobl::covariance_matrix (file_xi, rad, mean, cov_mat, doJK); 
  m_dataset->set_covariance(cov_mat);
}

