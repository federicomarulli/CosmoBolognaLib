/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli                          *
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
 *  @file NumberCounts/NumberCounts2D.cpp
 *
 *  @brief Methods of the class NumberCounts2D 
 *
 *  This file contains the implementation of the methods of the class
 *  NumberCounts2D, used to measure number counts
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "NumberCounts2D.h"
#include "Data2D_extra.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace measure::numbercounts;


// ============================================================================


cbl::measure::numbercounts::NumberCounts2D::NumberCounts2D(const catalogue::Var var1, const BinType bin_type1, const catalogue::Var var2, const BinType bin_type2, const catalogue::Catalogue data, const size_t nbins1, const size_t nbins2, const double minVar1, const double maxVar1, const double minVar2, const double maxVar2, const double shift1, const double shift2, const glob::HistogramType hist_type, const double fact)
{
  m_Var1 = var1;
  m_Var2 = var2;

  m_HistogramType = hist_type;
  m_fact = fact;

  set_data(data);

  m_histogram = make_shared<glob::Histogram2D> (glob::Histogram2D()); 

  double _minVar1 = (minVar1>par::defaultDouble) ? minVar1 : m_data->Min(m_Var1)*0.999;
  double _maxVar1 = (maxVar1>par::defaultDouble) ? maxVar1 : m_data->Max(m_Var1)*1.001;
  double _minVar2 = (minVar2>par::defaultDouble) ? minVar2 : m_data->Min(m_Var2)*0.999;
  double _maxVar2 = (maxVar2>par::defaultDouble) ? maxVar2 : m_data->Max(m_Var2)*1.001;

  m_histogram->set(nbins1, nbins2, _minVar1, _maxVar1, _minVar2, _maxVar2, shift1, shift2, bin_type1, bin_type2);
}


// ============================================================================


shared_ptr<data::Data> cbl::measure::numbercounts::NumberCounts2D::m_measurePoisson ()
{
  auto histogram =  make_shared<glob::Histogram2D> (glob::Histogram2D ());
  histogram->set(m_histogram->nbins1(), m_histogram->nbins2(), m_histogram->minVar1(), m_histogram->maxVar1(), m_histogram->minVar2(),  m_histogram->maxVar2(), m_histogram->shift1(), m_histogram->shift2(), m_histogram->bin_type1(), m_histogram->bin_type2());

  histogram->put(m_data->var(m_Var1), m_data->var(m_Var2), m_data->var(catalogue::Var::_Weight_));

  m_histogram = histogram;

  vector<double> bins1 = m_histogram->bins1();
  vector<double> edges1 = m_histogram->edges1();

  vector<double> bins2 = m_histogram->bins2();
  vector<double> edges2 = m_histogram->edges2();

  vector<double> hist(m_histogram->nbins1()*m_histogram->nbins2());
  vector<double> error(m_histogram->nbins1()*m_histogram->nbins2());

  vector<vector<double>> extra_info(5);

  for(size_t i=0; i<m_histogram->nbins1(); i++)
    for(size_t j=0; j<m_histogram->nbins2(); j++){
      hist[i*m_histogram->nbins2()+j] = m_histogram->operator()(i, j, m_HistogramType, m_fact);
      error[i*m_histogram->nbins2()+j] = m_histogram->poisson_error(i, j, m_HistogramType, m_fact);
      extra_info[0].push_back(edges1[i]);
      extra_info[1].push_back(edges1[i+1]);
      extra_info[2].push_back(edges2[j]);
      extra_info[3].push_back(edges2[j+1]);
      extra_info[4].push_back(m_histogram->weight(i, j));
    }

  auto dataset = make_shared<data::Data2D_extra> (data::Data2D_extra(bins1, bins2, hist, error, extra_info));

  return dataset;
}


// ============================================================================


shared_ptr<data::Data> cbl::measure::numbercounts::NumberCounts2D::m_measureJackknife (const string dir_output_resample)
{
  auto histogram =  make_shared<glob::Histogram2D> (glob::Histogram2D ());
  histogram->set(m_histogram->nbins1(), m_histogram->nbins2(), m_histogram->minVar1(), m_histogram->maxVar1(), m_histogram->minVar2(),  m_histogram->maxVar2(), m_histogram->shift1(), m_histogram->shift2(), m_histogram->bin_type1(), m_histogram->bin_type2());

  const int nRegions = m_data->nRegions();
  vector<shared_ptr<glob::Histogram>> histo_JK(nRegions);

  for (int i=0; i<nRegions; i++) {
    histo_JK[i] = make_shared<glob::Histogram2D>(glob::Histogram2D());
    histo_JK[i]->set(m_histogram->nbins1(), m_histogram->nbins2(), m_histogram->minVar1(), m_histogram->maxVar1(), m_histogram->minVar2(),  m_histogram->maxVar2(), m_histogram->shift1(), m_histogram->shift2(), m_histogram->bin_type1(), m_histogram->bin_type2());
  }

  vector<double> regions = m_data->var(catalogue::Var::_Region_);
  vector<double> var1 = m_data->var(m_Var1);
  vector<double> var2 = m_data->var(m_Var2);
  vector<double> weight = m_data->var(catalogue::Var::_Weight_);

  for (size_t i=0; i<m_data->nObjects(); i++) {
    vector<int> bin = histogram->digitize(var1[i], var2[i]);
    histogram->put(bin[0], bin[1], weight[i]);
    histo_JK[regions[i]]->put (bin[0], bin[1], weight[i]);
  }

  m_histogram = histogram;

  vector<double> bins1 = m_histogram->bins1();
  vector<double> edges1 = m_histogram->edges1();

  vector<double> bins2 = m_histogram->bins2();
  vector<double> edges2 = m_histogram->edges2();

  vector<double> hist(m_histogram->nbins1()*m_histogram->nbins2());
  vector<double> error(m_histogram->nbins1()*m_histogram->nbins2());

  vector<vector<double>> extra_info(5);

  for(size_t i=0; i<m_histogram->nbins1(); i++)
    for(size_t j=0; j<m_histogram->nbins2(); j++){
      hist[i*m_histogram->nbins2()+j] = m_histogram->operator()(i, j, m_HistogramType, m_fact);
      error[i*m_histogram->nbins2()+j] = m_histogram->poisson_error(i, j, m_HistogramType, m_fact);
      extra_info[0].push_back(edges1[i]);
      extra_info[1].push_back(edges1[i+1]);
      extra_info[2].push_back(edges2[i]);
      extra_info[3].push_back(edges2[i+1]);
      extra_info[4].push_back(m_histogram->weight(i, j));
    }

  compute_covariance(histo_JK, true);
  auto dataset = make_shared<data::Data2D_extra> (data::Data2D_extra(bins1, bins2, hist, m_dataset->covariance(), extra_info));

  // write jackknife outputs
  if (dir_output_resample!=par::defaultString)
    for (int i=0; i<nRegions; i++) 
      histo_JK[i]->write(dir_output_resample, "Jackknife_"+conv(i+1, par::fINT)+".dat", m_HistogramType, m_fact);

  return dataset;
}

// ============================================================================


shared_ptr<data::Data> cbl::measure::numbercounts::NumberCounts2D::m_measureBootstrap (const string dir_output_resample, const int nResamplings, const int seed) 
{
  auto histogram =  make_shared<glob::Histogram2D> (glob::Histogram2D ());
  histogram->set(m_histogram->nbins1(), m_histogram->nbins2(), m_histogram->minVar1(), m_histogram->maxVar1(), m_histogram->minVar2(),  m_histogram->maxVar2(), m_histogram->shift1(), m_histogram->shift2(), m_histogram->bin_type1(), m_histogram->bin_type2());

  const int nRegions = m_data->nRegions();
  vector<shared_ptr<glob::Histogram>> histo_BS(nResamplings);

  random::UniformRandomNumbers_Int ran(0., nRegions-1, seed);
  vector<vector<int>> region_weights (nResamplings, vector<int> (nRegions, 0));

  for (int i=0; i<nResamplings; i++) {
    histo_BS[i] = make_shared<glob::Histogram2D>(glob::Histogram2D());
    histo_BS[i]->set(m_histogram->nbins1(), m_histogram->nbins2(), m_histogram->minVar1(), m_histogram->maxVar1(), m_histogram->minVar2(),  m_histogram->maxVar2(), m_histogram->shift1(), m_histogram->shift2(), m_histogram->bin_type1(), m_histogram->bin_type2());
    for (int n=0; n<nRegions; n++)
      region_weights[i][ran()] ++;
  }

  vector<double> regions = m_data->var(catalogue::Var::_Region_);
  vector<double> var1 = m_data->var(m_Var1);
  vector<double> var2 = m_data->var(m_Var2);
  vector<double> weight = m_data->var(catalogue::Var::_Weight_);

  for (size_t i=0; i<m_data->nObjects(); i++) {
    vector<int> bin = m_histogram->digitize(var1[i], var2[i]);
    m_histogram->put(bin[0], bin[1], weight[i]);
    for (int j=0; j<nResamplings; j++) {
      histo_BS[j]->put (bin[0], bin[1], weight[i]*region_weights[j][regions[i]]);
    }
  }

  m_histogram = histogram;

  vector<double> bins1 = m_histogram->bins1();
  vector<double> edges1 = m_histogram->edges1();

  vector<double> bins2 = m_histogram->bins2();
  vector<double> edges2 = m_histogram->edges2();

  vector<double> hist(m_histogram->nbins1()*m_histogram->nbins2());
  vector<double> error(m_histogram->nbins1()*m_histogram->nbins2());

  vector<vector<double>> extra_info(5);

  for(size_t i=0; i<m_histogram->nbins1(); i++)
    for(size_t j=0; j<m_histogram->nbins2(); j++){
      hist[i*m_histogram->nbins2()+j] = m_histogram->operator()(i, j, m_HistogramType, m_fact);
      error[i*m_histogram->nbins2()+j] = m_histogram->poisson_error(i, j, m_HistogramType, m_fact);
      extra_info[0].push_back(edges1[i]);
      extra_info[1].push_back(edges1[i+1]);
      extra_info[2].push_back(edges2[i]);
      extra_info[3].push_back(edges2[i+1]);
      extra_info[4].push_back(m_histogram->weight(i, j));
    }


  compute_covariance(histo_BS, true);
  auto dataset = make_shared<data::Data2D_extra> (data::Data2D_extra(bins1, bins2, hist, m_dataset->covariance(), extra_info));

  // write bootstrap outputs
  if (dir_output_resample!=par::defaultString)
    for (int i=0; i<nResamplings; i++)
      histo_BS[i]->write(dir_output_resample, "Bootstraps_"+conv(i+1, par::fINT)+".dat", m_HistogramType, m_fact);

  return dataset;
}


// ============================================================================


void cbl::measure::numbercounts::NumberCounts2D::measure (const ErrorType errorType, const string dir_output_resample, const int nResamplings, const int seed)
{
  switch (errorType) {
    case (ErrorType::_Poisson_) :
      m_dataset = m_measurePoisson();
      break;
    case (ErrorType::_Jackknife_) :
      m_dataset = m_measureJackknife(dir_output_resample);
      break;
    case (ErrorType::_Bootstrap_) :
      m_dataset = m_measureBootstrap(dir_output_resample, nResamplings, seed);
      break;

    default:
      ErrorCBL("Error in measure() of NumberCounts2D.cpp, unknown type of error");
  }

}

// ============================================================================


void cbl::measure::numbercounts::NumberCounts2D::compute_covariance (const vector<shared_ptr<glob::Histogram>> histo, const bool JK)
{
  vector<vector<double>> measures(histo.size(), vector<double>(histo[0]->nbins1()*histo[0]->nbins2(),0));
  vector<vector<double>> cov;

  for(size_t i=0; i<histo.size(); i++)
    for(size_t j=0; j<histo[i]->nbins1(); j++)
      for(size_t k=0; k<histo[i]->nbins2(); k++)
	measures[i][j*histo[0]->nbins2()+k] = histo[i]->operator()(j, k, m_HistogramType, m_fact);

  covariance_matrix (measures, cov, JK);
  m_dataset->set_covariance(cov);
}


// ============================================================================


void cbl::measure::numbercounts::NumberCounts2D::write (const string dir, const string file, const int rank) const 
{
  string mkdir = "mkdir -p "+dir;
  if(system(mkdir.c_str())) {}

  string header = "# bin1_center bin2_center counts error r1_1 r2_1 r1_2 r2_2";
  m_dataset->write(dir, file, header, false, 8, rank);
}


// ============================================================================


void cbl::measure::numbercounts::NumberCounts2D::write_covariance (const string dir, const string file) const 
{ 
  string mkdir = "mkdir -p "+dir;
  if(system(mkdir.c_str())) {}
  m_dataset->write_covariance(dir, file, 8);
}
