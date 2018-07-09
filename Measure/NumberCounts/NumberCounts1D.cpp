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
 *  @file NumberCounts/NumberCounts1D.cpp
 *
 *  @brief Methods of the class NumberCounts1D 
 *
 *  This file contains the implementation of the methods of the class
 *  NumberCounts1D, used to measure number counts
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "NumberCounts1D.h"
#include "Data1D.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace measure::numbercounts;


// ============================================================================


cbl::measure::numbercounts::NumberCounts1D::NumberCounts1D(const catalogue::Var var, const BinType bin_type, const catalogue::Catalogue data, const size_t nbins, const double minVar, const double maxVar, const double shift, const glob::HistogramType hist_type, const double fact)
{
  m_Var = var;

  m_HistogramType = hist_type;
  m_fact = fact;

  set_data(data);

  m_histogram = make_shared<glob::Histogram1D> (glob::Histogram1D ()); 

  double _minVar = (minVar>par::defaultDouble) ? minVar : m_data->Min(m_Var)*0.999;
  double _maxVar = (maxVar>par::defaultDouble) ? maxVar : m_data->Max(m_Var)*1.001;

  m_histogram->set(nbins, _minVar, _maxVar, shift, bin_type);
}


// ============================================================================


shared_ptr<data::Data> cbl::measure::numbercounts::NumberCounts1D::m_measurePoisson ()
{
  auto histogram =  make_shared<glob::Histogram1D> (glob::Histogram1D ());
  histogram->set(m_histogram->nbins(), m_histogram->minVar(), m_histogram->maxVar(), m_histogram->shift(), m_histogram->bin_type());

  histogram->put(m_data->var(m_Var), m_data->var(catalogue::Var::_Weight_));

  m_histogram = histogram;

  vector<double> bins = m_histogram->bins();
  vector<double> edges = m_histogram->edges();
  vector<double> hist = m_histogram->operator()(m_HistogramType, m_fact);
  vector<double> weights = m_histogram->weights();
  vector<double> error = m_histogram->poisson_error(m_HistogramType, m_fact);

  vector<vector<double>> extra_info (3);
  extra_info[0] = slice(edges, 0, m_histogram->nbins());
  extra_info[1] = slice(edges, 1, m_histogram->nbins()+1);
  extra_info[2] = weights;

  auto dataset = make_shared<data::Data1D_extra> (data::Data1D_extra(bins, hist, error, extra_info));

  return dataset;
}


// ============================================================================


shared_ptr<data::Data> cbl::measure::numbercounts::NumberCounts1D::m_measureJackknife (const string dir_output_resample)
{
  auto histogram =  make_shared<glob::Histogram1D> (glob::Histogram1D ());
  histogram->set(m_histogram->nbins(), m_histogram->minVar(), m_histogram->maxVar(), m_histogram->shift(), m_histogram->bin_type());

  const int nRegions = m_data->nRegions();
  vector<shared_ptr<glob::Histogram>> histo_JK(nRegions);

  for (int i=0; i<nRegions; i++){
    histo_JK[i] = make_shared<glob::Histogram1D>(glob::Histogram1D());
    histo_JK[i]->set(m_histogram->nbins(), m_histogram->minVar(), m_histogram->maxVar(), m_histogram->shift(), m_histogram->bin_type());
  }

  vector<double> regions = m_data->var(catalogue::Var::_Region_);
  vector<double> var = m_data->var(m_Var);
  vector<double> weight = m_data->var(catalogue::Var::_Weight_);

  for (size_t i=0; i<m_data->nObjects(); i++) {
    int bin = histogram->digitize(var[i]);
    histogram->put(bin, weight[i]);
    histo_JK[regions[i]]->put (bin, weight[i]);
  }
  
  m_histogram = histogram;

  vector<double> bins = m_histogram->bins();
  vector<double> edges = m_histogram->edges();
  vector<double> hist = m_histogram->operator()(m_HistogramType, m_fact);
  vector<double> weights = m_histogram->weights();
  vector<double> error = m_histogram->poisson_error(m_HistogramType, m_fact);

  vector<vector<double>> extra_info (3);
  extra_info[0] = slice(edges, 0, m_histogram->nbins());
  extra_info[1] = slice(edges, 1, m_histogram->nbins()+1);
  extra_info[2] = weights;

  vector<vector<double>> measures(histo_JK.size(), vector<double>(histo_JK[0]->nbins()));
  vector<vector<double>> covariance;

  for (size_t i=0; i<histo_JK.size(); i++)
    for (size_t j=0; j<histo_JK[i]->nbins(); j++)
      measures[i][j] = histo_JK[i]->operator()(j, m_HistogramType, m_fact);

  covariance_matrix (measures, covariance, true);

  auto dataset = make_shared<data::Data1D_extra> (data::Data1D_extra(bins, hist, covariance, extra_info));

  //Write jackknife outputs

  if (dir_output_resample!=par::defaultString)
    for (int i=0; i<nRegions; i++) 
      histo_JK[i]->write(dir_output_resample, "Jackknife_"+conv(i+1, par::fINT)+".dat", m_HistogramType, m_fact);

  coutCBL << "Done!"<<endl;

  return dataset;
}

// ============================================================================


shared_ptr<data::Data> cbl::measure::numbercounts::NumberCounts1D::m_measureBootstrap (const string dir_output_resample, const int nResamplings, const int seed) 
{
  auto histogram =  make_shared<glob::Histogram1D> (glob::Histogram1D ());
  histogram->set(m_histogram->nbins(), m_histogram->minVar(), m_histogram->maxVar(), m_histogram->shift(), m_histogram->bin_type());

  const int nRegions = m_data->nRegions();
  vector<shared_ptr<glob::Histogram>> histo_BS(nResamplings);

  random::UniformRandomNumbers_Int ran(0., nRegions-1, seed);
  vector<vector<int>> region_weights (nResamplings, vector<int> (nRegions, 0));

  for (int i=0; i<nResamplings; i++) {
    histo_BS[i] = make_shared<glob::Histogram1D>(glob::Histogram1D());
    histo_BS[i]->set(m_histogram->nbins(), m_histogram->minVar(), m_histogram->maxVar(), m_histogram->shift(), m_histogram->bin_type());
    for (int n=0; n<nRegions; n++)
      region_weights[i][ran()] ++;
  }

  vector<double> regions = m_data->var(catalogue::Var::_Region_);
  vector<double> var = m_data->var(m_Var);
  vector<double> weight = m_data->var(catalogue::Var::_Weight_);

  for (size_t i=0; i<m_data->nObjects(); i++) {
    int bin = histogram->digitize(var[i]);
    histogram->put(bin, weight[i]);
    for (int j=0; j<nResamplings; j++) {
      histo_BS[j]->put (bin, weight[i]*region_weights[j][regions[i]]);
    }
  }

  m_histogram = histogram;

  vector<double> bins = m_histogram->bins();
  vector<double> edges = m_histogram->edges();
  vector<double> hist = m_histogram->operator()(m_HistogramType, m_fact);
  vector<double> weights = m_histogram->weights();
  vector<double> error = m_histogram->poisson_error(m_HistogramType, m_fact);

  vector<vector<double>> extra_info (3);
  extra_info[0] = slice(edges, 0, m_histogram->nbins()-1);
  extra_info[1] = slice(edges, 1, m_histogram->nbins());
  extra_info[2] = weights;

  vector<vector<double>> measures(histo_BS.size(), vector<double>(histo_BS[0]->nbins()));
  vector<vector<double>> covariance;

  for (size_t i=0; i<histo_BS.size(); i++)
    for (size_t j=0; j<histo_BS[i]->nbins(); j++)
      measures[i][j] = histo_BS[i]->operator()(j, m_HistogramType, m_fact);

  covariance_matrix (measures, covariance, false);

  auto dataset = make_shared<data::Data1D_extra> (data::Data1D_extra(bins, hist, covariance, extra_info));

  //Write bootstrap outputs

  if (dir_output_resample!=par::defaultString)
    for (int i=0; i<nResamplings; i++)
      histo_BS[i]->write(dir_output_resample, "Bootstraps_"+conv(i+1, par::fINT)+".dat", m_HistogramType, m_fact);

  return dataset;
}


// ============================================================================


void cbl::measure::numbercounts::NumberCounts1D::measure (const ErrorType errorType, const string dir_output_resample, const int nResamplings, const int seed)
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
      ErrorCBL("Error in measure() of NumberCounts1D.cpp, unknown type of error");
  }
}


// ============================================================================


void cbl::measure::numbercounts::NumberCounts1D::compute_covariance (const vector<shared_ptr<glob::Histogram>> histo, const bool JK)
{
  vector<vector<double>> measures(histo.size(), vector<double>(histo[0]->nbins()));
  vector<vector<double>> cov;

  for (size_t i=0; i<histo.size(); i++)
    for (size_t j=0; j<histo[i]->nbins(); j++)
      measures[i][j] = histo[i]->operator()(j, m_HistogramType, m_fact);

  covariance_matrix (measures, cov, JK);
  m_dataset->set_covariance(cov);
}


// ============================================================================


void cbl::measure::numbercounts::NumberCounts1D::write (const string dir, const string file, const int rank) const 
{
  string mkdir = "mkdir -p "+dir;
  if (system(mkdir.c_str())) {}

  string header = "# bin_center counts error  r1   r2";
  m_dataset->write(dir, file, header, 8, rank);
}


// ============================================================================


void cbl::measure::numbercounts::NumberCounts1D::write_covariance (const string dir, const string file) const 
{
  string mkdir = "mkdir -p "+dir;
  if (system(mkdir.c_str())) {}

  m_dataset->write_covariance(dir, file, 8);
}
