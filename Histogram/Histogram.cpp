/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo *
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
 *  @file Histogram/Histogram.cpp
 *
 *  @brief Methods of the class Histogram
 *
 *  This file contains the implementation of the methods of the class
 *  Histogram1D and Histogram_2D
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "Histogram.h"

using namespace std;

using namespace cbl;
using namespace glob;


// ============================================================================


cbl::glob::Histogram1D::Histogram1D (const std::vector<double> var, const std::vector<double> weight, const size_t nbins, const double minVar, const double maxVar, const double shift, const BinType bin_type)
{
  double _minVar = (minVar>par::defaultDouble) ? minVar : Min(var)*0.9999;
  double _maxVar = (maxVar>par::defaultDouble) ? maxVar : Max(var)*1.0001;

  set(nbins, _minVar, _maxVar, shift, bin_type);
  put(var, weight);
}


// ============================================================================


void cbl::glob::Histogram1D::set (const size_t nbins, const double minVar, const double maxVar, const double shift, const BinType bin_type)
{
  if (shift>1 || shift<0) 
    ErrorCBL("shift must be 0<shift<1!", "set", "Histogram.cpp");

  m_nbins = nbins;
  m_minVar = minVar;
  m_maxVar = maxVar;
  m_shift = shift;

  m_bins.resize(m_nbins, 0);
  m_edges.resize(m_nbins+1, 0);
  m_counts.resize(m_nbins, 0.);
  m_var.resize(m_nbins, 0.);
  m_var_err.resize(m_nbins, 0.);

  shared_ptr<gsl_histogram> histo1(gsl_histogram_alloc(m_nbins), gsl_histogram_free);
  shared_ptr<gsl_histogram> histo2(gsl_histogram_alloc(m_nbins), gsl_histogram_free);

  m_binType = bin_type;

  if (m_binType == BinType::_linear_) {
    m_binSize = (m_maxVar-m_minVar)/m_nbins;
    gsl_histogram_set_ranges_uniform(histo1.get(), m_minVar, m_maxVar);
    gsl_histogram_set_ranges_uniform(histo2.get(), m_minVar, m_maxVar);

    m_edges[0] = histo1.get()->range[0];
    for (size_t i=0; i<m_nbins; i++){
      m_edges[i+1] = histo1.get()->range[i+1];
      m_bins[i] = m_edges[i]+shift*m_binSize;
    }

  }
  else if (m_binType == BinType::_logarithmic_) {

    m_binSize = (log10(m_maxVar)-log10(m_minVar))/nbins;

    m_edges[0] = m_minVar;
    for (size_t i=0; i<m_nbins; i++){
      m_edges[i+1] = pow(10., log10(m_minVar)+(i+1)*m_binSize);
      m_bins[i] = pow(10., log10(m_edges[i])+m_shift*m_binSize);
    }

    gsl_histogram_set_ranges(histo1.get(), m_edges.data(), m_nbins+1);
    gsl_histogram_set_ranges(histo2.get(), m_edges.data(), m_nbins+1);
  }

  m_histo = histo1;
  m_histo_error =  histo2;
}


// ============================================================================


int cbl::glob::Histogram1D::digitize (const double var)
{
  gsl_set_error_handler_off();

  size_t i;
  int res = gsl_histogram_find(m_histo.get(), var, &i);

  if (res == GSL_EDOM)
    return -1;
  else
    return i;
}


// ============================================================================


std::vector<int> cbl::glob::Histogram1D::digitize (const std::vector<double> var)
{
  vector<int> bin(var.size());
  for (size_t i=0; i<var.size(); i++)
    bin[i] = digitize(var[i]);
  return bin;
}


// ============================================================================


void cbl::glob::Histogram1D::put (const double var, const double weight)
{
  const int bin = digitize(var);
  if (bin>-1) 
    put(bin, weight, var);
}


// ============================================================================


void cbl::glob::Histogram1D::put (const std::vector<double> var, const std::vector<double> weight)
{
  for (size_t i=0; i<var.size(); i++)
    put(var[i], weight[i]);
}


// ============================================================================


void cbl::glob::Histogram1D::put (const int bin, const double weight, const double var)
{
  m_histo.get()->bin[bin] += weight;
  m_histo_error.get()->bin[bin] += weight*weight;
  m_counts[bin] += 1;

  double Wn = m_histo.get()->bin[bin];
  double mean_prev = m_var[bin];
  double var_prev = m_var_err[bin];
  m_var[bin] = (Wn==0) ? 0 : mean_prev+weight/Wn*(var-mean_prev);
  m_var_err[bin] = (Wn==0) ? 0 : var_prev+weight*(var-mean_prev)*(var-m_var[bin]);

  if (m_var[bin]!= m_var[bin])
    cout << Wn << " " << mean_prev << " " << var_prev << " " << m_var[bin] << " " << m_var_err[bin] << " " << var << " " << weight << endl;
}


// ============================================================================


void cbl::glob::Histogram1D::put (const std::vector<int> bins, const std::vector<double> weight, const std::vector<double> var)
{
  for (size_t i=0; i<bins.size(); i++)
    put(bins[i], weight[i], var[i]);
}


// ============================================================================


double cbl::glob::Histogram1D::normalization ( const int i, const HistogramType hist_type, const double fact) const
{
  double norm;

  switch (hist_type) {

    case (HistogramType::_dn_dV_):
      norm = fact*(m_edges[i+1]-m_edges[i]) ;
      break;

    case (HistogramType::_dn_dlogV_):
      norm = fact*(log10(m_edges[i+1])-log10(m_edges[i])) ;
      break;

    case (HistogramType::_dn_dlnV_):
      norm = fact*(log(m_edges[i+1])-log(m_edges[i])) ;
      break;

    case (HistogramType::_N_V_):
      norm = 1.;
      break;

    case (HistogramType::_n_V_):
      norm = fact;
      break;

    default:
      ErrorCBL("no such a variable in the list!", "normalization", "Histogram1D.cpp");
  }

  return norm;
}


// ============================================================================


double cbl::glob::Histogram1D::operator() ( const int i, const HistogramType hist_type, const double fact) const
{
  return m_histo.get()->bin[i]/normalization(i, hist_type, fact);

}


// ============================================================================


std::vector<double> cbl::glob::Histogram1D::operator() ( const HistogramType hist_type, const double fact) const
{
  vector<double> hist(m_nbins, 0);

  for (size_t i=0; i<m_nbins; i++)
    hist[i] = m_histo.get()->bin[i]/normalization(i, hist_type, fact);

  return hist;
}


// ============================================================================


double cbl::glob::Histogram1D::error ( const int i, const HistogramType hist_type, const double fact) const
{
  return sqrt(m_histo_error.get()->bin[i])/normalization(i, hist_type, fact);
}


// ============================================================================


std::vector<double> cbl::glob::Histogram1D::error ( const HistogramType hist_type, const double fact) const
{
  vector<double> hist(m_nbins, 0);

  for (size_t i=0; i<m_nbins; i++)
    hist[i] = sqrt(m_histo_error.get()->bin[i])/normalization(i, hist_type, fact);

  return hist;
}

// ============================================================================


vector<double> cbl::glob::Histogram1D::error_bins () const
{
  std::vector<double> vv = m_var_err;
  for (size_t i=0; i<m_nbins; i++)
    vv[i] = sqrt(vv[i]/m_histo.get()->bin[i]);
  return vv;
}


// ============================================================================


void cbl::glob::Histogram1D::write (const string dir, const string file, const HistogramType hist_type, const double fact) const
{
  string mkdir = "mkdir -p "+dir;
  if (system(mkdir.c_str())) {}
  
  string output_file = dir+file;cout <<output_file<<endl;
  ofstream fout (output_file.c_str());

  fout << "# bin_center hist hist_err weighted_counts weighted_counts_err counts counts_err  lower_edge upper_edge histogram_normalization var var_err" << endl;

  for (size_t i=0; i<m_nbins; i++) {
    fout << bin(i) <<  " " << this->operator()(i, hist_type, fact) << " " << this->error(i, hist_type, fact);
    fout << " " << m_histo.get()->bin[i] << " " << sqrt(m_histo_error.get()->bin[i]);
    fout << " " << m_counts[i] << " " << sqrt(m_counts[i]);
    fout << " " << edge(i) << " " << edge(i+1) << " " << normalization(i, hist_type, fact);
    fout << " " << m_var[i] << " " << sqrt(m_var_err[i]/m_histo.get()->bin[i]) << endl;
  }

  fout.clear(); fout.close();
}


// ============================================================================


cbl::glob::Histogram2D::Histogram2D (const std::vector<double> var1, const std::vector<double> var2, const std::vector<double> weight, const size_t nbins1, const size_t nbins2, const double minVar1, const double maxVar1, const double minVar2, const double maxVar2, const double shift1, const double shift2, const BinType bin_type1, const BinType bin_type2) 
{
  double _minVar1 = (minVar1>par::defaultDouble) ? minVar1 : Min(var1)*0.9999;
  double _maxVar1 = (maxVar1>par::defaultDouble) ? maxVar1 : Max(var1)*1.0001;

  double _minVar2 = (minVar2>par::defaultDouble) ? minVar2 : Min(var2)*0.9999;
  double _maxVar2 = (maxVar2>par::defaultDouble) ? maxVar2 : Max(var2)*1.0001;

  set(nbins1, nbins2, _minVar1, _maxVar1, _minVar2, _maxVar2, shift1, shift2, bin_type1, bin_type2);
  put(var1, var2, weight);
}


// ============================================================================


void cbl::glob::Histogram2D::set (const size_t nbins1, const size_t nbins2, const double minVar1, const double maxVar1, const double minVar2, const double maxVar2, const double shift1, const double shift2, const BinType bin_type1, const BinType bin_type2)
{
  if (shift1>1 || shift1<0 || shift2>1 || shift2<0) 
    ErrorCBL("shift must be in the range (0,1) !", "set", "Histrogram.cpp");

  m_nbins1 = nbins1;
  m_minVar1 = minVar1;
  m_maxVar1 = maxVar1;
  m_shift1 = shift1;
  m_bins1.resize(m_nbins1, 0);
  m_edges1.resize(m_nbins1+1, 0);
  m_binType1 = bin_type1;

  m_nbins2 = nbins2;
  m_minVar2 = minVar2;
  m_maxVar2 = maxVar2;
  m_shift2 = shift2;
  m_bins2.resize(m_nbins2, 0);
  m_edges2.resize(m_nbins2+1, 0);
  m_binType2 = bin_type2;

  m_counts.resize(m_nbins1, vector<int>(m_nbins2, 0));
  m_var1.resize(m_nbins1, vector<double>(m_nbins2, 0));
  m_var1_err.resize(m_nbins1, vector<double>(m_nbins2, 0));
  m_var2.resize(m_nbins1, vector<double>(m_nbins2, 0));
  m_var2_err.resize(m_nbins1, vector<double>(m_nbins2, 0));

  shared_ptr<gsl_histogram2d> histo(gsl_histogram2d_alloc(m_nbins1, m_nbins2), gsl_histogram2d_free);
  shared_ptr<gsl_histogram2d> histo2(gsl_histogram2d_alloc(m_nbins1, m_nbins2), gsl_histogram2d_free);

  if (m_binType1 == BinType::_linear_ && m_binType2 == BinType::_linear_) {
    m_binSize1 = (m_maxVar1-m_minVar1)/m_nbins1;
    m_binSize2 = (m_maxVar2-m_minVar2)/m_nbins2;

    gsl_histogram2d_set_ranges_uniform(histo.get(), m_minVar1, m_maxVar1, m_minVar2, m_minVar2);
    gsl_histogram2d_set_ranges_uniform(histo2.get(), m_minVar1, m_maxVar1, m_minVar2, m_minVar2);

    m_edges1[0] = histo.get()->xrange[0];
    for (size_t i=0; i<m_nbins1; i++) {
      m_edges1[i+1] = histo.get()->xrange[i+1];
      m_bins1[i] = m_edges1[i]+m_shift1*m_binSize1;
    }

    m_edges2[0] = histo.get()->yrange[0];
    for (size_t i=0; i<m_nbins2; i++) {
      m_edges2[i+1] = histo.get()->yrange[i+1];
      m_bins2[i] = m_edges2[i]+m_shift2*m_binSize2;
    }

  }
  else {

    // set ranges for variable 1
    if ( m_binType1==BinType::_linear_) {
      m_binSize1 = (m_maxVar1-m_minVar1)/m_nbins1;
      m_edges1[0] = m_minVar1;
      for (size_t i=0; i<m_nbins1; i++) {
	m_edges1[i+1] = m_edges1[i]+m_binSize1;
	m_bins1[i] = m_edges1[i]+m_shift1*m_binSize1;
      }
    }
    else if (m_binType1==BinType::_logarithmic_) {
      m_binSize1 = (log10(m_maxVar1)-log10(m_minVar1))/m_nbins1;
      m_edges1[0] = m_minVar1;
      for (size_t i=0; i<m_nbins1; i++) {
	m_edges1[i+1] = pow(10., log10(m_edges1[i])+m_binSize1);
	m_bins1[i] = pow(10., log10(m_edges1[i])+m_shift1*m_binSize1);
      }
    }
    else 
      ErrorCBL("no such bin_type!", "set", "Histogram2D.cpp");

    // set ranges for variable 2
    if ( m_binType2 == BinType::_linear_) {
      m_binSize2 = (m_maxVar2-m_minVar2)/m_nbins2;
      m_edges2[0] = m_minVar2;
      for (size_t i=0; i<m_nbins2; i++){
	m_edges2[i+1] = m_edges2[i]+m_binSize2;
	m_bins2[i] = m_edges2[i]+m_shift2*m_binSize2;
      }
    }
    else if (m_binType2 == BinType::_logarithmic_) {
      m_binSize2 = (log10(m_maxVar2)-log10(m_minVar2))/m_nbins2;
      m_edges2[0] = m_minVar2;
      for (size_t i=0; i<m_nbins2; i++) {
	m_edges2[i+1] = pow(10., log10(m_edges2[i])+m_binSize2);
	m_bins2[i] = pow(10., log10(m_edges2[i])+m_shift2*m_binSize2);
      }
    }
    else 
      ErrorCBL("the input bin_type is not allowed!", "set", "Histrogram.cpp");
    
    gsl_histogram2d_set_ranges(histo.get(), m_edges1.data(), m_nbins1+1, m_edges2.data(), m_nbins2+1);
    gsl_histogram2d_set_ranges(histo2.get(), m_edges1.data(), m_nbins1+1, m_edges2.data(), m_nbins2+1);
  }

  m_histo = histo;
  m_histo_error = histo2;

}


// ============================================================================


std::vector<int> cbl::glob::Histogram2D::digitize (const double var1, const double var2)
{
  size_t i, j;
  int result = gsl_histogram2d_find(m_histo.get(), var1, var2, &i, &j);

  if (result==GSL_SUCCESS)
    return {(int)i, (int)j};
  return {-1, -1};
}


// ============================================================================


std::vector<std::vector<int>> cbl::glob::Histogram2D::digitize (const std::vector<double> var1, const std::vector<double> var2)
{
  vector<vector<int>> bins (var1.size(), vector<int>(2, 0));

  for (size_t i=0; i<var1.size(); i++)
    bins[i] = digitize(var1[i], var2[i]);

  return bins;
}


// ============================================================================


void cbl::glob::Histogram2D::put (const double var1, const double var2, const double weight)
{
  const vector<int> bins = digitize(var1, var2);
  if (bins[0]>-1 && bins[1]>-1){
    put(bins[0], bins[1], weight);
  }
}


// ============================================================================


void cbl::glob::Histogram2D::put (const std::vector<double> var1, const std::vector<double> var2, const std::vector<double> weight)
{
  for (size_t i=0; i<var1.size(); i++)
    put(var1[i], var2[i], weight[i]);
}


// ============================================================================


void cbl::glob::Histogram2D::put (const int i, const int j, const double weight, const double var1, const double var2)
{
  m_histo.get()->bin[i*m_nbins2+j] += weight;
  m_histo_error.get()->bin[i*m_nbins2+j] += weight*weight;
  m_counts[i][j] += 1;

  double Wn = m_histo.get()->bin[i*m_nbins2+j];
  double mean_prev = m_var1[i][j];
  double var_prev = m_var1_err[i][j];

  m_var1[i][j] = (Wn==0) ? 0 : mean_prev+weight/Wn*(var1-mean_prev);
  m_var1_err[i][j] = (Wn==0) ? 0 : var_prev+weight*(var1-mean_prev)*(var1-m_var1[i][j]);

  mean_prev = m_var2[i][j];
  var_prev = m_var2_err[i][j];

  m_var2[i][j] = (Wn==0) ? 0 :  mean_prev+weight/Wn*(var2-mean_prev);
  m_var2_err[i][j] = (Wn==0) ? 0 :  var_prev+weight*(var2-mean_prev)*(var2-m_var2[i][j]);
}


// ============================================================================


void cbl::glob::Histogram2D::put (const std::vector<std::vector<int>> bins, const std::vector<double> weight, const std::vector<double> var1, const std::vector<double> var2)
{
  for (size_t i=0; i<weight.size(); i++)
    put(bins[i][0], bins[i][1], weight[i], var1[i], var2[i]);
}


// ============================================================================


double cbl::glob::Histogram2D::normalization ( const int i, const int j, const HistogramType hist_type, const double fact) const
{
  double norm;

  switch (hist_type) {

    case (HistogramType::_dn_dV_):
      norm = fact*(m_edges1[i+1]-m_edges1[i])*(m_edges2[j+1]-m_edges2[j]) ;
      break;

    case (HistogramType::_dn_dlogV_):
      norm = fact*(log10(m_edges1[i+1])-log10(m_edges1[i]))*(log10(m_edges2[j+1])-log10(m_edges2[j])) ;
      break;

    case (HistogramType::_dn_dlnV_):
      norm = fact*(log(m_edges1[i+1])-log(m_edges1[i]))*(log(m_edges2[j+1])-log(m_edges2[j])) ;
      break;

    case (HistogramType::_N_V_):
      norm = 1.;
      break;

    case (HistogramType::_n_V_):
      norm = fact;
      break;

    default:
      ErrorCBL("no such a variable in the list!", "normalization", "Histogram.cpp");
  }

  return norm;
}


// ============================================================================


double cbl::glob::Histogram2D::operator() (const int i, const int j, const HistogramType hist_type, const double fact) const
{
  return m_histo.get()->bin[i*m_nbins2+j]/normalization(i, j, hist_type, fact);
}

// ============================================================================


vector<vector<double>> cbl::glob::Histogram2D::error_bins1 () const
{
  std::vector<vector<double>> vv = m_var1_err;
  for (size_t i=0; i<m_nbins1; i++)
    for (size_t j=0; j<m_nbins2; j++)
      vv[i][j] = sqrt(vv[i][j]/m_histo.get()->bin[i*m_nbins2+j]);
  return vv;
}


// ============================================================================


vector<vector<double>> cbl::glob::Histogram2D::error_bins2 () const
{
  std::vector<vector<double>> vv = m_var2_err;
  for (size_t i=0; i<m_nbins1; i++)
    for (size_t j=0; j<m_nbins2; j++)
      vv[i][j] = sqrt(vv[i][j]/m_histo.get()->bin[i*m_nbins2+j]);
  return vv;
}


// ============================================================================


double cbl::glob::Histogram2D::error (const int i, const int j, const HistogramType hist_type, const double fact) const
{
  return sqrt(m_histo_error.get()->bin[i*m_nbins2+j])/normalization(i, j, hist_type, fact);
}


// ============================================================================


void cbl::glob::Histogram2D::write (const string dir, const string file, const HistogramType hist_type, const double fact) const
{
  string mkdir = "mkdir -p "+dir;
  if (system(mkdir.c_str())) {}

  string output_file = dir+file;
  ofstream fout (output_file.c_str());

  fout << "# bin1_center bin2_center hist hist_err weighted_counts weighted_counts_err counts counts_err";
  fout << "lower_edge1 upper_edge1 lower_edge2 upper_edge2 histogram_normalization";
  fout << "var1 var1_err var2 var2_err" << endl;

  for (size_t i=0; i<m_nbins1; i++) {
    for (size_t j=0; j<m_nbins2; j++) {
      fout << bin(i) << " " << bin2(j) << " " <<  this->operator()(i, j, hist_type, fact) << " " << this->error(i, j, hist_type, fact);
      fout << " " << m_histo.get()->bin[i*m_nbins2+j] << " " << sqrt(m_histo_error.get()->bin[i*m_nbins2+j]);
      fout << " " << m_counts[i][j] << " " << sqrt(m_counts[i][j]);
      fout << " " << edge1(i) << " " << edge1(i+1) << " " << edge2(j) << " " << edge2(j+1) << " " << normalization(i, j, hist_type, fact) << endl;
      fout << " " << m_var1[i][j] << " " << sqrt(m_var1_err[i][j]/m_histo.get()->bin[i*m_nbins2+j]) << " " << m_var2[i][j] << " " << sqrt(m_var2_err[i][j]/m_histo.get()->bin[i*m_nbins2+j]) << endl;
    }
  }

  fout.clear(); fout.close();
}


