/********************************************************************
 *  Copyright (C) 2020 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Distribution/CombinedDistribution.cpp
 *
 *  @brief Methods of the class CombinedDistribution
 *
 *  This file contains the implementation of the methods of the class
 *  CombinedDistribution
 *
 *  @author Sofia Contarini
 *
 *  @author sofia.contarini3@unibo.it
 */

#include "CombinedDistribution.h"


// ======================================================================================


cbl::glob::CombinedDistribution::CombinedDistribution (const DistributionType distributionType, const std::vector<double> meanVec, const std::vector<std::vector<double>> covMat, const std::vector<double> xMinVec, const std::vector<double> xMaxVec, const int seed) 
{
  if (distributionType!=glob::DistributionType::_Gaussian_)
    ErrorCBL("this constructor only allows DistributionType::_Gaussian_!", "CombinedDistribution", "CombinedDistribution.cpp");

  m_xMinVec = xMinVec;
  m_xMaxVec = xMaxVec;
  std::vector<double> sigmaVec = cbl::wrapper::eigen::EigenToVector(cbl::wrapper::eigen::MatrixToEigen(covMat).diagonal().array().sqrt());

  glob::STR_multivariateGaussian parameters;
  parameters.MeanVec = cbl::wrapper::eigen::VectorToEigen(meanVec);
  parameters.CovMat = cbl::wrapper::eigen::MatrixToEigen(covMat);

  m_inputs = std::make_shared<glob::STR_multivariateGaussian>(parameters);  
  m_Func = &multivariateGaussian;

  m_distributionVec.erase(m_distributionVec.begin(), m_distributionVec.end());
  m_distributionVec.resize(meanVec.size(), NULL);
  
  for (size_t i=0; i<meanVec.size(); i++)
    m_distributionVec[i] = std::make_shared<cbl::glob::Distribution>(cbl::glob::Distribution(cbl::glob::DistributionType::_Gaussian_, {meanVec[i], sigmaVec[i]}, xMinVec[i], xMaxVec[i], seed+i));

}

// ======================================================================================


cbl::glob::CombinedDistribution::CombinedDistribution (const std::string filename, const std::string path, const std::vector<int> columns_to_read, const int skip_nlines, const int type_data, const bool normalize, const int distNum, const double rMAX, const double cell_size) 
{
  int nparams = columns_to_read.size()-1;
  std::vector<std::vector<double>> all_data = cbl::read_file(filename, path, columns_to_read, skip_nlines);
  std::vector<double> posterior;

  switch (type_data) {
  case 0:
    coutCBL << "Reading log(Posterior) distribution as column number " << columns_to_read[nparams] << std::endl;
    posterior = all_data[nparams];
    if (cbl::Min(posterior)<log(DBL_MIN)) ErrorCBL("please select a valid higher log(Posterior) values!", "CombinedDistribution", "CombinedDistribution.cpp");
    for (size_t i=0; i<posterior.size(); i++) posterior[i] = exp(posterior[i]);
    break;
 
  case 1:
    coutCBL << "Reading Posterior distribution as column number " << columns_to_read[nparams] << std::endl;
    posterior = all_data[nparams+1];
    if (cbl::Min(posterior)<0.) ErrorCBL("please select a valid positive Posterior values!", "CombinedDistribution", "CombinedDistribution.cpp");
    break;

  case 2:
    coutCBL << "Reading Chi2 distribution as column number " << columns_to_read[nparams] << std::endl;
    posterior = all_data[nparams];
    if (cbl::Min(posterior)<-2.*log(DBL_MIN)) ErrorCBL("please select a valid higher log(Posterior) values!", "CombinedDistribution", "CombinedDistribution.cpp");
    for (size_t i=0; i<posterior.size(); i++) posterior[i] = exp(-0.5*posterior[i]);
    break;
  
  default:
    ErrorCBL("please select a valid type_data!", "CombinedDistribution", "CombinedDistribution.cpp");
  }

  if (normalize) {
    const double Min_post = cbl::Min(posterior);
    const double Max_post = cbl::Max(posterior);
    const double Norm = fabs(Max_post-Min_post); 
    for (size_t i=0; i<posterior.size(); i++) posterior[i] = (posterior[i]-Min_post)/Norm;
  }
  
  m_xMinVec.erase(m_xMinVec.begin(), m_xMinVec.end());
  m_xMaxVec.erase(m_xMaxVec.begin(), m_xMaxVec.end());
  m_xMinVec.resize(nparams);
  m_xMaxVec.resize(nparams);
  
  std::vector<std::vector<double>> data(nparams);
  for (int i=0; i<nparams; i++) {
    data[i] = all_data[i];
    m_xMinVec[i] = cbl::Min(data[i])*0.99;
    m_xMaxVec[i] = cbl::Max(data[i])*1.01;
  }
  
  cbl::chainmesh::ChainMesh chMesh(cell_size, nparams);
  chMesh.normalize(data, posterior, rMAX);

  glob::STR_chainMeshInterpolate parameters;
  parameters.ChainMesh = chMesh;
  parameters.DistNum = distNum;
  
  m_inputs = std::make_shared<glob::STR_chainMeshInterpolate>(parameters);  
  m_Func = &chainMeshInterpolate;

  m_distributionVec.erase(m_distributionVec.begin(), m_distributionVec.end());
  m_distributionVec.resize(nparams, NULL);
  
  for (int i=0; i<nparams; i++)
    m_distributionVec[i] = std::make_shared<cbl::glob::Distribution>(cbl::glob::Distribution(cbl::glob::DistributionType::_Uniform_, m_xMinVec[i], m_xMaxVec[i]));

}
      
// ======================================================================================


double cbl::glob::CombinedDistribution::operator[] (std::vector<double> xx)
{   
  double fact = m_Func(xx, m_inputs);
  
  if (fact==0.) fact = DBL_MIN;
  
  else if (!std::isfinite(fact)) {
    fact = DBL_MIN;
    WarningMsgCBL("inf or nan values encountered", "operator[]", "CombinedDistribution.cpp");
  }
  
  return fact;
}
