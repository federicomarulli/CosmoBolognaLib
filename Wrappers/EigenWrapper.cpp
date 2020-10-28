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
 *  @file Wrappers/EigenWrapper.cpp
 *
 *  @brief functions that wrap Eigen routines for vector
 *  and matrix manipulation
 *
 *  This file contains the implementation of 
 *  wrappers of Eigen routines for vector
 *  and matrix manipulation
 *
 *  @author Alfonso Veropalumbo
 *
 *  @author alfonso.veropalumbo@uniroma3.it
 */

#include "EigenWrapper.h"

using namespace std;

using namespace cbl;
using namespace wrapper;
using namespace eigen;


// ============================================================================


std::vector<double> cbl::wrapper::eigen::EigenToVector (const Eigen::MatrixXd vec)
{
  vector<double> vector(vec.data(), vec.data()+vec.rows()*vec.cols());
  return vector;
}


// ============================================================================


vector<vector<double>> cbl::wrapper::eigen::EigenToMatrix (const Eigen::MatrixXd mat)
{
  vector<vector<double>> matrix;

  for (int i=0; i<mat.rows(); i++)
    matrix.push_back(EigenToVector(mat.row(i)));

  return matrix;
}


// ============================================================================


Eigen::MatrixXd cbl::wrapper::eigen::VectorToEigen (const std::vector<double> vec)
{ 
  return Eigen::VectorXd::Map(vec.data(), vec.size());
}


// ============================================================================


Eigen::MatrixXd cbl::wrapper::eigen::MatrixToEigen (const std::vector<std::vector<double>> mat)
{
  const int rows = static_cast<int>(mat.size());
  const int cols = static_cast<int>(mat[0].size());

  Eigen::MatrixXd matrix(rows, cols);
  
  for (int i=0; i<matrix.rows(); i++)
    matrix.row(i) = Eigen::VectorXd::Map(mat[i].data(), mat.size());

  return matrix;
}


// ============================================================================


Eigen::MatrixXd cbl::wrapper::eigen::SquareMatrixToEigen (const std::vector<double> mat)
{
  int order = sqrt(mat.size());
  return MatrixToEigen(reshape(mat, order, order));
}
