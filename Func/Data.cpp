/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Func/Data.cpp
 *
 *  @brief Methods of the class Data
 *
 *  This file contains the implementation of the methods of the class
 *  Data
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Data.h"
#include "Data1D.h"
#include "Data2D.h"
#include "Data1D_collection.h"

using namespace cosmobl;


// ======================================================================================


shared_ptr<Data> cosmobl::Data::Create (const dataType dType)
{
  if (dType==dataType::_1D_data_) return move(unique_ptr<Data1D>(new Data1D()));
  else if (dType==dataType::_1D_collection_data_) return move(unique_ptr<Data1D_collection>(new Data1D_collection()));
  else if (dType==dataType::_2D_data_) return move(unique_ptr<Data2D>(new Data2D()));
  else ErrorMsg("Error in cosmobl::Data::Create of Data.cpp: no such type of object, or error in the input parameters!!");

  return NULL;
}

// ======================================================================================


shared_ptr<Data> cosmobl::Data::Create (const string input_file, const double xmin, const double xmax)
{
  return move(unique_ptr<Data1D>(new Data1D(input_file, xmin, xmax)));
}


// ======================================================================================


shared_ptr<Data> cosmobl::Data::Create (const vector<double> x, const vector<double> fx, const double xmin, const double xmax)
{
  return move(unique_ptr<Data1D>(new Data1D(x, fx, xmin, xmax)));
}


// ======================================================================================


shared_ptr<Data> cosmobl::Data::Create (const vector<double> x, const vector<double> fx, const vector<double> error_fx, const double xmin, const double xmax)
{
  return move(unique_ptr<Data1D>(new Data1D(x, fx, error_fx, xmin, xmax)));
}


// ======================================================================================


shared_ptr<Data> cosmobl::Data::Create (const vector<double> x, const vector<double> fx, const vector<vector<double> > covariance_fx, const double xmin, const double xmax)
{
  return move(unique_ptr<Data1D>(new Data1D(x, fx, covariance_fx, xmin, xmax)));
}


// ======================================================================================


shared_ptr<Data> cosmobl::Data::Create (vector<double> x, vector<double> y, vector<vector<double> > fxy, double xmin, double xmax, double ymin, double ymax)
{
  return move(unique_ptr<Data2D>(new Data2D(x, y, fxy, xmin, xmax, ymin, ymax)));
}


// ======================================================================================


shared_ptr<Data> cosmobl::Data::Create (vector<double> x, vector<double> y, vector<vector<double> > fxy, vector<vector<double> > error_fxy, double xmin, double xmax, double ymin, double ymax)
{
  return move(unique_ptr<Data2D>(new Data2D(x, y, fxy, error_fxy, xmin, xmax, ymin, ymax)));
}
