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
 *  @file Func/Data1D_collection.cpp
 *
 *  @brief Methods of the class Data1D_collection
 *
 *  This file contains the implementation of the methods of the class
 *  Data1D_collection
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Data1D_collection.h"
using namespace cosmobl;


// ======================================================================================


cosmobl::Data1D_collection::Data1D_collection (const int n_data) {}


// ======================================================================================


cosmobl::Data1D_collection::Data1D_collection (const vector<double> x1, const vector<double> fx1, const vector<double> x2, const vector<double> fx2, const double x1_min, const double x1_max, const double x2_min, const double x2_max) {}


// ======================================================================================


cosmobl::Data1D_collection::Data1D_collection (const vector<double> x1, const vector<double> fx1, const vector<double> x2, const vector<double> fx2, const vector<double> x3, const vector<double> fx3, const double x1_min, const double x1_max, const double x2_min, const double x2_max , const double x3_min, const double x3_max) {}


// ======================================================================================


cosmobl::Data1D_collection::Data1D_collection (const vector< vector< double > > x, const vector< vector< double > > fx, const vector<double> x_min, const vector<double> x_max) {}


// ======================================================================================


cosmobl::Data1D_collection::Data1D_collection (const vector<double> x1, const vector<double> fx1, const vector<double> error_fx1, const vector<double> x2, const vector<double> fx2, const vector<double> error_fx2, const double x1_min, const double x1_max, const double x2_min, const double x2_max ) {}


// ======================================================================================


cosmobl::Data1D_collection::Data1D_collection (const vector<double> x1, const vector<double> fx1, const vector<double> error_fx1, const vector<double> x2, const vector<double> fx2, const vector<double> error_fx2, const vector<double> x3, const vector<double> fx3, const vector<double> error_fx3, const double x1_min, const double x1_max, const double x2_min, const double x2_max, const double x3_min, const double x3_max) {}


// ======================================================================================


cosmobl::Data1D_collection::Data1D_collection (const vector< vector< double > > x, const vector< vector< double > > fx, const vector< vector< double > > error_fx, const vector<double> x_min, const vector<double> x_max) {}


// ======================================================================================


cosmobl::Data1D_collection::Data1D_collection (const vector<double> x1, const vector<double> fx1, const vector<vector<double>> covariance_fx1, const vector<double> x2, const vector<double> fx2, const vector<vector<double>> covariance_fx2, const vector<vector<double> > covariance_12 , const double x1_min, const double x1_max, const double x2_min, const double x2_max ) {}


// ======================================================================================


cosmobl::Data1D_collection::Data1D_collection (const vector<double> x1, const vector<double> fx1, const vector<vector<double>> covariance_fx1, const vector<double> x2, const vector<double> fx2, const vector<vector<double>> covariance_fx2, const vector<double> x3, const vector<double> fx3, const vector<vector<double>> covariance_fx3, const vector<vector<double> > covariance_12 , const vector<vector<double> > covariance_13 , const vector<vector<double> > covariance_23 ,double x1_min, const double x1_max, const double x2_min, const double x2_max , const double x3_min, const double x3_max) {}


// ======================================================================================


cosmobl::Data1D_collection::Data1D_collection (const vector< Data1D > data) {}


// ======================================================================================


cosmobl::Data1D_collection::Data1D_collection (const vector<string> input_file, const vector<double> x_min, const vector<double> x_max) {}


// ======================================================================================


double cosmobl::Data1D_collection::covariance_fx (const int i, const int j, const int ax1, const int ax2) const { return 1; }


// ======================================================================================


double cosmobl::Data1D_collection::inverse_covariance_fx (const int i, const int j, const int ax1, const int ax2) const { return 1; }


// ======================================================================================


void cosmobl::Data1D_collection::set_limits (const int i, const double xmin, const double xmax) {}


// ======================================================================================


void cosmobl::Data1D_collection::set_covariance_fx (const int i, const string filename) {}


// ======================================================================================


void cosmobl::Data1D_collection::set_covariance_fx (const int i, const vector<vector<double> > covariance_fx) {}


// ======================================================================================


void cosmobl::Data1D_collection::set_covariance_fx (const int i, const int j, const vector<vector<double> > covariance_fx) {} 


// ======================================================================================


int cosmobl::Data1D_collection::ndata_eff () const { return 1; } 


// ======================================================================================


int cosmobl::Data1D_collection::ndata () const { return 1; }  
