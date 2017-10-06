/*******************************************************************
 *  Copyright (C) 2015 by Federico Marulli                         *
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
 *  @file Measure/TwoPointCorrelation/Pair.cpp
 *
 *  @brief Methods of the class Pair  
 *
 *  This file contains the implementation of all the methods of the
 *  classes Pair*, used to handle pairs of objects of any kind
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Pair1D_extra.h"
#include "Pair2D_extra.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;


// ============================================================================================


shared_ptr<Pair> cosmobl::pairs::Pair::Create (const PairType type, const PairInfo info, const double Min, const double Max, const int nbins, const double shift, const CoordUnits angularUnits, function<double(double)> angularWeight)
{
  if (info==_standard_) {
    
    if (type==_angular_lin_)
      return move(unique_ptr<Pair1D_angular_lin>{new Pair1D_angular_lin(Min, Max, nbins, shift, angularUnits, angularWeight)});

    else if (type==_angular_log_)
      return move(unique_ptr<Pair1D_angular_log>{new Pair1D_angular_log(Min, Max, nbins, shift, angularUnits, angularWeight)});

    else if (type==_comoving_lin_)
      return move(unique_ptr<Pair1D_comoving_lin>{new Pair1D_comoving_lin(Min, Max, nbins, shift, angularUnits, angularWeight)});

    else if (type==_comoving_log_) 
      return move(unique_ptr<Pair1D_comoving_log>{new Pair1D_comoving_log(Min, Max, nbins, shift, angularUnits, angularWeight)});
  
    else ErrorCBL("Error in cosmobl::pairs::Create of Pairs.cpp: no such type of object!");

  }

  else if (info==_extra_) {

    if (type==_angular_lin_)
      return move(unique_ptr<Pair1D_angular_lin_extra>{new Pair1D_angular_lin_extra(Min, Max, nbins, shift, angularUnits, angularWeight)});
    
    else if (type==_angular_log_)
      return move(unique_ptr<Pair1D_angular_log_extra>{new Pair1D_angular_log_extra(Min, Max, nbins, shift, angularUnits, angularWeight)});
    
    else if (type==_comoving_lin_)
      return move(unique_ptr<Pair1D_comoving_lin_extra>{new Pair1D_comoving_lin_extra(Min, Max, nbins, shift, angularUnits, angularWeight)});
    
    else if (type==_comoving_log_)
      return move(unique_ptr<Pair1D_comoving_log_extra>{new Pair1D_comoving_log_extra(Min, Max, nbins, shift, angularUnits, angularWeight)});

    else ErrorCBL("Error in cosmobl::pairs::Create of Pairs.cpp: no such type of object!");

  }
  
  else ErrorCBL("Error in cosmobl::pairs::Create of Pairs.cpp: no such type of object!");
   
  return NULL;
}


// ============================================================================================


shared_ptr<Pair> cosmobl::pairs::Pair::Create (const PairType type, const PairInfo info, const double Min, const double Max, const double binSize, const double shift, const CoordUnits angularUnits, function<double(double)> angularWeight)
{
  if (info==_standard_) {
  
    if (type==_angular_lin_)
      return move(unique_ptr<Pair1D_angular_lin>{new Pair1D_angular_lin(Min, Max, binSize, shift, angularUnits, angularWeight)});

    else if (type==_angular_log_)
      return move(unique_ptr<Pair1D_angular_log>{new Pair1D_angular_log(Min, Max, binSize, shift, angularUnits, angularWeight)});

    else if (type==_comoving_lin_)
      return move(unique_ptr<Pair1D_comoving_lin>{new Pair1D_comoving_lin(Min, Max, binSize, shift, angularUnits, angularWeight)});

    else if (type==_comoving_log_)
      return move(unique_ptr<Pair1D_comoving_log>{new Pair1D_comoving_log(Min, Max, binSize, shift, angularUnits, angularWeight)});
    
    else ErrorCBL("Error in cosmobl::pairs::Create of Pairs.cpp: no such type of object!");

  }
  
  else if (info==_extra_) {

    if (type==_angular_lin_)
      return move(unique_ptr<Pair1D_angular_lin_extra>{new Pair1D_angular_lin_extra(Min, Max, binSize, shift, angularUnits, angularWeight)});

    else if (type==_angular_log_)
      return move(unique_ptr<Pair1D_angular_log_extra>{new Pair1D_angular_log_extra(Min, Max, binSize, shift, angularUnits, angularWeight)});

    else if (type==_comoving_lin_)
      return move(unique_ptr<Pair1D_comoving_lin_extra>{new Pair1D_comoving_lin_extra(Min, Max, binSize, shift, angularUnits, angularWeight)});

    else if (type==_comoving_log_) 
      return move(unique_ptr<Pair1D_comoving_log_extra>{new Pair1D_comoving_log_extra(Min, Max, binSize, shift, angularUnits, angularWeight)});
  
    else ErrorCBL("Error in cosmobl::pairs::Create of Pairs.cpp: no such type of object!");
    
  }

  else ErrorCBL("Error in cosmobl::pairs::Create of Pairs.cpp: no such type of object!");
  
  return NULL;
}


// ============================================================================================


shared_ptr<Pair> cosmobl::pairs::Pair::Create (const PairType type, const PairInfo info, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const double Min_D2, const double Max_D2, const int nbins_D2, const double shift_D2, const CoordUnits angularUnits, function<double(double)> angularWeight)
{
  if (info==_standard_) {
  
    if (type==_comovingCartesian_linlin_)
      return move(unique_ptr<Pair2D_comovingCartesian_linlin>{new Pair2D_comovingCartesian_linlin(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});
  
    else if (type==_comovingCartesian_linlog_)
      return move(unique_ptr<Pair2D_comovingCartesian_linlog>{new Pair2D_comovingCartesian_linlog(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingCartesian_loglin_)
      return move(unique_ptr<Pair2D_comovingCartesian_loglin>{new Pair2D_comovingCartesian_loglin(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingCartesian_loglog_)
      return move(unique_ptr<Pair2D_comovingCartesian_loglog>{new Pair2D_comovingCartesian_loglog(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});
    
    else if (type==_comovingPolar_linlin_)
      return move(unique_ptr<Pair2D_comovingPolar_linlin>{new Pair2D_comovingPolar_linlin(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});
  
    else if (type==_comovingPolar_linlog_)
      return move(unique_ptr<Pair2D_comovingPolar_linlog>{new Pair2D_comovingPolar_linlog(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingPolar_loglin_)
      return move(unique_ptr<Pair2D_comovingPolar_loglin>{new Pair2D_comovingPolar_loglin(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingPolar_loglog_)
      return move(unique_ptr<Pair2D_comovingPolar_loglog>{new Pair2D_comovingPolar_loglog(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});
  
    else ErrorCBL("Error in cosmobl::pairs::Create of Pairs.cpp: no such type of object!");

  }

  else if (info==_extra_) {
    
    if (type==_comovingCartesian_linlin_)
      return move(unique_ptr<Pair2D_comovingCartesian_linlin_extra>{new Pair2D_comovingCartesian_linlin_extra(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingCartesian_linlog_)
      return move(unique_ptr<Pair2D_comovingCartesian_linlog_extra>{new Pair2D_comovingCartesian_linlog_extra(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingCartesian_loglin_)
      return move(unique_ptr<Pair2D_comovingCartesian_loglin_extra>{new Pair2D_comovingCartesian_loglin_extra(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingCartesian_loglog_)
      return move(unique_ptr<Pair2D_comovingCartesian_loglog_extra>{new Pair2D_comovingCartesian_loglog_extra(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingPolar_linlin_)
      return move(unique_ptr<Pair2D_comovingPolar_linlin_extra>{new Pair2D_comovingPolar_linlin_extra(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingPolar_linlog_)
      return move(unique_ptr<Pair2D_comovingPolar_linlog_extra>{new Pair2D_comovingPolar_linlog_extra(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingPolar_loglin_)
      return move(unique_ptr<Pair2D_comovingPolar_loglin_extra>{new Pair2D_comovingPolar_loglin_extra(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingPolar_loglog_)
      return move(unique_ptr<Pair2D_comovingPolar_loglog_extra>{new Pair2D_comovingPolar_loglog_extra(Min_D1, Max_D1, nbins_D1, shift_D1, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight)});
  
    else ErrorCBL("Error in cosmobl::pairs::Create of Pairs.cpp: no such type of object!");
	  
  }

  else ErrorCBL("Error in cosmobl::pairs::Create of Pairs.cpp: no such type of object!");
  
  return NULL;
}


// ============================================================================================


shared_ptr<Pair> cosmobl::pairs::Pair::Create (const PairType type, const PairInfo info, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const double Min_D2, const double Max_D2, const double binSize_D2, const double shift_D2, const CoordUnits angularUnits, function<double(double)> angularWeight)
{
  if (info==_standard_) {
    
    if (type==_comovingCartesian_linlin_)
      return move(unique_ptr<Pair2D_comovingCartesian_linlin>{new Pair2D_comovingCartesian_linlin(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingCartesian_linlog_)
      return move(unique_ptr<Pair2D_comovingCartesian_linlog>{new Pair2D_comovingCartesian_linlog(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingCartesian_loglin_)
      return move(unique_ptr<Pair2D_comovingCartesian_loglin>{new Pair2D_comovingCartesian_loglin(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingCartesian_loglog_)
      return move(unique_ptr<Pair2D_comovingCartesian_loglog>{new Pair2D_comovingCartesian_loglog(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingPolar_linlin_)
      return move(unique_ptr<Pair2D_comovingPolar_linlin>{new Pair2D_comovingPolar_linlin(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});
  
    else if (type==_comovingPolar_linlog_)
      return move(unique_ptr<Pair2D_comovingPolar_linlog>{new Pair2D_comovingPolar_linlog(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingPolar_loglin_)
      return move(unique_ptr<Pair2D_comovingPolar_loglin>{new Pair2D_comovingPolar_loglin(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingPolar_loglog_)
      return move(unique_ptr<Pair2D_comovingPolar_loglog>{new Pair2D_comovingPolar_loglog(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

  }

  else if (info==_extra_) {

    if (type==_comovingCartesian_linlin_)
      return move(unique_ptr<Pair2D_comovingCartesian_linlin_extra>{new Pair2D_comovingCartesian_linlin_extra(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingCartesian_linlog_)
      return move(unique_ptr<Pair2D_comovingCartesian_linlog_extra>{new Pair2D_comovingCartesian_linlog_extra(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingCartesian_loglin_)
      return move(unique_ptr<Pair2D_comovingCartesian_loglin_extra>{new Pair2D_comovingCartesian_loglin_extra(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingCartesian_loglog_)
      return move(unique_ptr<Pair2D_comovingCartesian_loglog_extra>{new Pair2D_comovingCartesian_loglog_extra(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingPolar_linlin_)
      return move(unique_ptr<Pair2D_comovingPolar_linlin_extra>{new Pair2D_comovingPolar_linlin_extra(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingPolar_linlog_)
      return move(unique_ptr<Pair2D_comovingPolar_linlog_extra>{new Pair2D_comovingPolar_linlog_extra(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingPolar_loglin_)
      return move(unique_ptr<Pair2D_comovingPolar_loglin_extra>{new Pair2D_comovingPolar_loglin_extra(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});

    else if (type==_comovingPolar_loglog_)
      return move(unique_ptr<Pair2D_comovingPolar_loglog_extra>{new Pair2D_comovingPolar_loglog_extra(Min_D1, Max_D1, binSize_D1, shift_D1, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight)});
    
  }
  
  else ErrorCBL("Error in cosmobl::pairs::Create of Pairs.cpp: no such type of object!");
 
  return NULL;
}



