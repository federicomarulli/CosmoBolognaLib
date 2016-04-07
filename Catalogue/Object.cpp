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
 *  @file Catalogue/Object.cpp
 *
 *  @brief Methods of the class Object 
 *
 *  This file contains the implementation of the methods of the class
 *  Object, used to handle astronomical sources
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */


#include "Catalogue.h"

using namespace cosmobl;
using namespace catalogue;


// ============================================================================


shared_ptr<Object> cosmobl::catalogue::Object::Create (const ObjType type, const double xx, const double yy, const double zz, const double weight)
{
  if (type==_RandomObject_) return move(unique_ptr<RandomObject>(new RandomObject(xx, yy, zz, weight)));
  else if (type==_Mock_)    return move(unique_ptr<Mock>(new Mock(xx, yy, zz, weight)));
  else if (type==_Halo_)    return move(unique_ptr<Halo>(new Halo(xx, yy, zz, weight)));
  else if (type==_Galaxy_)  return move(unique_ptr<Galaxy>(new Galaxy(xx, yy, zz, weight)));
  else if (type==_Cluster_) return move(unique_ptr<Cluster>(new Cluster(xx, yy, zz, weight)));
  else if (type==_Void_)    return move(unique_ptr<Void>(new Void(xx, yy, zz, weight)));
  else ErrorMsg("Error in Create of Object.cpp: no such type of object!");
  return NULL;
}


// ============================================================================


shared_ptr<Object> cosmobl::catalogue::Object::Create(const ObjType type, const double ra, const double dec, const double redshift, const Cosmology &cosm, const double weight)
{
  if (type==_RandomObject_) return move(unique_ptr<RandomObject>(new RandomObject(ra, dec, redshift, cosm, weight)));
  else if (type==_Mock_)    return move(unique_ptr<Mock>(new Mock(ra, dec, redshift, cosm, weight)));
  else if (type==_Halo_)    return move(unique_ptr<Halo>(new Halo(ra, dec, redshift, cosm, weight)));
  else if (type==_Galaxy_)  return move(unique_ptr<Galaxy>(new Galaxy(ra, dec, redshift, cosm, weight)));
  else if (type==_Cluster_) return move(unique_ptr<Cluster>(new Cluster(ra, dec, redshift, cosm, weight)));
  else if (type==_Void_)    return move(unique_ptr<Void>(new Void(ra, dec, redshift, cosm, weight)));
  else ErrorMsg("Error in Create of Object.cpp: no such type of object!");
  return NULL;
}


