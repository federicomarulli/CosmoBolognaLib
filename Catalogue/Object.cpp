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


shared_ptr<Object> cosmobl::catalogue::Object::Create (const ObjType objType, const comovingCoordinates coord, const double weight)
{
  if (objType==_RandomObject_) return move(unique_ptr<RandomObject>(new RandomObject(coord, weight)));
  else if (objType==_Mock_)    return move(unique_ptr<Mock>(new Mock(coord, weight)));
  else if (objType==_Halo_)    return move(unique_ptr<Halo>(new Halo(coord, weight)));
  else if (objType==_Galaxy_)  return move(unique_ptr<Galaxy>(new Galaxy(coord, weight)));
  else if (objType==_Cluster_) return move(unique_ptr<Cluster>(new Cluster(coord, weight)));
  else if (objType==_Void_)    return move(unique_ptr<Void>(new Void(coord, weight)));
  else ErrorMsg("Error in Create of Object.cpp: no such type of object!");
  return NULL;
}


// ============================================================================


shared_ptr<Object> cosmobl::catalogue::Object::Create (const ObjType objType, const comovingCoordinates coord, const cosmology::Cosmology &cosm, const double z1_guess, const double z2_guess, const double weight)
{
  if (objType==_RandomObject_) return move(unique_ptr<RandomObject>(new RandomObject(coord, cosm, z1_guess, z2_guess, weight)));
  else if (objType==_Mock_)    return move(unique_ptr<Mock>(new Mock(coord, cosm, z1_guess, z2_guess, weight)));
  else if (objType==_Halo_)    return move(unique_ptr<Halo>(new Halo(coord, cosm, z1_guess, z2_guess, weight)));
  else if (objType==_Galaxy_)  return move(unique_ptr<Galaxy>(new Galaxy(coord, cosm, z1_guess, z2_guess, weight)));
  else if (objType==_Cluster_) return move(unique_ptr<Cluster>(new Cluster(coord, cosm, z1_guess, z2_guess, weight)));
  else if (objType==_Void_)    return move(unique_ptr<Void>(new Void(coord, cosm, z1_guess, z2_guess, weight)));
  else ErrorMsg("Error in Create of Object.cpp: no such type of object!");
  return NULL;
}

// ============================================================================


shared_ptr<Object> cosmobl::catalogue::Object::Create (const ObjType objType, const observedCoordinates coord, const double weight)
{
  if (objType==_RandomObject_) return move(unique_ptr<RandomObject>(new RandomObject(coord, weight)));
  else if (objType==_Mock_)    return move(unique_ptr<Mock>(new Mock(coord, weight)));
  else if (objType==_Halo_)    return move(unique_ptr<Halo>(new Halo(coord, weight)));
  else if (objType==_Galaxy_)  return move(unique_ptr<Galaxy>(new Galaxy(coord, weight)));
  else if (objType==_Cluster_) return move(unique_ptr<Cluster>(new Cluster(coord, weight)));
  else if (objType==_Void_)    return move(unique_ptr<Void>(new Void(coord, weight)));
  else ErrorMsg("Error in Create of Object.cpp: no such type of object!");
  return NULL;
}

// ============================================================================


shared_ptr<Object> cosmobl::catalogue::Object::Create (const ObjType objType, const observedCoordinates coord, const CoordUnits inputUnits, const double weight)
{
  if (objType==_RandomObject_) return move(unique_ptr<RandomObject>(new RandomObject(coord, inputUnits, weight)));
  else if (objType==_Mock_)    return move(unique_ptr<Mock>(new Mock(coord, inputUnits, weight)));
  else if (objType==_Halo_)    return move(unique_ptr<Halo>(new Halo(coord, inputUnits, weight)));
  else if (objType==_Galaxy_)  return move(unique_ptr<Galaxy>(new Galaxy(coord, inputUnits, weight)));
  else if (objType==_Cluster_) return move(unique_ptr<Cluster>(new Cluster(coord, inputUnits, weight)));
  else if (objType==_Void_)    return move(unique_ptr<Void>(new Void(coord, inputUnits, weight)));
  else ErrorMsg("Error in Create of Object.cpp: no such type of object!");
  return NULL;
}

// ============================================================================


shared_ptr<Object> cosmobl::catalogue::Object::Create (const ObjType objType, const observedCoordinates coord, const cosmology::Cosmology &cosm, const double weight)
{
  if (objType==_RandomObject_) return move(unique_ptr<RandomObject>(new RandomObject(coord, cosm, weight)));
  else if (objType==_Mock_)    return move(unique_ptr<Mock>(new Mock(coord, cosm, weight)));
  else if (objType==_Halo_)    return move(unique_ptr<Halo>(new Halo(coord, cosm, weight)));
  else if (objType==_Galaxy_)  return move(unique_ptr<Galaxy>(new Galaxy(coord, cosm, weight)));
  else if (objType==_Cluster_) return move(unique_ptr<Cluster>(new Cluster(coord, cosm, weight)));
  else if (objType==_Void_)    return move(unique_ptr<Void>(new Void(coord, cosm, weight)));
  else ErrorMsg("Error in Create of Object.cpp: no such type of object!");
  return NULL;
}

// ============================================================================


shared_ptr<Object> cosmobl::catalogue::Object::Create (const ObjType objType, const observedCoordinates coord, const CoordUnits inputUnits, const cosmology::Cosmology &cosm, const double weight)
{
  if (objType==_RandomObject_) return move(unique_ptr<RandomObject>(new RandomObject(coord, inputUnits, cosm, weight)));
  else if (objType==_Mock_)    return move(unique_ptr<Mock>(new Mock(coord, inputUnits, cosm, weight)));
  else if (objType==_Halo_)    return move(unique_ptr<Halo>(new Halo(coord, inputUnits, cosm, weight)));
  else if (objType==_Galaxy_)  return move(unique_ptr<Galaxy>(new Galaxy(coord, inputUnits, cosm, weight)));
  else if (objType==_Cluster_) return move(unique_ptr<Cluster>(new Cluster(coord, inputUnits, cosm, weight)));
  else if (objType==_Void_)    return move(unique_ptr<Void>(new Void(coord, inputUnits, cosm, weight)));
  else ErrorMsg("Error in Create of Object.cpp: no such type of object!");
  return NULL;
}

// ============================================================================


shared_ptr<Object> cosmobl::catalogue::Object::Create (const ObjType objType, const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight)
{
  if (objType==_RandomObject_) return move(unique_ptr<RandomObject>(new RandomObject(xx, yy, zz, ra, dec, redshift, weight)));
  else if (objType==_Mock_)    return move(unique_ptr<Mock>(new Mock(xx, yy, zz, ra, dec, redshift, weight)));
  else if (objType==_Halo_)    return move(unique_ptr<Halo>(new Halo(xx, yy, zz, ra, dec, redshift, weight)));
  else if (objType==_Galaxy_)  return move(unique_ptr<Galaxy>(new Galaxy(xx, yy, zz, ra, dec, redshift, weight)));
  else if (objType==_Cluster_) return move(unique_ptr<Cluster>(new Cluster(xx, yy, zz, ra, dec, redshift, weight)));
  else if (objType==_Void_)    return move(unique_ptr<Void>(new Void(xx, yy, zz, ra, dec, redshift, weight)));
  else ErrorMsg("Error in Create of Object.cpp: no such type of object!");
  return NULL;
}


