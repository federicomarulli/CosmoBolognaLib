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
 *  @author federico.marulli3@unibo.it
 */


#include "Catalogue.h"

using namespace std;

using namespace cbl;
using namespace catalogue;


// ============================================================================


std::shared_ptr<Object> cbl::catalogue::Object::Create (const ObjectType objType)
{
  if (objType==ObjectType::_Random_) return move(unique_ptr<RandomObject>(new RandomObject));
  else if (objType==ObjectType::_Mock_) return move(unique_ptr<Mock>(new Mock()));
  else if (objType==ObjectType::_Halo_) return move(unique_ptr<Halo>(new Halo()));
  else if (objType==ObjectType::_Galaxy_) return move(unique_ptr<Galaxy>(new Galaxy()));
  else if (objType==ObjectType::_Cluster_) return move(unique_ptr<Cluster>(new Cluster()));
  else if (objType==ObjectType::_Void_) return move(unique_ptr<Void>(new Void()));
  else if (objType==ObjectType::_HostHalo_) return move(unique_ptr<HostHalo>(new HostHalo()));
  else if (objType==ObjectType::_ChainMeshCell_) return move(unique_ptr<ChainMeshCell>(new ChainMeshCell()));
  else ErrorCBL("no such type of object!", "Create", "Object.cpp");
  return NULL;
}


// ============================================================================


std::shared_ptr<Object> cbl::catalogue::Object::Create (const ObjectType objType, const comovingCoordinates coord, const double weight, const long region, const int ID, const std::string field, const double x_displacement, const double y_displacement, const double z_displacement)
{
  if (objType==ObjectType::_Random_) return move(unique_ptr<RandomObject>(new RandomObject(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Mock_) return move(unique_ptr<Mock>(new Mock(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Halo_) return move(unique_ptr<Halo>(new Halo(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Galaxy_) return move(unique_ptr<Galaxy>(new Galaxy(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Cluster_) return move(unique_ptr<Cluster>(new Cluster(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Void_) return move(unique_ptr<Void>(new Void(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_HostHalo_) return move(unique_ptr<HostHalo>(new HostHalo(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else ErrorCBL("no such type of object!", "Create", "Object.cpp");
  return NULL;
}


// ============================================================================


std::shared_ptr<Object> cbl::catalogue::Object::Create (const ObjectType objType, const comovingCoordinates coord, const cosmology::Cosmology &cosm, const double z1_guess, const double z2_guess, const double weight, const long region, const int ID, const std::string field, const double x_displacement, const double y_displacement, const double z_displacement)
{
  if (objType==ObjectType::_Random_) return move(unique_ptr<RandomObject>(new RandomObject(coord, cosm, z1_guess, z2_guess, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Mock_) return move(unique_ptr<Mock>(new Mock(coord, cosm, z1_guess, z2_guess, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Halo_) return move(unique_ptr<Halo>(new Halo(coord, cosm, z1_guess, z2_guess, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Galaxy_) return move(unique_ptr<Galaxy>(new Galaxy(coord, cosm, z1_guess, z2_guess, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Cluster_) return move(unique_ptr<Cluster>(new Cluster(coord, cosm, z1_guess, z2_guess, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Void_) return move(unique_ptr<Void>(new Void(coord, cosm, z1_guess, z2_guess, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_HostHalo_) return move(unique_ptr<HostHalo>(new HostHalo(coord, cosm, z1_guess, z2_guess, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else ErrorCBL("no such type of object!", "Create", "Object.cpp");
  return NULL;
}

// ============================================================================


std::shared_ptr<Object> cbl::catalogue::Object::Create (const ObjectType objType, const observedCoordinates coord, const double weight, const long region, const int ID, const std::string field, const double x_displacement, const double y_displacement, const double z_displacement)
{
  if (objType==ObjectType::_Random_) return move(unique_ptr<RandomObject>(new RandomObject(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Mock_) return move(unique_ptr<Mock>(new Mock(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Halo_) return move(unique_ptr<Halo>(new Halo(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Galaxy_) return move(unique_ptr<Galaxy>(new Galaxy(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Cluster_) return move(unique_ptr<Cluster>(new Cluster(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Void_) return move(unique_ptr<Void>(new Void(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_HostHalo_) return move(unique_ptr<HostHalo>(new HostHalo(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else ErrorCBL("no such type of object!", "Create", "Object.cpp");
  return NULL;
}

// ============================================================================


std::shared_ptr<Object> cbl::catalogue::Object::Create (const ObjectType objType, const observedCoordinates coord, const CoordinateUnits inputUnits, const double weight, const long region, const int ID, const std::string field, const double x_displacement, const double y_displacement, const double z_displacement)
{
  if (objType==ObjectType::_Random_) return move(unique_ptr<RandomObject>(new RandomObject(coord, inputUnits, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Mock_) return move(unique_ptr<Mock>(new Mock(coord, inputUnits, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Halo_) return move(unique_ptr<Halo>(new Halo(coord, inputUnits, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Galaxy_) return move(unique_ptr<Galaxy>(new Galaxy(coord, inputUnits, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Cluster_) return move(unique_ptr<Cluster>(new Cluster(coord, inputUnits, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Void_) return move(unique_ptr<Void>(new Void(coord, inputUnits, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_HostHalo_) return move(unique_ptr<HostHalo>(new HostHalo(coord, inputUnits, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else ErrorCBL("no such type of object!", "Create", "Object.cpp");
  return NULL;
}

// ============================================================================


std::shared_ptr<Object> cbl::catalogue::Object::Create (const ObjectType objType, const observedCoordinates coord, const cosmology::Cosmology &cosm, const double weight, const long region, const int ID, const std::string field, const double x_displacement, const double y_displacement, const double z_displacement)
{
  if (objType==ObjectType::_Random_) return move(unique_ptr<RandomObject>(new RandomObject(coord, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Mock_) return move(unique_ptr<Mock>(new Mock(coord, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Halo_) return move(unique_ptr<Halo>(new Halo(coord, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Galaxy_) return move(unique_ptr<Galaxy>(new Galaxy(coord, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Cluster_) return move(unique_ptr<Cluster>(new Cluster(coord, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Void_) return move(unique_ptr<Void>(new Void(coord, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_HostHalo_) return move(unique_ptr<HostHalo>(new HostHalo(coord, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else ErrorCBL("no such type of object!", "Create", "Object.cpp");
  return NULL;
}

// ============================================================================


std::shared_ptr<Object> cbl::catalogue::Object::Create (const ObjectType objType, const observedCoordinates coord, const CoordinateUnits inputUnits, const cosmology::Cosmology &cosm, const double weight, const long region, const int ID, const std::string field, const double x_displacement, const double y_displacement, const double z_displacement)
{
  if (objType==ObjectType::_Random_) return move(unique_ptr<RandomObject>(new RandomObject(coord, inputUnits, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Mock_) return move(unique_ptr<Mock>(new Mock(coord, inputUnits, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Halo_) return move(unique_ptr<Halo>(new Halo(coord, inputUnits, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Galaxy_) return move(unique_ptr<Galaxy>(new Galaxy(coord, inputUnits, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Cluster_) return move(unique_ptr<Cluster>(new Cluster(coord, inputUnits, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Void_) return move(unique_ptr<Void>(new Void(coord, inputUnits, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_HostHalo_) return move(unique_ptr<HostHalo>(new HostHalo(coord, inputUnits, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else ErrorCBL("no such type of object!", "Create", "Object.cpp");
  return NULL;
}

// ============================================================================


std::shared_ptr<Object> cbl::catalogue::Object::Create (const ObjectType objType, const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight, const long region, const int ID, const std::string field, const double x_displacement, const double y_displacement, const double z_displacement)
{
  if (objType==ObjectType::_Random_) return move(unique_ptr<RandomObject>(new RandomObject(xx, yy, zz, ra, dec, redshift, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Mock_) return move(unique_ptr<Mock>(new Mock(xx, yy, zz, ra, dec, redshift, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Halo_) return move(unique_ptr<Halo>(new Halo(xx, yy, zz, ra, dec, redshift, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Galaxy_) return move(unique_ptr<Galaxy>(new Galaxy(xx, yy, zz, ra, dec, redshift, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Cluster_) return move(unique_ptr<Cluster>(new Cluster(xx, yy, zz, ra, dec, redshift, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_Void_) return move(unique_ptr<Void>(new Void(xx, yy, zz, ra, dec, redshift, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else if (objType==ObjectType::_HostHalo_) return move(unique_ptr<HostHalo>(new HostHalo(xx, yy, zz, ra, dec, redshift, weight, region, ID, field, x_displacement, y_displacement, z_displacement)));
  else ErrorCBL("no such type of object!", "Create", "Object.cpp");
  return NULL;
}

// ============================================================================
std::shared_ptr<Object> cbl::catalogue::Object::Create (const comovingCoordinates coord, const int ID, const std::vector<unsigned int> part, std::vector<std::vector<unsigned int>> nearCells)
{
  return move(unique_ptr<ChainMeshCell>(new ChainMeshCell(coord, ID, part, nearCells)));
}


