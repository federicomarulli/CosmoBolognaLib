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
 *  @file Headers/Object.h
 *
 *  @brief The class Object
 *
 *  This file defines the interface of the class Object, that defines
 *  the members and methods of the derived classes used to handle
 *  astronomical objects, and the static factories used to construct
 *  the objects
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __OBJECT__
#define __OBJECT__


// ============================================================================================


namespace cbl {

  namespace catalogue {

    /**
     * @enum ObjectType
     * @brief the object types
     */
    enum class ObjectType {

      /// random object
      _Random_,

      /// mock object
      _Mock_,

      /// dark matter halo
      _Halo_,

      /// galaxy
      _Galaxy_,

      /// galaxy cluster
      _Cluster_,

      /// cosmic void
      _Void_,

      /// host halo
      _HostHalo_

    };

    /**
     * @brief return a vector containing the
     * ObjectType names
     * @return a vector containing the
     * ObjectType names
     */
    inline std::vector<std::string> ObjectTypeNames ()
    { return {"Random", "Mock", "Halo", "Galaxy", "Cluster", "Void", "HostHalo"}; }

    /**
     * @brief cast an enum of type ObjectType
     * from its index
     * @param objectTypeIndex the objectType index
     * @return object of class ObjectType
     */
    inline ObjectType ObjectTypeCast (const int objectTypeIndex)
    { return castFromValue<ObjectType>(objectTypeIndex); }

    /**
     * @brief cast an enum of type ObjectType
     * from its name
     * @param objectTypeName the objectType name
     * @return object of class ObjectType
     */
    inline ObjectType ObjectTypeCast (const std::string objectTypeName)
    { return castFromName<ObjectType>(objectTypeName, ObjectTypeNames()); }

    /**
     * @brief cast an enum of type ObjectType
     * from indeces
     * @param objectTypeIndeces the objectType indeces
     * @return object of class ObjectType
     */
    inline std::vector<ObjectType> ObjectTypeCast (const std::vector<int> objectTypeIndeces)
    { return castFromValues<ObjectType>(objectTypeIndeces); }

    /**
     * @brief cast an enum of type ObjectType
     * from thier names
     * @param objectTypeNames the objectType names
     * @return vector of ObjectType enums
     */
    inline std::vector<ObjectType> ObjectTypeCast (const std::vector<std::string> objectTypeNames)
    { return castFromNames<ObjectType>(objectTypeNames, ObjectTypeNames()); }

    /**
     *  @class Object Object.h "Headers/Object.h"
     *
     *  @brief The class Object
     *
     *  This class is used to handle objects of type <EM> object
     *  </EM>
     */
    class Object {

    protected:

      /// comoving coordinate x
      double m_xx = cbl::par::defaultDouble;

      /// comoving coordinate y
      double m_yy = cbl::par::defaultDouble;

      /// comoving coordinate z
      double m_zz = cbl::par::defaultDouble;

      /// Right Ascension
      double m_ra = cbl::par::defaultDouble;

      /// Declination
      double m_dec = cbl::par::defaultDouble;

      /// redshift
      double m_redshift = cbl::par::defaultDouble;

      /// comoving distance
      double m_dc = cbl::par::defaultDouble;

      /// weight
      double m_weight = 1.;

      /// region (used for jackknife/bootstrap)
      long m_region = cbl::par::defaultLong;

      /// ID
      int m_ID = cbl::par::defaultInt;

      /// observed field
      std::string m_field = cbl::par::defaultString;

      /// displacement along the x-axis
      double m_x_displacement = cbl::par::defaultDouble;

      /// displacement along the y-axis
      double m_y_displacement = cbl::par::defaultDouble;

      /// displacement along the z-axis
      double m_z_displacement = cbl::par::defaultDouble;


    public :

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Object
       */
      Object () = default;

      /**
       *  @brief constructor that uses comoving coordinates
       *
       *  @param coord structure containing the comoving coordinates
       *  {x, y, z}
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @return object of class Object
       */
      Object (const comovingCoordinates coord, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble)
	: m_xx(coord.xx), m_yy(coord.yy), m_zz(coord.zz), m_ra(par::defaultDouble), m_dec(par::defaultDouble), m_redshift(par::defaultDouble), m_dc(par::defaultDouble), m_weight(weight), m_region(region), m_ID(ID), m_field(field), m_x_displacement(x_displacement), m_y_displacement(y_displacement), m_z_displacement(z_displacement) {}

      /**
       *  @brief constructor that uses comoving coordinates and a
       *  cosmological model to estimate the redshift
       *
       *  @param coord structure containing the comoving coordinates
       *  {x, y, z}
       *
       *  @param cosm object of class Cosmology, used to estimate
       *  comoving distances
       *
       *  @param z1_guess minimum prior on the redshift
       *
       *  @param z2_guess maximum prior on the redshift
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @return object of class Object
       */
      Object (const comovingCoordinates coord, const cosmology::Cosmology &cosm, const double z1_guess=0., const double z2_guess=10., const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble)
	: m_xx(coord.xx), m_yy(coord.yy), m_zz(coord.zz), m_weight(weight), m_region(region), m_ID(ID), m_field(field), m_x_displacement(x_displacement), m_y_displacement(y_displacement), m_z_displacement(z_displacement)
	{
	  cbl::polar_coord(m_xx, m_yy, m_zz, m_ra, m_dec, m_dc);
	  m_redshift = cosm.Redshift(m_dc, z1_guess, z2_guess);
	}

      /**
       *  @brief constructor that uses observed coordinates in radians
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshift}
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @return object of class Object
       */
      Object (const observedCoordinates coord, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble)
	: m_ra(coord.ra), m_dec(coord.dec), m_redshift(coord.redshift), m_dc(par::defaultDouble), m_weight(weight), m_region(region), m_ID(ID), m_field(field), m_x_displacement(x_displacement), m_y_displacement(y_displacement), m_z_displacement(z_displacement) {}

      /**
       *  @brief constructor that uses observed coordinates in any
       *  angular units
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshift}
       *
       *  @param inputUnits the units of the input coordinates
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @return object of class Object
       */
      Object (const observedCoordinates coord, const CoordinateUnits inputUnits, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble)
	: m_ra(radians(coord.ra, inputUnits)), m_dec(radians(coord.dec, inputUnits)), m_redshift(coord.redshift), m_dc(par::defaultDouble), m_weight(weight), m_region(region), m_ID(ID), m_field(field), m_x_displacement(x_displacement), m_y_displacement(y_displacement), m_z_displacement(z_displacement) {}

      /**
       *  @brief constructor that uses observed coordinates in radians
       *  and a cosmological model to estimate the comoving
       *  coordinates
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshitf}
       *
       *  @param cosm object of class Cosmology, used to estimate
       *  comoving distances
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @return object of class Object
       */
      Object (const observedCoordinates coord, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble)
	: m_ra(coord.ra), m_dec(coord.dec), m_redshift(coord.redshift), m_dc(par::defaultDouble), m_weight(weight), m_region(region), m_ID(ID), m_field(field), m_x_displacement(x_displacement), m_y_displacement(y_displacement), m_z_displacement(z_displacement)
	{
	  m_dc = cosm.D_C(m_redshift);
	  cbl::cartesian_coord(m_ra, m_dec, m_dc, m_xx, m_yy, m_zz);
	}

      /**
       *  @brief constructor that uses observed coordinates and a
       *  cosmological model to estimate the comoving coordinates
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshift}
       *
       *  @param inputUnits the units of the input coordinates
       *
       *  @param cosm object of class Cosmology, used to estimate
       *  comoving distances
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @return object of class Object
       */
      Object (const observedCoordinates coord, const CoordinateUnits inputUnits, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble)
	: m_ra(radians(coord.ra, inputUnits)), m_dec(radians(coord.dec, inputUnits)), m_redshift(coord.redshift), m_dc(par::defaultDouble), m_weight(weight), m_region(region), m_ID(ID), m_field(field), m_x_displacement(x_displacement), m_y_displacement(y_displacement), m_z_displacement(z_displacement)
	{
	  m_dc = cosm.D_C(m_redshift);
	  cbl::cartesian_coord(m_ra, m_dec, m_dc, m_xx, m_yy, m_zz);
	}

      /**
       *  @brief constructor that uses both comoving and observed
       *  coordinates
       *  @param xx comoving coordinate
       *  @param yy comoving coordinate
       *  @param zz comoving coordinate
       *  @param ra Right Ascension
       *  @param dec Declination
       *  @param redshift redshift
       *  @param weight weight
       *  @param region region, used e.g. for jackknife and bootstrap
       *  @param ID the object ID
       *  @param field the field where the object has been observed
       *  @param x_displacement the displacement along the x-axis
       *  @param y_displacement the displacement along the y-axis
       *  @param z_displacement the displacement along the z-axis
       *  @return object of class Object
       */
      Object (const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble)
      : m_xx(xx), m_yy(yy), m_zz(zz), m_ra(ra), m_dec(dec), m_redshift(redshift), m_dc(sqrt(xx*xx+yy*yy+zz*zz)), m_weight(weight), m_region(region), m_ID(ID), m_field(field), m_x_displacement(x_displacement), m_y_displacement(y_displacement), m_z_displacement(z_displacement)
      {}

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Object () = default;

      ///@}


      /**
       *  @name Static factories
       */
      ///@{

      /**
       *  @brief static factory used to construct objects of any type,
       *  providing in input comoving coordinates
       *
       *  @param ObjectType the object type; it can be: GenericObject,
       *  RandomObject, Mock, Halo, Galaxy, Cluster, Void, HostHalo
       *
       *  @return object of a given type
       */
      static std::shared_ptr<Object> Create (const ObjectType ObjectType );

      /**
       *  @brief static factory used to construct objects of any type,
       *  providing in input comoving coordinates
       *
       *  @param ObjectType the object type; it can be: GenericObject,
       *  RandomObject, Mock, Halo, Galaxy, Cluster, Void, HostHalo
       *
       *  @param coord structure containing the comoving coordinates
       *  {x, y, z}
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @return object of a given type
       */
      static std::shared_ptr<Object> Create (const ObjectType ObjectType, const comovingCoordinates coord, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble);

      /**
       *  @brief static factory used to construct objects of any type,
       *  providing in input comoving coordinates and a cosmological
       *  model to estimate the redshift
       *
       *  @param ObjectType the object type; it can be: GenericObject,
       *  RandomObject, Mock, Halo, Galaxy, Cluster, Void, HostHalo
       *
       *  @param coord structure containing the comoving coordinates
       *  {x, y, z}
       *
       *  @param cosm object of class Cosmology, used to estimate
       *  comoving distances
       *
       *  @param z1_guess minimum prior on the redshift
       *
       *  @param z2_guess maximum prior on the redshift
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @return object of a given type
       */
      static std::shared_ptr<Object> Create (const ObjectType ObjectType, const comovingCoordinates coord, const cosmology::Cosmology &cosm, const double z1_guess=0., const double z2_guess=10., const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble);

      /**
       *  @brief static factory used to construct objects of any kind,
       *  providing in input observed coordinates in radians
       *
       *  @param ObjectType the object type; it can be: GenericObject,
       *  RandomObject, Mock, Halo, Galaxy, Cluster, Void, HostHalo
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshitf}
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @return object of a given type
       */
      static std::shared_ptr<Object> Create (const ObjectType ObjectType, const observedCoordinates coord, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble);

      /**
       *  @brief static factory used to construct objects of any kind,
       *  providing in input observed coordinates in any angular units
       *
       *  @param ObjectType the object type; it can be: GenericObject,
       *  RandomObject, Mock, Halo, Galaxy, Cluster, Void, HostHalo
       *
       *  @param coord structure containing the observed coordinates
       *  {ra Right Ascension}
       *
       *  @param inputUnits the units of the input coordinates
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @return object of a given type
       */
      static std::shared_ptr<Object> Create (const ObjectType ObjectType, const observedCoordinates coord, const CoordinateUnits inputUnits, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble);

      /**
       *  @brief static factory used to construct objects of any kind,
       *  providing in input observed coordinates in radians and a
       *  cosmological model to estimate the comoving coordinates
       *
       *  @param ObjectType the object type; it can be: GenericObject,
       *  RandomObject, Mock, Halo, Galaxy, Cluster, Void, HostHalo
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshitf}
       *
       *  @param cosm object of class Cosmology, used to estimate
       *  comoving distances
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @return object of a given type
       */
      static std::shared_ptr<Object> Create (const ObjectType ObjectType, const observedCoordinates coord, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble);

      /**
       *  @brief static factory used to construct objects of any kind,
       *  providing in input observed coordinates and a cosmological
       *  model to estimate the comoving coordinates
       *
       *  @param ObjectType the object type; it can be: GenericObject,
       *  RandomObject, Mock, Halo, Galaxy, Cluster, Void, HostHalo
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshitf}
       *
       *  @param inputUnits the units of the input coordinates
       *
       *  @param cosm object of class Cosmology, used to estimate
       *  comoving distances
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @return object of a given type
       */
      static std::shared_ptr<Object> Create (const ObjectType ObjectType, const observedCoordinates coord, const CoordinateUnits inputUnits, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble);

      /**
       *  @brief static factory used to construct objects of any kind,
       *  providing in input both comoving and observed coordinates
       *
       *  @param ObjectType the object type; it can be: GenericObject,
       *  RandomObject, Mock, Halo, Galaxy, Cluster, Void, HostHalo
       *
       *  @param xx comoving coordinate
       *
       *  @param yy comoving coordinate
       *
       *  @param zz comoving coordinate
       *
       *  @param ra Right Ascension
       *
       *  @param dec Declination
       *
       *  @param redshift redshift
       *
       *  @param weight weight
       *
       *  @param region the object region, used e.g. for jackknife and
       *  bootstrap
       *
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @return object of a given type
       *
       */
      static std::shared_ptr<Object> Create (const ObjectType ObjectType, const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble);

      ///@}


      /**
       *  @name Member functions used to get the protected members
       */
      ///@{

      /**
       *  @brief get the member \e m_xx
       *  @return the coordinate x of the derived object
       */
      double xx () const
      { return (cbl::isSet(m_xx)) ? m_xx : ErrorCBL("the m_xx variable is not defined!", "xx", "Object.h"); }

      /**
       *  @brief get the member \e m_yy
       *  @return the coordinate y of the derived object
       */
      double yy () const
      { return (cbl::isSet(m_yy)) ? m_yy : ErrorCBL("the m_yy variable is not defined!", "yy", "Object.h"); }

      /**
       *  @brief get the member \e m_zz
       *  @return the coordinate z of the derived object
       */
      double zz () const
      { return (cbl::isSet(m_zz)) ? m_zz : ErrorCBL("the m_zz variable is not defined!", "zz", "Object.h"); }

      /**
       *  @brief get the member \e m_dc
       *  @return the comoving distance of the object
       */
      double dc () const
      { return (cbl::isSet(m_dc)) ? m_dc : ErrorCBL("the m_dc variable is not defined!", "dc", "Object.h"); }

      /**
       *  @brief get the member \e m_ra
       *  @return the Right Ascension of the object
       */
      double ra () const
      { return (cbl::isSet(m_ra)) ? m_ra : ErrorCBL("the m_ra variable is not defined!", "ra", "Object.h"); }

      /**
       *  @brief get the member \e m_dec
       *  @return the Declination of the object
       */
      double dec () const
      { return (cbl::isSet(m_dec)) ? m_dec : ErrorCBL("the m_dec variable is not defined!", "dec", "Object.h"); }

      /**
       *  @brief get the member \e m_redshift
       *  @return the redshift of the object
       */
      double redshift () const
      { return (cbl::isSet(m_redshift)) ? m_redshift : ErrorCBL("the m_redshift variable is not defined!", "redshift", "Object.h"); }

      /**
       *  @brief get the member \e m_weight
       *  @return the weight of the object
       */
      double weight () const
      { return (cbl::isSet(m_weight)) ? m_weight : ErrorCBL("the m_region variable is not defined!", "weight", "Object.h"); }

      /**
       *  @brief get the member \e m_region
       *  @return the index of the subRegion in which the object is located
       */
      long region () const
      { return (cbl::isSet(m_region)) ? m_region : ErrorCBL("the m_region variable is not defined!", "region", "Object.h"); }

      /**
       *  @brief get the member \e m_radius
       *  @return the ID of the derived object, or an
       *  error message if the derived object does not have this member
       */
      int ID () const
      { return (cbl::isSet(m_ID)) ? m_ID : cbl::ErrorCBL("the m_ID variable is not defined!", "ID", "Object.h"); }

      /**
       *  @brief get the member \e m_field
       *  @return the field where the object has been observed
       */
      std::string field () const
	{
	  if (!cbl::isSet(m_field)) ErrorCBL("the m_field variable is not defined!", "field", "Object.h");
	  return m_field;
	}

      /**
       *   @brief get the member \e m_x_displacement
       *   @return the displacement along the x axis
       */
      double x_displacement () const
      { return (cbl::isSet(m_x_displacement)) ? m_x_displacement : ErrorCBL("the m_x_displacement variable is not defined!", "x_displacement", "Object.h"); }

      /**
       *   @brief get the member \e m_y_displacement
       *   @return the displacement along the x axis
       */
      double y_displacement () const
      { return (cbl::isSet(m_y_displacement)) ? m_y_displacement : ErrorCBL("the m_y_displacement variable is not defined!", "y_displacement", "Object.h"); }

      /**
       *   @brief get the member \e m_z_displacement
       *   @return the displacement along the z axis
       */
      double z_displacement () const
      { return (cbl::isSet(m_z_displacement)) ? m_z_displacement : ErrorCBL("the m_z_displacement variable is not defined!", "z_displacement", "Object.h"); }

      /**
       *  @brief get the object coordinates
       *  @return a vector containing the object coordinates
       */
      std::vector<double> coords () const
	{
	  if (!cbl::isSet(m_xx) || !cbl::isSet(m_yy) || !cbl::isSet(m_zz))
	    ErrorCBL("one or more of the m_xx, m_yy, m_zz variables is not defined!", "coords", "Object.h");
	  return {m_xx, m_yy, m_zz};
	}

  /**
  *  @brief get the object coordinates
  *  @return a vector containing the object coordinates
  */
  Vector4D eigen_coords () const
  {
    return {m_xx, m_yy, m_zz, 0};
  }

      /**
       *  @brief get the member \e m_vx
       *  @return the peculiar velocity along the x direction of the
       *  derived object, or an error message if the derived object does
       *  not have this member
       */
      virtual double vx () const
      { return cbl::ErrorCBL("", "vx", "Object.h"); }

      /**
       *  @brief get the member \e m_vy
       *  @return the peculiar velocity along the y direction of the
       *  derived object, or an error message if the derived object does
       *  not have this member
       */
      virtual double vy () const
      { return cbl::ErrorCBL("", "vy", "Object.h"); }

      /**
       *  @brief get the member \e m_vz
       *  @return the peculiar velocity along the z direction of the
       *  derived object, or an message if the derived object does not
       *  have this member
       */
      virtual double vz () const
      { return cbl::ErrorCBL("", "vz", "Object.h"); }

      /**
       *  @brief get the member \e m_mass
       *  @return the mass of the derived object, or an error message if
       *  the derived object does not have this member
       */
      virtual double mass () const
      { return cbl::ErrorCBL("", "mass", "Object.h"); }

      /**
       *  @brief get the member \e m_magnitude
       *  @return the magnitude of the derived object, or an error
       *  message if the derived object does not have this member
       */
      virtual double magnitude () const
      { return cbl::ErrorCBL("", "magnitude", "Object.h"); }

      /**
       *  @brief get the private member Galaxy::m_SFR
       *  @return the star formation rate of the galaxy
       */
      virtual double SFR () const
      { return cbl::ErrorCBL("", "SFR", "Object.h"); }

      /**
       *  @brief get the private member Galaxy::m_sSFR
       *  @return the specific star formation rate of the galaxy
       */
      virtual double sSFR () const
      { return cbl::ErrorCBL("", "sSFR", "Object.h"); }

      /**
       *  @brief get the member \e m_richness
       *  @return the richness of the derived object, or an error
       *  message if the derived object does not have this member
       */
      virtual double richness () const
      { return cbl::ErrorCBL("", "richness", "Object.h"); }

      /**
       *  @brief get the member \e m_richness_error
       *  @return the richness_error of the derived object, or an
       *  error message if the derived object does not have this
       *  member
       */
      virtual double richness_error () const
      { return cbl::ErrorCBL("", "richness_error", "Object.h"); }

      /**
       *  @brief get the member \e m_bias
       *  @return the bias of the derived object, or an error
       *  message if the derived object does not have this member
       */
      virtual double bias () const
      { return cbl::ErrorCBL("", "bias", "Object.h"); }

      /**
       *  @brief get the member \e m_generic
       *  @return the generic variable of the derived object, or an
       *  error message if the derived object does not have this member
       */
      virtual double generic () const
      { return cbl::ErrorCBL("", "generic", "Object.h"); }

      /**
       *  @brief get the member \e m_radius
       *  @return the radius of the derived object, or an
       *  error message if the derived object does not have this member
       */
      virtual double radius () const
      { return cbl::ErrorCBL("", "radius", "Object.h"); }

      /**
       *  @brief get the member \e m_radius
       *  @return the density contrast of the derived object, or an
       *  error message if the derived object does not have this member
       */
      virtual double densityContrast () const
      { return cbl::ErrorCBL("", "densityContrast", "Object.h"); }

      /**
       *  @brief get the member \e m_radius
       *  @return the central density of the derived object, or an
       *  error message if the derived object does not have this member
       */
      virtual double centralDensity () const
      { return cbl::ErrorCBL("", "centralDensity", "Object.h"); }

      /**
       *  @brief get the private member \e m_mass_estimate
       *  @return the mass estimate of the group, or an
       *  error message if the derived object does not have this member
       */
      virtual double mass_estimate () const
      { return cbl::ErrorCBL("", "mass_estimate", "Object.h"); }

      /**
       *  @brief get the private member \e m_radius_estimate
       *  @return the radius estimate of the group, or an
       *  error message if the derived object does not have this member
       */
      virtual double radius_estimate () const
      { return cbl::ErrorCBL("", "radius_estimate", "Object.h"); }

      /**
       *  @brief get the private member \e m_veldisp_estimate
       *  @return the velocity dispersion estimate of the group, or an
       *  error message if the derived object does not have this member
       */
      virtual double veldisp_estimate () const
      { return cbl::ErrorCBL("", "veldisp_estimate", "Object.h"); }

      /**
       *  @brief get the private member \e m_xcm
       *  @return the x-axis coordinate of the centre of mass, or an
       *  error message if the derived object does not have this member
       */
      virtual double xcm () const
      { return cbl::ErrorCBL("", "xcm", "Object.h"); }

      /**
       *  @brief get the private member \e m_ycm
       *  @return the y-axis coordinate of the centre of mass, or an
       *  error message if the derived object does not have this member
       */
      virtual double ycm () const
      { return cbl::ErrorCBL("", "ycm", "Object.h"); }

      /**
       *  @brief get the private member \e m_zcm
       *  @return the z-axis coordinate of the centre of mass, or an
       *  error message if the derived object does not have this member
       */
      virtual double zcm () const
      { return cbl::ErrorCBL("", "zcm", "Object.h"); }

      /**
       *  @brief get the private member \e m_spin_x
       *  @return the x-axis component of the spin, or an
       *  error message if the derived object does not have this member
       */
      virtual double spin_x () const
      { return cbl::ErrorCBL("", "spin_x", "Object.h"); }

      /**
       *  @brief get the private member \e m_spin_y
       *  @return the y-axis component of the spin, or an
       *  error message if the derived object does not have this member
       */
      virtual double spin_y () const
      { return cbl::ErrorCBL("", "spin_y", "Object.h"); }

      /**
       *  @brief get the private member \e m_spin_z
       *  @return the z-axis component of the spin, or an
       *  error message if the derived object does not have this member
       */
      virtual double spin_z () const
      { return cbl::ErrorCBL("", "spin_z", "Object.h"); }

      /**
       *  @brief get the private member \e m_veldisp
       *  @return the velocity dispersion of the sub-group, or an
       *  error message if the derived object does not have this member
       */
      virtual double veldisp () const
      { return cbl::ErrorCBL("", "veldisp", "Object.h"); }

      /**
       *  @brief get the private member \e m_vmax
       *  @return the maximum total velocity of the sub-group, or an
       *  error message if the derived object does not have this member
       */
      virtual double vmax () const
      { return cbl::ErrorCBL("", "vmax", "Object.h"); }

      /**
       *  @brief get the private member \e m_vmax_rad
       *  @return the maximum radial velocity of the sub-group, or an
       *  error message if the derived object does not have this member
       */
      virtual double vmax_rad () const
      { return cbl::ErrorCBL("", "vmax_rad", "Object.h"); }

      /**
       *  @brief get the private member \e m_tot_mass
       *  @return the total mass of the host halo (sum over all contributions), or an
       *  error message if the derived object does not have this member
       */
      virtual double tot_mass () const
      { return cbl::ErrorCBL("", "tot_mass", "Object.h"); }

      /**
       *  @brief get the private member \e m_parent
       *  @return the id of the parent group, or an
       *  error message if the derived object does not have this member
       */
      virtual int parent () const { return cbl::ErrorCBL("", "parent", "Object.h"); }

      /**
       *  @brief get the private member \e m_nsub
       *  @return the number of sub-groups in the group, or an
       *  error message if the derived object does not have this member
       */
      virtual int nsub () const
      { return cbl::ErrorCBL("", "nsub", "Object.h"); }

      /**
       *  @brief get the member \e m_satellites
       *  @return a vector of pointers to satellite objects of the derived object, or an
       *  error message if the derived object does not have this member
       */
      virtual std::vector<std::shared_ptr<Object>> satellites () const
      { cbl::ErrorCBL("", "satellites", "Object.h"); return {}; }

      ///@}


      /**
       *  @name Member functions used to set the protected members
       */
      ///@{

      /**
       *  @brief set the member \e m_xx
       *  @param xx the coordinate x of the object
       *  @return none
       */
      void set_xx (const double xx)
      { m_xx = xx; }

      /**
       *  @brief set the member \e m_yy
       *  @param yy the coordinate y of the object
       *  @return none
       */
      void set_yy (const double yy)
      { m_yy = yy; }

      /**
       *  @brief set the member \e m_zz
       *  @param zz the coordinate z of the object
       *  @return none
       */
      void set_zz (const double zz)
      { m_zz = zz; }

      /**
       *  @brief set the member \e m_ra, updating the
       *  comoving coordinates accordingly (if already set)
       *  @param ra the Right Ascension of the object
       *  @param inputUnits the units of the input coordinates
       *  @return none
       */
      void set_ra (const double ra, const CoordinateUnits inputUnits=CoordinateUnits::_radians_)
      {
	m_ra = radians(ra, inputUnits);
	if (m_dc>par::defaultDouble) cbl::cartesian_coord(m_ra, m_dec, m_dc, m_xx, m_yy, m_zz);
      }

      /**
       *  @brief set the member \e m_dec, updating the
       *  comoving coordinates accordingly (if already set)
       *  @param dec the Declination of the object
       *  @param inputUnits the units of the input coordinates
       *  @return none
       */
      void set_dec (const double dec, const CoordinateUnits inputUnits=CoordinateUnits::_radians_)
      {
	m_dec = radians(dec, inputUnits);
	if (m_dc>par::defaultDouble) cbl::cartesian_coord(m_ra, m_dec, m_dc, m_xx, m_yy, m_zz);
      }

      /**
       *  @brief set the member \e m_redshift, updating
       *  the comoving coordinates accordingly (if already set)
       *
       *  @param redshift the redshift of the object
       *
       *  @param cosmology object of class Cosmology, used to estimate
       *  the comoving distance from the given redshift
       *
       *  @return none
       */
      void set_redshift (const double redshift, const cosmology::Cosmology cosmology)
      {
	m_redshift = redshift;
	m_dc = cosmology.D_C(m_redshift);
	cbl::cartesian_coord(m_ra, m_dec, m_dc, m_xx, m_yy, m_zz);
      }

      /**
       *  @brief set the member \e m_dc, updating the
       *  comoving coordinates accordingly
       *  @param dc the comoving distance of the object
       *  @return none
       */
      void set_dc (const double dc)
      {
	m_dc = dc;
	cbl::cartesian_coord(m_ra, m_dec, m_dc, m_xx, m_yy, m_zz);
      }

      /**
       *  @brief set the member \e m_weight
       *  @param weight the weight of the object
       *  @return none
       */
      void set_weight (const double weight)
      { m_weight = weight; }

      /**
       *  @brief set the member \e m_region
       *  @param region the index of the subRegion in which the object
       *  is located
       *  @return none
       */
      void set_region (const long region)
      { if (region<0) ErrorCBL("Error in Object.h: region must be >0 !", "", "Object.h"); m_region = region; }

      /**
       *  @brief set the member \e m_ID
       *  @param ID the ID
       *  @return none
       */
      void set_ID (const int ID)
      { m_ID = ID; }

      /**
       *  @brief set the member \e m_field
       *  @param field the field were the object has been observed
       *  @return none
       */
      void set_field (const std::string field)
      { m_field = field; }

      /**
       *  @brief set the member \e m_x_displacement
       *  @param x_displacement the displacement (in Mpc) of the x coordinate
       *  @return none
       */
      void set_x_displacement (const double x_displacement)
      { m_x_displacement = x_displacement; }

      /**
       *  @brief set the member \e m_y_displacement
       *  @param y_displacement the displacement (in Mpc) of the y coordinate
       *  @return none
       */
      void set_y_displacement (const double y_displacement)
      { m_y_displacement = y_displacement; }

      /**
       *  @brief set the member \e m_z_displacement
       *  @param z_displacement the displacement (in Mpc) of the x coordinate
       *  @return none
       */
      void set_z_displacement (const double z_displacement)
      { m_z_displacement = z_displacement; }

      /**
       *  @brief set the member \e m_vx
       *  @param vx the peculiar velocity along the x direction
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_vx (const double vx)
      { (void)vx; cbl::ErrorCBL("", "set_vx", "Object.h"); }

      /**
       *  @brief set the member \e m_vy
       *  @param vy the peculiar velocity along the y direction
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_vy (const double vy)
      { (void)vy; cbl::ErrorCBL("", "set_vy", "Object.h"); }

      /**
       *  @brief set the member \e m_vz
       *  @param vz the peculiar velocity along the z direction
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_vz (const double vz)
      { (void)vz; cbl::ErrorCBL("", "set_vz", "Object.h"); }

      /**
       *  @brief set the member \e m_mass
       *  @param mass the mass
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_mass (const double mass)
      { (void)mass; cbl::ErrorCBL("", "set_mass", "Object.h"); }

      /**
       *  @brief set the member \e m_magnitude
       *  @param magnitude the magnitude
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_magnitude (const double magnitude)
      { (void)magnitude; cbl::ErrorCBL("", "set_magnitude", "Object.h"); }

      /**
       *  @brief set the private member Galaxy::m_SFR
       *  @param SFR the star formation rate of the galaxy
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_SFR (const double SFR)
      { (void)SFR; cbl::ErrorCBL("", "set_SFR", "Object.h"); }

      /**
       *  @brief set the private member Galaxy::m_sSFR
       *  @param sSFR the specific star formation rate of the galaxy
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_sSFR (const double sSFR)
      { (void)sSFR; cbl::ErrorCBL("", "set_sSFR", "Object.h"); }

      /**
       *  @brief set the member \e m_richness
       *  @param richness the richness
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_richness (const double richness)
      { (void)richness; cbl::ErrorCBL("", "richness", "Object.h"); }

      /**
       *  @brief set the member \e m_richness_error
       *  @param richness_error the richness error
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_richness_error (const double richness_error)
      { (void)richness_error; cbl::ErrorCBL("", "set_richness_error", "Object.h"); }

      /**
       *  @brief set the member \e m_bias
       *  @param bias the bias
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_bias (const double bias)
      { (void)bias; cbl::ErrorCBL("", "set_bias", "Object.h"); }

      /**
       *  @brief set the member \e m_generic
       *  @param generic the generic variable
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_generic (const double generic)
      { (void)generic; cbl::ErrorCBL("", "set_generic", "Object.h"); }

      /**
       *  @brief set the member \e m_radius
       *  @param radius the radius
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_radius (const double radius)
      { (void)radius; cbl::ErrorCBL("", "set_radius", "Object.h"); }

      /**
       *  @brief set the member \e m_densityContrast
       *  @param densityContrast the density contrast
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_densityContrast (const double densityContrast)
      { (void)densityContrast; cbl::ErrorCBL("", "set_densityContrast", "Object.h"); }

      /**
       *  @brief set the member \e m_centralDensity
       *  @param centralDensity the central density
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_centralDensity (const double centralDensity)
      { (void)centralDensity; cbl::ErrorCBL("", "set_centralDensity", "Object.h"); }

      /**
       *  @brief set the private member \e m_mass_estimate
       *  @param mass_estimate the mass estimate of the group
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_mass_estimate (const double mass_estimate)
      { (void)mass_estimate; cbl::ErrorCBL("", "set_mass_estimate", "Object.h"); }

      /**
       *  @brief set the private member \e m_radius_estimate
       *  @param radius_estimate the radius estimate of the group
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_radius_estimate (const double radius_estimate)
      { (void)radius_estimate; cbl::ErrorCBL("", "set_radius_estimate", "Object.h"); }


      /**
       *  @brief set the private member \e m_veldisp_estimate
       *  @param veldisp_estimate the velocity dispersion estimate of the group
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_veldisp_estimate (const double veldisp_estimate)
      { (void)veldisp_estimate; cbl::ErrorCBL("", "set_veldisp_estimate", "Object.h"); }

      /**
       *  @brief set the private member \e m_xcm
       *  @param xcm the x-axis coordinate of the centre of mass
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_xcm (const double xcm)
      { (void)xcm; cbl::ErrorCBL("", "set_xcm", "Object.h"); }

      /**
       *  @brief set the private member \e m_ycm
       *  @param ycm the y-axis coordinate of the centre of mass
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_ycm (const double ycm)
      { (void)ycm; cbl::ErrorCBL("", "set_ycm", "Object.h"); }

      /**
       *  @brief set the private member \e m_zcm
       *  @param zcm the z-axis coordinate of the centre of mass
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_zcm (const double zcm)
      { (void)zcm; cbl::ErrorCBL("", "set_zcm", "Object.h"); }

      /**
       *  @brief set the private member \e m_spin_x
       *  @param spin_x the x-axis component of the spin
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_spin_x (const double spin_x)
      { (void)spin_x; cbl::ErrorCBL("", "spin_x", "Object.h"); }

      /**
       *  @brief set the private member \e m_spin_y
       *  @param spin_y the y-axis component of the spin
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_spin_y (const double spin_y)
      { (void)spin_y; cbl::ErrorCBL("", "set_spin_y", "Object.h"); }

      /**
       *  @brief set the private member \e m_spin_z
       *  @param spin_z the z-axis component of the spin
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_spin_z (const double spin_z)
      { (void)spin_z; cbl::ErrorCBL("", "set_spin_z", "Object.h"); }

      /**
       *  @brief set the private member \e m_veldisp
       *  @param veldisp the velocity dispersion of the sub-group
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_veldisp (const double veldisp)
      { (void)veldisp; cbl::ErrorCBL("", "set_veldisp", "Object.h"); }

      /**
       *  @brief set the private member \e m_vmax
       *  @param vmax the maximum total velocity of the sub-group
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_vmax (const double vmax)
      { (void)vmax; cbl::ErrorCBL("", "set_vmax", "Object.h"); }

      /**
       *  @brief set the private member \e m_vmax_rad
       *  @param vmax_rad the maximum radial velocity of the sub-group
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_vmax_rad (const double vmax_rad)
      { (void)vmax_rad; cbl::ErrorCBL("", "set_vmax_rad", "Object.h"); }

      /**
       *  @brief set the private member \e m_tot_mass
       *  @param tot_mass the total mass of the parent halo
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_tot_mass (const double tot_mass)
      { (void)tot_mass; cbl::ErrorCBL("", "set_tot_mass", "Object.h"); }

      /**
       *  @brief set the private member \e m_parent
       *  @param parent the id of the parent group
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_parent (const int parent)
      { (void)parent; cbl::ErrorCBL("", "set_parent", "Object.h"); }

      /**
       *  @brief set the private member \e m_nsub
       *  @param nsub the number of sub-groups in the group
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_nsub (const int nsub)
      { (void)nsub; cbl::ErrorCBL("", "set_nsub", "Object.h"); }

      /**
       *  @brief set the private member \em m_satellites
       *  @param satellite the shared pointers to be added to the satellite objects pointer vector
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_satellite (const std::shared_ptr<Object> satellite)
      {	(void)satellite; cbl::ErrorCBL("", "set_satellite", "Object.h"); }

      /**
       *  @brief set the private member \em m_satellites
       *  @param satellites the vector of shared pointers to satellite objects
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_satellites (const std::vector<std::shared_ptr<Object>> satellites)
      {	(void)satellites; cbl::ErrorCBL("", "set_satellites", "Object.h"); }

      ///@}


      /**
       *  @name Member functions used to check if the protected members are set
       */
      ///@{

      /**
       *  @brief check if the member \e m_xx is set
       *
       *  @return true if the coordinate x is set; false otherwise
       */
      bool isSet_xx () const
      { return (cbl::isSet(m_xx)) ? true : false; }

      /**
       *  @brief check if the member \e m_yy is set
       *
       *  @return true if the coordinate y is set; false otherwise
       */
      bool isSet_yy () const
      { return (cbl::isSet(m_yy)) ? true : false; }

      /**
       *  @brief check if the member \e m_zz is set
       *
       *  @return true if the coordinate z is set; false otherwise
       */
      bool isSet_zz ()
      { return (cbl::isSet(m_zz)) ? true : false; }

      /**
       *  @brief check if the member \e m_ra is set
       *
       *  @return true if the coordinate RA is set; false otherwise
       */
      bool isSet_ra ()
      { return (cbl::isSet(m_ra)) ? true : false; }

      /**
       *  @brief check if the member \e m_dec is set
       *
       *  @return true if the coordinate Dec is set; false otherwise
       */
      bool isSet_dec ()
      { return (cbl::isSet(m_dec)) ? true : false; }

      /**
       *  @brief check if the member \e m_redshift is
       *  set
       *
       *  @return true if the redshift is set; false otherwise
       */
      bool isSet_redshift ()
      { return (cbl::isSet(m_redshift)) ? true : false;	}

      /**
       *  @brief check if the member \e m_dc is set
       *
       *  @return true if the comoving distance is set; false
       *  otherwise
       */
      bool isSet_dc ()
      { return (cbl::isSet(m_dc)) ? true : false; }

      /**
       *  @brief check if the member \e m_weight is set
       *
       *  @return true if the weight is set; false otherwise
       */
      bool isSet_weight ()
      { return (cbl::isSet(m_weight)) ? true : false; }

      /**
       *  @brief check if the member \e m_region is set
       *
       *  @return true if the region is set; false otherwise
       */
      bool isSet_region ()
      { return (cbl::isSet(m_region)) ? true : false; }

      /**
       *  @brief check if the member \e m_field is set
       *
       *  @return true if the field is set; false otherwise
       */
      bool isSet_field ()
      { return (cbl::isSet(m_field)) ? true : false; }

      /**
       *  @brief check if the member
       *  \e m_x_displacement is set
       *
       *  @return true if the displacement along the x direction is
       *  set; false otherwise
       */
      bool isSet_x_displacement ()
      { return (cbl::isSet(m_x_displacement)) ? true : false; }

      /**
       *  @brief check if the member
       *  \e m_y_displacement is set
       *
       *  @return true if the displacement along the y direction is
       *  set; false otherwise
       */
      bool isSet_y_displacement ()
      { return (cbl::isSet(m_y_displacement)) ? true : false; }

      /**
       *  @brief check if the member
       *  \e m_z_displacement is set
       *
       *  @return true if the displacement along the z direction is
       *  set; false otherwise
       */
      bool isSet_z_displacement ()
      { return (cbl::isSet(m_z_displacement)) ? true : false; }

      /**
       *  @brief check if the member \e m_vx is set
       *
       *  @return true if the velocity component Vx is set; false
       *  otherwise, or an error message if the derived object does
       *  not have this member
       */
      virtual bool isSet_vx ()
      { return cbl::ErrorCBL("", "isSet_vx", "Object.h"); }

      /**
       *  @brief check if the member \e m_vy is set
       *
       *  @return true if the the velocity component Vy is set; false
       *  otherwise, or an error message if the derived object does
       *  not have this member
       */
      virtual bool isSet_vy ()
      { return cbl::ErrorCBL("", "isSet_vy", "Object.h"); }

      /**
       *  @brief check if the member \e m_vz is set
       *
       *  @return true if the the velocity component Vz is set; false
       *  otherwise, or an error message if the derived object does
       *  not have this member
       */
      virtual bool isSet_vz ()
      { return cbl::ErrorCBL("", "isSet_vz", "Object.h"); }

      /**
       *  @brief check if the member \e m_mass is set
       *
       *  @return true if the mass is set; false otherwise, or an
       *  error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_mass ()
      { return cbl::ErrorCBL("", "isSet_mass", "Object.h"); }

      /**
       *  @brief check if the member \e m_magnitude is set
       *
       *  @return true if the magnitue is set; false otherwise, or an
       *  error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_magnitude ()
      { return cbl::ErrorCBL("", "isSet_magnitude", "Object.h"); }

      /**
       *  @brief check if the member \em m_SFR is set
       *
       *  @return true if the SFR is set; false otherwise, or an error
       *  message if the derived object does not have this member
       */
      virtual bool isSet_SFR ()
      { return cbl::ErrorCBL("", "isSet_SFR", "Object.h"); }

      /**
       *  @brief check if the private member Galaxy::m_sSFR is set
       *
       *  @return true if the sSFR is set; false otherwise, or an
       *  error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_sSFR ()
      { return cbl::ErrorCBL("", "isSet_sSFR", "Object.h"); }

      /**
       *  @brief check if the member \e m_richness is set
       *
       *  @return true if the richness is set; false otherwise, or an
       *  error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_richness ()
      { return cbl::ErrorCBL("", "isSet_richness", "Object.h"); }

      /**
       *  @brief check if the member \e m_richness_error is set
       *
       *  @return true if the error on the richness is set; false
       *  otherwise, or an error message if the derived object does
       *  not have this member
       */
      virtual bool isSet_richness_error ()
      { return cbl::ErrorCBL("", "isSet_richness_error", "Object.h"); }

      /**
       *  @brief check if the member \e m_bias is set
       *
       *  @return true if the bias is set; false otherwise, or an
       *  error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_bias ()
      { return cbl::ErrorCBL("", "isSet_bias", "Object.h"); }

      /**
       *  @brief check if the member \e m_generic is set
       *
       *  @return true if the generic properties is set; false
       *  otherwise, or an error message if the derived object does
       *  not have this member
       */
      virtual bool isSet_generic ()
      { return cbl::ErrorCBL("", "isSet_generic", "Object.h"); }

      /**
       *  @brief check if the member \e m_radius is set
       *
       *  @return true if the radius is set; false otherwise, or an
       *  error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_radius ()
      { return cbl::ErrorCBL("", "isSet_radius", "Object.h"); }

      /**
       *  @brief check if the member \e m_densityContrast is set
       *
       *  @return true if the density contrast is set; false
       *  otherwise, or an error message if the derived object does
       *  not have this member
       */
      virtual bool isSet_densityContrast ()
      { return cbl::ErrorCBL("", "isSet_densityContrast", "Object.h"); }

      /**
       *  @brief check if the member \e m_centralDensity is set
       *
       *  @return true if the central density is set; false otherwise,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_centralDensity ()
      { return cbl::ErrorCBL("", "isSet_centralDensity", "Object.h"); }

      /**
       *  @brief check if the member \e m_mass_estimate is set
       *
       *  @return true if the  is set; false otherwise,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_mass_estimate ()
      { return cbl::ErrorCBL("", "isSet_mass_estimate", "Object.h"); }

      /**
       *  @brief check if the member \e m_radius_estimate is set
       *
       *  @return true if the  is set; false otherwise,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_radius_estimate ()
      { return cbl::ErrorCBL("", "isSet_radius_estimate", "Object.h"); }

      /**
       *  @brief check if the member \e m_veldisp_estimate is set
       *
       *  @return true if the  is set; false otherwise,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_veldisp_estimate ()
      { return cbl::ErrorCBL("", "isSet_veldisp_estimate", "Object.h"); }

      /**
       *  @brief check if the member \e m_xcm is set
       *
       *  @return true if the  is set; false otherwise,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_xcm ()
      { return cbl::ErrorCBL("", "isSet_xcm", "Object.h"); }

      /**
       *  @brief check if the member \e m_ycm is set
       *
       *  @return true if the  is set; false otherwise,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_ycm ()
      { return cbl::ErrorCBL("", "isSet_ycm", "Object.h"); }

      /**
       *  @brief check if the member \e m_zcm is set
       *
       *  @return true if the  is set; false otherwise,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_zcm ()
      { return cbl::ErrorCBL("", "isSet_zcm", "Object.h"); }

      /**
       *  @brief check if the member \e m_spin_x is set
       *
       *  @return true if the  is set; false otherwise,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_spin_x ()
      { return cbl::ErrorCBL("", "isSet_spin_x", "Object.h"); }

      /**
       *  @brief check if the member \e m_spin_y is set
       *
       *  @return true if the  is set; false otherwise,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_spin_y ()
      { return cbl::ErrorCBL("", "isSet_spin_y", "Object.h"); }

      /**
       *  @brief check if the member \e m_spin_z is set
       *
       *  @return true if the  is set; false otherwise,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_spin_z ()
      { return cbl::ErrorCBL("", "isSet_spin_z", "Object.h"); }

      /**
       *  @brief check if the member \e m_veldisp is set
       *
       *  @return true if the  is set; false otherwise,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_veldisp ()
      { return cbl::ErrorCBL("", "isSet_veldisp", "Object.h"); }

      /**
       *  @brief check if the member \e m_vmax is set
       *
       *  @return true if the  is set; false otherwise,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_vmax ()
      { return cbl::ErrorCBL("", "isSet_vmax", "Object.h"); }

      /**
       *  @brief check if the member \e m_vmax_rad is set
       *
       *  @return true if the  is set; false otherwise,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_vmax_rad ()
      { return cbl::ErrorCBL("", "isSet_vmax_rad", "Object.h"); }

      /**
       *  @brief check if the member \e m_tot_mass is set
       *
       *  @return true if the  is set; false otherwise,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_tot_mass ()
      { return cbl::ErrorCBL("", "isSet_tot_mass", "Object.h"); }

      /**
       *  @brief check if the member \e m_parent is set
       *
       *  @return true if the  is set; false otherwise,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_parent ()
      { return cbl::ErrorCBL("", "isSet_parent", "Object.h"); }

      /**
       *  @brief check if the member \e m_nsub is set
       *
       *  @return true if the  is set; false otherwise,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_nsub ()
      { return cbl::ErrorCBL("", "isSet_nsub", "Object.h"); }

      /**
       *  @brief check if the member \e m_ID is set
       *
       *  @return true if the object ID is set; false otherwise, or an
       *  error message if the derived object does not have this
       *  member
       */
      virtual bool isSet_ID ()
      { return cbl::ErrorCBL("", "isSet_ID", "Object.h"); }

      ///@}

    };
  }
}

#endif
