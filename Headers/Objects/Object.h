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
 *  @file Headers/Objects/Object.h
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
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __OBJECT__
#define __OBJECT__ 


// ============================================================================================


namespace cosmobl {

  namespace catalogue {
    
    /**
     * @enum ObjType
     * @brief the object types
     */
    enum ObjType {
    
      /// generic object
      _GenericObject_,
    
      /// random object
      _RandomObject_,
    
      /// mock object
      _Mock_,
    
      /// dark matter halo
      _Halo_,
    
      /// galaxy
      _Galaxy_,
    
      /// galaxy cluster
      _Cluster_,
    
      /// cosmic void
      _Void_
    
    };
  
    /**
     *  @class Object Object.h "Headers/Lib/Object.h"
     *
     *  @brief The class Object
     *
     *  This class is used to handle objects of type <EM> object
     *  </EM>
     */
    class Object {

    protected:
    
      /// comoving coordinate x
      double m_xx;
    
      /// comoving coordinate y
      double m_yy;

      /// comoving coordinate z
      double m_zz;

      /// Right Ascension 
      double m_ra; 

      /// Declination 
      double m_dec; 

      /// redshift
      double m_redshift;

      /// comoving distance
      double m_dc;

      /// weight
      double m_weight;
  
      /// region
      long m_region;

    
    public :
    
      /**
       *  @name Constructors/destructors
       */
      ///@{
    
      /**
       *  @brief default constructor
       *  @return object of class Object
       */
      Object () 
	: m_xx(par::defaultDouble), m_yy(par::defaultDouble), m_zz(par::defaultDouble), m_ra(par::defaultDouble), m_dec(par::defaultDouble), m_redshift(par::defaultDouble), m_dc(par::defaultDouble), m_weight(1.), m_region(par::defaultInt) {}
    
      /**
       *  @brief constructor that uses comoving coordinates
       *  @param xx comoving coordinate
       *  @param yy comoving coordinate
       *  @param zz comoving coordinate 
       *  @param weight weight
       *  @return object of class GenericObject
       */
      Object (const double xx, const double yy, const double zz, const double weight=1.) 
	: m_xx(xx), m_yy(yy), m_zz(zz), m_ra(par::defaultDouble), m_dec(par::defaultDouble), m_redshift(par::defaultDouble), m_dc(par::defaultDouble), m_weight(weight), m_region(-1) {}
    
      /**
       *  @brief constructor that uses comoving coordinates
       *  @param xx comoving coordinate
       *  @param yy comoving coordinate
       *  @param zz comoving coordinate 
       *  @param cosm object of class Cosmology, used to estimate
       *  comoving distances
       *  @param z1_guess minimum redshift used to search the redshift
       *  @param z2_guess maximum redshift used to search the redshift
       *  @param weight weight
       *  @return object of class GenericObject
       */
      Object (const double xx, const double yy, const double zz, const Cosmology &cosm, const double z1_guess, const double z2_guess, const double weight=1.) 
	: m_xx(xx), m_yy(yy), m_zz(zz), m_dc(-1.), m_weight(weight), m_region(-1)
	{
	  cosmobl::polar_coord(m_xx, m_yy, m_zz, m_ra, m_dec, m_dc);
	  m_redshift = cosm.Redshift(m_dc, z1_guess, z2_guess);
	}
    
      /**
       *  @brief constructor that uses observed coordinates
       *  @param ra Right Ascension
       *  @param dec Declination
       *  @param redshift redshift
       *  @param cosm object of class Cosmology, used to estimate comoving distances
       *  @param weight weight
       *  @return object of class GenericObject
       */
      Object (const double ra, const double dec, const double redshift, const Cosmology &cosm, const double weight=1.) 
	: m_ra(ra), m_dec(dec), m_redshift(redshift), m_dc(-1.), m_weight(weight), m_region(-1)
	{ 
	  m_dc = cosm.D_C(m_redshift); 
	  cosmobl::cartesian_coord(m_ra, m_dec, m_dc, m_xx, m_yy, m_zz);
	}

      /**
       *  @brief constructor that uses both comoving and observed coordinates
       *  @param xx comoving coordinate
       *  @param yy comoving coordinate
       *  @param zz comoving coordinate 
       *  @param ra Right Ascension
       *  @param dec Declination
       *  @param redshift redshift
       *  @param weight weight
       *  @return object of class GenericObject
       */
      Object (const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight=1.) 
	: m_xx(xx), m_yy(yy), m_zz(zz), m_ra(ra), m_dec(dec), m_redshift(redshift), m_dc(sqrt(xx*xx+yy*yy+zz*zz)), m_weight(weight), m_region(-1) {}
    
      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Object () {}

      ///@}

    
      /**
       *  @name Static factories
       */
      ///@{
    
      /**
       * @brief static factory used to construct objects of any type 
       * @param type the object type; it can be: GenericObject,
       * RandomObject, Mock, Halo, Galaxy, Cluster, Void
       * @param xx comoving coordinate
       * @param yy comoving coordinate
       * @param zz comoving coordinate
       * @param weight weight
       * @return object of a given type
       */
      static shared_ptr<Object> Create (const ObjType type, const double xx, const double yy, const double zz, const double weight=1);

      /**
       * @brief static factory used to construct objects of any kind 
       * @param type the object type; it can be: GenericObject,
       * RandomObject, Mock, Halo, Galaxy, Cluster, Void
       * @param ra Right Ascension
       * @param dec Declination
       * @param redshift redshift
       * @param cosm object of class Cosmology, used to estimate comoving distances 
       * @param weight weight
       * @return object of a given type
       */
      static shared_ptr<Object> Create (const ObjType type, const double ra, const double dec, const double redshift, const Cosmology &cosm, const double weight=1);
    
      ///@}

    
      /**
       *  @name Member functions used to get the protected members
       */
      ///@{
    
      /**
       *  @brief get the member \e m_xx
       *  @return the coordinate x of the derived object
       */
      double xx () const { return m_xx; }

      /**
       *  @brief get the member \e m_yy
       *  @return the coordinate y of the derived object
       */
      double yy () const { return m_yy; }

      /**
       *  @brief get the member \e m_zz
       *  @return the coordinate z of the derived object
       */
      double zz () const { return m_zz; }

      /**
       *  @brief get the protected member GenericObject::m_dc
       *  @return the comoving distance of the object
       */
      double dc () const { return m_dc; }
    
      /**
       *  @brief get the protected member GenericObject::m_ra
       *  @return the Right Ascension of the object
       */
      double ra () const { return m_ra; } 
    
      /**
       *  @brief get the protected member GenericObject::m_dec
       *  @return the Declination of the object
       */
      double dec () const { return m_dec; }
    
      /**
       *  @brief get the protected member GenericObject::m_redshift
       *  @return the redshift of the object
       */
      double redshift () const { return m_redshift; }
    
      /**
       *  @brief get the protected member GenericObject::m_weight
       *  @return the weight of the object
       */
      double weight () const { return m_weight; }

      /**
       *  @brief get the protected member GenericObject::m_region
       *  @return the index of the subRegion in which the object is located
       */
      long region () const { return m_region; }

      /**
       *  @brief get the object coordinates
       *  @return a vector containing the object coordinates
       */
      vector<double> coords () const { return {m_xx, m_yy, m_zz}; }
    
      /**
       *  @brief get the member \e m_vx
       *  @return the peculiar velocity along the x direction of the
       *  derived object, or an error message if the derived object does
       *  not have this member
       */
      virtual double vx () const { cosmobl::ErrorMsg("Error in vx() of Objech.h!"); return 0; }
    
      /**
       *  @brief get the member \e m_vy
       *  @return the peculiar velocity along the y direction of the
       *  derived object, or an error message if the derived object does
       *  not have this member
       */
      virtual double vy () const { cosmobl::ErrorMsg("Error in vy() of Objech.h!"); return 0; }
    
      /**
       *  @brief get the member \e m_vz
       *  @return the peculiar velocity along the z direction of the
       *  derived object, or an message if the derived object does not
       *  have this member
       */
      virtual double vz () const { cosmobl::ErrorMsg("Error in vz() of Objech.h!"); return 0; }
  
      /**
       *  @brief get the member \e m_mass
       *  @return the mass of the derived object, or an error message if
       *  the derived object does not have this member
       */
      virtual double mass () const { cosmobl::ErrorMsg("Error in mass() of Objech.h!"); return 0; }
    
      /**
       *  @brief get the member \e m_magnitude
       *  @return the magnitude of the derived object, or an error
       *  message if the derived object does not have this member
       */
      virtual double magnitude () const { cosmobl::ErrorMsg("Error in magnitude() of Objech.h!"); return 0; }  
    
      /**
       *  @brief get the member \e m_richness
       *  @return the richness of the derived object, or an error
       *  message if the derived object does not have this member
       */
      virtual double richness () const { cosmobl::ErrorMsg("Error in richness() of Objech.h!"); return 0; }  
    
      /**
       *  @brief get the member \e m_generic
       *  @return the generic variable of the derived object, or an
       *  error message if the derived object does not have this member
       */
      virtual double generic () const { cosmobl::ErrorMsg("Error in generic() of Objech.h!"); return 0; }  
    
      /**
       *  @brief get the member \e m_radius
       *  @return the radius of the derived object, or an
       *  error message if the derived object does not have this member
       */
      virtual double radius () const { cosmobl::ErrorMsg("Error in radius() of Objech.h!"); return 0; }

      ///@}

    
      /**
       *  @name Member functions used to set the protected members
       */
      ///@{ 
    
      /**
       *  @brief set the protected member GenericObject::m_xx
       *  @param xx the coordinate x of the object
       *  @return none
       */
      void set_xx (const double xx=par::defaultDouble) { m_xx = xx; }
 
      /**
       *  @brief set the protected member GenericObject::m_yy
       *  @param yy the coordinate y of the object
       *  @return none
       */
      void set_yy (const double yy=par::defaultDouble) { m_yy = yy; }
    
      /**
       *  @brief set the protected member GenericObject::m_zz
       *  @param zz the coordinate z of the object
       *  @return none
       */
      void set_zz (const double zz=par::defaultDouble) { m_zz = zz; }
    
      /**
       *  @brief set the protected member GenericObject::m_ra
       *  @param ra the Right Ascension of the object
       *  @return none
       */
      void set_ra (const double ra=par::defaultDouble) { m_ra = ra; }
    
      /**
       *  @brief set the protected member GenericObject::m_dec
       *  @param dec the Declination of the object
       *  @return none
       */
      void set_dec (const double dec=par::defaultDouble) { m_dec = dec; }
    
      /**
       *  @brief set the protected member GenericObject::m_redshift
       *  @param redshift the redshift of the object
       *  @return none
       */
      void set_redshift (const double redshift=par::defaultDouble) { m_redshift = redshift; }
    
      /**
       *  @brief set the protected member GenericObject::m_dc
       *  @param dc the comoving distance of the object
       *  @return none
       */
      void set_dc (const double dc=par::defaultDouble) { m_dc = dc; }
    
      /**
       *  @brief set the protected member GenericObject::m_weight
       *  @param weight the weight of the object
       *  @return none
       */
      void set_weight (const double weight=par::defaultDouble) { m_weight = weight; }

      /**
       *  @brief set the protected member GenericObject::m_region
       *  @param region the index of the subRegion in which the object is located
       *  @return none
       */
      void set_region (const long region=par::defaultLong) { m_region = region; }
    
      /**
       *  @brief set the member \e m_vx
       *  @param vx the peculiar velocity along the x direction
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_vx (const double vx=par::defaultDouble) { cosmobl::ErrorMsg("Error in set_vx() of Objech.h!"); }
  
      /**
       *  @brief set the member \e m_vy
       *  @param vy the peculiar velocity along the y direction 
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_vy (const double vy=par::defaultDouble) { cosmobl::ErrorMsg("Error in set_vy() of Objech.h!"); }
    
      /**
       *  @brief set the member \e m_vz
       *  @param vz the peculiar velocity along the z direction
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_vz (const double vz=par::defaultDouble) { cosmobl::ErrorMsg("Error in set_vz() of Objech.h!"); }

      /**
       *  @brief set the member \e m_mass
       *  @param mass the mass
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_mass (const double mass=par::defaultDouble) { cosmobl::ErrorMsg("Error in set_mass() of Objech.h!"); }
    
      /**
       *  @brief set the member \e m_magnitude
       *  @param magnitude the magnitude
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_magnitude (const double magnitude=par::defaultDouble) { cosmobl::ErrorMsg("Error in set_magnitude() of Objech.h!"); }  
    
      /**
       *  @brief set the member \e m_richness
       *  @param richness the richness 
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_richness (const double richness=par::defaultDouble) { cosmobl::ErrorMsg("Error in set_richness() of Objech.h!"); }  
    
      /**
       *  @brief set the member \e m_generic
       *  @param generic the generic variable
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_generic (const double generic=par::defaultDouble) { cosmobl::ErrorMsg("Error in set_generic() of Objech.h!"); }  
    
      /**
       *  @brief set the member \e m_radius
       *  @param radius the radius
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_radius (const double radius=par::defaultDouble) { cosmobl::ErrorMsg("Error in set_radius() of Objech.h!"); }

      ///@}
    
    };
  }
}

#endif
