/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Tommaso Ronconi	    *
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
 *  @file Headers/Void.h
 *
 *  @brief The class Void
 *
 *  This file defines the interface of the class Void, used to
 *  handle objects of type <EM> void </EM>
 *
 *  @authors Federico Marulli, Tommaso Ronconi
 *
 *  @authors federico.marulli3@unibo.it, tommaso.ronconi@studio.unibo.it
 */

#ifndef __VOID__
#define __VOID__


// ===================================================================================================

namespace cbl {

  namespace catalogue {
    
    /**
     *  @class Void Void.h "Headers/Void.h"
     *
     *  @brief The class Void
     *
     *  This class is used to handle objects of type <EM> void
     *  </EM>
     */
    class Void : public Object {

    private:

      /// void radius [Mpc/h]
      double m_radius;

      /// void density contrast (\f${\rho_v}/{\rho_m}\f$) 
      double m_densityContrast;

      /// void central density
      double m_centralDensity;
      
      /// generic variable
      double m_generic;

    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{
      
      /**
       * @brief function that allows copying private variables of the class 
       * when an object of class Catalogue is copied
       * 
       * @return a shared pointer to the Object
       *
       */
      std::shared_ptr<Object> getShared() {
        return std::make_shared<Void>(*this);
      }
      
      /**
       *  @brief default constructor
       */
      Void () 
      	: Object(), m_radius(par::defaultDouble), m_densityContrast(par::defaultDouble), m_centralDensity(par::defaultDouble) {}
      
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
       *  @param radius radius of the sphere with volume equivalent to that of the void
       *
       *  @param densityContrast ratio between the central density and the density at the border of the void
       *
       *  @param centralDensity central density
       *
       *  
       */
      Void (const comovingCoordinates coord, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double radius=par::defaultDouble, const double densityContrast=par::defaultDouble, const double centralDensity=par::defaultDouble) 
	: Object(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement), m_radius(radius), m_densityContrast(densityContrast), m_centralDensity(centralDensity) {}

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
       *  @param radius radius of the sphere with volume equivalent to that of the void
       *
       *  @param densityContrast ratio between the central density and the density at the border of the void
       *
       *  @param centralDensity central density
       *
       *  
       */
      Void (const comovingCoordinates coord, const cosmology::Cosmology &cosm, const double z1_guess=0., const double z2_guess=10., const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double radius=par::defaultDouble, const double densityContrast=par::defaultDouble, const double centralDensity=par::defaultDouble) 
	: Object(coord, cosm, z1_guess, z2_guess, weight, region, ID, field, x_displacement, y_displacement, z_displacement), m_radius(radius), m_densityContrast(densityContrast), m_centralDensity(centralDensity) {}

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
       *  @param radius radius of the sphere with volume equivalent to that of the void
       *
       *  @param densityContrast ratio between the central density and the density at the border of the void
       *
       *  @param centralDensity central density
       *
       *  
       */
      Void (const observedCoordinates coord, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double radius=par::defaultDouble, const double densityContrast=par::defaultDouble, const double centralDensity=par::defaultDouble) 
	: Object(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement), m_radius(radius), m_densityContrast(densityContrast), m_centralDensity(centralDensity) {}
      
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
       *  @param radius radius of the sphere with volume equivalent to that of the void
       *
       *  @param densityContrast ratio between the central density and the density at the border of the void
       *
       *  @param centralDensity central density
       *
       *  
       */
      Void (const observedCoordinates coord, const CoordinateUnits inputUnits, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double radius=par::defaultDouble, const double densityContrast=par::defaultDouble, const double centralDensity=par::defaultDouble) 
	: Object(coord, inputUnits, weight, region, ID, field, x_displacement, y_displacement, z_displacement), m_radius(radius), m_densityContrast(densityContrast), m_centralDensity(centralDensity) {}
      
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
       *  @param radius radius of the sphere with volume equivalent to that of the void
       *
       *  @param densityContrast ratio between the central density and the density at the border of the void
       *
       *  @param centralDensity central density
       *
       *  
       */
      Void (const observedCoordinates coord, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double radius=par::defaultDouble, const double densityContrast=par::defaultDouble, const double centralDensity=par::defaultDouble) 
	: Object(coord, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement), m_radius(radius), m_densityContrast(densityContrast), m_centralDensity(centralDensity) {}

      /**
       *  @brief constructor that uses observed coordinates and a
       *  cosmological model to estimate the comoving coordinates
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshift}
       *
       *  @param inputUnits the units of the input coordinates
       *
       *  @param cosm object of class Cosmology, used to estimate comoving distances
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
       *  @param radius radius of the sphere with volume equivalent to that of the void
       *
       *  @param densityContrast ratio between the central density and the density at the border of the void
       *
       *  @param centralDensity central density
       *
       *  
       */
      Void (const observedCoordinates coord, const CoordinateUnits inputUnits, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double radius=par::defaultDouble, const double densityContrast=par::defaultDouble, const double centralDensity=par::defaultDouble) 
	: Object(coord, inputUnits, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement), m_radius(radius), m_densityContrast(densityContrast), m_centralDensity(centralDensity) {}

      /**
       *  @brief constructor that uses both comoving and observed coordinates
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
       *  @param radius radius 
       *
       *  @param densityContrast density contrast
       *
       *  @param centralDensity central density
       *
       *  
       */
      Void (const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double radius=par::defaultDouble, const double densityContrast=par::defaultDouble, const double centralDensity=par::defaultDouble) 
	: Object(xx, yy, zz, ra, dec, redshift, weight, region, ID, field, x_displacement, y_displacement, z_displacement), m_radius(radius), m_densityContrast(densityContrast), m_centralDensity(centralDensity) {}
      
      /**
       *  @brief default destructor
       */
      ~Void () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to get the private members 
       */
      ///@{
    
      /**
       *  @brief get the private member \e m_radius
       *  @return the radius of the sphere with volume equivalent to that of the void
       */
      double radius () const override
      { return m_radius; }

      /**
       *  @brief get the private member \e m_densityContrast
       *  @return the ratio between the central density of the void and density at its border
       */
      double densityContrast () const
      { return m_densityContrast; }

      /**
       *  @brief get the private member \e m_centralDensity
       *  @return the central density of the void
       */
      double centralDensity () const
      { return m_centralDensity; }
      
      /**
       *  @brief get the private member \e m_generic
       *  @return the generic variable of the void object
       */
      double generic () const override { return m_generic; }

      ///@}
  
  
      /**
       *  @name Member functions used to set the private members 
       */
      ///@{
    	
      /**
       *  @brief get the private member \e m_radius
       *
       *  @param radius the radius of the sphere with volume
       *  equivalent to that of the void
       */
      void set_radius (const double radius=par::defaultDouble) override
      { m_radius = radius; }
	
      /**
       *  @brief get the private member \e m_densityContrast
       *
       *  @param densityContrast ratio between the central density of
       *  the void and the density at its border
       */
      void set_densityContrast (const double densityContrast=par::defaultDouble)
      { m_densityContrast = densityContrast; }
	
      /**
       *  @brief get the private member \e m_centralDensity
       *  @param centralDensity the central density of the void
       */
      void set_centralDensity (const double centralDensity=par::defaultDouble)
      { m_centralDensity = centralDensity; }
      
      /**
       *  @brief set the private member \e m_generic
       *  @param generic the generic variable of the mock object
       */
      void set_generic (const double generic=par::defaultDouble) override { m_generic = generic; }
    
      ///@}
      
      /**
       *  @name Member functions used to check if the private members are set 
       */
      ///@{
    	
      /**
       *  @brief get the private member \e m_radius
       *  
       *  @return true if the radius is set; false otherwise
       */
      bool isSet_radius () override
      { return (cbl::isSet(m_radius)) ? true : false; }
	
      /**
       *  @brief get the private member \e m_densityContrast
       *  
       *  @return true if the density contrast is set; false otherwise
       */
      bool isSet_densityContrast ()
      { return (cbl::isSet(m_densityContrast)) ? true : false; }
	
      /**
       *  @brief get the private member \e m_centralDensity
       *  
       *  @return true if the central density is set; false otherwise
       */
      bool isSet_centralDensity ()
      { return (cbl::isSet(m_centralDensity)) ? true : false; }
      
      /**
       *  @brief set the private member \e m_generic
       *  
       *  @return true if the generic variable is set; false otherwise
       */
      bool isSet_generic () override
      { return (cbl::isSet(m_generic)) ? true : false; }
    
      ///@}
      
    };
  }
}

#endif
