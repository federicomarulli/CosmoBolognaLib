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
 *  @file Headers/Objects/Void.h
 *
 *  @brief The class Void
 *
 *  This file defines the interface of the class Void, used to
 *  handle objects of type <EM> void </EM>
 *
 *  @authors Federico Marulli, Tommaso Ronconi
 *
 *  @authors federico.marulli3@unbo.it, tommaso.ronconi@studio.unibo.it
 */

#ifndef __VOID__
#define __VOID__


// ===================================================================================================

namespace cosmobl {

  namespace catalogue {
    
    /**
     *  @class Void Void.h "Headers/Lib/Void.h"
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

      /// void ID
      int m_ID;

    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{
      
      /**
       *  @brief default constructor
       *  @return object of class Void
       */
      Void () 
      	: Object(), m_radius(par::defaultDouble), m_densityContrast(par::defaultDouble), m_centralDensity(par::defaultDouble), m_ID(par::defaultInt) {}
      
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
       *  @param ID identification number
       *
       *  @return object of class Void
       */
      Void (const comovingCoordinates coord, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double radius=par::defaultDouble, const double densityContrast=par::defaultDouble, const double centralDensity=par::defaultDouble, const int ID=par::defaultInt) 
	: Object(coord, weight, region, field, x_displacement, y_displacement, z_displacement), m_radius(radius), m_densityContrast(densityContrast), m_centralDensity(centralDensity), m_ID(ID) {}

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
       *  @param ID identification number
       *
       *  @return object of class Void
       */
      Void (const comovingCoordinates coord, const cosmology::Cosmology &cosm, const double z1_guess=0., const double z2_guess=10., const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double radius=par::defaultDouble, const double densityContrast=par::defaultDouble, const double centralDensity=par::defaultDouble, const int ID=par::defaultInt) 
	: Object(coord, cosm, z1_guess, z2_guess, weight, region, field, x_displacement, y_displacement, z_displacement), m_radius(radius), m_densityContrast(densityContrast), m_centralDensity(centralDensity), m_ID(ID) {}

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
       *  @param ID identification number
       *
       *  @return object of class Void
       */
      Void (const observedCoordinates coord, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double radius=par::defaultDouble, const double densityContrast=par::defaultDouble, const double centralDensity=par::defaultDouble, const int ID=par::defaultInt) 
	: Object(coord, weight, region, field, x_displacement, y_displacement, z_displacement), m_radius(radius), m_densityContrast(densityContrast), m_centralDensity(centralDensity), m_ID(ID) {}
      
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
       *  @param ID identification number
       *
       *  @return object of class Void
       */
      Void (const observedCoordinates coord, const CoordUnits inputUnits, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double radius=par::defaultDouble, const double densityContrast=par::defaultDouble, const double centralDensity=par::defaultDouble, const int ID=par::defaultInt) 
	: Object(coord, inputUnits, weight, region, field, x_displacement, y_displacement, z_displacement), m_radius(radius), m_densityContrast(densityContrast), m_centralDensity(centralDensity), m_ID(ID) {}
      
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
       *  @param ID identification number
       *
       *  @return object of class Void
       */
      Void (const observedCoordinates coord, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double radius=par::defaultDouble, const double densityContrast=par::defaultDouble, const double centralDensity=par::defaultDouble, const int ID=par::defaultInt) 
	: Object(coord, cosm, weight, region, field, x_displacement, y_displacement, z_displacement), m_radius(radius), m_densityContrast(densityContrast), m_centralDensity(centralDensity), m_ID(ID) {}

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
       *  @param ID identification number
       *
       *  @return object of class Void
       */
      Void (const observedCoordinates coord, const CoordUnits inputUnits, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double radius=par::defaultDouble, const double densityContrast=par::defaultDouble, const double centralDensity=par::defaultDouble, const int ID=par::defaultInt) 
	: Object(coord, inputUnits, cosm, weight, region, field, x_displacement, y_displacement, z_displacement), m_radius(radius), m_densityContrast(densityContrast), m_centralDensity(centralDensity), m_ID(ID) {}

      /**
       *  @brief constructor that uses both comoving and observed coordinates
       *  @param xx comoving coordinate
       *  @param yy comoving coordinate
       *  @param zz comoving coordinate 
       *  @param ra Right Ascension
       *  @param dec Declination
       *  @param redshift redshift
       *  @param weight weight
       *  @param region region, used e.g. for jackknife and bootstrap
       *  @param field the field where the object has been observed
       *  @param x_displacement the displacement along the x-axis
       *  @param y_displacement the displacement along the y-axis
       *  @param z_displacement the displacement along the z-axis
       *  @param radius radius 
       *  @param densityContrast density contrast
       *  @param centralDensity central density
       *  @param ID identification number
       *
       *  @return object of class Void
       */
      Void (const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double radius=par::defaultDouble, const double densityContrast=par::defaultDouble, const double centralDensity=par::defaultDouble, const int ID=par::defaultInt) 
	: Object(xx, yy, zz, ra, dec, redshift, weight, region, field, x_displacement, y_displacement, z_displacement), m_radius(radius), m_densityContrast(densityContrast), m_centralDensity(centralDensity), m_ID(ID) {}
      
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Void () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to get the private members 
       */
      ///@{
    
      /**
       *  @brief get the private member Void::m_radius
       *  @return the radius of the sphere with volume equivalent to that of the void
       */
      double radius () const override { return m_radius; }

      /**
       *  @brief get the private member Void::m_densityContrast
       *  @return the ratio between the central density of the void and density at its border
       */
      double densityContrast () const { return m_densityContrast; }

      /**
       *  @brief get the private member Void::m_centralDensity
       *  @return the central density of the void
       */
      double centralDensity () const { return m_centralDensity; }

      /**
       *  @brief get the private member Void::m_ID
       *  @return the identification number of the void
       */
      int ID () const { return m_ID; }

      ///@}
  
  
      /**
       *  @name Member functions used to set the private members 
       */
      ///@{
    	
      /**
       *  @brief get the private member Void::m_radius
       *  @param radius the radius of the sphere with volume equivalent to that of the void
       *  @return none
       */
      void set_radius (const double radius=par::defaultDouble) override { m_radius = radius; }
	
      /**
       *  @brief get the private member Void::m_densityContrast
       *  @param densityContrast ratio between the central density of the void and the density at its border
       *  @return none
       */
      void set_densityContrast (const double densityContrast=par::defaultDouble) { m_densityContrast = densityContrast; }
	
      /**
       *  @brief get the private member Void::m_centralDensity
       *  @param centralDensity the central density of the void
       *  @return none
       */
      void set_centralDensity (const double centralDensity=par::defaultDouble) { m_centralDensity = centralDensity; }
	
      /**
       *  @brief get the private member Void::m_ID
       *  @param ID the identification number of the void
       *  @return none
       */
      void set_ID (const int ID=par::defaultInt) { m_ID = ID; } 
    
      ///@}
    
    };
  }
}

#endif
