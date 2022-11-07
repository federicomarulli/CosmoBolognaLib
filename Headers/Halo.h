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
 *  @file Headers/Halo.h
 *
 *  @brief The class Halo 
 *
 *  This file defines the interface of the class Halo, used to handle
 *  objects of type <EM> dark matter halo </EM>
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __HALO__
#define __HALO__ 


// ===================================================================================================


namespace cbl {

  namespace catalogue {
    
    /**
     *  @class Halo Halo.h "Headers/Halo.h"
     *
     *  @brief The class Halo
     *
     *  This class is used to handle objects of type <EM> halo
     *  </EM>
     */
    class Halo : public Object {

    protected :

      /// halo peculiar velocity along the x direction
      double m_vx;

      /// halo peculiar velocity along the y direction
      double m_vy;

      /// halo peculiar velocity along the z direction
      double m_vz;

      /// halo mass
      double m_mass; 


    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  
       */
      Halo ()
	: Object(), m_vx(par::defaultDouble), m_vy(par::defaultDouble), m_vz(par::defaultDouble), m_mass(par::defaultDouble) {}
  
      /**
       * @brief function that allows copying private variables of the class 
       * when an object of class Catalogue is copied
       * 
       * @return a shared pointer to the Object
       *
       */
      std::shared_ptr<Object> getShared() {
        return std::make_shared<Halo>(*this);
      }
       
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
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the halo mass
       *
       *  
       */
      Halo (const comovingCoordinates coord, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble) 
	: Object(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement), m_vx(vx), m_vy(vy), m_vz(vz), m_mass(mass) {}
      
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
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the halo mass
       *
       *  
       */
      Halo (const comovingCoordinates coord, const cosmology::Cosmology &cosm, const double z1_guess=0., const double z2_guess=10., const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble) 
	: Object(coord, cosm, z1_guess, z2_guess, weight, region, ID, field, x_displacement, y_displacement, z_displacement), m_vx(vx), m_vy(vy), m_vz(vz), m_mass(mass) {}

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
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the halo mass
       *
       *  
       */
      Halo (const observedCoordinates coord, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble) 
	: Object(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement), m_vx(vx), m_vy(vy), m_vz(vz), m_mass(mass) {}
      
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
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the halo mass
       *
       *  
       */
      Halo (const observedCoordinates coord, const CoordinateUnits inputUnits, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble) 
	: Object(coord, inputUnits, weight, region, ID, field, x_displacement, y_displacement, z_displacement), m_vx(vx), m_vy(vy), m_vz(vz), m_mass(mass) {}
      
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
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the halo mass
       *
       *  
       */
      Halo (const observedCoordinates coord, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble) 
	: Object(coord, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement), m_vx(vx), m_vy(vy), m_vz(vz), m_mass(mass) {}

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
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the halo mass
       *
       *  
       */
      Halo (const observedCoordinates coord, const CoordinateUnits inputUnits, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble) 
	: Object(coord, inputUnits, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement), m_vx(vx), m_vy(vy), m_vz(vz), m_mass(mass) {}

      /**
       *  @brief constructor that uses both comoving and observed
       *  coordinates
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
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the halo mass
       *
       *  
       */
      Halo (const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble) 
	: Object(xx, yy, zz, ra, dec, redshift, weight, region, ID, field, x_displacement, y_displacement, z_displacement), m_vx(vx), m_vy(vy), m_vz(vz), m_mass(mass) {}
      
      /**
       *  @brief default destructor
       */
      ~Halo () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to get the private members 
       */
      ///@{
    
      /**
       *  @brief get the protected member \e m_vx
       *  @return the halo peculiar velocity along the x direction
       */
      double vx () const override
      { return m_vx; }

      /**
       *  @brief get the protected member \e m_vy
       *  @return the halo peculiar velocity along the y direction
       */
      double vy () const override
      { return m_vy; }

      /**
       *  @brief get the protected member \e m_vz
       *  @return the halo peculiar velocity along the z direction
       */
      double vz () const override
      { return m_vz; }

      /**
       *  @brief get the protected member \e m_mass
       *  @return the mass of the halo
       */
      double mass () const override
      { return m_mass; }

      ///@}
  
  
      /**
       *  @name Member functions used to set the private members 
       */
      ///@{
    
      /**
       *  @brief set the protected member \e m_vx
       *  @param vx the the halo peculiar velocity along the x direction
       */
      void set_vx (const double vx=par::defaultDouble) override
      { m_vx = vx; }
    
      /**
       *  @brief set the protected member \e m_vy
       *  @param vy the the halo peculiar velocity along the y direction
       */
      void set_vy (const double vy=par::defaultDouble) override
      { m_vy = vy; }

      /**
       *  @brief set the protected member \e m_vz
       *  @param vz the the halo peculiar velocity along the z direction
       */
      void set_vz (const double vz=par::defaultDouble) override
      { m_vz = vz; }

      /**
       *  @brief set the protected member \e m_mass
       *  @param mass the mass of the halo
       */
      void set_mass (const double mass=par::defaultDouble) override
      { m_mass = mass; }

      ///@}


      /**
       *  @name Member functions used to check if the private members are set 
       */
      ///@{
    
      /**
       *  @brief check if the protected member \e m_vx is set
       *
       *  @return true if the peculiar velocity along the x direction
       *  is set; false otherwise
       */
      bool isSet_vx () override
      { return (cbl::isSet(m_vx)) ? true : false; }
    
      /**
       *  @brief check if the protected member \e m_vy is set
       *  
       *  @return true if the peculiar velocity along the y direction
       *  is set; false otherwise
       */
      bool isSet_vy () override
      { return (cbl::isSet(m_vy)) ? true : false; }

      /**
       *  @brief check if the protected member \e m_vz is set
       *  
       *  @return true if the peculiar velocity along the z direction
       *  is set; false otherwise
       */
      bool isSet_vz () override
      { return (cbl::isSet(m_vz)) ? true : false; }

      /**
       *  @brief check if the protected member \e m_mass is set
       *  
       *  @return true if the mass is set; false otherwise
       */
      bool isSet_mass () override
      { return (cbl::isSet(m_mass)) ? true : false; }

      ///@}
    
    };
  }
}

#endif
