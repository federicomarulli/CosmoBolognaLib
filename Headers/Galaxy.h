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
 *  @file Headers/Galaxy.h
 *
 *  @brief The class Galaxy 
 *
 *  This file defines the interface of the class Galaxy, used to
 *  handle objects of type <EM> galaxies </EM>
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#ifndef __GALAXY__
#define __GALAXY__ 


// ===================================================================================================


namespace cbl {

  namespace catalogue {
    
    /**
     *  @class Galaxy Galaxy.h "Headers/Galaxy.h"
     *
     *  @brief The class Galaxy
     *
     *  This class is used to handle objects of type <EM> galaxy
     *  </EM>
     */
    class Galaxy : public Object {
    
    private :

      /// mass
      double m_mass;      

      /// magnitude
      double m_magnitude; 

      /// star formation rate
      double m_SFR; 

      /// specific star formation rate
      double m_sSFR; 


    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{
      
      /**
       *  @brief default constructor
       *  @return object of class Galaxy
       */
      Galaxy ()
	: Object(), m_mass(par::defaultDouble), m_magnitude(par::defaultDouble), m_SFR(par::defaultDouble), m_sSFR(par::defaultDouble) {}

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
       *  @param mass the galaxy mass
       *
       *  @param magnitude the galaxy magnitude
       *
       *  @param SFR the galaxy star formation rate
       *
       *  @param sSFR the galaxy specific star formation rate
       *
       *  @return object of class Galaxy
       */
      Galaxy (const comovingCoordinates coord, const double weight=1., const long region=par::defaultLong, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double mass=par::defaultDouble, const double magnitude=par::defaultDouble, const double SFR=par::defaultDouble, const double sSFR=par::defaultDouble) 
	: Object(coord, weight, region, field, x_displacement, y_displacement, z_displacement), m_mass(mass), m_magnitude(magnitude), m_SFR(SFR), m_sSFR(sSFR) {}

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
       *  @param mass the galaxy mass
       *
       *  @param magnitude the galaxy magnitude
       *
       *  @param SFR the galaxy star formation rate
       *
       *  @param sSFR the galaxy specific star formation rate
       *
       *  @return object of class Galaxy
       */
      Galaxy (const comovingCoordinates coord, const cosmology::Cosmology &cosm, const double z1_guess=0., const double z2_guess=10., const double weight=1., const long region=par::defaultLong, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double mass=par::defaultDouble, const double magnitude=par::defaultDouble, const double SFR=par::defaultDouble, const double sSFR=par::defaultDouble) 
	: Object(coord, cosm, z1_guess, z2_guess, weight, region, field, x_displacement, y_displacement, z_displacement), m_mass(mass), m_magnitude(magnitude), m_SFR(SFR), m_sSFR(sSFR) {}

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
       *  @param mass the galaxy mass
       *
       *  @param magnitude the galaxy magnitude
       *
       *  @param SFR the galaxy star formation rate
       *
       *  @param sSFR the galaxy specific star formation rate
       *
       *  @return object of class Galaxy
       */
      Galaxy (const observedCoordinates coord, const double weight=1., const long region=par::defaultLong, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double mass=par::defaultDouble, const double magnitude=par::defaultDouble, const double SFR=par::defaultDouble, const double sSFR=par::defaultDouble) 
	: Object(coord, weight, region, field, x_displacement, y_displacement, z_displacement), m_mass(mass), m_magnitude(magnitude), m_SFR(SFR), m_sSFR(sSFR) {}
      
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
       *  @param mass the galaxy mass
       *
       *  @param magnitude the galaxy magnitude
       *
       *  @param SFR the galaxy star formation rate
       *
       *  @param sSFR the galaxy specific star formation rate
       *
       *  @return object of class Galaxy
       */
      Galaxy (const observedCoordinates coord, const CoordinateUnits inputUnits, const double weight=1., const long region=par::defaultLong, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double mass=par::defaultDouble, const double magnitude=par::defaultDouble, const double SFR=par::defaultDouble, const double sSFR=par::defaultDouble) 
	: Object(coord, inputUnits, weight, region, field, x_displacement, y_displacement, z_displacement), m_mass(mass), m_magnitude(magnitude), m_SFR(SFR), m_sSFR(sSFR) {}
      
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
       *  @param mass the galaxy mass
       *
       *  @param magnitude the galaxy magnitude
       *
       *  @param SFR the galaxy star formation rate
       *
       *  @param sSFR the galaxy specific star formation rate
       *
       *  @return object of class Galaxy
       */
      Galaxy (const observedCoordinates coord, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double mass=par::defaultDouble, const double magnitude=par::defaultDouble, const double SFR=par::defaultDouble, const double sSFR=par::defaultDouble) 
	: Object(coord, cosm, weight, region, field, x_displacement, y_displacement, z_displacement), m_mass(mass), m_magnitude(magnitude), m_SFR(SFR), m_sSFR(sSFR) {}

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
       *  @param mass the galaxy mass
       *
       *  @param magnitude the galaxy magnitude
       *
       *  @param SFR the galaxy star formation rate
       *
       *  @param sSFR the galaxy specific star formation rate
       *
       *  @return object of class Galaxy
       */
      Galaxy (const observedCoordinates coord, const CoordinateUnits inputUnits, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double mass=par::defaultDouble, const double magnitude=par::defaultDouble, const double SFR=par::defaultDouble, const double sSFR=par::defaultDouble) 
	: Object(coord, inputUnits, cosm, weight, region, field, x_displacement, y_displacement, z_displacement), m_mass(mass), m_magnitude(magnitude), m_SFR(SFR), m_sSFR(sSFR) {}

      /**
       *  @brief constructor that uses both comoving and observed
       *  coordinates
       *
       *  @param xx comoving coordinate
       *  @param yy comoving coordinate
       *  @param zz comoving coordinate 
       *  @param ra Right Ascension
       *  @param dec Declination
       *  @param redshift redshift
       *  @param weight weight  
       *  @param region region, used e.g. for jackknife and boostrap 
       *  @param field the field where the object has been observed
       *  @param x_displacement the displacement along the x-axis
       *  @param y_displacement the displacement along the y-axis
       *  @param z_displacement the displacement along the z-axis
       *  @param mass the galaxy mass
       *  @param magnitude the galaxy magnitude
       *  @param SFR the galaxy star formation rate
       *  @param sSFR the galaxy specific star formation rate
       *
       *  @return object of class Galaxy
       */
      Galaxy (const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight=1., const long region=par::defaultLong, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double mass=par::defaultDouble, const double magnitude=par::defaultDouble, const double SFR=par::defaultDouble, const double sSFR=par::defaultDouble) 
	: Object(xx, yy, zz, ra, dec, redshift, weight, region, field, x_displacement, y_displacement, z_displacement), m_mass(mass), m_magnitude(magnitude), m_SFR(SFR), m_sSFR(sSFR) {}
      
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Galaxy () = default;
    
      ///@}
  
  
      /**
       *  @name Member functions used to get the private members 
       */
      ///@{
    
      /**
       *  @brief get the private member \e m_mass
       *  @return the mass of the galaxy
       */
      double mass () const override
      { return m_mass; }

      /**
       *  @brief get the private member \e m_magnitude
       *  @return the magnitude of the galaxy
       */
      double magnitude () const override
      { return m_magnitude; }

      /**
       *  @brief get the private member \e m_SFR
       *  @return the star formation rate of the galaxy
       */
      double SFR () const override
      { return m_SFR; }

      /**
       *  @brief get the private member \e m_sSFR
       *  @return the specific star formation rate of the galaxy
       */
      double sSFR () const override
      { return m_sSFR; }

      ///@}
  
  
      /**
       *  @name Member functions used to set the private members 
       */
      ///@{
    
      /**
       *  @brief set the private member \e m_mass
       *  @param mass the mass of the galaxy
       *  @return none
       */
      void set_mass (const double mass=par::defaultDouble) override
      { m_mass = mass; } 

      /**
       *  @brief set the private member \e m_magnitude
       *  @param magnitude the magnitude of the galaxy
       *  @return none
       */
      void set_magnitude (const double magnitude=par::defaultDouble) override
      { m_magnitude = magnitude; }

      /**
       *  @brief set the private member \e m_SFR
       *  @param SFR the star formation rate of the galaxy
       *  @return none
       */
      void set_SFR (const double SFR=par::defaultDouble) override
      { m_SFR = SFR; }
    
      /**
       *  @brief set the private member \e m_sSFR
       *  @param sSFR the specific star formation rate of the galaxy
       *  @return none
       */
      void set_sSFR (const double sSFR=par::defaultDouble) override
      { m_sSFR = sSFR; }
    
      ///@}


      /**
       *  @name Member functions used to check if the private members are set 
       */
      ///@{
    
      /**
       *  @brief check if the private member \e m_mass is set
       *
       *  @return none
       */
      bool isSet_mass () override
      { return (cbl::isSet(m_mass)) ? true : false; } 

      /**
       *  @brief check if the private member \e m_magnitude is set
       *  
       *  @return true if the magnitude is set; false otherwise
       */
      bool isSet_magnitude () override
      { return (cbl::isSet(m_magnitude)) ? true : false; }

      /**
       *  @brief check if the private member \e m_SFR is set
       *  
       *  @return true if the SFR is set; false otherwise
       */
      bool isSet_SFR () override
      { return (cbl::isSet(m_SFR)) ? true : false; }
    
      /**
       *  @brief check if the private member \e m_sSFR is set
       *  
       *  @return true if the sSFR is set; false otherwise
       */
      bool isSet_sSFR () override
      { return (cbl::isSet(m_sSFR)) ? true : false; }
    
      ///@}
      
    };
  }
}

#endif
