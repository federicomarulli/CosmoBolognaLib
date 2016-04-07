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
 *  @file Headers/Objects/Galaxy.h
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


namespace cosmobl {

  namespace catalogue {
    
    /**
     *  @class Galaxy Galaxy.h "Headers/Lib/Galaxy.h"
     *
     *  @brief The class Galaxy
     *
     *  This class is used to handle objects of type <EM> galaxy
     *  </EM>
     */
    class Galaxy : public Object {
    
    private :

      /// galaxy mass
      double m_mass;      

      /// galaxy magnitude
      double m_magnitude; 


    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{
      
      /**
       *  @brief default constructor
       *  @return object of class Galaxy
       */
      Galaxy () {}

      /**
       *  @brief constructor that uses comoving coordinates
       *  @param xx comoving coordinate
       *  @param yy comoving coordinate
       *  @param zz comoving coordinate
       *  @param weight weight
       *  @param mass mass
       *  @param magnitude magnitude
       *  @return object of class Galaxy
       */
      Galaxy (const double xx, const double yy, const double zz, const double weight=1., const double mass=par::defaultDouble, const double magnitude=par::defaultDouble) 
	: Object(xx, yy, zz, weight), m_mass(mass), m_magnitude(magnitude) {}

      /**
       *  @brief constructor that uses observed coordinates
       *  @param ra Right Ascension
       *  @param dec Declination
       *  @param redshift redshift
       *  @param cosm object of class Cosmology, used to estimate
       *  comoving distances
       *  @param weight weight
       *  @param mass mass
       *  @param magnitude magnitude
       *  @return object of class Galaxy
       */
      Galaxy (const double ra, const double dec, const double redshift, const Cosmology &cosm, const double weight=1., const double mass=par::defaultDouble, const double magnitude=par::defaultDouble) 
	: Object(ra, dec, redshift, cosm, weight), m_mass(mass), m_magnitude(magnitude) {}

      /**
       *  @brief default destructor
       *  @return none
       */
      ~Galaxy () {}
    
      ///@}
  
  
      /**
       *  @name Member functions used to get the private members 
       */
      ///@{
    
      /**
       *  @brief get the private member Galaxy::m_mass
       *  @return the mass of the galaxy
       */
      double mass () const override { return m_mass; }

      /**
       *  @brief get the private member Galaxy::m_magnitude
       *  @return the magnitude of the galaxy
       */
      double magnitude () const override { return m_magnitude; }

      ///@}
  
  
      /**
       *  @name Member functions used to set the private members 
       */
      ///@{
    
      /**
       *  @brief set the private member Galaxy::m_mass
       *  @param mass the mass of the galaxy
       *  @return none
       */
      void set_mass (const double mass=par::defaultDouble) override { m_mass = mass; } 

      /**
       *  @brief set the private member Galaxy::m_magnitude
       *  @param magnitude the magnitude of the galaxy
       *  @return none
       */
      void set_magnitude (const double magnitude=par::defaultDouble) override { m_magnitude = magnitude; }
    
      ///@}
    
    };
  }
}

#endif
