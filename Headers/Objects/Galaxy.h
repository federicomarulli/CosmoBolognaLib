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

/** @file Headers/Objects/Galaxy.h
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

  /** @class Galaxy Galaxy.h "Headers/Lib/Galaxy.h"
   *
   *  @brief The class Galaxy
   *
   *  This class is used to handle objects of type <EM> galaxy
   *  </EM>
   */
  class Galaxy : public GenericObject {
    
  protected :

    /// galaxy mass
    double m_mass;      

    /// galaxy magnitude
    double m_magnitude; 


  public:

    /** @brief default constructor
     *  @return object of class Galaxy, with protected members set to
     *  -1.e30
     */
    Galaxy () 
      : GenericObject(), m_mass(-1.e30), m_magnitude(-1.e30) {}
 
    /** @brief constructor that uses comoving coordinates
     *  @param xx comoving coordinate
     *  @param yy comoving coordinate
     *  @param zz comoving coordinate
     *  @param mass mass
     *  @param magnitude magnitude
     *  @param weight weight
     *  @return object of class Galaxy
     */
    Galaxy (const double xx, const double yy, const double zz, const double mass=-1.e30, const double magnitude=-1.e30, const double weight=1.) 
      : GenericObject(xx, yy, zz, weight), m_mass(mass), m_magnitude(magnitude) {}

    /** @brief constructor that uses observed coordinates
     *  @param ra Right Ascension
     *  @param dec Declination
     *  @param redshift redshift
     *  @param cosm object of class Cosmology, used to estimate comoving distances
     *  @param mass mass
     *  @param magnitude magnitude
     *  @param weight weight
     *  @return object of class Galaxy
     */
    Galaxy (const double ra, const double dec, const double redshift, const Cosmology &cosm, const double mass=-1.e30, const double magnitude=-1.e30, const double weight=1.) 
      : GenericObject(ra, dec, redshift, cosm, weight), m_mass(mass), m_magnitude(magnitude) {}

    /** @brief default destructor
     *  @return none
     */
    ~Galaxy () {}


    /** @brief get the protected member Galaxy::m_mass
     *  @return the mass of the galaxy
     */
    double mass () const override { return m_mass; }

    /** @brief get the protected member Galaxy::m_magnitude
     *  @return the magnitude of the galaxy
     */
    double magnitude () const override { return m_magnitude; }

    /** @brief set the protected member Galaxy::m_mass
     *  @param mass the mass of the galaxy
     *  @return none
     */
    void set_mass (const double mass) override { m_mass = mass; } 

    /** @brief set the protected member Galaxy::m_magnitude
     *  @param magnitude the magnitude of the galaxy
     *  @return none
     */
    void set_magnitude (const double magnitude) override { m_magnitude = magnitude; }

  };
}

#endif
