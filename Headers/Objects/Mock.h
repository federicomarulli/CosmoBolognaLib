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

/** @file Headers/Objects/Mock.h
 *
 *  @brief The class Mock 
 *
 *  This file defines the interface of the class Mock, used to handle
 *  objects of type <EM> mock object </EM>
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MOCK__
#define __MOCK__ 


// ===================================================================================================


namespace cosmobl {

  /** @class Mock Mock.h "Headers/Lib/Mock.h"
   *
   *  @brief The class Mock
   *
   *  This class is used to handle objects of type <EM> mock object
   *  </EM>
   */
  class Mock : public Halo {

  protected :

    /// generic variable of the mock object
    double m_generic; 


  public:

    /** @brief default constructor
     *  @return object of class Mock, with protected members set to
     *  -1.e30
     */
    Mock () 
      : Halo(), m_generic(-1.e30) {}
    
    /** @brief constructor that uses comoving coordinates
     *  @param xx comoving coordinate
     *  @param yy comoving coordinate
     *  @param zz comoving coordinate
     *  @param vx peculiar velocity
     *  @param vy peculiar velocity
     *  @param vz peculiar velocity    
     *  @param mass mass
     *  @param generic generic variable
     *  @param weight weight
     *  @return object of class Mock
     */
    Mock (double xx, double yy, double zz, double vx=-1.e30, double vy=-1.e30, double vz=-1.e30, double mass=-1.e30, double generic=-1.e30, double weight=1.) 
      : Halo(xx, yy, zz, vx, vy, vz, mass, weight), m_generic(generic) {}

    /** @brief constructor that uses observed coordinates
     *  @param ra Right Ascension
     *  @param dec Declination
     *  @param redshift redshift
     *  @param cosm object of class Cosmology, used to estimate comoving distances
     *  @param vx peculiar velocity
     *  @param vy peculiar velocity
     *  @param vz peculiar velocity    
     *  @param mass mass
     *  @param generic generic variable
     *  @param weight weight
     *  @return object of class Mock
     */
    Mock (double ra, double dec, double redshift, Cosmology &cosm, double vx=-1.e30, double vy=-1.e30, double vz=-1.e30, double mass=-1.e30, double generic=-1.e30, double weight=1.) 
      : Halo(ra, dec, redshift, cosm, vx, vy, vz, mass, weight), m_generic(generic) {}

    /** @brief default destructor
     *  @return none
     */
    ~Mock () {}

    /** @brief get the protected member Mock::m_generic
     *  @return the generic variable of the mock object
     */
    double generic ()  { return m_generic; }

    /** @brief set the protected member Mock::m_generic
     *  @param generic the generic variable of the mock object
     *  @return none
     */
    void set_generic (double generic)  { m_generic = generic; }

  };
}

#endif
