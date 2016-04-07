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
 *  @file Headers/Objects/RandomObject.h
 *
 *  @brief The class RandomObject 
 *
 *  This file defines the interface of the class RandomObject, used to
 *  handle objects of type <EM> galaxies </EM>
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __RANDOMOBJ__
#define __RANDOMOBJ__ 

#include "Object.h"


// ===================================================================================================


namespace cosmobl {

  namespace catalogue {
    
    /**
     *  @class RandomObject RandomObject.h "Headers/Lib/RandomObject.h"
     *
     *  @brief The class RandomObject
     *
     *  This class is used to handle objects of type <EM> RandomObject
     *  </EM>
     */
    class RandomObject : public Object {
    
    public :
    
      /**
       *  @brief default constructor
       *  @return object of class RandomObject
       */
      RandomObject () {}
    
      /**
       *  @brief constructor that uses comoving coordinates
       *  @param xx comoving coordinate
       *  @param yy comoving coordinate
       *  @param zz comoving coordinate
       *  @param weight weight
       *  @return object of class RandomObject
       */
      RandomObject (const double xx, const double yy, const double zz, const double weight=1.) 
	: Object(xx, yy, zz, weight) {}

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
       *  @return object of class RandomObject
       */
      RandomObject (const double xx, const double yy, const double zz, const Cosmology &cosm, const double z1_guess, const double z2_guess, const double weight=1.) 
	: Object(xx, yy, zz, cosm, z1_guess, z2_guess, weight) {}
    
      /**
       *  @brief constructor that uses both comoving and observed coordinates
       *  @param xx comoving coordinate
       *  @param yy comoving coordinate
       *  @param zz comoving coordinate
       *  @param ra Right Ascension
       *  @param dec Declination
       *  @param redshift redshift
       *  @param weight weight
       *  @return object of class RandomObject
       */
      RandomObject (const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight=1.) 
	: Object(xx, yy, zz, ra, dec, redshift, weight) {}

      /**
       *  @brief constructor that uses observed coordinates
       *  @param ra Right Ascension
       *  @param dec Declination
       *  @param redshift redshift
       *  @param cosm object of class Cosmology, used to estimate comoving distances
       *  @param weight weight
       *  @return object of class RandomObject
       */
      RandomObject (const double ra, const double dec, const double redshift, const Cosmology &cosm, const double weight=1.) 
	: Object(ra, dec, redshift, cosm, weight) {} 
    
    };  
  }
}

#endif
