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

/** @file Headers/Objects/GenericObject.h
 *
 *  @brief The class GenericObject 
 *
 *  This file defines the interface of the class GenericObject, used
 *  to handle generic astronomical objects
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __GENERICOBJ__
#define __GENERICOBJ__ 


// ===================================================================================================


namespace cosmobl {

  /** @class GenericObject GenericObject.h "Headers/Lib/GenericObject.h"
   *
   *  @brief The class GenericObject
   *
   *  This class is used to handle objects of type <EM> generic object
   *  </EM>
   */
  class GenericObject : public Object {

  protected :

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

  public:

    /** @brief default constructor
     *  @return object of class GenericObject, with protected members
     *  set to -1.e30, exept m_weight=1.
     */
    GenericObject () 
      : m_xx(-1.e30), m_yy(-1.e30), m_zz(-1.e30), m_ra(-1.e30), m_dec(-1.e30), m_redshift(-1.e30), m_dc(-1.e30), m_weight(1.), m_region(-1) {}
    
    /** @brief constructor that uses comoving coordinates
     *  @param xx comoving coordinate
     *  @param yy comoving coordinate
     *  @param zz comoving coordinate
     *  @param weight weight
     *  @return object of class GenericObject
     */
    GenericObject (double xx, double yy, double zz, double weight=1.) 
      : m_xx(xx), m_yy(yy), m_zz(zz), m_ra(-1.e30), m_dec(-1.e30), m_redshift(-1.e30), m_dc(-1.e30), m_weight(weight), m_region(-1) {}

    /** @brief constructor that uses comoving coordinates
     *  @param xx comoving coordinate
     *  @param yy comoving coordinate
     *  @param zz comoving coordinate 
     *  @param cosm object of class Cosmology, used to estimate comoving distances
     *  @param z1_guess minimum redshift used to search the redshift
     *  @param z2_guess maximum redshift used to search the redshift
     *  @param weight weight
     *  @return object of class GenericObject
     */
    GenericObject (double xx, double yy, double zz, Cosmology &cosm, double &z1_guess, double &z2_guess, double weight=1.) 
      : m_xx(xx), m_yy(yy), m_zz(zz), m_dc(-1.), m_weight(weight), m_region(-1)
      {
	 cosmobl::polar_coord(m_xx, m_yy, m_zz, &m_ra, &m_dec, &m_dc);
	 m_redshift = cosm.Redshift(m_dc, z1_guess, z2_guess);
      }
    
    /** @brief constructor that uses observed coordinates
     *  @param ra Right Ascension
     *  @param dec Declination
     *  @param redshift redshift
     *  @param cosm object of class Cosmology, used to estimate comoving distances
     *  @param weight weight
     *  @return object of class GenericObject
     */
    GenericObject (double ra, double dec, double redshift, Cosmology &cosm, double weight=1.) 
      : m_ra(ra), m_dec(dec), m_redshift(redshift), m_dc(-1.), m_weight(weight), m_region(-1)
    { 
      m_dc = cosm.D_C(m_redshift);
      cosmobl::cartesian_coord(m_ra, m_dec, m_dc, &m_xx, &m_yy, &m_zz);
    }

    /** @brief constructor that uses both comoving and observed coordinates
     *  @param xx comoving coordinate
     *  @param yy comoving coordinate
     *  @param zz comoving coordinate 
     *  @param ra Right Ascension
     *  @param dec Declination
     *  @param redshift redshift
     *  @param weight weight
     *  @return object of class GenericObject
     */
    GenericObject (double xx, double yy, double zz, double ra, double dec, double redshift, double weight=1.) 
      : m_xx(xx), m_yy(yy), m_zz(zz), m_ra(ra), m_dec(dec), m_redshift(redshift), m_dc(sqrt(xx*xx+yy*yy+zz*zz)), m_weight(weight), m_region(-1) {}
    
    /** @brief default destructor
     *  @return none
     */
    ~GenericObject () {}


    /** @brief get the protected member GenericObject::m_xx
     *  @return the coordinate x of the object
     */
    double xx () { return m_xx; }
    
    /** @brief get the protected member GenericObject::m_yy
     *  @return the coordinate y of the object
     */
    double yy () { return m_yy; }
    
    /** @brief get the protected member GenericObject::m_zz
     *  @return the coordinate z of the object
     */
    double zz () { return m_zz; }
    
    /** @brief get the protected member GenericObject::m_dc
     *  @return the comoving distance of the object
     */
    double dc () { return m_dc; }
    
    /** @brief get the protected member GenericObject::m_ra
     *  @return the Right Ascension of the object
     */
    double ra () { return m_ra; } 
    
    /** @brief get the protected member GenericObject::m_dec
     *  @return the Declination of the object
     */
    double dec () { return m_dec; }
    
    /** @brief get the protected member GenericObject::m_redshift
     *  @return the redshift of the object
     */
    double redshift () { return m_redshift; }
    
    /** @brief get the protected member GenericObject::m_weight
     *  @return the weight of the object
     */
    double weight () { return m_weight; }

    /** @brief get the protected member GenericObject::m_region
     *  @return the index of the subRegion in which the object is located
     */
    long region () { return m_region; }

    /** @brief get the object coordinates
     *  @return a vector containing the object coordinates
     */
    vector<double> coords () { return {m_xx, m_yy, m_zz}; }
    
    /** @brief set the protected member GenericObject::m_xx
     *  @param xx the coordinate x of the object
     *  @return none
     */
    void set_xx (double xx) { m_xx = xx; }
 
    /** @brief set the protected member GenericObject::m_yy
     *  @param yy the coordinate y of the object
     *  @return none
     */
    void set_yy (double yy) { m_yy = yy; }
    
    /** @brief set the protected member GenericObject::m_zz
     *  @param zz the coordinate z of the object
     *  @return none
     */
    void set_zz (double zz) { m_zz = zz; }
    
    /** @brief set the protected member GenericObject::m_ra
     *  @param ra the Right Ascension of the object
     *  @return none
     */
    void set_ra (double ra) { m_ra = ra; }
    
    /** @brief set the protected member GenericObject::m_dec
     *  @param dec the Declination of the object
     *  @return none
     */
    void set_dec (double dec) { m_dec = dec; }
    
    /** @brief set the protected member GenericObject::m_redshift
     *  @param redshift the redshift of the object
     *  @return none
     */
    void set_redshift (double redshift) { m_redshift = redshift; }
    
    /** @brief set the protected member GenericObject::m_dc
     *  @param dc the comoving distance of the object
     *  @return none
     */
    void set_dc (double dc) { m_dc = dc; }
    
    /** @brief set the protected member GenericObject::m_weight
     *  @param weight the weight of the object
     *  @return none
     */
    void set_weight (double weight) { m_weight = weight; }

    /** @brief set the protected member GenericObject::m_region
     *  @param region the index of the subRegion in which the object is located
     *  @return none
     */
    void set_region (long region) { m_region = region; }

  };

  /** @class RandomObject GenericObject.h "Headers/Lib/GenericObject.h"
   *
   *  @brief The class RandomObject
   *
   *  This class is used to handle objects of type <EM> random object
   *  </EM>
   */
  class RandomObject : public GenericObject {
    
  public :
    
    /** @brief default constructor
     *  @return object of class RandomObject, with protected members
     *  set to -1.e30, exept m_weight=1.
     */
    RandomObject() 
      : GenericObject() {}
    
    /** @brief constructor that uses comoving coordinates
     *  @param xx comoving coordinate
     *  @param yy comoving coordinate
     *  @param zz comoving coordinate
     *  @param weight weight
     *  @return object of class RandomObject
     */
    RandomObject(double xx, double yy, double zz, double weight=1.) 
      : GenericObject(xx, yy, zz, weight) {}

    /** @brief constructor that uses comoving coordinates
     *  @param xx comoving coordinate
     *  @param yy comoving coordinate
     *  @param zz comoving coordinate
     *  @param cosm object of class Cosmology, used to estimate comoving distances
     *  @param z1_guess minimum redshift used to search the redshift
     *  @param z2_guess maximum redshift used to search the redshift
     *  @param weight weight
     *  @return object of class RandomObject
     */
    RandomObject(double xx, double yy, double zz, Cosmology &cosm, double &z1_guess, double &z2_guess, double weight=1.) 
      : GenericObject(xx, yy, zz, cosm, z1_guess, z2_guess, weight) {}
    
    /** @brief constructor that uses both comoving and observed coordinates
     *  @param xx comoving coordinate
     *  @param yy comoving coordinate
     *  @param zz comoving coordinate
     *  @param ra Right Ascension
     *  @param dec Declination
     *  @param redshift redshift
     *  @param weight weight
     *  @return object of class RandomObject
     */
    RandomObject(double xx, double yy, double zz, double ra, double dec, double redshift, double weight=1.) 
      : GenericObject(xx, yy, zz, ra, dec, redshift, weight) {}

    /** @brief constructor that uses observed coordinates
     *  @param ra Right Ascension
     *  @param dec Declination
     *  @param redshift redshift
     *  @param cosm object of class Cosmology, used to estimate comoving distances
     *  @param weight weight
     *  @return object of class GenericObject
     */
    RandomObject (double ra, double dec, double redshift, Cosmology &cosm, double weight=1.) 
      : GenericObject(ra, dec, redshift, cosm, weight) {} 
    
  };
  
}

#endif
