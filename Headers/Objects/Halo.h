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

/** @file Headers/Objects/Halo.h
 *
 *  @brief The class Halo 
 *
 *  This file defines the interface of the class Halo, used to handle
 *  objects of type <EM> dark matter halo </EM>
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __HALO__
#define __HALO__ 


// ===================================================================================================


namespace cosmobl {

  /** @class Halo Halo.h "Headers/Lib/Halo.h"
   *
   *  @brief The class Halo
   *
   *  This class is used to handle objects of type <EM> halo
   *  </EM>
   */
  class Halo : public GenericObject {

  protected :

    /// halo peculiar velocity along the x direction
    double m_vx;

    /// halo peculiar velocity along the y direction
    double m_vy;

    /// halo peculiar velocity along the y direction
    double m_vz;

    /// halo mass
    double m_mass; 


  public:

    /** @brief default constructor
     *  @return object of class Halo, with protected members set to
     *  -1.e30
     */
    Halo () 
      : GenericObject(), m_vx(-1.e30), m_vy(-1.e30), m_vz(-1.e30), m_mass(-1.e30) {}

    /** @brief constructor that uses comoving coordinates
     *  @param xx comoving coordinate
     *  @param yy comoving coordinate
     *  @param zz comoving coordinate
     *  @param vx peculiar velocity
     *  @param vy peculiar velocity
     *  @param vz peculiar velocity    
     *  @param mass mass
     *  @param weight weight
     *  @return object of class Halo
     */
    Halo (const double xx, const double yy, const double zz, const double vx=-1.e30, const double vy=-1.e30, const double vz=-1.e30, const double mass=-1.e30, const double weight=1.) 
      : GenericObject(xx, yy, zz, weight), m_vx(vx), m_vy(vy), m_vz(vz), m_mass(mass) {}
    
    /** @brief constructor that uses observed coordinates
     *  @param ra Right Ascension
     *  @param dec Declination
     *  @param redshift redshift
     *  @param cosm object of class Cosmology, used to estimate comoving distances
     *  @param vx peculiar velocity
     *  @param vy peculiar velocity
     *  @param vz peculiar velocity    
     *  @param mass mass
     *  @param weight weight
     *  @return object of class Halo
     */
    Halo (const double ra, const double dec, const double redshift, const Cosmology &cosm, const double vx=-1.e30, const double vy=-1.e30, const double vz=-1.e30, const double mass=-1.e30, const double weight=1.) 
      : GenericObject(ra, dec, redshift, cosm, weight), m_vx(vx), m_vy(vy), m_vz(vz), m_mass(mass) {}

    /** @brief default destructor
     *  @return none
     */
    ~Halo () {}

    /** @brief get the protected member Halo::m_vx
     *  @return the halo peculiar velocity along the x direction
     */
    double vx () const override { return m_vx; }

    /** @brief get the protected member Halo::m_vy
     *  @return the halo peculiar velocity along the y direction
     */
    double vy () const override { return m_vy; }

    /** @brief get the protected member Halo::m_vz
     *  @return the halo peculiar velocity along the z direction
     */
    double vz () const override { return m_vz; }

    /** @brief get the protected member Halo::m_mass
     *  @return the mass of the halo
     */
    double mass () const override { return m_mass; }

    /** @brief set the protected member Halo::m_vx
     *  @param vx the the halo peculiar velocity along the x direction
     *  @return none
     */
    void set_vx (const double vx) override { m_vx = vx; }
    
    /** @brief set the protected member Halo::m_vy
     *  @param vy the the halo peculiar velocity along the y direction
     *  @return none
     */
    void set_vy (const double vy) override { m_vy = vy; }

    /** @brief set the protected member Halo::m_vz
     *  @param vz the the halo peculiar velocity along the z direction
     *  @return none
     */
    void set_vz (const double vz) override { m_vz = vz; }

    /** @brief set the protected member Halo::m_mass
     *  @param mass the mass of the halo
     *  @return none
     */
    void set_mass (const double mass) override { m_mass = mass; }

  };
}

#endif
