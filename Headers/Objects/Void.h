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

/** @file Headers/Objects/Void.h
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

  /** @class Void Void.h "Headers/Lib/Void.h"
   *
   *  @brief The class Void
   *
   *  This class is used to handle objects of type <EM> void
   *  </EM>
   */
  class Void : public GenericObject {

  protected:

    /// void normalized volume
    double m_volNorm;

    /// void radius [Mpc/h]
    double m_radius;

    /// void volume [(Mpc/h)^3]
    double m_volume;

    /// void ID
    int m_ID;

    /// void density contrast (\f${\rho_v}/{\rho_m}\f$) 
    double m_densContr;

    /// ID of the void parent (-1 for main voids)
    int m_parentID;

    /// void hierarchy level (0 for main voids)
    int m_treeLevel;

    /// number of children
    int m_child;

    /// void central density
    double m_rho0;

    /// void core density (normalized)
    double m_rho0Norm;

  public:

    /** @brief default constructor
     *  @return object of class Void, with protected members set to
     *  -1.e30
     */		
  Void()
    : GenericObject(), m_volNorm(-1.e30), m_radius(-1.e30), m_volume(-1.e30), m_ID(-1e4), m_densContr(-1.e30), m_parentID(-1e4), m_treeLevel(-1e4), m_child(-1e4), m_rho0(-1.e30), m_rho0Norm(-1.e30) {} 

    /**
     *  @brief constructor that uses comoving coordinates
     *  @param xx comoving coordinate
     *  @param yy comoving coordinate
     *  @param zz comoving coordinate
     *  @param volNorm normalized volume
     *  @param radius radius 
     *  @param volume volume 
     *  @param ID identification number
     *  @param densContr density contrast
     *  @param parentID identification number of the parent void
     *  @param treeLevel hierarchy level
     *  @param child number of childrens
     *  @param rho0 central density
     *  @param rho0Norm normalized core density
     *  @param weight weight
     *  @return object of class Void
     */
    Void (const double xx, const double yy, const double zz, const double volNorm=-1.e30, const double radius=-1.e30, const double volume=-1.e30, const int ID=-1e4, const double densContr=-1.e30, const int parentID=-1e4, const int treeLevel=-1e4, const int child=-1e4, const double rho0=-1.e30, const double rho0Norm=-1.e30, const double weight = 1.)
      : GenericObject(xx, yy, zz, weight), m_volNorm(volNorm), m_radius(radius), m_volume(volume), m_ID(ID), m_densContr(densContr), m_parentID(parentID), m_treeLevel(treeLevel), m_child(child), m_rho0(rho0), m_rho0Norm(rho0Norm) {}

    /**
     *  @brief default destructor
     *  @return none
     */
    ~Void () {}


    /** @brief get the protected member Void::m_volNorm
     *  @return the normalized volume of the void
     */
    double volNorm () const { return m_volNorm; }

    /** @brief get the protected member Void::m_radius
     *  @return the radius of the spherical void 
     */
    double radius () const { return m_radius; }

    /** @brief get the protected member Void::m_volume
     *  @return the volume of the void 
     */
    double volume () const { return m_volume; }

    /** @brief get the protected member Void::m_ID
     *  @return the identification number of the void
     */
    int ID () const { return m_ID; }

    /** @brief get the protected member Void::m_densContr
     *  @return the density contrast between the void and the Universe
     */
    double densContr () const { return m_densContr; }

    /** @brief get the protected member Void::m_parentID
     *  @return the identification number of the parent void
     */
    int parentID () const { return m_parentID; }

    /** @brief get the protected member Void::m_treeLevel
     *  @return the hierarchy level of the void
     */
    int treeLevel () const { return m_treeLevel; }

    /** @brief get the protected member Void::m_child
     *  @return the number of children of the void
     */
    int child () const { return m_child; }

    /** @brief get the protected member Void::m_rho0
     *  @return the central density of the void
     */
    double rho0 () const { return m_rho0; }

    /** @brief get the protected member Void::m_rho0Norm
     *  @return the normalized core density of the void
     */
    double rho0Norm () const { return m_rho0Norm; }

    /** @brief set the protected member Void::m_volNorm
     *  @param volNorm the normalized volume of the void
     *  @return none
     */
    void set_volNorm (const double volNorm) { m_volNorm = volNorm; }
	
    /** @brief get the protected member Void::m_radius
     *  @param radius the radius of the spherical void [Mpc/h]
     *  @return none
     */
    void set_radius (const double radius) { m_radius = radius; }
	
    /** @brief get the protected member Void::m_volume
     *  @param volume the volume of the void [(Mpc/h)^3]
     *  @return none
     */
    void set_volume (const double volume) { m_volume = volume; }
	
    /** @brief get the protected member Void::m_ID
     *  @param ID the identification number of the void
     *  @return none
     */
    void set_ID (const int ID) { m_ID = ID; } 
	
    /** @brief get the protected member Void::m_densContr
     *  @param densContr the density contrast between the void and the Universe
     *  @return none
     */
    void set_densContr (const double densContr) { m_densContr = densContr; }
	
    /** @brief get the protected member Void::m_parentID
     *  @param parentID the identification number of the parent void
     *  @return none
     */
    void set_parentID (int parentID) { m_parentID = parentID; }
	
    /** @brief get the protected member Void::m_treeLevel
     *  @param treeLevel the hierarchy level of the void
     *  @return none
     */
    void set_treeLevel (int treeLevel) { m_treeLevel = treeLevel; }
	
    /** @brief get the protected member Void::m_child
     *  @param child the number of children of the void
     *  @return none
     */
    void set_child (int child) { m_child = child; }
	
    /** @brief get the protected member Void::m_rho0
     *  @param rho0 the central density of the void
     *  @return none
     */
    void set_rho0 (const double rho0) { m_rho0 = rho0; }    
	
    /** @brief get the protected member Void::m_rho0Norm
     *  @param rho0Norm the normalized core density of the void
     *  @return none
     */
    void set_rho0Norm (const double rho0Norm) { m_rho0Norm = rho0Norm; }

  };
}

#endif
