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

/** @file Headers/Objects/Cluster.h
 *
 *  @brief The class Cluster 
 *
 *  This file defines the interface of the class Cluster, used to
 *  handle objects of type <EM> cluster of galaxies </EM>
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __CLUSTER__
#define __CLUSTER__ 


// ============================================================================================


namespace cosmobl {

  /** @class Cluster Cluster.h "Headers/Lib/Cluster.h"
   *
   *  @brief The class Cluster
   *
   *  This class is used to handle objects of type <EM> cluster of
   *  galaxies </EM>
   */
  class Cluster : public GenericObject { 

  protected :

    /// cluster mass
    double m_mass;
    
    /// cluster richness
    double m_richness;


  public:

    /** @brief default constructor
     *  @return object of class Cluster, with protected members set to
     *  -1.e30
     */
    Cluster () 
      : GenericObject(), m_mass(-1.e30), m_richness(-1.e30) {}
  
    /** @brief constructor that uses comoving coordinates
     *  @param xx comoving coordinate
     *  @param yy comoving coordinate
     *  @param zz comoving coordinate
     *  @param mass mass
     *  @param richness richness
     *  @param weight weight
     *  @return object of class Cluster
     */
    Cluster (double xx, double yy, double zz, double mass=-1.e30, double richness=-1.e30, double weight=1.) 
      : GenericObject(xx, yy, zz, weight), m_mass(mass), m_richness(richness) {}

    /** @brief constructor that uses observed coordinates
     *  @param ra Right Ascension
     *  @param dec Declination
     *  @param redshift redshift
     *  @param cosm object of class Cosmology, used to estimate comoving distances
     *  @param mass mass
     *  @param richness magnitude
     *  @param weight weight
     *  @return object of class Cluster
     */
    Cluster (double ra, double dec, double redshift, Cosmology &cosm, double mass=-1.e30, double richness=-1.e30, double weight=1.) 
      : GenericObject(ra, dec, redshift, cosm, weight), m_mass(mass), m_richness(richness) {}
  
    /** @brief default destructor
     *  @return none
     */
    ~Cluster () {}


    /** @brief get the protected member Cluster::m_mass
     *  @return the mass of the cluster
     */
    double mass () { return m_mass; }

    /** @brief get the protected member Cluster::m_richness
     *  @return the richness of the cluster
     */
    double richness () { return m_richness; }


    /** @brief set the protected member Cluster::m_mass
     *  @param mass the mass of the cluster
     *  @return none
     */
    void set_mass (double mass) { m_mass = mass; }

    /** @brief set the protected member Cluster::m_richness
     *  @param richness the richness of the cluster
     *  @return none
     */
    void set_richness (double richness) { m_richness = richness; }

  };
}

#endif

