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
 *  @file Headers/Objects/Cluster.h
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
  
  namespace catalogue {

    /**
     *  @class Cluster Cluster.h "Headers/Lib/Cluster.h"
     *
     *  @brief The class Cluster
     *
     *  This class is used to handle objects of type <EM> cluster of
     *  galaxies </EM>
     */
    class Cluster : public Object { 

    private :

      /// cluster mass
      double m_mass;
    
      /// cluster richness
      double m_richness;


    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Cluster
       */
      Cluster () {}
      
      /**
       *  @brief constructor that uses comoving coordinates
       *  @param xx comoving coordinate
       *  @param yy comoving coordinate
       *  @param zz comoving coo
       *  @param weight weightrdinate
       *  @param mass mass
       *  @param richness richness
       *  @return object of class Cluster
       */
      Cluster (const double xx, const double yy, const double zz, const double weight=1., const double mass=par::defaultDouble, const double richness=par::defaultDouble) 
	: Object(xx, yy, zz, weight), m_mass(mass), m_richness(richness) {}

      /**
       *  @brief constructor that uses observed coordinates
       *  @param ra Right Ascension
       *  @param dec Declination
       *  @param redshift redshift
       *  @param cosm object of class Cosmology, used to estimate
       *  comoving distances
       *  @param weight weight
       *  @param mass mass
       *  @param richness magnitude
       *  @return object of class Cluster
       */
      Cluster (const double ra, const double dec, const double redshift, const Cosmology &cosm, const double weight=1., const double mass=par::defaultDouble, const double richness=par::defaultDouble) 
	: Object(ra, dec, redshift, cosm, weight), m_mass(mass), m_richness(richness) {}
  
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Cluster () {}

      ///@}
  
  
      /**
       *  @name Member functions used to get the private members 
       */
      ///@{
    
      /**
       *  @brief get the private member Cluster::m_mass
       *  @return the mass of the cluster
       */
      double mass () const override { return m_mass; }

      /**
       *  @brief get the private member Cluster::m_richness
       *  @return the richness of the cluster
       */
      double richness () const override { return m_richness; }

      ///@}
  
  
      /**
       *  @name Member functions used to set the private members 
       */
      ///@{
    
      /**
       *  @brief set the private member Cluster::m_mass
       *  @param mass the mass of the cluster
       *  @return none
       */
      void set_mass (const double mass=par::defaultDouble) override { m_mass = mass; }

      /**
       *  @brief set the private member Cluster::m_richness
       *  @param richness the richness of the cluster
       *  @return none
       */
      void set_richness (const double richness=par::defaultDouble) override { m_richness = richness; }

      ///@}
    
    };
  }
}

#endif

