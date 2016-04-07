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

/**
 *  @file Headers/Objects/Void.h
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

  namespace catalogue {
    
    /**
     *  @class Void Void.h "Headers/Lib/Void.h"
     *
     *  @brief The class Void
     *
     *  This class is used to handle objects of type <EM> void
     *  </EM>
     */
    class Void : public Object {

    private:

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
    
      /**
       *  @name Constructors/destructors
       */
      ///@{
      
      /**
       *  @brief default constructor
       *  @return object of class Void
       */
      Void () {}
      
      /**
       *  @brief constructor that uses comoving coordinates
       *  @param xx comoving coordinate
       *  @param yy comoving coordinate
       *  @param zz comoving coordinate
       *  @param weight weight
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
       *  @return object of class Void
       */
      Void (const double xx, const double yy, const double zz, const double weight=1., const double volNorm=par::defaultDouble, const double radius=par::defaultDouble, const double volume=par::defaultDouble, const int ID=par::defaultInt, const double densContr=par::defaultDouble, const int parentID=par::defaultInt, const int treeLevel=par::defaultInt, const int child=par::defaultInt, const double rho0=par::defaultDouble, const double rho0Norm=par::defaultDouble)
	: Object(xx, yy, zz, weight), m_volNorm(volNorm), m_radius(radius), m_volume(volume), m_ID(ID), m_densContr(densContr), m_parentID(parentID), m_treeLevel(treeLevel), m_child(child), m_rho0(rho0), m_rho0Norm(rho0Norm) {}

      /**
       *  @brief constructor that uses observed coordinates
       *  @param ra Right Ascension
       *  @param dec Declination
       *  @param redshift comoving coordinate
       *  @param cosm object of class Cosmology, used to estimate
       *  comoving distances
       *  @param weight weight
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
       *  @return object of class Void
       */
      Void (const double ra, const double dec, const double redshift, const Cosmology &cosm, const double weight=1., const double volNorm=par::defaultDouble, const double radius=par::defaultDouble, const double volume=par::defaultDouble, const int ID=par::defaultInt, const double densContr=par::defaultDouble, const int parentID=par::defaultInt, const int treeLevel=par::defaultInt, const int child=par::defaultInt, const double rho0=par::defaultDouble, const double rho0Norm=par::defaultDouble)
	: Object(ra, dec, redshift, cosm, weight), m_volNorm(volNorm), m_radius(radius), m_volume(volume), m_ID(ID), m_densContr(densContr), m_parentID(parentID), m_treeLevel(treeLevel), m_child(child), m_rho0(rho0), m_rho0Norm(rho0Norm) {}
    
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Void () {}

      ///@}
  
  
      /**
       *  @name Member functions used to get the private members 
       */
      ///@{
    
      /**
       *  @brief get the private member Void::m_radius
       *  @return the radius of the spherical void 
       */
      double radius () const override { return m_radius; }
    
      /**
       *  @brief get the private member Void::m_volNorm
       *  @return the normalized volume of the void
       */
      double volNorm () const { return m_volNorm; }

      /**
       *  @brief get the private member Void::m_volume
       *  @return the volume of the void 
       */
      double volume () const { return m_volume; }

      /**
       *  @brief get the private member Void::m_ID
       *  @return the identification number of the void
       */
      int ID () const { return m_ID; }

      /**
       *  @brief get the private member Void::m_densContr
       *  @return the density contrast between the void and the Universe
       */
      double densContr () const { return m_densContr; }

      /**
       *  @brief get the private member Void::m_parentID
       *  @return the identification number of the parent void
       */
      int parentID () const { return m_parentID; }

      /**
       *  @brief get the private member Void::m_treeLevel
       *  @return the hierarchy level of the void
       */
      int treeLevel () const { return m_treeLevel; }

      /**
       *  @brief get the private member Void::m_child
       *  @return the number of children of the void
       */
      int child () const { return m_child; }

      /**
       *  @brief get the private member Void::m_rho0
       *  @return the central density of the void
       */
      double rho0 () const { return m_rho0; }

      /**
       *  @brief get the private member Void::m_rho0Norm
       *  @return the normalized core density of the void
       */
      double rho0Norm () const  { return m_rho0Norm; }

      ///@}
  
  
      /**
       *  @name Member functions used to set the private members 
       */
      ///@{
    	
      /**
       *  @brief get the private member Void::m_radius
       *  @param radius the radius of the spherical void [Mpc/h]
       *  @return none
       */
      void set_radius (const double radius=par::defaultDouble) override { m_radius = radius; }
    
      /**
       *  @brief set the private member Void::m_volNorm
       *  @param volNorm the normalized volume of the void
       *  @return none
       */
      void set_volNorm (const double volNorm=par::defaultDouble) { m_volNorm = volNorm; }
	
      /**
       *  @brief get the private member Void::m_volume
       *  @param volume the volume of the void [(Mpc/h)^3]
       *  @return none
       */
      void set_volume (const double volume=par::defaultDouble) { m_volume = volume; }
	
      /**
       *  @brief get the private member Void::m_ID
       *  @param ID the identification number of the void
       *  @return none
       */
      void set_ID (const int ID=par::defaultInt) { m_ID = ID; } 
	
      /**
       *  @brief get the private member Void::m_densContr
       *  @param densContr the density contrast between the void and the Universe
       *  @return none
       */
      void set_densContr (const double densContr=par::defaultDouble) { m_densContr = densContr; }
	
      /**
       *  @brief get the private member Void::m_parentID
       *  @param parentID the identification number of the parent void
       *  @return none
       */
      void set_parentID (const int parentID=par::defaultInt) { m_parentID = parentID; }
	
      /**
       *  @brief get the private member Void::m_treeLevel
       *  @param treeLevel the hierarchy level of the void
       *  @return none
       */
      void set_treeLevel (const int treeLevel=par::defaultInt) { m_treeLevel = treeLevel; }
	
      /**
       *  @brief get the private member Void::m_child
       *  @param child the number of children of the void
       *  @return none
       */
      void set_child (const int child=par::defaultInt) { m_child = child; }
	
      /**
       *  @brief get the private member Void::m_rho0
       *  @param rho0 the central density of the void
       *  @return none
       */
      void set_rho0 (const double rho0=par::defaultDouble) { m_rho0 = rho0; }    
	
      /**
       *  @brief get the private member Void::m_rho0Norm
       *  @param rho0Norm the normalized core density of the void
       *  @return none
       */
      void set_rho0Norm (const double rho0Norm=par::defaultDouble) { m_rho0Norm = rho0Norm; }
    
      ///@}
    
    };
  }
}

#endif
