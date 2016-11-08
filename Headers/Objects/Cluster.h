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
      Cluster ()
	: Object(), m_mass(par::defaultDouble), m_richness(par::defaultDouble) {}

      /**
       *  @brief constructor that uses comoving coordinates
       *
       *  @param coord structure containing the comoving coordinates
       *  {x, y, z}
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param mass the cluster mass
       *
       *  @param richness the cluster richness
       *
       *  @return object of class Cluster
       */
      Cluster (const comovingCoordinates coord, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double mass=par::defaultDouble, const double richness=par::defaultDouble) 
	: Object(coord, weight, region, field, x_displacement, y_displacement, z_displacement), m_mass(mass), m_richness(richness) {}

      /**
       *  @brief constructor that uses comoving coordinates and a
       *  cosmological model to estimate the redshift
       *
       *  @param coord structure containing the comoving coordinates
       *  {x, y, z}
       *
       *  @param cosm object of class Cosmology, used to estimate
       *  comoving distances
       *
       *  @param z1_guess minimum prior on the redshift
       *
       *  @param z2_guess maximum prior on the redshift 
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *   
       *  @param mass the cluster mass
       *
       *  @param richness the cluster richness
       *
       *  @return object of class Cluster
       */
      Cluster (const comovingCoordinates coord, const cosmology::Cosmology &cosm, const double z1_guess=0., const double z2_guess=10., const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double mass=par::defaultDouble, const double richness=par::defaultDouble) 
	: Object(coord, cosm, z1_guess, z2_guess, weight, region, field, x_displacement, y_displacement, z_displacement), m_mass(mass), m_richness(richness) {}

      /**
       *  @brief constructor that uses observed coordinates in radians
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshift}
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param mass the cluster mass
       *
       *  @param richness the cluster richness
       *
       *  @return object of class Cluster
       */
      Cluster (const observedCoordinates coord, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double mass=par::defaultDouble, const double richness=par::defaultDouble) 
	: Object(coord, weight, region, field, x_displacement, y_displacement, z_displacement), m_mass(mass), m_richness(richness) {}
      
      /**
       *  @brief constructor that uses observed coordinates in any
       *  angular units
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshift}
       *
       *  @param inputUnits the units of the input coordinates
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param mass the cluster mass
       *
       *  @param richness the cluster richness
       *
       *  @return object of class Cluster
       */
      Cluster (const observedCoordinates coord, const CoordUnits inputUnits, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double mass=par::defaultDouble, const double richness=par::defaultDouble) 
	: Object(coord, inputUnits, weight, region, field, x_displacement, y_displacement, z_displacement), m_mass(mass), m_richness(richness) {}
      
      /**
       *  @brief constructor that uses observed coordinates in radians
       *  and a cosmological model to estimate the comoving
       *  coordinates
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshitf}
       *
       *  @param cosm object of class Cosmology, used to estimate
       *  comoving distances
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param mass the cluster mass
       *
       *  @param richness the cluster richness
       *
       *  @return object of class Cluster
       */
      Cluster (const observedCoordinates coord, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double mass=par::defaultDouble, const double richness=par::defaultDouble) 
	: Object(coord, cosm, weight, region, field, x_displacement, y_displacement, z_displacement), m_mass(mass), m_richness(richness) {}

      /**
       *  @brief constructor that uses observed coordinates and a
       *  cosmological model to estimate the comoving coordinates
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshift}
       *
       *  @param inputUnits the units of the input coordinates
       *
       *  @param cosm object of class Cosmology, used to estimate comoving distances
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param mass the cluster mass
       *
       *  @param richness the cluster richness
       *
       *  @return object of class Cluster
       */
      Cluster (const observedCoordinates coord, const CoordUnits inputUnits, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double mass=par::defaultDouble, const double richness=par::defaultDouble) 
	: Object(coord, inputUnits, cosm, weight, region, field, x_displacement, y_displacement, z_displacement), m_mass(mass), m_richness(richness) {}

      /**
       *  @brief constructor that uses both comoving and observed coordinates
       *  @param xx comoving coordinate
       *  @param yy comoving coordinate
       *  @param zz comoving coordinate 
       *  @param ra Right Ascension
       *  @param dec Declination
       *  @param redshift redshift
       *  @param weight weight  
       *  @param region region, used e.g. for jackknife and bootstrap
       *  @param field the field where the object has been observed
       *  @param x_displacement the displacement along the x-axis
       *  @param y_displacement the displacement along the y-axis
       *  @param z_displacement the displacement along the z-axis
       *  @param mass the cluster mass
       *  @param richness the cluster richness
       *
       *  @return object of class Cluster
       */
      Cluster (const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double mass=par::defaultDouble, const double richness=par::defaultDouble) 
	: Object(xx, yy, zz, ra, dec, redshift, weight, region, field, x_displacement, y_displacement, z_displacement), m_mass(mass), m_richness(richness) {}
      
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Cluster () = default;

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

