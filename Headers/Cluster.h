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
 *  @file Headers/Cluster.h
 *
 *  @brief The class Cluster 
 *
 *  This file defines the interface of the class Cluster, used to
 *  handle objects of type <EM> cluster of galaxies </EM>
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __CLUSTER__
#define __CLUSTER__ 


// ============================================================================================


namespace cbl {
  
  namespace catalogue {

    /**
     *  @class Cluster Cluster.h "Headers/Cluster.h"
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
    
      /// cluster mass proxy
      double m_mass_proxy;

      /// cluster proxy error
      double m_mass_proxy_error;
      
      /// cluster linear bias
      double m_bias;
      
    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  
       */
      Cluster ()
	: Object(), m_mass(par::defaultDouble), m_mass_proxy(par::defaultDouble), m_mass_proxy_error(par::defaultDouble), m_bias(par::defaultDouble) {}

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
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param redshiftMin minimum redshift
       *
       *  @param redshiftMax maximum redshift
       *
       *  @param sn signal-to-noise
       *
       *  @param mass the cluster mass
       *
       *  @param mass_proxy the cluster mass proxy
       *
       *  @param mass_proxy_error the cluster mass proxy error
       *
       *  @param bias the cluster linear bias
       *
       *  
       */
      Cluster (const comovingCoordinates coord, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double mass_proxy=par::defaultDouble, const double mass_proxy_error=par::defaultDouble, const double bias=par::defaultDouble) 
	: Object(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_mass_proxy(mass_proxy), m_mass_proxy_error(mass_proxy_error), m_bias(bias) {}

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
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param redshiftMin minimum redshift
       *
       *  @param redshiftMax maximum redshift
       *
       *  @param sn signal-to-noise
       *   
       *  @param mass the cluster mass
       *
       *  @param mass_proxy the cluster mass proxy
       *
       *  @param mass_proxy_error the error on the cluster mass proxy
       *
       *  @param bias the cluster linear bias
       *
       *  
       */
      Cluster (const comovingCoordinates coord, const cosmology::Cosmology &cosm, const double z1_guess=0., const double z2_guess=10., const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double mass_proxy=par::defaultDouble, const double mass_proxy_error=par::defaultDouble, const double bias=par::defaultDouble) 
	: Object(coord, cosm, z1_guess, z2_guess, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_mass_proxy(mass_proxy), m_mass_proxy_error(mass_proxy_error), m_bias(bias) {}

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
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param redshiftMin minimum redshift
       *
       *  @param redshiftMax maximum redshift
       *
       *  @param sn signal-to-noise
       *
       *  @param mass the cluster mass
       *
       *  @param mass_proxy the cluster mass proxy
       *
       *  @param mass_proxy_error the error on the cluster mass proxy
       *
       *  @param bias the cluster linear bias
       *
       *  
       */
      Cluster (const observedCoordinates coord, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double mass_proxy=par::defaultDouble, const double mass_proxy_error=par::defaultDouble, const double bias=par::defaultDouble) 
	: Object(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_mass_proxy(mass_proxy), m_mass_proxy_error(mass_proxy_error), m_bias(bias) {}
      
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
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param redshiftMin minimum redshift
       *
       *  @param redshiftMax maximum redshift
       *
       *  @param sn signal-to-noise
       *
       *  @param mass the cluster mass
       *
       *  @param mass_proxy the cluster mass proxy
       *
       *  @param mass_proxy_error the error on the cluster mass proxy
       *
       *  @param bias the cluster linear bias
       *
       *  
       */
      Cluster (const observedCoordinates coord, const CoordinateUnits inputUnits, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double mass_proxy=par::defaultDouble, const double mass_proxy_error=par::defaultDouble, const double bias=par::defaultDouble) 
	: Object(coord, inputUnits, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_mass_proxy(mass_proxy), m_mass_proxy_error(mass_proxy_error), m_bias(bias) {}
      
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
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param redshiftMin minimum redshift
       *
       *  @param redshiftMax maximum redshift
       *
       *  @param sn signal-to-noise
       *
       *  @param mass the cluster mass
       *
       *  @param mass_proxy the cluster mass proxy
       *
       *  @param mass_proxy_error the error on the cluster mass proxy
       *
       *  @param bias the cluster linear bias
       *
       *  
       */
      Cluster (const observedCoordinates coord, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double mass_proxy=par::defaultDouble, const double mass_proxy_error=par::defaultDouble, const double bias=par::defaultDouble) 
	: Object(coord, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_mass_proxy(mass_proxy), m_mass_proxy_error(mass_proxy_error), m_bias(bias) {}

      /**
       *  @brief constructor that uses observed coordinates and a
       *  cosmological model to estimate the comoving coordinates
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshift}
       *
       *  @param inputUnits the units of the input coordinates
       *
       *  @param cosm object of class Cosmology, used to estimate
       *  comoving distances
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param redshiftMin minimum redshift
       *
       *  @param redshiftMax maximum redshift
       *
       *  @param sn signal-to-noise
       *
       *  @param mass the cluster mass
       *
       *  @param mass_proxy the cluster mass proxy
       *
       *  @param mass_proxy_error the error on the cluster mass proxy
       *
       *  @param bias the cluster linear bias
       *
       *  
       */
      Cluster (const observedCoordinates coord, const CoordinateUnits inputUnits, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double mass_proxy=par::defaultDouble, const double mass_proxy_error=par::defaultDouble, const double bias=par::defaultDouble) 
	: Object(coord, inputUnits, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_mass_proxy(mass_proxy), m_mass_proxy_error(mass_proxy_error), m_bias(bias) {}

      /**
       *  @brief constructor that uses both comoving and observed
       *  coordinates
       *
       *  @param xx comoving coordinate
       *
       *  @param yy comoving coordinate
       *
       *  @param zz comoving coordinate 
       *
       *  @param ra Right Ascension
       *
       *  @param dec Declination
       *
       *  @param redshift redshift
       *
       *  @param weight weight  
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param ID the object ID
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param redshiftMin minimum redshift
       *
       *  @param redshiftMax maximum redshift
       *
       *  @param sn signal-to-noise
       *
       *  @param mass the cluster mass
       *
       *  @param mass_proxy the cluster mass proxy
       *
       *  @param mass_proxy_error the error on the cluster mass proxy
       *
       *  @param bias the cluster bias
       *
       *  
       */
      Cluster (const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double mass_proxy=par::defaultDouble, const double mass_proxy_error=par::defaultDouble, const double bias=par::defaultDouble) 
	: Object(xx, yy, zz, ra, dec, redshift, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_mass_proxy(mass_proxy), m_mass_proxy_error(mass_proxy_error), m_bias(bias) {}
      
      /**
       *  @brief default destructor
       */
      ~Cluster () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to get the private members 
       */
      ///@{
    
      /**
       *  @brief get the private member \e m_mass
       *  @return the mass of the cluster
       */
      double mass () const override
      { return m_mass; }

      /**
       *  @brief get the private member \e m_mass_proxy
       *  @return the mass proxy of the cluster
       */
      double mass_proxy () const override
      { return m_mass_proxy; }

      /**
       *  @brief get the private member \e m_mass_proxy_error
       *  @return the mass proxy error of the cluster
       */
      double mass_proxy_error () const override
      { return m_mass_proxy_error; }

      /**
       *  @brief get the private member \e m_bias
       *  @return the linear bias of the cluster
       */
      double bias () const override
      { return m_bias; }

      ///@}
  
  
      /**
       *  @name Member functions used to set the private members 
       */
      ///@{
    
      /**
       *  @brief set the private member \e m_mass
       *  @param mass the mass of the cluster
       */
      void set_mass (const double mass=par::defaultDouble) override
      { m_mass = mass; }

      /**
       *  @brief set the private member \e m_mass_proxy
       *  @param mass_proxy the mass proxy of the cluster
       */
      void set_mass_proxy (const double mass_proxy=par::defaultDouble) override
      { m_mass_proxy = mass_proxy; }

      /**
       *  @brief set the private member \e m_mass_proxy_error
       *  @param mass_proxy_error the mass proxy of the cluster
       */
      void set_mass_proxy_error (const double mass_proxy_error=par::defaultDouble) override
      { m_mass_proxy_error = mass_proxy_error; }

      /**
       *  @brief set the private member \e m_bias
       *  @param bias the linear bias of the cluster
       */
      void set_bias (const double bias=par::defaultDouble) override
      { m_bias = bias; }

      ///@}


      /**
       *  @name Member functions used to check if the private members are set 
       */
      ///@{
    
      /**
       *  @brief check if the private member \e m_mass is set
       *
       *  @return true if the mass is set; false otherwise
       */
      bool isSet_mass () override
      { return (cbl::isSet(m_mass)) ? true : false; }
      
      /**
       *  @brief check if the private member \e m_mass_proxy is set
       *  
       *  @return true if the proxy is set; false otherwise
       */
      bool isSet_mass_proxy () override
      { return (cbl::isSet(m_mass_proxy)) ? true : false; }

      /**
       *  @brief check if the private member \e m_mass_proxy_error is set
       *  
       *  @return true if the error on the mass proxy error is set; false
       *  otherwise
       */
      bool isSet_mass_proxy_error () override
      { return (cbl::isSet(m_mass_proxy_error)) ? true : false; }

      /**
       *  @brief check if the private member \e m_bias is set
       *  
       *  @return true if the bias is set; false otherwise
       */
      bool isSet_bias () override
      { return (cbl::isSet(m_bias)) ? true : false; }

      ///@}
      
    
    };
  }
}

#endif

