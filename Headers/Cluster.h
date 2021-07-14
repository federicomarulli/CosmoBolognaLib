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
      
      /// cluster mass logarithm (undefined base here)
      double m_logM; 
    
      /// cluster mass proxy
      double m_mass_proxy;

      /// cluster proxy error
      double m_mass_proxy_error;
      
      /// cluster linear bias
      double m_bias;
      
      /// cluster concentration
      double m_concentration;
      
      /// fraction of miscentered cluster population
      double m_f_off;
      
      /// rms of the miscentered distribution
      double m_sigma_off;
      
      /// normalization of the mass-observable scaling relation
      double m_alpha_scaling_rel;
      
      /// slope of the mass-observable scaling relation
      double m_beta_scaling_rel;
      
      /// z evolution factor of the mass-observable scaling relation
      double m_gamma_scaling_rel;
      
      /// constant term of the intrinsic scatter of the mass-observable scaling relation
      double m_scatter0_scaling_rel;
      
      /// multiplicative factor in the mass/mass proxy dependent term in the intrinsic scatter of the mass-observable scaling relation
      double m_scatterM_scaling_rel;
      
      /// exponent in the mass/mass proxy dependent term in the intrinsic scatter of the mass-observable scaling relation
      double m_scatterM_exponent_scaling_rel;
      
      /// multiplicative factor in the redshift dependent term in the intrinsic scatter of the mass-observable scaling relation
      double m_scatterz_scaling_rel;
      
      /// exponent in the redshift dependent term in the intrinsic scatter of the mass-observable scaling relation
      double m_scatterz_exponent_scaling_rel;
      
      /// cluster redshift bias
      double m_zbias;
      
      /// cluster mass proxy bias
      double m_proxybias;
      
      /// cluster redshift uncertainty
      double m_zerror;
      
      /// cluster mass proxy uncertainty
      double m_proxyerror;
      
      /// \f$ a \f$ term in the function describing the cluster abundance, i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
      double m_Plambda_a;
      
      /// \f$ b \f$ term in the function describing the cluster abundance, i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
      double m_Plambda_b;
      
      /// \f$ c \f$ term in the function describing the cluster abundance, i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
      double m_Plambda_c;
      
      /// 
      
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
	: Object(), m_mass(par::defaultDouble), m_logM(par::defaultDouble), m_mass_proxy(par::defaultDouble), m_mass_proxy_error(par::defaultDouble), m_bias(par::defaultDouble), m_concentration(par::defaultDouble), m_f_off(par::defaultDouble), m_sigma_off(par::defaultDouble), m_alpha_scaling_rel(par::defaultDouble), m_beta_scaling_rel(par::defaultDouble), m_gamma_scaling_rel(par::defaultDouble), m_scatter0_scaling_rel(par::defaultDouble), m_scatterM_scaling_rel(par::defaultDouble), m_scatterM_exponent_scaling_rel(par::defaultDouble), m_scatterz_scaling_rel(par::defaultDouble), m_scatterz_exponent_scaling_rel(par::defaultDouble), m_zbias(par::defaultDouble), m_proxybias(par::defaultDouble), m_zerror(par::defaultDouble), m_proxyerror(par::defaultDouble), m_Plambda_a(par::defaultDouble), m_Plambda_b(par::defaultDouble), m_Plambda_c(par::defaultDouble) {}

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
       *  @param bias cluster linear bias
       *
       *  @param logM the cluster mass logarithm
       *
       *  @param conc concentration
       *
       *  @param f_off fraction of miscentered cluster population
       * 
       *  @param sigma_off rms of the miscentered distribution
       *
       *  @param alpha_scaling_rel normalization of the mass-observable scaling relation
       *
       *  @param beta_scaling_rel slope of the mass-observable scaling relation
       *
       *  @param gamma_scaling_rel z evolution factor of the mass-observable scaling relation
       *
       *  @param scatter0_scaling_rel constant term of the intrinsic scatter 
       *  of the mass-observable scaling relation
       *
       *  @param scatterM_scaling_rel multiplicative factor in the mass/mass proxy 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterM_exponent_scaling_rel exponent in the mass/mass proxy 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterz_scaling_rel multiplicative factor in the redshift
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterz_exponent_scaling_rel exponent in the redshift 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param zbias cluster redshift bias
       *
       *  @param proxybias cluster mass proxy bias
       *
       *  @param zerror cluster redshift error
       *
       *  @param proxyerror cluster redshift bias
       *
       *  @param Plambda_a \f$ a \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  @param Plambda_b \f$ b \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  @param Plambda_c \f$ c \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *  
       */
      Cluster (const comovingCoordinates coord, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double mass_proxy=par::defaultDouble, const double mass_proxy_error=par::defaultDouble, const double bias=par::defaultDouble, const double logM=par::defaultDouble, const double conc=par::defaultDouble, const double f_off=par::defaultDouble, const double sigma_off=par::defaultDouble, const double alpha_scaling_rel=par::defaultDouble, const double beta_scaling_rel=par::defaultDouble, const double gamma_scaling_rel=par::defaultDouble, const double scatter0_scaling_rel=par::defaultDouble, const double scatterM_scaling_rel=par::defaultDouble, const double scatterM_exponent_scaling_rel=par::defaultDouble, const double scatterz_scaling_rel=par::defaultDouble, const double scatterz_exponent_scaling_rel=par::defaultDouble, const double zbias=par::defaultDouble, const double proxybias=par::defaultDouble, const double zerror=par::defaultDouble, const double proxyerror=par::defaultDouble, const double Plambda_a=par::defaultDouble, const double Plambda_b=par::defaultDouble, const double Plambda_c=par::defaultDouble) 
      : Object(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_logM(logM), m_mass_proxy(mass_proxy), m_mass_proxy_error(mass_proxy_error), m_bias(bias), m_concentration(conc), m_f_off(f_off), m_sigma_off(sigma_off), m_alpha_scaling_rel(alpha_scaling_rel), m_beta_scaling_rel(beta_scaling_rel), m_gamma_scaling_rel(gamma_scaling_rel), m_scatter0_scaling_rel(scatter0_scaling_rel), m_scatterM_scaling_rel(scatterM_scaling_rel), m_scatterM_exponent_scaling_rel(scatterM_exponent_scaling_rel), m_scatterz_scaling_rel(scatterz_scaling_rel), m_scatterz_exponent_scaling_rel(scatterz_exponent_scaling_rel), m_zbias(zbias), m_proxybias(proxybias), m_zerror(zerror), m_proxyerror(proxyerror), m_Plambda_a(Plambda_a), m_Plambda_b(Plambda_b), m_Plambda_c(Plambda_c) {}

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
       *  @param logM the cluster mass logarithm
       *
       *  @param conc concentration
       *
       *  @param f_off fraction of miscentered cluster population
       * 
       *  @param sigma_off rms of the miscentered distribution
       *
       *  @param alpha_scaling_rel normalization of the mass-observable scaling relation
       *
       *  @param beta_scaling_rel slope of the mass-observable scaling relation
       *
       *  @param gamma_scaling_rel z evolution factor of the mass-observable scaling relation
       *
       *  @param scatter0_scaling_rel constant term of the intrinsic scatter 
       *  of the mass-observable scaling relation
       *
       *  @param scatterM_scaling_rel multiplicative factor in the mass/mass proxy 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterM_exponent_scaling_rel exponent in the mass/mass proxy 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterz_scaling_rel multiplicative factor in the redshift
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterz_exponent_scaling_rel exponent in the redshift 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param zbias cluster redshift bias
       *
       *  @param proxybias cluster mass proxy bias
       *
       *  @param zerror cluster redshift error
       *
       *  @param proxyerror cluster redshift bias
       *
       *  @param Plambda_a \f$ a \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  @param Plambda_b \f$ b \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  @param Plambda_c \f$ c \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  
       */
      Cluster (const comovingCoordinates coord, const cosmology::Cosmology &cosm, const double z1_guess=0., const double z2_guess=10., const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double mass_proxy=par::defaultDouble, const double mass_proxy_error=par::defaultDouble, const double bias=par::defaultDouble, const double logM=par::defaultDouble, const double conc=par::defaultDouble, const double f_off=par::defaultDouble, const double sigma_off=par::defaultDouble, const double alpha_scaling_rel=par::defaultDouble, const double beta_scaling_rel=par::defaultDouble, const double gamma_scaling_rel=par::defaultDouble, const double scatter0_scaling_rel=par::defaultDouble, const double scatterM_scaling_rel=par::defaultDouble, const double scatterM_exponent_scaling_rel=par::defaultDouble, const double scatterz_scaling_rel=par::defaultDouble, const double scatterz_exponent_scaling_rel=par::defaultDouble, const double zbias=par::defaultDouble, const double proxybias=par::defaultDouble, const double zerror=par::defaultDouble, const double proxyerror=par::defaultDouble, const double Plambda_a=par::defaultDouble, const double Plambda_b=par::defaultDouble, const double Plambda_c=par::defaultDouble) 
      : Object(coord, cosm, z1_guess, z2_guess, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_logM(logM), m_mass_proxy(mass_proxy), m_mass_proxy_error(mass_proxy_error), m_bias(bias), m_concentration(conc), m_f_off(f_off), m_sigma_off(sigma_off), m_alpha_scaling_rel(alpha_scaling_rel), m_beta_scaling_rel(beta_scaling_rel), m_gamma_scaling_rel(gamma_scaling_rel), m_scatter0_scaling_rel(scatter0_scaling_rel), m_scatterM_scaling_rel(scatterM_scaling_rel), m_scatterM_exponent_scaling_rel(scatterM_exponent_scaling_rel), m_scatterz_scaling_rel(scatterz_scaling_rel), m_scatterz_exponent_scaling_rel(scatterz_exponent_scaling_rel), m_zbias(zbias), m_proxybias(proxybias), m_zerror(zerror), m_proxyerror(proxyerror), m_Plambda_a(Plambda_a), m_Plambda_b(Plambda_b), m_Plambda_c(Plambda_c) {}

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
       *  @param logM the cluster mass logarithm
       *
       *  @param conc concentration
       *
       *  @param f_off fraction of miscentered cluster population
       * 
       *  @param sigma_off rms of the miscentered distribution
       *
       *  @param alpha_scaling_rel normalization of the mass-observable scaling relation
       *
       *  @param beta_scaling_rel slope of the mass-observable scaling relation
       *
       *  @param gamma_scaling_rel z evolution factor of the mass-observable scaling relation
       *
       *  @param scatter0_scaling_rel constant term of the intrinsic scatter 
       *  of the mass-observable scaling relation
       *
       *  @param scatterM_scaling_rel multiplicative factor in the mass/mass proxy 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterM_exponent_scaling_rel exponent in the mass/mass proxy 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterz_scaling_rel multiplicative factor in the redshift
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterz_exponent_scaling_rel exponent in the redshift 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param zbias cluster redshift bias
       *
       *  @param proxybias cluster mass proxy bias
       *
       *  @param zerror cluster redshift error
       *
       *  @param proxyerror cluster redshift bias
       *
       *  @param Plambda_a \f$ a \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  @param Plambda_b \f$ b \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  @param Plambda_c \f$ c \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  
       */
      Cluster (const observedCoordinates coord, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double mass_proxy=par::defaultDouble, const double mass_proxy_error=par::defaultDouble, const double bias=par::defaultDouble, const double logM=par::defaultDouble, const double conc=par::defaultDouble, const double f_off=par::defaultDouble, const double sigma_off=par::defaultDouble, const double alpha_scaling_rel=par::defaultDouble, const double beta_scaling_rel=par::defaultDouble, const double gamma_scaling_rel=par::defaultDouble, const double scatter0_scaling_rel=par::defaultDouble, const double scatterM_scaling_rel=par::defaultDouble, const double scatterM_exponent_scaling_rel=par::defaultDouble, const double scatterz_scaling_rel=par::defaultDouble, const double scatterz_exponent_scaling_rel=par::defaultDouble, const double zbias=par::defaultDouble, const double proxybias=par::defaultDouble, const double zerror=par::defaultDouble, const double proxyerror=par::defaultDouble, const double Plambda_a=par::defaultDouble, const double Plambda_b=par::defaultDouble, const double Plambda_c=par::defaultDouble) 
      : Object(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_logM(logM), m_mass_proxy(mass_proxy), m_mass_proxy_error(mass_proxy_error), m_bias(bias), m_concentration(conc), m_f_off(f_off), m_sigma_off(sigma_off), m_alpha_scaling_rel(alpha_scaling_rel), m_beta_scaling_rel(beta_scaling_rel), m_gamma_scaling_rel(gamma_scaling_rel), m_scatter0_scaling_rel(scatter0_scaling_rel), m_scatterM_scaling_rel(scatterM_scaling_rel), m_scatterM_exponent_scaling_rel(scatterM_exponent_scaling_rel), m_scatterz_scaling_rel(scatterz_scaling_rel), m_scatterz_exponent_scaling_rel(scatterz_exponent_scaling_rel), m_zbias(zbias), m_proxybias(proxybias), m_zerror(zerror), m_proxyerror(proxyerror), m_Plambda_a(Plambda_a), m_Plambda_b(Plambda_b), m_Plambda_c(Plambda_c) {}
      
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
       *  @param logM the cluster mass logarithm
       *
       *  @param conc concentration
       *
       *  @param f_off fraction of miscentered cluster population
       * 
       *  @param sigma_off rms of the miscentered distribution
       *
       *  @param alpha_scaling_rel normalization of the mass-observable scaling relation
       *
       *  @param beta_scaling_rel slope of the mass-observable scaling relation
       *
       *  @param gamma_scaling_rel z evolution factor of the mass-observable scaling relation
       *
       *  @param scatter0_scaling_rel constant term of the intrinsic scatter 
       *  of the mass-observable scaling relation
       *
       *  @param scatterM_scaling_rel multiplicative factor in the mass/mass proxy 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterM_exponent_scaling_rel exponent in the mass/mass proxy 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterz_scaling_rel multiplicative factor in the redshift
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterz_exponent_scaling_rel exponent in the redshift 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param zbias cluster redshift bias
       *
       *  @param proxybias cluster mass proxy bias
       *
       *  @param zerror cluster redshift error
       *
       *  @param proxyerror cluster redshift bias
       *
       *  @param Plambda_a \f$ a \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  @param Plambda_b \f$ b \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  @param Plambda_c \f$ c \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  
       */
      Cluster (const observedCoordinates coord, const CoordinateUnits inputUnits, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double mass_proxy=par::defaultDouble, const double mass_proxy_error=par::defaultDouble, const double bias=par::defaultDouble, const double logM=par::defaultDouble, const double conc=par::defaultDouble, const double f_off=par::defaultDouble, const double sigma_off=par::defaultDouble, const double alpha_scaling_rel=par::defaultDouble, const double beta_scaling_rel=par::defaultDouble, const double gamma_scaling_rel=par::defaultDouble, const double scatter0_scaling_rel=par::defaultDouble, const double scatterM_scaling_rel=par::defaultDouble, const double scatterM_exponent_scaling_rel=par::defaultDouble, const double scatterz_scaling_rel=par::defaultDouble, const double scatterz_exponent_scaling_rel=par::defaultDouble, const double zbias=par::defaultDouble, const double proxybias=par::defaultDouble, const double zerror=par::defaultDouble, const double proxyerror=par::defaultDouble, const double Plambda_a=par::defaultDouble, const double Plambda_b=par::defaultDouble, const double Plambda_c=par::defaultDouble) 
      : Object(coord, inputUnits, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_logM(logM), m_mass_proxy(mass_proxy), m_mass_proxy_error(mass_proxy_error), m_bias(bias), m_concentration(conc), m_f_off(f_off), m_sigma_off(sigma_off), m_alpha_scaling_rel(alpha_scaling_rel), m_beta_scaling_rel(beta_scaling_rel), m_gamma_scaling_rel(gamma_scaling_rel), m_scatter0_scaling_rel(scatter0_scaling_rel), m_scatterM_scaling_rel(scatterM_scaling_rel), m_scatterM_exponent_scaling_rel(scatterM_exponent_scaling_rel), m_scatterz_scaling_rel(scatterz_scaling_rel), m_scatterz_exponent_scaling_rel(scatterz_exponent_scaling_rel), m_zbias(zbias), m_proxybias(proxybias), m_zerror(zerror), m_proxyerror(proxyerror), m_Plambda_a(Plambda_a), m_Plambda_b(Plambda_b), m_Plambda_c(Plambda_c) {}
      
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
       *  @param logM the cluster mass logarithm
       *
       *  @param conc concentration
       *
       *  @param f_off fraction of miscentered cluster population
       * 
       *  @param sigma_off rms of the miscentered distribution
       *
       *  @param alpha_scaling_rel normalization of the mass-observable scaling relation
       *
       *  @param beta_scaling_rel slope of the mass-observable scaling relation
       *
       *  @param gamma_scaling_rel z evolution factor of the mass-observable scaling relation
       *
       *  @param scatter0_scaling_rel constant term of the intrinsic scatter 
       *  of the mass-observable scaling relation
       *
       *  @param scatterM_scaling_rel multiplicative factor in the mass/mass proxy 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterM_exponent_scaling_rel exponent in the mass/mass proxy 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterz_scaling_rel multiplicative factor in the redshift
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterz_exponent_scaling_rel exponent in the redshift 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param zbias cluster redshift bias
       *
       *  @param proxybias cluster mass proxy bias
       *
       *  @param zerror cluster redshift error
       *
       *  @param proxyerror cluster redshift bias
       *
       *  @param Plambda_a \f$ a \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  @param Plambda_b \f$ b \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  @param Plambda_c \f$ c \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  
       */
      Cluster (const observedCoordinates coord, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double mass_proxy=par::defaultDouble, const double mass_proxy_error=par::defaultDouble, const double bias=par::defaultDouble, const double logM=par::defaultDouble, const double conc=par::defaultDouble, const double f_off=par::defaultDouble, const double sigma_off=par::defaultDouble, const double alpha_scaling_rel=par::defaultDouble, const double beta_scaling_rel=par::defaultDouble, const double gamma_scaling_rel=par::defaultDouble, const double scatter0_scaling_rel=par::defaultDouble, const double scatterM_scaling_rel=par::defaultDouble, const double scatterM_exponent_scaling_rel=par::defaultDouble, const double scatterz_scaling_rel=par::defaultDouble, const double scatterz_exponent_scaling_rel=par::defaultDouble, const double zbias=par::defaultDouble, const double proxybias=par::defaultDouble, const double zerror=par::defaultDouble, const double proxyerror=par::defaultDouble, const double Plambda_a=par::defaultDouble, const double Plambda_b=par::defaultDouble, const double Plambda_c=par::defaultDouble) 
      : Object(coord, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_logM(logM), m_mass_proxy(mass_proxy), m_mass_proxy_error(mass_proxy_error), m_bias(bias), m_concentration(conc), m_f_off(f_off), m_sigma_off(sigma_off), m_alpha_scaling_rel(alpha_scaling_rel), m_beta_scaling_rel(beta_scaling_rel), m_gamma_scaling_rel(gamma_scaling_rel), m_scatter0_scaling_rel(scatter0_scaling_rel), m_scatterM_scaling_rel(scatterM_scaling_rel), m_scatterM_exponent_scaling_rel(scatterM_exponent_scaling_rel), m_scatterz_scaling_rel(scatterz_scaling_rel), m_scatterz_exponent_scaling_rel(scatterz_exponent_scaling_rel), m_zbias(zbias), m_proxybias(proxybias), m_zerror(zerror), m_proxyerror(proxyerror), m_Plambda_a(Plambda_a), m_Plambda_b(Plambda_b), m_Plambda_c(Plambda_c) {}

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
       *  @param logM the cluster mass logarithm
       *
       *  @param conc concentration
       *
       *  @param f_off fraction of miscentered cluster population
       * 
       *  @param sigma_off rms of the miscentered distribution
       *
       *  @param alpha_scaling_rel normalization of the mass-observable scaling relation
       *
       *  @param beta_scaling_rel slope of the mass-observable scaling relation
       *
       *  @param gamma_scaling_rel z evolution factor of the mass-observable scaling relation
       *
       *  @param scatter0_scaling_rel constant term of the intrinsic scatter 
       *  of the mass-observable scaling relation
       *
       *  @param scatterM_scaling_rel multiplicative factor in the mass/mass proxy 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterM_exponent_scaling_rel exponent in the mass/mass proxy 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterz_scaling_rel multiplicative factor in the redshift
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterz_exponent_scaling_rel exponent in the redshift 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param zbias cluster redshift bias
       *
       *  @param proxybias cluster mass proxy bias
       *
       *  @param zerror cluster redshift error
       *
       *  @param proxyerror cluster redshift bias
       *
       *  @param Plambda_a \f$ a \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  @param Plambda_b \f$ b \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  @param Plambda_c \f$ c \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  
       */
      Cluster (const observedCoordinates coord, const CoordinateUnits inputUnits, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double mass_proxy=par::defaultDouble, const double mass_proxy_error=par::defaultDouble, const double bias=par::defaultDouble, const double logM=par::defaultDouble, const double conc=par::defaultDouble, const double f_off=par::defaultDouble, const double sigma_off=par::defaultDouble, const double alpha_scaling_rel=par::defaultDouble, const double beta_scaling_rel=par::defaultDouble, const double gamma_scaling_rel=par::defaultDouble, const double scatter0_scaling_rel=par::defaultDouble, const double scatterM_scaling_rel=par::defaultDouble, const double scatterM_exponent_scaling_rel=par::defaultDouble, const double scatterz_scaling_rel=par::defaultDouble, const double scatterz_exponent_scaling_rel=par::defaultDouble, const double zbias=par::defaultDouble, const double proxybias=par::defaultDouble, const double zerror=par::defaultDouble, const double proxyerror=par::defaultDouble, const double Plambda_a=par::defaultDouble, const double Plambda_b=par::defaultDouble, const double Plambda_c=par::defaultDouble) 
      : Object(coord, inputUnits, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_logM(logM), m_mass_proxy(mass_proxy), m_mass_proxy_error(mass_proxy_error), m_bias(bias), m_concentration(conc), m_f_off(f_off), m_sigma_off(sigma_off), m_alpha_scaling_rel(alpha_scaling_rel), m_beta_scaling_rel(beta_scaling_rel), m_gamma_scaling_rel(gamma_scaling_rel), m_scatter0_scaling_rel(scatter0_scaling_rel), m_scatterM_scaling_rel(scatterM_scaling_rel), m_scatterM_exponent_scaling_rel(scatterM_exponent_scaling_rel), m_scatterz_scaling_rel(scatterz_scaling_rel), m_scatterz_exponent_scaling_rel(scatterz_exponent_scaling_rel), m_zbias(zbias), m_proxybias(proxybias), m_zerror(zerror), m_proxyerror(proxyerror), m_Plambda_a(Plambda_a), m_Plambda_b(Plambda_b), m_Plambda_c(Plambda_c) {}

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
       *  @param logM the cluster mass logarithm
       *
       *  @param conc concentration
       *
       *  @param f_off fraction of miscentered cluster population
       * 
       *  @param sigma_off rms of the miscentered distribution
       *
       *  @param alpha_scaling_rel normalization of the mass-observable scaling relation
       *
       *  @param beta_scaling_rel slope of the mass-observable scaling relation
       *
       *  @param gamma_scaling_rel z evolution factor of the mass-observable scaling relation
       *
       *  @param scatter0_scaling_rel constant term of the intrinsic scatter 
       *  of the mass-observable scaling relation
       *
       *  @param scatterM_scaling_rel multiplicative factor in the mass/mass proxy 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterM_exponent_scaling_rel exponent in the mass/mass proxy 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterz_scaling_rel multiplicative factor in the redshift
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param scatterz_exponent_scaling_rel exponent in the redshift 
       *  dependent term in the intrinsic scatter of the mass-observable scaling relation
       *
       *  @param zbias cluster redshift bias
       *
       *  @param proxybias cluster mass proxy bias
       *
       *  @param zerror cluster redshift error
       *
       *  @param proxyerror cluster redshift bias
       *
       *  @param Plambda_a \f$ a \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  @param Plambda_b \f$ b \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  @param Plambda_c \f$ c \f$ term in the function describing the cluster abundance, 
       *  i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
       *
       *  
       */
      Cluster (const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double mass_proxy=par::defaultDouble, const double mass_proxy_error=par::defaultDouble, const double bias=par::defaultDouble, const double logM=par::defaultDouble, const double conc=par::defaultDouble, const double f_off=par::defaultDouble, const double sigma_off=par::defaultDouble, const double alpha_scaling_rel=par::defaultDouble, const double beta_scaling_rel=par::defaultDouble, const double gamma_scaling_rel=par::defaultDouble, const double scatter0_scaling_rel=par::defaultDouble, const double scatterM_scaling_rel=par::defaultDouble, const double scatterM_exponent_scaling_rel=par::defaultDouble, const double scatterz_scaling_rel=par::defaultDouble, const double scatterz_exponent_scaling_rel=par::defaultDouble, const double zbias=par::defaultDouble, const double proxybias=par::defaultDouble, const double zerror=par::defaultDouble, const double proxyerror=par::defaultDouble, const double Plambda_a=par::defaultDouble, const double Plambda_b=par::defaultDouble, const double Plambda_c=par::defaultDouble) 
      : Object(xx, yy, zz, ra, dec, redshift, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_logM(logM), m_mass_proxy(mass_proxy), m_mass_proxy_error(mass_proxy_error), m_bias(bias), m_concentration(conc), m_f_off(f_off), m_sigma_off(sigma_off), m_alpha_scaling_rel(alpha_scaling_rel), m_beta_scaling_rel(beta_scaling_rel), m_gamma_scaling_rel(gamma_scaling_rel), m_scatter0_scaling_rel(scatter0_scaling_rel), m_scatterM_scaling_rel(scatterM_scaling_rel), m_scatterM_exponent_scaling_rel(scatterM_exponent_scaling_rel), m_scatterz_scaling_rel(scatterz_scaling_rel), m_scatterz_exponent_scaling_rel(scatterz_exponent_scaling_rel), m_zbias(zbias), m_proxybias(proxybias), m_zerror(zerror), m_proxyerror(proxyerror), m_Plambda_a(Plambda_a), m_Plambda_b(Plambda_b), m_Plambda_c(Plambda_c) {}
      
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
       *  @brief get the private member \e m_logM
       *  @return the logarithm of the cluster mass
       */
      double logM () const override
      { return m_logM; }

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
      
      /**
       *  @brief get the private member \e m_concentration
       *  @return the concentration of the cluster
       */
      double concentration () const override
      { return m_concentration; }
      
      /**
       *  @brief get the private member \e m_f_off
       *  @return the miscentered cluster fraction
       */
      double f_off () const override
      { return m_f_off; }
      
      /**
       *  @brief get the private member \e m_sigma_off
       *  @return the rms of the miscentered cluster population
       */
      double sigma_off () const override
      { return m_sigma_off; }
      
      /**
       *  @brief get the private member \e m_alpha_scaling_rel
       *  @return normalization of the mass-observable scaling relation
       */
      double alpha_scaling_rel () const override
      { return m_alpha_scaling_rel; }
      
      /**
       *  @brief get the private member \e m_beta_scaling_rel
       *  @return slope of the mass-observable scaling relation
       */
      double  beta_scaling_rel () const override
      { return m_beta_scaling_rel; }
      
      /**
       *  @brief get the private member \e m_gamma_scaling_rel
       *  @return z evolution factor of the mass-observable scaling relation
       */
      double gamma_scaling_rel () const override
      { return m_gamma_scaling_rel; }

      /**
       *  @brief get the private member \e m_scatter0_scaling_rel
       *  @return scatter0
       */
      double scatter0_scaling_rel () const override
      { return m_scatter0_scaling_rel; }

      /**
       *  @brief get the private member \e m_scatterM_scaling_rel
       *  @return scatterM
       */
      double scatterM_scaling_rel () const override
      { return m_scatterM_scaling_rel; }

      /**
       *  @brief get the private member \e m_scatterM_exponent_scaling_rel
       *  @return scatterM_exponent
       */
      double scatterM_exponent_scaling_rel () const override
      { return m_scatterM_exponent_scaling_rel; }

      /**
       *  @brief get the private member \e m_scatterz_scaling_rel
       *  @return scatterz
       */
      double scatterz_scaling_rel () const override
      { return m_scatterz_scaling_rel; }

      /**
       *  @brief get the private member \e m_scatterz_exponent_scaling_rel
       *  @return scatterz_exponent
       */
      double scatterz_exponent_scaling_rel () const override
      { return m_scatterz_exponent_scaling_rel; }

      /**
       *  @brief get the private member \e m_zbias
       *  @return zbias
       */
      double zbias () const override
      { return m_zbias; }
      
      /**
       *  @brief get the private member \e m_proxybias
       *  @return proxybias
       */
      double proxybias () const override
      { return m_proxybias; }
      
      /**
       *  @brief get the private member \e m_zerror
       *  @return zerror
       */
      double zerror () const override
      { return m_zerror; }
      
      /**
       *  @brief get the private member \e m_proxyerror
       *  @return proxyerror
       */
      double proxyerror () const override
      { return m_proxyerror; }

      /**
       *  @brief get the private member \e m_Plambda_a
       *  @return Plambda_a
       */
      double Plambda_a () const override
      { return m_Plambda_a; }

      /**
       *  @brief get the private member \e m_Plambda_b
       *  @return Plambda_b
       */
      double Plambda_b () const override
      { return m_Plambda_b; }

      /**
       *  @brief get the private member \e m_Plambda_c
       *  @return Plambda_c
       */
      double Plambda_c () const override
      { return m_Plambda_c; }
      

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
       *  @brief set the private member \e m_logM
       *  @param logM the log-mass of the cluster
       */
      void set_logM (const double logM=par::defaultDouble) override
      { m_logM = logM; }

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
      
      /**
       *  @brief set the private member \e m_concentration
       *  @param conc concentration
       */
      void set_concentration (const double conc=par::defaultDouble) override
      { m_concentration = conc; }
      
      /**
       *  @brief set the private member \e m_f_off
       *  @param f_off f_off
       */
      void set_f_off (const double f_off=par::defaultDouble) override
      { m_f_off = f_off; }
      
      /**
       *  @brief set the private member \e m_sigma_off
       *  @param sigma_off sigma_off
       */
      void set_sigma_off (const double sigma_off=par::defaultDouble) override
      { m_sigma_off = sigma_off; }
      
      /**
       *  @brief set the private member \e m_alpha_scaling_rel
       *  @param alpha normalization of the cluster mass-observable scaling relation
       */
      void set_alpha_scaling_rel (const double alpha=par::defaultDouble) override
      { m_alpha_scaling_rel = alpha; }
      
      /**
       *  @brief set the private member \e m_beta_scaling_rel
       *  @param beta slope of the cluster mass-observable scaling relation
       */
      void set_beta_scaling_rel (const double beta=par::defaultDouble) override
      { m_beta_scaling_rel = beta; }
      
      /**
       *  @brief set the private member \e m_gamma_scaling_rel
       *  @param gamma z evolution factor of the cluster mass-observable scaling relation
       */
      void set_gamma_scaling_rel (const double gamma=par::defaultDouble) override
      { m_gamma_scaling_rel = gamma; }

      /**
       *  @brief set the private member \e m_scatter0_scaling_rel
       *  @param scatter0 scatter0
       */
      void set_scatter0_scaling_rel (const double scatter0=par::defaultDouble) override
      { m_scatter0_scaling_rel = scatter0; }

      /**
       *  @brief set the private member \e m_scatterM_scaling_rel
       *  @param scatterM scatterM
       */
      void set_scatterM_scaling_rel (const double scatterM=par::defaultDouble) override
      { m_scatterM_scaling_rel = scatterM; }

      /**
       *  @brief set the private member \e m_scatterM_exponent_scaling_rel
       *  @param scatterM_exponent scatterM_exponent
       */
      void set_scatterM_exponent_scaling_rel (const double scatterM_exponent=par::defaultDouble) override
      { m_scatterM_exponent_scaling_rel = scatterM_exponent; }

      /**
       *  @brief set the private member \e m_scatterz_scaling_rel
       *  @param scatterz scatterz
       */
      void set_scatterz_scaling_rel (const double scatterz=par::defaultDouble) override
      { m_scatterz_scaling_rel = scatterz; }

      /**
       *  @brief set the private member \e m_scatterz_exponent_scaling_rel
       *  @param scatterz_exponent scatterz_exponent
       */
      void set_scatterz_exponent_scaling_rel (const double scatterz_exponent=par::defaultDouble) override
      { m_scatterz_exponent_scaling_rel = scatterz_exponent; }

      /**
       *  @brief set the private member \e m_zbias
       *  @param zbias zbias
       */
      void set_zbias (const double zbias=par::defaultDouble) override
      { m_zbias = zbias; }
      
      /**
       *  @brief set the private member \e m_proxybias
       *  @param proxybias proxybias
       */
      void set_proxybias (const double proxybias=par::defaultDouble) override
      { m_proxybias = proxybias; }
      
      /**
       *  @brief set the private member \e m_zerror
       *  @param zerror zerror
       */
      void set_zerror (const double zerror=par::defaultDouble) override
      { m_zerror = zerror; }
      
      /**
       *  @brief set the private member \e m_proxyerror
       *  @param proxyerror proxyerror
       */
      void set_proxyerror (const double proxyerror=par::defaultDouble) override
      { m_proxyerror = proxyerror; }

      /**
       *  @brief set the private member \e m_Plambda_a
       *  @param Plambda_a Plambda_a
       */
      void set_Plambda_a (const double Plambda_a=par::defaultDouble) override
      { m_Plambda_a = Plambda_a; }

      /**
       *  @brief set the private member \e m_Plambda_b
       *  @param Plambda_b Plambda_b
       */
      void set_Plambda_b (const double Plambda_b=par::defaultDouble) override
      { m_Plambda_b = Plambda_b; }

      /**
       *  @brief set the private member \e m_Plambda_c
       *  @param Plambda_c Plambda_c
       */
      void set_Plambda_c (const double Plambda_c=par::defaultDouble) override
      { m_Plambda_c = Plambda_c; }

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
       *  @brief check if the private member \e m_logM is set
       *
       *  @return true if set; false otherwise
       */
      bool isSet_logM () override
      { return (cbl::isSet(m_logM)) ? true : false; }
      
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
      
      /**
       *  @brief check if the private member \e m_concentration is set
       *
       *  @return true if set; false otherwise
       */
      bool isSet_concentration () override
      { return (cbl::isSet(m_concentration)) ? true : false; }
      
      /**
       *  @brief check if the private member \e m_f_off is set
       *
       *  @return true if set; false otherwise
       */
      bool isSet_f_off () override
      { return (cbl::isSet(m_f_off)) ? true : false; }
      
      /**
       *  @brief check if the private member \e m_sigma_off is set
       *
       *  @return true if set; false otherwise
       */
      bool isSet_sigma_off () override
      { return (cbl::isSet(m_sigma_off)) ? true : false; }
      
      /**
       *  @brief check if the private member \e m_alpha_scaling_rel is set
       *
       *  @return true if set; false otherwise
       */
      bool isSet_alpha_scaling_rel () override
      { return (cbl::isSet(m_alpha_scaling_rel)) ? true : false; }
      
      /**
       *  @brief check if the private member \e m_beta_scaling_rel is set
       *
       *  @return true if set; false otherwise
       */
      bool isSet_beta_scaling_rel () override
      { return (cbl::isSet(m_beta_scaling_rel)) ? true : false; }
      
      /**
       *  @brief check if the private member \e m_gamma_scaling_rel is set
       *
       *  @return true if set; false otherwise
       */
      bool isSet_gamma_scaling_rel () override
      { return (cbl::isSet(m_gamma_scaling_rel)) ? true : false; }

      ///@}
      
    
    };
  }
}

#endif

