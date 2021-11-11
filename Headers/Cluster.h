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

#include "Cosmology.h"
#include "Object.h"


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
    
      /// pointer to the input cosmology
      std::shared_ptr<cosmology::Cosmology> m_cosmology = NULL;
      
      /// cluster density profile author
      std::string m_profile_author = par::defaultString;
      
      /// cluster density profile function
      double (Cluster::*m_profile) (const double, const double, std::vector<double>, const std::vector<cbl::cosmology::CosmologicalParameter>);
    
      /// if true, the concentration-mass relation is set
      bool m_isSet_cM_relation = false;
      
      /// if true, the density profile is set
      bool m_isSet_profile = false;
      
      /// pointer to the concentration-mass relation function
      double (Cluster::*m_concentration_from_mass) (std::vector<double>, const std::vector<cbl::cosmology::CosmologicalParameter>);
      
      /// parameter A in the c-M relation
      double m_A = par::defaultDouble;
      
      /// parameter B in the c-M relation
      double m_B = par::defaultDouble;
      
      /// parameter C in the c-M relation
      double m_C = par::defaultDouble;

      /// cluster mass
      double m_mass = par::defaultDouble;
      
      /// cluster mass logarithm (undefined base here)
      double m_logM = par::defaultDouble; 
    
      /// cluster mass proxy
      double m_mass_proxy = par::defaultDouble;

      /// cluster proxy error
      double m_mass_proxy_error = par::defaultDouble;
      
      /// cluster linear bias
      double m_bias = par::defaultDouble;
      
      /// cluster concentration
      double m_concentration = par::defaultDouble;
      
      /// fraction of miscentered cluster population
      double m_f_off = par::defaultDouble;
      
      /// rms of the miscentered distribution
      double m_sigma_off = par::defaultDouble;
      
      /// normalization of the mass-observable scaling relation
      double m_alpha_scaling_rel = par::defaultDouble;
      
      /// slope of the mass-observable scaling relation
      double m_beta_scaling_rel = par::defaultDouble;
      
      /// z evolution factor of the mass-observable scaling relation
      double m_gamma_scaling_rel = par::defaultDouble;
      
      /// constant term of the intrinsic scatter of the mass-observable scaling relation
      double m_scatter0_scaling_rel = par::defaultDouble;
      
      /// multiplicative factor in the mass/mass proxy dependent term in the intrinsic scatter of the mass-observable scaling relation
      double m_scatterM_scaling_rel = par::defaultDouble;
      
      /// exponent in the mass/mass proxy dependent term in the intrinsic scatter of the mass-observable scaling relation
      double m_scatterM_exponent_scaling_rel = par::defaultDouble;
      
      /// multiplicative factor in the redshift dependent term in the intrinsic scatter of the mass-observable scaling relation
      double m_scatterz_scaling_rel = par::defaultDouble;
      
      /// exponent in the redshift dependent term in the intrinsic scatter of the mass-observable scaling relation
      double m_scatterz_exponent_scaling_rel = par::defaultDouble;
      
      /// cluster redshift bias
      double m_zbias = par::defaultDouble;
      
      /// cluster mass proxy bias
      double m_proxybias = par::defaultDouble;
      
      /// cluster redshift uncertainty
      double m_zerror = par::defaultDouble;
      
      /// cluster mass proxy uncertainty
      double m_proxyerror = par::defaultDouble;
      
      /// \f$ a \f$ term in the function describing the cluster abundance, i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
      double m_Plambda_a = par::defaultDouble;
      
      /// \f$ b \f$ term in the function describing the cluster abundance, i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
      double m_Plambda_b = par::defaultDouble;
      
      /// \f$ c \f$ term in the function describing the cluster abundance, i.e. \f$ P(\lambda|z) = a \, \lambda^{-b} \, e^{-c\lambda} \f$, where \f$\lambda\f$ is a mass proxy
      double m_Plambda_c = par::defaultDouble;
      
      /**
       *  @name Protected member functions 
       */
      ///@{
	
      /**
       *  @brief The concentration-mass relation by Duffy et
       *  al. (2008): \f[c(M_h, z) = A(M_h/M_{pivot})^B\,(1+z)^C\f]
       *
       *  @param cosmo_par cosmological parameters
       *
       *  @param cosmo_par_names names of the cosmological parameters
       *
       *  @return the concentration
       *
       */
      double m_concentration_Duffy (const std::vector<double> cosmo_par={}, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_par_names={});
      
      /**
       *  @brief The NFW density profile
       *
       *  @param conc the concentration
       *
       *  @param rad the radius
       *
       *  @param cosmo_par cosmological parameters
       *
       *  @param cosmo_par_names names of the cosmological parameters
       *
       *  @return the NFW density profile
       *
       */
      double m_NFW (const double conc, const double rad, const std::vector<double> cosmo_par={}, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_par_names={});

      ///@}
      
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
	: Object() {}

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
       */
      Cluster (const comovingCoordinates coord, const double weight, const long region, const int ID, const std::string field, const double x_displacement, const double y_displacement, const double z_displacement) 
      : Object(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement) {}

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
       */
      Cluster (const comovingCoordinates coord, const cosmology::Cosmology &cosm, const double z1_guess, const double z2_guess, const double weight, const long region, const int ID, const std::string field, const double x_displacement, const double y_displacement, const double z_displacement) 
      : Object(coord, cosm, z1_guess, z2_guess, weight, region, ID, field, x_displacement, y_displacement, z_displacement) {}

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
       *  
       */
      Cluster (const observedCoordinates coord, const double weight, const long region, const int ID, const std::string field, const double x_displacement, const double y_displacement, const double z_displacement) 
      : Object(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement) {}
      
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
       */
      Cluster (const observedCoordinates coord, const CoordinateUnits inputUnits, const double weight, const long region, const int ID, const std::string field, const double x_displacement, const double y_displacement, const double z_displacement) 
      : Object(coord, inputUnits, weight, region, ID, field, x_displacement, y_displacement, z_displacement) {}
      
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
       */
      Cluster (const observedCoordinates coord, const cosmology::Cosmology &cosm, const double weight, const long region, const int ID, const std::string field, const double x_displacement, const double y_displacement, const double z_displacement) 
      : Object(coord, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement) {}

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
       */
      Cluster (const observedCoordinates coord, const CoordinateUnits inputUnits, const cosmology::Cosmology &cosm, const double weight, const long region, const int ID, const std::string field, const double x_displacement, const double y_displacement, const double z_displacement) 
      : Object(coord, inputUnits, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement) {}

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
       */
      Cluster (const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight, const long region, const int ID, const std::string field, const double x_displacement, const double y_displacement, const double z_displacement) 
      : Object(xx, yy, zz, ra, dec, redshift, weight, region, ID, field, x_displacement, y_displacement, z_displacement) {}
      
      /**
       *  @brief default destructor
       */
      ~Cluster () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to get the private members or
       *  to compute the cluster halo profiles
       */
      ///@{
      
      /**
       *  @brief The concentration obtained
       *  from a concentration-mass relation. If cosmo_par is provided,
       *  it is computed by updating the selected cosmological parameters
       *
       *  @param cosmo_par input cosmological parameters
       *
       *  @param cosmo_par_names names of the 
       *  input cosmological parameters
       *
       *  @return the concentration
       *
       */
      double concentration_from_mass(const std::vector<double> cosmo_par={}, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_par_names={});
	
      /**
       *  @brief compute the halo concentration
       *  @author Carlo Giocoli
       *  @author cgiocoli@gmail.com
       *  @param Vmax V<SUB>max</SUB>
       *  @param Rmax R<SUB>max</SUB>
       *  @return the halo concentration
       */
      double concentration2 (const double Vmax, const double Rmax) const;
      
      /**
       *  @brief the normalised halo density profile, computed by
       *  assuming a concentration-mass relation. 
       *  If cosmo_par is provided, it is computed
       *  by updating the selected cosmological parameters
       *
       *  this function computes the normalised density distribution
       *  of dark matter haloes; the Navarro-Frenk-White profile is
       *  the only one currently implemented (see e.g. eq. 74 of
       *  Cooray & Sheth 2002, eq. 68 of van den Bosch et
       *  al. 2012eq. 38 of Coe 2010):
       *
       *  \f[u_h(r, M_h, z) = \frac{\rho_h(r, M_h, z)}{M_h} =
       *  \frac{\rho_s}{(r/r_s)(1+r/r_s)^2}\f]
       *
       *  where
       *
       *  \f[\rho_s =
       *  \frac{\rho_{crit}\Delta_c}{3}\frac{c^3}{\ln(1+c)-c/(1+c)}\f]
       *
       *  the relation between the halo concentration,
       *  \f$c=c_{vir}=r_{vir}/r_s\f$, and halo mass, \f$M_h\f$, is
       *  computed by cbl::modelling::twopt::concentration;
       *  \f$\Delta_c(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::Delta_c and
       *  \f$\rho_{crit}(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::rho_crit; \f$r_{vir}(M_h,
       *  z)\f$ is computed by cbl::cosmology::Cosmology::r_vir
       *
       *  @param rad the scale
       *
       *  @param cosmo_par input cosmological parameters
       *
       *  @param cosmo_par_names names of the 
       *  input cosmological parameters
       *
       *  @return the halo density profile
       */
      double density_profile (const double rad, std::vector<double> cosmo_par={}, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_par_names={});
      
      /**
       *  @brief the normalised halo density profile. 
       *  If cosmo_par is provided, it is computed
       *  by updating the selected cosmological parameters
       *
       *  this function computes the normalised density distribution
       *  of dark matter haloes; the Navarro-Frenk-White profile is
       *  the only one currently implemented (see e.g. eq. 74 of
       *  Cooray & Sheth 2002, eq. 68 of van den Bosch et
       *  al. 2012eq. 38 of Coe 2010):
       *
       *  \f[u_h(r, M_h, z) = \frac{\rho_h(r, M_h, z)}{M_h} =
       *  \frac{\rho_s}{(r/r_s)(1+r/r_s)^2}\f]
       *
       *  where
       *
       *  \f[\rho_s =
       *  \frac{\rho_{crit}\Delta_c}{3}\frac{c^3}{\ln(1+c)-c/(1+c)}\f]
       *
       *  the relation between the halo concentration,
       *  \f$c=c_{vir}=r_{vir}/r_s\f$, and halo mass, \f$M_h\f$, is
       *  computed by cbl::modelling::twopt::concentration;
       *  \f$\Delta_c(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::Delta_c and
       *  \f$\rho_{crit}(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::rho_crit; \f$r_{vir}(M_h,
       *  z)\f$ is computed by cbl::cosmology::Cosmology::r_vir
       *
       *  @param rad the scale
       *
       *  @param conc the concentration
       *
       *  @param cosmo_par input cosmological parameters
       *
       *  @param cosmo_par_names names of the 
       *  input cosmological parameters
       *
       *  @return the halo density profile
       */
      double density_profile (const double rad, const double conc, std::vector<double> cosmo_par={}, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_par_names={});
      
      /**
       *  @brief the Fourier transform of the normalised halo density
       *  profile. If cosmo_par is provided, it is computed by
       *  updating the selected cosmological parameters
       *
       *  this function computes the Fourier transform of the
       *  normalised density distribution of dark matter haloes; the
       *  Navarro-Frenk-White profile is the only one currently
       *  implemented (see e.g. eq. 81 of Cooray & Sheth 2002 and 70
       *  of van den Bosch et al. 2012)
       *
       *  \f[\tilde{u}_h(k, M_h, z) = \frac{4\pi\rho_sr_s^3}{M_h}
       *  \left[\cos\mu\left[{\rm Ci}(\mu+\mu c) - {\rm
       *  Ci}(\mu)\right] + \sin\mu\left[{\rm Si}(\mu+\mu c) - {\rm
       *  Si}(\mu)\right] - \frac{\sin\mu c}{\mu+\mu c}\right]\f]
       *
       *  where 
       *
       *  \f[\mu\equiv kr_s\,,\f]
       *
       *  \f[\rho_s = \frac{\rho_{crit}\Delta_c}{3}
       *  \frac{c^3}{\ln(1+c)-c/(1+c)}\,,\f]
       *
       *  \f[{\rm Ci}(x)=-\int_x^\infty\frac{\cos t}{t}\,{\rm d}t\,,\f]
       *
       *  \f[{\rm Si}(x)=-\int_x^\infty\frac{\sin t}{t}\,{\rm d}t\f]
       *
       *  the relation between the halo concentration,
       *  \f$c=c_{vir}=r_{vir}/r_s\f$, and halo mass, \f$M_h\f$, is
       *  computed by cbl::modelling::twopt::concentration;
       *  \f$\Delta_c(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::Delta_c and
       *  \f$\rho_{crit}(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::rho_crit; \f$r_{vir}(M_h,
       *  z)\f$ is computed by cbl::cosmology::Cosmology::r_vir
       *
       *  @param kk the wave vector module at which the model is
       *  computed
       *
       *  @param cosmo_par input cosmological parameters
       *
       *  @param cosmo_par_names names of the 
       *  input cosmological parameters
       *
       *  @return the halo density profile
       */
      double density_profile_FourierSpace (const double kk, std::vector<double> cosmo_par={}, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_par_names={});
      
      /**
       *  @brief The halo mass converted to a different value of
       *  \f$\Delta\f$, assuming the Navarro-Frenk-White density
       *  profile
       *
       *  This function converts a given input mass \f$M_\Delta\f$ to
       *  \f$M_{\Delta^{new}}\f$ (e.g. \f$M_{500} \rightarrow
       *  M_{200}\f$).
       *
       *  Specifically, the algorithm currently implemented can be
       *  derived as follows. Given the Navarro-Frenk-White profile:
       *
       *  \f[ \rho(r) = \frac{\rho_0}{\frac{r}{R_s}
       *  \left(1+\frac{r}{R_s} \right)^2} \f]
       *
       *  the total halo mass contained within a radius \f$R_\Delta\f$
       *  is:
       *  
       *  \f[ M_\Delta = \int_0^{R_\Delta} 4\pi r^2\rho(r)dr = 4\pi
       *  \rho_0 R_s^3 \left[
       *  \ln(1+c_\Delta)-\frac{c_\Delta}{1+c_\Delta} \right] \f]
       *
       *  where the concentration is defined as \f$c_\Delta\equiv
       *  R_\Delta/R_s\f$. Thus, we can write:
       *
       *  \f[ \frac{\ln(1+c_\Delta) - \frac{c_\Delta}{1+c_\Delta}}
       *  {M_\Delta} = \frac{\ln(1+c_{\Delta^{new}}) -
       *  \frac{c_{\Delta^{new}}}{1+c_{\Delta^{new}}}}
       *  {M_{\Delta^{new}}} \f]
       *
       *  \f$ M_\Delta \f$ can be written as follows:
       *
       *  \f[ M_\Delta = \frac{4}{3}\pi\Delta\rho_{crit}R_\Delta^3 =
       *  \frac{\Delta}{\Delta^{new}}
       *  \left(\frac{R_\Delta}{R_{\Delta^{new}}} \right)^3
       *  M_{\Delta^{new}} = \frac{\Delta}{\Delta^{new}} x^3
       *  M_{\Delta^{new}} \f]
       *
       *  where \f$ x\equiv R_\Delta/R_{\Delta^{new}} \f$. Thus we
       *  have:
       *
       *  \f[ \frac{M_\Delta}{M_{\Delta^{new}}} \left[
       *  \ln(1+c_{\Delta^{new}}) -
       *  \frac{c_{\Delta^{new}}}{1+c_{\Delta^{new}}} \right] - \left[
       *  \ln(1+c_\Delta) - \frac{c_\Delta}{1+c_\Delta} \right] = 0
       *  \f]
       *
       *  where \f$ c_{\Delta^{new}} \equiv R_{\Delta^{new}}/R_s =
       *  c_\Delta R_{\Delta^{new}}/R_\Delta = c_\Delta/x \f$. The
       *  algorithm solves the above equation as a function of \f$ x
       *  \f$, providing in output:
       *
       *  \f[ M_{\Delta^{new}} =
       *  \frac{1}{x^3}\frac{\Delta^{new}}{\Delta} M_\Delta \f]
       *
       *  @param Mass the input mass \f$M_\Delta\f$
       *  (e.g. \f$M_{500}\f$)
       *
       *  @param Delta_in the input \f$\Delta\f$
       *  (e.g. \f$\Delta=500\f$)
       *
       *  @param Delta_out the output \f$\Delta^{new}\f$
       *  (e.g. \f$\Delta^{new}=200\f$)
       *
       *  @param conc the concentration related to either the input
       *  \f$\Delta\f$, or the output \f$\Delta^{new}\f$ (this is
       *  specified by is_input_conc)
       *
       *  @param is_input_conc true \f$\rightarrow\f$ the given
       *  concentration is related to the input \f$\Delta\f$; false
       *  \f$\rightarrow\f$ the given concentration is related to the
       *  output \f$\Delta^{new}\f$
       * 
       *  @param rRmin_guess the minimum guess value of
       *  \f$r_\Delta/r_{\Delta^{new}}\f$ used by the gsl minimisation
       *  function
       *
       *  @param rRmax_guess the maximum guess value of
       *  \f$r_\Delta/r_{\Delta^{new}}\f$ used by the gsl minimisation
       *  function
       *
       *  @return \f$M_{\Delta^{new}}\f$
       *
       *  @warning the current implementation assumes the
       *  Navarro-Frenk-White profile
       */
      double Mass_Delta (const double Mass, const double Delta_in, const double Delta_out, const double conc, const bool is_input_conc, const double rRmin_guess=1.e-3, const double rRmax_guess=10.) const;
    
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
       *  @brief set the concentration-mass relation. 
       *  A function is set, computing the concentration of a dark matter
       *  halo of a given a mass, at a given redshift; the models
       *  implemented are the following:
       *
       *  - Duffy et al. 2008:
       *  \f[c(M_h, z) = A(M_h/M_{pivot})^B\,(1+z)^C\f]
       *
       *
       *  @param cosmology the cosmology
       *
       *  @param Mass the halo mass
       *
       *  @param redshift the redshift
       *
       *  @param cM_author author(s) who proposed the 
       *  concentration-mass relation. Possibilities are:
       *  "Duffy" (Duffy et al. 2008)
       *
       *  @param profile_author the density profile author(s); available options are:
       *  "NFW" \f$\rightarrow\f$ Navarro-Frenk-White profile;
       *  "Einasto" \f$\rightarrow\f$ Einasto profile
       *
       *  @param halo_def the halo definition; available options are:
       *  "vir" \f$\rightarrow\f$ all matter within the radius
       *  \f$r_{vir}\f$ for which the mean internal density is
       *  \f$\Delta\f$ times the critical density
       *  \f$\rho_{crit}=3H^2/8\pi G\f$; "200" \f$\rightarrow\f$ all
       *  matter within the radius \f$r_{200}\f$ for which the mean
       *  internal density is 200 times the critical density; "mean"
       *  \f$\rightarrow\f$ all matter withing the radius
       *  \f$r_{200}\f$ for which the mean internal density is 200
       *  times the critical mean background density
       *
       *  @warning the Duffy et al. concentrantion-mass relation
       *  refers to the 0<z<2 redshift range, obtained from their full
       *  samples (see Table 1 of Duffy et al. 2008).
       *
       */
      void set_profile (const cbl::cosmology::Cosmology cosmology, const double Mass, const double redshift, const std::string cM_author="Duffy", const std::string profile_author="NFW", const std::string halo_def="vir");
      
      /**
       *  @brief set the cosmological model
       *
       *  @param cosmology the cosmology
       *
       */
      void set_cosmology (const cbl::cosmology::Cosmology cosmology);
      
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

