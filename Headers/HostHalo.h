/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Tommaso Ronconi      *
 *  federico.marulli3@sissa.it tronconi@sissa.it                    *
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
 *  @file Headers/HostHalo.h
 *
 *  @brief The class HostHalo 
 *
 *  This file defines the interface of the class HostHalo, used to handle
 *  objects of type <EM> dark matter particles group </EM>
 *
 *  @authors Federico Marulli, Tommaso Ronconi
 *
 *  @authors federico.marulli3@unbo.it, tronconi@sissa.it
 */

#ifndef __HOSTHALO__
#define __HOSTHALO__


// ===================================================================================================


namespace cbl {

  namespace catalogue {
    
    /**
     *  @class HostHalo HostHalo.h "Headers/Lib/HostHalo.h"
     *
     *  @brief The class HostHalo
     *
     *  This class is used to handle objects of type <EM> group
     *  </EM>
     */
    class HostHalo : public Halo {

    private :

      /// gas component of the total group mass
      double m_tot_mass;

      /// star formation rate
      double m_SFR;

      /// mass estimate
      double m_mass_estimate;

      /// radius estimate
      double m_radius_estimate;

      /// velocity dispersion estimate
      double m_veldisp_estimate;

      /// x-coordinate of the centre of mass
      double m_xcm;

      /// y-coordinate of the centre of mass
      double m_ycm;

      /// z-coordinate of the centre of mass
      double m_zcm;

      /// x-axis component of the spin
      double m_spin_x;

      /// y-axis component of the spin
      double m_spin_y;

      /// z-axis component of the spin
      double m_spin_z;

      /// velocity dispersion
      double m_veldisp;

      /// maximum velocity
      double m_vmax;

      /// maximum radial velocity
      double m_vmax_rad;

      /// half-mass radius
      double m_radius;

      /// unambiguous group ID
      int m_ID;

      /// ID of the parent halo
      int m_parent;

      /// number of sub-haloes within the main group
      int m_nsub;

      /// vector of pointers to satellites
      std::vector<std::shared_ptr<Object>> m_satellites;
      
    public :
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class HostHalo
       */
      HostHalo ()
	: Halo(), m_tot_mass(par::defaultDouble), m_mass_estimate(par::defaultDouble), m_radius_estimate(par::defaultDouble), m_veldisp_estimate(par::defaultDouble), m_xcm(par::defaultDouble), m_ycm(par::defaultDouble), m_zcm(par::defaultDouble), m_spin_x(par::defaultDouble), m_spin_y(par::defaultDouble), m_spin_z(par::defaultDouble), m_veldisp(par::defaultDouble), m_vmax(par::defaultDouble), m_vmax_rad(par::defaultDouble), m_radius(par::defaultDouble), m_ID(par::defaultInt), m_parent(par::defaultInt), m_nsub(par::defaultInt), m_satellites({}) {}
      
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
       *  @param vx velocity along the x-axis
       *
       *  @param vy velocity along the y-axis
       *
       *  @param vz velocity along the z-axis
       *
       *  @param mass halo mass
       *
       *  @param tot_mass gas component of the total group mass
       *
       *  @param mass_estimate mass estimate
       *
       *  @param radius_estimate radius estimate
       *
       *  @param veldisp_estimate velocity dispersion estimate
       *
       *  @param cm_coord x-coordinate of the centre of mass
       *
       *  @param spin_x x-axis component of the spin
       *
       *  @param spin_y y-axis component of the spin
       *
       *  @param spin_z z-axis component of the spin
       *
       *  @param veldisp velocity dispersion
       *
       *  @param vmax maximum velocity
       *
       *  @param vmax_rad maximum radial velocity
       *
       *  @param radius half-mass radius
       *
       *  @param ID unambiguous group ID
       *
       *  @param parent ID of the parent halo
       *
       *  @param nsub number of sub-haloes within the main group
       *
       *  @param satellites vector of pointers to satellites
       *
       *  @return object of class HostHalo
       */
      HostHalo (const comovingCoordinates coord, const double weight=1., const long region=par::defaultLong, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble, const double tot_mass=par::defaultDouble, const double mass_estimate=par::defaultDouble, const double radius_estimate=par::defaultDouble, const double veldisp_estimate=par::defaultDouble, const comovingCoordinates cm_coord={par::defaultDouble, par::defaultDouble, par::defaultDouble}, const double spin_x=par::defaultDouble, const double spin_y=par::defaultDouble, const double spin_z=par::defaultDouble, const double veldisp=par::defaultDouble, const double vmax=par::defaultDouble, const double vmax_rad=par::defaultDouble, const double radius=par::defaultDouble, const int ID=par::defaultInt, const int parent=par::defaultInt, const int nsub=par::defaultInt, std::vector<std::shared_ptr<Object>> satellites = {})
	: Halo(coord, weight, region, field, x_displacement, y_displacement, z_displacement, vx, vy, vz, mass), m_tot_mass(tot_mass), m_mass_estimate(mass_estimate), m_radius_estimate(radius_estimate), m_veldisp_estimate(veldisp_estimate), m_xcm(cm_coord.xx), m_ycm(cm_coord.yy), m_zcm(cm_coord.zz), m_spin_x(spin_x), m_spin_y(spin_y), m_spin_z(spin_z), m_veldisp(veldisp), m_vmax(vmax), m_vmax_rad(vmax_rad), m_radius(radius), m_ID(ID), m_parent(parent), m_nsub(nsub), m_satellites(satellites) {}

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
       *  @param vx velocity along the x-axis
       *
       *  @param vy velocity along the y-axis
       *
       *  @param vz velocity along the z-axis
       *
       *  @param mass halo mass
       *
       *  @param tot_mass gas component of the total group mass
       *
       *  @param mass_estimate mass estimate
       *
       *  @param radius_estimate radius estimate
       *
       *  @param veldisp_estimate velocity dispersion estimate
       *
       *  @param cm_coord x-coordinate of the centre of mass
       *
       *  @param spin_x x-axis component of the spin
       *
       *  @param spin_y y-axis component of the spin
       *
       *  @param spin_z z-axis component of the spin
       *
       *  @param veldisp velocity dispersion
       *
       *  @param vmax maximum velocity
       *
       *  @param vmax_rad maximum radial velocity
       *
       *  @param radius half-mass radius
       *
       *  @param ID unambiguous group ID
       *
       *  @param parent ID of the parent halo
       *
       *  @param nsub number of sub-haloes within the main group
       *
       *  @param satellites vector of pointers to satellites
       *
       *  @return object of class HostHalo
       */
      HostHalo (const comovingCoordinates coord, const cosmology::Cosmology &cosm, const double z1_guess=0., const double z2_guess=10., const double weight=1., const long region=par::defaultLong, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble, const double tot_mass=par::defaultDouble, const double mass_estimate=par::defaultDouble, const double radius_estimate=par::defaultDouble, const double veldisp_estimate=par::defaultDouble, const comovingCoordinates cm_coord={par::defaultDouble, par::defaultDouble, par::defaultDouble}, const double spin_x=par::defaultDouble, const double spin_y=par::defaultDouble, const double spin_z=par::defaultDouble, const double veldisp=par::defaultDouble, const double vmax=par::defaultDouble, const double vmax_rad=par::defaultDouble, const double radius=par::defaultDouble, const int ID=par::defaultInt, const int parent=par::defaultInt, const int nsub=par::defaultInt, std::vector<std::shared_ptr<Object>> satellites = {}) 
	: Halo(coord, cosm, z1_guess, z2_guess, weight, region, field, x_displacement, y_displacement, z_displacement, vx, vy, vz, mass), m_tot_mass(tot_mass), m_mass_estimate(mass_estimate), m_radius_estimate(radius_estimate), m_veldisp_estimate(veldisp_estimate), m_xcm(cm_coord.xx), m_ycm(cm_coord.yy), m_zcm(cm_coord.zz), m_spin_x(spin_x), m_spin_y(spin_y), m_spin_z(spin_z), m_veldisp(veldisp), m_vmax(vmax), m_vmax_rad(vmax_rad), m_radius(radius), m_ID(ID), m_parent(parent), m_nsub(nsub), m_satellites(satellites) {}

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
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the mock mass
       *
       *  @param tot_mass gas component of the total group mass
       *
       *  @param mass_estimate mass estimate
       *
       *  @param radius_estimate radius estimate
       *
       *  @param veldisp_estimate velocity dispersion estimate
       *
       *  @param cm_coord x-coordinate of the centre of mass
       *
       *  @param spin_x x-axis component of the spin
       *
       *  @param spin_y y-axis component of the spin
       *
       *  @param spin_z z-axis component of the spin
       *
       *  @param veldisp velocity dispersion
       *
       *  @param vmax maximum velocity
       *
       *  @param vmax_rad maximum radial velocity
       *
       *  @param radius half-mass radius
       *
       *  @param ID unambiguous group ID
       *
       *  @param parent ID of the parent halo
       *
       *  @param nsub number of sub-haloes within the main group
       *
       *  @param satellites vector of pointers to satellites
       *
       *  @return object of class HostHalo
       */
      HostHalo (const observedCoordinates coord, const double weight=1., const long region=par::defaultLong, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble, const double tot_mass=par::defaultDouble, const double mass_estimate=par::defaultDouble, const double radius_estimate=par::defaultDouble, const double veldisp_estimate=par::defaultDouble, const comovingCoordinates cm_coord={par::defaultDouble, par::defaultDouble, par::defaultDouble}, const double spin_x=par::defaultDouble, const double spin_y=par::defaultDouble, const double spin_z=par::defaultDouble, const double veldisp=par::defaultDouble, const double vmax=par::defaultDouble, const double vmax_rad=par::defaultDouble, const double radius=par::defaultDouble, const int ID=par::defaultInt, const int parent=par::defaultInt, const int nsub=par::defaultInt, std::vector<std::shared_ptr<Object>> satellites = {}) 
	: Halo(coord, weight, region, field, x_displacement, y_displacement, z_displacement, vx, vy, vz, mass), m_tot_mass(tot_mass), m_mass_estimate(mass_estimate), m_radius_estimate(radius_estimate), m_veldisp_estimate(veldisp_estimate), m_xcm(cm_coord.xx), m_ycm(cm_coord.yy), m_zcm(cm_coord.zz), m_spin_x(spin_x), m_spin_y(spin_y), m_spin_z(spin_z), m_veldisp(veldisp), m_vmax(vmax), m_vmax_rad(vmax_rad), m_radius(radius), m_ID(ID), m_parent(parent), m_nsub(nsub), m_satellites(satellites) {}
      
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
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the mock mass
       *
       *  @param tot_mass gas component of the total group mass
       *
       *  @param mass_estimate mass estimate
       *
       *  @param radius_estimate radius estimate
       *
       *  @param veldisp_estimate velocity dispersion estimate
       *
       *  @param cm_coord x-coordinate of the centre of mass
       *
       *  @param spin_x x-axis component of the spin
       *
       *  @param spin_y y-axis component of the spin
       *
       *  @param spin_z z-axis component of the spin
       *
       *  @param veldisp velocity dispersion
       *
       *  @param vmax maximum velocity
       *
       *  @param vmax_rad maximum radial velocity
       *
       *  @param radius half-mass radius
       *
       *  @param ID unambiguous group ID
       *
       *  @param parent ID of the parent halo
       *
       *  @param nsub number of sub-haloes within the main group
       *
       *  @param satellites vector of pointers to satellites
       *
       *  @return object of class HostHalo
       */
      HostHalo (const observedCoordinates coord, const CoordinateUnits inputUnits, const double weight=1., const long region=par::defaultLong, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble, const double tot_mass=par::defaultDouble, const double mass_estimate=par::defaultDouble, const double radius_estimate=par::defaultDouble, const double veldisp_estimate=par::defaultDouble, const comovingCoordinates cm_coord={par::defaultDouble, par::defaultDouble, par::defaultDouble}, const double spin_x=par::defaultDouble, const double spin_y=par::defaultDouble, const double spin_z=par::defaultDouble, const double veldisp=par::defaultDouble, const double vmax=par::defaultDouble, const double vmax_rad=par::defaultDouble, const double radius=par::defaultDouble, const int ID=par::defaultInt, const int parent=par::defaultInt, const int nsub=par::defaultInt, std::vector<std::shared_ptr<Object>> satellites = {}) 
	: Halo(coord, inputUnits, weight, region, field, x_displacement, y_displacement, z_displacement, vx, vy, vz, mass), m_tot_mass(tot_mass), m_mass_estimate(mass_estimate), m_radius_estimate(radius_estimate), m_veldisp_estimate(veldisp_estimate), m_xcm(cm_coord.xx), m_ycm(cm_coord.yy), m_zcm(cm_coord.zz), m_spin_x(spin_x), m_spin_y(spin_y), m_spin_z(spin_z), m_veldisp(veldisp), m_vmax(vmax), m_vmax_rad(vmax_rad), m_radius(radius), m_ID(ID), m_parent(parent), m_nsub(nsub), m_satellites(satellites) {}
      
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
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the mock mass
       *
       *  @param tot_mass gas component of the total group mass
       *
       *  @param mass_estimate mass estimate
       *
       *  @param radius_estimate radius estimate
       *
       *  @param veldisp_estimate velocity dispersion estimate
       *
       *  @param cm_coord x-coordinate of the centre of mass
       *
       *  @param spin_x x-axis component of the spin
       *
       *  @param spin_y y-axis component of the spin
       *
       *  @param spin_z z-axis component of the spin
       *
       *  @param veldisp velocity dispersion
       *
       *  @param vmax maximum velocity
       *
       *  @param vmax_rad maximum radial velocity
       *
       *  @param radius half-mass radius
       *
       *  @param ID unambiguous group ID
       *
       *  @param parent ID of the parent halo
       *
       *  @param nsub number of sub-haloes within the main group
       *
       *  @param satellites vector of pointers to satellites
       *
       *  @return object of class HostHalo
       */
      HostHalo (const observedCoordinates coord, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble, const double tot_mass=par::defaultDouble, const double mass_estimate=par::defaultDouble, const double radius_estimate=par::defaultDouble, const double veldisp_estimate=par::defaultDouble, const comovingCoordinates cm_coord={par::defaultDouble, par::defaultDouble, par::defaultDouble}, const double spin_x=par::defaultDouble, const double spin_y=par::defaultDouble, const double spin_z=par::defaultDouble, const double veldisp=par::defaultDouble, const double vmax=par::defaultDouble, const double vmax_rad=par::defaultDouble, const double radius=par::defaultDouble, const int ID=par::defaultInt, const int parent=par::defaultInt, const int nsub=par::defaultInt, std::vector<std::shared_ptr<Object>> satellites = {}) 
	: Halo(coord, cosm, weight, region, field, x_displacement, y_displacement, z_displacement, vx, vy, vz, mass), m_tot_mass(tot_mass), m_mass_estimate(mass_estimate), m_radius_estimate(radius_estimate), m_veldisp_estimate(veldisp_estimate), m_xcm(cm_coord.xx), m_ycm(cm_coord.yy), m_zcm(cm_coord.zz), m_spin_x(spin_x), m_spin_y(spin_y), m_spin_z(spin_z), m_veldisp(veldisp), m_vmax(vmax), m_vmax_rad(vmax_rad), m_radius(radius), m_ID(ID), m_parent(parent), m_nsub(nsub), m_satellites(satellites) {}

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
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the mock mass
       *
       *  @param tot_mass gas component of the total group mass
       *
       *  @param mass_estimate mass estimate
       *
       *  @param radius_estimate radius estimate
       *
       *  @param veldisp_estimate velocity dispersion estimate
       *
       *  @param cm_coord x-coordinate of the centre of mass
       *
       *  @param spin_x x-axis component of the spin
       *
       *  @param spin_y y-axis component of the spin
       *
       *  @param spin_z z-axis component of the spin
       *
       *  @param veldisp velocity dispersion
       *
       *  @param vmax maximum velocity
       *
       *  @param vmax_rad maximum radial velocity
       *
       *  @param radius half-mass radius
       *
       *  @param ID unambiguous group ID
       *
       *  @param parent ID of the parent halo
       *
       *  @param nsub number of sub-haloes within the main group
       *
       *  @param satellites vector of pointers to satellites
       *
       *  @return object of class HostHalo
       */
      HostHalo (const observedCoordinates coord, const CoordinateUnits inputUnits, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble, const double tot_mass=par::defaultDouble, const double mass_estimate=par::defaultDouble, const double radius_estimate=par::defaultDouble, const double veldisp_estimate=par::defaultDouble, const comovingCoordinates cm_coord={par::defaultDouble, par::defaultDouble, par::defaultDouble}, const double spin_x=par::defaultDouble, const double spin_y=par::defaultDouble, const double spin_z=par::defaultDouble, const double veldisp=par::defaultDouble, const double vmax=par::defaultDouble, const double vmax_rad=par::defaultDouble, const double radius=par::defaultDouble, const int ID=par::defaultInt, const int parent=par::defaultInt, const int nsub=par::defaultInt, std::vector<std::shared_ptr<Object>> satellites = {}) 
	: Halo(coord, inputUnits, cosm, weight, region, field, x_displacement, y_displacement, z_displacement, vx, vy, vz, mass), m_tot_mass(tot_mass), m_mass_estimate(mass_estimate), m_radius_estimate(radius_estimate), m_veldisp_estimate(veldisp_estimate), m_xcm(cm_coord.xx), m_ycm(cm_coord.yy), m_zcm(cm_coord.zz), m_spin_x(spin_x), m_spin_y(spin_y), m_spin_z(spin_z), m_veldisp(veldisp), m_vmax(vmax), m_vmax_rad(vmax_rad), m_radius(radius), m_ID(ID), m_parent(parent), m_nsub(nsub), m_satellites(satellites) {}

      /**
       *  @brief constructor that uses both comoving and observed coordinates
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
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the mock mass
       *
       *  @param tot_mass gas component of the total group mass
       *
       *  @param mass_estimate mass estimate
       *
       *  @param radius_estimate radius estimate
       *
       *  @param veldisp_estimate velocity dispersion estimate
       *
       *  @param cm_coord x-coordinate of the centre of mass
       *
       *  @param spin_x x-axis component of the spin
       *
       *  @param spin_y y-axis component of the spin
       *
       *  @param spin_z z-axis component of the spin
       *
       *  @param veldisp velocity dispersion
       *
       *  @param vmax maximum velocity
       *
       *  @param vmax_rad maximum radial velocity
       *
       *  @param radius half-mass radius
       *
       *  @param ID unambiguous group ID
       *
       *  @param parent ID of the parent halo
       *
       *  @param nsub number of sub-haloes within the main group
       *
       *  @param satellites vector of pointers to satellites
       *
       *  @return object of class HostHalo
       */
      HostHalo (const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight=1., const long region=par::defaultLong, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble, const double tot_mass=par::defaultDouble, const double mass_estimate=par::defaultDouble, const double radius_estimate=par::defaultDouble, const double veldisp_estimate=par::defaultDouble, const comovingCoordinates cm_coord={par::defaultDouble, par::defaultDouble, par::defaultDouble}, const double spin_x=par::defaultDouble, const double spin_y=par::defaultDouble, const double spin_z=par::defaultDouble, const double veldisp=par::defaultDouble, const double vmax=par::defaultDouble, const double vmax_rad=par::defaultDouble, const double radius=par::defaultDouble, const int ID=par::defaultInt, const int parent=par::defaultInt, const int nsub=par::defaultInt, std::vector<std::shared_ptr<Object>> satellites = {}) 
	: Halo(xx, yy, zz, ra, dec, redshift, weight, region, field, x_displacement, y_displacement, z_displacement, vx, vy, vz, mass), m_tot_mass(tot_mass), m_mass_estimate(mass_estimate), m_radius_estimate(radius_estimate), m_veldisp_estimate(veldisp_estimate), m_xcm(cm_coord.xx), m_ycm(cm_coord.yy), m_zcm(cm_coord.zz), m_spin_x(spin_x), m_spin_y(spin_y), m_spin_z(spin_z), m_veldisp(veldisp), m_vmax(vmax), m_vmax_rad(vmax_rad), m_radius(radius), m_ID(ID), m_parent(parent), m_nsub(nsub), m_satellites(satellites) {}
      
      /**
       *  @brief default destructor
       *  @return none
       */
      ~HostHalo () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to get the private members 
       */
      ///@{
    
      /**
       *  @brief get the private member HostHalo::m_tot_mass
       *  @return the total mass of the host halo (sum over all contributions)
       */
      double tot_mass () const override { return m_tot_mass; }
    
      /**
       *  @brief get the private member HostHalo::m_mass_estimate
       *  @return the mass estimate of the host halo
       */
      double mass_estimate () const override { return m_mass_estimate; }
    
      /**
       *  @brief get the private member HostHalo::m_radius_estimate
       *  @return the radius estimate of the host halo
       */
      double radius_estimate () const override { return m_radius_estimate; }
    
      /**
       *  @brief get the private member HostHalo::m_veldisp_estimate
       *  @return the velocity dispersion estimate of the host halo
       */
      double veldisp_estimate () const override { return m_veldisp_estimate; }
    
      /**
       *  @brief get the private member HostHalo::m_xcm
       *  @return the x-axis coordinate of the centre of mass
       */
      double xcm () const override { return m_xcm; }
    
      /**
       *  @brief get the private member HostHalo::m_ycm
       *  @return the y-axis coordinate of the centre of mass
       */
      double ycm () const override { return m_ycm; }
    
      /**
       *  @brief get the private member HostHalo::m_zcm
       *  @return the z-axis coordinate of the centre of mass
       */
      double zcm () const override { return m_zcm; }
    
      /**
       *  @brief get the private member HostHalo::m_spin_x
       *  @return the x-axis component of the spin
       */
      double spin_x () const override { return m_spin_x; }
    
      /**
       *  @brief get the private member HostHalo::m_spin_y
       *  @return the y-axis component of the spin
       */
      double spin_y () const override { return m_spin_y; }
    
      /**
       *  @brief get the private member HostHalo::m_spin_z
       *  @return the z-axis component of the spin
       */
      double spin_z () const override { return m_spin_z; }
    
      /**
       *  @brief get the private member HostHalo::m_veldisp
       *  @return the velocity dispersion of the host halo
       */
      double veldisp () const override { return m_veldisp; }
    
      /**
       *  @brief get the private member HostHalo::m_vmax
       *  @return the maximum total velocity of the host halo
       */
      double vmax () const override { return m_vmax; }
    
      /**
       *  @brief get the private member HostHalo::m_vmax_rad
       *  @return the maximum radial velocity of the host halo
       */
      double vmax_rad () const override { return m_vmax_rad; }
    
      /**
       *  @brief get the private member HostHalo::m_radius
       *  @return the radius of the host halo
       */
      double radius () const override { return m_radius; }
    
      /**
       *  @brief get the private member HostHalo::m_ID
       *  @return the unambiguous host halo ID
       */
      int ID () const override { return m_ID; }
    
      /**
       *  @brief get the private member HostHalo::m_parent
       *  @return the id of the parent host halo
       */
      int parent () const override { return m_parent; }
    
      /**
       *  @brief get the private member HostHalo::m_nsub
       *  @return the number of sub-groups in the host halo
       */
      int nsub () const override { return m_nsub; }
    
      /**
       *  @brief get the private member HostHalo::m_satellites
       *  @return the vector of pointers to the satellite objects
       */
      std::vector<std::shared_ptr<Object>> satellites () const override { return m_satellites; }

      ///@}
  
  
      /**
       *  @name Member functions used to set the private members 
       */
      ///@{
    
      /**
       *  @brief set the private member HostHalo::m_tot_mass
       *  @param tot_mass the total mass of the parent halo
       *  @return none
       */
      void set_tot_mass (const double tot_mass=par::defaultDouble) override { m_tot_mass = tot_mass; }
    
      /**
       *  @brief set the private member HostHalo::m_mass_estimate
       *  @param mass_estimate the mass estimate of the host halo
       *  @return none
       */
      void set_mass_estimate (const double mass_estimate=par::defaultDouble) override { m_mass_estimate = mass_estimate; }
    
      /**
       *  @brief set the private member HostHalo::m_radius_estimate
       *  @param radius_estimate the radius estimate of the host halo
       *  @return none
       */
      void set_radius_estimate (const double radius_estimate=par::defaultDouble) override { m_radius_estimate = radius_estimate; }
    
    
      /**
       *  @brief set the private member HostHalo::m_veldisp_estimate
       *  @param veldisp_estimate the velocity dispersion estimate of the host halo
       *  @return none
       */
      void set_veldisp_estimate (const double veldisp_estimate=par::defaultDouble) override { m_veldisp_estimate = veldisp_estimate; }
    
      /**
       *  @brief set the private member HostHalo::m_xcm
       *  @param xcm the x-axis coordinate of the centre of mass
       *  @return none
       */
      void set_xcm (const double xcm=par::defaultDouble) override { m_xcm = xcm; }
    
      /**
       *  @brief set the private member HostHalo::m_ycm
       *  @param ycm the y-axis coordinate of the centre of mass
       *  @return none
       */
      void set_ycm (const double ycm=par::defaultDouble) override { m_ycm = ycm; }
    
      /**
       *  @brief set the private member HostHalo::m_zcm
       *  @param zcm the z-axis coordinate of the centre of mass
       *  @return none
       */
      void set_zcm (const double zcm=par::defaultDouble) override { m_zcm = zcm; }
    
      /**
       *  @brief set the private member HostHalo::m_spin_x
       *  @param spin_x the x-axis component of the spin
       *  @return none
       */
      void set_spin_x (const double spin_x=par::defaultDouble) override { m_spin_x = spin_x; }
    
      /**
       *  @brief set the private member HostHalo::m_spin_y
       *  @param spin_y the y-axis component of the spin
       *  @return none
       */
      void set_spin_y (const double spin_y=par::defaultDouble) override { m_spin_y = spin_y; }
    
      /**
       *  @brief set the private member HostHalo::m_spin_z
       *  @param spin_z the z-axis component of the spin
       *  @return none
       */
      void set_spin_z (const double spin_z=par::defaultDouble) override { m_spin_z = spin_z; }
    
      /**
       *  @brief set the private member HostHalo::m_veldisp
       *  @param veldisp the velocity dispersion of the host halo
       *  @return none
       */
      void set_veldisp (const double veldisp=par::defaultDouble) override { m_veldisp = veldisp; }
    
      /**
       *  @brief set the private member HostHalo::m_vmax
       *  @param vmax the maximum total velocity of the host halo
       *  @return none
       */
      void set_vmax (const double vmax=par::defaultDouble) override { m_vmax = vmax; }
    
      /**
       *  @brief set the private member HostHalo::m_vmax_rad
       *  @param vmax_rad the maximum radial velocity of the host halo
       *  @return none
       */
      void set_vmax_rad (const double vmax_rad=par::defaultDouble) override { m_vmax_rad = vmax_rad; }
    
      /**
       *  @brief set the private member HostHalo::m_radius
       *  @param radius the maximum radial velocity of the host halo
       *  @return none
       */
      void set_radius (const double radius=par::defaultDouble) override { m_radius = radius; }
      
      /**
       *  @brief set the private member HostHalo::m_ID
       *  @param ID the unambiguous host halo ID
       *  @return none
       */
      void set_ID (const int ID=par::defaultInt) override { m_ID = ID; }
    
      /**
       *  @brief set the private member HostHalo::m_parent
       *  @param parent the id of the parent host halo
       *  @return none
       */
      void set_parent (const int parent=par::defaultInt) override { m_parent = parent; }

      /**
       *  @brief set the private member HostHalo::m_nsub
       *  @param nsub the number of sub-groups in the host halo
       *  @return none
       */
      void set_nsub (const int nsub=par::defaultInt) override { m_nsub = nsub; }

      /**
       *  @brief set the private member HostHalo::m_satellites
       *  @param satellite the vector of shared pointers to satellite
       *  objects
       *  @return none
       */
      void set_satellite (const std::shared_ptr<Object> satellite={}) override {
	m_satellites.push_back(satellite);
      }

      /**
       *  @brief set the private member HostHalo::m_satellites
       *  @param satellites the vector of shared pointers to satellite objects
       *  @return none
       */
      void set_satellites (const std::vector<std::shared_ptr<Object>> satellites={}) override {
	for (size_t ii = 0; ii<satellites.size(); ii++)	m_satellites.push_back(satellites[ii]);
      }
    
      ///@}


      /**
       *  @name Member functions used to check if the private members are set 
       */
      ///@{
    	
      /**
       *  @brief get the private member \e m_tot_mass
       *  
       *  @return true if the total mass is set; false otherwise
       */
      bool isSet_tot_mass () override
      { return (cbl::isSet(m_tot_mass)) ? true : false; }
    	
      /**
       *  @brief get the private member \e m_mass_estimate
       *  
       *  @return true if the mass estimate is set; false otherwise
       */
      bool isSet_mass_estimate () override
      { return (cbl::isSet(m_mass_estimate)) ? true : false; }
    	
      /**
       *  @brief get the private member \e m_radius_estimate
       *  
       *  @return true if the radius estimate is set; false otherwise
       */
      bool isSet_radius_estimate () override
      { return (cbl::isSet(m_radius_estimate)) ? true : false; }
    	
      /**
       *  @brief get the private member \e m_veldisp_estimate
       *  
       *  @return true if the velocity dispersion estimate is set; false otherwise
       */
      bool isSet_veldisp_estimate () override
      { return (cbl::isSet(m_veldisp_estimate)) ? true : false; }
    	
      /**
       *  @brief get the private member \e m_xcm
       *  
       *  @return true if the x-axis coordinate of the centre of mass is set; false otherwise
       */
      bool isSet_xcm () override
      { return (cbl::isSet(m_xcm)) ? true : false; }
    	
      /**
       *  @brief get the private member \e m_ycm
       *  
       *  @return true if the y-axis coordinate of the centre of mass is set; false otherwise
       */
      bool isSet_ycm () override
      { return (cbl::isSet(m_ycm)) ? true : false; }
    	
      /**
       *  @brief get the private member \e m_zcm
       *  
       *  @return true if the z-axis coordinate of the centre of mass is set; false otherwise
       */
      bool isSet_zcm () override
      { return (cbl::isSet(m_zcm)) ? true : false; }
    	
      /**
       *  @brief get the private member \e m_spin_x
       *  
       *  @return true if the x-axis component of the spin is set; false otherwise
       */
      bool isSet_spin_x () override
      { return (cbl::isSet(m_spin_x)) ? true : false; }
    	
      /**
       *  @brief get the private member \e m_spin_y
       *  
       *  @return true if the y-axis component of the spin is set; false otherwise
       */
      bool isSet_spin_y () override
      { return (cbl::isSet(m_spin_y)) ? true : false; }
    	
      /**
       *  @brief get the private member \e m_spin_z
       *  
       *  @return true if the z-axis component of the spin is set; false otherwise
       */
      bool isSet_spin_z () override
      { return (cbl::isSet(m_spin_z)) ? true : false; }
    	
      /**
       *  @brief get the private member \e m_veldisp
       *  
       *  @return true if the velocity dispersion is set; false otherwise
       */
      bool isSet_veldisp () override
      { return (cbl::isSet(m_veldisp)) ? true : false; }
    	
      /**
       *  @brief get the private member \e m_vmax
       *  
       *  @return true if the maximum total velocity is set; false otherwise
       */
      bool isSet_vmax () override
      { return (cbl::isSet(m_vmax)) ? true : false; }
    	
      /**
       *  @brief get the private member \e m_vmax_rad
       *  
       *  @return true if the maximum radial velocity is set; false otherwise
       */
      bool isSet_vmax_rad () override
      { return (cbl::isSet(m_vmax_rad)) ? true : false; }
    	
      /**
       *  @brief get the private member \e m_radius
       *  
       *  @return true if the radius is set; false otherwise
       */
      bool isSet_radius () override
      { return (cbl::isSet(m_radius)) ? true : false; }
    	
      /**
       *  @brief get the private member \e m_ID
       *  
       *  @return true if the unambiguous host halo ID is set; false otherwise
       */
      bool isSet_ID () override
      { return (cbl::isSet(m_ID)) ? true : false; }
    	
      /**
       *  @brief get the private member \e m_parent
       *  
       *  @return true if the id of the parent is set; false otherwise
       */
      bool isSet_parent () override
      { return (cbl::isSet(m_parent)) ? true : false; }
    	
      /**
       *  @brief get the private member \e m_nsub
       *  
       *  @return true if the the number of sub-groups in the host halo is set; false otherwise
       */
      bool isSet_nsub () override
      { return (cbl::isSet(m_nsub)) ? true : false; }
    
      ///@}
  
    }; //class HostHalo

  } //namespace catalogue
  
} //namespace cbl

#endif //__HOSTHALO__
