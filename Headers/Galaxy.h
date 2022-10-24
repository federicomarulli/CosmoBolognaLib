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
 *  @file Headers/Galaxy.h
 *
 *  @brief The class Galaxy 
 *
 *  This file defines the interface of the class Galaxy, used to
 *  handle objects of type <EM> galaxies </EM>
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */


#ifndef __GALAXY__
#define __GALAXY__ 


// ===================================================================================================


namespace cbl {

  namespace catalogue {
    
    /**
     *  @class Galaxy Galaxy.h "Headers/Galaxy.h"
     *
     *  @brief The class Galaxy
     *
     *  This class is used to handle objects of type <EM> galaxy
     *  </EM>
     */
    class Galaxy : public Object {
    
    private :

      /// mass
      double m_mass = par::defaultDouble;      

      /// magnitude
      double m_magnitude = par::defaultDouble; 
      
      /// u-band magnitude
      double m_magnitudeU = par::defaultDouble; 
      
      /// g-band magnitude
      double m_magnitudeG = par::defaultDouble; 
      
      /// r-band magnitude
      double m_magnitudeR = par::defaultDouble; 
      
      /// i-band magnitude
      double m_magnitudeI = par::defaultDouble; 

      /// star formation rate
      double m_SFR = par::defaultDouble; 

      /// specific star formation rate
      double m_sSFR = par::defaultDouble; 
      
      /// first shear component
      double m_shear1 = par::defaultDouble;
      
      /// second shear component
      double m_shear2 = par::defaultDouble;
      
      /// odds
      double m_odds = par::defaultDouble;
      
      /// lensing weight
      double m_lensingWeight = par::defaultDouble;
      
      /// lensing calibration factor
      double m_lensingCalib = par::defaultDouble;

      /// Id_halo
      int m_IDHost = par::defaultInt;
     
      /// Tag of a galaxy "central" or "satellite"
      double m_galaxyTag = par::defaultDouble;

      /// stellar mass
      double m_mstar = par::defaultDouble;

      /// infall mass of substructure
      double m_massinfall = par::defaultDouble;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{
      
      /**
       *  @brief default constructor
       *  
       */
      Galaxy ()
	: Object() {}

      /**
       * @brief function that allows copying private variables of the class 
       * when an object of class Catalogue is copied
       * 
       * @return a shared pointer to the Object
       *
       */	
      std::shared_ptr<Object> getShared() {
        return std::make_shared<Galaxy>(*this);
      }

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
       *  @param mass the galaxy mass
       *
       *  @param magnitude the galaxy magnitude
       *
       *  @param magnitudeU the u-band galaxy magnitude
       *
       *  @param magnitudeG the g-band galaxy magnitude
       *
       *  @param magnitudeR the r-band galaxy magnitude
       *
       *  @param magnitudeI the i-band galaxy magnitude
       *
       *  @param SFR the galaxy star formation rate
       *
       *  @param sSFR the galaxy specific star formation rate
       *
       *  @param odds the lensing odds
       *
       *  @param shear1 first shear component
       *
       *  @param shear2 second shear component
       *
       *  @param lensingWeight lensing weight
       *
       *  @param lensingCalib lensing calibration factor
       *
       *  @param IDHost the Id of halo that host galaxy
       *
       *  @param galaxyTag identify the tag of a galaxy "central" or "satellite"
       *
       *  @param massinfall the infall mass of the substructure
       *
       *  @param mstar the galaxy stellar mass
       *  
       */
      Galaxy (const comovingCoordinates coord, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double magnitude=par::defaultDouble, const double magnitudeU=par::defaultDouble, const double magnitudeG=par::defaultDouble, const double magnitudeR=par::defaultDouble, const double magnitudeI=par::defaultDouble, const double SFR=par::defaultDouble, const double sSFR=par::defaultDouble, const double odds=par::defaultDouble, const double shear1=par::defaultDouble, const double shear2=par::defaultDouble, const double lensingWeight=par::defaultDouble, const double lensingCalib=par::defaultDouble, const int IDHost=par::defaultInt, const double galaxyTag= par::defaultDouble, const double mstar=par::defaultDouble, const double massinfall=par::defaultDouble) 
	: Object(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_magnitude(magnitude), m_magnitudeU(magnitudeU), m_magnitudeG(magnitudeG), m_magnitudeR(magnitudeR), m_magnitudeI(magnitudeI), m_SFR(SFR), m_sSFR(sSFR), m_shear1(shear1), m_shear2(shear2), m_odds(odds), m_lensingWeight(lensingWeight), m_lensingCalib(lensingCalib), m_IDHost(IDHost), m_galaxyTag(galaxyTag), m_mstar(mstar), m_massinfall(massinfall) {}

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
       *  @param mass the galaxy mass
       *
       *  @param magnitude the galaxy magnitude
       *
       *  @param magnitudeU the u-band galaxy magnitude
       *
       *  @param magnitudeG the g-band galaxy magnitude
       *
       *  @param magnitudeR the r-band galaxy magnitude
       *
       *  @param magnitudeI the i-band galaxy magnitude
       *
       *  @param SFR the galaxy star formation rate
       *
       *  @param sSFR the galaxy specific star formation rate
       *
       *  @param odds the lensing odds
       *
       *  @param shear1 first shear component
       *
       *  @param shear2 second shear component
       *
       *  @param lensingWeight lensing weight
       *
       *  @param lensingCalib lensing calibration factor
       *
       *  @param IDHost the Id of halo that host galaxy
       *
       *  @param galaxyTag identify the tag of a galaxy "central" or "satellite"
       *
       *  @param massinfall the infall mass of the substructure
       *
       *  @param mstar the galaxy stellar mass
       *  
       */
      Galaxy (const comovingCoordinates coord, const cosmology::Cosmology &cosm, const double z1_guess=0., const double z2_guess=10., const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double magnitude=par::defaultDouble, const double magnitudeU=par::defaultDouble, const double magnitudeG=par::defaultDouble, const double magnitudeR=par::defaultDouble, const double magnitudeI=par::defaultDouble, const double SFR=par::defaultDouble, const double sSFR=par::defaultDouble, const double odds=par::defaultDouble, const double shear1=par::defaultDouble, const double shear2=par::defaultDouble, const double lensingWeight=par::defaultDouble, const double lensingCalib=par::defaultDouble, const int IDHost=par::defaultInt, const double galaxyTag= par::defaultDouble, const double mstar=par::defaultDouble, const double massinfall=par::defaultDouble) 
	: Object(coord, cosm, z1_guess, z2_guess, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_magnitude(magnitude), m_magnitudeU(magnitudeU), m_magnitudeG(magnitudeG), m_magnitudeR(magnitudeR), m_magnitudeI(magnitudeI), m_SFR(SFR), m_sSFR(sSFR), m_shear1(shear1), m_shear2(shear2), m_odds(odds), m_lensingWeight(lensingWeight), m_lensingCalib(lensingCalib), m_IDHost(IDHost), m_galaxyTag(galaxyTag), m_mstar(mstar), m_massinfall(massinfall) {}

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
       *  @param mass the galaxy mass
       *
       *  @param magnitude the galaxy magnitude
       *
       *  @param magnitudeU the u-band galaxy magnitude
       *
       *  @param magnitudeG the g-band galaxy magnitude
       *
       *  @param magnitudeR the r-band galaxy magnitude
       *
       *  @param magnitudeI the i-band galaxy magnitude
       *
       *  @param SFR the galaxy star formation rate
       *
       *  @param sSFR the galaxy specific star formation rate
       *
       *  @param odds the lensing odds
       *
       *  @param shear1 first shear component
       *
       *  @param shear2 second shear component
       *
       *  @param lensingWeight lensing weight
       *
       *  @param lensingCalib lensing calibration factor
       *
       *  @param IDHost the Id of halo that host galaxy
       *
       *  @param galaxyTag identify the tag of a galaxy "central" or "satellite"
       *
       *  @param massinfall the infall mass of the substructure
       *
       *  @param mstar the galaxy stellar mass
       *  
       */
      Galaxy (const observedCoordinates coord, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double magnitude=par::defaultDouble, const double magnitudeU=par::defaultDouble, const double magnitudeG=par::defaultDouble, const double magnitudeR=par::defaultDouble, const double magnitudeI=par::defaultDouble, const double SFR=par::defaultDouble, const double sSFR=par::defaultDouble, const double odds=par::defaultDouble, const double shear1=par::defaultDouble, const double shear2=par::defaultDouble, const double lensingWeight=par::defaultDouble, const double lensingCalib=par::defaultDouble, const int IDHost=par::defaultInt, const double galaxyTag= par::defaultDouble, const double mstar=par::defaultDouble, const double massinfall=par::defaultDouble) 
	: Object(coord, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_magnitude(magnitude), m_magnitudeU(magnitudeU), m_magnitudeG(magnitudeG), m_magnitudeR(magnitudeR), m_magnitudeI(magnitudeI), m_SFR(SFR), m_sSFR(sSFR), m_shear1(shear1), m_shear2(shear2), m_odds(odds), m_lensingWeight(lensingWeight), m_lensingCalib(lensingCalib), m_IDHost(IDHost), m_galaxyTag(galaxyTag), m_mstar(mstar), m_massinfall(massinfall) {}
      
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
       *  @param mass the galaxy mass
       *
       *  @param magnitude the galaxy magnitude
       *
       *  @param magnitudeU the u-band galaxy magnitude
       *
       *  @param magnitudeG the g-band galaxy magnitude
       *
       *  @param magnitudeR the r-band galaxy magnitude
       *
       *  @param magnitudeI the i-band galaxy magnitude
       *
       *  @param SFR the galaxy star formation rate
       *
       *  @param sSFR the galaxy specific star formation rate
       *
       *  @param odds the lensing odds
       *
       *  @param shear1 first shear component
       *
       *  @param shear2 second shear component
       *
       *  @param lensingWeight lensing weight
       *
       *  @param lensingCalib lensing calibration factor
       *
       *  @param IDHost the Id of halo that host galaxy
       *
       *  @param galaxyTag identify the tag of a galaxy "central" or "satellite"
       *
       *  @param massinfall the infall mass of the substructure
       *
       *  @param mstar the galaxy stellar mass
       *
       */
      Galaxy (const observedCoordinates coord, const CoordinateUnits inputUnits, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double magnitude=par::defaultDouble, const double magnitudeU=par::defaultDouble, const double magnitudeG=par::defaultDouble, const double magnitudeR=par::defaultDouble, const double magnitudeI=par::defaultDouble, const double SFR=par::defaultDouble, const double sSFR=par::defaultDouble, const double odds=par::defaultDouble, const double shear1=par::defaultDouble, const double shear2=par::defaultDouble, const double lensingWeight=par::defaultDouble, const double lensingCalib=par::defaultDouble, const int IDHost=par::defaultInt, const double galaxyTag= par::defaultDouble, const double mstar=par::defaultDouble, const double massinfall=par::defaultDouble) 
	: Object(coord, inputUnits, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_magnitude(magnitude), m_magnitudeU(magnitudeU), m_magnitudeG(magnitudeG), m_magnitudeR(magnitudeR), m_magnitudeI(magnitudeI), m_SFR(SFR), m_sSFR(sSFR), m_shear1(shear1), m_shear2(shear2), m_odds(odds), m_lensingWeight(lensingWeight), m_lensingCalib(lensingCalib), m_IDHost(IDHost), m_galaxyTag(galaxyTag), m_mstar(mstar), m_massinfall(massinfall) {}
      
      /**
       *  @brief constructor that uses observed coordinates in radians
       *  and a cosmological model to estimate the comoving
       *  coordinates
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshift}
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
       *  @param mass the galaxy mass
       *
       *  @param magnitude the galaxy magnitude
       *
       *  @param magnitudeU the u-band galaxy magnitude
       *
       *  @param magnitudeG the g-band galaxy magnitude
       *
       *  @param magnitudeR the r-band galaxy magnitude
       *
       *  @param magnitudeI the i-band galaxy magnitude
       *
       *  @param SFR the galaxy star formation rate
       *
       *  @param sSFR the galaxy specific star formation rate
       *
       *  @param odds the lensing odds
       *
       *  @param shear1 first shear component
       *
       *  @param shear2 second shear component
       *
       *  @param lensingWeight lensing weight
       *
       *  @param lensingCalib lensing calibration factor
       *
       *  @param IDHost the Id of halo that host galaxy
       *
       *  @param galaxyTag identify the tag of a galaxy "central" or "satellite"
       *
       *  @param massinfall the infall mass of the substructure
       *
       *  @param mstar the galaxy stellar mass
       *
       */
      Galaxy (const observedCoordinates coord, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double magnitude=par::defaultDouble, const double magnitudeU=par::defaultDouble, const double magnitudeG=par::defaultDouble, const double magnitudeR=par::defaultDouble, const double magnitudeI=par::defaultDouble, const double SFR=par::defaultDouble, const double sSFR=par::defaultDouble, const double odds=par::defaultDouble, const double shear1=par::defaultDouble, const double shear2=par::defaultDouble, const double lensingWeight=par::defaultDouble, const double lensingCalib=par::defaultDouble, const int IDHost=par::defaultInt, const double galaxyTag= par::defaultDouble, const double mstar=par::defaultDouble, const double massinfall=par::defaultDouble) 
	: Object(coord, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_magnitude(magnitude), m_magnitudeU(magnitudeU), m_magnitudeG(magnitudeG), m_magnitudeR(magnitudeR), m_magnitudeI(magnitudeI), m_SFR(SFR), m_sSFR(sSFR), m_shear1(shear1), m_shear2(shear2), m_odds(odds), m_lensingWeight(lensingWeight), m_lensingCalib(lensingCalib), m_IDHost(IDHost), m_galaxyTag(galaxyTag), m_mstar(mstar), m_massinfall(massinfall) {}

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
       *  @param mass the galaxy mass
       *
       *  @param magnitude the galaxy magnitude
       *
       *  @param magnitudeU the u-band galaxy magnitude
       *
       *  @param magnitudeG the g-band galaxy magnitude
       *
       *  @param magnitudeR the r-band galaxy magnitude
       *
       *  @param magnitudeI the i-band galaxy magnitude
       *
       *  @param SFR the galaxy star formation rate
       *
       *  @param sSFR the galaxy specific star formation rate
       *
       *  @param odds the lensing odds
       *
       *  @param shear1 first shear component
       *
       *  @param shear2 second shear component
       *
       *  @param lensingWeight lensing weight
       *
       *  @param lensingCalib lensing calibration factor
       *  
       *  @param IDHost the Id of halo that host galaxy
       *
       *  @param galaxyTag identify the tag of a galaxy "central" or "satellite"
       *
       *  @param massinfall the infall mass of the substructure
       *
       *  @param mstar the galaxy stellar mass
       *
       */
      Galaxy (const observedCoordinates coord, const CoordinateUnits inputUnits, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double magnitude=par::defaultDouble, const double magnitudeU=par::defaultDouble, const double magnitudeG=par::defaultDouble, const double magnitudeR=par::defaultDouble, const double magnitudeI=par::defaultDouble, const double SFR=par::defaultDouble, const double sSFR=par::defaultDouble, const double odds=par::defaultDouble, const double shear1=par::defaultDouble, const double shear2=par::defaultDouble, const double lensingWeight=par::defaultDouble, const double lensingCalib=par::defaultDouble, const int IDHost=par::defaultInt, const double galaxyTag= par::defaultDouble, const double mstar=par::defaultDouble, const double massinfall=par::defaultDouble) 
	: Object(coord, inputUnits, cosm, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_magnitude(magnitude), m_magnitudeU(magnitudeU), m_magnitudeG(magnitudeG), m_magnitudeR(magnitudeR), m_magnitudeI(magnitudeI), m_SFR(SFR), m_sSFR(sSFR), m_shear1(shear1), m_shear2(shear2), m_odds(odds), m_lensingWeight(lensingWeight), m_lensingCalib(lensingCalib), m_IDHost(IDHost), m_galaxyTag(galaxyTag), m_mstar(mstar), m_massinfall(massinfall) {}

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
       *  @param region region, used e.g. for jackknife and boostrap 
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
       *  @param mass the galaxy mass
       *
       *  @param magnitude the galaxy magnitude
       *
       *  @param magnitudeU the u-band galaxy magnitude
       *
       *  @param magnitudeG the g-band galaxy magnitude
       *
       *  @param magnitudeR the r-band galaxy magnitude
       *
       *  @param magnitudeI the i-band galaxy magnitude
       *
       *  @param SFR the galaxy star formation rate
       *
       *  @param sSFR the galaxy specific star formation rate
       *
       *  @param odds the lensing odds
       *
       *  @param shear1 first shear component
       *
       *  @param shear2 second shear component
       *
       *  @param lensingWeight lensing weight
       *
       *  @param lensingCalib lensing calibration factor
       *  
       *  @param IDHost the Id of halo that host galaxy
       *
       *  @param galaxyTag identify the tag of a galaxy "central" or "satellite"
       *
       *  @param massinfall the infall mass of the substructure
       *
       *  @param mstar the galaxy stellar mass
       *
       */
      Galaxy (const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight=1., const long region=par::defaultLong, const int ID=par::defaultInt, const std::string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double redshiftMin=par::defaultDouble, const double redshiftMax=par::defaultDouble, const double sn=par::defaultDouble, const double mass=par::defaultDouble, const double magnitude=par::defaultDouble, const double magnitudeU=par::defaultDouble, const double magnitudeG=par::defaultDouble, const double magnitudeR=par::defaultDouble, const double magnitudeI=par::defaultDouble, const double SFR=par::defaultDouble, const double sSFR=par::defaultDouble, const double odds=par::defaultDouble, const double shear1=par::defaultDouble, const double shear2=par::defaultDouble, const double lensingWeight=par::defaultDouble, const double lensingCalib=par::defaultDouble, const int IDHost=par::defaultInt, const double galaxyTag= par::defaultDouble, const double mstar=par::defaultDouble, const double massinfall=par::defaultDouble) 
	: Object(xx, yy, zz, ra, dec, redshift, weight, region, ID, field, x_displacement, y_displacement, z_displacement, redshiftMin, redshiftMax, sn), m_mass(mass), m_magnitude(magnitude), m_magnitudeU(magnitudeU), m_magnitudeG(magnitudeG), m_magnitudeR(magnitudeR), m_magnitudeI(magnitudeI), m_SFR(SFR), m_sSFR(sSFR), m_shear1(shear1), m_shear2(shear2), m_odds(odds), m_lensingWeight(lensingWeight), m_lensingCalib(lensingCalib), m_IDHost(IDHost), m_galaxyTag(galaxyTag), m_mstar(mstar), m_massinfall(massinfall) {}
      
      /**
       *  @brief default destructor
       */
      ~Galaxy () = default;
    
      ///@}
  
  
      /**
       *  @name Member functions used to get the private members 
       */
      ///@{
    
      /**
       *  @brief get the private member \e m_mass
       *  @return the mass of the galaxy
       */
      double mass () const override
      { return m_mass; }

      /**
       *  @brief get the private member \e m_magnitude
       *  @return the magnitude of the galaxy
       */
      double magnitude () const override
      { return m_magnitude; }
      
      /**
       *  @brief get the private member \e m_magnitudeU
       *  @return the u-band magnitude of the galaxy
       */
      double magnitudeU () const override
      { return m_magnitudeU; }
      
      /**
       *  @brief get the private member \e m_magnitudeG
       *  @return the g-band magnitude of the galaxy
       */
      double magnitudeG () const override
      { return m_magnitudeG; }
      
      /**
       *  @brief get the private member \e m_magnitudeR
       *  @return the r-band magnitude of the galaxy
       */
      double magnitudeR () const override
      { return m_magnitudeR; }
      
      /**
       *  @brief get the private member \e m_magnitudeI
       *  @return the i-band magnitude of the galaxy
       */
      double magnitudeI () const override
      { return m_magnitudeI; }
      
      /**
       *  @brief get the private member \e m_odds
       *  @return the odds of the galaxy
       */
      double odds () const override
      { return m_odds; }
      
      /**
       *  @brief get the private member \e m_shear1
       *  @return the first shear component of the galaxy
       */
      double shear1 () const override
      { return m_shear1; }
      
      /**
       *  @brief get the private member \e m_shear2
       *  @return the second shear component of the galaxy
       */
      double shear2 () const override
      { return m_shear2; }
      
      /**
       *  @brief get the private member \e m_lensingWeight
       *  @return the lensing weight of the galaxy
       */
      double lensingWeight () const override
      { return m_lensingWeight; }
      
      /**
       *  @brief get the private member \e m_lensingCalib
       *  @return the lensing calibration factor of the galaxy
       */
      double lensingCalib () const override
      { return m_lensingCalib; }

      /**
       *  @brief get the private member \e m_SFR
       *  @return the star formation rate of the galaxy
       */
      double SFR () const override
      { return m_SFR; }

      /**
       *  @brief get the private member \e m_sSFR
       *  @return the specific star formation rate of the galaxy
       */
      double sSFR () const override
      { return m_sSFR; }

      /**
       *  @brief get the private member \e m_IDHost
       *  @return the Id of the halo that host galaxy
       */
      int IDHost () const override
      { return m_IDHost; }

      /**
       *  @brief get the private member \e m_galaxyTag
       *  @return the tag of the galaxy
       */
      double galaxyTag () const override
      { return m_galaxyTag; }

      /**
       *  @brief get the private member \e m_mstar
       *  @return the stellar mass of the galaxy
       */
      double mstar () const override
      { return m_mstar; }

      /**
       *  @brief get the private member \e m_massinfall
       *  @return the infall mass of the substructure
       */
      double massinfall () const override
      { return m_massinfall; }

      ///@}
  
  
      /**
       *  @name Member functions used to set the private members 
       */
      ///@{
    
      /**
       *  @brief set the private member \e m_mass
       *  @param mass the mass of the galaxy
       */
      void set_mass (const double mass=par::defaultDouble) override
      { m_mass = mass; } 

      /**
       *  @brief set the private member \e m_magnitude
       *  @param magnitude the magnitude of the galaxy
       */
      void set_magnitude (const double magnitude=par::defaultDouble) override
      { m_magnitude = magnitude; }
      
      /**
       *  @brief set the private member \e m_magnitudeU
       *  @param magnitudeU u-band magnitude the magnitude of the galaxy
       */
      void set_magnitudeU (const double magnitudeU=par::defaultDouble) override
      { m_magnitudeU = magnitudeU; }
      
      /**
       *  @brief set the private member \e m_magnitudeG
       *  @param magnitudeG g-band magnitude the magnitude of the galaxy
       */
      void set_magnitudeG (const double magnitudeG=par::defaultDouble) override
      { m_magnitudeG = magnitudeG; }
      
      /**
       *  @brief set the private member \e m_magnitudeR
       *  @param magnitudeR r-band magnitude the magnitude of the galaxy
       */
      void set_magnitudeR (const double magnitudeR=par::defaultDouble) override
      { m_magnitudeR = magnitudeR; }
      
      /**
       *  @brief set the private member \e m_magnitudeI
       *  @param magnitudeI i-band magnitude the magnitude of the galaxy
       */
      void set_magnitudeI (const double magnitudeI=par::defaultDouble) override
      { m_magnitudeI = magnitudeI; }
      
      /**
       *  @brief set the private member \e m_odds
       *  @param odds odds of the galaxy
       */
      void set_odds (const double odds=par::defaultDouble) override
      { m_odds = odds; }
      
      /**
       *  @brief set the private member \e m_shear1
       *  @param shear1 first shear component of the galaxy
       */
      void set_shear1 (const double shear1=par::defaultDouble) override
      { m_shear1 = shear1; }
      
      /**
       *  @brief set the private member \e m_shear2
       *  @param shear2 second shear component of the galaxy
       */
      void set_shear2 (const double shear2=par::defaultDouble) override
      { m_shear2 = shear2; }
      
      /**
       *  @brief set the private member \e m_lensingWeight
       *  @param lensingWeight lensing weight of the galaxy
       */
      void set_lensingWeight (const double lensingWeight=par::defaultDouble) override
      { m_lensingWeight = lensingWeight; }
      
      /**
       *  @brief set the private member \e m_lensingCalib
       *  @param lensingCalib lensing calibration factor of the galaxy
       */
      void set_lensingCalib (const double lensingCalib=par::defaultDouble) override
      { m_lensingCalib = lensingCalib; }

      /**
       *  @brief set the private member \e m_SFR
       *  @param SFR the star formation rate of the galaxy
       */
      void set_SFR (const double SFR=par::defaultDouble) override
      { m_SFR = SFR; }
    
      /**
       *  @brief set the private member \e m_sSFR
       *  @param sSFR the specific star formation rate of the galaxy
       */
      void set_sSFR (const double sSFR=par::defaultDouble) override
      { m_sSFR = sSFR; }

      /**
       *  @brief set the private member \e m_IDHost
       *  @param IDHost the Id of the halo that host galaxy
       */
      void set_IDHost (const int IDHost=par::defaultInt) override
      { m_IDHost = IDHost; }
      
      /**
       *  @brief set the private member \e m_galaxyTag
       *  @param galaxyTag the tag of a galaxy
       */
      void set_galaxyTag (const double galaxyTag) override
      { m_galaxyTag = galaxyTag; }

      /**
       *  @brief set the private member \e m_mstar
       *  @param mstar the stellar mass of the galaxy
       */
      void set_mstar (const double mstar=par::defaultDouble) override
      { m_mstar = mstar; }

      /**
       *  @brief set the private member \e m_massinfall
       *  @param massinfall the infall mass of the substructure
       */
      void set_massinfall (const double massinfall=par::defaultDouble) override
      { m_massinfall = massinfall; }
      
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
       *  @brief check if the private member \e m_magnitude is set
       *  
       *  @return true if the magnitude is set; false otherwise
       */
      bool isSet_magnitude () override
      { return (cbl::isSet(m_magnitude)) ? true : false; }
      
      /**
       *  @brief check if the private member \e m_magnitudeU is set
       *  
       *  @return true if the u magnitude is set; false otherwise
       */
      bool isSet_magnitudeU () override
      { return (cbl::isSet(m_magnitudeU)) ? true : false; }
      
      /**
       *  @brief check if the private member \e m_magnitudeG is set
       *  
       *  @return true if the g magnitude is set; false otherwise
       */
      bool isSet_magnitudeG () override
      { return (cbl::isSet(m_magnitudeG)) ? true : false; }
      
      /**
       *  @brief check if the private member \e m_magnitudeR is set
       *  
       *  @return true if the r magnitude is set; false otherwise
       */
      bool isSet_magnitudeR () override
      { return (cbl::isSet(m_magnitudeR)) ? true : false; }
      
      /**
       *  @brief check if the private member \e m_magnitudeI is set
       *  
       *  @return true if the i magnitude is set; false otherwise
       */
      bool isSet_magnitudeI () override
      { return (cbl::isSet(m_magnitudeI)) ? true : false; }

      /**
       *  @brief check if the private member \e m_odds is set
       *  
       *  @return true if the odds is set; false otherwise
       */
      bool isSet_odds () override
      { return (cbl::isSet(m_odds)) ? true : false; }
      
      /**
       *  @brief check if the private member \e m_shear1 is set
       *  
       *  @return true if the shear is set; false otherwise
       */
      bool isSet_shear1 () override
      { return (cbl::isSet(m_shear1)) ? true : false; }
      
      /**
       *  @brief check if the private member \e m_shear2 is set
       *  
       *  @return true if the shear is set; false otherwise
       */
      bool isSet_shear2 () override
      { return (cbl::isSet(m_shear2)) ? true : false; }
      
      /**
       *  @brief check if the private member \e m_lensingWeight is set
       *  
       *  @return true if the lensing weight is set; false otherwise
       */
      bool isSet_lensingWeight () override
      { return (cbl::isSet(m_lensingWeight)) ? true : false; }
      
      /**
       *  @brief check if the private member \e m_lensingCalib is set
       *  
       *  @return true if the lensing calibration factor is set; false otherwise
       */
      bool isSet_lensingCalib () override
      { return (cbl::isSet(m_lensingCalib)) ? true : false; }
      
      /**
       *  @brief check if the private member \e m_SFR is set
       *  
       *  @return true if the SFR is set; false otherwise
       */
      bool isSet_SFR () override
      { return (cbl::isSet(m_SFR)) ? true : false; }
    
      /**
       *  @brief check if the private member \e m_sSFR is set
       *  
       *  @return true if the sSFR is set; false otherwise
       */
      bool isSet_sSFR () override
      { return (cbl::isSet(m_sSFR)) ? true : false; }

      /**
       *  @brief check if the private member \e m_IDHost is set
       *  
       *  @return true if the IDHost is set; false otherwise
       */
      bool isSet_IDHost () override
      { return (cbl::isSet(m_IDHost)) ? true : false; }
    
      /**
       *  @brief check if the private member \e m_galaxyTag is set
       *  
       *  @return true if the galaxyTag is set; false otherwise
       */
      bool isSet_galaxyTag () override
      { return (cbl::isSet(m_galaxyTag)) ? true : false; }

      /**
       *  @brief check if the private member \e m_mstar is set
       *  
       *  @return true if the mstar is set; false otherwise
       */
      bool isSet_mstar () override
      { return (cbl::isSet(m_mstar)) ? true : false; }

      /**
       *  @brief check if the private member \e m_massinfall is set
       *  
       *  @return true if the massinfall is set; false otherwise
       */
      bool isSet_massinfall () override
      { return (cbl::isSet(m_massinfall)) ? true : false; }

      ///@}
      
    };
  }
}

#endif
