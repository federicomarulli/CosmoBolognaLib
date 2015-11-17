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
 *  @file Headers/Objects/Object.h
 *
 *  @brief The class Object 
 *
 *  This file defines the interface of the class Object, with only
 *  virtual methods that are implemented in the derived classes
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __OBJECT__
#define __OBJECT__ 


// ============================================================================================


namespace cosmobl {

  /**
   *  @class Object Object.h "Headers/Lib/Object.h"
   *
   *  @brief The class Object
   *
   *  This class is used to handle objects of type <EM> object
   *  </EM>
   */
  class Object {

  public :

    /**
     *  @brief default constructor
     *  @return object of class Object
     */
    Object () {}

    /**
     *  @brief default destructor
     *  @return none
     */
    virtual ~Object () {}

    /**
     *  @brief get the member \e m_xx
     *  @return the coordinate x of the derived object, or an error
     *  message if the derived object does not have this member
     */
    virtual double xx () const { cosmobl::ErrorMsg ("Error in xx() of Objech.h!"); return 0; }

    /**
     *  @brief get the member \e m_yy
     *  @return the coordinate y of the derived object, or an error
     *  message if the derived object does not have this member
     */
    virtual double yy () const { cosmobl::ErrorMsg ("Error in yy() of Objech.h!"); return 0; }

    /**
     *  @brief get the member \e m_zz
     *  @return the coordinate z of the derived object, or an error
     *  message if the derived object does not have this member
     */
    virtual double zz () const { cosmobl::ErrorMsg ("Error in zz() of Objech.h!"); return 0; }

    /**
     *  @brief get the member \e m_vx
     *  @return the peculiar velocity along the x direction of the
     *  derived object, or an error message if the derived object does
     *  not have this member
     */
    virtual double vx () const { cosmobl::ErrorMsg ("Error in vx() of Objech.h!"); return 0; }
    
    /**
     *  @brief get the member \e m_vy
     *  @return the peculiar velocity along the y direction of the
     *  derived object, or an error message if the derived object does
     *  not have this member
     */
    virtual double vy () const { cosmobl::ErrorMsg ("Error in vy() of Objech.h!"); return 0; }
    
    /**
     *  @brief get the member \e m_vz
     *  @return the peculiar velocity along the z direction of the
     *  derived object, or an message if the derived object does not
     *  have this member
     */
    virtual double vz () const { cosmobl::ErrorMsg ("Error in vz() of Objech.h!"); return 0; }
    
    /**
     *  @brief get the member \e m_dc
     *  @return the comoving distance of the derived object, or an
     *  error message if the derived object does not have this member
     */
    virtual double dc () const { cosmobl::ErrorMsg ("Error in dc() of Objech.h!"); return 0; }
    
    /**
     *  @brief get the member \e m_ra
     *  @return the Right Ascension of the derived object, or an error
     *  message if the derived object does not have this member
     */
    virtual double ra () const { cosmobl::ErrorMsg ("Error in ra() of Objech.h!"); return 0; }
    
    /**
     *  @brief get the member \e m_dec
     *  @return the Declination of the derived object, or an error
     *  message if the derived object does not have this member
     */
    virtual double dec () const { cosmobl::ErrorMsg ("Error in dec() of Objech.h!"); return 0; }
    
    /**
     *  @brief get the member \e m_redshift
     *  @return the redshift of the derived object, or an error
     *  message if the derived object does not have this member
     */
    virtual double redshift () const { cosmobl::ErrorMsg ("Error in redshift() of Objech.h!"); return 0; }
    
    /**
     *  @brief get the member \e m_weight
     *  @return the weight of the derived object, or an error message
     *  if the derived object does not have this member
     */
    virtual double weight () const { cosmobl::ErrorMsg ("Error in weight() of Objech.h!"); return 0; }

    /**
     *  @brief get the member \e m_region
     *  @return the index of the subRegion in which the object is located
     */
    virtual long region () const { cosmobl::ErrorMsg ("Error in region() of Objech.h!"); return 0; }

    /**
     *  @brief get the member \e m_mass
     *  @return the mass of the derived object, or an error message if
     *  the derived object does not have this member
     */
    virtual double mass () const { cosmobl::ErrorMsg ("Error in mass() of Objech.h!"); return 0; }
    
    /**
     *  @brief get the member \e m_magnitude
     *  @return the magnitude of the derived object, or an error
     *  message if the derived object does not have this member
     */
    virtual double magnitude () const { cosmobl::ErrorMsg ("Error in magnitude() of Objech.h!"); return 0; }  
    
    /**
     *  @brief get the member \e m_richness
     *  @return the richness of the derived object, or an error
     *  message if the derived object does not have this member
     */
    virtual double richness () const { cosmobl::ErrorMsg ("Error in richness() of Objech.h!"); return 0; }  
    
    /**
     *  @brief get the member \e m_generic
     *  @return the generic variable of the derived object, or an
     *  error message if the derived object does not have this member
     */
    virtual double generic () const { cosmobl::ErrorMsg ("Error in generic() of Objech.h!"); return 0; }  
    
    /**
     *  @brief get the member \e m_radius
     *  @return the radius of the derived object, or an
     *  error message if the derived object does not have this member
     */
    virtual double radius () const { cosmobl::ErrorMsg ("Error in radius() of Objech.h!"); return 0; }

    /**
     *  @brief get the object coordinates
     *  @return a vector containing the object coordinates
     */
    virtual vector<double> coords () const { cosmobl::ErrorMsg ("Error in generic() of Objech.h!"); vector<double> cc; return cc;}

    /**
     *  @brief set the member \e m_xx
     *  @param xx the coordinate x 
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_xx (const double xx) { cosmobl::ErrorMsg ("Error in set_xx() of Objech.h!"); }

    /**
     *  @brief set the member \e m_yy
     *  @param yy the coordinate y
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_yy (const double yy) { cosmobl::ErrorMsg ("Error in set_yy() of Objech.h!"); }

    /**
     *  @brief set the member \e m_zz
     *  @param zz the coordinate z 
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_zz (const double zz) { cosmobl::ErrorMsg ("Error in set_zz() of Objech.h!"); }

    /**
     *  @brief set the member \e m_vx
     *  @param vx the peculiar velocity along the x direction
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_vx (const double vx) { cosmobl::ErrorMsg ("Error in set_vx() of Objech.h!"); }
  
    /**
     *  @brief set the member \e m_vy
     *  @param vy the peculiar velocity along the y direction 
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_vy (const double vy) { cosmobl::ErrorMsg ("Error in set_vy() of Objech.h!"); }
    
    /**
     *  @brief set the member \e m_vz
     *  @param vz the peculiar velocity along the z direction
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_vz (const double vz) { cosmobl::ErrorMsg ("Error in set_vz() of Objech.h!"); }
    
    /**
     *  @brief set the member \e m_dc
     *  @param dc the comoving distance
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_dc (const double dc) { cosmobl::ErrorMsg ("Error in set_dc() of Objech.h!"); }
    
    /**
     *  @brief set the member \e m_ra
     *  @param ra the Right Ascension
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_ra (const double ra) { cosmobl::ErrorMsg ("Error in set_ra() of Objech.h!"); }
    
    /**
     *  @brief set the member \e m_dec
     *  @param dec the Declination
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_dec (const double dec) { cosmobl::ErrorMsg ("Error in set_dec() of Objech.h!"); }
    
    /**
     *  @brief set the member \e m_redshift
     *  @param redshift the redshift
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_redshift (const double redshift) { cosmobl::ErrorMsg ("Error in set_redshift() of Objech.h!"); }
    
    /**
     *  @brief set the member \e m_weight
     *  @param weight the weight
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_weight (const double weight) { cosmobl::ErrorMsg ("Error in set_weight() of Objech.h!"); }

    /**
     *  @brief set the member \e m_region
     *  @param region the index of the subRegion in which the object is located
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_region (const long region) { cosmobl::ErrorMsg ("Error in set_region() of Objech.h!"); }

    /**
     *  @brief set the member \e m_mass
     *  @param mass the mass
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_mass (const double mass) { cosmobl::ErrorMsg ("Error in set_mass() of Objech.h!"); }
    
    /**
     *  @brief set the member \e m_magnitude
     *  @param magnitude the magnitude
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_magnitude (const double magnitude) { cosmobl::ErrorMsg ("Error in set_magnitude() of Objech.h!"); }  
    
    /**
     *  @brief set the member \e m_richness
     *  @param richness the richness 
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_richness (const double richness) { cosmobl::ErrorMsg ("Error in set_richness() of Objech.h!"); }  
    
    /**
     *  @brief set the member \e m_generic
     *  @param generic the generic variable
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_generic (const double generic) { cosmobl::ErrorMsg ("Error in set_generic() of Objech.h!"); }  
    
    /**
     *  @brief set the member \e m_radius
     *  @param radius the radius
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_radius (const double radius) { cosmobl::ErrorMsg ("Error in set_radius() of Objech.h!"); }

  };
}

#endif
