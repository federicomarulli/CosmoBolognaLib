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
 *  @file Headers/Objects/Mock.h
 *
 *  @brief The class Mock 
 *
 *  This file defines the interface of the class Mock, used to handle
 *  objects of type <EM> mock object </EM>
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MOCK__
#define __MOCK__ 


// ===================================================================================================


namespace cosmobl {

  namespace catalogue {
    
    /**
     *  @class Mock Mock.h "Headers/Lib/Mock.h"
     *
     *  @brief The class Mock
     *
     *  This class is used to handle objects of type <EM> mock object
     *  </EM>
     */
    class Mock : public Halo {

    private :

      /// generic variable of the mock object
      double m_generic; 

    
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Mock
       */
      Mock () {}
    
      /**
       *  @brief constructor that uses comoving coordinates
       *  @param xx comoving coordinate
       *  @param yy comoving coordinate
       *  @param zz comoving coordinate
       *  @param weight weight
       *  @param vx peculiar velocity
       *  @param vy peculiar velocity
       *  @param vz peculiar velocity    
       *  @param mass mass
       *  @param generic generic variable
       *  @return object of class Mock
       */
      Mock (const double xx, const double yy, const double zz, const double weight=1., const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble, const double generic=par::defaultDouble) 
	: Halo(xx, yy, zz, weight, vx, vy, vz, mass), m_generic(generic) {}

      /**
       *  @brief constructor that uses observed coordinates
       *  @param ra Right Ascension
       *  @param dec Declination
       *  @param redshift redshift
       *  @param cosm object of class Cosmology, used to estimate
       *  comoving distances
       *  @param weight weight
       *  @param vx peculiar velocity
       *  @param vy peculiar velocity
       *  @param vz peculiar velocity    
       *  @param mass mass
       *  @param generic generic variable
       *  @return object of class Mock
       */
      Mock (const double ra, const double dec, const double redshift, const Cosmology &cosm, const double weight=1., const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble, const double generic=par::defaultDouble) 
	: Halo(ra, dec, redshift, cosm, weight, vx, vy, vz, mass), m_generic(generic) {}

      /**
       *  @brief default destructor
       *  @return none
       */
      ~Mock () {}

      ///@}
  
  
      /**
       *  @name Member functions used to get the private members 
       */
      ///@{
    
      /**
       *  @brief get the private member Mock::m_generic
       *  @return the generic variable of the mock object
       */
      double generic () const override { return m_generic; }

      ///@}
  
  
      /**
       *  @name Member functions used to set the private members 
       */
      ///@{
    
      /**
       *  @brief set the private member Mock::m_generic
       *  @param generic the generic variable of the mock object
       *  @return none
       */
      void set_generic (const double generic=par::defaultDouble) override { m_generic = generic; }
    
      ///@}

    };
  }
}

#endif
