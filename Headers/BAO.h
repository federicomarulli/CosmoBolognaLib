/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo  *
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
 *******************************************************************/

/**
 *  @file Headers/BAO.h
 *
 *  @brief The class BAO
 *
 *  This file defines the interface of the class BAO, used to model the features of the baryon acoustic oscillations
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __BAO__
#define __BAO__

#include "PkXi.h"


// ===================================================================================================


namespace cbl {
  
  namespace cosmology {

    /**
     *  @class BAO BAO.h "Headers/BAO.h"
     *
     *  @brief The class BAO
     *
     *  This class is used to handle objects of type <EM> BAO </EM>.
     */
    class BAO {

    private:

      /// pointer to the input cosmology
      std::shared_ptr<cosmology::Cosmology> m_cosmology = NULL;
      
    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      BAO () = default;
      
      /**
       *  @brief constructor
       *
       *  @param cosmology pointer to an object of class Cosmology
       */
      BAO (std::shared_ptr<cbl::cosmology::Cosmology> cosmology)
      : m_cosmology(std::move(cosmology)) {}
      
      /**
       *  @brief default destructor
       */
      ~BAO () = default;

      ///@}

      
      /**
       *  @name Functions to get the private members of the class
       */
      ///@{
      
      /**
       *  @brief Get the private member m_cosmology
       *
       *  @return the cosmological model
       */
      std::shared_ptr<cosmology::Cosmology> cosmology ()
      { return move(m_cosmology); }

      ///@}
      
      
      /**
       *  @name Functions used to model baryon acoustic oscillations
       */
      ///@{
      
      /**
       *  @brief redshift at wich occurs baryon photon decoupling
       *
       *  see Hu & Sugiyama (1996)
       *
       *  @return z<SUB>dec</SUB>
       */
      double z_decoupling () const;

      /**
       *  @brief redshift of drag epoch 
       *
       *  see Hu & Sugiyama (1996).
       *
       *  @return z<SUB>dec</SUB>
       */
      double z_drag () const;   

      /**
       *  @brief the sound horizon at the drag epoch
       *  r<SUB>s</SUB>(z<SUB>d</SUB>), valid choices for method_Pk
       *  are: EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html],
       *  CAMB [http://camb.info/]
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param method_Pk the method to compute the sound horizon
       *  
       *  @param T_CMB T<SUB>CMB</SUB>: the present day CMB temperature [K]
       *     
       *  @return r<SUB>s</SUB>
       */
      double rs (const std::string method_Pk, const double T_CMB=par::TCMB) const;

      /**
       *  @brief the sound horizon
       *
       *  @param redshift the redshift
       *
       *  @param T_CMB the temperature of the Cosmic Microwave
       *  Background
       *
       *  @return the sound horizon
       */
      double rs (const double redshift, const double T_CMB=2.7255) const;

      /**
       *  @brief the fiducial cosmology independent ratio
       *  r<SUB>s</SUB>/D<SUB>V</SUB>,  valid choices for method_Pk are: EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html],
       *  CAMB [http://camb.info/]
       *
       *  both r<SUB>s</SUB> and D<SUB>V</SUB> are in Mpc
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param redshift the redshift
       *
       *  @param method_Pk method used to compute the sound horizon
       *  (i.e. the Boltzmann solver)
       *
       *  @param T_CMB CMB temperature
       *
       *  @return y<SUB>s</SUB>
       */
      double ys (const double redshift, const std::string method_Pk, const double T_CMB=par::TCMB) const;
  
      /**
       *  @brief the acoustic parameter 
       *
       *  see Eisenstein 2005 
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *  @param redshift the redshift
       *  @return the acoustic parameter 
       */
      double Az (const double redshift) const;

      /**
       *  @brief the sound speed
       *
       *  @param redshift the redshift
       *
       *  @param T_CMB the temperature of the Cosmic Microwave
       *  Background
       *
       *  @return the sound speed
       */
      double sound_speed (const double redshift, const double T_CMB=2.7255) const;

      /**
       *  @brief the sound horizon integrand
       *
       *  @param redshift the redshift
       *
       *  @param T_CMB the temperature of the Cosmic Microwave
       *  Background
       *
       *  @return the sound horizon integrand
       */
      double rs_integrand (const double redshift, const double T_CMB=2.7255) const;

      /**
       *  @brief the linear point  
       *
       *  see Anselmi et al. 2016 
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param redshift the redshift
       *
       *  @param rmin the minimum scale
       *
       *  @param rmax the maximum scale
       *
       *  @param nbinr the number of scale bins
       *
       *  @param interpType the interpolation type
       *
       *  @return vector containing the linear point, the dip 
       *  and the BAO peak for the correlation function at the redshift
       *  provided
       *
       */
      std::vector<double> linear_point (const double redshift, const double rmin=60., const double rmax=150., const int nbinr=100, const std::string interpType="Spline");
      
      ///@}

    };
  }
}

#endif
