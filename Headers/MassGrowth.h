/********************************************************************
 *  Copyright (C) 2023 by Federico Marulli                          *
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
 *  @file Headers/MassGrowth.h
 *
 *  @brief The class MassGrowth
 *
 *  This file defines the interface of the class MassGrowth, used to
 *  model the cosmic mass accretion history
 *
 *  @author Federico Marulli, Carlo Giocoli  
 *
 *  @author federico.marulli3@unibo.it, carlo.giocoli@inaf.it
 */

#ifndef __MASSGROWTH__
#define __MASSGROWTH__

#include "Sigma.h"


// ===================================================================================================


namespace cbl {

  namespace cosmology {

    /**
     *  @class MassGrowth MassGrowth.h "Headers/MassGrowth.h"
     *
     *  @brief The class MassGrowth
     *
     *  This class is used to handle objects of type <EM> MassGrowth
     *  </EM>. It is used to model the cosmic mass accretion history
     */
    class MassGrowth {

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
      MassGrowth () = default;

      /**
       *  @brief constructor
       *
       *  @param cosmology pointer to an object of class Cosmology
       */
      MassGrowth (std::shared_ptr<cbl::cosmology::Cosmology> cosmology)
      : m_cosmology(std::move(cosmology)) {}
      
      /**
       *  @brief default destructor
       */
      ~MassGrowth () = default;

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
       *  @name Functions to estimate the cosmic mass accretion
       *  history
       */
      ///@{

      /**
       *  @brief The differential distribution 
       *
       *  This member function provides the differential rescaled and
       *  generalized formation redshift distributions
       *
       *  @author Carlo Giocoli
       *  @author carlo.giocoli@inaf.it
       *
       *  @param ww rescaled variable w as in Lacey and Coles 1993
       *
       *  @param ff assembled fraction
       *
       *  @param author valid authors are: NS (Nusser and Sheth), GTS
       *  (Giocoli et al. 2012)
       *
       *  @return p(w): the differential distribution
       */
      double pw (const double ww, const double ff, const std::string author) const; 

      /**
       *  @brief The formation probability
       *
       *  This member function provides the probability that a halo of
       *  a given mass m0 at redshift z0 makes a mass fraction f at
       *  redshift z
       *
       *  @author Carlo Giocoli
       *  @author carlo.giocoli@inaf.it
       *
       *  @param m0 halo mass
       *
       *  @param z0 redshift when the halo has a mass m0
       *
       *  @param frac mass fraction
       *
       *  @param redshift the redshift
       *
       *  @param model_model valid authors are: NS (Nusser and Sheth),
       *  GTS (Giocoli et al. 2012)
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  CLASS [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *     
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @return p(z): the formation probability
       */
      double pz (const double m0, const double z0, const double frac, const double redshift, const std::string model_model, const std::string method_SS, const bool store_output=true, const std::string output_root="test") const; 

      /**
       *  @brief The cumulative distribution 
       *
       *  This member function provides the cumulative rescaled and generalized
       *  formation redshift distributions
       *
       *  @author Carlo Giocoli
       *  @author carlo.giocoli@inaf.it
       *
       *  @param ww rescaled variable w as in Lacey and Coles 1993
       *
       *  @param ff assembled fraction
       *
       *  @param author valid authors are: NS (Nusser and Sheth), GTS
       *  (Giocoli et al. 2012)
       *
       *  @return P(w): the cumulative distribution
       */
      double cumPw (const double ww, const double ff, const std::string author) const; 
    
      /**
       *  @brief The median formation w   
       *  
       *  @author Carlo Giocoli
       *  @author carlo.giocoli@inaf.it
       *
       *  @param [in] ff assembled fraction
       *
       *  @param [in] model_model valid authors are: NS (Nusser and
       *  Sheth), GTS (Giocoli et al. 2012)
       *
       *  @param [out] wf vector of w(f)
       */
      void medianwf (const double ff, const std::string model_model, std::vector<double> &wf) const; 

      /**
       *  @brief The median formation z 
       *
       *  @author Carlo Giocoli
       *  @author carlo.giocoli@inaf.it
       *
       *  @param [in] ff assembled fraction
       *
       *  @param [in] mass halo mass
       *
       *  @param [in] z0 redshift when the halo has a mass mass
       *
       *  @param [in] model_model valid authors are: NS (Nusser and
       *  Sheth), GTS (Giocoli et al. 2012)
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  CLASS [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param [out] zf vector of z(f)
       *      
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       */
      void medianzf (const double ff, const double mass, const double z0, const std::string model_model, const std::string method_SS, std::vector<double> &zf, const bool store_output=true, const std::string output_root="test") const; 
  
      /**
       *  @brief The rescaled variable w 
       *
       *  This member functions provides the conditional variable
       *  w=[&delta;<SUB>c</SUB>(zf) - &delta;<SUB>c</SUB>(z)] /
       *  &radic; [s(fm)-s(m)] where&delta;<SUB>c</SUB>(z) =
       *  &delta;<SUB>c0</SUB>(z)/D+(z), as in Lacey and Coles 1993
       *
       *  @author Carlo Giocoli
       *  @author carlo.giocoli@inaf.it
       *
       *  @param mm halo mass
       *
       *  @param redshift the redshift 
       *
       *  @param ff assembled fraction
       *
       *  @param zf redshift at which the mass is accreted
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  CLASS [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *     
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @return the conditional variable w
       */
      double wf (const double mm, const double redshift, const double ff, const double zf, const std::string method_SS, const bool store_output=true, const std::string output_root="test") const; 

      /**
       *  @brief The redshift at a given wf
       *
       *  This member function estimates the redshift from wf, given
       *  the parent halo mass at z=z', z', and its assembled fraction
       *  f
       *
       *  @author Carlo Giocoli
       *  @author cgiocoli@gmail.com
       *
       *  @param mm mass
       *
       *  @param redshift the redshift
       *
       *  @param ff assembled fraction
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass) const; valid method_SS are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param wwf rescaled variable w as in Lacey and Coles 1993
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @return the redshift
       */
      double Redshift (const double mm, const double redshift, const double ff, const std::string method_SS, const double wwf, const bool store_output=true, const std::string output_root="test") const;
      
      ///@}
      
    };
  }
}

#endif
