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
 *  @file Headers/Velocities.h
 *
 *  @brief The class Velocities
 *
 *  This file defines the interface of the class Velocities, used for
 *  calculations of the statistics of the cosmic velocity field
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __VELOCITIES__
#define __VELOCITIES__

#include "PkXi.h"


// ===================================================================================================


namespace cbl {

  namespace cosmology {

    /**
     *  @class Velocities Velocities.h "Headers/Velocities.h"
     *
     *  @brief The class Velocities
     *
     *  This class is used to handle objects of type <EM> velocities
     *  </EM>.
     */
    class Velocities {

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
      Velocities () = default;
      
      /**
       *  @brief constructor
       *
       *  @param cosmology pointer to an object of class Cosmology
       */
      Velocities (std::shared_ptr<cbl::cosmology::Cosmology> cosmology)
      : m_cosmology(std::move(cosmology)) {}
      
      /**
       *  @brief default destructor
       */
      ~Velocities () = default;

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
       *  @name functions used to compute the statistics of the cosmic
       *  velocity field
       */
      ///@{
      
      /**
       *  @brief the mean square bulk flow
       *
       *  @param rr comoving radius 
       *
       *  @param k_int_min minimum wave vector module up to which the
       *  integral is computed
       *
       *  @param method_Pk method used to compute the power spectrum
       *  and &sigma;(mass) (i.e. the Boltzmann solver); valid choices
       *  for method_Pk are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
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
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return the mean square bulk flow
       */
      double square_bulk_flow (const double rr, const double k_int_min, const std::string method_Pk, const double redshift, const bool store_output=true, const std::string output_root="test", const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString);
                
      /**
       *  @brief the mean square bulk flow
       *
       *  @param rr comoving radius 
       *
       *  @param k_int_min minimum wave vector module up to which the
       *  integral is computed
       *
       *  @param lgkk vector of log(k)
       *
       *  @param lgPk vector of log(P(k))
       *
       *  @param redshift the redshift
       *
       *  @return the mean square bulk flow
       */
      double square_bulk_flow_Table (const double rr, const double k_int_min, const std::vector<double> lgkk, const std::vector<double> lgPk, const double redshift) const; 
                
      /**
       *  @brief the mean square velocity dispersion
       *
       *  @param rr comoving radius 
       *
       *  @param k_int_min minimum wave vector module up to which the
       *  integral is computed
       *
       *  @param method_Pk method used to compute the power spectrum
       *  and &sigma;(mass) (i.e. the Boltzmann solver); valid choices
       *  for method_Pk are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
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
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return the mean square velocity dispersion
       */
      double square_velocity_dispersion (const double rr, const double k_int_min, const std::string method_Pk, const double redshift, const bool store_output=true, const std::string output_root="test", const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString);
                
      /**
       *  @brief the Cosmic Mach Number
       *
       *  @param rr comoving radius 
       *
       *  @param k_int_min minimum wave vector module up to which the
       *  integral is computed
       *
       *  @param method_Pk method used to compute the power spectrum
       *  and &sigma;(mass) (i.e. the Boltzmann solver); valid choices
       *  for method_Pk are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
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
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return the Cosmic Mach Number
       */
      double CMN (const double rr, const double k_int_min, const std::string method_Pk, const double redshift, const bool store_output=true, const std::string output_root="test", const double k_max=100., const std::string file_par=par::defaultString) const;

      ///@}
      
    };
  }
}

#endif
