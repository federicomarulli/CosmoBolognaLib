/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation_monopole.h
 *
 *  @brief The class Modelling_TwoPointCorrelation_monopole
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation_monopole, used for modelling monopole of 2pcf
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __MODELLINGMONO__
#define __MODELLINGMONO__


#include "Modelling_TwoPointCorrelation.h"


// ===================================================================================================


namespace cosmobl {
  
  namespace modelling {
    
    /**
     *  @class Modelling_TwoPointCorrelation_monopole Modelling_TwoPointCorrelation_monopole.h
     *  "Headers/Lib/Modelling_TwoPointCorrelation_monopole.h"
     *
     *  @brief The class Modelling_TwoPointCorrelation_monopole
     *
     *  This file defines the interface of the base class Modelling_TwoPointCorrelation_monopole,
     *  used for modelling the monopole of 2pcf
     *
     */
    class Modelling_TwoPointCorrelation_monopole : public Modelling_TwoPointCorrelation {

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Modelling_TwoPointCorrelation_monopole
       */
      Modelling_TwoPointCorrelation_monopole () = default;

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Modelling_TwoPointCorrelation_monopole () = default;

      /**
       *  @brief constructor of the ModellingTwoPointCorrelation_monopole
       *  
       *  @param twop the two-point correlation function to model
       *
       *  @return object of type Modelling_TwoPointCorrelation_monopole
       */
      Modelling_TwoPointCorrelation_monopole (const shared_ptr<cosmobl::twopt::TwoPointCorrelation> twop);
	
      ///@}

      /**
       * @brief set the fiducial model for dark matter 
       * two point correlation function
       *
       *  @return none
       */
      void set_fiducial_twop () override;

      /**
       * @brief fit the monopole of the two-point correlation function
       * taking into accout geometric distortions (i.e. the
       * Alcock-Paczynski effect). The model used is \f$\xi(s)= B^2
       * \xi_{DM}(\alpha s)\ + A_0 + A_1/s +A_2/s^2\f$, where
       * \f$\xi_{DM}\f$ is computed at the fiducial cosmology, and {B,
       * A<SUB>0</SUB>, A<SUB>1</SUB>, A<SUB>2</SUB>} are considered
       * as nuisance parameters
       *
       * @param alpha_prior the prior for &alpha;
       *
       * @param B_prior the prior for B
       *
       * @param A0_prior the prior for A0
       *
       * @param A1_prior the prior for A1
       *
       * @param A2_prior the prior for A2
       *
       * @param pT_alpha the parameter type of &alpha;: it can be
       * either statistics::_free_ or _fixed_
       *
       * @param pT_B the parameter type of B: it can be either
       * statistics::_free_ or _fixed_
       *
       * @param pT_A0 the parameter type of A0: it can be either
       * statistics::_free_ or _fixed_
       *
       * @param pT_A1 the parameter type of A1: it can be either
       * statistics::_free_ or _fixed_
       *
       * @param pT_A2 the parameter type of A2: it can be either
       * statistics::_free_ or _fixed_
       *
       * @return none
       */
      void set_model_AP_isotropic (const statistics::Prior alpha_prior, const statistics::Prior B_prior, const statistics::Prior A0_prior, const statistics::Prior A1_prior, const statistics::Prior A2_prior, const statistics::ParameterType pT_alpha=statistics::_free_, const statistics::ParameterType pT_B=statistics::_free_, const statistics::ParameterType pT_A0=statistics::_free_, const statistics::ParameterType pT_A1=statistics::_free_, const statistics::ParameterType pT_A2=statistics::_free_) override;

      /**
       * @brief compute and write the model using the stored 
       * parameter values
       *
       * @param xx vector of point at which the model
       * is computed
       * @param dir_model the output directory of the model
       * @param file_model the name of the file
       *
       * @return none
       */
       void write_model(const vector<double> xx, const string dir_model, const string file_model)
       { m_model->write_model(xx, dir_model, file_model); }

      /**
       * @brief compute and write the model using the stored 
       * parameter values
       *
       * @param xx vector of point at which the model
       * is computed
       * @param parameters vector of parameters values
       * at which the model is computed
       * @param dir_model the output directory of the model
       * @param file_model the name of the file
       *
       * @return none
       */
       virtual void write_model_parameters(const vector<double> xx, const vector<double> parameters, const string dir_model, const string file_model)
       { m_model->write_model(xx, parameters, dir_model, file_model); }
    };
  }
}

#endif
