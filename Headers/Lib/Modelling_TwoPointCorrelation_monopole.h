/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  Modelling_TwoPointCorrelation_monopole, used to model the monopole
 *  of two-point correlation function
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGMONO__
#define __MODELLINGMONO__


#include "Modelling_TwoPointCorrelation1D.h"


// ===================================================================================================


namespace cosmobl {
  
  namespace modelling {
    
    /**
     *  @class Modelling_TwoPointCorrelation_monopole
     *  Modelling_TwoPointCorrelation_monopole.h
     *  "Headers/Lib/Modelling_TwoPointCorrelation_monopole.h"
     *
     *  @brief The class Modelling_TwoPointCorrelation_monopole
     *
     *  This file defines the interface of the base class
     *  Modelling_TwoPointCorrelation_monopole, used for modelling the
     *  monopole of two-point correlation function
     *
     */
    class Modelling_TwoPointCorrelation_monopole : public Modelling_TwoPointCorrelation1D {

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class
       *  Modelling_TwoPointCorrelation_monopole
       */
      Modelling_TwoPointCorrelation_monopole () = default;

      /**
       *  @brief constructor
       *  
       *  @param twop the two-point correlation function to model
       *
       *  @return object of type
       *  Modelling_TwoPointCorrelation_monopole
       */
      Modelling_TwoPointCorrelation_monopole (const shared_ptr<cosmobl::twopt::TwoPointCorrelation> twop)
	: Modelling_TwoPointCorrelation1D(twop) {}
	
      /**
       *  @brief constructor
       *  
       *  @param twop_dataset the dataset containing the two-point
       *  correlation function to model
       *
       *  @return object of type
       *  Modelling_TwoPointCorrelation_monopole
       */
      Modelling_TwoPointCorrelation_monopole (const shared_ptr<data::Data> twop_dataset)
	: Modelling_TwoPointCorrelation1D() { set_data(twop_dataset); }

      
      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Modelling_TwoPointCorrelation_monopole () = default;
	
      ///@}

      
      /**
       *  @name Member functions used to get the best-fit values of
       *  model parameters
       */
      ///@{

      /**
       * @brief return the best-fit value of \f$\alpha\f$
       *
       * @return the best-fit value of \f$\alpha\f$
       */
      double alpha_bestfit () const { return m_model->parameter(0)->value(); }

      /**
       * @brief return the best-fit value of \f$f(z)\sigma_8(z)\f$
       *
       * @return the best-fit value of \f$f(z)\sigma_8(z)\f$
       */
      double fsigma8_bestfit () const { return m_model->parameter(1)->value(); }

      /**
       * @brief return the best-fit value of \f$b(z)\sigma_8(z)\f$
       *
       * @return the best-fit value of \f$b(z)\sigma_8(z)\f$
       */
      double bsigma8_bestfit () const { return m_model->parameter(2)->value(); }
             
      /**
       * @brief return the best-fit value of \f$A_0\f$
       *
       * @return the best-fit value of \f$A_0\f$
       */
      double A0_bestfit () const { return m_model->parameter(3)->value(); }

      /**
       * @brief return the best-fit value of \f$A_1\f$
       *
       * @return the best-fit value of \f$A_1\f$
       */
      double A1_bestfit () const { return m_model->parameter(4)->value(); }

      /**
       * @brief return the best-fit value of \f$A_2\f$
       *
       * @return the best-fit value of \f$A_2\f$
       */
      double A2_bestfit () const { return m_model->parameter(5)->value(); }
      
      ///@}

      
      /**
       * @brief set the fiducial model for dark matter two-point
       * correlation function
       *
       *  @return none
       */
      virtual void set_fiducial_xiDM () override;


      /**
       *  @brief set the parameters to model the monopole of the
       *  two-point correlation function in redshift space
       * 
       *  redshift-space distorsions are modelled in the Kaiser limit,
       *  that is neglecting non-linearities in dynamics and bias;
       *  specifically, the model considered is the following:
       *  
       *  \f$\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3} f\sigma_8
       *  \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2 \right] \cdot
       *  \xi_{\rm DM}(s)/\sigma_8^2\f$
       *
       *  @param fsigma8_guess guess value for the parameter
       *  \f$f(z)\sigma_8(z)\f$
       * 
       *  @param fsigma8_prior prior for the parameter
       *  \f$f(z)\sigma_8(z)\f$
       *
       *  @param fsigma8_type type of the parameter
       *  \f$f(z)\sigma_8(z)\f$: it can be either _free_ or _fixed_
       *
       *  @param bsigma8_guess guess value for the parameter
       *  \f$b(z)\sigma_8(z)\f$
       *
       *  @param bsigma8_prior prior for the parameter
       *  \f$b(z)\sigma_8(z)\f$
       *
       *  @param bsigma8_type type of the parameter
       *  \f$b(z)\sigma_8(z)\f$: it can be either _free_ or _fixed_
       *
       *  @return none
       */
      void set_model_Kaiser (const double fsigma8_guess, const statistics::Prior fsigma8_prior, const statistics::ParameterType fsigma8_type, const double bsigma8_guess, const statistics::Prior bsigma8_prior, const statistics::ParameterType bsigma8_type)
      { set_model_monopole(1., {}, statistics::_fixed_, fsigma8_guess, fsigma8_prior, fsigma8_type, bsigma8_guess, bsigma8_prior, bsigma8_type, 0., {}, statistics::_fixed_, 0., {}, statistics::_fixed_, 0., {}, statistics::_fixed_); }

      
      /**
       *  @brief overloading of the function used to set the
       *  parameters to model the monopole of the two-point
       *  correlation function in redshift space
       * 
       *  redshift-space distorsions are modelled in the Kaiser limit,
       *  that is neglecting non-linearities in dynamics and bias;
       *  specifically, the model considered is the following:
       *  
       *  \f$\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3} f\sigma_8
       *  \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2 \right] \cdot
       *  \xi_{\rm DM}(s)/\sigma_8^2\f$
       *
       *  inizial guess values for \f$f(z)\sigma_8(z)\f$ and
       *  \f$b(z)\sigma_8(z)\f$ are extracted from the prior
       *  distributions; both the two parameters are free
       *
       *  @param fsigma8_prior prior for the parameter
       *  \f$f(z)\sigma_8(z)\f$
       *
       *  @param bsigma8_prior prior for the parameter
       *  \f$b(z)\sigma_8(z)\f$
       *
       *  @return none
       */
      void set_model_Kaiser (const statistics::Prior fsigma8_prior, const statistics::Prior bsigma8_prior)
      { set_model_Kaiser(par::defaultDouble, fsigma8_prior, statistics::_free_, par::defaultDouble, bsigma8_prior, statistics::_free_); }
      
      
      /**
       *  @brief set the parameter to model the monopole of the
       *  two-point correlation function in real space, taking into
       *  accout geometric distortions (that is the Alcock-Paczynski
       *  effect)
       *
       *  the model used is the following:
       *
       *  \f$\xi(s)= b^2 \xi_{DM}(\alpha s)\ + A_0 + A_1/s +A_2/s^2\f$
       *
       *  where \f$\xi_{DM}\f$ is computed at the fiducial (fixed)
       *  cosmology, and {\f$b\sigma_8\f$, \f$A_0\f$, \f$A_1\f$,
       *  \f$A_2\f$} are considered as nuisance parameters
       *
       *  @param alpha_guess guess value for the parameter \f$\alpha\f$
       *
       *  @param alpha_prior prior for the parameter \f$\alpha\f$
       *
       *  @param alpha_type type of the parameter \f$\alpha\f$: it
       *  can be either _free_ or _fixed_
       *
       *  @param bsigma8_guess guess value for the parameter
       *  \f$b(z)\sigma_8(z)\f$
       *
       *  @param bsigma8_prior prior for the parameter
       *  \f$b(z)\sigma_8(z)\f$
       *
       *  @param bsigma8_type type of the parameter
       *  \f$b(z)\sigma_8(z)\f$: it can be either _free_ or _fixed_
       *
       *  @param A0_guess guess value for the parameter \f$A_0\f$
       *
       *  @param A0_prior prior for the parameter \f$A_0\f$
       *
       *  @param A0_type type of the parameter \f$A_0\f$: it can be
       *  either _free_ or _fixed_
       *
       *  @param A1_guess guess value for the parameter \f$A_1\f$
       *
       *  @param A1_prior prior for the parameter \f$A_1\f$
       *
       *  @param A1_type type of the parameter \f$A_1\f$: it can be
       *  either _free_ or _fixed_
       *
       *  @param A2_guess guess value for the parameter \f$A_2\f$
       *
       *  @param A2_prior prior for the parameter \f$A_2\f$
       *
       *  @param A2_type type of the parameter \f$A_2\f$: it can be
       *  either _free_ or _fixed_
       *
       *  @return none
       */
      void set_model_BAO (const double alpha_guess=par::defaultDouble, const statistics::Prior alpha_prior={}, const statistics::ParameterType alpha_type=statistics::_free_, const double bsigma8_guess=par::defaultDouble, const statistics::Prior bsigma8_prior={}, const statistics::ParameterType bsigma8_type=statistics::_free_, const double A0_guess=par::defaultDouble, const statistics::Prior A0_prior={}, const statistics::ParameterType A0_type=statistics::_free_, const double A1_guess=par::defaultDouble, const statistics::Prior A1_prior={}, const statistics::ParameterType A1_type=statistics::_free_, const double A2_guess=par::defaultDouble, const statistics::Prior A2_prior={}, const statistics::ParameterType A2_type=statistics::_free_)
      {
	set_model_monopole(alpha_guess, alpha_prior, alpha_type, 0., {}, statistics::_fixed_, bsigma8_guess, bsigma8_prior, bsigma8_type, A0_guess, A0_prior, A0_type, A1_guess, A1_prior, A1_type, A2_guess, A2_prior, A2_type);
      }


      /**
       *  @brief overloading of the function used to set the
       *  parameters to model the monopole of the two-point
       *  correlation function in real space, taking into accout
       *  geometric distortions (that is the Alcock-Paczynski effect)
       *
       *  the model used is the following:
       *
       *  \f$\xi(s)= b^2 \xi_{DM}(\alpha s)\ + A_0 + A_1/s +A_2/s^2\f$
       *
       *  where \f$\xi_{DM}\f$ is computed at the fiducial (fixed)
       *  cosmology, and {\f$b\sigma_8\f$, \f$A_0\f$, \f$A_1\f$,
       *  \f$A_2\f$} are considered as nuisance parameters
       *
       *  @param alpha_prior prior for the parameter \f$\alpha\f$
       *
       *  @param bsigma8_prior prior for the parameter
       *  \f$b(z)\sigma_8(z)\f$
       *
       *  @param A0_prior prior for the parameter \f$A_0\f$
       *
       *  @param A1_prior prior for the parameter \f$A_1\f$
       *
       *  @param A2_prior prior for the parameter \f$A_2\f$
       *
       *  @return none
       */
      void set_model_BAO (const statistics::Prior alpha_prior, const statistics::Prior bsigma8_prior, const statistics::Prior A0_prior, const statistics::Prior A1_prior, const statistics::Prior A2_prior)
      { set_model_monopole(par::defaultDouble, alpha_prior, statistics::_free_, 0., {}, statistics::_fixed_, par::defaultDouble, bsigma8_prior, statistics::_free_, par::defaultDouble, A0_prior, statistics::_free_, par::defaultDouble, A1_prior, statistics::_free_, par::defaultDouble, A2_prior, statistics::_free_); }
      
    };
  }
}

#endif
