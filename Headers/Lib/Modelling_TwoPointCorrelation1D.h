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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation1D.h
 *
 *  @brief The class Modelling_TwoPointCorrelation1D
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation1D, that contains all the common
 *  methods to model 1D two-point correlation functions
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLING1D__
#define __MODELLING1D__


#include "Modelling_TwoPointCorrelation.h"


// ===================================================================================================


namespace cosmobl {
  
  namespace modelling {
    
    /**
     *  @class Modelling_TwoPointCorrelation1D
     *  Modelling_TwoPointCorrelation1D.h
     *  "Headers/Lib/Modelling_TwoPointCorrelation1D.h"
     *
     *  @brief The class Modelling_TwoPointCorrelation1D
     *
     *  This file defines the interface of the base class
     *  Modelling_TwoPointCorrelation1D, that contains all the common
     *  methods to model 1D two-point correlation functions
     *
     */
    class Modelling_TwoPointCorrelation1D : public Modelling_TwoPointCorrelation {

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class
       *  Modelling_TwoPointCorrelation1D
       */
      Modelling_TwoPointCorrelation1D () = default;

      /**
       *  @brief constructor
       *  
       *  @param twop the two-point correlation function to model
       *
       *  @return object of type
       *  Modelling_TwoPointCorrelation1D
       */
      Modelling_TwoPointCorrelation1D (const shared_ptr<cosmobl::twopt::TwoPointCorrelation> twop)
	: Modelling_TwoPointCorrelation(twop) {}
	
	/**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Modelling_TwoPointCorrelation1D () = default;
	
      ///@}

      
      /**
       *  @brief set the fiducial model for dark matter two-point
       *  correlation function
       *
       *  @return none
       */
      virtual void set_fiducial_xiDM () = 0;

      
      /**
       *  @brief set the parameters used to model the monopole of
       *  the two-point correlation function 
       *
       *  the model is the following:
       *
       *  \f$\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3} f\sigma_8
       *  \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2 \right] \cdot
       *  \xi_{\rm DM}(\alpha\cdot s)/\sigma_8^2 + A_0 + A_1/s +
       *  A_2/s^2\f$
       *
       *  the model has 6 parameters: 
       *    - \f$\alpha\f$
       *    - \f$f(z)\sigma_8(z)\f$
       *    - \f$b(z)\sigma_8(z)\f$
       *    - \f$A_0\f$
       *    - \f$A_1\f$
       *    - \f$A_2\f$ 
       *
       *  the dark matter two-point correlation function is computed
       *  using the input cosmological parameters
       *
       *  @param alpha_guess guess value for the parameter \f$\alpha\f$
       *
       *  @param alpha_prior prior for the parameter \f$\alpha\f$
       *
       *  @param alpha_type type of the parameter \f$\alpha\f$: it
       *  can be either _free_ or _fixed_
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
      void set_model_monopole (const double alpha_guess=par::defaultDouble, const statistics::Prior alpha_prior={}, const statistics::ParameterType alpha_type=statistics::_free_, const double fsigma8_guess=par::defaultDouble, const statistics::Prior fsigma8_prior={}, const statistics::ParameterType fsigma8_type=statistics::_free_, const double bsigma8_guess=par::defaultDouble, const statistics::Prior bsigma8_prior={}, const statistics::ParameterType bsigma8_type=statistics::_free_, const double A0_guess=par::defaultDouble, const statistics::Prior A0_prior={}, const statistics::ParameterType A0_type=statistics::_free_, const double A1_guess=par::defaultDouble, const statistics::Prior A1_prior={}, const statistics::ParameterType A1_type=statistics::_free_, const double A2_guess=par::defaultDouble, const statistics::Prior A2_prior={}, const statistics::ParameterType A2_type=statistics::_free_);

      
      /**
       *  @brief set the parameters to model the monopole of the
       *  two-point correlation function in real space assuming a
       *  linear bias
       *
       *  @param bsigma8_guess guess value for the parameter
       *  \f$b(z)\sigma_8(z)\f$
       *
       *  @param bsigma8_prior prior for the parameter
       *  \f$b(z)\sigma_8(z)\f$
       *
       *  @return none
       */
      void set_model_linearBias (const double bsigma8_guess, const statistics::Prior bsigma8_prior)
      {
	set_model_monopole(1., {}, statistics::_fixed_, 0., {}, statistics::_fixed_, bsigma8_guess, bsigma8_prior, statistics::_free_, 0., {}, statistics::_fixed_, 0., {}, statistics::_fixed_, 0., {}, statistics::_fixed_);
      }

      
      /**
       *  @brief set the parameters to model the monopole of the
       *  two-point correlation function in real space assuming a
       *  linear bias
       *
       *  @param bsigma8_prior prior for the parameter
       *  \f$b(z)\sigma_8(z)\f$
       *
       *  @return none
       */
      void set_model_linearBias (const statistics::Prior bsigma8_prior)
      {
	set_model_monopole(1., {}, statistics::_fixed_, 0., {}, statistics::_fixed_, par::defaultDouble, bsigma8_prior, statistics::_free_, 0., {}, statistics::_fixed_, 0., {}, statistics::_fixed_, 0., {}, statistics::_fixed_);
      }
      
      
      /**
       *  @name Member functions used to write the outputs
       */
      ///@{
      
      /**
       *  @brief compute and write the model
       *
       *  @param xx vector of points at which the model is computed
       *  @param dir the output directory
       *  @param file the output file
       *
       *  @param parameter vector containing the input parameters used
       *  to compute the model; if it is not provided, internal
       *  parameters will be used
       *
       *  @return none
       */
      void write_model (const vector<double> xx, const string dir, const string file, const vector<double> parameter={}) const override
      { m_model->write_model(xx, dir, file, parameter); }
      
      ///@}
      
    };
  }
}

#endif
