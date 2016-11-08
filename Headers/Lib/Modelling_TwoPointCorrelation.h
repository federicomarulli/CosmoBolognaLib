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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation.h
 *
 *  @brief The class Modelling_TwoPointCorrelation
 *
 *  This file defines the interface of the class Modelling, used to
 *  model two-point correlation functions of any kind
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLING2P__
#define __MODELLING2P__


#include "Cosmology.h"
#include "TwoPointCorrelation.h"
#include "Modelling.h"
#include "ModelFunction.h"


// ===================================================================================================


namespace cosmobl {

  namespace modelling {

    /**
     *  @class Modelling_TwoPointCorrelation Modelling_TwoPointCorrelation.h
     *  "Headers/Lib/Modelling_TwoPointCorrelation.h"
     *
     *  @brief The class Modelling_TwoPointCorrelation
     *
     *  This file defines the interface of the base class Modelling_TwoPointCorrelation,
     *  used for modelling any kind of measurements
     *
     */
    class Modelling_TwoPointCorrelation : public Modelling 
    {
    protected:
	
      /// the two-point correlation function type
      twopt::TwoPType m_twoPType;

      /// container of parameters for two-point correlation function model computation
      modelling::STR_twop_model m_twop_parameters;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class ModellingTwoPointCorrelation
       */
      Modelling_TwoPointCorrelation () = default;
	
      /**
       *  @brief constuctor
       *  @param twop the two-point correlation function to model
       *  @return object of class Modelling_TwoPointCorrelation
       */
      Modelling_TwoPointCorrelation (const shared_ptr<cosmobl::twopt::TwoPointCorrelation> twop)
	{ m_data = twop->dataset(); }
	
      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Modelling_TwoPointCorrelation () = default;

      /**
       *  @brief static factory used to construct modelling of 
       *  two-point correlation functions of any type
       *
       *  @param twop the two-point correlation function to model
       *
       *  @return a pointer to an object of class 
       *  Modelling_TwoPointCorrelation of a given type
       */
      static shared_ptr<Modelling_TwoPointCorrelation> Create (const shared_ptr<twopt::TwoPointCorrelation> twop);

      /**
       *  @brief static factory used to construct modelling of 
       *  two-point correlation functions of any type
       *
       *  @param twoPType type of the two-point correlation function
       *
       *  @param twop_dataset the dataset containing the two-point
       *  correlation function to model
       *
       *  @return a pointer to an object of class 
       *  Modelling_TwoPointCorrelation of a given type
       */
      static shared_ptr<Modelling_TwoPointCorrelation> Create (const twopt::TwoPType twoPType, const shared_ptr<data::Data> twop_dataset);

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
      virtual double alpha_bestfit () const { return m_model->parameter(0)->value(); }

      /**
       * @brief return the best-fit value of \f$f(z)\sigma_8(z)\f$
       *
       * @return the best-fit value of \f$f(z)\sigma_8(z)\f$
       */
      virtual double fsigma8_bestfit () const { return m_model->parameter(1)->value(); }

      /**
       * @brief return the best-fit value of \f$b(z)\sigma_8(z)\f$
       *
       * @return the best-fit value of \f$b(z)\sigma_8(z)\f$
       */
      virtual double bsigma8_bestfit () const { return m_model->parameter(2)->value(); }
             
      /**
       * @brief return the best-fit value of \f$A_0\f$
       *
       * @return the best-fit value of \f$A_0\f$
       */
      virtual double A0_bestfit () const { return m_model->parameter(3)->value(); }

      /**
       * @brief return the best-fit value of \f$A_1\f$
       *
       * @return the best-fit value of \f$A_1\f$
       */
      virtual double A1_bestfit () const { return m_model->parameter(4)->value(); }

      /**
       * @brief return the best-fit value of \f$A_2\f$
       *
       * @return the best-fit value of \f$A_2\f$
       */
      virtual double A2_bestfit () const { return m_model->parameter(5)->value(); }
	
      ///@}

	
      /**
       * @brief return the type of correlation function
       * @return the type of correlation function
       */
      twopt::TwoPType twoPType () { return m_twoPType; }

	
      /**
       *  @brief set the parameters for the computation of the dark
       *  matter two-point correlation function
       *
       *  @param fiducial_radDM scales at wich the fiducal model for
       *  &xi;<SUB>DM</SUB> is computed
       *
       *  @param cosmology the cosmological model used to compute
       *  &xi;<SUB>DM</SUB>
       *
       *  @param redshift redshift
       *
       *  @param method_Pk method used to compute the power
       *  spectrum; valid choices for method_Pk are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param sigmaNL damping of the wiggles in the linear power spectrum
       *
       *  @param NL 0 &rarr; linear power spectrum; 1 &rarr;
       *  non-linear power spectrum
       *
       *  @param pimax the upper limit of the line of sight integration
       *
       *  @param r_min minimum separation up to which the
       *  correlation function is computed
       *
       *  @param r_max maximum separation up to which the
       *  correlation function is computed
       *  
       *  @param output_root output_root of the parameter file used
       *  to compute the power spectrum and &sigma;(mass); it can be
       *  any name
       *  
       *  @param norm 0 &rarr; don't normalize the power spectrum; 1
       *  &rarr; normalize the power spectrum
       *  
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *  
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *  
       *  @param aa parameter \e a of Eq. 24 of Anderson et al. 2012
       *  
       *  @param GSL 0 &rarr; the FFTlog libraries are used; 1
       *  &rarr; the GSL libraries are used
       *  
       *  @param prec accuracy of the GSL integration
       *  
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return none
       */
      void set_parameters_xiDM (const vector<double> fiducial_radDM, const cosmology::Cosmology cosmology, const double redshift, const string method_Pk="CAMB", const double sigmaNL=0, const bool NL=true, const double pimax=40, const double r_min=1.e-3, const double r_max=350., const string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const double aa=0, const bool GSL=true, const double prec=1.e-3, const string file_par=par::defaultString);
	

      /**
       *  @brief set the parameters for the computation of the dark
       *  matter two-point correlation function
       *
       *  @param fiducial_radDM scales at wich the fiducal model for
       *  &xi;<SUB>DM</SUB> is computed
       *  
       *  @param cosmology the cosmology used
       *
       *  @param redshift redshift
       *
       *  @param method_Pk method used to compute the power spectrum
       *  and &sigma;(mass); valid choices for method_Pk are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param sigmaNL damping of the wiggles in the linear power
       *  spectrum
       *
       *  @param NL 0 &rarr; linear power spectrum; 1 &rarr;
       *  non-linear power spectrum
       *
       *  @param FV 0 &rarr; exponential form for f(v); 1 &rarr;
       *  Gaussian form for f(v); where f(v) is the velocity
       *  distribution function
       *
       *  @param Xi vector of &xi;(r), the two-point correlation
       *  function of dark matter
       *
       *  @param Xi_ vector of barred &xi;(r),
       *
       *  @param Xi__ vector of double-barred &xi;(r)
       *
       *  @param output_root output_root of the parameter file used
       *  to compute the power spectrum and &sigma;(mass); it can be
       *  any name
       *
       *  @param bias_nl 0 &rarr; linear bias; 1 &rarr; non-linear bias 
       *
       *  @param bA b<SUB>a</SUB> non-linear bias parameter
       *
       *  @param xiType 0 &rarr; standard; 1 &rarr; Chuang & Wang model
       *
       *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
       *
       *  @param xiNL 0 &rarr; linear power spectrum; 1 &rarr;
       *  non-linear power spectrum
       *
       *  @param v_min minimum velocity used in the convolution of
       *  the correlation function
       *
       *  @param v_max maximum velocity used in the convolution of
       *  the correlation function
       *
       *  @param step_v number of steps used in the convolution of
       *  the correlation function
       *
       *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
       *  normalize the power spectrum
       *
       *  @param r_min minimum separation up to which the
       *  correlation function is computed
       *
       *  @param r_max maximum separation up to which the
       *  correlation function is computed
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param aa parameter \e a of Eq. 24 of Anderson et al. 2012
       *
       *  @param GSL 0 &rarr; the Numerical libraries are used; 1
       *  &rarr; the GSL libraries are used
       *
       *  @param prec accuracy of the GSL integration
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return none
       */
      void set_parameters_xi2D_DM (const vector<double> fiducial_radDM, const cosmology::Cosmology cosmology, const double redshift, const string method_Pk="CAMB", const double sigmaNL=0, const bool NL=true, const int FV=0, const vector<double> Xi={}, const vector<double> Xi_={}, const vector<double> Xi__={}, const string output_root="test", const bool bias_nl=false, const double bA=-1., const bool xiType=false, const double k_star=-1., const bool xiNL=false, const double v_min=-5000., const double v_max=5000., const int step_v=500, const int norm=-1, const double r_min=0.1, const double r_max=150., const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=true, const double prec=1.e-2, const string file_par=par::defaultString);

	
      /**
       *  @brief set the fiducial model for dark matter two-point
       *  correlation function
       *
       *  @return none
       */
      virtual void set_fiducial_xiDM () = 0;

    };
  }
}

#endif
