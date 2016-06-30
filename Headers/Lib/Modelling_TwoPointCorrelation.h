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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation.h
 *
 *  @brief The class Modelling_TwoPointCorrelation
 *
 *  This file defines the interface of the class
 *  Modelling, used for modelling any kind of 
 *  2pcf measurements
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
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
	
	/// two-point correlation function type
	twopt::TwoPType m_twoPType;

	/// container of parameters for two point model computation
	glob::STR_twop_model m_twop_parameters;

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class ModellingTwoPointCorrelation
	 */
	Modelling_TwoPointCorrelation () {}

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_TwoPointCorrelation () {}

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
	 *  @param twop_dataset the dataset containing
	 *  the two-point correlation function to model
	 *
	 *  @param twoPType the type of two-point correlation function
	 *
	 *  @return a pointer to an object of class 
	 *  Modelling_TwoPointCorrelation of a given type
	 */
	static shared_ptr<Modelling_TwoPointCorrelation> Create (const shared_ptr<data::Data> twop_dataset, const twopt::TwoPType twoPType);

	///@}

	/**
	 * @brief return the type of correlation function
	 * @return the type of correlation function
	 */
	twopt::TwoPType twoPType () {return m_twoPType;}


	/**
	 * @brief set the parameters for the computation
	 * of the dark matter two point correlation function
	 *
	 *  @param model_scales scales at wich the fiducal model for &xi;<SUB>DM</SUB>
	 *  @param cosmology the cosmology used
	 *  @param redshift redshift
	 *  @param method method used to compute the power spectrum; valid
	 *  choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
	 *  [http://class-code.net/], MPTbreeze-v1
	 *  [http://arxiv.org/abs/1207.1465], EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *  @param sigmaNL damping of the wiggles in the linear power spectrum
	 *  @param NL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power
	 *  spectrum
	 *  @param pimax the upper limit of the line of sight integration
	 *  @param r_min minimum separation up to which the
	 *  correlation function is computed
	 *  @param r_max maximum separation up to which the
	 *  correlation function is computed
	 *  @param output_root output_root of the parameter file used to
	 *  compute the power spectrum and &sigma;(mass); it can be any
	 *  name
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *  @param k_min minimum wave vector module up to which the power
	 *  spectrum is computed
	 *  @param k_max maximum wave vector module up to which the power
	 *  spectrum is computed
	 *  @param aa parameter \e a of Eq. 24 of Anderson et al. 2012
	 *  @param GSL 0 &rarr; the FFTlog libraries are used; 1 &rarr;
	 *  the GSL libraries are used
	 *  @param prec accuracy of the GSL integration 
	 *  @param file_par name of the parameter file; if a parameter
	 *  file is provided (i.e. file_par!=NULL), it will be used,
	 *  ignoring the cosmological parameters of the object
	 *
	 *  @return none
	 */
	void set_parameters_twop_DM(const vector<double> model_scales, const cosmology::Cosmology cosmology, const double redshift, const string method="CAMB", const double sigmaNL=0, const bool NL=1, const double pimax=40, const double r_min=1.e-3, const double r_max=350., const string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const double aa=0, const bool GSL=1, const double prec=1.e-3, const string file_par=par::defaultString);

	/**
	 * @brief set the fiducial model for dark matter 
	 * two point correlation function
	 *
	 *  @return none
	 */
	virtual void set_fiducial_twop ()
	{ ErrorMsg("Error in set_fiducial_twop() of Modelling_TwoPointCorrelation.h!"); }

	/**
	 * @brief fit the bias of the measured
	 * two point correlation function
	 *
	 * @param bias_prior prior for the bias
	 *
	 * @param pT_bias the parameter type of the bias: it can be
	 * either _free_ or _fixed_
	 *
	 * @return none
	 */
	virtual void set_model_bias (const statistics::Prior bias_prior, const statistics::ParameterType pT_bias=statistics::_free_)
	{ ErrorMsg("Error in fit_bias of Modelling_TwoPointCorrelation.h!"); }

	/**
	 * @brief fit the Alcock-Paczynski effect for the two point
	 * correlation function. 
	 * \f$\xi(s)= b^2  \xi_{DM}(\alpha s)\f$.
	 * &xi<SUB>DM</SUB> is computed at the fiducial cosmology.
	 *
	 * @param bias_prior prior for the bias
	 *
	 * @param alpha_prior prior for &alpha;
	 *
	 * @param pT_bias the parameter type of the bias: it can be
	 * either _free_ or _fixed_
	 *
	 * @param pT_alpha the parameter type of &alpha;: it can be
	 * either _free_ or _fixed_
	 *
	 * @return none
	 */
	virtual void set_model_bias_AP_isotropic(const statistics::Prior bias_prior, const statistics::Prior alpha_prior, const statistics::ParameterType pT_bias=statistics::_free_, const statistics::ParameterType pT_alpha=statistics::_free_)
	{ ErrorMsg("Error in fit_bias_AP_isotropic of Modelling_TwoPointCorrelation.h!"); }

	/**
	 * @brief fit the Alcock-Paczynski effect for the two point
	 * correlation function. 
	 * \f$\xi(s)= B^2  \xi_{DM}(\alpha s)\ + A_0 + A_1/s +A_2/s^2\f$.
	 * &xi<SUB>DM</SUB> is computed at the fiducial cosmology.
	 * B, A<SUB>0</SUB>, A<SUB>1</SUB>, A<SUB>2</SUB> are nuisance
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
	 * either _free_ or _fixed_
	 *
	 * @param pT_B the parameter type of B: it can be either _free_
	 * or _fixed_
	 *
	 * @param pT_A0 the parameter type of A0: it can be either
	 * _free_ or _fixed_
	 *
	 * @param pT_A1 the parameter type of A1: it can be either
	 * _free_ or _fixed_
	 *
	 * @param pT_A2 the parameter type of A2: it can be either
	 * _free_ or _fixed_
	 *
	 * @return none
	 */
	virtual void set_model_AP_isotropic (const statistics::Prior alpha_prior, const statistics::Prior B_prior, const statistics::Prior A0_prior, const statistics::Prior A1_prior, const statistics::Prior A2_prior, const statistics::ParameterType pT_alpha=statistics::_free_, const statistics::ParameterType pT_B=statistics::_free_, const statistics::ParameterType pT_A0=statistics::_free_, const statistics::ParameterType pT_A1=statistics::_free_, const statistics::ParameterType pT_A2=statistics::_free_)
	{ ErrorMsg("Error in fit_AP_isotropic of Modelling_TwoPointCorrelation.h!"); }


    };
  }
}

#endif
