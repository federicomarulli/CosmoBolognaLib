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
 *  @file Headers/Modelling_NumberCounts.h
 *
 *  @brief The class Modelling_NumberCounts
 *
 *  This file defines the interface of the class Modelling, used to
 *  model number counts of any kind
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGNC__
#define __MODELLINGNC__


#include "NumberCounts.h"
#include "Modelling.h"
#include "ModelFunction_NumberCounts.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    /**
     *  @brief The namespace of the <B> number counts
     *  modelling </B>
     *  
     *  The \e modelling::numbercounts namespace contains all the
     *  functions and classes to model number counts
     */
    namespace numbercounts {
    
      /**
       *  @class Modelling_NumberCounts
       *  Modelling_NumberCounts.h
       *  "Headers/Modelling_NumberCounts.h"
       *
       *  @brief The class Modelling_NumberCounts
       *
       *  This file defines the interface of the base class
       *  Modelling_NumberCounts, used for modelling any kind of
       *  number counts measurements
       */
      class Modelling_NumberCounts : public Modelling {
      
      protected:

	/// the histogram type
	glob::HistogramType m_HistogramType;

	/// the normalization factor
	double m_fact;

	/// the container of parameters for number counts model computation
	modelling::numbercounts::STR_NC_data_model m_data_model;

	/// the container of parameters for size number counts model computation
	modelling::numbercounts::STR_NCSF_data_model m_data_model_SF;

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  _NumberCounts
	 */
	Modelling_NumberCounts () = default;
	
	/**
	 *  @brief constuctor
	 *  @param nc the number counts to model
	 *  _NumberCounts
	 */
	Modelling_NumberCounts (const std::shared_ptr<cbl::measure::numbercounts::NumberCounts> nc)
	{ m_HistogramType = nc->HistogramType(); m_fact = nc->fact(); }
	
	/**
	 *  @brief constuctor
	 *  @param hist_type the histogram type
	 *  @param fact the normalization factor
	 *  _NumberCounts
	 */
	Modelling_NumberCounts (glob::HistogramType hist_type, double fact)
	{ m_HistogramType = hist_type; m_fact = fact;}
	
	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~Modelling_NumberCounts () = default;

	///@}

	/**
	 *  @name Member functions used to get the protected members of the class
	 */
	///@{

	/**
	 * @brief get the member \e m_data_model
	 * @return the container of parameters for number counts
         * model computation
	 */
	modelling::numbercounts::STR_NC_data_model data_model () { return m_data_model; }


	///@{

	/**
	 * @brief get the member \e m_data_model_SF
	 * @return the container of parameters for size number counts
         * model computation
	 */
	modelling::numbercounts::STR_NCSF_data_model data_model_SF () { return m_data_model_SF; }
	
	///@}

	
	/**
	 *  @name Member functions used to set the model parameters
	 */
	///@{

	/**
	 *  @brief set the data used to construct generic models of
	 *  number counts
         *  
	 *  @param cosmology the cosmological model used to compute
	 *  &xi;<SUB>DM</SUB>
	 *
	 *  @param redshift redshift
	 *
	 *  @param method_Pk method used to compute the power spectrum
	 *  (i.e. the Boltzmann solver); valid choices for method_Pk
	 *  are: CAMB [http://camb.info/], CLASS
	 *  [http://class-code.net/], MPTbreeze-v1
	 *  [http://arxiv.org/abs/1207.1465], EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *    
	 *  @param k_min minimum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *  
	 *  @param k_max maximum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *
	 *  @param step number of steps used to compute the binned
	 *  power spectrum
	 *
	 *  @param store_output if true the output files created by
	 *  the Boltzmann solver are stored; if false the output files
	 *  are removed
	 *
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *
	 *  @param Delta \f$\Delta\f$: the overdensity, defined as the
	 *  mean interior density relative to the background
	 *
	 *  @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
	 *  virial overdensity
	 *
	 *  @param model_MF author(s) who proposed the mass function
         *
	 *  @param selection_function_file input file with the
	 *  selection function
	 *
	 *  @param selection_function_column vector containing the
	 *  columns with {mass, redshift, selection function}
	 *
	 *  @param z_min minimum redshift
	 * 
	 *  @param z_max maximum redshift
	 *
	 *  @param z_step the number of bins
	 *  for redshift vector
	 *  
	 *  @param Mass_min minimum halo mass
	 *
	 *  @param Mass_max maximum halo mass
	 *
	 *  @param Mass_step the number of bins for mass vector
	 *
	 *  @param area_degrees the area in degrees
	 *
	 *  @param prec the precision
	 *
	 *  
	 */
	void set_data_model (const cbl::cosmology::Cosmology cosmology={}, const double redshift=0., const std::string method_Pk="CAMB", const double k_min=1.e-4, const double k_max=100., const int step=500, const bool store_output=true, const int norm=-1, const double Delta=200., const bool isDelta_vir=true, const std::string model_MF="Tinker", const std::string selection_function_file=par::defaultString, const std::vector<int> selection_function_column={}, const double z_min=par::defaultDouble, const double z_max=par::defaultDouble, const int z_step=50, const double Mass_min=par::defaultDouble, const double Mass_max=par::defaultDouble, const int Mass_step=100, const double area_degrees=par::defaultDouble, const double prec=1.e-4);
	
	/**
	 *  @brief set the data used to construct a model of
	 *  number counts as a function of a mass proxy, here 
	 *  expressed as \f$\lambda\f$, with the following
	 *  functional form:
	 *  
	 *  \f$ \langle N(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle 
	 *  = w(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\,\,\Omega 
	 *  \int_{0}^{\infty} {\rm d} z_{\rm tr}\,\,
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{0}^{\infty}{\rm d}\lambda_{\rm tr}\,\,
	 *  P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\,
	 *  \int_{\Delta z_{\text{ob},j}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  \int_{\Delta\lambda_{\text{ob},i}}{\rm d} \lambda_{\rm ob} 
	 *  \,\,P(\lambda_{\rm ob}|\lambda_{\rm tr}), \f$
	 *
	 *  where \f$ w(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j}) \f$ is
	 *  the weight derived from the selection function (see e.g. Lesci et al. 2021),
	 *  and \f$\Omega\f$ is the survey effective area.
	 *
	 *  The distirbution \f$P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\f$ is expressed as:
	 *
	 *  \f$P(\lambda_{\rm tr}|M_{\rm tr},z_{\rm tr})= 
	 *  P(M_{\rm tr}|\lambda_{\rm tr},z_{\rm tr})\,
	 *  P(\lambda_{\rm tr}|z_{\rm tr})\,/\,P( M_{\rm tr}|z_{\rm tr}),\f$
	 *
	 *  where \f$P(M_{\rm tr}|\lambda_{\rm tr},z_{\rm tr})\f$ is
	 *  a log-normal whose mean is given by the mass-mass proxy
	 *  relation, i.e. 
	 *
	 *  \f$\log M = \alpha + \beta 
	 *  \log (\lambda/\lambda_{\rm piv}) + \gamma \log (f(z)),\f$
	 *
	 *  and whose mean is given by the intrinsic scatter, \f$ \sigma_{\rm intr} \f$, expressed as:
	 *
	 *  \f$ \sigma_{\rm intr} = \sigma_0 + \sigma_{\lambda} 
	 *  \log (\lambda/\lambda_{\rm piv})^{e_{\lambda}} + \sigma_z \log (f(z))^{e_z}.\f$
	 *
	 *  The distribution \f$ P(\lambda_{\rm tr}|z_{\rm tr}) \f$ is derived
	 *  from mock catalogues, and it has the functional form of a power-law
	 *  with an exponential cut-off, i.e.:
	 *  \f$ P(\lambda_{\rm tr}|z_{\rm tr}) = a \, \lambda_{\rm tr}^{-b} \, e^{-c\lambda_{\rm tr}} \f$.
	 *
	 *  Furthermore, \f$P( M_{\rm tr}|z_{\rm tr}) = \int_{0}^{\infty}
	 *  {\rm d} \lambda_{\rm tr}\, 
	 *  P(M_{\rm tr}|\lambda_{\rm tr},z_{\rm tr})\,P(\lambda_{\rm tr}|z_{\rm tr}).\f$
	 *
	 *  Finally, \f$P(z_{\rm ob}|z_{\rm tr})\f$ and \f$P(\lambda_{\rm ob}|\lambda_{\rm tr})\f$
	 *  are Gaussian distributions.
	 *  
	 *  @param cosmology the cosmological model
	 *
	 *  @param cluster object of the class Cluster
	 *
	 *  @param SF_weights the counts weights derived from the selection
	 *  function (see e.g. Lesci et al. 2021)
	 *
	 *  @param z_error the absolute error on the redshifts, defining
	 *  the width of \f$ P(z_{\rm ob}|z_{\rm tr}) \f$
	 *
	 *  @param proxy_relative_error the relative error on the mass proxies, defining
	 *  the width of \f$ P(\lambda_{\rm ob}|\lambda_{\rm tr}) \f$
	 *
	 *  @param z_pivot the redshift pivot in the scaling relation
	 *
	 *  @param proxy_pivot the mass proxy pivot in the scaling relation
	 *
	 *  @param mass_pivot the mass pivot; for example, if the scaling relation
	 *  is written as \f$ \log [M/(10^{14}M_\odot /h)] = 
	 *  \alpha+\beta \log(\lambda/\lambda_{piv})+\gamma\log f(z) \f$,
	 *  then mass_pivot = 1.e14
	 *
	 *  @param scalrel_z_evo functional form of the redshift evolution
	 *  function in the scaling relation: "E_z" \f$\rightarrow\f$ 
	 *  \f$ f(z)=E(z)/E(z_{piv}) \f$, "direct" \f$\rightarrow\f$ 
	 *  \f$ f(z)=(1+z)/(1+z_{piv}) \f$
	 *
	 *  @param log_base the base of the logarithm used in the 
	 *  mass-mass proxy scaling relation
	 *
	 *  @param method_Pk method used to compute the power spectrum
	 *  (i.e. the Boltzmann solver); valid choices for method_Pk
	 *  are: CAMB [http://camb.info/], CLASS
	 *  [http://class-code.net/], MPTbreeze-v1
	 *  [http://arxiv.org/abs/1207.1465], EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *
	 *  @param store_output if true the output files created by
	 *  the Boltzmann solver are stored; if false the output files
	 *  are removed
	 *
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *
	 *  @param Delta \f$\Delta\f$: the overdensity, defined as the
	 *  mean interior density relative to the background
	 *
	 *  @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
	 *  virial overdensity
	 *
	 *  @param model_MF author(s) who proposed the mass function
	 *
	 *  @param z_min minimum redshift
	 * 
	 *  @param z_max maximum redshift
	 *
	 *  @param area_degrees the area in degrees
	 *
	 *  @param prec the precision
	 *
	 *  
	 */
	void set_data_model (const cbl::cosmology::Cosmology cosmology, const cbl::catalogue::Cluster cluster, const std::vector<double> SF_weights, const double z_error, const double proxy_relative_error, const double z_pivot, const double proxy_pivot, const double mass_pivot, const std::string scalrel_z_evo, const double log_base, const std::string method_Pk, const bool store_output=true, const int norm=-1, const double Delta=200., const bool isDelta_vir=true, const std::string model_MF="Tinker", const double z_min=par::defaultDouble, const double z_max=par::defaultDouble, const double area_degrees=par::defaultDouble, const double prec=1.e-4);

	///@}
		
	/**
	 *  @brief Member functions used to set the model parameters
	 *
	 *  @param cosmology the cosmological model used to compute
	 *  the void size function
	 *
	 *  @param radii the void radii
	 *
	 *  @param redshift the redshift
	 *
	 *  @param model size function model name; valid choices for
	 *  model name are SvdW (Sheth and van de Weygaert, 2004),
	 *  linear and Vdn (Jennings et al., 2013)
	 *
	 *  @param b_eff the effective bias of the sample
	 *
	 *  @param slope first coefficent to convert the effective bias
	 *  (default value set to \f$0.854\f$)
	 *
	 *  @param offset second coefficent to convert the effective
	 *  bias (default value set to \f$0.420\f$)
	 *
	 *  @param deltav_NL the non linear density contrast:
	 *  \f$\rho_v/\rho_m\f$ (default value set to \f$-0.795\f$)
	 *
	 *  @param del_c critical value of the linear density field
	 *  (default value set to \f$1.06\f$)
	 *
	 *  @param method_Pk method used to compute the power spectrum
	 *  (i.e. the Boltzmann solver); valid choices for method_Pk
	 *  are: CAMB [http://camb.info/], CLASS
	 *  [http://class-code.net/], MPTbreeze-v1
	 *  [http://arxiv.org/abs/1207.1465], EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *
	 *  @param k_Pk_ratio wave vector module required to compute
	 *  the growth factor (cbl::cosmology::Cosmology::DD_norm())
	 *  with the method "Pk_ratio". This parameter represents
	 *  wavenumber at which the ratio between the the power
	 *  spectra at different redshift is computed. It is
	 *  recommended to use this method when dealing with
	 *  cosmologies alternative to the LCDM. To avoid to compute
	 *  the growth factor with this method, set k_Pk_ratio=-1.
	 *
	 *  @param store_output if true the output files created by
	 *  the Boltzmann solver are stored; if false the output files
	 *  are removed
	 *
	 *  @param output_root output_root of the parameter file used
	 *  to compute the power spectrum and &sigma;(mass); it can be
	 *  any name
	 *
	 *  @param interpType method to interpolate the power spectrum
	 *
	 *  @param k_max maximum wave vector module up to which the power
	 *  spectrum is computed
	 *           
	 *  @param input_file either the parameter file or the power
	 *  spectrum file; if a parameter file is provided,
	 *  i.e. input_file!=NULL and is_parameter_file=true, it will
	 *  be used to compute the power spectrum; if a power spectrum
	 *  file is provided, i.e. input_file!=NULL and
	 *  is_parameter_file=false, then the provided power spectrum
	 *  will be used directly; in both cases
	 *  &sigma;<SUP>2</SUP>(M) is computed by integrating the
	 *  computed/provided power spectrum ignoring the cosmological
	 *  parameters of the object
	 *
	 *  @param is_parameter_file true \f$\rightarrow\f$ the
	 *  input_file is a parameter file, used to compute the power
	 *  spectrum with the method specified by method_Pk; false
	 *  \f$\rightarrow\f$ the input_file is a file containing the
	 *  power spectrum
	 *
	 *  
	 */
	
	void set_data_model_SF (const cosmology::Cosmology cosmology, const std::vector<double> radii, const double redshift, const std::string model, const double b_eff, double slope=0.854, double offset=0.420, const double deltav_NL=-0.795, const double del_c=1.69, const std::string method_Pk="EisensteinHu", const double k_Pk_ratio=-1., const bool store_output=true, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true); 


	/**
	 *  @brief set the data used to construct mass
	 *  number counts of simulation snapshots
         *  
	 *  @param cosmology the cosmological model used to compute
	 *  &xi;<SUB>DM</SUB>
	 *
	 *  @param redshift redshift
	 *
	 *  @param method_Pk method used to compute the power spectrum
	 *  (i.e. the Boltzmann solver); valid choices for method_Pk
	 *  are: CAMB [http://camb.info/], CLASS
	 *  [http://class-code.net/], MPTbreeze-v1
	 *  [http://arxiv.org/abs/1207.1465], EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *    
	 *  @param k_min minimum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *  
	 *  @param k_max maximum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *
	 *  @param step number of steps used to compute the binned
	 *  power spectrum
	 *
	 *  @param store_output if true the output files created by
	 *  the Boltzmann solver are stored; if false the output files
	 *  are removed
	 *
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *
	 *  @param Delta \f$\Delta\f$: the overdensity, defined as the
	 *  mean interior density relative to the background
	 *
	 *  @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
	 *  virial overdensity
	 *
	 *  @param model_MF author(s) who proposed the mass function
	 *
	 *  @param Volume the snapshot volume
	 *  
	 *  @param Mass_min minimum halo mass
	 *
	 *  @param Mass_max maximum halo mass
	 *
	 *  @param Mass_step the number of bins for mass vector
       	 *
	 *  @param prec the precision
	 *
	 *  @return none, or an error message if the derived object
	 *  does not have this member
	 */
	virtual void set_data_model_snapshot (const cbl::cosmology::Cosmology cosmology={}, const double redshift=0., const std::string method_Pk="CAMB", const double k_min=1.e-4, const double k_max=100., const int step=500, const bool store_output=true, const int norm=-1, const double Delta=200., const bool isDelta_vir=true, const std::string model_MF="Tinker", const double Volume=par::defaultDouble, const double Mass_min=par::defaultDouble, const double Mass_max=par::defaultDouble, const int Mass_step=100, const double prec=1.e-4) 
	{ 
	  (void)cosmology; (void)redshift; (void)method_Pk; (void)k_min; (void)k_max; (void)step;
	  (void)store_output; (void)norm; (void)Delta; (void)isDelta_vir; (void)model_MF;
	  (void)Volume; (void)prec; (void)Mass_min; (void)Mass_max; (void)Mass_step;
	  cbl::ErrorCBL("", "set_data_model_snapshot", "Modelling_NumberCounts.h");
	}


	///@}
	
      };
    }
  }
}

#endif
