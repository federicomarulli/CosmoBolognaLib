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

#ifndef __MODELLINGTWOPOINT1D__
#define __MODELLINGTWOPOINT1D__


#include "Modelling_TwoPointCorrelation.h"
#include "ModelFunction_TwoPointCorrelation1D.h"


// ===================================================================================================


namespace cosmobl {
  
  namespace modelling {
    
    namespace twopt {
      
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

      protected:

	/// the container of parameters for the HOD modelling of the two-point correlation function
	modelling::twopt::STR_data_HOD m_data_HOD;
	
	/**
	 *  @name free parameters to model the two-point correlation function
	 */
	///@{

	/// \f$\sigma_8\f$
	shared_ptr<statistics::Parameter> m_sigma8;

	/// the linear bias, \f$b\f$
	shared_ptr<statistics::Parameter> m_bias;
	
	/// \f$b\sigma_8\f$
	shared_ptr<statistics::Parameter> m_bsigma8;

	/// \f$f\sigma_8\f$
	shared_ptr<statistics::Parameter> m_fsigma8;

	/// polynomial 
	vector<shared_ptr<statistics::Parameter>> m_polynomial;

	/// position of the BAO peak
	shared_ptr<statistics::Parameter> m_peak;

	/// position of the BAO dip
	shared_ptr<statistics::Parameter> m_dip;

	/// position of the BAO linear point
	shared_ptr<statistics::Parameter> m_linear_point;

        /// damping of the BAO peak
	shared_ptr<statistics::Parameter> m_sigmaNL;

	///@}

	
	/**
	 *  @name HOD free parameters
	 */
	///@{
	
	/// \f$M_{min}\f$: the mass scale at which 50% of haloes host a satellite galaxy
	shared_ptr<statistics::Parameter> m_Mmin;

	/// \f$\sigma_{\log M_h}\f$: transition width reflecting the scatter in the luminosity-halo mass relation
	shared_ptr<statistics::Parameter> m_sigmalgM;

	/// \f$M_0\f$: the cutoff mass 
	shared_ptr<statistics::Parameter> m_M0;

	/// \f$M_1\f$: the amplitude of the power law
	shared_ptr<statistics::Parameter> m_M1;
	
	/// \f$\alpha\f$: the slope of the power law	  
	shared_ptr<statistics::Parameter> m_alpha;
	
	///@}

	
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
	 *  @return object of type Modelling_TwoPointCorrelation1D
	 */
	Modelling_TwoPointCorrelation1D (const shared_ptr<cosmobl::measure::twopt::TwoPointCorrelation> twop)
	  : Modelling_TwoPointCorrelation(twop) {}
	
	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_TwoPointCorrelation1D () = default;
	
	///@}


	/**
	 *  @name Member functions used to get the protected members of the class
	 */
	///@{
	
	/**
	 * @brief get the member \e m_data_HOD	 
	 * @return the container of parameters for the HOD modelling
	 * of the two-point correlation function
	 */
	modelling::twopt::STR_data_HOD data_HOD () { return m_data_HOD; }
	
	///@}
	
	
	/**
	 *  @name Member functions used to get the model parameters
	 */
	///@{

	/**
	 *  @brief get the private member m_sigma8
	 *
	 *  @return the parameter \f$\sigma_8\f$
	 */
	shared_ptr<statistics::Parameter> sigma8 () { return m_sigma8; }
	
	/**
	 *  @brief get the private member m_bias
	 *
	 *  @return the linear bias, \f$b\f$
	 */
	shared_ptr<statistics::Parameter> bias () { return m_bias; }
	
	/**
	 *  @brief get the private member m_bsigma8
	 *
	 *  @return the parameter \f$b\sigma_8\f$
	 */
	shared_ptr<statistics::Parameter> bsigma8 () { return m_bsigma8; }

	/**
	 *  @brief get the private member m_fsigma8
	 *
	 *  @return the parameter \f$f\sigma_8\f$
	 */  
	shared_ptr<statistics::Parameter> fsigma8 () { return m_fsigma8; }
	  
	/**
	 *  @brief get the private member m_polynomial
	 *
	 *  @return m_polynomial
	 */ 
	vector<shared_ptr<statistics::Parameter>> polynomial () { return m_polynomial; }
	  
	/**
	 *  @brief get the private member m_polynomial
	 *
	 *  @param i the i-th index
	 *  @return the i-th element of m_polynomial
	 */ 
	shared_ptr<statistics::Parameter> polynomial (const int i) {return m_polynomial[i];}
	  
	/**
	 *  @brief get the private member m_peak
	 * 
	 *  @return the BAO peak position
	 */
	shared_ptr<statistics::Parameter> peak () { return m_peak; }
	  	  
	/**
	 *  @brief get the private member m_dip
	 * 
	 *  @return the BAO dip position
	 */
	shared_ptr<statistics::Parameter> dip () { return m_dip; }
	  	  
	/**
	 *  @brief get the private member m_linear_point
	 * 
	 *  @return the BAO linear_point
	 */
	shared_ptr<statistics::Parameter> linear_point () { return m_linear_point; }
	  	  
	/**
	 *  @brief get the private member m_sigmaNL
	 * 
	 *  @return the BAO damping \f$\Sigma_{NL}\f$
	 */
	shared_ptr<statistics::Parameter> sigmaNL () { return m_sigmaNL; }

	/**
	 *  @brief get the private member m_Mmin
	 * 
	 *  @return the parameter \f$M_{min}\f$
	 */
	shared_ptr<statistics::Parameter> Mmin () { return m_Mmin; }

	/**
	 *  @brief get the private member m_sigmalgM
	 * 
	 *  @return the parameter \f$sigma_{\log M_h}\f$
	 */
	shared_ptr<statistics::Parameter> sigmalgM () { return m_sigmalgM; }

	/**
	 *  @brief get the private member m_M0
	 * 
	 *  @return the parameter \f$M_0\f$
	 */
	shared_ptr<statistics::Parameter> M0 () { return m_M0; }

	/**
	 *  @brief get the private member m_M1
	 * 
	 *  @return the parameter \f$M_1\f$
	 */
	shared_ptr<statistics::Parameter> M1 () { return m_M1; }

	/**
	 *  @brief get the private member m_alpha
	 * 
	 *  @return the parameter \f$\alpha\f$
	 */
	shared_ptr<statistics::Parameter> alpha () { return m_alpha; }
	
	///@}
	

	/**
	 *  @name Member functions used to write the outputs
	 */
	///@{
      
	/**
	 *  @brief compute and write the model
	 *
	 *  @param dir the output directory
	 *
	 *  @param file the output file
	 *
	 *  @param xx vector of points at which the model is computed
	 *
	 *  @param parameter vector containing the input parameters
	 *  used to compute the model; if this vector is not provided,
	 *  the model will be computed using the best-fit parameters
	 *
	 *  @param start the starting position for each chain
	 *
	 *  @param thin the position step
	 *
	 *  @return none
	 */
        void write_model (const string dir, const string file, const vector<double> xx={}, const vector<double> parameter={}, const int start=0, const int thin=10) override;

        /**
         *  @brief compute the best-fit model from MCMC chains
         *
         *  @param xx vector of points at which the model is computed
	 *  @param start the starting position for each chain
	 *  @param thin the position step
         *  @param median_model the median model
         *  @param low_model the model at 16th percentile
         *  @param up_model the model at 84th percentile
         *  @return none
         */
        void compute_model_from_chains (const vector<double> xx, const int start, const int thin, vector<double> &median_model, vector<double> &low_model, vector<double> &up_model) override;

	///@}

	
	/**
	 *  @name Member functions used to set the model parameters
	 */
	///@{

	/**
	 *  @brief set the scale range used for the fit
	 *
	 *  @param xmin the minimum x value
	 *  @param xmax the maximum x value
	 *  @return none
	 */
	void set_fit_range (const double xmin, const double xmax)
        { m_data_fit = m_data->cut(xmin, xmax); m_fit_range = true; }
	
	/**
	 *  @brief set the scale range used for the fit
	 *
	 *  @param xmin the minimum x value
	 *  @param xmax the maximum x value
	 *  @param ndataset the number of dataset
	 *  @return none
	 */
	void set_fit_range (const double xmin, const double xmax, const int ndataset)
	{(void)xmin; (void)xmax; (void)ndataset; ErrorCBL("Error in set_fit_range of Modelling_TwoPointCorrelation1D!");}

	/**
	 *  @brief set the data used to construct generic models of
	 *  the two-point correlation function
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
	 *  @param sigmaNL_perp damping of the wiggles in the linear
	 *  power spectrum perpendicular to the line of sight
	 *
	 *  @param sigmaNL_par damping of the wiggles in the linear
	 *  power spectrum parallel to the line of sight
	 *
	 *  @param NL 0 &rarr; linear power spectrum; 1 &rarr;
	 *  non-linear power spectrum
	 *
	 *  @param bias the linear bias
	 *
	 *  @param pimax the upper limit of the line of sight
	 *  integration
	 *
	 *  @param r_min minimum separation up to which the binned
	 *  dark matter correlation function is computed
	 *
	 *  @param r_max maximum separation up to which the binned
	 *  dark matter correlation function is computed
	 *    
	 *  @param k_min minimum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *  
	 *  @param k_max maximum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *
	 *  @param step number of steps used to compute the binned
	 *  dark matter correlation function
	 *
	 *  @param output_dir the output_dir directory
	 *  where the output of external codes are written
	 *
	 *  @param output_root output_root of the parameter file used
	 *  to compute the power spectrum and &sigma;(mass); it can be
	 *  any name
	 *  
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
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
	 *  @param Delta \f$\Delta\f$: the overdensity, defined as the
	 *  mean interior density relative to the background
	 *
	 *  @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
	 *  virial overdensity
	 *
	 *  @param cluster_redshift vector containing the cluster
	 *  redshifts
	 *
	 *  @param cluster_mass_proxy vector containing the cluster
	 *  mass proxies
	 *
	 *  @param cluster_mass_proxy_error vector containing the
	 *  errors on the cluster mass proxies
	 *
	 *  @param model_bias author(s) who proposed the bias; valid
	 *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
	 *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
	 *  correction of Warren 2004), Tinker (Tinker et al. 2010)
	 *
	 *  @param meanType meanType="mean_bias" \f$\rightarrow\f$ the
	 *  effective bias is the mean bias; meanType="mean_pair_bias"
	 *  \f$\rightarrow\f$ the effective bias is pair mean bias
	 *
	 *  @param seed seed for scatter in scaling relation
	 *
	 *  @return none
	 */
	void set_data_model (const cosmology::Cosmology cosmology={}, const double redshift=0., const string method_Pk="CAMB", const double sigmaNL_perp=0., const double sigmaNL_par=0., const bool NL=true, const double bias=1., const double pimax=40., const double r_min=1., const double r_max=350., const double k_min=1.e-4, const double k_max=100., const int step=500,  const string output_dir=par::defaultString, const string output_root="test", const int norm=-1, const double aa=0., const bool GSL=true, const double prec=1.e-3, const string file_par=par::defaultString, const double Delta=200., const bool isDelta_vir=true, const vector<double> cluster_redshift={}, const vector<double> cluster_mass_proxy={}, const vector<double> cluster_mass_proxy_error={}, const string model_bias="Tinker", const string meanType="mean_bias", const int seed=744524123);
	
	/**
	 *  @brief set the data used to construct the HOD model
	 *
	 *  @param cosmology the cosmological model used to compute
	 *  &xi;<SUB>DM</SUB>
	 *
	 *  @param redshift redshift
	 *
	 *  @param model_MF author(s) who proposed the mass function;
	 *  valid authors are: PS (Press & Schechter), ST (Sheth &
	 *  Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
	 *  al. 2006), Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH
	 *  (halo MF by Shen et al. 2006), ShenF (filament MF by Shen et
	 *  al. 2006), ShenS (sheet MF by Shen et al. 2006), Tinker
	 *  (Tinker et al. 2008), Crocce (Crocce et al. 2010),
	 *  Angulo_FOF (FoF MF by Angulo et al. 2012), Angulo_Sub
	 *  (SUBFIND MF by Angulo et al. 2012), Watson_FOF (FoF MF by
	 *  Watson et al. 2012), Watson_SOH (Spherical Overdensity halo
	 *  MF by Watson et al. 2012), Manera (Manera et al. 2010),
	 *  Bhattacharya (Bhattacharya et al. 2011), Courtin (Courtin et
	 *  al. 2010), Peacock (by Peacock at al. 2007)
	 *	 
	 *  @param model_bias author(s) who proposed the bias; valid
	 *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
	 *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
	 *  correction of Warren 2004), Tinker (Tinker et al. 2010)
	 *
	 *  @param Mh_min minimum halo mass
	 *
	 *  @param Mh_max maximum halo mass
	 *
	 *  @param pi_max the upper limit of the line of sight
	 *  integration
	 *
	 *  @param r_max_int the maximum separation used to count
	 *  pairs; it is used to compute the upper limit of the
	 *  line-of-sight integration in approximate projected
	 *  clustering estimators
	 *
	 *  @param r_min minimum separation up to which the binned
	 *  dark matter correlation function is computed
	 *
	 *  @param r_max maximum separation up to which the binned
	 *  dark matter correlation function is computed
	 *  
	 *  @param k_min minimum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *  
	 *  @param k_max maximum wave vector module up to which the
	 *  binned dark matter power spectrum is computed

	 *  @param step number of steps used to compute the binned
	 *  dark matter correlation function and power spectrum
	 *
	 *  @param method_Pk method used to compute the power
	 *  spectrum; valid choices for method_Pk are: CAMB
	 *  [http://camb.info/], classgal_v1 [http://class-code.net/],
	 *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465],
	 *  EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *  
	 *  @param NL false \f$rightarrow\f$ linear power spectrum;
	 *  true \f$rightarrow\f$ non-linear power spectrum
	 * 
	 *  @param output_root output_root of the parameter file used
	 *  to compute the power spectrum and &sigma;(mass); it can be
	 *  any name
	 *
	 *  @param Delta \f$\Delta\f$ the overdensity, defined as the
	 *  mean interior density relative to the background
	 *
	 *  @param kk wave vector module
	 *
	 *  @param interpType method to interpolate the power spectrum
	 *  
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *  
	 *  @param prec accuracy of the GSL integration
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
	 *  @param is_parameter_file true \f$rightarrow\f$ the
	 *  input_file is a parameter file, used to compute the power
	 *  spectrum with the method specified by method_Pk; false
	 *  \f$rightarrow\f$ the input_file is a file containing the
	 *  power spectrum
	 *
	 *  @param model_cM the author(s) of the concentration-mass
	 *  relation (see cosmobl::modelling::twopt::concentration)
	 *  
	 *  @param profile the density profile (see
	 *  cosmobl::modelling::twopt::concentration)
	 *
	 *  @param halo_def the halo definition (see
	 *  cosmobl::modelling::twopt::concentration)
	 *
	 *  @return none
	 */
	void set_data_HOD (const cosmology::Cosmology cosmology={}, const double redshift=0., const string model_MF="Tinker", const string model_bias="Tinker", const double Mh_min=0., const double Mh_max=1.e16, const double pi_max=100., const double r_max_int=100., const double r_min=1.e-3, const double r_max=350., const double k_min=0., const double k_max=100., const int step=200, const string method_Pk="CAMB", const bool NL=true, const string output_root="test", const double Delta=200., const double kk=0., const string interpType="Linear", const int norm=-1, const double prec=1.e-2, const string input_file=par::defaultString, const bool is_parameter_file=true, const string model_cM="Duffy", const string profile="NFW", const string halo_def="vir");
	
	/**
	 *  @brief set the data used to construct the model for
	 *  clustering analyses using the selection function
	 *
	 *  @param cosmology the cosmological model used to compute
	 *  &xi;<SUB>DM</SUB>
	 *	
	 *  @param test_cosmology the cosmological model used to
	 *  measure &xi;<SUB>DM</SUB>
	 *
	 *  @param mean_redshift mean redshift of the sample
	 *
	 *  @param model_MF author(s) who proposed the mass function;
	 *  valid authors are: PS (Press & Schechter), ST (Sheth &
	 *  Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
	 *  al. 2006), Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH
	 *  (halo MF by Shen et al. 2006), ShenF (filament MF by Shen et
	 *  al. 2006), ShenS (sheet MF by Shen et al. 2006), Tinker
	 *  (Tinker et al. 2008), Crocce (Crocce et al. 2010),
	 *  Angulo_FOF (FoF MF by Angulo et al. 2012), Angulo_Sub
	 *  (SUBFIND MF by Angulo et al. 2012), Watson_FOF (FoF MF by
	 *  Watson et al. 2012), Watson_SOH (Spherical Overdensity halo
	 *  MF by Watson et al. 2012), Manera (Manera et al. 2010),
	 *  Bhattacharya (Bhattacharya et al. 2011), Courtin (Courtin et
	 *  al. 2010), Peacock (by Peacock at al. 2007)
	 *	 
	 *  @param model_bias author(s) who proposed the bias; valid
	 *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
	 *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
	 *  correction of Warren 2004), Tinker (Tinker et al. 2010)
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
	 *  @param Mass_min minimum halo mass
	 *
	 *  @param Mass_max maximum halo mass
	 *
	 *  @param Delta \f$\Delta\f$ the overdensity
	 *
	 *  @param isDelta_vir true \f$\rightarray\f$ \f$\Delta\f$ is
	 *  the virial overdensity
	 *
	 *  @param method_Pk method used to compute the power
	 *  spectrum; valid choices for method_Pk are: CAMB
	 *  [http://camb.info/], classgal_v1 [http://class-code.net/],
	 *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465],
	 *  EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *  
	 *  @param output_dir
	 *  
	 *  @param prec accuracy of the GSL integration
	 *  
	 *  @param k_min minimum wave vector module up to which the
	 *  power spectrum is computed
	 *
	 *  @param k_max maximum wave vector module up to which the
	 *  power spectrum is computed
	 *
	 *  @param prec accuracy of the GSL integration 
	 *
	 *  @param step
	 * 
	 *  @param mass_step
	 *
	 *  @return none
	 */
	void set_data_model_cluster_selection_function (const cosmology::Cosmology cosmology, const cosmology::Cosmology test_cosmology, const double mean_redshift, const string model_MF, const string model_bias, const string selection_function_file, const vector<int> selection_function_column={}, const double z_min=par::defaultDouble, const double z_max=par::defaultDouble, const double Mass_min=par::defaultDouble, const double Mass_max=par::defaultDouble, const double Delta=200, const bool isDelta_vir=false, const string method_Pk="CAMB", const string output_dir=par::defaultString, const double k_min=1.e-4, const double k_max=100, const double prec=1.e-2, const int step=200, const int mass_step=50);

	///@}
	
      };
    }
  }
}

#endif
