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
 *  @file Headers/Modelling_TwoPointCorrelation1D.h
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


#include "Modelling1D.h"
#include "Modelling_TwoPointCorrelation.h"
#include "ModelFunction_TwoPointCorrelation1D.h"


// ===================================================================================================


namespace cbl {
  
  namespace modelling {
    
    namespace twopt {
      
      /**
       *  @class Modelling_TwoPointCorrelation1D
       *  Modelling_TwoPointCorrelation1D.h
       *  "Headers/Modelling_TwoPointCorrelation1D.h"
       *
       *  @brief The class Modelling_TwoPointCorrelation1D
       *
       *  This file defines the interface of the base class
       *  Modelling_TwoPointCorrelation1D, that contains all the common
       *  methods to model 1D two-point correlation functions
       *
       */
      class Modelling_TwoPointCorrelation1D : public Modelling1D, public Modelling_TwoPointCorrelation 
      {

      protected:

	/// the container of parameters for the HOD modelling of the two-point correlation function
	modelling::twopt::STR_data_HOD m_data_HOD;
	
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
	Modelling_TwoPointCorrelation1D (const std::shared_ptr<cbl::measure::twopt::TwoPointCorrelation> twop);

	/**
	 *  @brief constructor
	 *  
	 *  @param dataset the two-point correlation dataset
	 *  
	 *  @param twoPType the two-point correlation type
	 *
	 *  @return object of type Modelling_TwoPointCorrelation1D
	 */
	Modelling_TwoPointCorrelation1D (const std::shared_ptr<cbl::data::Data> dataset, const measure::twopt::TwoPType twoPType);
	
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
	 *  @param cosmology_mass cosmology used to measure the cluster
	 *  masses
	 *
	 *  @param redshift_source vector containing the redshifts of the
	 *  source galaxies, in case the cluster masses are estimated from
	 *  weak lensing
	 *
	 *  @return none
	 */
	void set_data_model (const cbl::cosmology::Cosmology cosmology={}, const double redshift=0., const std::string method_Pk="CAMB", const double sigmaNL_perp=0., const double sigmaNL_par=0., const bool NL=true, const double bias=1., const double pimax=40., const double r_min=1., const double r_max=350., const double k_min=1.e-4, const double k_max=100., const int step=500,  const std::string output_dir=par::defaultString, const std::string output_root="test", const int norm=-1, const double aa=0., const bool GSL=true, const double prec=1.e-3, const std::string file_par=par::defaultString, const double Delta=200., const bool isDelta_vir=true, const std::vector<double> cluster_redshift={}, const std::vector<double> cluster_mass_proxy={}, const std::vector<double> cluster_mass_proxy_error={}, const std::string model_bias="Tinker", const std::string meanType="mean_bias", const int seed=666, const cbl::cosmology::Cosmology cosmology_mass={}, const std::vector<double> redshift_source={});
	
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
	 *  relation (see cbl::modelling::twopt::concentration)
	 *  
	 *  @param profile the density profile (see
	 *  cbl::modelling::twopt::concentration)
	 *
	 *  @param halo_def the halo definition (see
	 *  cbl::modelling::twopt::concentration)
	 *
	 *  @return none
	 */
	void set_data_HOD (const cbl::cosmology::Cosmology cosmology={}, const double redshift=0., const std::string model_MF="Tinker", const std::string model_bias="Tinker", const double Mh_min=0., const double Mh_max=1.e16, const double pi_max=100., const double r_max_int=100., const double r_min=1.e-3, const double r_max=350., const double k_min=1.e-4, const double k_max=100., const int step=200, const std::string method_Pk="CAMB", const bool NL=true, const std::string output_root="test", const double Delta=200., const double kk=0., const std::string interpType="Linear", const int norm=-1, const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true, const std::string model_cM="Duffy", const std::string profile="NFW", const std::string halo_def="vir");
	
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
	 *  @param output_dir the output_dir directory
	 *  where the output of external codes are written
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
	void set_data_model_cluster_selection_function (const cbl::cosmology::Cosmology cosmology, const cbl::cosmology::Cosmology test_cosmology, const double mean_redshift, const std::string model_MF, const std::string model_bias, const std::string selection_function_file, const std::vector<int> selection_function_column={}, const double z_min=par::defaultDouble, const double z_max=par::defaultDouble, const double Mass_min=par::defaultDouble, const double Mass_max=par::defaultDouble, const double Delta=200, const bool isDelta_vir=false, const std::string method_Pk="CAMB", const std::string output_dir=par::defaultString, const double k_min=1.e-4, const double k_max=100, const double prec=1.e-2, const int step=200, const int mass_step=50);

	///@}
	
      };
    }
  }
}

#endif
