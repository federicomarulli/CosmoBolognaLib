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
 *  @file Headers/Modelling_TwoPointCorrelation1D_angular.h
 *
 *  @brief The class Modelling_TwoPointCorrelation1D_angular
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation1D_angular, used to model the angular
 *  of two-point correlation function
 *
 *  @author Federico Marulli, Alfonso Veropalumbo, Massimiliano Romanello
 *
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it, massimilia.romanell2@unibo.it
 */

#ifndef __MODELLINGTWOPOINTANG__
#define __MODELLINGTWOPOINTANG__

#include "PkXi.h"
#include "Data1D_extra.h"
#include "TwoPointCorrelation1D_angular.h"
#include "Modelling_TwoPointCorrelation1D.h"
#include "ModelFunction_TwoPointCorrelation1D_angular.h"


// ===================================================================================================


namespace cbl {
  
  namespace modelling {

    namespace twopt {
    
      /**
       *  @class Modelling_TwoPointCorrelation1D_angular
       *  Modelling_TwoPointCorrelation1D_angular.h
       *  "Headers/Modelling_TwoPointCorrelation1D_angular.h"
       *
       *  @brief The class Modelling_TwoPointCorrelation1D_angular
       *
       *  This file defines the interface of the base class
       *  Modelling_TwoPointCorrelation1D_angular, used for modelling
       *  the angular two-point correlation function
       *
       */
      class Modelling_TwoPointCorrelation1D_angular : public Modelling_TwoPointCorrelation1D {
      protected:
	
	
      public:
	
	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 */
	Modelling_TwoPointCorrelation1D_angular () = default;

	/**
	 *  @brief constructor
	 *  
	 *  @param twop the two-point correlation function to model
	 */
	Modelling_TwoPointCorrelation1D_angular (const std::shared_ptr<cbl::measure::twopt::TwoPointCorrelation> twop)
	{m_data = twop->dataset();}
	
	/**
	 *  @brief constructor
	 *  
	 *  @param twop_dataset the dataset containing the two-point
	 *  correlation function to model
	 */
	Modelling_TwoPointCorrelation1D_angular (const std::shared_ptr<data::Data> twop_dataset)
	{m_data=twop_dataset;}
	
	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~Modelling_TwoPointCorrelation1D_angular () = default;
	
	///@}

	/**
	 *  @brief Set the data used to construct models of
	 *  the angular correlation function.
	 *
	 *  @param z_abs_err the redshift absolute error
	 *
	 *  @param cosmology the cosmological model used to measure angular correlation function
	 *
	 *  @param z_min the minimum redshift of the first redshift slice or of the full integral in one bin case
	 *
	 *  @param z_max the maximum redshift of the first redshift slice or of the full integral in one bin case
	 *
	 *  @param theta_units the coordinate units
	 *
	 *  @param method_Pk method used to compute the power spectrum
	 *  (i.e. the Boltzmann solver); valid choices for method_Pk
	 *  are: CAMB [http://camb.info/], CLASS
	 *  [http://class-code.net/], MPTbreeze-v1
	 *  [http://arxiv.org/abs/1207.1465], EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *
	 *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
	 *  non-linear power spectrum
	 *
	 *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
	 *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
	 *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
	 *
	 *  @param k_min minimum wave vector module up to which the
	 *  power spectrum is computed in order to estimate the power
	 *  spectrum normalisation; this parameter is used only if
	 *  either norm=1, or norm=-1 and sigma8 is set
	 *
	 *  @param k_max maximum wave vector module up to which the
	 *  power spectrum is computed to estimate the power spectrum
	 *  normalisation; this parameter is used only if norm=1
	 *
	 *  @param dN_par parameters of the dN/dz normalized distribution for the first redshift slice
	 *
	 *  @param distribution_type the shape of the dN/dz distribution.
	 *  Can be "polynomial", for a fit of the photometric redshift
	 *  distribution, "gaussian" or "skewnormal", in the correct
	 *  case of convolution with P(z true|z_phot)
	 *
	 *  @param redshift_space true->redshift space, include
	 *  redshift space distorsion, false->real space
	 *
	 *  @param interpType interpolation type
	 *
	 */
	void set_data_model (const double z_abs_err, const std::shared_ptr<cosmology::Cosmology> cosmology, const double z_min=0., const double z_max=10., const cbl::CoordinateUnits theta_units=cbl::CoordinateUnits::_arcminutes_, const std::string method_Pk="CAMB", const bool NL=false, const int norm=-1, const double k_min=0.001, const double k_max=5., const std::vector<double> dN_par={}, const std::string distribution_type="gaussian", const bool redshift_space=false, const std::string interpType="Spline");

	/**
	 *  @brief Set the data used to construct models of
	 *  the angular correlation function, in which the redshift selection
	 *  function is computed through a theoretical mass function.
	 *
	 *  @param z_abs_err the redshift absolute error. Zero for spectroscopic survey, >0 for photometric survey.
	 *
	 *  @param cosmology the cosmological model used to compute the power spectrum
	 *
	 *  @param z_min the minimum redshift of the photometric redshift slice
	 *
	 *  @param z_max the maximum redshift of the photometric redshift slice
	 *
	 *  @param theta_units the coordinate units
	 *
	 *  @param method_Pk method used to compute the power
	 *  spectrum; valid choices for method_Pk are: CAMB
	 *  [http://camb.info/], CLASS [http://class-code.net/],
	 *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465],
	 *  EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *
	 *  @param NL 0 &rarr; linear power spectrum; 1 &rarr;
	 *  non-linear power spectrum
	 *
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *
	 *  @param k_min minimum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *  
	 *  @param k_max maximum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *
	 *  @param area_degrees the area in deg. covered by the survey
	 *
	 *  @param redshift_space true->redshift space, false->real space
	 *
	 *  @param interpType the interpolation type 
	 *
	 *  @param Delta \f$\Delta\f$, the overdensity
	 *
	 *  @param isDelta_critical \f$\rightarrow\f$ \f$\Delta\f$ is the
	 *  overdensity defined with respect to the critical density
	 *
	 *  @param model_MF author(s) who proposed the mass function;
	 *  valid authors are: PS (Press & Schechter), ST(Sheth & Tormen),
	 *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
	 *  ShenH (halo MF, Shen et al. 2006), ShenF
	 *  (filament MF, Shen et al. 2006), ShenS
	 *  (sheet MF, Shen et al. 2006), Pan (Pan 2007), Peacock
	 *  (Peacock at al. 2007), Reed (Reed et al. 2007), Tinker
	 *  (Tinker et al. 2008), Crocce (Crocce et al. 2010), Manera
	 *  (Manera et al. 2010), Courtin (Courtin et al. 2010),
	 *  Bhattacharya (Bhattacharya et al. 2011), Angulo_FOF
	 *  (FoF MF, Angulo et al. 2012), Angulo_Sub
	 *  (SUBFIND MF, Angulo et al. 2012), Watson_FOF
	 *  (FoF MF, Watson et al. 2012),
	 *  Watson_SOH (Spherical Overdensity halo MF, Watson et al. 2012)
	 *, Despali_Z0, Despali_AllZ, Despali_AllZAllCosmo,
	 *  Despali_HighM (Despali et al. 2016)
	 *
	 *  @param Mass_min minimum mass
	 *  
	 *  @param Mass_max maximum mass
	 *
	 *  @param Mass_step the number of mass step
	 *
	 */
	void set_data_model_MF_theory (const double z_abs_err, const std::shared_ptr<cosmology::Cosmology> cosmology, const double z_min=par::defaultDouble, const double z_max=par::defaultDouble, const cbl::CoordinateUnits theta_units=cbl::CoordinateUnits::_arcminutes_, const std::string method_Pk="CAMB", const bool NL=false, const int norm=-1, const double k_min=0.001, const double k_max=5., const double area_degrees=cbl::par::defaultDouble, const bool redshift_space=false, const std::string interpType="Spline", const double Delta = 200., const bool isDelta_critical = true, const std::string model_MF = "Tinker", const double Mass_min = par::defaultDouble, const double Mass_max = par::defaultDouble, const int Mass_step = 100);

	/**
	 *  @brief Set the data used to construct models of
	 *  the angular correlation function, in which the bias and the redshift selection
	 *  function are computed theorically.
	 *
	 *  @param z_abs_err the redshift absolute error. Zero for spectroscopic survey, >0 for photometric survey.
	 *
	 *  @param cosmology the cosmological model used to compute the power spectrum
	 *
	 *  @param cluster_mass vector conteining the masses of the catalogue
	 *
	 *  @param cluster_redshift vector conteining the redshifts of the catalogue
	 *
	 *  @param z_min the minimum redshift of the photometric redshift slice
	 *
	 *  @param z_max the maximum redshift of the photometric redshift slice
	 *
	 *  @param theta_units the coordinate units
	 *
	 *  @param method_Pk method used to compute the power
	 *  spectrum; valid choices for method_Pk are: CAMB
	 *  [http://camb.info/], CLASS [http://class-code.net/],
	 *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465],
	 *  EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *
	 *  @param NL 0 &rarr; linear power spectrum; 1 &rarr;
	 *  non-linear power spectrum
	 *
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *
	 *  @param k_min minimum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *  
	 *  @param k_max maximum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *
	 *  @param output_dir the output_dir directory
	 *  where the output of external codes are written
	 *
	 *  @param output_root output_root of the parameter file used
	 *  to compute the power spectrum and &sigma;(mass); it can be
	 *  any name
	 *
	 *  @param area_degrees the area in deg. covered by the survey
	 *
	 *  @param redshift_space true->redshift space, false->real space
	 *
	 *  @param interpType the interpolation type 
	 *
	 *  @param Delta \f$\Delta\f$, the overdensity
	 *
	 *  @param isDelta_critical \f$\rightarrow\f$ \f$\Delta\f$ is the
	 *  overdensity defined with respect to the critical density
	 *
	 *  @param model_MF author(s) who proposed the mass function;
	 *  valid authors are: PS (Press & Schechter), ST(Sheth & Tormen),
	 *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
	 *  ShenH (halo MF, Shen et al. 2006), ShenF
	 *  (filament MF, Shen et al. 2006), ShenS
	 *  (sheet MF, Shen et al. 2006), Pan (Pan 2007), Peacock
	 *  (Peacock at al. 2007), Reed (Reed et al. 2007), Tinker
	 *  (Tinker et al. 2008), Crocce (Crocce et al. 2010), Manera
	 *  (Manera et al. 2010), Courtin (Courtin et al. 2010),
	 *  Bhattacharya (Bhattacharya et al. 2011), Angulo_FOF
	 *  (FoF MF, Angulo et al. 2012), Angulo_Sub
	 *  (SUBFIND MF, Angulo et al. 2012), Watson_FOF
	 *  (FoF MF, Watson et al. 2012),
	 *  Watson_SOH (Spherical Overdensity halo MF, Watson et al. 2012)
	 *, Despali_Z0, Despali_AllZ, Despali_AllZAllCosmo,
	 *  Despali_HighM (Despali et al. 2016)
	 *
	 *  @param model_bias author(s) who proposed the bias; valid
	 *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
	 *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
	 *  correction of Warren 2004), Tinker (Tinker et al. 2010)
	 *
	 *  @param Mass_step the number of mass step
	 *
	 *  @param bias_from_sample if true compute the bias as an average over the clusters bias in the catalogue
	 *
	 *  @param include_quadrupole if true, include the quadrupole of the 3D correlation function in the model. It works only in redshift space
	 */
	void set_data_model_bias_MF_theory (const double z_abs_err, const std::shared_ptr<cosmology::Cosmology> cosmology, const std::vector<double> cluster_mass, const std::vector<double> cluster_redshift, const double z_min=par::defaultDouble, const double z_max=par::defaultDouble, const cbl::CoordinateUnits theta_units=cbl::CoordinateUnits::_arcminutes_, const std::string method_Pk="CAMB", const bool NL=false, const int norm=-1, const double k_min=0.001, const double k_max=100., const std::string output_dir=par::defaultString, const std::string output_root="test", const double area_degrees=cbl::par::defaultDouble, const bool redshift_space=false, const std::string interpType="Spline", const double Delta = 200., const bool isDelta_critical = true, const std::string model_MF = "Tinker", const std::string model_bias = "Tinker", const int Mass_step = 100, const bool bias_from_sample=false, const bool include_quadrupole=false);
	
	/**
	 *  @brief Set the data used to construct models of
	 *  the two-point correlation function, in which the effective bias
	 *  is computed through a mass-observable scaling relation.
	 *
	 *  In particular, the cluster masses are derived from the following
	 *  scaling relation:
	 *
	 *  \f$ \log (M/M_{\rm piv}) = \alpha + \beta\log(\lambda/\lambda_{\rm piv})+
	 *  \gamma \log f(z;z_{\rm piv}) \f$
	 *
	 *  See cbl::modelling::massobsrel::Modelling_MassObservableRelation for details.
	 *
	 *  @param z_abs_err absolute redshift error 
	 *
	 *  @param proxy_rel_err proxy relative error
	 *
	 *  @param cosmology the cosmological model used to measure
	 *  &xi;<SUB>DM</SUB>
	 *
	 *  @param redshift redshift
	 *
	 *  @param cluster_redshift vector containing the cluster
	 *  redshifts
	 *
	 *  @param cluster_mass_proxy vector containing the cluster
	 *  mass proxies
	 *
	 *  @param redshift_pivot redshift pivot in the scaling relation
	 *
	 *  @param proxy_pivot proxy pivot in the scaling relation
	 *
	 *  @param mass_pivot mass pivot in the scaling relation
	 *
	 *  @param log_base base of the logarithms in the scaling relation
	 *
	 *  @param method_Pk method used to compute the power
	 *  spectrum; valid choices for method_Pk are: CAMB
	 *  [http://camb.info/], CLASS [http://class-code.net/],
	 *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465],
	 *  EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *
	 *  @param NL 0 &rarr; linear power spectrum; 1 &rarr;
	 *  non-linear power spectrum
	 *    
	 *  @param k_min minimum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *  
	 *  @param k_max maximum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *  
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *
	 *  @param redshift_space true->redshift space, false->real space
	 *
	 *  @param dN_par parameters of the dN/dz normalized distribution for the first redshift slice
	 *
	 *  @param distribution_type the shape of the dN/dz distribution.
	 *  Can be "polynomial", for a fit of the photometric redshift
	 *  distribution, "gaussian" or "skewnormal", in the correct
	 *  case of convolution with P(z true|z_phot)
	 *
	 *  @param z_min the minimum redshift of the photometric redshift slice
	 *
	 *  @param z_max the maximum redshift of the photometric redshift slice
	 *  @param theta_units angular units for \f$w(\theta)\f$
	 *
	 *  @param output_dir the output_dir directory
	 *  where the output of external codes are written
	 *
	 *  @param output_root output_root of the parameter file used
	 *  to compute the power spectrum and &sigma;(mass); it can be
	 *  any name
	 *
	 *  @param prec accuracy of the GSL integration
	 *  
	 *  @param file_par name of the parameter file; if a parameter
	 *  file is provided (i.e. file_par!=NULL), it will be used,
	 *  ignoring the cosmological parameters of the object
	 *
	 *  @param Delta \f$\Delta\f$, the overdensity
	 *
	 *  @param isDelta_critical \f$\rightarrow\f$ \f$\Delta\f$ is the
	 *  overdensity defined with respect to the critical density
	 *
	 *  @param model_bias author(s) who proposed the bias; valid
	 *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
	 *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
	 *  correction of Warren 2004), Tinker (Tinker et al. 2010)
	 *
	 *  @param interpType the interpolation type 
	 *
	 */
	void set_data_model_scaling_relation_cosmology (const double z_abs_err, const double proxy_rel_err, const std::shared_ptr<cosmology::Cosmology> cosmology, const double redshift, const std::vector<double> cluster_redshift, const std::vector<double> cluster_mass_proxy, const double redshift_pivot, const double proxy_pivot, double mass_pivot, const double log_base, const std::string method_Pk="CAMB", const bool NL=true, const double k_min=1.e-4, const double k_max=100., const int norm=-1, const bool redshift_space=false, std::vector<double> dN_par={}, const std::string distribution_type="gaussian", const double z_min=par::defaultDouble, const double z_max=par::defaultDouble, const cbl::CoordinateUnits theta_units=cbl::CoordinateUnits::_arcminutes_, const std::string output_dir=par::defaultString, const std::string output_root="test", const double prec=1.e-3, const std::string file_par=par::defaultString, const double Delta=200., const bool isDelta_critical=true, const std::string model_bias="Tinker", const std::string interpType="Spline");
	
	/**
	 *  @brief Set the data used to construct models of
	 *  the two-point correlation function, in which the effective bias
	 *  is computed through a mass-observable scaling relation.
	 *
	 *  In particular, the cluster masses are derived from the following
	 *  scaling relation:
	 *
	 *  \f$ \log (M/M_{\rm piv}) = \alpha + \beta\log(\lambda/\lambda_{\rm piv})+
	 *  \gamma \log f(z;z_{\rm piv}) \f$
	 *
	 *  See cbl::modelling::massobsrel::Modelling_MassObservableRelation for details.
	 *
	 *  @param z_abs_err the redshift absolute error
	 *
	 *  @param proxy_rel_err the proxy relative error
	 *
	 *  @param cosmology the cosmological model used to compute the power spectrum
	 *
	 *  @param redshift redshift
	 *
	 *  @param cluster_redshift vector containing the cluster
	 *  redshifts
	 *
	 *  @param cluster_mass_proxy vector containing the cluster
	 *  mass proxies
	 *
	 *  @param redshift_pivot redshift pivot in the scaling relation
	 *
	 *  @param proxy_pivot proxy pivot in the scaling relation
	 *
	 *  @param mass_pivot mass pivot in the scaling relation
	 *
	 *  @param log_base base of the logarithms in the scaling relation
	 *
	 *  @param method_Pk method used to compute the power
	 *  spectrum; valid choices for method_Pk are: CAMB
	 *  [http://camb.info/], CLASS [http://class-code.net/],
	 *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465],
	 *  EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *
	 *  @param NL 0 &rarr; linear power spectrum; 1 &rarr;
	 *  non-linear power spectrum
	 *    
	 *  @param k_min minimum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *  
	 *  @param k_max maximum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *
	 *  @param redshift_space true->redshift space, false->real space
	 *
	 *  @param area_degrees the area in deg. covered by the survey
	 *
	 *  @param z_min the minimum redshift of the photometric redshift slice
	 *
	 *  @param z_max the maximum redshift of the photometric redshift slice
	 *
	 *  @param proxy_min the minimum proxy of the photometric redshift slice
	 *
	 *  @param proxy_max the maximum proxy of the photometric redshift slice
	 *
	 *  @param theta_units angular units for \f$w(\theta)\f$
	 *
	 *  @param output_dir the output_dir directory
	 *  where the output of external codes are written
	 *
	 *  @param output_root output_root of the parameter file used
	 *  to compute the power spectrum and &sigma;(mass); it can be
	 *  any name
	 *
	 *  @param prec accuracy of the GSL integration
	 *  
	 *  @param file_par name of the parameter file; if a parameter
	 *  file is provided (i.e. file_par!=NULL), it will be used,
	 *  ignoring the cosmological parameters of the object
	 *
	 *  @param Delta \f$\Delta\f$, the overdensity
	 *
	 *  @param isDelta_critical \f$\rightarrow\f$ \f$\Delta\f$ is the
	 *  overdensity defined with respect to the critical density
	 *
	 *  @param model_bias author(s) who proposed the bias; valid
	 *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
	 *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
	 *  correction of Warren 2004), Tinker (Tinker et al. 2010)
	 *
	 *  @param interpType interpolation type
	 *
	 *  @param model_MF author(s) who proposed the mass function;
	 *  valid authors are: PS (Press & Schechter), ST(Sheth & Tormen),
	 *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
	 *  ShenH (halo MF, Shen et al. 2006), ShenF
	 *  (filament MF, Shen et al. 2006), ShenS
	 *  (sheet MF, Shen et al. 2006), Pan (Pan 2007), Peacock
	 *  (Peacock at al. 2007), Reed (Reed et al. 2007), Tinker
	 *  (Tinker et al. 2008), Crocce (Crocce et al. 2010), Manera
	 *  (Manera et al. 2010), Courtin (Courtin et al. 2010),
	 *  Bhattacharya (Bhattacharya et al. 2011), Angulo_FOF
	 *  (FoF MF, Angulo et al. 2012), Angulo_Sub
	 *  (SUBFIND MF, Angulo et al. 2012), Watson_FOF
	 *  (FoF MF, Watson et al. 2012),
	 *  Watson_SOH (Spherical Overdensity halo MF, Watson et al. 2012)
	 *, Despali_Z0, Despali_AllZ, Despali_AllZAllCosmo,
	 *  Despali_HighM (Despali et al. 2016)
	 *
	 *  @param SF_weights weights for completeness and purity of number counts in lambda* redshift bins
	 *
	 */
	void set_data_model_scaling_relation_cosmology_MF_theory (const double z_abs_err, const double proxy_rel_err, const std::shared_ptr<cosmology::Cosmology> cosmology, const double redshift, const std::vector<double> cluster_redshift, const std::vector<double> cluster_mass_proxy, const double redshift_pivot, const double proxy_pivot, double mass_pivot, const double log_base, const std::string method_Pk="CAMB", const bool NL=true, const double k_min=1.e-4, const double k_max=100., const int norm=-1, const bool redshift_space=false, const double area_degrees=cbl::par::defaultDouble, const double z_min=par::defaultDouble, const double z_max=par::defaultDouble, const double proxy_min=par::defaultDouble, const double proxy_max=par::defaultDouble, const cbl::CoordinateUnits theta_units=cbl::CoordinateUnits::_arcminutes_, const std::string output_dir=par::defaultString, const std::string output_root="test", const double prec=1.e-3, const std::string file_par=par::defaultString, const double Delta=200., const bool isDelta_critical=true, const std::string model_bias="Tinker", const std::string interpType="Spline", const std::string model_MF = "Tinker", const std::vector<double> SF_weights={1.});
		
	/**
	 *  @brief set the model to fit the angular power spectrum
	 *
	 *  The model is based on the redshift distribution obtained directly from the data (e.g. galaxy catalogues). The true redshift distribution in a given redshift bin is obtained as: 

	 \f[\frac{dN}{dz_{true}} = \int_{0}^{\infty} dz_{phot} W(z_{phot}) \frac{dN}{dz_{phot}} P(z_{true}|z_{phot})  \f]
	 *
	 *  where \f$ W(z_{phot}) \f$ is a top-hat window function in the given photometric bin and \f$ P(z_{true}|z_{phot}) \f$ is assumed as gaussian with zero mean and std dev given by \f$ \sigma_{0,z}(1+z_{phot}) \f$. See Balaguera et al. (2018)
	 *
	 *  The model computes:
	 *
	 *  \f[w(\theta) = \int\int dz_1 dz_2 \phi(z_1) \phi(z_2) \xi(s), \f]
	 *
	 *  where,
	 *
	 *  \f$ s=\sqrt{r^2(z_1)+r^2(z_2)-2r_1r_2\cos(\theta)} \f$
	 *  and: 
	 *  \f[ \phi(z)=\frac{\frac{dN(z)}{dz}}{\int dz \frac{dN(z)}{dz}}\f]
	 *  
	 *  Redshift-space distorsions are modelled in the Kaiser
	 *  limit, that is neglecting non-linearities in dynamics
	 *  and bias. Specifically:
	 * 
	 *  \f[\xi(s) = \left[ b^2 + \frac{2}{3}
	 *  f \cdot b + \frac{1}{5}f^2
	 *  \right] \cdot \xi_{\rm DM}(s) \f]
	 *
	 *  the model has N cosmological parameters, plus 1 parameter: 
	 *    - \f$ b \f$
	 *    
	 *  the parameters of the true dN/dz distribution ("gaussian", "skewnormal" or "polynomial") in a given photometric redshift bin are fixed
	 *  by the user.
	 *
	 *  @param cosmo_param list of cosmological parameter for the
	 *  model
	 *
	 *  @param cosmo_prior prior for the cosmological parameters
	 *
	 *  @param bias_prior prior for the parameter \f$ b \f$
	 */
	void set_model (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution bias_prior);
	
	/**
	 *  @brief set the model to fit the angular correlation function,
	 *  with theoretical mass function.
	 *  The model computes:
	 *
	 *  \f[w(\theta) = \int\int dz_1 dz_2 \phi(z_1) \phi(z_2) \xi(s), \f]
	 *
	 *  where,
	 *
	 *  \f$ s=\sqrt{r^2(z_1)+r^2(z_2)-2r_1r_2\cos(\theta)} \f$
	 *  and: 
	 *  \f[ \phi(z)=\frac{\frac{dN(z)}{dz}}{\int dz \frac{dN(z)}{dz}}\f]
	 *  
	 *  Redshift-space distorsions are modelled in the Kaiser
	 *  limit, that is neglecting non-linearities in dynamics
	 *  and bias. Specifically:
	 * 
	 *  \f[\xi(s) = \left[ b^2 + \frac{2}{3}
	 *  f \cdot b + \frac{1}{5}f^2
	 *  \right] \cdot \xi_{\rm DM}(s) \f]
	 *
	 *  where \f$ b \f$ is the bias.
	 *
	 *  *  the model has N cosmological parameters plus 1 parameter: 
	 *    - \f$ b \f$
	 *    
	 *  the dN/dz distribution in computed from the theoretical mass
	 *  function:
	 *
	 *  \f[ \frac{{\rm d}N }{{\rm d} z_{\rm tr}}(\Delta z_{\text{ob}})=\,\,\Omega 
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{\Delta z_{\text{ob}}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  , \f]
	 *  the angular correlation function is computed
	 *  using the input cosmological parameters
	 *
	 *  @param cosmo_param list of cosmological parameter for the model
	 *
	 *  @param cosmo_prior prior for the cosmological parameters
	 *
	 *  @param bias_prior prior for the parameter \f$ b \f$
	 *
	 */
	void set_model_MF_theory (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution bias_prior);
	
	/**
	 *  @brief set the model to fit the angular correlation function,
	 *  with theoretical mass function.
	 *  The model computes:
	 *
	 *  \f[w(\theta) = \int\int dz_1 dz_2 \phi(z_1) \phi(z_2) \xi(s), \f]
	 *
	 *  where,
	 *
	 *  \f$ s=\sqrt{r^2(z_1)+r^2(z_2)-2r_1r_2\cos(\theta)} \f$
	 *  and: 
	 *  \f[ \phi(z)=\frac{\frac{dN(z)}{dz}}{\int dz \frac{dN(z)}{dz}}\f]
	 *  
	 *  Redshift-space distorsions are modelled in the Kaiser
	 *  limit, that is neglecting non-linearities in dynamics
	 *  and bias, only if the user choose redshift space in set_data_model_bias. 
	 *  Specifically:
	 * 
	 *  \f[\xi(s) = \left[ b^2 + \frac{2}{3}
	 *  f \cdot b + \frac{1}{5}f^2
	 *  \right] \cdot \xi_{\rm DM}(s) \f]
	 *
	 *  where \f$ b \f$ is the linear effective bias, computed theoretically for each cosmology in two ways: 
	 *  using the provided halo masses with cbl::cosmology::bias_eff_mass, 
	 *  or with cbl::cosmology::bias_eff_selection_function The dN/dz distribution is computed 
	 *  from the theoretical mass function:
	 *
	 *  \f[ \frac{{\rm d}N }{{\rm d} z_{\rm tr}}(\Delta z_{\text{ob}})=\,\,\Omega 
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{\Delta z_{\text{ob}}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  , \f]
	 *  the angular correlation function is computed
	 *  using the input cosmological parameters
	 *
	 *  The model has N cosmological parameters. The bias is considered a derived parameter, from cosmology
	 *  @param cosmo_param list of cosmological parameter for the model
	 *
	 *  @param cosmo_prior prior for the cosmological parameters
	 *
	 */
	void set_model_bias_MF_theory (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior);
	
	/**
	 *  @brief set the model to fit the angular correlation function,
	 *  with theoretical mass function.
	 *  The model computes:
	 *
	 *  \f[w(\theta) = \int\int dz_1 dz_2 \phi(z_1) \phi(z_2) \xi(s), \f]
	 *
	 *  where,
	 *
	 *  \f$ s=\sqrt{r^2(z_1)+r^2(z_2)-2r_1r_2\cos(\theta)} \f$
	 *  and: 
	 *  \f[ \phi(z)=\frac{\frac{dN(z)}{dz}}{\int dz \frac{dN(z)}{dz}}\f]
	 *  
	 *  Redshift-space distorsions are modelled in the Kaiser
	 *  limit, that is neglecting non-linearities in dynamics
	 *  and bias, only if the user choose redshift space in set_data_model_bias. Specifically:
	 * 
	 *  \f[\xi_0(s) = b^2 \left[ 1 + \frac{2 \beta}{3}
	 *  \cdot + \frac{\beta^2}{5} \right] \cdot
	 *  \xi_{\rm DM}(s, \Sigma_{NL})\ \f]
	 *
	 *  where \f$ b \f$ is the linear effective bias, computed theoretically for each cosmology in two ways: 
	 *  using the provided halo masses with cbl::cosmology::bias_eff_mass, 
	 *  or with cbl::cosmology::bias_eff_selection_function The dN/dz distribution is computed 
	 *  from the theoretical mass function:
	 *
	 *  \f[ \frac{{\rm d}N }{{\rm d} z_{\rm tr}}(\Delta z_{\text{ob}})=\,\,\Omega 
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{\Delta z_{\text{ob}}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  , \f]
	 *  the angular correlation function is computed
	 *  using the input cosmological parameters
	 *
	 *  The model has N cosmological parameters, plus: 
	 *    - \f$\Sigma_{NL}\f$ damping at BAO
	 *
	 *  The bias is considered a derived parameter, from cosmology
	 *
	 *  @param cosmo_param list of cosmological parameter for the model
	 *
	 *  @param cosmo_prior prior for the cosmological parameters
	 *
	 *  @param SigmaNL_prior prior for the parameter \f$\Sigma_{NL}\f$
	 *
	 */
	void set_model_bias_MF_theory_BAO (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution SigmaNL_prior);
	
	/**
	 *  @brief set the model to fit the angular correlation function,
	 *  with theoretical mass function.
	 *  The model computes:
	 *
	 *  \f[w(\theta) = \int\int dz_1 dz_2 \phi(z_1) \phi(z_2) \xi(s), \f]
	 *
	 *  where,
	 *
	 *  \f$ s=\sqrt{r^2(z_1)+r^2(z_2)-2r_1r_2\cos(\theta)} \f$
	 *  and: 
	 *  \f[ \phi(z)=\frac{\frac{dN(z)}{dz}}{\int dz \frac{dN(z)}{dz}}\f]
	 *  
	 *  Redshift-space distorsions are modelled in the Kaiser
	 *  limit, that is neglecting non-linearities in dynamics
	 *  and bias, only if the user choose redshift space in set_data_model_bias. Specifically:
	 * 
	 *  \f[\xi_0(s) = b^2 \left[ 1 + \frac{2 \beta}{3}
	 *  \cdot + \frac{\beta^2}{5} \right] \cdot
	 *  \xi_{\rm DM}(s, \Sigma_{NL})\ + A_0 + A_1/s
	 *  +A_2/s^2\f]
	 *
	 *  where \f$ b \f$ is the linear effective bias, computed theoretically for each cosmology in two ways: 
	 *  using the provided halo masses with cbl::cosmology::bias_eff_mass, 
	 *  or with cbl::cosmology::bias_eff_selection_function The dN/dz distribution is computed 
	 *  from the theoretical mass function:
	 *
	 *  \f[ \frac{{\rm d}N }{{\rm d} z_{\rm tr}}(\Delta z_{\text{ob}})=\,\,\Omega 
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{\Delta z_{\text{ob}}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  , \f]
	 *  the angular correlation function is computed
	 *  using the input cosmological parameters
	 *
	 *  The model has N cosmological parameters, plus: 
	 *    - \f$A_0, A_1, A2\f$ parameters of the polynomial
	 *    - \f$\Sigma_{NL}\f$ damping at BAO
	 *
	 *  The bias is considered a derived parameter, from cosmology
	 *
	 *  @param cosmo_param list of cosmological parameter for the model
	 *
	 *  @param cosmo_prior prior for the cosmological parameters
	 *
	 *  @param SigmaNL_prior prior for the parameter \f$\Sigma_{NL}\f$
	 *
	 *  @param A0_prior prior for the parameter \f$A_0\f$
	 *
	 *  @param A1_prior prior for the parameter \f$A_1\f$
	 *
	 *  @param A2_prior prior for the parameter \f$A_2\f$
	 *
	 */
	void set_model_bias_MF_theory_BAO_poly (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution SigmaNL_prior, const statistics::PriorDistribution A0_prior, const statistics::PriorDistribution A1_prior, const statistics::PriorDistribution A2_prior);
	
	/**
	 *  @brief Set the parameters to model the monopole of the
	 *  two-point correlation function in redshift space, where
	 *  the effective bias can be computed in two ways, depending
	 *  on the cbl::modelling::twopt::Modelling_TwoPointCorrelation1D::set\_data\_model used.
	 *  In one case, the masses are given by a mass-observable scaling relation
	 *  with the following functional form:
	 *
	 *  \f[ \log \frac{M}{M_{\rm piv}} = \alpha+
	 *  \beta  \log \frac{\lambda}{\lambda_{\rm piv}}+
	 *  \gamma  \log f(z;z_{\rm piv}) + \sigma_{\rm intr}, \f]
	 *
	 *  where \f$\lambda\f$ is the mass proxy and
	 *  the intrinsic scatter, \f$\sigma_{\rm intr}\f$, should
	 *  be set as both negative and positive and it has the 
	 *  following functional form:
	 *
	 *  \f$ \sigma_{\rm intr} = \sigma_0 + \sigma_{M} 
	 *  \log (M/M_{\rm piv})^{e_{M}} + \sigma_z \log (f(z))^{e_z}. \f$
	 *
	 *  The bias is then computed for each object in correspondence of
	 *  such masses, and then averaged to obtain the effective bias.
	 *
	 *  In the second case, the effective bias is derived as follows:
	 *
	 *  \f$ b_{\rm eff} = \frac{1}{N}\sum_{i=1}^N b(z_{{\rm ob},i},\lambda_{{\rm ob},i}) = 
	 *  \frac{1}{N}\sum_{i=1}^N\int_0^\infty {\rm d}z\,\int_0^\infty{\rm d}\lambda\,\int_0^\infty{\rm d}M\,\, b(M,z) \, 
	 *  P(M|\lambda,z) \, P(z|z_{{\rm ob},i}) \, P(\lambda|\lambda_{{\rm ob},i}), \f$
	 *
	 *  where \f$N\f$ is the number of objects, \f$z_{{\rm ob},i}\f$ and \f$\lambda_{{\rm ob},i}\f$
	 *  are the observed redshift and mass proxy of the \f$i\f$th object, respectively,
	 *  \f$P(M|\lambda,z)\f$ is a log-normal whose mean is the mass proxy - mass relation and whose
	 *  rms is the intrinsic scatter of such relation (for details, see e.g. 
	 *  cbl::modelling::massobsrel::Modelling_MassObservableRelation), while \f$P(z|z_{{\rm ob},i})\f$
	 *  and \f$P(\lambda|\lambda_{{\rm ob},i})\f$ are Gaussian distributions whose standard deviations
	 *  are given by the errors on redshift and proxy, respectively.
	 *
	 *  The model computes:
	 *
	 *  \f[w(\theta) = \int\int dz_1 dz_2 \phi(z_1) \phi(z_2) \xi(s), \f]
	 *
	 *  where,
	 *
	 *  \f$ s=\sqrt{r^2(z_1)+r^2(z_2)-2r_1r_2\cos(\theta)} \f$
	 *  and: 
	 *  \f[ \phi(z)=\frac{\frac{dN(z)}{dz}}{\int dz \frac{dN(z)}{dz}}\f]
	 *  
	 *  Redshift-space distorsions are modelled in the Kaiser
	 *  limit, that is neglecting non-linearities in dynamics
	 *  and bias. Specifically:
	 * 
	 *  \f[\xi(s) = \left[ b_{eff}^2 + \frac{2}{3}
	 *  f \cdot b_{eff} + \frac{1}{5}f^2
	 *  \right] \cdot \xi_{\rm DM}(s) \f]
	 *
	 *  where \f$ b_{eff} \f$ is the effective bias.
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_prior vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @param alpha_prior prior on the scaling relation normalization
	 *
	 *  @param beta_prior prior on the scaling relation slope
	 *
	 *  @param gamma_prior prior on the redshift evolution factor of the scaling relation
	 *
	 *  @param scatter0_prior prior on the 
	 *  constant term of the intrinsic scatter, \f$ \sigma_0 \f$
	 *
	 *  @param scatterM_prior prior on the factor in the
	 *  proxy-dependent term of the intrinsic scatter, \f$ \sigma_{\lambda} \f$
	 *
	 *  @param scatterM_exponent_prior prior on the exponent in the
	 *  proxy-dependent term of the intrinsic scatter, \f$ e_{\lambda} \f$
	 *
	 *  @param scatterz_prior prior on the factor in the
	 *  redshift-dependent term of the intrinsic scatter, \f$ \sigma_z \f$
	 *
	 *  @param scatterz_exponent_prior prior on the exponent in the
	 *  redshift-dependent term of the intrinsic scatter, \f$ e_z \f$
	 *
	 *  @param z_evo functional form of the redshift evolution
	 *  function in the scaling relation. See the documentation
	 *  of cbl::modelling::massobsrel::Modelling_MassObservableRelation
	 *  
	 */
	void set_model_scaling_relation_cosmology (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const std::string z_evo);

	/**
	 *  @brief Set the parameters to model the monopole of the
	 *  two-point correlation function in redshift space, where
	 *  the effective bias can be computed in two ways, depending
	 *  on the
	 *  cbl::modelling::twopt::Modelling_TwoPointCorrelation1D::set\_data\_model
	 *  used.  In one case, the masses are given by a
	 *  mass-observable scaling relation with the following
	 *  functional form:
	 *
	 *  \f[ \log \frac{M}{M_{\rm piv}} = \alpha+
	 *  \beta  \log \frac{\lambda}{\lambda_{\rm piv}}+
	 *  \gamma  \log f(z;z_{\rm piv}) + \sigma_{\rm intr}, \f]
	 *
	 *  where \f$\lambda\f$ is the mass proxy and
	 *  the intrinsic scatter, \f$\sigma_{\rm intr}\f$, should
	 *  be set as both negative and positive and it has the 
	 *  following functional form:
	 *
	 *  \f$ \sigma_{\rm intr} = \sigma_0 + \sigma_{M} \log
	 *  (M/M_{\rm piv})^{e_{M}} + \sigma_z \log (f(z))^{e_z}. \f$
	 *
	 *  The bias is then computed for each object in correspondence of
	 *  such masses, and then averaged to obtain the effective bias.
	 *
	 *  In the second case, the effective bias is derived as follows:
	 *
	 *  \f$ b_{\rm eff} = \frac{1}{N}\sum_{i=1}^N b(z_{{\rm
	 *  ob},i},\lambda_{{\rm ob},i}) =
	 *  \frac{1}{N}\sum_{i=1}^N\int_0^\infty {\rm
	 *  d}z\,\int_0^\infty{\rm d}\lambda\,\int_0^\infty{\rm
	 *  d}M\,\, b(M,z) \, P(M|\lambda,z) \, P(z|z_{{\rm ob},i}) \,
	 *  P(\lambda|\lambda_{{\rm ob},i}), \f$
	 *
	 *  where \f$N\f$ is the number of objects, \f$z_{{\rm
	 *  ob},i}\f$ and \f$\lambda_{{\rm ob},i}\f$ are the observed
	 *  redshift and mass proxy of the \f$i\f$th object,
	 *  respectively, \f$P(M|\lambda,z)\f$ is a log-normal whose
	 *  mean is the mass proxy - mass relation and whose rms is
	 *  the intrinsic scatter of such relation (for details, see
	 *  e.g.
	 *  cbl::modelling::massobsrel::Modelling_MassObservableRelation),
	 *  while \f$P(z|z_{{\rm ob},i})\f$ and
	 *  \f$P(\lambda|\lambda_{{\rm ob},i})\f$ are Gaussian
	 *  distributions whose standard deviations are given by the
	 *  errors on redshift and proxy, respectively.
	 * 
	 *  The model computes:
	 *
	 *  \f[w(\theta) = \int\int dz_1 dz_2 \phi(z_1) \phi(z_2)
	 *  \xi(s), \f]
	 *
	 *  where,
	 *
	 *  \f$ s=\sqrt{r^2(z_1)+r^2(z_2)-2r_1r_2\cos(\theta)} \f$
	 *  and: \f[ \phi(z)=\frac{\frac{dN(z)}{dz}}{\int dz
	 *  \frac{dN(z)}{dz}}\f]
	 *  
	 *  Finally, \f[ \frac{{\rm d}N }{{\rm d} z_{\rm tr}}(\Delta
	 *  z_{\text{ob}, j})=\,\,\Omega \frac{{\rm d} V}{{\rm d}
	 *  z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} {\rm d} M_{\rm
	 *  tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d}
	 *  M_{\rm tr}}\,\, \int_{0}^{\infty}{\rm d}\lambda_{\rm
	 *  tr}\,\, P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\,
	 *  \int_{\Delta z_{\text{ob},j}}{\rm d} z_{\rm ob}
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  \int_{\Delta\lambda_{\text{ob},i}}{\rm d} \lambda_{\rm ob}
	 *  \,\,P(\lambda_{\rm ob}|\lambda_{\rm
	 *  tr})w(\lambda_{ob},\Delta z_{\text{ob},j}), \f]
	 *
	 *  as expressed in
	 *  cbl::modelling::numbercounts::Modelling_NumberCounts1D_MassProxy::set_model_NumberCounts_cosmology
	 *
	 *  Redshift-space distorsions are modelled in the Kaiser
	 *  limit, that is neglecting non-linearities in dynamics
	 *  and bias. Specifically:
	 * 
	 *  \f[\xi(s) = \left[ b_{eff}^2 + \frac{2}{3}
	 *  f \cdot b_{eff} + \frac{1}{5}f^2
	 *  \right] \cdot \xi_{\rm DM}(s) \f]
	 *
	 *  where \f$ b_{eff} \f$ is the effective bias.
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_prior vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @param alpha_prior prior on the scaling relation
	 *  normalization
	 *
	 *  @param beta_prior prior on the scaling relation slope
	 *
	 *  @param gamma_prior prior on the redshift evolution factor
	 *  of the scaling relation
	 *
	 *  @param scatter0_prior prior on the 
	 *  constant term of the intrinsic scatter, \f$ \sigma_0 \f$
	 *
	 *  @param scatterM_prior prior on the factor in the
	 *  proxy-dependent term of the intrinsic scatter, \f$
	 *  \sigma_{\lambda} \f$
	 *
	 *  @param scatterM_exponent_prior prior on the exponent in the
	 *  proxy-dependent term of the intrinsic scatter, \f$ e_{\lambda} \f$
	 *
	 *  @param scatterz_prior prior on the factor in the
	 *  redshift-dependent term of the intrinsic scatter, \f$
	 *  \sigma_z \f$
	 *
	 *  @param scatterz_exponent_prior prior on the exponent in the
	 *  redshift-dependent term of the intrinsic scatter, \f$ e_z \f$
	 *
	 *  @param Plambda_prior priors on the three parameters,
	 *  respectively a, b, c, defining the distribution
	 *  \f$P(\lambda_{tr}|z_{tr})\f$
	 *
	 *  @param z_evo functional form of the redshift evolution
	 *  function in the scaling relation. See the documentation of
	 *  cbl::modelling::massobsrel::Modelling_MassObservableRelation
	 */
	void set_model_scaling_relation_cosmology_MF_theory (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const std::vector<statistics::PriorDistribution> Plambda_prior, const std::string z_evo);

	/**
	 *  @brief set the model to fit the covariance of the angular correlation function
	 *
	 *  The model is the following:
	 *
	 *  \f[ \sigma_{w(\theta)} = \frac{(1+\alpha)^2}{4\pi^2 f_{sky} \overline{n}^2 (\mu_- - \mu_+)} + \sum_\ell \frac{2\ell+1}{8 \pi^2 f_{sky}} \hat{L_\ell}(\mu) \hat{L_\ell}(\mu) \left(\beta^2 C_\ell^2 + \frac{2 \beta C_\ell (1+\alpha)}{\overline{n}} \right) \f]
	 *
	 *  where \f$\overline{n}=N_{obj}/\Omega_{sky}\f$ and \f$\hat{L_\ell}(\mu)\f$ are
	 *  the band-averaged Legendre polynomials of degree \f$\ell\f$, computed as:
	 *
	 *  \f[ \hat{L_\ell} = \frac{L_{\ell+1}(\mu_+) - L_{\ell+1}(\mu_-) - L_{\ell-1}(\mu_+) + L_{\ell-1}(\mu_-)} {(2\ell+1) (\mu_- - \mu_+)}, \f]
	 *
	 *   \f$ \theta_\pm = \theta \pm \Delta\theta/2 \f$ represent the upper and the lower limits
	 *   of the angular bins, while \f$ \mu_\pm \f$ are the corresponding cosines.
	 *   The angular power spectrum model is computed with set_model_bias_MF_theory and can model
	 *   photometric errors and RSDs.
	 *
	 *  The parameter of the model are:
	 *
	 *  - N cosmological parameters
	 *  - \f$\beta\f$ halo bias inaccuracy
	 *  - \f$\alpha\f$ deviations from poissonian shot noise
	 *
	 *  @param cosmo_param list of cosmological parameter for the model
	 *
	 *  @param cosmo_prior prior for the cosmological parameters
	 *  
	 *  @param beta_prior prior for the beta parameter
	 *
	 *  @param alpha_prior prior for the alpha parameter
 	 *
	 */
	void set_model_covariance_bias_MF_theory (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution alpha_prior);
	
      };
    }
  }
}

#endif
