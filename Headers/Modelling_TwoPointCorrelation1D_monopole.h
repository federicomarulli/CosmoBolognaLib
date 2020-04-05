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
 *  @file Headers/Modelling_TwoPointCorrelation1D_monopole.h
 *
 *  @brief The class Modelling_TwoPointCorrelation1D_monopole
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation1D_monopole, used to model the monopole
 *  of two-point correlation function
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGTWOPOINTMON__
#define __MODELLINGTWOPOINTMON__


#include "Modelling_TwoPointCorrelation1D.h"
#include "ModelFunction_TwoPointCorrelation1D_monopole.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    namespace twopt {
      
      /**
       *  @class Modelling_TwoPointCorrelation1D_monopole
       *  Modelling_TwoPointCorrelation1D_monopole.h
       *  "Headers/Modelling_TwoPointCorrelation1D_monopole.h"
       *
       *  @brief The class Modelling_TwoPointCorrelation1D_monopole
       *
       *  This file defines the interface of the base class
       *  Modelling_TwoPointCorrelation1D_monopole, used for modelling the
       *  monopole of two-point correlation function
       *
       */
      class Modelling_TwoPointCorrelation1D_monopole : public Modelling_TwoPointCorrelation1D {

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class
	 *  Modelling_TwoPointCorrelation1D_monopole
	 */
	Modelling_TwoPointCorrelation1D_monopole () = default;

	/**
	 *  @brief constructor
	 *  
	 *  @param twop the two-point correlation function to model
	 *
	 *  @return object of type
	 *  Modelling_TwoPointCorrelation1D_monopole
	 */
	Modelling_TwoPointCorrelation1D_monopole (const std::shared_ptr<cbl::measure::twopt::TwoPointCorrelation> twop)
	  : Modelling_TwoPointCorrelation1D(twop) {}

	/**
	 *  @brief constructor
	 *  
	 *  @param twop_dataset the dataset containing the two-point
	 *  correlation function to model
	 *
	 *  @return object of type
	 *  Modelling_TwoPointCorrelation1D_monopole
	 */
	Modelling_TwoPointCorrelation1D_monopole (const std::shared_ptr<data::Data> twop_dataset)
	  : Modelling_TwoPointCorrelation1D(twop_dataset, cbl::measure::twopt::TwoPType::_monopole_) {}

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_TwoPointCorrelation1D_monopole () = default;

	///@}
	

	/**
	 *  @name Member functions used to set the model parameters
	 */
	///@{

	/**
	 *  @brief set the fiducial model for the dark matter
	 *  two-point correlation function and associated quantities
	 *
	 *  @return none
	 */
	void set_fiducial_xiDM ();

	/**
	 *  @brief set the fiducial model for the dark matter power
	 *  spectrum
	 *
	 *  @return none
	 */
	void set_fiducial_PkDM ();
		
	/**
	 *  @brief set the fiducial model for the variance
	 *  \f$\sigma(M)\f$ 
	 *
	 *  @return none
	 */
	void set_fiducial_sigma_data_model ();

	/**
	 *  @brief set the fiducial model for the variance
	 *  \f$\sigma(M)\f$ and its derivative \f$d\ln\sigma(M)/d\ln
	 *  M\f$
	 *
	 *  @return none
	 */
	void set_fiducial_sigma ();

	/**
	 *  @brief set a grid with effective bias values estimating
	 *  from a set of masses
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param min_par the minimum value for the parameter where
	 *  the effective bias is computed
	 *
	 *  @param max_par the maximum value for the parameter where
	 *  the effective bias is computed
	 *
	 *  @param nbins_par the number of points for the
	 *  parameter where the effective bias is computed
	 *
	 *  @param dir the directory where is the grid of effective
	 *  bias is stored
	 *
	 *  @param file_grid_bias the file where is the grid
	 *  of effective bias is stored
	 *
	 *  @return none
	 */
	void set_bias_eff_grid (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<double> min_par, const std::vector<double> max_par, const std::vector<int> nbins_par, const std::string dir, const std::string file_grid_bias);

	/**
	 *  @brief set a grid with effective bias values estimating
	 *  from a selection function and a theoretical mass function
	 *
	 *  @param file_selection_function input file with the
	 *  selection function
	 *	 
	 *  @param column vector containing the columns with {mass,
	 *  redshift, selection function}
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param min_par the minimum value for the parameter where
	 *  the effective bias is computed
	 *  
	 *  @param max_par the maximum value for the parameter where
	 *  the effective bias is computed
	 *
	 *  @param nbins_par the number of points for the
	 *  parameter where the effective bias is computed
	 *
	 *  @param dir the directory where is the grid of effective
	 *  bias is stored
	 *
	 *  @param file_grid_bias the file where is the grid
	 *  of effective bias is stored
	 *
	 *  @return none
	 */
	void set_bias_eff_grid (const std::string file_selection_function, const std::vector<int> column, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<double> min_par, const std::vector<double> max_par, const std::vector<int> nbins_par, const std::string dir, const std::string file_grid_bias);
	
	/**
	 *  @brief set the model to fit the full shape of the monopole
	 *  of the two-point correlation function
	 *
	 *  the model is the following:
	 *
	 *  \f[\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3}
	 *  f\sigma_8 \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2
	 *  \right] \cdot \xi_{\rm DM}(\alpha\cdot s)/\sigma_8^2 +
	 *  \sum_{i=0}^N \frac{A_i}{r^{i}}\f]
	 *
	 *  the model has 3+N parameters: 
	 *    - \f$\alpha\f$
	 *    - \f$f(z)\sigma_8(z)\f$
	 *    - \f$b(z)\sigma_8(z)\f$
	 *    - \f$A_i\f$
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters
	 *
	 *  @param alpha_prior prior for the parameter \f$\alpha\f$
	 *
	 *  @param fsigma8_prior prior for the parameter
	 *  \f$f(z)\sigma_8(z)\f$
	 *
	 *  @param bsigma8_prior prior for the parameter
	 *  \f$b(z)\sigma_8(z)\f$
	 *
	 *  @param polynomial_prior vector containing the priors
	 *  for the polynomial part: order of the polynomial is
	 *  the size of the vector \f$-1\f$
	 *
	 *  @return none
	 */
	void set_model_linear (const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const std::vector<statistics::PriorDistribution> polynomial_prior);

	/**
	 *  @brief set the model to fit the full shape of the monopole
	 *  of the two-point correlation function, providing in output
	 *  the positions of the BAO peak, deap and linear point
	 *
	 *  the model is the following:
	 *
	 *  \f[\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3}
	 *  f\sigma_8 \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2
	 *  \right] \cdot \xi_{\rm DM}(\alpha\cdot s)/\sigma_8^2 +
	 *  \sum_{i=0}^N \frac{A_i}{r^{i}}\f]
	 *
	 *  the model has 3+N parameters: 
	 *    - \f$\alpha\f$
	 *    - \f$f(z)\sigma_8(z)\f$
	 *    - \f$b(z)\sigma_8(z)\f$
	 *    - \f$A_i\f$
	 *
	 *  the positions of the BAO peak, dip, and linear point are
	 *  provided in output as derived parameters
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters
	 *
	 *  @param alpha_prior prior for the parameter \f$\alpha\f$
	 *
	 *  @param fsigma8_prior prior for the parameter
	 *  \f$f(z)\sigma_8(z)\f$
	 *
	 *  @param bsigma8_prior prior for the parameter
	 *  \f$b(z)\sigma_8(z)\f$
	 *
	 *  @param polynomial_prior vector containing the priors
	 *  for the polynomial part: order of the polynomial is
	 *  the size of the vector \f$-1\f$
	 *
	 *  @return none
	 */
	void set_model_linear_LinearPoint (const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const std::vector<statistics::PriorDistribution> polynomial_prior);

	/**
	 *  @brief set the model to fit the full shape of the monopole
	 *  of the two-point correlation function as a polynomial,
	 *  providing in output the positions of the BAO peak, deap
	 *  and linear point
	 *
	 *  @param polynomial_prior vector containing the priors
	 *  for the polynomial part: order of the polynomial is
	 *  the size of the vector \f$-1\f$
	 *
	 *  @return none
	 */
	void set_model_polynomial_LinearPoint (const std::vector<statistics::PriorDistribution> polynomial_prior);

	/**
	 *  @brief set the parameters to model the monopole of the
	 *  two-point correlation function in redshift space
	 * 
	 *  redshift-space distorsions are modelled in the Kaiser
	 *  limit, that is neglecting non-linearities in dynamics
	 *  and bias; specifically, the model considered is the
	 *  following:
	 *  
	 *  \f[\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3} f\sigma_8
	 *  \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2 \right] \cdot
	 *  \xi_{\rm DM}(s)/\sigma_8^2\f]
	 * 
	 *  @param fsigma8_prior prior for the parameter
	 *  \f$f(z)\sigma_8(z)\f$
	 *
	 *  @param bsigma8_prior prior for the parameter
	 *  \f$b(z)\sigma_8(z)\f$
	 *
	 *  @return none
	 */
	void set_model_Kaiser (const statistics::PriorDistribution fsigma8_prior={}, const statistics::PriorDistribution bsigma8_prior={});

	/**
	 *  @brief set the parameters to model the monopole of the
	 *  two-point correlation function in redshift space
	 * 
	 *  redshift-space distorsions are modelled in the Kaiser
	 *  limit, that is neglecting non-linearities in dynamics
	 *  and bias; specifically, the model considered is the
	 *  following:
	 *         
	 *  \f[\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3}
	 *  f(z)\sigma_8(z) \cdot b(z)\sigma_8(z) +
	 *  \frac{1}{5}(f(z)^2\sigma_8(z)^2) \right] \cdot \xi_{\rm
	 *  DM}(s)\frac{\sigma_8^2}{\sigma_8^2(z)} \f]
	 *
	 *  @param sigma8_prior prior for the parameter
	 *  \f$\sigma_8(z)\f$
	 *
	 *  @param bias_prior prior for the parameter bias
	 *  \f$b(z)\f$
	 *
	 *  @return none
	 */
	void set_model_sigma8_bias (const statistics::PriorDistribution sigma8_prior={}, const statistics::PriorDistribution bias_prior={});
	
	/**
	 *  @brief set the model to fit the full shape of the monopole
	 *  of the two-point correlation function
	 *
	 *  the model is the following:
	 *
	 *  \f[\xi_0(s) = \left[ b^2 + \frac{2}{3}f \cdot b +
	 *  \frac{1}{5}f^2 \right] \cdot \xi_{\rm
	 *  DM}\left(\frac{D_V(z)}{D_V^{fid}(z)}\cdot s\right) \f]
	 *
	 *  the model has 1+N parameters: 
	 *    - bias
	 *    - cosmological paramters
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters
	 *
	 *  @param bias_prior prior for the parameter
	 *  \f$b(z)\sigma_8(z)\f$
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_param_prior vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @return none
	 */
	void set_model_linear_bias_cosmology (const statistics::PriorDistribution bias_prior={}, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param={}, const std::vector<statistics::PriorDistribution> cosmo_param_prior={});
	
        /**
	 *  @brief set the model to fit the full shape of the monopole
	 *  of the two-point correlation function with cluster masses
	 *  provided in input, with only \f$sigma_8\f$ as a free
	 *  parameter
	 *
	 *  the model is the following:
	 *
	 *  \f[\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3} f\sigma_8
	 *  \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2 \right] \cdot
	 *  \xi_{\rm DM}(s)/\sigma_8^2\f]
	 *
	 *  the model has 1 parameter: \f$sigma_8\f$
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters; the linear
	 *  effective bias is computed for each cosmology using the
	 *  provided halo masses by cbl::cosmology::bias_eff
	 *	 
	 *  @param sigma8_prior prior for the parameter
	 *  \f$\sigma_8(z)\f$
	 *
	 *  @return none
	 */
	void set_model_linear_sigma8_clusters (const statistics::PriorDistribution sigma8_prior={});

	/**
	 *  @brief set the model to fit the full shape of the monopole
	 *  of the two-point correlation function with either the
	 *  cluster masses provided in input, or using a selection
	 *  function to estimate the bias
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters; the linear
	 *  effective bias is computed for each cosmology using the
	 *  provided halo masses by cbl::cosmology::bias_eff
	 *
	 *  @param cosmo_param the cosmological
	 *  parameter
	 *
	 *  @param cosmo_param_prior the prior for
	 *  the cosmological parameter 
	 *
	 *  @param dir the directory where is the grid
	 *  of effective bias is stored
	 *
	 *  @param file_grid_bias the file where is the grid
	 *  of effective bias is stored
	 *
	 *  @param min_par the minimum value for the
	 *  parameter where the effective bias is computed
	 *  
	 *  @param max_par the maximum value for the
	 *  parameter where the effective bias is computed
	 *
	 *  @param nbins_par the number of points for the
	 *  parameter where the effective bias is computed
	 *
	 *  @param file_selection_function input file with the
	 *  selection function
	 *
	 *  @param column vector containing the columns with {mass,
	 *  redshift, selection function}
	 *
	 *  @return none
	 */
	void set_model_linear_cosmology_clusters_grid (const cbl::cosmology::CosmologicalParameter cosmo_param, const statistics::PriorDistribution cosmo_param_prior, const std::string dir, const std::string file_grid_bias, const double min_par, const double max_par, const int nbins_par, const std::string file_selection_function=par::defaultString, const std::vector<int> column={0, 1, 2});
		
	/**
	 *  @brief set the model to fit the full shape of the monopole
	 *  of the two-point correlation function with either the
	 *  cluster masses provided in input, or using a selection
	 *  function to estimate the bias
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters; the linear
	 *  effective bias is computed for each cosmology using the
	 *  provided halo masses by cbl::cosmology::bias_eff
	 *
	 *  @param cosmo_param1 the first cosmological
	 *  parameter
	 *
	 *  @param cosmo_param_prior1 the prior for
	 *  the first cosmological parameter 
	 *
	 *  @param cosmo_param2 the second cosmological
	 *  parameter
	 *
	 *  @param cosmo_param_prior2 the prior for
	 *  the second cosmological parameter 
	 *
	 *  @param dir the directory where is the grid
	 *  of effective bias is stored
	 *
	 *  @param file_grid_bias the file where is the grid
	 *  of effective bias is stored
	 *
	 *  @param min_par1 the minimum value for the first
	 *  parameter where the effective bias is computed
	 *  
	 *  @param max_par1 the maximum value for the first
	 *  parameter where the effective bias is computed
	 *
	 *  @param nbins_par1 the number of points for the first
	 *  parameter where the effective bias is computed
	 *
	 *  @param min_par2 the minimum value for the second
	 *  parameter where the effective bias is computed
	 *  
	 *  @param max_par2 the maximum value for the second
	 *  parameter where the effective bias is computed
	 *
	 *  @param nbins_par2 the number of points for the second
	 *  parameter where the effective bias is computed
	 *	
	 *  @param file_selection_function input file with the
	 *  selection function
	 *	 
	 *  @param column vector containing the columns with {mass,
	 *  redshift, selection function}
	 *
	 *  @return none
	 */
	void set_model_linear_cosmology_clusters_grid (const cbl::cosmology::CosmologicalParameter cosmo_param1, const statistics::PriorDistribution cosmo_param_prior1, const cbl::cosmology::CosmologicalParameter cosmo_param2, const statistics::PriorDistribution cosmo_param_prior2, const std::string dir, const std::string file_grid_bias, const double min_par1, const double max_par1, const int nbins_par1, const double min_par2, const double max_par2, const int nbins_par2, const std::string file_selection_function=par::defaultString, const std::vector<int> column={0, 1, 2});
	
	/**
	 *  @brief set the model to fit the full shape of the monopole
	 *  of the two-point correlation function in redshift space,
	 *  using a given selection function to estimate the bias
	 * 
	 *  the monopole of the two-point correlation function in
	 *  redshift space is computed by
	 *  cbl::modelling::twopt::xi0_linear_cosmology_clusters_selection_function
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters; the linear
	 *  effective bias is a derived parameter, that is computed
	 *  for each cosmology, using the provided halo masses, by
	 *  cbl::cosmology::bias_eff
	 *
	 *  @param alpha_prior prior for the \f$\alpha\f$ parameter of
	 *  the cluster mass scaling relation
	 *
	 *  @param cosmo_param the model cosmological parameters
	 *
	 *  @param cosmo_param_prior the prior for the model
	 *  cosmological parameters
	 *
	 *  @return none
	 */
	void set_model_linear_cosmology_cluster_selection_function (const statistics::PriorDistribution alpha_prior, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_param_prior);

	
	/**
	 *  @brief set the model to fit the full shape of the monopole
	 *  of the two-point correlation function with cluster masses
	 *  provided in input
	 *
	 *  the model is the following:
	 *
	 *  \f[\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3} f\sigma_8
	 *  \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2 \right] \cdot
	 *  \xi_{\rm DM}(s)/\sigma_8^2\f]
	 *
	 *  the model has N cosmological parameters
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters; the linear
	 *  effective bias is computed for each cosmology using the
	 *  provided halo masses by cbl::cosmology::bias_eff
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_param_prior vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @return none
	 */
	void set_model_linear_cosmology_clusters (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param={}, const std::vector<statistics::PriorDistribution> cosmo_param_prior={});
	
	/**
	 *  @brief set the function to model the monopole of the
	 *  two-point correlation function, taking into accout
	 *  geometric distortions (that is the Alcock-Paczynski
	 *  effect), and using a second order polynomial. The
	 *  redshift-space distortions are accounted by the polynomial
	 *  factors, and not used to constrain the linear growth rate
	 *  parameter.
	 *
	 *  The model is the following:
	 *
	 *  \f[ \xi(s)= \left(b\sigma_8\right)^2 \frac{\xi_{DM}(\alpha
	 *  r)}{\sigma_8}\ + A_0 + \frac{A_1}{r} + \frac{A_2}{r^2} \f]
	 *
	 *  where the dark matter two-point correlation function,
	 *  \f$\xi_{DM}\f$, is computed at the fiducial (fixed)
	 *  cosmology, and {\f$A_0\f$, \f$A_1\f$, \f$A_2\f$} are
	 *  considered as nuisance parameters
	 *
	 *  @param alpha_prior prior for the parameter \f$\alpha\f$
	 *
	 *  @param bs8_prior prior for the parameter
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
	void set_model_BAO (const statistics::PriorDistribution alpha_prior={}, const statistics::PriorDistribution bs8_prior={}, const statistics::PriorDistribution A0_prior={}, const statistics::PriorDistribution A1_prior={}, const statistics::PriorDistribution A2_prior={});

	/**
	 *  @brief set the parameter to model the monopole of the
	 *  two-point correlation function in real space, taking into
	 *  accout geometric distortions (that is the Alcock-Paczynski
	 *  effect), and using a second order polynomial
	 *
	 *  the model used is the following:
	 *
	 *  \f[\xi(s)= b^2 \xi_{DM}(\alpha r, \Sigma_{NL})\ + A_0 + A_1/r
	 *  +A_2/r^2\f]
	 *
	 *  where \f$\xi_{DM}\f$ is computed at the fiducial (fixed)
	 *  cosmology, with damping of the BAO peak
	 *
	 *  @param sigmaNL_prior prior for the parameter \f$\Sigma_{NL}\f$
	 *
	 *  @param alpha_prior prior for the parameter \f$\alpha\f$
	 *
	 *  @param BB_prior prior for the parameter
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
	void set_model_BAO_sigmaNL (const statistics::PriorDistribution sigmaNL_prior={}, const statistics::PriorDistribution alpha_prior={}, const statistics::PriorDistribution BB_prior={}, const statistics::PriorDistribution A0_prior={}, const statistics::PriorDistribution A1_prior={}, const statistics::PriorDistribution A2_prior={});

	/**
	 *  @brief set the parameter to model the monopole of the
	 *  two-point correlation function in real space, taking into
	 *  accout geometric distortions (that is the Alcock-Paczynski
	 *  effect), and using a second order polynomial, providing in
	 *  output the positions of the BAO peak, deap and linear
	 *  point
	 *
	 *  the model used is the following:
	 *
	 *  \f[\xi(s)= b^2 \xi_{DM}(\alpha r)\ + A_0 + A_1/r
	 *  +A_2/r^2\f]
	 *
	 *  the positions of the BAO peak, deap and linear point are
	 *  provided in output as derived parameters
	 *
	 *  where \f$\xi_{DM}\f$ is computed at the fiducial (fixed)
	 *  cosmology, and {\f$b\sigma_8\f$, \f$A_0\f$, \f$A_1\f$,
	 *  \f$A_2\f$} are considered as nuisance parameters
	 *
	 *  @param alpha_prior prior for the parameter \f$\alpha\f$
	 *
	 *  @param BB_prior prior for the parameter
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
	void set_model_BAO_LinearPoint (const statistics::PriorDistribution alpha_prior={}, const statistics::PriorDistribution BB_prior={}, const statistics::PriorDistribution A0_prior={}, const statistics::PriorDistribution A1_prior={}, const statistics::PriorDistribution A2_prior={});

	/**
	 *  @brief set the HOD parameters used to model the full shape
	 *  of the monopole of the two-point correlation function
	 *
	 *  the HOD model for \f$\xi(r)\f$ is the one implemented by
	 *  cbl::modelling::twopt::xi_HOD
	 *
	 *  the model has 5 free parameters:
	 *
	 *  - \f$M_{min}\f$: the mass scale at which 50% of haloes
	 *   host a satellite galaxy
	 *
	 *  - \f$\sigma_{\log M_h}\f$: transition width reflecting
	 *   the scatter in the luminosity-halo mass relation
	 *
	 *  - \f$M_0\f$: the cutoff mass
	 *
	 *  - \f$M_1\f$: the amplitude of the power law
	 *
	 *  - \f$\alpha\f$: the slope of the power law
	 *
	 *  @param Mmin_prior \f$M_{min}\f$ prior
	 *
	 *  @param sigmalgM_prior \f$\sigma_{\log M_h}\f$ prior
	 *
	 *  @param M0_prior \f$M_0\f$ prior
	 *
	 *  @param M1_prior \f$\alpha\f$ prior
	 *
	 *  @param alpha_prior \f$\alpha\f$ prior
	 *
	 *  @return none
	 */
	void set_model_HOD (const statistics::PriorDistribution Mmin_prior={}, const statistics::PriorDistribution sigmalgM_prior={}, const statistics::PriorDistribution M0_prior={}, const statistics::PriorDistribution M1_prior={}, const statistics::PriorDistribution alpha_prior={});

	/**
	 *  @brief set the parameters to model the monopole of the
	 *  two-point correlation function in redshift space
	 * 
	 *  redshift-space distorsions are modelled in the Kaiser
	 *  limit, that is neglecting non-linearities in dynamics
	 *  and bias; specifically, the model considered is the
	 *  following:
	 * 
	 *  \f$\xi(s) = b^2 \xi'(r) + b \xi''(r) + \xi'''(r) \, ;\f$
	 *
	 *  where b is the linear bias and the terms \f$\xi'(r)\f$,
	 *  \f$\xi''(r)\f$, \f$\xi'''(r)\f$ are the Fourier
	 *  anti-transform of the power spectrum terms obtained
	 *  integrating the redshift space 2D power spectrum along
	 *  \f$\mu\f$ (see cbl::modelling::twopt.:damped_Pk_terms)
	 *
	 *  @param bias_prior prior for the parameter bias
	 *  \f$b(z)\f$
	 *
	 *  @param sigmaz_prior prior for the parameter
	 *  \f$\sigma_z(z)\f$
	 *
	 *  @return none
	 */
	void set_model_bias_sigmaz (const statistics::PriorDistribution bias_prior={}, const statistics::PriorDistribution sigmaz_prior={});

	/**
	 *  @brief set the parameters to model the monopole of the
	 *  two-point correlation function in redshift space
	 * 
	 *  redshift-space distorsions are modelled in the Kaiser
	 *  limit, that is neglecting non-linearities in dynamics
	 *  and bias; specifically, the model considered is the
	 *  following:
	 * 
	 *  \f$\xi(s) = b^2 \xi'(r) + b \xi''(r) + \xi'''(r) \, ;\f$
	 *
	 *  where b is the linear bias and the terms \f$\xi'(r)\f$,
	 *  \f$\xi''(r)\f$, \f$\xi'''(r)\f$ are the Fourier
	 *  anti-transform of the power spectrum terms obtained
	 *  integrating the redshift space 2D power spectrum along
	 *  \f$\mu\f$ (see cbl::modelling::twopt.:damped_Pk_terms)
	 *
	 *  @param M0_prior prior for the parameter \f$M_0\f$, the
	 *  intercept of the scaling relation
	 *
	 *  @param slope_prior prior for the slope of the scaling
	 *  relation
	 *
	 *  @param scatter_prior for the scatter of the scaling
	 *  relatino
	 *
	 *  @param sigmaz_prior prior for the parameter
	 *  \f$\sigma_z(z)\f$
	 *
	 *  @return none
	 */
	void set_model_scaling_relation_sigmaz (const statistics::PriorDistribution M0_prior={}, const statistics::PriorDistribution slope_prior={}, const statistics::PriorDistribution scatter_prior={}, const statistics::PriorDistribution sigmaz_prior={});

	///@}

      };
    }
  }
}

#endif
