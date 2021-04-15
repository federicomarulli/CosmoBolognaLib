/********************************************************************
 *  Copyright (C) 2021 by Federico Marulli and Giorgio Lesci        *
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
 *  @file Headers/CombinedModelling.h
 *
 *  @brief The class CombinedModelling
 *
 *  This file defines the interface of the class CombinedModelling, used for
 *  jointly modelling any kind of measurements
 *
 *  @author Giorgio Lesci, Sofia Contarini (and Federico Marulli)
 *
 *  @author giorgio.lesci2@unibo.it, sofia.contarini3@unibo.it (and federico.marulli3@unibo.it)
 */

#ifndef __COMBMODELLING__
#define __COMBMODELLING__


#include "CombinedPosterior.h"
#include "Modelling.h"


// ===================================================================================================


namespace cbl {

  /**
   *  @brief The namespace of the functions and classes used for <B>
   *  modelling </B>
   *  
   *  The \e modelling namespace contains all the functions and
   *  classes used to model any kind of measurements
   */
  namespace modelling {

    /**
     *  @class CombinedModelling CombinedModelling.h "Headers/CombinedModelling.h"
     *
     *  @brief The class CombinedModelling
     *
     *  This file defines the interface of the base class CombinedModelling,
     *  used for combining any kind of modelling
     *
     */
    class CombinedModelling : public Modelling {

    protected:

      /// combined posterior
      std::shared_ptr<statistics::CombinedPosterior> m_combined_posterior;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{
      
      /**
       *  @brief default constuctor 
       */
      CombinedModelling () = default;
      
      /**
       *  @brief Constuctor used to combine statistically independent probes
       *
       *  @param modelling vector of pointers to the single Modelling objects
       *
       *  @param repeated_par parameters shared by different probes, for which
       *  the user wants different posteriors for each probe. For example, if
       *  the probes A and B depend on the same parameter \f$p\f$, whose identification
       *  string is "par", then "par" must be given in input if the user desires to have
       *  different posteriors of \f$p\f$ for the two probes A and B. Every probe depending
       *  on \f$p\f$ will provide a different posterior on \f$p\f$. This is useful
       *  in the case of astrophysical parameters, e.g. the parameters describing
       *  galaxy cluster profiles. 
       */
      CombinedModelling (std::vector<std::shared_ptr<modelling::Modelling>> modelling, std::vector<std::string> repeated_par={});
      
      /**
       *  @brief Constructor used to set the modelling of
       *  statistically dependent probes. For each vector of Modelling
       *  objects in the \f$modelling\f$ argument, a CovarianceMatrix
       *  object must be defined. The probes in each vector within
       *  \f$modelling\f$ are described by a Gaussian likelihood.
       *  The final likelihood is given by the sum of the logarithms
       *  of each Gaussian likelihood, where each Gaussian likelihood
       *  describes a set of dependent probes.
       *
       *  If additional independent probes are provided through \f$independent\_modelling\f$,
       *  then the relative likelihood functions are the ones defined by the user for
       *  each of those Modelling objects. The logarithms of these likelihoods
       *  are added to the previous Gaussian log-likelihoods.
       *
       *  @param modelling vector of vectors of pointers 
       *  to Modelling objects. In each vector a set of dependent
       *  probes is contained, with a covariance matrix
       *  defined by the corresponding CovarianceMatrix object
       *  given as input in the second argument of this constructor
       *
       *  @param covariance objects defining the covariance
       *  matrices for the Modelling objects given in input through
       *  the modelling argument. 
       *  
       *  @param independent_modelling set of Modelling objects
       *  independent from each other and from all the others Modelling
       *  objects.
       *
       *  @param repeated_par parameters shared by different probes, for which
       *  the user wants different posteriors for each probe. For example, if
       *  the probes A and B depend on the same parameter \f$p\f$, whose identification
       *  string is "par", then "par" must be given in input if the user desires to have
       *  different posteriors of \f$p\f$ for the two probes A and B. Every probe depending
       *  on \f$p\f$ will provide a different posterior on \f$p\f$. This is useful
       *  in the case of astrophysical parameters, e.g. the parameters describing
       *  galaxy cluster profiles.
       *
       */
      CombinedModelling (std::vector<std::vector<std::shared_ptr<modelling::Modelling>>> modelling, const std::vector<std::shared_ptr<data::CovarianceMatrix>> covariance, std::vector<std::shared_ptr<modelling::Modelling>> independent_modelling={}, std::vector<std::string> repeated_par={});

      /**
       *  @brief default destructor 
       */
      virtual ~CombinedModelling () = default;

      ///@}
      
      /**
       *  @name Member functions used to manage likelihood/posterior
       *  distributions
       */
      ///@{

      /**
       *  @brief function that maximizes the combined posterior, finds the
       *  best-fit parameters and stores them in the model
       *
       *  this function exploits the Nelder-Mead method
       *  https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
       *
       *  the algorithm defines a simplex (i.e a k-dimensional
       *  polytope which is the convex hull of its k+1 vertices) in
       *  the parameter space. At each step, it identifies the
       *  simplex vertex at which the function to be minimised
       *  (i.e. the negative posterior in this case) has the
       *  greatest value, and moves it, via reflections and scaling,
       *  to a new position in which the function has a lower
       *  value. This iteration stops when the simplex area becomes
       *  lower than the tolerance. For instance, in 2D, the
       *  starting vertices of the simplex (a triangle in 2D) are
       *  the following: (start[0], start[1]) ; (start[0]+epsilon,
       *  start[1]) ; (start[0], start[1]+epsilon)
       *
       *  @param start vector containing initial values for
       *  the posterior maximization
       *
       *  @param max_iter the maximum number of iterations
       *
       *  @param tol the tolerance to find convergence
       *
       *  @param epsilon the relative fraction of the initial
       *  simplex size
       *
       */
      void maximize_combined_posterior (const std::vector<double> start, const unsigned int max_iter=10000, const double tol=1.e-6, const double epsilon=1.e-3);
      
      /**
       *  @brief sample the posterior, initializing the chains by
       *  drawing from the prior distributions
       *
       *  the starting values of the chain are extracted from the
       *  (possibly different) distributions of the priors
       * 
       *  @param chain_size the chain lenght
       *
       *  @param nwalkers the number of parallel chains
       *
       *  @param aa the parameter of the \f$g(z)\f$ distribution
       *
       *  @param parallel false \f$\rightarrow\f$ non-parallel
       *  sampler; true \f$\rightarrow\f$ parallel sampler
       */
      void sample_combined_posterior (const int chain_size, const int nwalkers, const double aa=2, const bool parallel=true);
      
      /**
       *  @brief sample the posterior, initializing the chains in a
       *  ball around the posterior best-fit parameters values
       *
       *  the starting values of the chain are extracted from
       *  uniform distributions in the range [parameter-radius,
       *  parameter+radius] (for each likelihood parameter)
       *
       *  this function first maximizes the posterior, starting the
       *  computation at the values of the input vector 'start',
       *  then it inizializes the chain
       *
       *  @param chain_size the chain lenght
       *
       *  @param nwalkers the number of parallel chains
       *
       *  @param radius radius of the ball in parameter space
       *
       *  @param start std::vector containing initial values for the
       *  posterior maximization
       *
       *  @param max_iter the maximum number of iterations
       *
       *  @param tol the tolerance in finding convergence 
       *
       *  @param epsilon the simplex side
       *
       *  @param aa the parameter of the \f$g(z)\f$ distribution
       *
       *  @param parallel false \f$\rightarrow\f$ non-parallel
       *  sampler; true \f$\rightarrow\f$ parallel sampler
       */
      void sample_combined_posterior (const int chain_size, const int nwalkers, const double radius, const std::vector<double> start, const unsigned int max_iter=10000, const double tol=1.e-6, const double epsilon=1.e-3, const double aa=2, const bool parallel=true);
      
      ///@}
      
      
      /**
       *  @name Member functions used for Input/Output
       */
      ///@{      
      /**
       *  @brief write the results of the MCMC sampling to file
       *  
       *  this function stores to file the posterior mean, the
       *  posterior standard deviation, the posterior median, 18th and
       *  82th posterior percentiles, and, optionally, the posterior
       *  mode.
       *
       *  If the covariance matrix has been estimated from a set of
       *  mock catalogues, and the input parameters ns (\f$n_s\f$, the
       *  number of samples used to estimate the covariance matrix) is
       *  provided (>0), then the parameter errors (\f$\sigma_p\f$)
       *  will be corrected to take into account the uncertainities in
       *  the covariance estimate (Percival et al. 2014):
       *
       *  \f[ \sigma_p = \sqrt{\frac{1+B(n_b-n_p)}{1+A+B(n_p+1)}} \f]
       *
       *  where 
       *
       *  \f[ A = \frac{2}{(n_s-n_b-1)(n_s-n_b-4)} \,, \f]
       *
       *  \f[ B = \frac{(n_s-n_b-2)}{(n_s-n_b-1)(n_s-n_b-4)} \,, \f]
       *
       *  where \f$n_b\f$ is number of data measurements (e.g. the
       *  bins of the dataset).
       *
       *  This correction can be applied only
       *  if the likelihood is Gaussian. Morever, the inverce
       *  covariance matrix estimator has to be corrected to take into
       *  account the inverse Wishart distribution (Hartlap, Simon and
       *  Schneider 2006).
       *
       *  @param output_dir the output director
       *
       *  @param root_file the root of the output files: -
       *  file_root_parameters.dat file containing the output of the
       *  MCMC sampling for each parameter -
       *  file_root_covariance.dat file containing the covariance of
       *  the parameters - file_root_chain file containing the
       *  chains: the extention can be .dat or .fits
       *
       *  @param start the minimum chain position to be written
       *
       *  @param thin the step used for dilution on screen
       *
       *  @param nbins the number of bins
       *
       *  @param fits false \f$\rightarrow\f$ ascii file; true
       *  \f$\rightarrow\f$ fits file
       *
       *  @param compute_mode true \f$\rightarrow\f$ compute the
       *  posterior mode; false \f$\rightarrow\f$ do not compute the
       *  posterior mode
       *
       *  @param ns number of samples used to estimate the covariance
       *  matrix
       */
      void write_combined_results (const std::string output_dir, const std::string root_file, const int start=0, const int thin=1, const int nbins=50, const bool fits=false, const bool compute_mode=false, const int ns=-1);
      
      /**
       * @brief write the model computing 16th, 50th and
       * 84th percentiles from the MCMC
       *
       * @param output_dir the output directory
       *
       * @param output_file the output file
       *
       * @param start the minimum chain position to be written
       *
       * @param thin the step used for dilution on screen
       */
      void write_model_from_combined_chain (const std::string output_dir, const std::string output_file, const int start, const int thin);
      
      ///@}
      
    };
  }
}

#endif
