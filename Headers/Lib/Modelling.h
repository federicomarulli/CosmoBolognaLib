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
 *  @file Headers/Lib/Modelling.h
 *
 *  @brief The class Modelling
 *
 *  This file defines the interface of the class Modelling, used for
 *  modelling any kind of measurements
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLING__
#define __MODELLING__


#include "Likelihood.h"


// ===================================================================================================


namespace cosmobl {
  
  /**
   *  @brief The namespace of the functions and classes used for <B>
   *  modelling </B>
   *  
   * The \e modelling namespace contains all the functions and classes
   * used to model any kind of measurements
   */
  namespace modelling {
    
    /**
     *  @class Modelling Modelling.h
     *  "Headers/Lib/Modelling.h"
     *
     *  @brief The class Modelling
     *
     *  This file defines the interface of the base class Modelling,
     *  used for modelling any kind of measurements
     *
     */
    class Modelling {
      
    protected:
      
      /// input data to be modelled
      shared_ptr<data::Data> m_data;
      
      /// check if fit range has been set
      bool m_fit_range = false;

      /// input data to be modelled
      shared_ptr<data::Data> m_data_fit;

      /// input model
      shared_ptr<statistics::Model> m_model;

      /// likelihood parameters
      shared_ptr<statistics::LikelihoodParameters> m_parameters;

      /// likelihood
      shared_ptr<statistics::Likelihood> m_likelihood;

      
    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Modelling
       */
      Modelling () = default;

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Modelling () = default;

      ///@}

      
      /**
       *  @name Member functions used to get the protected members
       */
      ///@{
      
      /**
       * @brief return the dataset
       * @return pointer to the current dataset
       */
      shared_ptr<data::Data> data () { return m_data; }

      /**
       * @brief return the model
       * @return pointer to the current model
       */    
      shared_ptr<statistics::Model> model () { return m_model; }
      
      /**
       * @brief return the parameter
       *
       * @param pp the parameter index
       *
       * @return pointer to the parameters
       */    
      shared_ptr<statistics::Parameter> parameter (const int pp) { return m_parameters->parameter(pp); }

      /**
       * @brief return the parameters
       * @return pointer to the likelihood parameters
       */    
      shared_ptr<statistics::LikelihoodParameters> parameters () { return m_parameters; }

      /**
       * @brief return the likelihood
       * @return pointer to the current likelihood
       */    
      shared_ptr<statistics::Likelihood> likelihood () { return m_likelihood; }

      ///@}

      
      /**
       *  @name Member functions used to set internal parameters
       */
      ///@{
      
      /**
       *  @brief reset the fit range 
       *
       *  set m_fit_range = false, that means that the fit range is
       *  unset
       *
       *  @return none
       */
      void reset_fit_range () { m_fit_range = false; }

      /**
       *  @brief set the fit range 
       *
       *  @param xmin minimum x value used for the fit
       *
       *  @param xmax maximum x value used for the fit
       *
       *  @return none
       */
      virtual void set_fit_range (const double xmin, const double xmax)
      { (void)xmin; (void)xmax; ErrorCBL("Error in set_fit_range of Modelling.h."); }

      /**
       *  @brief set the fit range 
       *
       *  @param xmin minimum x value used for the fit
       *
       *  @param xmax maximum x value used for the fit
       *
       *  @param ymin minimum y value used for the fit
       *
       *  @param ymax maximum y value used for the fit
       *
       *  @return none
       */
      virtual void set_fit_range (const double xmin, const double xmax, const double ymin, const double ymax)
      { (void)xmin; (void)xmax; (void)ymin; (void)ymax; ErrorCBL("Error in set_fit_range of Modelling.h."); }

      /**
       * @brief set the dataset
       *
       * @param dataset the dataset 
       *
       * @return none
       */
      void set_data (const shared_ptr<data::Data> dataset) { m_data = move(dataset); }
      
      /**
       * @brief set the model
       *
       * @param model the model
       *
       * @return none
       */    
      void set_model (const shared_ptr<statistics::Model> model) { m_model = move(model); }

      /**
       * @brief set the parameters
       *
       * @param parameters the likelihood parameters
       *
       * @return none
       */    
      void set_parameters (const vector<shared_ptr<statistics::Parameter>> parameters);

      /**
       * @brief set the likelihood
       *
       * @param likelihood_type the likelihood type
       *
       * @return none
       */
      void set_likelihood (const statistics::LikelihoodType likelihood_type);

      /**
       * @brief set the likelihood
       *
       * @param likelihood_func the likelihood function
       *
       * @return none
       */
      void set_likelihood (const statistics::LogLikelihood_function likelihood_func);

      /**
       *  @brief function that maximize the likelihood, find the
       *  best-fit parameters and store them in model
       *
       *  @param guess vector containing initial guess values and
       *  result of the likelihood maximization
       *
       *  @param usePriors false &rarr; minimize \f$
       *  -\log(\mathcal{L}) \f$; true &rarr; minimize \f$
       *  -\log(\mathcal{L})- \Sum\log(\mathcal{Prior}) \f$:
       *
       *  @param max_iter the maximum number of iterations
       *
       *  @param tol the tolerance in finding convergence 
       *
       *  @param epsilon the relative fraction of the interval size
       *
       *  @return none
       */
      void maximize_likelihood (vector<double> &guess, const bool usePriors=false, const unsigned int max_iter=100, const double tol=1.e-6, const double epsilon=1.e-3); 

      /**
       *  @brief function that maximize the likelihood, find the
       *  best-fit parameters and store them in model, firstly
       *  extracting prior values to find the initial guess
       * 
       *  @param guess vector containing the results of the likelihood
       *  maximization
       *
       *  @param ntry number of prior extractions to find the initial
       *  guess
       *
       *  @param prior_seed base prior seed
       *
       *  @param usePriors false &rarr; minimize \f$
       *  -\log(\mathcal{L}) \f$; true &rarr; minimize \f$
       *  -\log(\mathcal{L})- \Sum\log(\mathcal{Prior}) \f$:
       *
       *  @param max_iter maximum number of iterations 
       *
       *  @param tol the tolerance in finding convergence
       *
       *  @param epsilon the relative fraction of the interval size
       *
       *  @return none
       */
      void maximize_likelihood (vector<double> &guess, const int ntry, const int prior_seed=431412, const bool usePriors=false, const unsigned int max_iter=100, const double tol=1.e-6, const double epsilon=1.e-3);

      /**
       * @brief sample the likelihood
       *
       * @param chain_size size the chain lenght
       *
       * @param nwalkers the number of parallel walkers
       * 
       * @param seed the base seed for initialization
       *
       * @param aa parameter for the stretch-move distribution
       *
       * @return none
       */
      void sample_likelihood (const int chain_size, const int nwalkers, const int seed, const double aa=2);
      
      /**
       * @brief sample the likelihood
       *
       * @param chain_size size the chain lenght
       *
       * @param nwalkers the number of parallel walkers
       * 
       * @param seed the base seed for initialization
       *
       * @param start the means  for
       * normal distribution random extraction
       *
       * @param radius the stardand distribution for
       * normal distribution random extraction
       *
       * @param aa parameter for the stretch-move distribution
       *
       * @return none
       */
      void sample_likelihood (const int chain_size, const int nwalkers, const int seed, vector<double> &start, const double radius, const double aa=2);

      /**
       *  @brief function that initialize chains
       *
       *  @param chain_size the chain size
       *
       *  @param chain_values the value of the chain 
       *
       *  @param seed the base seed for initialization
       *
       *  @param aa parameter for the stretch-move distribution
       *
       *  @return none
       */
      void sample_likelihood (const int chain_size, const vector<vector<double>> chain_values, const int seed, const double aa=2);

      /**
       * @brief sample the likelihood
       *
       * @param chain_size size the chain lenght
       *
       * @param nwalkers the number of parallel walkers
       * 
       * @param seed the base seed for initialization
       *
       * @param input_dir directory of input for chains 
       *
       * @param input_file file of input for chains 
       *
       * @param aa parameter for the stretch-move distribution
       *
       * @return none
       */
      void sample_likelihood (const int chain_size, const int nwalkers, const int seed, const string input_dir, const string input_file, const double aa=2);


      ///@}

      
      /**
       *  @name Member functions used to write the outputs
       */
      ///@{
      
      /**
       *  @brief show fit results on screen
       *
       *  @param start the starting position for each chain
       *  @param thin the position step
       *  @param seed the base seed for initialization
       *
       *  @return none
       */
      void show_results (const int start, const int thin=1, const int seed=434213);

      /**
       *  @brief write fit results
       *
       *  @param dir the output directory
       *  @param file the output file
       *  @param start the starting position for each chain
       *  @param thin the position step
       *  @param seed the base seed for initialization
       *  @return none
       */
      void write_results (const string dir, const string file, const int start, const int thin=1, const int seed=434213);

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
      virtual void write_model (const string dir, const string file, const vector<double> xx, const vector<double> parameter, const int start, const int thin)
      { (void)xx; (void)parameter; (void)dir; (void)file; (void)start; (void)thin; cosmobl::ErrorCBL("Error in write_model() of Modelling.h!"); }

      /**
       *  @brief compute the best-fit model from MCMC chains
       *
       *  @param xx vector of points at which the model is computed
       *  @param start the starting position for each chain
       *  @param thin the position step
       *  @param median_model the median model
       *  @param low_model the model at 16th percentile
       *  @param up_model the model at 84th percentile
       *
       *  @return none
       */
      virtual void compute_model_from_chains (const vector<double> xx, const int start, const int thin, vector<double> &median_model, vector<double> &low_model, vector<double> &up_model)
      { (void)xx; (void)median_model; (void)low_model; (void)up_model; (void)start; (void)thin; cosmobl::ErrorCBL("Error in model_from_chains of Modelling.h!"); }

      /**
       *  @brief compute and write the model 
       *
       *  @param dir the output directory
       *
       *  @param file the output file
       *
       *  @param xx vector of points at which the model is computed,
       *  first axis
       *
       *  @param yy vector of points at which the model is computed,
       *  second axis
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
      virtual void write_model (const string dir, const string file, const vector<double> xx, const vector<double> yy, const vector<double> parameter, const int start, const int thin)
      { (void)xx; (void)yy; (void)parameter; (void)dir; (void)file; (void)start; (void)thin; cosmobl::ErrorCBL("Error in write_model() of Modelling.h!"); }

      /**
       *  @brief compute the best-fit model from MCMC chains
       *
       *  @param xx vector of points at which the model is computed
       *  @param yy vector of points at which the model is computed
       *  @param start the starting position for each chain
       *  @param thin the position step
       *  @param median_model the median model
       *  @param low_model the model at 16th percentile
       *  @param up_model the model at 84th percentile
       *
       *  @return none
       */
      virtual void compute_model_from_chains (const vector<double> xx, const vector<double> yy, const int start, const int thin, vector<vector<double>> &median_model, vector<vector<double>> &low_model, vector<vector<double>> &up_model)
      { (void)xx; (void)yy; (void)start; (void)thin; (void)median_model; (void)low_model; (void)up_model; cosmobl::ErrorCBL("Error in model_from_chains of Modelling.h!"); }
      
      /**
       *  @brief read MCMC chains
       *
       *  @param dir input directory
       *  @param file input file
       *  @param nwalkers the number of parallel walkers
       *  @param skip_header number of header lines skipped
       *
       *  @return none
       */
      virtual void read_chain (const string dir, const string file, const int nwalkers, const int skip_header=0);
      
      /**
       *  @brief write MCMC chains in a FITS file
       *
       *  @param dir the output directory
       *
       *  @param file the output file
       *
       *  @param table_name the table name
       *
       *  @param burn_in the number of the first chain steps to be
       *  discarded
       *
       *  @param thin take one chain step every \e thin steps
       *
       *  @return none
       */
      void write_chain_FITS (const string dir, const string file, const string table_name, const int burn_in=0, const int thin=1)
      {
	m_likelihood->write_chain_FITS(dir, file, table_name, burn_in, thin);
      }

      /**
       *  @brief read MCMC chains in a FITS file
       *
       *  @param dir the input directory
       *
       *  @param file the input file
       *
       *  @param table_name the table name
       *
       *  @return none
       */
      virtual void read_chain_FITS (const string dir, const string file, const string table_name)
      {
	m_likelihood->read_chain_FITS(dir, file, table_name);
      }

      ///@}
      
    };
  }
}

#endif
