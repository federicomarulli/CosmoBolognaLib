/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli                         *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file Headers/Lib/ModelFunction.h
 *
 *  @brief Functions to model data 
 *
 *  This file contains the prototypes of a set of functions to model
 *  data
 *  
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __MODFUNC__
#define __MODFUNC__


// ============================================================================

#include "Cosmology.h"


namespace cosmobl {

  namespace glob{

    /**
     * @struct STR_params
     * @brief the struct STR_params
     *
     * This struct contains the data
     * and the model for the &chi;&sup2; analysis
     */
    struct STR_twop_model{

      /// cosmology
      shared_ptr<Cosmology> cosmology;

      /// cosmological parameters
      vector<CosmoPar> Cpar;

      /// redshift
      double redshift;

      /// growth rate
      double f;

      /// method to compute the dark matter power spectrum
      string method;

      string output_root;

      bool NL;
      int norm;
      double k_min;
      double k_max;
      double aa;
      bool GSL;
      double prec;
      string file_par;
      double r_min;
      double r_max;

      shared_ptr<classfunc::func_grid_GSL> func_xi;

      double pi_max;

      STR_twop_model() {}

    };

    /**
     * @brief model for the projected correlation function.
     * Model parameter is the bias, the dark matter clustering 
     * is computed for a given cosmology:
     * \f$w_p(r_p)= b^2 w_p^{DM}(r)\f$.
     *
     * @param r the scale at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the value of the two point correlation function model
     */
    double wp_bias(double r, shared_ptr<void> parameters, vector<double> model_parameters); 

    /**
     * @brief model for the projected correlation function.
     * Model parameters are the bias and the cosmological parameters that enters in the 
     * theoretical clustering model for the dark matter.
     * \f$w_p(r_p)= b^2 w_p^{DM}(r)\f$.
     *
     * @param r the scale at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the value of the two point correlation function model
     */
    double wp_bias_cosmology(double r, shared_ptr<void> parameters, vector<double> model_parameters); 

    /**
     * @brief model for the correlation function in real space.
     * Model parameter is the bias, the dark matter clustering 
     * is computed for a given cosmology:
     * \f$\xi(r)= b^2 \xi_{DM}(r)\f$.
     *
     * @param r the scale at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the value of the two point correlation function model
     */
    double xi_bias(double r, shared_ptr<void> parameters, vector<double> model_parameters); 

    /**
     * @brief model for the correlation function in real space.
     * Model parameters are the bias and the cosmological parameters
     * that enters in the theoretical clustering model for the dark matter.
     * \f$\xi(r)= b^2 \xi_{DM}(r)\f$.
     *
     * @param r the scale at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the value of the two point correlation function model
     */
    double xi_bias_cosmology(double r, shared_ptr<void> parameters, vector<double> model_parameters); 

    /**
     * @brief model for the monopole of the correlation 
     * function in redshift space. Model parameter is the bias,
     * the dark matter clustering is computed for a given cosmology:
     * \f$\xi(s)= b^2 (1+2*\beta/3+\beta^2/5) \xi_{DM}(s)\f$.
     *
     * @param r the scale at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the value of the two point correlation function model
     */
    double xi0_bias(double r, shared_ptr<void> parameters, vector<double> model_parameters); 

    /**
     * @brief model for the monopole of the correlation 
     * function in redshift space. 
     * Model parameters are the bias and the cosmological parameters
     * that enters in the theoretical clustering model for the dark matter.
     * \f$\xi(s)= b^2 (1+2*\beta/3+\beta^2/5) \xi_{DM}(s)\f$.
     *
     * @param r the scale at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the value of the two point correlation function model
     */
    double xi0_bias_cosmology(double r, shared_ptr<void> parameters, vector<double> model_parameters); 

    /**
     * @brief model for the monopole of the correlation function. 
     * Model parameters are the bias and the &alpha; &rarr; the shift of
     * the measured two point correlation function relative
     * to the clustering model for the fiducial cosmology.
     * For BAO measure only. (see Anderson et al. 2015, and references therein)
     * The dark matter clustering is computed at the fiducial cosmology:
     * \f$\xi(s)= B^2  \xi_{DM}(\alpha s)\f$.
     *
     * @param r the scale at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the value of the two point correlation function model
     */
    double xi_alpha_B(double r, shared_ptr<void> parameters, vector<double> model_parameters);

    /**
     * @brief model for the monopole of the correlation function. 
     * Model parameters are the bias, &alpha; &rarr; the shift of
     * the measured two point correlation function relative
     * to the clustering model for the fiducial cosmology and coefficient
     * of a polynomial that model systematics in the data. 
     * For BAO measure only. (see Anderson et al. 2015, and references therein)
     * The dark matter clustering is computed at the fiducial cosmology:
     * \f$\xi(s)= B^2  \xi_{DM}(\alpha s)\ + A_0 + A_1/s +A_2/s^2\f$.
     *
     * @param r the scale at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the value of the two point correlation function model
     */
    double xi_alpha_B_poly(double r, shared_ptr<void> parameters, vector<double> model_parameters);

    /**
     * @brief model for the monopole of the correlation 
     * function in redshift space. 
     * Model parameters are the bias and &alpha; &rarr; the shift of
     * the measured two point correlation function relative
     * to the clustering model for the fiducial cosmology.
     * \f$\xi(s)= b^2 (1+2*\beta/3+\beta^2/5) \xi_{DM}(\alpha s)\f$.
     *
     * @param r the scale at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the value of the two point correlation function model
     */
    double xi0_alpha_bias(double r, shared_ptr<void> parameters, vector<double> model_parameters);

    /**
     * @brief model for the monopole of the correlation 
     * function in redshift space. 
     * Model parameters are the bias, &alpha; &rarr; the shift of
     * the measured two point correlation function relative
     * to the clustering model for the fiducial cosmology 
     * and the cosmological parameters that enters in the 
     * theoretical clustering model for the dark matter.
     * \f$\xi(s)= b^2 (1+2*\beta/3+\beta^2/5) \xi_{DM}(\alpha s)\f$.
     *
     * @param r the scale at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the value of the two point correlation function model
     */
    double xi0_alpha_bias_cosmology(double r, shared_ptr<void> parameters, vector<double> model_parameters);

    /**
     * @brief model for the projected correlation function.
     * Model parameter is the bias, the dark matter clustering 
     * is computed for a given cosmology:
     * \f$w_p(r_p)= b^2 w_p^{DM}(r)\f$.
     *
     * @param r the scales at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the values of the two point correlation function model
     */
    vector<double> wp_bias_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters); 

    /**
     * @brief model for the projected correlation function.
     * Model parameters are the bias and the cosmological parameters that enters in the 
     * theoretical clustering model for the dark matter.
     * \f$w_p(r_p)= b^2 w_p^{DM}(r)\f$.
     *
     * @param r the scales at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the values of the two point correlation function model
     */
    vector<double> wp_bias_cosmology_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters); 

    /**
     * @brief model for the correlation function in real space.
     * Model parameter is the bias, the dark matter clustering 
     * is computed for a given cosmology:
     * \f$\xi(r)= b^2 \xi_{DM}(r)\f$.
     *
     * @param r the scales at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the values of the two point correlation function model
     */
    vector<double> xi_bias_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters); 

    /**
     * @brief model for the correlation function in real space.
     * Model parameters are the bias and the cosmological parameters
     * that enters in the theoretical clustering model for the dark matter.
     * \f$\xi(r)= b^2 \xi_{DM}(r)\f$.
     *
     * @param r the scales at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the values of the two point correlation function model
     */
    vector<double> xi_bias_cosmology_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters); 

    /**
     * @brief model for the monopole of the correlation 
     * function in redshift space. Model parameter is the bias,
     * the dark matter clustering is computed for a given cosmology:
     * \f$\xi(s)= b^2 (1+2*\beta/3+\beta^2/5) \xi_{DM}(s)\f$.
     *
     * @param r the scales at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the values of the two point correlation function model
     */
    vector<double> xi0_bias_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters); 

    /**
     * @brief model for the monopole of the correlation 
     * function in redshift space. 
     * Model parameters are the bias and the cosmological parameters
     * that enters in the theoretical clustering model for the dark matter.
     * \f$\xi(s)= b^2 (1+2*\beta/3+\beta^2/5) \xi_{DM}(s)\f$.
     *
     * @param r the scales at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the values of the two point correlation function model
     */
    vector<double> xi0_bias_cosmology_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters); 

    /**
     * @brief model for the monopole of the correlation function. 
     * Model parameters are the bias and the &alpha; &rarr; the shift of
     * the measured two point correlation function relative
     * to the clustering model for the fiducial cosmology.
     * For BAO measure only. (see Anderson et al. 2015, and references therein)
     * The dark matter clustering is computed at the fiducial cosmology:
     * \f$\xi(s)= B^2  \xi_{DM}(\alpha s)\f$.
     *
     *
     * @param r the scales at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the values of the two point correlation function model
     */
    vector<double> xi_alpha_B_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters);

    /**
     * @brief model for the monopole of the correlation function. 
     * Model parameters are the bias, &alpha; &rarr; the shift of
     * the measured two point correlation function relative
     * to the clustering model for the fiducial cosmology and coefficient
     * of a polynomial that model systematics in the data. 
     * For BAO measure only. (see Anderson et al. 2015, and references therein)
     * The dark matter clustering is computed at the fiducial cosmology:
     * \f$\xi(s)= B^2  \xi_{DM}(\alpha s)\ + A_0 + A_1/s +A_2/s^2\f$.
     *
     * @param r the scales at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the values of the two point correlation function model
     */
    vector<double> xi_alpha_B_poly_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters);

    /**
     * @brief model for the monopole of the correlation 
     * function in redshift space. 
     * Model parameters are the bias and &alpha; &rarr; the shift of
     * the measured two point correlation function relative
     * to the clustering model for the fiducial cosmology.
     * \f$\xi(s)= b^2 (1+2*\beta/3+\beta^2/5) \xi_{DM}(\alpha s)\f$.
     *
     * @param r the scales at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the values of the two point correlation function model
     */
    vector<double> xi0_alpha_bias_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters);

    /**
     * @brief model for the monopole of the correlation 
     * function in redshift space. 
     * Model parameters are the bias, &alpha; &rarr; the shift of
     * the measured two point correlation function relative
     * to the clustering model for the fiducial cosmology 
     * and the cosmological parameters that enters in the 
     * theoretical clustering model for the dark matter.
     * \f$\xi(s)= b^2 (1+2*\beta/3+\beta^2/5) \xi_{DM}(\alpha s)\f$.
     *
     * @param r the scales at which the model is computed
     * @param parameters pointer to an object of type STR_twop_model
     * @param model_parameters free parameters of the model
     *
     * @return the values of the two point correlation function model
     */
    vector<double> xi0_alpha_bias_cosmology_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters);

  }
}

#endif
