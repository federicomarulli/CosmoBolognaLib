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
 *  @file Headers/Lib/ModelFunction.h
 *
 *  @brief Functions to model data 
 *
 *  This file contains all the prototypes of the functions used to
 *  model any kind of data
 *  
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODFUNC__
#define __MODFUNC__


// ============================================================================


#include "Cosmology.h"


namespace cosmobl {

  namespace modelling {

    /**
     *  @struct STR_twop_model
     *  @brief the structure STR_twop_model
     *
     *  This structure contains the data used for statistical analyses
     *  of the two-point correlation function
     */
    struct STR_twop_model {

      /// cosmology
      shared_ptr<cosmology::Cosmology> cosmology;

      /// redshift
      double redshift;

      /// method to compute the dark matter power spectrum
      string method_Pk;
      
      /// output root of the parameter file used to compute the dark matter power spectrum
      string output_root;

      /// 0 &rarr; linear power spectrum; 1 &rarr; non-linear power spectrum
      bool NL;

      /// sigmaNL damping of the wiggles in the linear power spectrum
      double sigmaNL;

      /// 0 &rarr; don't normalize the power spectrum; 1 &rarr; normalize the power spectrum
      int norm;

      /// minimum wave vector module up to which the power spectrum is computed
      double k_min;

      /// maximum wave vector module up to which the power spectrum is computed
      double k_max;

      /// parameter \e a of Eq. 24 of Anderson et al. 2012
      double aa;

      /// 0 &rarr; FFTlog is used; 1 &rarr; the GSL libraries are used
      bool GSL;

      /// accuracy of the GSL integration
      double prec;

      /// name of the parameter file
      string file_par;

      /// scales at wich the fiducal dark matter two-point correlation function is computed
      vector<double> fiducial_radDM;

      /// values of the fiducial dark matter two-point correlation function
      vector<double> fiducial_xiDM;

      /// pointer to a function of func_grid_GSL class, used to interpolate of the two-point correlation function
      shared_ptr<glob::FuncGrid> func_xi;

      /// upper limit of integration for the projected correlation function
      double pi_max;
      
      /// minimum separation up to which the correlation function is computed
      double r_min;

      /// maximum separation up to which the correlation function is computed
      double r_max;

      /// the linear growth rate at redshift z
      double linear_growth_rate_z;

      /// &sigma<SUB>8</SUB>; at redshift z
      double sigma8_z;

      /// (1+z)/HH(z)
      double var;

      /// cosmological parameters
      vector<cosmology::CosmoPar> Cpar;

      /// FV 0 &rarr; exponential form for f(v); 1 &rarr; Gaussian form for f(v); where f(v) is the velocity distribution function
      int FV;

      /// vector of &xi;(r), 
      vector<double> Xi;

      /// vector of barred &xi;(r), 
      vector<double> Xi_;

      /// vector of double-barred &xi;(r), 
      vector<double> Xi__;

      /// &xi;(r) as pointer to an interpolation function
      shared_ptr<glob::FuncGrid> funcXiR;

      /// barred &xi;(r) as pointer to an interpolation function
      shared_ptr<glob::FuncGrid> funcXiR_;

      /// double-barred &xi;(r) as pointer to an interpolation function
      shared_ptr<glob::FuncGrid> funcXiR__;

      /// 0 &rarr; linear bias; 1 &rarr; non-linear bias 
      int bias_nl;
       
      /// non-linear bias parameter
      double bA;
       
      /// 0 &rarr; standard; 1 &rarr; Chuang & Wang model
      int xiType;
       
      /// k<SUB>*</SUB> of the Chuang & Wang model
      double k_star;
       
      /// 0 &rarr; linear two-point correlation function; 1 &rarr; non-linear two-point correlation function
      int xiNL;
       
      /// v_min minimum velocity used in the convolution of the two-point correlation function
      double v_min;
       
      /// v_max maximum velocity used in the convolution of the two-point correlation function
      double v_max;
       
      /// number of steps used in the convolution of the two-point correlation function
      int step_v;
      
      /// index for pre-computed two-point correlation function
      int xi_real_index;

      /**
       * @brief default constructor
       * @return object of type STR_twop_model
       */
      STR_twop_model () = default;
    };

    
    /**
     *  @brief model for the monopole of the two-point correlation
     *  function
     *
     *  the function computes:
     *
     *  \f$\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3} f\sigma_8
     *  \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2 \right] \cdot
     *  \xi_{\rm DM}(\alpha\cdot s)/\sigma_8^2 + A_0 + A_1/s +
     *  A_2/s^2\f$
     *
     *  the model has 6 parameters: 
     *    - \f$\alpha\f$
     *    - \f$f(z)\sigma_8(z)\f$
     *    - \f$b(z)\sigma_8(z)\f$
     *    - \f$A_0\f$
     *    - \f$A_1\f$
     *    - \f$A_2\f$ 
     *
     *  the dark matter two-point correlation function is fixed and
     *  provided in input
     *
     *  @param rad the scale at which the model is computed
     *
     *  @param inputs pointer to the structure that contains the dark
     *  matter two-point correlation function and \f$sigma_8(z)\f$,
     *  computed at a given (fixed) cosmology
     *
     *  @param parameter 6D vector containing the input parameters
     *
     *  @return the monopole of the two-point correlation function
     */
    double xi0_linear (const double rad, const shared_ptr<void> inputs, vector<double> parameter); 

    
    /**
     *  @brief model for the monopole of the two-point correlation
     *  function
     *
     *  the function computes:
     *
     *  \f$\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3} f\sigma_8
     *  \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2 \right] \cdot
     *  \xi_{\rm DM}(\alpha\cdot s)/\sigma_8^2 + A_0 + A_1/s +
     *  A_2/s^2\f$
     *
     *  the model has 6 parameters: 
     *    - \f$\alpha\f$
     *    - \f$f(z)\sigma_8(z)\f$
     *    - \f$b(z)\sigma_8(z)\f$
     *    - \f$A_0\f$
     *    - \f$A_1\f$
     *    - \f$A_2\f$ 
     *
     *  the dark matter two-point correlation function is computed
     *  using the input cosmological parameters
     *
     *  @param rad the scale at which the model is computed
     *
     *  @param inputs pointer to the structure that contains the
     *  cosmological paramters used to compute the dark matter
     *  two-point correlation function
     *
     *  @param parameter 1D vector containing the linear bias
     *
     *  @return the monopole of the two-point correlation function
     */
    double xi0_linear_cosmology (const double rad, const shared_ptr<void> inputs, vector<double> parameter); 

    
    /**
     *  @brief model for the 2D two-point correlation function, in
     *  Cartesian coordinates
     *
     *  the function computes \f$\xi(r_p,\pi)\f$ with the dispersion
     *  model (see e.g. http://arxiv.org/abs/1203.1002)
     *
     *  the model has 5 parameters: 
     *    - the two Alcock-Paczynski parameters
     *      (\f$\alpha_\perp=\frac{D_{\rm A,1}(z)}{D_{\rm A,2}(z)}\f$
     *      and \f$\alpha_\parallel=\frac{H_2(z)}{H_1(z)}\f$)
     *    - \f$f(z)\sigma_8(z)\f$
     *    - \f$b(z)\sigma_8(z)\f$ 
     *    - \f$\sigma_{12}(z)\f$
     *
     *  @param rp the scale perpendicular to the line of sight at
     *  which the model is computed
     *
     *  @param pi the scale parallel to the line of sight at which the
     *  model is computed
     *
     *  @param inputs pointer to the structure that contains all the
     *  data required to implement the dispersion model
     *
     *  @param parameter 5D vector containing the input parameters
     *
     *  @return the 2D two-point correlation function in redshift
     *  space
     */
    double xi2D_dispersionModel (const double rp, const double pi, const shared_ptr<void> inputs, vector<double> parameter); 


    /**
     *  @brief Halo Occupation Distribution model for the projected
     *  two-point correlation function function
     *
     *  the function computes:
     *
     *  \f$\w_p(r_p) = ...\f$
     *
     *  the model has * parameters: 
     *    - \f$\M_{min}\f$
     *    - \f$sigma_log_M\f$
     *    - \f$\alpha\f$
     *    - \f$M_1\f$
     *    - \f$M_2\f$
     *
     *  the dark matter two-point correlation function is fixed and
     *  provided in input
     *
     *  @param rp the scale perpendicular to the line of sight at
     *  which the model is computed
     *
     *  @param inputs pointer to the structure that contains the dark
     *  matter two-point correlation function 
     *
     *  @param parameter *D vector containing the input parameters
     *
     *  @return the model projected two-point correlation function
     */
    double HOD (const double rp, const shared_ptr<void> inputs, vector<double> parameter); 
  }
}

#endif
