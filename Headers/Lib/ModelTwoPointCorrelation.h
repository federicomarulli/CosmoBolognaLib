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
 *  @file Headers/Lib/ModelTwoPointCorrelation.h
 *
 *  @brief The class ModelTwoPointCorrelation
 *
 *  This file defines the interface of the class ModelTwoPointCorrelation, 
 *  used to model the two-point correlation function
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __MTWOPOINT__
#define __MTWOPOINT__ 

#include "TwoPointCorrelation.h"


// ============================================================================================


namespace cosmobl {

  /**
   *  @class ModelTwoPointCorrelation ModelTwoPointCorrelation.h
   * "Headers/Lib/ModelTwoPointCorrelation.h"
   *
   *  @brief The class ModelTwoPointCorrelation
   *
   *  This class is used to handle objects of type <EM>
   *  ModelTwoPointCorrelation </EM>. It is used to model the
   *  two-point correlation function and extract cosmological
   *  constraints.
   */
  class ModelTwoPointCorrelation {

  private:
  
    // two-point correlation functions to be modelled
    shared_ptr<TwoPointCorrelation> m_TwoP;

    // constrained parameters 
    double m_beta_best, m_error_beta, m_sigma12_best, m_error_sigma12, m_bias_best, m_bias_sigma8_best, m_f_sigma8_best, m_bA_best, m_E_min;

    // limits of the fitting region
    vector<int> m_lim_index_fit;
    vector<double> m_lim_fit;


  public:

    // ==========================================================
    // ============== constructors and destructors ==============
    // ==========================================================

    ModelTwoPointCorrelation () 
      : m_beta_best(0.), m_error_beta(0.), m_sigma12_best(0.), m_error_sigma12(0.), m_bias_best(0.), m_E_min(0.) {} 
  
    ModelTwoPointCorrelation (shared_ptr<TwoPointCorrelation> TwoP)
      : m_TwoP(move(TwoP)) {}
      
    ~ModelTwoPointCorrelation () {}
    

    // ====================================================================
    // ============== methods to set the internal parameters ==============
    // ====================================================================
 
    void setLimit (double, double);
    void setLimitLog (double, double);
    void setLimit (double, double, double, double);


    // ======================================================================
    // ============== methods to model two-point correlations  ==============
    // ======================================================================

    // fit the z-space correlation function with a double power-law function
    void fit_xi_dpl (vector<double>, vector<double> &, string &, string file_cov="NULL", string dir_in="NULL");

    // fit the projected correlation function
    void fit_projected_xi (vector<double>, vector<double> &, string &, Cosmology &, double &, bool &, vector<string>, string file_cov="NULL", string dir_in="NULL", double median_redshift=-1, double median_magnitude=-1, double median_lgMass=-1);


    // ======================================================================
    // ============== methods to model clustering anisotropies ==============
    // ======================================================================

    // compute the (linear or non-linear) model xi(rp,pi) using the measured xi(r)
    void compute_xi2D_modelXiMeasured (vector< vector<double> > &, double &, double &, Cosmology &, double &, bool &, int &, double &, string &, double &, int &, bool bias_nl=0, double bA=0., double v_min=-3000., double v_max=3000., int step_v=500);

    // compute the (linear or non-linear) model xi(rp,pi) 
    void compute_xi2D_model (vector< vector<double> > &, double &, double &, double &, Cosmology &, double &, bool &, string &, string &, int &, bool bias_nl=0, double bA=0., bool xiType=0, bool xiNL=0, double v_min=-3000., double v_max=3000., int step_v=500, int norm=-1, double r_min=0.1, double r_max=150., double k_min=0., double k_max=100., double aa=0., bool GSL=1, double prec=1.e-2, string file_par="NULL");

    // compute the multipoles of the model (non-)linear xi(rp,pi) using the measured xi(r)
    void compute_multipoles_modelXiMeasured (vector<double> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &, double &, double &, Cosmology &, double &, bool &, int &, double &, string &, double &, int &, bool bias_nl=0, double bA=0., double v_min=-3000., double v_max=3000., int step_v=500, int step_cos=3000);

    // compute the effective multipoles of the model (non-)linear xi(rp,pi) using the measured xi(r) 
    void compute_effective_multipoles_modelXiMeasured (vector<double> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &, double &, double &, Cosmology &, double &, bool &, int &, double &, string &, double &, int &, bool bias_nl=0, double bA=0., double v_min=-3000., double v_max=3000., int step_v=500);
  
    // measure the chi2, given beta, modelling xi(s)/xi(r) at large scales
    double chi2_beta_KaiserLimit (vector<double>, double &, Cosmology &, double &, int &, double &, string &, int &, double &, string &); 

    // measure the chi2, given beta, modelling the linear xi(rp,pi)
    double chi2_beta (vector<double>, double &, int &, Cosmology &, double &, int &, double &, string &, int &, double &, string &); 

    // measure the chi2, given beta and sigma12, modelling the non-linear xi(rp,pi)
    double chi2_beta_sigma12 (double &, int &, int &, Cosmology &, double &, int &, double &, string &, int &, double &, double &, string &, string &);

    // write the chi2 map
    void write_chi2 (string &, double &, double &, vector<int> &, int &, double &);

    // measure beta, modelling xi(s)/xi(r) at large scales
    void measure_beta_KaiserLimit_XiMeasured (int &, double &, string &, double &, double beta_guess=1.);

    // measure beta, modelling xi(rp,pi) with the linear Kaiser model using the measured xi(r)
    void measure_beta_XiMeasured (int &, int &, double &, string &, double &, double beta_guess=1.);

    // measure beta and sigma12, modelling xi(rp,pi) with the dispersion model and using the measured xi(r)
    void measure_beta_sigma12_DispersionModel_XiMeasured (int &, int &, Cosmology &, double &, int &, double &, string &, double &, string fileMCMC="NULL", double v_min=-3000., double v_max=3000., int step_v=500);

    // measure beta and sigma12, modelling the multipoles of xi(rp,pi) with the dispersion model and using the measured xi(r)
    void measure_beta_sigma12_DispersionModel_XiMeasured_multipoles (int &, int &, Cosmology &, double &, int &, double &, string &, double &, double v_min=-3000., double v_max=3000., int step_v=500);

    // measure beta, modelling xi(s) at large scales
    void measure_fsigma8_KaiserLimit (Cosmology &, double &, string &, double &, string &, int &, double fsigma8_guess=1., bool xiType=0, double k_star=-1., bool xiNL=0, int norm=-1, double r_min=0.1, double r_max=150., double k_min=0., double k_max=100., double aa=0., bool GSL=1, double prec=1.e-2, string file_par="NULL", int rank=0);

    // measure f*sigma8 and bias*sigma8 (and sigma12) modelling xi(rp,pi) with the dispersion model
    void measure_fsigma8_bsigma8 (Cosmology &, string &, double &, string &, int &, int &, bool &, vector<vector<double>>, vector<double>, bool bias_nl=0, bool xiType=0, double k_star=-1., bool xiNL=0, double v_min=-3000., double v_max=3000., int step_v=500, int norm=-1, double r_min=0.1, double r_max=150., double k_min=0., double k_max=100., double aa=0., bool GSL=1, double prec=1.e-2, string file_par="NULL", int rank=0);

    // expected error on beta (from Bianchi et al. 2012)
    double relative_error_beta_catalogue (double &, double &, double &, bool proj=0, bool NL=0); 
    double relative_error_beta_catalogue_box (double &, double &, bool proj=0, bool NL=0); 
  
    // store xi(rp,pi) estimated with the dispersion model with measure xi(r) (for SM and gnuplot)
    void write_xi2D_DispersionModelXiMeasured (string &, double &, double &, Cosmology &, double &, bool &, int &, double &, string &, double &, int &, bool bias_nl=0, double bA=0., double v_min=-3000., double v_max=3000., int step_v=500, int rank=0); 
  
    // store the effective multipoles, estimated with the dispersion model with measure xi(r)
    void write_multipoles_DispersionModelXiMeasured (string &, string &, double &, double &, Cosmology &, double &, bool &, int &, double &, string &, double &, int &, bool effective=1, double v_min=-3000., double v_max=3000., int step_v=500, int step_cos=3000, int rank=0); 

    // store xi(rp,pi) estimated with the dispersion model with theoretical xi(r) (for SM and gnuplot)
    void write_xi2D_DispersionModel (string &, double &, double &, double &, Cosmology &, double &, bool &, string &, string &, int &, bool bias_nl=0, double bA=0., bool xiType=0, bool xiNL=0, double v_min=-3000., double v_max=3000., int step_v=500, int norm=-1, double rmin=3., double rmax=200., double k_min=0., double k_max=100., double aa=0., bool GSL=1, double prec=1.e-2, string file_par="NULL", int rank=0);
  
    // store xi(rp,pi) estimated with the Chuang&Wang model (for SM and gnuplot)
    void write_xi2D_CWModel (string &, double &, double &, double &, double &, double &, double &, double &, Cosmology &, double &, string &, bool BAO=1, bool xiType=0, bool xiNL=0, double r_min=0.1, double r_max=150., double v_min=-3000., double v_max=3000., int step_v=500, double k_min=0., double k_max=100., double x_min=-3000., double x_max=3000., int step_x=500, double aa=0., bool GSL=1, double prec=1.e-2, string file_par="NULL", int rank=0);

    // store a table for SM and gnuplot
    void write_map (vector< vector<double> >, string &, int rank=0);
  

    // minimum mass of the host DM haloes, estimed from the mean 2PCF bias
    double MhaloMin (Cosmology &, double &, double &, double &, string &, string &, string &, string &, double &, double &, double &, double &, bool proj=0, bool NL=0);


    // =========================================================
    // ============== methods to access variables ============== 
    // =========================================================

    double beta_best ()        { return m_beta_best; }
    double error_beta ()       { return m_error_beta; }
    double sigma12_best ()     { return m_sigma12_best; }
    double error_sigma12 ()    { return m_error_sigma12; }
    double bias_best ()        { return m_bias_best; }
    double bias_sigma8_best () { return m_bias_sigma8_best; }
    double f_sigma8_best ()    { return m_f_sigma8_best; }
    double bA_best()           { return m_bA_best; }
    double E_min ()            { return m_E_min; }

  };


  // =====================================================================================


  namespace glob {

    struct STR_xi0K
    {
      Cosmology cosm;
      double bias_sigma8;
      string author;
      double redshift;
      string Model;
      bool xiType; 
      double k_star; 
      bool xiNL;
      double v_min;
      double v_max;
      int step_v;
      int norm;
      double r_min;
      double r_max;
      double k_min;
      double k_max;
      double aa;
      bool GSL;
      double prec;
      string file_par;
    };
  }

}

#include "ModelTwoPointClassFunc.h"

#endif
