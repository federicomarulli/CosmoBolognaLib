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
 *  @file Headers/Lib/TwoPointCorrelation.h
 *
 *  @brief The class TwoPointCorrelation
 *
 *  This file defines the interface of the class TwoPointCorrelation, 
 *  used to measure the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINT__
#define __TWOPOINT__ 


#include "ChainMesh_Catalogue.h"
#include "Pairs.h"


// ===================================================================================================


namespace cosmobl {

  /**
   *  @class TwoPointCorrelation TwoPointCorrelation.h
   * "Headers/Lib/TwoPointCorrelation.h"
   *
   *  @brief The class TwoPointCorrelation
   *
   *  This class is used to handle objects of type <EM>
   *  TwoPointCorrelation </EM>. It is used to measure the two-point
   *  correlation function.
   */
  class TwoPointCorrelation {

  protected :

    // =================================================
    // ============== internal parameters ==============
    // =================================================

    shared_ptr<Catalogue> m_data;
    shared_ptr<Catalogue> m_random;

    // number of real and random objects
    double m_nGal, m_nRan;

    // parameters of interal use
    bool m_read_only, m_do_all;

    // number of pairs
    vector<double> m_gg_log, m_rr_log, m_gr_log, m_gg_lin, m_rr_lin, m_gr_lin;
    vector< vector<double> > m_gg_2d, m_rr_2d, m_gr_2d, m_gg_slog, m_rr_slog, m_gr_slog,
      m_gg_coslog, m_rr_coslog, m_gr_coslog, m_gg_coslin, m_rr_coslin, m_gr_coslin;
     
    // measured correlation functions
    vector<double> m_rad_log, m_xi_log, m_error_xi_log, m_rad_lin, m_xi_lin, m_error_xi_lin, m_cos_lin; // 1D
    vector< vector<double> > m_xi_2d_lin, m_error_xi_2d_lin, m_xi_2d_loglin, m_error_xi_2d_loglin,
      m_xi_coslog, m_error_xi_coslog, m_xi_coslin, m_error_xi_coslin; // 2D

    // projected correlation function
    vector<double> m_rp_proj, m_xi_proj, m_error_xi_proj;

    // real space correlation function
    vector<double> m_rad_real_log, m_xi_real_log, m_rad_real_lin, m_xi_real_lin, m_error_xi_real_lin, m_xi_real_lin_extr;

    // real space correlation function and integral forms interpolated in the grid used for the model xi(rp,pi)
    vector<double> m_xi_real_lin_interp, m_xi_interp, m_xi__interp;

    // multipoles of xi(r,mu)
    vector<double> m_xi0_log, m_xi2_log, m_xi4_log, m_error_xi0_log, m_error_xi2_log, m_error_xi4_log, 
      m_xi0_lin, m_xi2_lin, m_xi4_lin, m_error_xi0_lin, m_error_xi2_lin, m_error_xi4_lin, m_quad, m_error_quad;

    // bias
    vector<double> m_rr_bias_lin_xi, m_bias_lin_xi, m_error_bias_lin_xi, m_rr_bias_nl_xi, m_bias_nl_xi, m_error_bias_nl_xi, 
      m_rr_bias_lin_wp, m_bias_lin_wp, m_error_bias_lin_wp, m_rr_bias_nl_wp, m_bias_nl_wp, m_error_bias_nl_wp;
    double m_bias, m_bias_min, m_bias_max;

    // parameters used to measure the correlation function 
    double m_rMIN, m_rMAX, m_rMAX_eff, m_N_R, m_BOX, 
      m_shift_log, m_shift_lin, m_shift_cos, m_logbinsz, m_linbinsz, m_logbinSize, m_binSize, m_cosSize, m_cosbinsz;
    int m_nlogbins, m_nlinbins, m_ncosbins, m_Nboot, m_nCosm;
    
    // parameters of the linear extrapolation of xi(r) at small scales
    double m_r0_linextr, m_gamma_linextr;


    /* ======== Alfonso Veropalumbo ======== */

    // number of pairs for the angular correlation function
    vector<double> m_gg_theta_log, m_rr_theta_log, m_gr_theta_log, m_gg_theta_lin, m_rr_theta_lin, m_gr_theta_lin;

    // measured angular correlation function
    vector<double> m_theta_log, m_wtheta_log, m_error_wtheta_log, m_theta_lin, m_wtheta_lin, m_error_wtheta_lin;

    // parameters used to measure the angular correlation function 
    double m_thetaMIN, m_thetaMAX, m_thetaMAX_eff, m_shift_theta_log, m_shift_theta_lin, m_logthetabinsz, m_linthetabinsz, m_logthetabinSize, m_linthetabinSize;
    int m_nlogthetabins, m_nlinthetabins;

    vector<shared_ptr<TwoPointCorrelation> > m_twop_mock;


    // ==================================================================
    // ============== methods to count the number of pairs ==============
    // ==================================================================

    // count the number of pairs
    void count_pairs (shared_ptr<Catalogue>, ChainMesh_Catalogue &, Pairs &, bool, bool tcount=0);

    // simple count of object pairs (to test performances)
    void count_pairs_direct (Catalogue &, Catalogue &); 

    // create a chain mesh
    void create_chain_mesh (Catalogue &, vector<double>, double &, double &, double &, vector<double> &, vector< vector< vector<double> > > &, bool cross=0, bool order=0);


    // ==========================================
    // ============== I/O methods  ==============
    // ==========================================

    // write/read the number of pairs
    void write_pairs (vector<double> &, vector<double> &, vector< vector<double> > &, vector< vector<double> > &, vector< vector<double> > &, vector< vector<double> > &, string &, string &);
    void read_pairs (vector<double> &, vector<double> &, vector< vector<double> > &, vector< vector<double> > &, vector< vector<double> > &, vector< vector<double> > &, vector<string> &, string &);


    /* ======== Alfonso Veropalumbo ======== */
    // write/read the number of pairs for the angular correlation function
    void write_pairs (vector<double> &, vector<double> &, string &, string &);
    void read_pairs (vector<double> &, vector<double> &,  vector<string> &, string &);

  
  public :

    // ==========================================================
    // ============== constructors and destructors ==============
    // ==========================================================

    TwoPointCorrelation () { m_read_only = 1; m_do_all = 1; }

    ~TwoPointCorrelation () {}
    
    TwoPointCorrelation (shared_ptr<Catalogue>, shared_ptr<Catalogue>, bool do_all=1); 
  

    // ====================================================================
    // ============== methods to set the internal parameters ==============
    // ====================================================================

    // set the parameters used to measure the correlation function
    void setParameters (double &, double &, double &, double &, double &, bool ANG=0);

    // vector allocation and inizialization
    void allocate_vectors_xi (bool &);
    void allocate_vectors_ACF (bool &);

    // erase the vectors with multipoles
    void erase_multipoles ();


    // ===============================================================================================
    // ============== methods to measure two-point correlations and related quantities  ==============
    // ===============================================================================================

    // measure the correlation function (xi(r) and xi(rp,pi))
    void measure_xi (string, vector<string>, int, int, int, bool, bool tcount=0);
    void measure_xi (string dir_output_pairs, bool tcount=0)
    {
      measure_xi(dir_output_pairs, {}, 1, 1, 1, 1, tcount);
    };
    
    void measure_xi (vector<string> dir_input_pairs)
    {
      measure_xi("NULL", dir_input_pairs, 0, 0, 0, 1, 0);
    };
    
    void measure_xi (vector<string> dir_input_pairs, int count_gg, int count_rr, int count_gr, bool doGR, bool tcount=0)
    {
      measure_xi("NULL", dir_input_pairs, count_gg, count_rr, count_gr, doGR, tcount);
    };
    
    void measure_xi (string dir_output, vector<string> dir_input_pairs)
    {
      string file_in = dir_output+"nObjects";
      ifstream fin (file_in.c_str()); checkIO(file_in, 1);
      fin >>m_nGal>>m_nRan;
      fin.clear(); fin.close();
      measure_xi("NULL", dir_input_pairs, 0, 0, 0, 1, 0);
    };
    
    // write the outputs
    void write_xi (string &, int rank=0);
    void write_xi (vector< vector<double> > &, string &, int rank=0);

    
    /* ======== Alfonso Veropalumbo ======== */

    // ================================================================================
    // ======================== methods for Internal Error   ==========================
    // ================================================================================

    // count pairs when Catalogue is subsampld in regions
    void count_pairs_regions (shared_ptr<Catalogue> , ChainMesh_Catalogue &, vector<shared_ptr<Pairs> >, bool, bool tcount=0);

    // measure the error with Catalogue subsampling
    void measure_xi (string, string, string, bool, int nSamples=0., vector<string> dir_input_pairs={}, int count_gg=1, int count_rr=1, int count_gr=1, bool doGR=1, string suffix="NULL", bool tcount=0);
    
    void measure_xi (string dir_output_pairs, string dir_output_pairs_subSamples, string dir_output_covariance, string suffix="NULL", bool tcount=0) // to measure the 2pt correlation function and the covariance matrix with jackknife
    {
      measure_xi(dir_output_pairs, dir_output_pairs_subSamples, dir_output_covariance, 1, 0., {}, 1, 1, 1, 1, suffix, tcount);
    }
    
    void measure_xi (string dir_output_pairs, string dir_output_covariance, bool tcount=0) // to measure the 2pt correlation function and the covariance matrix with jackknife without storing object pairs
    {
      measure_xi(dir_output_pairs, "NULL", dir_output_covariance, 1, 0., {}, 1, 1, 1, 1, "NULL", tcount);
    }

    void measure_xi (vector<string> dir_input_pairs, string dir_output_covariance, string suffix="NULL") // to read the 2pt correlation function and the covariance matrix with jackknife
    {
      measure_xi("NULL", "NULL", dir_output_covariance, 1, 0., dir_input_pairs, 0, 0, 0, 1, suffix, 0);
    }

    void measure_xi (string dir_output_pairs, string dir_output_pairs_subSamples, string dir_output_covariance, int nSamples, string suffix="NULL", bool tcount=0) // to measure the 2pt correlation function and the covariance matrix with bootstrap
    {
      measure_xi(dir_output_pairs,  dir_output_pairs_subSamples, dir_output_covariance, 0, nSamples, {}, 1, 1, 1, 1, suffix, tcount);
    }

    void measure_xi (vector<string> dir_input_pairs, string dir_output_covariance, int nSamples, string suffix="NULL") // to read the 2pt correlation function and the covariance matrix with bootstrap
    {
      measure_xi("NULL", "NULL", dir_output_covariance, 0, nSamples, dir_input_pairs, 0, 0, 0, 1, suffix, 0);
    }

    
    // measure the covariance matrix given a collection of Two Point Correlation Function
    void get_covariance (string &, bool &, string suffix="NULL");
    
    // Pairs
    void set_gg_pairs (vector<double>, vector<double>, vector< vector<double> >, vector< vector<double> >, vector< vector<double> >, vector< vector<double> >);
    void set_rr_pairs (vector<double>, vector<double>, vector< vector<double> >, vector< vector<double> >, vector< vector<double> >, vector< vector<double> >);
    void set_gr_pairs (vector<double>, vector<double>, vector< vector<double> >, vector< vector<double> >, vector< vector<double> >, vector< vector<double> >);

    void write_pairs_subSamples(vector<shared_ptr<Pairs> > , int, string &, string &);
     
    void read_pairs_subSamples (vector<shared_ptr<Pairs> > , int, vector<string> &, string &);
    void write_xi_Mocks (string);
    void measure_projected_xi_Mocks (double, string);

    // measure the angular correlation function
    void measure_wtheta (string dir_output_pairs="NULL", vector<string> dir_input_pairs={}, int count_gg=1, int count_rr=1, int count_gr=1, bool doGR=1, bool tcount=0);
    void measure_wtheta (string dir_output_pairs, bool tcount)
    {
      measure_wtheta(dir_output_pairs, {}, 1, 1, 1, 1, tcount);
    };
    void measure_wtheta (vector<string> dir_input_pairs, bool tcount)
    {
      measure_wtheta("NULL", dir_input_pairs, 0, 0, 0, 1, tcount);
    };
    
    
    // write the outputs
    void write_wtheta (string &);

    // expected error on the correlation function
    double Error (double &, double &);
    double Error (double &, double &, double &);

    // estimeate the full covariance matrix for the 1D correlation function
    void measure_covariance_1D (vector<string> &, vector<double> &, vector< vector<double> > &, string file_out="NULL");

    // measure the projected correlation function
    void measure_projected_xi (double &, string &);

    // measure the multipoles of xi(rp,pi)
    void measure_multipoles_xirmu (double rApprox=0.);
    void measure_multipoles_xirppi (int step_cos=3000);
    void measure_effective_multipoles_xirppi (); // from Chuang&Wang 2012/2013 (check the propagated errors!)

    // write/read the multipoles of the correlation function
    void write_multipoles (string &, string &); 
    void read_multipoles (string &, string &);

    // measure the normalized quadrupole
    void measure_normalized_quadrupole (double rApprox=3.);

    // rms mass fluctuation
    double sigmaR_obj (double, int); 

    // measure the bias 
    void measure_bias (Cosmology &, double &, string &, string &, string file_bias="NULL", bool proj=0);
    double mean_bias (double &, double &, bool proj=0, bool NL=0); 
    double error_mean_bias (double &, double &, bool proj=0, bool NL=0); 
    void fit_bias (double &, double &, bool proj=0, bool NL=0); 
    void read_bias (string &, bool proj=0);
 

    // ================================================================================
    // ============== methods for the real-space two-point correlations  ==============
    // ================================================================================

    // compute the real space correlation function with the deprojection method
    void compute_real_space_xi_deprojected (double &, string &);

    // general method to derive the real space correlation function
    void derive_real_xi (int &, double &, string &, double &);

  
    // ==============================================================
    // ============== methods to get private variables ============== 
    // ==============================================================

    double N_R () {return m_N_R;};
    int nlinbins () {return m_nlinbins;};
    int nlogbins () {return m_nlogbins;};
    int ncosbins () {return m_ncosbins;};
    double linbinsz () {return m_linbinsz;};
    double logbinsz () {return m_logbinsz;};
    double cosbinsz () {return m_cosbinsz;};
    double shift_lin () {return m_shift_lin;};
    double shift_log () {return m_shift_log;};
    double shift_cos () {return m_shift_cos;};
 
    double rad_log (int i) {cosmobl::checkDim(m_rad_log, i, "rad_log"); return m_rad_log[i];};
    double xi_log (int i) {cosmobl::checkDim(m_xi_log, i, "xi_log"); return m_xi_log[i];};
    double error_xi_log (int i) {cosmobl::checkDim(m_error_xi_log, i, "error_xi_log"); return m_error_xi_log[i];};
 
    double rad_lin (int i) {cosmobl::checkDim(m_rad_lin, i, "rad_lin"); return m_rad_lin[i];};
    double xi_lin (int i) {cosmobl::checkDim(m_xi_lin, i, "xi_lin"); return m_xi_lin[i];};
    double error_xi_lin (int i) {cosmobl::checkDim(m_error_xi_lin, i, "error_xi_lin"); return m_error_xi_lin[i];};
  
    double xi_2d_lin (int i, int j) {cosmobl::checkDim(m_xi_2d_lin, i, j,"xi_2d_lin"); return m_xi_2d_lin[i][j];};
    double error_xi_2d_lin (int i, int j) {cosmobl::checkDim(m_error_xi_2d_lin, i, j,"xi_2d_lin"); return m_error_xi_2d_lin[i][j];};
  
    double rad_real_lin (int i) {cosmobl::checkDim(m_rad_real_lin, i, "rad_real_lin"); return m_rad_real_lin[i];};
    double xi_real_lin (int i) {cosmobl::checkDim(m_xi_real_lin, i, "xi_real_lin"); return m_xi_real_lin[i];};
    double xi_real_lin_extr (int i) {cosmobl::checkDim(m_xi_real_lin_extr, i, "xi_real_lin_extr"); return m_xi_real_lin_extr[i];};
    double xi_real_lin_interp (int i) {cosmobl::checkDim(m_xi_real_lin_interp, i, "xi_real_lin_interp"); return m_xi_real_lin_interp[i];};
    double error_xi_real_lin (int i) {cosmobl::checkDim(m_error_xi_real_lin, i, "error_xi_real_lin"); return m_error_xi_real_lin[i];};

    double rad_real_log (int i) {cosmobl::checkDim(m_rad_real_log, i, "rad_real_log"); return m_rad_real_log[i];};
    double xi_real_log (int i) {cosmobl::checkDim(m_xi_real_log, i, "xi_real_log"); return m_xi_real_log[i];};

    double rp_proj (int i) {cosmobl::checkDim(m_rp_proj, i, "rp_proj"); return m_rp_proj[i];};
    double xi_proj (int i) {cosmobl::checkDim(m_xi_proj, i, "xi_proj"); return m_xi_proj[i];};
    double error_xi_proj (int i) {cosmobl::checkDim(m_error_xi_proj, i, "error_xi_proj"); return m_error_xi_proj[i];};
  
    double xi0_log (int i) {cosmobl::checkDim(m_xi0_log, i, "xi0_log"); return m_xi0_log[i];};
    double xi2_log (int i) {cosmobl::checkDim(m_xi2_log, i, "xi2_log"); return m_xi2_log[i];};
    double xi4_log (int i) {cosmobl::checkDim(m_xi4_log, i, "xi4_log"); return m_xi4_log[i];};
    double xi0_lin (int i) {cosmobl::checkDim(m_xi0_lin, i, "xi0_lin"); return m_xi0_lin[i];};
    double xi2_lin (int i) {cosmobl::checkDim(m_xi2_lin, i, "xi2_lin"); return m_xi2_lin[i];};
    double xi4_lin (int i) {cosmobl::checkDim(m_xi4_lin, i, "xi4_lin"); return m_xi4_lin[i];};
    double error_xi0_log (int i) {cosmobl::checkDim(m_error_xi0_log, i, "error_xi0_log"); return m_error_xi0_log[i];};
    double error_xi2_log (int i) {cosmobl::checkDim(m_error_xi2_log, i, "error_xi2_log"); return m_error_xi2_log[i];};
    double error_xi4_log (int i) {cosmobl::checkDim(m_error_xi4_log, i, "error_xi4_log"); return m_error_xi4_log[i];};
    double error_xi0_lin (int i) {cosmobl::checkDim(m_error_xi0_lin, i, "error_xi0_lin"); return m_error_xi0_lin[i];};
    double error_xi2_lin (int i) {cosmobl::checkDim(m_error_xi2_lin, i, "error_xi2_lin"); return m_error_xi2_lin[i];};
    double error_xi4_lin (int i) {cosmobl::checkDim(m_error_xi4_lin, i, "error_xi4_lin"); return m_error_xi4_lin[i];};
    double quad (int i) {cosmobl::checkDim(m_quad, i, "quad"); return m_quad[i];};
    double error_quad (int i) {cosmobl::checkDim(m_error_quad, i, "error_quad"); return m_error_quad[i];};
 
    double rMIN () {return m_rMIN;};
    double rMAX () {return m_rMAX;};
    double rMAX_eff () {return m_rMAX_eff;};
 
    double rr_bias_lin_xi (int i) {cosmobl::checkDim(m_rr_bias_lin_xi, i, "rr_bias_lin_xi"); return m_rr_bias_lin_xi[i];};
    double bias_lin_xi (int i) {cosmobl::checkDim(m_bias_lin_xi, i, "bias_lin_xi"); return m_bias_lin_xi[i];};
    double error_bias_lin_xi (int i) {cosmobl::checkDim(m_error_bias_lin_xi, i, "error_bias_lin_xi"); return m_error_bias_lin_xi[i];};
    double rr_bias_nl_xi (int i) {cosmobl::checkDim(m_rr_bias_nl_xi, i, "rr_bias_nl_xi"); return m_rr_bias_nl_xi[i];};
    double bias_nl_xi (int i) {cosmobl::checkDim(m_bias_nl_xi, i, "bias_nl_xi"); return m_bias_nl_xi[i];};
    double error_bias_nl_xi (int i) {cosmobl::checkDim(m_error_bias_nl_xi, i, "error_bias_nl_xi"); return m_error_bias_nl_xi[i];};
    double rr_bias_lin_wp (int i) {cosmobl::checkDim(m_rr_bias_lin_wp, i, "rr_bias_lin_wp"); return m_rr_bias_lin_wp[i];};
    double bias_lin_wp (int i) {cosmobl::checkDim(m_bias_lin_wp, i, "bias_lin_wp"); return m_bias_lin_wp[i];};
    double error_bias_lin_wp (int i) {cosmobl::checkDim(m_error_bias_lin_wp, i, "error_bias_lin_wp"); return m_error_bias_lin_wp[i];};
    double rr_bias_nl_wp (int i) {cosmobl::checkDim(m_rr_bias_nl_wp, i, "rr_bias_nl_wp"); return m_rr_bias_nl_wp[i];};
    double bias_nl_wp (int i) {cosmobl::checkDim(m_bias_nl_wp, i, "bias_nl_wp"); return m_bias_nl_wp[i];};
    double error_bias_nl_wp (int i) {cosmobl::checkDim(m_error_bias_nl_wp, i, "error_bias_nl_wp"); return m_error_bias_nl_wp[i];};

    double bias() {return m_bias;};
    double bias_min() {return m_bias_min;};
    double bias_max() {return m_bias_max;};
  
    double r0_linextr () {return m_r0_linextr;};
    double gamma_linextr () {return m_gamma_linextr;};

    int sizeof_xi_log () {return m_xi_log.size();};
    int sizeof_xi_lin () {return m_xi_lin.size();};
    int sizeof_xi_real_lin () {return m_xi_real_lin.size();};
    int sizeof_xi_real_lin_extr () {return m_xi_real_lin_extr.size();};
    int sizeof_xi_real_lin_interp () {return m_xi_real_lin_interp.size();};
    int sizeof_xi_proj () {return m_xi_proj.size();};
    int sizeof_xi0_log () {return m_xi0_log.size();};
    int sizeof_xi2_log () {return m_xi2_log.size();};
    int sizeof_xi4_log () {return m_xi4_log.size();};
    int sizeof_xi0_lin () {return m_xi0_lin.size();};
    int sizeof_xi2_lin () {return m_xi2_lin.size();};
    int sizeof_xi4_lin () {return m_xi4_lin.size();};
    int sizeof_quad () {return m_quad.size();};
    int sizeof_bias_lin_xi() {return m_bias_lin_xi.size();};
    int sizeof_bias_nl_xi() {return m_bias_nl_xi.size();};
    int sizeof_bias_lin_wp() {return m_bias_lin_wp.size();};
    int sizeof_bias_nl_wp() {return m_bias_nl_wp.size();};

    int bin_proj () {return m_rp_proj.size();};

    int nObjects () {return m_data->nObjects();};
  
    double Xmin () { vector<double> xx; m_data->MinMax_var(Var::_XX_, xx); return xx[0]; }; 
    double Xmax () { vector<double> xx; m_data->MinMax_var(Var::_XX_, xx); return xx[1]; };
    double Ymin () { vector<double> yy; m_data->MinMax_var(Var::_YY_, yy); return yy[0]; }; 
    double Ymax () { vector<double> yy; m_data->MinMax_var(Var::_YY_, yy); return yy[1]; };
    double Zmin () { vector<double> zz; m_data->MinMax_var(Var::_ZZ_, zz); return zz[0]; }; 
    double Zmax () { vector<double> zz; m_data->MinMax_var(Var::_ZZ_, zz); return zz[1]; };



    // ======================================================
    // ============== methods to set variables ============== 
    // ======================================================

    
    void set_nGal (int nGal) { m_nGal = nGal; }
    void set_nRan (int nRan) { m_nRan = nRan; }
    
    void set_coordinates (vector<double>, vector<double>, vector<double>);
    void set_weight (vector<double>);
    void set_rad_log (vector<double>);
    void set_xi_log (vector<double>);
    void set_error_xi_log (vector<double>);
    void set_rad_lin (vector<double>);
    void set_xi_lin (vector<double>);
    void set_error_xi_lin (vector<double>);
    void set_error_xi_2d_lin (vector< vector<double> >);
    void set_error_xi_2d_loglin (vector< vector<double> >);
    void set_error_xi_coslog (vector< vector<double> >);
    void set_error_xi_coslin (vector< vector<double> >);
    void set_error_proj (vector<double>);
    void set_correlations (vector<double>, vector<double>, vector< vector<double> >, vector< vector<double> >, vector< vector<double> >, vector< vector<double> >);
 
    void set_xi_real_log (vector<double>);
    void set_xi_real_lin (vector<double>);
    void set_xi_real_lin_interp (vector<double>);
    void set_xi_real_lin_extr (vector<double>);
    void set_xi_proj (vector<double>, vector<double>, vector<double>);

  };
}

#endif
