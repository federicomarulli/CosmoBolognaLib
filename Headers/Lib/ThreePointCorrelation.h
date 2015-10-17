/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli, Michele Moresco         *
 *  and Alfonso Veropalumbo                                         *
 *                                                                  *
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
 *  @file Headers/Lib/ThreePointCorrelation.h
 *
 *  @brief The class ThreePointCorrelation
 *
 *  This file defines the interface of the class ThreePointCorrelation,
 *  used to measure the three-point correlation function
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#ifndef __THREEPOINT__
#define __THREEPOINT__ 

#include "Triplets.h"


// ===================================================================================================


namespace cosmobl {

  /**
   *  @class ThreePointCorrelation ThreePointCorrelation.h
   * "Headers/Lib/ThreePointCorrelation.h"
   *
   *  @brief The class ThreePointCorrelation
   *
   *  This class is used to handle objects of type <EM>
   *  ThreePointCorrelation </EM>. It is used to measure the
   *  three-point correlation function.
   */
  class ThreePointCorrelation : TwoPointCorrelation {

  private :

    // =================================================
    // ============== internal parameters ==============
    // =================================================

    // number of triplets
    vector<double> m_ggg, m_rrr, m_ggr, m_grr;

    // measured correlation functions
    vector<double> m_zeta, m_error_zeta, m_zeta_red, m_error_zeta_red;

    // parameters used to measure the 3-point correlation function
    double m_N_R3p, m_BOX3p, m_binsize, m_side_s, m_side_u, m_perc_increase;
    int m_nbins, m_Nboot3p, m_nCosm3p;
    string m_type_binning;


    // =====================================================================
    // ============== methods to count the number of triplets ==============
    // =====================================================================

    // count the number of triplets
    void count_triplets (shared_ptr<Catalogue>, ChainMesh_Catalogue &, ChainMesh_Catalogue &, Triplets &, bool do_3D, bool tcount=0);

    void count_triplets_gg (Triplets &, Triplets &, ChainMesh_Catalogue &, ChainMesh_Catalogue &, ChainMesh_Catalogue &, bool do_3D, bool tcount=0);
    void count_triplets_rr (Triplets &, Triplets &, ChainMesh_Catalogue &, ChainMesh_Catalogue &, ChainMesh_Catalogue &, bool do_3D, bool tcount=0);
    void count_triplets_gr (Triplets &, Triplets &, ChainMesh_Catalogue &, ChainMesh_Catalogue &, ChainMesh_Catalogue &, bool do_3D, bool tcount=0);
    void count_triplets_rg (Triplets &, Triplets &, ChainMesh_Catalogue &, ChainMesh_Catalogue &, ChainMesh_Catalogue &, bool do_3D, bool tcount=0);
    
    // simple count of object triplets (to test performances)
    void count_triplets_direct (vector<double> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &);

    // write and read the number of triplets
    void write_triplets (vector<double> &, string &, string);
    void read_triplets (vector<double> &, vector<string> &, string);


  public :


    // ==========================================================
    // ============== constructors and destructors ==============
    // ==========================================================

    ThreePointCorrelation () {};
    ~ThreePointCorrelation () {};

    ThreePointCorrelation (shared_ptr<Catalogue> real, shared_ptr<Catalogue> random) 
      : TwoPointCorrelation (real, random, 0) {} ;

  

    // ====================================================================
    // ============== methods to set the internal parameters ==============
    // ====================================================================

    // set the parameters 
    void setParameters3p (double &, double &, double &, string &, int &);

    // vector allocation and inizialization
    void allocate_vectors_zeta ();



    // ==========================================================================
    // ============== methods to measure three-point correlations  ==============
    // ==========================================================================

    // measure the 3-point correlation function
    void measure_Q (string dir_output_triplets="NULL", string dir_output_2pt="NULL", vector<string> dir_input_triplets={}, int count_ggg=1, int count_rrr=1, int count_ggr=1, int count_grr=1, bool tcount=0);
    void measure_Q (string dir_output_triplets, string dir_output_2pt, bool tcount)
    {
      measure_Q(dir_output_triplets, dir_output_2pt, {}, 1, 1, 1, 1, tcount);
    };
    void measure_Q (vector<string> dir_input_triplets, bool tcount)
    {
      measure_Q("NULL", "NULL", dir_input_triplets, 0, 0, 0, 0, tcount);
    };

    void measure_Q_TEST (string dir_output_triplets="NULL", string dir_output_2pt="NULL", vector<string> dir_input_triplets={}, bool count=1, bool tcount=0);
    
    
    void measure_Q_ang (string dir_output_triplets="NULL", string dir_output_2pt="NULL", vector<string> dir_input_triplets={}, int count_ggg=1, int count_rrr=1, int count_ggr=1, int count_grr=1, bool tcount=0);
    void measure_Q_ang (string dir_output_triplets, string dir_output_2pt, bool tcount)
    {
      measure_Q_ang(dir_output_triplets, dir_output_2pt, {}, 1, 1, 1, 1, tcount);
    };
    void measure_Q_ang (vector<string> dir_input_triplets, bool tcount)
    {
      measure_Q_ang("NULL", "NULL", dir_input_triplets, 0, 0, 0, 0, tcount);
    };
    
    
    // write the outputs
    void write_Q (string &);

  };
}

#endif
