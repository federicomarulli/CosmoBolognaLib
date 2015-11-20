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
  class ThreePointCorrelation : public TwoPointCorrelation {

  private :

    // =================================================
    // ============== internal parameters ==============
    // =================================================

    /// number of object-object-object triplets
    vector<double> m_ggg;

    /// number of random-random-random triplets
    vector<double> m_rrr;

    /// number of object-object-random triplets
    vector<double> m_ggr;
    
    /// number of object-random-random triplets
    vector<double> m_grr;

    /// connected three-point correlation function
    vector<double> m_zeta;
    
    /// error on the connected three-point correlation function
    vector<double> m_error_zeta;
    
    /// reduced three-point correlation function
    vector<double> m_zeta_red;
    
    /// error on the reduced three-point correlation function
    vector<double> m_error_zeta_red;

    /// size of the bins
    double m_binsize;
    
    /// lenght of one of the triangle sides
    double m_side_s;
    
    /// lenght of one of the triangle sides
    double m_side_u;

    /// percentage difference between two of the triangle sides
    double m_perc_increase;

    /// number of bins
    int m_nbins;

    /// type of binning: lin &rarr; linear binning; ang &rarr; angular binning
    string m_type_binning;


    // =====================================================================
    // ============== methods to count the number of triplets ==============
    // =====================================================================

    /**
     * @brief method to count the number of triplets
     *
     * @param cat1 object of class \e Catalogue: the input catalogue
     * to analyse
     *
     * @param ChainMesh_rMAX1 object of class \e ChainMesh_Catalogue:
     * the chain mesh used to count the pairs relative to the first
     * side of the triangle (1-2)
     *
     * @param ChainMesh_rMAX2 object of class \e ChainMesh_Catalogue:
     * the chain mesh used to count the pairs relative to the second
     * side of the triangle (1-3)
     *
     * @param tt object of class \e Triplets
     *
     * @param do3D 1 &rarr; compute the 3D three-point correlation
     * function; 0 &rarr; compute the angular three-point correlation
     * function
     *
     * @param tcount 1 &rarr; activate the CPU time counter; 0 &rarr; no time counter
     * @return none
     *
     * @warning the angular three-point correlation function is not
     * implemented yet
     */
    void count_triplets (const shared_ptr<Catalogue>, const ChainMesh_Catalogue &, const ChainMesh_Catalogue &, Triplets &, const bool, const bool tcount=0);

    /**
     * @brief method to count the number of triplets of type:
     * object-object-*
     *
     * @param GGG object of class \e Triplets: number of
     * object-object-object-triplets
     *
     * @param GGR object of class \e Triplets: number of
     * object-object-random-triplets
     *
     * @param ChainMesh_rMAX1 object of class \e ChainMesh_Catalogue:
     * the chain mesh used to count the pairs relative to the first
     * side of the triangle (1-2)
     *
     * @param ChainMesh1_rMAX2 object of class \e ChainMesh_Catalogue:
     * first chain mesh used to count the pairs relative to the second
     * side of the triangle (1-3)
     *
     * @param ChainMesh2_rMAX2 object of class \e ChainMesh_Catalogue:
     * second chain mesh used to count the pairs relative to the second
     * side of the triangle (1-3)
     *
     * @param do3D 1 &rarr; compute the 3D three-point correlation
     * function; 0 &rarr; compute the angular three-point correlation
     * function
     *
     * @param tcount 1 &rarr; activate the CPU time counter; 0 &rarr; no time counter
     * @return none
     *
     * @warning the angular three-point correlation function is not
     * implemented yet
     */
    void count_triplets_gg (Triplets &, Triplets &, const ChainMesh_Catalogue &, const ChainMesh_Catalogue &, const ChainMesh_Catalogue &, const bool, const bool tcount=0);

     /**
      * @brief method to count the number of triplets of type:
     * random-random-*
     *
     * @param RRG object of class \e Triplets: number of
     * random-random-object-triplets
     *
     * @param RRR object of class \e Triplets: number of
     * random-random-random-triplets
     *
     * @param ChainMesh_rMAX1 object of class \e ChainMesh_Catalogue:
     * the chain mesh used to count the pairs relative to the first
     * side of the triangle (1-2)
     *
     * @param ChainMesh1_rMAX2 object of class \e ChainMesh_Catalogue:
     * first chain mesh used to count the pairs relative to the second
     * side of the triangle (1-3)
     *
     * @param ChainMesh2_rMAX2 object of class \e ChainMesh_Catalogue:
     * second chain mesh used to count the pairs relative to the second
     * side of the triangle (1-3)
     *
     * @param do3D 1 &rarr; compute the 3D three-point correlation
     * function; 0 &rarr; compute the angular three-point correlation
     * function
     *
     * @param tcount 1 &rarr; activate the CPU time counter; 0 &rarr; no time counter
     * @return none
     *
     * @warning the angular three-point correlation function is not
     * implemented yet
     */
    void count_triplets_rr (Triplets &, Triplets &, const ChainMesh_Catalogue &, const ChainMesh_Catalogue &, const ChainMesh_Catalogue &, const bool, const bool tcount=0);

     /**
     * @brief method to count the number of triplets of type:
     * object-random-*
     *
     * @param GRG object of class \e Triplets: number of
     * object-random-object-triplets
     *
     * @param GRR object of class \e Triplets: number of
     * object-random-random-triplets
     *
     * @param ChainMesh_rMAX1 object of class \e ChainMesh_Catalogue:
     * the chain mesh used to count the pairs relative to the first
     * side of the triangle (1-2)
     *
     * @param ChainMesh1_rMAX2 object of class \e ChainMesh_Catalogue:
     * first chain mesh used to count the pairs relative to the second
     * side of the triangle (1-3)
     *
     * @param ChainMesh2_rMAX2 object of class \e ChainMesh_Catalogue:
     * second chain mesh used to count the pairs relative to the second
     * side of the triangle (1-3)
     *
     * @param do3D 1 &rarr; compute the 3D three-point correlation
     * function; 0 &rarr; compute the angular three-point correlation
     * function
     *
     * @param tcount 1 &rarr; activate the CPU time counter; 0 &rarr; no time counter
     * @return none
     *
     * @warning the angular three-point correlation function is not
     * implemented yet
     */
    void count_triplets_gr (Triplets &, Triplets &, const ChainMesh_Catalogue &, const ChainMesh_Catalogue &, const ChainMesh_Catalogue &, const bool, const bool tcount=0);

     /**
     * @brief method to count the number of triplets of type:
     * random-object-*
     *
     * @param RGG object of class \e Triplets: number of
     * random-object-object-triplets
     *
     * @param RGR object of class \e Triplets: number of
     * rahdom-object-random-triplets
     *
     * @param ChainMesh_rMAX1 object of class \e ChainMesh_Catalogue:
     * the chain mesh used to count the pairs relative to the first
     * side of the triangle (1-2)
     *
     * @param ChainMesh1_rMAX2 object of class \e ChainMesh_Catalogue:
     * first chain mesh used to count the pairs relative to the second
     * side of the triangle (1-3)
     *
     * @param ChainMesh2_rMAX2 object of class \e ChainMesh_Catalogue:
     * second chain mesh used to count the pairs relative to the second
     * side of the triangle (1-3)
     *
     * @param do3D 1 &rarr; compute the 3D three-point correlation
     * function; 0 &rarr; compute the angular three-point correlation
     * function
     *
     * @param tcount 1 &rarr; activate the CPU time counter; 0 &rarr;
     * no time counter
     *
     * @return none
     *
     * @warning the angular three-point correlation function is not
     * implemented yet
     */
    void count_triplets_rg (Triplets &, Triplets &, const ChainMesh_Catalogue &, const ChainMesh_Catalogue &, const ChainMesh_Catalogue &, const bool, const bool tcount=0);
    
    /**
     * @brief write the binned number of triplets to a file
     * @param TT vector containing the number of triplets
     * @param dir vector containing the input directories
     * @param file name of the input files
     * @return none
     */
    void write_triplets (const vector<double>, const string, const string);

    /**
     * @brief read the files where the binned number of triplets are
     * stored
     * @param TT vector containing the number of triplets
     * @param dir vector containing the input directories
     * @param file name of the input files
     * @return none
     */
    void read_triplets (vector<double>, const vector<string>, const string);


  public :

    // ======================================================
    // ============== constructors/destructors ==============
    // ======================================================

    /**
     * @brief default constructor
     * @return object of class ThreePointCorrelation
     */
    ThreePointCorrelation () {};

    /**
     * @brief constructor 
     * @param real vector of objects of type \e Catalogue containing
     * the input data catalogue
     * @param random vector of objects of type \e Catalogue containing
     * the random data catalogue
     * @return object of class Catalogue
     */
    ThreePointCorrelation (const shared_ptr<Catalogue> real, const shared_ptr<Catalogue> random) 
      : TwoPointCorrelation (real, random, 0) {} ;

    /**
     * @brief default destructor
     * @return none
     */
    ~ThreePointCorrelation () {};
    

    // ====================================================================
    // ============== methods to set the internal parameters ==============
    // ====================================================================

    /**
     * @brief default destructor
     * @return none
     */
    void setParameters3p (const double, const double, const double, const string, const int);

    /**
     * @brief vector allocation and inizialization
     * @return none
     */
    void allocate_vectors_zeta ();



    // ==========================================================================
    // ============== methods to measure three-point correlations  ==============
    // ==========================================================================

    /**
     * @brief method to measure the three-point correlation function
     *
     * @param dir_output_triplets name of the output directory used to
     * store the number of triplets
     *
     * @param dir_output_2pt name of the output directory used to
     * store the two-point correlation functions
     * 
     * @param dir_input_triplets name of the input directories containing the number of triplets
     *
     * @param count_ggg 1 &rarr; count the object-object-object
     * triplets; 0 &rarr; read the object-object-object triplets from
     * a file
     *
     * @param count_rrr 1 &rarr; count the random-random-random
     * triplets; 0 &rarr; read the random-random-random triplets from
     * a file
     *
     * @param count_ggr 1 &rarr; count the object-object-random
     * triplets; 0 &rarr; read the object-object-random triplets from
     * a file
     *
     * @param count_grr 1 &rarr; count the object-random-random
     * triplets; 0 &rarr; read the object-random-random triplets from
     * a file
     *
     * @param tcount 1 &rarr; activate the CPU time counter; 0 &rarr;
     * no time counter
     *
     * @return none
     */
    void measure_Q (const string dir_output_triplets="NULL", const string dir_output_2pt="NULL", const vector<string>& dir_input_triplets=vector<string>(), const int count_ggg=1, const int count_rrr=1, const int count_ggr=1, const int count_grr=1, const bool tcount=0);


    /**
     * @brief overloading of the method to measure the three-point
     * correlation function
     *
     * @param dir_output_triplets name of the output directory used to
     * store the number of triplets
     *
     * @param dir_output_2pt name of the output directory used to
     * store the two-point correlation functions
     *
     * @param tcount 1 &rarr; activate the CPU time counter; 0 &rarr;
     * no time counter
     *
     * @return none
     */
    void measure_Q (const string dir_output_triplets, const string dir_output_2pt, const bool tcount)
    {
      measure_Q(dir_output_triplets, dir_output_2pt, {}, 1, 1, 1, 1, tcount);
    };

     /**
     * @brief overloading of the method to measure the three-point
     * correlation function
     *
     * @param dir_input_triplets name of the input directories
     * containing the number of triplets
     *
     * @param tcount 1 &rarr; activate the CPU time counter; 0 &rarr;
     * no time counter
     *
     * @return none
     */
    void measure_Q (vector<string> dir_input_triplets, const bool tcount)
    {
      measure_Q("NULL", "NULL", dir_input_triplets, 0, 0, 0, 0, tcount);
    };

    
    /// @cond TEMP_TEST
    
    void measure_Q_TEST (const string dir_output_triplets="NULL", const string dir_output_2pt="NULL", const vector<string>& dir_input_triplets=vector<string>(), const bool count=1, const bool tcount=0);
    
    void measure_Q_ang (const string dir_output_triplets="NULL", const string dir_output_2pt="NULL", const vector<string>& dir_input_triplets=vector<string>(), const int count_ggg=1, const int count_rrr=1, const int count_ggr=1, const int count_grr=1, const bool tcount=0);

    void measure_Q_ang (const string dir_output_triplets, const string dir_output_2pt, const bool tcount)
    {
      measure_Q_ang(dir_output_triplets, dir_output_2pt, {}, 1, 1, 1, 1, tcount);
    };

    void measure_Q_ang (const vector<string> dir_input_triplets, const bool tcount)
    {
      measure_Q_ang("NULL", "NULL", dir_input_triplets, 0, 0, 0, 0, tcount);
    };

    /// @endcond

    
     /**
     * @brief write the measured three-point correlation functions
     * @param dir name of the output directory
     * @return none
     */
    void write_Q (const string);

  };
}

#endif
