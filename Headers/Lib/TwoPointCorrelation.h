/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Headers/Lib/TwoPointCorrelation.h
 *
 *  @brief The class TwoPointCorrelation
 *
 *  This file defines the interface of the class TwoPointCorrelation,
 *  used to measure the two-point correlation function
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __TWOPOINT__
#define __TWOPOINT__


#include "ChainMesh_Catalogue.h"
#include "Pair1D_extra.h"
#include "Pair2D_extra.h"


// ===================================================================================================


namespace cosmobl {
  
  /**
   *  @brief The namespace of the <B> two-point correlation function
   *  </B>
   *  
   *  The \e twopt namespace contains all the functions and classes to
   *  measure the two-point correlation function
   */
  namespace twopt {

    /**
     * @enum ErrorType
     * @brief the two-point correlation function error type
     */
    enum ErrorType { 

      /// Poissonian error
      _Poisson_,

      /// Jackknife resampling
      _Jackknife_,

      /// Bootstrap resampling
      _Bootstrap_
      
    };

    /**
     * @enum TwoPType
     * @brief the two-point correlation function type
     */
    enum TwoPType { 

      /// the angle-averaged two-point correlation function, i.e. the monopole, &xi;(r)
      _1D_monopole_,

      /// the projected two-point correlation function, w(r<SUB>p</SUB>)
      _1D_projected_,
    
      /// the deprojected two-point correlation function, &xi;(r)
      _1D_deprojected_,
    
      /// the multipoles of the two-point correlation function, &xi;<SUB>i</SUB>(r)
      _1D_multipoles_,

      /// the wedges of the two-point correlation function, &xi;<SUB>i</SUB>(r)
      _1D_wedges_,
      
      /// filtered two-point correlation function
      _1D_filtered_,

      /// angular two-point correlation function
      _1D_angular_,
    
      /// 2D two-point correlation function in Cartesian coordinates, &xi;(r<SUB>p</SUB>,&pi;)
      _2D_Cartesian_,

      /// 2D two-point correlation function in polar coordinates, &xi;(r,&mu;)
      _2D_polar_
    
    };


    /**
     *  @enum Estimator
     *  @brief the estimator used to measure the two-point correlation
     *  function
     */
    enum Estimator { 

      /// natural estimator
      _natural_,

      /// Landy&Szalay estimator
      _LandySzalay_

    };
    
    
    /**
     *  @class TwoPointCorrelation TwoPointCorrelation.h
     *  "Headers/Lib/TwoPointCorrelation.h"
     *
     *  @brief The class TwoPointCorrelation
     *
     *  This is the base class used to measure the two-point
     *  correlation function
     *
     */
    class TwoPointCorrelation {

    protected :

      /**
       *  @name Two-point correlation function data
       */
      ///@{
      
      /// two-point correlation function type
      TwoPType m_twoPType;
      
      /// the dataset of the two-point correlation function
      shared_ptr<data::Data> m_dataset;
      
      ///@}

      
      /**
       *  @name Input and random catalogues
       */
      ///@{
    
      /// input data catalogue
      shared_ptr<catalogue::Catalogue> m_data;

      /// output data catalogue
      shared_ptr<catalogue::Catalogue> m_random;
    
      ///@}


      /**
       *  @name Object pairs
       */
      ///@{
    
      /// number of data-data pairs
      shared_ptr<pairs::Pair> m_dd;

      /// number of random-random pairs
      shared_ptr<pairs::Pair> m_rr;

      /// number of data-random pairs
      shared_ptr<pairs::Pair> m_dr;

      ///@}

      
      /**
       *  @name Other parameters
       */
      ///@{
	 
      /// true &rarr; compute extra information related to the pairs, such as the mean pair separation and redshift
      bool m_compute_extra_info;
      
      ///@}
      
      /**
       *  @name Internal input/output member functions (customized in all the derived classes)
       */
      ///@{
    
      /**
       *  @brief write the number of pairs
       *  @param PP pointer to an object of class Pair
       *  @param dir output directory
       *  @param file output file
       *  @return none
       */
      virtual void write_pairs (const shared_ptr<pairs::Pair> PP, const string dir, const string file) const = 0;
 
      /**
       *  @brief read the number of pairs
       *  @param [out] PP pointer to an object of class Pair
       *  @param [in] dir input directory
       *  @param [in] file input file
       *  @return none
       */
      virtual void read_pairs (shared_ptr<pairs::Pair> PP, const vector<string> dir, const string file) const = 0;

      /**
       *  @brief write the number of pairs
       *  @param PP pointer to a vector of objects of class Pair
       *  @param dir output directory
       *  @param file output file
       *  @return none
       */
      virtual void write_pairs (const vector<shared_ptr<pairs::Pair> > PP, const string dir, const string file) const = 0;
 
      /**
       *  @brief read the number of pairs
       *  @param [out] PP pointer to a vector of objects of class Pair
       *  @param [in] dir input directory
       *  @param [in] file input file
       *  @return none
       */
      virtual void read_pairs (vector<shared_ptr<pairs::Pair> > PP, const vector<string> dir, const string file) const = 0;

      ///@}


      /**
       *  @name Member functions to count the number of pairs 
       */
      ///@{

      /**
       *  @brief count the number of pairs
       *
       *  @param cat1 pointer to an object of class Catalogue,
       *  containing the first catalogue
       *
       *  @param ChM object of class ChainMesh_Catalogue, used to
       *  construct the chain-mesh
       *
       *  @param pp pointer to an object of class Pair
       *
       *  @param cross true &rarr; count the number of cross pairs
       *  (e.g. data-random pairs); false &rarr; count the number of
       *  pairs of objects of the same type (e.g. data-data and
       *  random-random pairs)
       *
       *  @param tcount true &rarr; activate the time counter; false
       *  &rarr; no time counter
       * 
       *  @return none
       */
      void count_pairs (const shared_ptr<catalogue::Catalogue> cat1, const chainmesh::ChainMesh_Catalogue &ChM, shared_ptr<pairs::Pair> pp, const bool cross=true, const bool tcount=false);

      /**
       *  @brief count the number of pairs, used for Jackknife/Bootstrap
       *  methods
       *
       *  @param cat1 pointer to an object of class Catalogue,
       *  containing the first catalogue
       *
       *  @param ChM object of class ChainMesh_Catalogue, used to
       *  construct the chain-mesh
       *
       *  @param pp pointer to an object of class Pair
       *
       *  @param pp_regions vector containing pointers to object of class Pair
       *
       *  @param cross true &rarr; count the number of cross pairs
       *  (e.g. data-random pairs); false &rarr; count the number of
       *  pairs of objects of the same type (e.g. data-data and
       *  random-random pairs)
       *
       *  @param tcount true &rarr; activate the time counter; false
       *  &rarr; no time counter
       * 
       *  @return none
       */
      void count_pairs_region (const shared_ptr<catalogue::Catalogue> cat1, const chainmesh::ChainMesh_Catalogue &ChM, shared_ptr<pairs::Pair> pp, vector< shared_ptr<pairs::Pair> > pp_regions, const bool cross=true, const bool tcount=false);
      
      /**
       *  @brief count the data-data, random-random and data-random
       *  pairs, used to construct the estimator of the two-point
       *  correlation function
       *
       *  @param type the type of two-point correlation function; it
       *  can be: \_1D_monopole\_, \_1D_projected\_,
       *  \_1D_deprojected\_, \_1D_multipoles\_, \_1D_angular\_,
       *  \_2D_Cartesian\_, \_2D_polar\_
       *
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to
       *  store the number of pairs (if the pairs are read from files)
       *
       *  @param count_dd true &rarr; count the number of data-data
       *  pairs; false &rarr; read the number of data-data pairs from
       *  file
       *
       *  @param count_rr true &rarr; count the number of random-random
       *  pairs; false &rarr; read the number of random-random pairs from
       *  file
       *
       *  @param count_dr true &rarr; count the number of data-random
       *  pairs; false &rarr; read the number of data-random pairs from
       *  file
       *
       *  @param tcount true &rarr; activate the time counter; false
       *  &rarr; no time counter
       *
       *  @param estimator the estimator used to measure the two-point
       *  correlation function
       *
       *  @return none
       */
      void count_allPairs (const TwoPType type, const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={}, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_);

      /**
       *  @brief count the data-data, random-random and data-random
       *  pairs, used to construct the estimator of the two-point
       *  correlation function
       *
       *  @param dd_regions data-data pairs in the sub-regions
       *
       *  @param rr_regions random-random pairs in the sub-regions
       *
       *  @param dr_regions data-random pairs in the sub-regions
       *
       *  @param type the type of two-point correlation function; it
       *  can be: \_1D_monopole\_, \_1D_projected\_,
       *  \_1D_deprojected\_, \_1D_multipoles\_, \_1D_angular\_,
       *  \_2D_Cartesian\_, \_2D_polar\_
       *
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to
       *  store the number of pairs (if the pairs are read from files)
       *
       *  @param count_dd true &rarr; count the number of data-data
       *  pairs; false &rarr; read the number of data-data pairs from
       *  file
       *
       *  @param count_rr true &rarr; count the number of random-random
       *  pairs; false &rarr; read the number of random-random pairs from
       *  file
       *
       *  @param count_dr true &rarr; count the number of data-random
       *  pairs; false &rarr; read the number of data-random pairs from
       *  file
       *
       *  @param tcount true &rarr; activate the time counter; false
       *  &rarr; no time counter
       *
       *  @param estimator the estimator used to measure the two-point
       *  correlation function
       *
       *  @return none
       */
       void count_allPairs_region (vector<shared_ptr<pairs::Pair> > &dd_regions, vector<shared_ptr<pairs::Pair> > &rr_regions, vector<shared_ptr<pairs::Pair> > &dr_regions, const TwoPType type, const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={}, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_);

      ///@}

       /**
       *  @name Member functions to compute the two-point correlation
       *  function
       */
      ///@{

       /**
       *  @brief measure the two-point correlation function using the
       *  Natural Estimator
       *  
       *  @param dd pointer to an object of type Pair containing the
       *  data-data pairs
       *
       *  @param rr pointer to an object of type Pair containing the
       *  random-random pairs
       *
       *  @param nData number of objects in the data catalogue
       *
       *  @param nRandom number of objects in the random catalogue
       *
       *  @return pointer to an object of type Data
       */
      virtual shared_ptr<data::Data> NaturalEstimator (const shared_ptr<pairs::Pair> dd, const shared_ptr<pairs::Pair> rr, const int nData, const int nRandom) = 0;

      /**
       *  @brief measure the two-point correlation function using the
       *  Landy-Szalay estimator
       *  
       *  @param dd pointer to an object of type Pair containing the
       *  data-data pairs
       *
       *  @param rr pointer to an object of type Pair containing the
       *  random-random pairs
       *
       *  @param dr pointer to an object of type Pair containing the
       *  data-random pairs
       *
       *  @param nData number of objects in the data catalogue
       *
       *  @param nRandom number of objects in the random catalogue
       *
       *  @return pointer to an object of type Data
       */
      virtual shared_ptr<data::Data> LandySzalayEstimator (const shared_ptr<pairs::Pair> dd, const shared_ptr<pairs::Pair> rr, const shared_ptr<pairs::Pair> dr, const int nData, const int nRandom) = 0;

      /**
       *  @brief measure the filtered two-point correlation function,
       *  \f$w(r_c)=2\pi \int dr \xi(r) W(r,r_c) r^2\f$, where \f$W(x)
       *  = (2x)^2(1-x)(0.5-x)r_c^{-3}\f$
       *  
       *  @param data pointer to an object of type Pair containing the
       *  data-data pairs
       *
       *  @return pointer to an object of type Data
       */
      virtual shared_ptr<data::Data> Filtered (const shared_ptr<data::Data> data)
      { (void)data; ErrorCBL("Error in Filtered() of TwoPointCorrelation.h"); shared_ptr<data::Data> dd; return dd; }
      
      /**
       *  @brief measure the projected two-point correlation function
       *
       *  @param rp projected separation
       *
       *  @param pi line of sight separation
       *
       *  @param xi the 2D two-point correlation function
       *
       *  @param error_xi errors on the 2D two-point correlation
       *  function
       *
       *  @return pointer to an object of type Data
       */
      virtual shared_ptr<data::Data> Projected (const vector<double> rp, const vector<double> pi, const vector<vector<double> > xi, const vector<vector<double> > error_xi)
      { (void)rp; (void)pi; (void)xi; (void)error_xi; ErrorCBL("Error in Projected() of TwoPointCorrelation.h"); shared_ptr<data::Data> data; return data; }

      /**
       *  @brief measure the deprojected two-point correlation
       *  function
       *  
       *  @param rp perpendicular separation
       *
       *  @param ww projected two-point correlation function
       *
       *  @param error_ww error on the projected two-point correlation
       *  function
       *
       *  @return pointer to an object of type Data
       */
      virtual shared_ptr<data::Data> Deprojected (const vector<double> rp, const vector<double> ww, const vector<double> error_ww)
      { (void)rp; (void)ww; (void)error_ww; ErrorCBL("Error in Deprojected() of TwoPointCorrelation.h"); shared_ptr<data::Data> data; return data; }

      /**
       *  @brief measure the multipoles of the two-point correlation
       *  function
       *  
       *  @param rr absolute separation 
       *
       *  @param mu angular separation
       *
       *  @param xi the 2D two-point correlation function
       *
       *  @param error_xi errors on the 2D two-point correlation
       *  function
       *
       *  @return pointer to an object of type Data
       */
      virtual shared_ptr<data::Data> Multipoles(const vector<double> rr, const vector<double> mu, const vector<vector<double>> xi, const vector<vector<double>> error_xi)
      { (void)rr; (void)mu; (void)xi; (void)error_xi; ErrorCBL("Error in Multipoles() of TwoPointCorrelation.h"); shared_ptr<data::Data> data; return data; }

      /**
       *  @brief measure the wedges of the two-poinr correlation
       *  function
       *  
       *  @param rr absolute separation 
       *
       *  @param mu angular separation
       *
       *  @param xi the 2D two-point correlation function
       *
       *  @param error_xi errors on the 2D two-point correlation
       *  function
       *
       *  @return pointer to an object of type Data
       */
      virtual shared_ptr<data::Data> Wedges(const vector<double> rr, const vector<double> mu, const vector<vector<double>> xi, const vector<vector<double> > error_xi)
      { (void)rr; (void)mu; (void)xi; (void)error_xi; ErrorCBL("Error in Wedges() of TwoPointCorrelation.h"); shared_ptr<data::Data> data; return data; }


      ///@}

       /**
       *  @name Member functions to count pairs and compute the
       *  two-point correlation function
       */
      ///@{


      /**
       *  @brief measure the two-point correlation function with
       *  Poisson errors
       *
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to
       *  store the number of pairs (if the pairs are read from files)
       *
       *  @param count_dd true &rarr; count the number of data-data
       *  pairs; false &rarr; read the number of data-data pairs from
       *  file
       *
       *  @param count_rr true &rarr; count the number of
       *  random-random pairs; false &rarr; read the number of
       *  random-random pairs from file
       *
       *  @param count_dr true &rarr; count the number of data-random
       *  pairs; false &rarr; read the number of data-random pairs
       *
       *  @param tcount true &rarr; activate the time counter; false
       *  &rarr; no time counter
       *
       *  @param estimator the estimator used to measure the two-point
       *  correlation function
       *
       *  @return none
       */
      virtual void measurePoisson (const string dir_output_pairs = par::defaultString, const vector<string> dir_input_pairs={}, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_)
      { (void)dir_output_pairs; (void)dir_input_pairs; (void)count_dd; (void)count_rr; (void)count_dr; (void)tcount; (void)estimator; cosmobl::ErrorCBL("Error inmeasurePoisson() of TwoPointCorrelation.h!"); }

      /**
       *  @brief measure the two-point correlation function estimating
       *  the covariance with Jackknife resampling
       *
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to
       *  store the number of pairs (if the pairs are read from files)
       *
       *  @param dir_output_resample output directory used to store
       *  the Jackknife resampling correlation functions, with Poisson
       *  errors
       *
       *  @param count_dd true &rarr; count the number of data-data
       *  pairs; false &rarr; read the number of data-data pairs from
       *  file
       *
       *  @param count_rr true &rarr; count the number of
       *  random-random pairs; false &rarr; read the number of
       *  random-random pairs from file
       *
       *  @param count_dr true &rarr; count the number of data-random
       *  pairs; false &rarr; read the number of data-random pairs
       *
       *  @param tcount true &rarr; activate the time counter; false
       *  &rarr; no time counter
       *
       *  @param estimator the estimator used to measure the two-point
       *  correlation function
       *
       *  @return none
       */
      virtual void measureJackknife (const string dir_output_pairs = par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_resample=par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_)
      { (void)dir_output_pairs; (void)dir_input_pairs; (void)dir_output_resample; (void)count_dd; (void)count_rr; (void)count_dr; (void)tcount; (void)estimator; cosmobl::ErrorCBL("Error in measureJackknife() of TwoPointCorrelation.h!"); }

      /**
       *  @brief measure the two-point correlation function estimating
       *  the covariance with Bootstrap resampling
       *
       *  @param nMocks number of mocks to be generated with bootstrap
       *  resampling
       *
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to
       *  store the number of pairs (if the pairs are read from files)
       *
       *  @param dir_output_resample output directory used to store
       *  the Bootstrap resampling correlation function, with Poisson
       *  errors
       *
       *  @param count_dd true &rarr; count the number of data-data
       *  pairs; false &rarr; read the number of data-data pairs from
       *  file
       *
       *  @param count_rr true &rarr; count the number of
       *  random-random pairs; false &rarr; read the number of
       *  random-random pairs from file
       *
       *  @param count_dr true &rarr; count the number of data-random
       *  pairs; false &rarr; read the number of data-random pairs
       *
       *  @param tcount true &rarr; activate the time counter; false
       *  &rarr; no time counter
       *
       *  @param estimator the estimator used to measure the two-point
       *  correlation function
       *
       *  @return none
       */
      virtual void measureBootstrap (const int nMocks, const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_resample=par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_)
      { (void)nMocks; (void)dir_output_pairs; (void)dir_input_pairs; (void)dir_output_resample; (void)count_dd; (void)count_rr; (void)count_dr; (void)tcount; (void)estimator; cosmobl::ErrorCBL("Error in measureBootstrap() of TwoPointCorrelation.h!"); }

      /**
       *  @brief measure the jackknife resampling of the two-point
       *  correlation function, &xi;(r)
       *
       *  @param dd vector of data-data pairs, divided per regions
       *
       *  @param rr vector of random-random pairs, divided per regions
       *
       *  @return none
       */
      virtual vector<shared_ptr<data::Data> > XiJackknife (const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr)
      { (void)dd; (void)rr; cosmobl::ErrorCBL("Error in vector<shared_ptr<data::Data> > XiJackknife of TwoPointCorrelation.h!"); vector<shared_ptr<data::Data> > data; return data; }

      /**
       *  @brief measure the jackknife resampling of the two-point
       *  correlation function, &xi;(r)
       *
       *  @param dd vector of data-data pairs, divided per regions
       *
       *  @param rr vector of random-random pairs, divided per regions
       *
       *  @param dr vector of random-random pairs, divided per regions   
       *
       *  @return none
       */
      virtual vector<shared_ptr<data::Data> > XiJackknife (const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr)
      { (void)dd; (void)rr; (void)dr; cosmobl::ErrorCBL("Error in vector<shared_ptr<data::Data> > XiJackknife of TwoPointCorrelation.h!"); vector<shared_ptr<data::Data> > data; return data; }

      /**
       *  @brief measure the bootstrap resampling of the two-point
       *  correlation function, &xi;(r)
       *
       *  @param nMocks number of bootstrap resampling
       *
       *  @param dd vector of data-data pairs, divided per regions
       *
       *  @param rr vector of random-random pairs, divided per regions      
       *
       *  @return none
       */
      virtual vector<shared_ptr<data::Data> > XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr)
      { (void)nMocks; (void)dd; (void)rr; cosmobl::ErrorCBL("Error in vector<shared_ptr<data::Data> > XiJackknife of TwoPointCorrelation.h!"); vector<shared_ptr<data::Data> > data; return data; }

      /**
       *  @brief measure the bootstrap resampling of the two-point
       *  correlation function, &xi;(r)
       *
       *  @param nMocks number of bootstrap resampling
       *
       *  @param dd vector of data-data pairs, divided per regions
       *
       *  @param rr vector of random-random pairs, divided per regions 
       *
       *  @param dr vector of random-random pairs, divided per regions  
       *
       *  @return none
       */
      virtual vector<shared_ptr<data::Data> > XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr)
      { (void)nMocks; (void)dd; (void)rr; (void)dr; cosmobl::ErrorCBL("Error in vector<shared_ptr<data::Data> > XiBootstrap of TwoPointCorrelation.h!"); vector<shared_ptr<data::Data> > data; return data; }

      ///@}

    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class TwoPointCorrelation
       */
      TwoPointCorrelation () = default;

      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random object of class Catalogue containing the random
       *  data catalogue
       *  @param compute_extra_info true &rarr; compute extra
       *  information related to the pairs, such as the mean pair
       *  separation and redshift
       *  @return object of class TwoPointCorrelation
       */
      TwoPointCorrelation (const catalogue::Catalogue data, const catalogue::Catalogue random, const bool compute_extra_info=false) 
	: m_data(make_shared<catalogue::Catalogue>(catalogue::Catalogue(move(data)))), m_random(make_shared<catalogue::Catalogue>(catalogue::Catalogue(move(random)))), m_compute_extra_info(compute_extra_info) {}
    
      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~TwoPointCorrelation () = default;

      /**
       *  @brief static factory used to construct two-point
       *  correlation functions of any type
       *
       *  @param type the type of two-point correlation function; it
       *  can be: \_1D_monopole\_, \_1D_angular\_
       *
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param binType binning type: 0 &rarr; linear; 1 &rarr;
       *  logarithmic
       *  @param Min minimum separation used to count the pairs
       *  @param Max maximum separation used to count the pairs
       *  @param nbins number of bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @param compute_extra_info true &rarr; compute extra
       *  information related to the pairs, such as the mean pair
       *  separation and redshift
       *  @return a pointer to an object of class TwoPointCorrelation of
       *  a given type
       */
      static shared_ptr<TwoPointCorrelation> Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const binType binType, const double Min, const double Max, const int nbins, const double shift, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false);

      /**
       *  @brief static factory used to construct two-point
       *  correlation functions of any type
       *
       *  @param type the type of two-point correlation function; it
       *  can be: \_1D_monopole\_, \_1D_angular\_
       *
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param binType binning type: 0 &rarr; linear; 1 &rarr;
       *  logarithmic
       *  @param Min minimum separation used to count the pairs
       *  @param Max maximum separation used to count the pairs
       *  @param binSize the bin size
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @param compute_extra_info true &rarr; compute extra
       *  information related to the pairs, such as the mean pair
       *  separation and redshift
       *  @return a pointer to an object of class TwoPointCorrelation of
       *  a given type
       */
      static shared_ptr<TwoPointCorrelation> Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const binType binType, const double Min, const double Max, const double binSize, const double shift, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false);

       /**
       *  @brief static factory used to construct two-point
       *  correlation functions of any type
       *
       *  @param type the type of two-point correlation function; it
       *  can be: \_2D_Cartesian\_, \_2D_polar\_
       *
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *
       *  @param binType_D1 binning type in the first dimension: 0
       *  &rarr; linear; 1 &rarr; logarithmic
       *
       *  @param Min_D1 minimum separation in the first dimensionused
       *  to count the pairs
       *
       *  @param Max_D1 maximum separation in the first dimension used
       *  to count the pairs
       *
       *  @param nbins_D1 number of bins in the first dimension
       *
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *
       *  @param binType_D2 binning type in the second dimension: 0
       *  &rarr; linear; 1 &rarr; logarithmic
       *
       *  @param Min_D2 minimum separation in the second dimensionused
       *  to count the pairs
       *
       *  @param Max_D2 maximum separation in the second dimension used
       *  to count the pairs
       *
       *  @param nbins_D2 number of bins in the second dimension
       *
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       * 
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *
       *  @param compute_extra_info true &rarr; compute extra
       *  information related to the pairs, such as the mean pair
       *  separation and redshift
       *
       *  @return a pointer to an object of class TwoPointCorrelation of
       *  a given type
       */
      static shared_ptr<TwoPointCorrelation> Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const binType binType_D1, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const binType binType_D2, const double Min_D2, const double Max_D2, const int nbins_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false);

      /**
       *  @brief static factory used to construct two-point
       *  correlation functions of any type
       *
       *  @param type the type of two-point correlation function; it
       *  can be: \_2D_Cartesian\_, \_2D_polar\_
       *
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *
       *  @param binType_D1 binning type in the first dimension: 0
       *  &rarr; linear; 1 &rarr; logarithmic
       *
       *  @param Min_D1 minimum separation in the first dimension used
       *  to count the pairs
       *
       *  @param Max_D1 maximum separation in the first dimension used
       *  to count the pairs
       *
       *  @param binSize_D1 the bin size in the first dimension
       *
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *
       *  @param binType_D2 binning type in the second dimension: 0
       *  &rarr; linear; 1 &rarr; logarithmic
       *
       *  @param Min_D2 minimum separation in the second dimension used
       *  to count the pairs
       *
       *  @param Max_D2 maximum separation in the second dimension used
       *  to count the pairs
       *
       *  @param binSize_D2 the bin size in the second dimension
       *
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *
       *  @param compute_extra_info true &rarr; compute extra
       *  information related to the pairs, such as the mean pair
       *  separation and redshift
       *
       *  @return a pointer to an object of class TwoPointCorrelation of
       *  a given type
       */
      static shared_ptr<TwoPointCorrelation> Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const binType binType_D1, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const binType binType_D2, const double Min_D2, const double Max_D2, const double binSize_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false);
      
       /**
       *  @brief static factory used to construct two-point
       *  correlation functions of any type
       *
       *  @param type the type of two-point correlation function; it
       *  can be: \_1D_projected\_, \_1D_deprojected\_,
       *  \_1D_multipoles\_
       *
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *
       *  @param binType_D1 binning type in the first dimension: 0
       *  &rarr; linear; 1 &rarr; logarithmic
       *
       *  @param Min_D1 minimum separation in the first dimensionused
       *  to count the pairs
       *
       *  @param Max_D1 maximum separation in the first dimension used
       *  to count the pairs
       *
       *  @param nbins_D1 number of bins in the first dimension
       *
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *
       *  @param Min_D2 minimum separation in the second dimensionused
       *  to count the pairs
       *
       *  @param Max_D2 maximum separation in the second dimension used
       *  to count the pairs
       *
       *  @param nbins_D2 number of bins in the second dimension
       *
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *
       *  @param piMax_integral upper limits of the integral
       *
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *
       *  @param compute_extra_info true &rarr; compute extra
       *  information related to the pairs, such as the mean pair
       *  separation and redshift
       *
       *  @return a pointer to an object of class TwoPointCorrelation of
       *  a given type
       */
      static shared_ptr<TwoPointCorrelation> Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const binType binType_D1, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const double Min_D2, const double Max_D2, const int nbins_D2, const double shift_D2, const double piMax_integral=50., const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false);

      /**
       *  @brief static factory used to construct two-point
       *  correlation functions of any type
       *
       *  @param type the type of two-point correlation function; it
       *  can be: \_1D_projected\_, \_1D_deprojected\_,
       *  \_1D_multipoles\_
       *
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *
       *  @param binType_D1 binning type in the first dimension: 0
       *  &rarr; linear; 1 &rarr; logarithmic
       *
       *  @param Min_D1 minimum separation in the first dimension used
       *  to count the pairs
       *
       *  @param Max_D1 maximum separation in the first dimension used
       *  to count the pairs
       *
       *  @param binSize_D1 the bin size in the first dimension
       *
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *
       *  @param Min_D2 minimum separation in the second dimension used
       *  to count the pairs
       *
       *  @param Max_D2 maximum separation in the second dimension used
       *  to count the pairs
       *
       *  @param binSize_D2 the bin size in the second dimension
       *
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *
       *  @param piMax_integral upper limits of the integral
       *
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *
       *  @param compute_extra_info true &rarr; compute extra
       *  information related to the pairs, such as the mean pair
       *  separation and redshift
       *
       *  @return a pointer to an object of class TwoPointCorrelation of
       *  a given type
       */
      static shared_ptr<TwoPointCorrelation> Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const binType binType_D1, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const double Min_D2, const double Max_D2, const double binSize_D2, const double shift_D2, const double piMax_integral=50., const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false);

      ///@}

    
      /**
       *  @name Member functions to get the private/protected members
       */
      ///@{
      
      /**
       *  @brief get the protected member m_twoPType
       *  @return the two-point correlation function type
       */
      TwoPType twoPType () const { return m_twoPType; }
      
      /**
       *  @brief get the protected member dataset
       *  @return a shared pointer to the dataset
       */
      virtual shared_ptr<data::Data> dataset () const { return m_dataset; }
      
      /**
       *  @brief get the protected member m_data
       *  @return the input data catalogue
       */
      shared_ptr<catalogue::Catalogue> data () const { return m_data; }
   
      /**
       *  @brief get the protected member m_random
       *  @return the input random catalogue
       */
      shared_ptr<catalogue::Catalogue> random () const { return m_random; }

      /**
       *  @brief get the protected member m_dd
       *  @return the number of data-data pairs
       */
      shared_ptr<pairs::Pair> dd () const { return m_dd; }

      /**
       *  @brief get the protected member m_rr
       *  @return the number of random-random pairs
       */
      shared_ptr<pairs::Pair> rr () const { return m_rr; }

      /**
       *  @brief get the protected member m_dr
       *  @return the number of data-random pairs
       */
      shared_ptr<pairs::Pair> dr () const { return m_dr; }

      /**
       *  @brief get the protected member m_dr
       *  @return the number of data-random pairs
       */
      bool compute_extra_info () const { return m_compute_extra_info; }

      /**
       *  @brief get the x coordinates
       *  @return the x coordinates
       */
      virtual vector<double> xx () const = 0;

      /**
       *  @brief get the y coordinates
       *  @return the y coordinates
       */
      virtual vector<double> yy () const { cosmobl::ErrorCBL("Error in yy() of TwoPointCorrelation.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the the binned correlation function 
       *  @return the binned correlation function 
       */
      virtual vector<double> xi1D () const { cosmobl::ErrorCBL("Error in xi() of TwoPointCorrelation.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the error on the binned correlation function
       *  function
       *  @return the error on the binned correlation function
       *  function
       */
      virtual vector<double> error1D () const { cosmobl::ErrorCBL("Error in error() of TwoPointCorrelation.h!"); vector<double> vv; return vv; }
      
      /**
       *  @brief get the the binned correlation function 
       *  @return the binned correlation function 
       */
      virtual vector<vector<double> > xi2D () const { cosmobl::ErrorCBL("Error in xi() of TwoPointCorrelation.h!"); vector<vector<double> > vv; return vv; }

      /**
       *  @brief get the error on the binned correlation function
       *  @return the error on the binned correlation function
       */
      virtual vector<vector<double> > error2D () const { cosmobl::ErrorCBL("Error in error() of TwoPointCorrelation.h!"); vector<vector<double> > vv; return vv; }

      /**
       *  @brief get the monopole of the polar xi
       *  @return the xiMonopole
       */
      virtual vector<double> xiMonopole () const 
      { cosmobl::ErrorCBL("Error in xiMonopole() of TwoPointCorrelation.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the error on the monopole of the polar xi
       *  @return the xiMonopole
       */
      virtual vector<double> errorMonopole () const 
      { cosmobl::ErrorCBL("Error in errorMonopole() of TwoPointCorrelation.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the quadrupole of the polar xi
       *  @return the xiQuadrupole
       */
      virtual vector<double> xiQuadrupole () const 
      { cosmobl::ErrorCBL("Error in xiQuadrupole() of TwoPointCorrelation.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the error on the quadrupole of the polar xi
       *  @return the xiQuadrupole
       */
      virtual vector<double> errorQuadrupole () const 
      { cosmobl::ErrorCBL("Error in errorQuadrupole() of TwoPointCorrelation.h!"); vector<double> vv; return vv; }
      /**
       *  @brief get the octupole of the polar xi
       *  @return the xiOctupole
       */
      virtual vector<double> xiOctupole () const 
      { cosmobl::ErrorCBL("Error in xiOctupole() of TwoPointCorrelation.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the error on the octupole of the polar xi
       *  @return the error on xiOctupole
       */
      virtual vector<double> errorOctupole () const 
      { cosmobl::ErrorCBL("Error in xiOctupole() of TwoPointCorrelation.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the perpendicular wedge of the polar xi
       *  @return the perpendicular wedge of the polar xi
       */
      virtual vector<double> xiPerpendicular () const 
      { cosmobl::ErrorCBL("Error in xiPerpendicular() of TwoPointCorrelation.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the perpendicular wedge of the polar xi
       *  @return the error on the perpendicular wedge of the polar xi
       */
      virtual vector<double> errorPerpendicular () const 
      { cosmobl::ErrorCBL("Error in errorPerpendicular() of TwoPointCorrelation.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the parallel wedge of the polar xi
       *  @return the xiPar coordinates
       */
      virtual vector<double> xiParallel () const 
      { cosmobl::ErrorCBL("Error in xiParallel() of TwoPointCorrelation.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the error on the parallel wedge of the polar xi
       *  @return the error on the parallel wedge of the polar xi
       */
      virtual vector<double> errorParallel () const 
      { cosmobl::ErrorCBL("Error in errorParallel() of TwoPointCorrelation.h!"); vector<double> vv; return vv; }

      ///@}
      

      /**
       *  @name Member functions to set protected members
       */
      ///@{
    
      /**
       *  @brief add a data catalogue
       *  @param data object of class Catalogue 
       *  @return none
       */
      void set_data (const catalogue::Catalogue data) { m_data = make_shared<catalogue::Catalogue>(catalogue::Catalogue(move(data))); }

      /**
       *  @brief add a random catalogue
       *  @param random object of class Catalogue 
       *  @return none
       */
      void set_random (const catalogue::Catalogue random) { m_random = make_shared<catalogue::Catalogue>(catalogue::Catalogue(move(random))); }

      ///@}

    
      /**
       *  @name Member functions to count measure the two-point correlation function
       */
      ///@{

      /**
       *  @brief measure the two-point correlation function
       *
       *  @param errorType type of error
       *  
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to
       *  store the number of pairs (if the pairs are read from files)
       *
       *  @param dir_output_resample output directory of the
       *  resampled correlation function
       *
       *  @param nMocks number of resampling used for bootstrap
       *
       *  @param count_dd true &rarr; count the number of data-data
       *  pairs; false &rarr; read the number of data-random pairs from
       *  file
       *
       *  @param count_rr true &rarr; count the number of random-random
       *  pairs; false &rarr; read the number of random-random pairs
       *
       *  @param count_dr true &rarr; count the number of data-random
       *  pairs; false &rarr; read the number of data-random pairs
       *
       *  @param tcount true &rarr; activate the time counter; false
       *  &rarr; no time counter
       *
       *  @param estimator the estimator used to measure the two-point
       *  correlation function
       *
       *  @return none
       */
      virtual void measure (const ErrorType errorType=ErrorType::_Poisson_, const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_resample=par::defaultString, const int nMocks=0, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_) = 0;

      ///@}

  
      /**
       *  @name Input/Output member functions (customized in all the derived classes)
       */
      ///@{

      /**
       *  @brief read the measured two-point correlation
       *  @param dir input directory
       *  @param file input file
       *  @return none
       */
      virtual void read (const string dir, const string file) = 0;

      /**
       *  @brief write the measured two-point correlation
       *  @param dir output directory
       *  @param file output file
       *  @param rank cpu index (for MPI usage)
       *  @return none
       */
      virtual void write (const string dir=par::defaultString, const string file=par::defaultString, const int rank=0) const = 0;

      /**
       *  @brief write the measured two-point correlation
       *  @param dir output directory
       *  @param file output file
       *  @param full false &rarr; simply store the data; true &rarr;
       *  duplicate the data in the other three quadrands (usefull
       *  e.g. when storing the 2D correlation function)
       *  @param rank cpu index (for MPI usage)
       *  @return none
       */
      virtual void write (const string dir, const string file, const bool full, const int rank=0) const
      { (void)dir; (void)file; (void)full; (void)rank; cosmobl::ErrorCBL("Error in write() of TwoPointCorrelation.h!"); }
      
      ///@}
      

      /**
       *  @name Member functions to compute, read and write covariance matrix
       */
      ///@{ 

      /**
       *  @brief read the measured covariance matrix
       *  @param dir input directory
       *  @param file input file
       *  @return none
       */
      virtual void read_covariance_matrix (const string dir, const string file)
      { (void)dir; (void)file; cosmobl::ErrorCBL("Error in read_covariance_matrix() of TwoPointCorrelation.h!"); }

      /**
       *  @brief write the measured two-point correlation
       *  @param dir output directory
       *  @param file output file
       *  @return none
       */
      virtual void write_covariance_matrix (const string dir, const string file) const
      { (void)dir; (void)file; cosmobl::ErrorCBL("Error in write_covariance_matrix() of TwoPointCorrelation.h!");}

      /**
       *  @brief compute the covariance matrix
       *  @param xi_collection vector containing the xi to compute the covariance matrix
       *  @param doJK 1 &rarr; compute jackknife covariance matrix; 0 compute standard covariance matrix
       *  @return none
       */
      virtual void compute_covariance_matrix (vector<shared_ptr<data::Data>> xi_collection, bool doJK)
      { (void)xi_collection; (void)doJK; cosmobl::ErrorCBL("Error in compute_covariance_matrix() of TwoPointCorrelation.h!");}

      /**
       *  @brief compute the covariance matrix
       *  @param file_xi vector containing the path to the xi to compute the covariance matrix
       *  @param doJK 1 &rarr; compute jackknife covariance matrix; 0 compute standard covariance matrix
       *  @return none
       */
      virtual void compute_covariance_matrix (vector<string> file_xi, bool doJK)
      { (void)file_xi; (void)doJK; cosmobl::ErrorCBL("Error in compute_covariance_matrix() of TwoPointCorrelation.h!");}

      ///@} 

      /**
       *  @name Member functions to estimate the errors
       */
      ///@{
    
      /**
       *  @brief estimate the Poisson error 
       *  @param dd number of data-data pairs
       *  @param rr number of random-random pairs
       *  @param dr number of data-random pairs
       *  @param nData number of data points
       *  @param nRandom number of random points
       *  @return the Poisson error
       */
      double PoissonError (const double dd, const double rr, const double dr, const int nData, const int nRandom) const;
    
      ///@}

    };
  }
}

#endif
