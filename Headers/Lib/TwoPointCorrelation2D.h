/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/Lib/TwoPointCorrelation2D.h
 *
 *  @brief The class TwoPointCorrelation2D
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation2D, used to measure the 2D two-point
 *  correlation function in Cartesian and polar coordinates
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINT2D__
#define __TWOPOINT2D__


#include "TwoPointCorrelation.h"


// ===================================================================================================


namespace cosmobl {

  namespace twopt {
    
    /**
     *  @class TwoPointCorrelation2D TwoPointCorrelation2D.h
     *  "Headers/Lib/TwoPointCorrelation2D.h"
     *
     *  @brief The class TwoPointCorrelation2D
     *
     *  This class is used to handle objects of type <EM>
     *  TwoPointCorrelation2D </EM>. It is used to measure the 2D
     *  two-point correlation function, in Cartesian and polar
     *  coordinates, \f$\xi(r_p,\pi)\f$ or \f$\xi(r,\mu)\f$, that is as
     *  a function of perpendicular, \f$r_p\f$, and parallel, \f$\pi\f$,
     *  line-of-sight separations, and as a function of absolute
     *  separation, \f$r=\sqrt{r_p^2+\pi^2}\f$, and the cosine of the
     *  angle between the separation vector and the line of sight,
     *  \f$\mu\equiv\cos\theta=s_\parallel/s\f$, respectively.
     */
    class TwoPointCorrelation2D : public TwoPointCorrelation {

    protected :
    
      /**
       *  @name Internal input/output methods
       */
      ///@{
    
      /**
       *  @brief write the number of pairs
       *  @param PP pointer to an object of class Pair
       *  @param dir output directory
       *  @param file output file
       *  @return none
       */
      void write_pairs (const shared_ptr<pairs::Pair> PP, const string dir, const string file) const override;

      /**
       *  @brief read the number of pairs
       *  @param [out] PP pointer to an object of class Pair
       *  @param [in] dir vector of input directories
       *  @param [in] file input file
       *  @return none
       */
      void read_pairs (shared_ptr<pairs::Pair> PP, const vector<string> dir, const string file) const override;

      /**
       *  @brief write the number of pairs
       *  @param PP pointer to a vector of objects of class Pair
       *  @param dir output directory
       *  @param file output file
       *  @return none
       */
      void write_pairs (const vector<shared_ptr<pairs::Pair> >  PP, const string dir, const string file) const override;

      /**
       *  @brief read the number of pairs
       *  @param [out] PP pointer to a vector of objects of class Pair
       *  @param [in] dir vector of input directories
       *  @param [in] file input file
       *  @return none
       */
      void read_pairs (vector<shared_ptr<pairs::Pair> > PP, const vector<string> dir, const string file) const override;

      ///@}

      /**
       *  @name Member functions to measure the two-point correlation function
       */
      ///@{

      /**
       *  @brief measure the xi with Poisson error using measured pairs
       *  using the Natural Estimator
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
      shared_ptr<Data> NaturalEstimatorTwoP (const shared_ptr<pairs::Pair> dd, const shared_ptr<pairs::Pair> rr, const int nData, const int nRandom) override;

      /**
       *  @brief measure the xi with Poisson error using measured pairs
       *  using the Landy-Szalay estimator
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
      shared_ptr<Data> LandySzalayEstimatorTwoP (const shared_ptr<pairs::Pair> dd, const shared_ptr<pairs::Pair> rr, const shared_ptr<pairs::Pair> dr, const int nData, const int nRandom) override;

      /**
       *  @brief measure the jackknife resampling of the two-point correlation
       *  function, &xi;(r) 
       *
       *  @param dd vector of data-data pairs, divider per regions
       *
       *  @param rr vector of random-random pairs, divider per regions
       *
       *  @return none
       */
      vector<shared_ptr<Data> > XiJackknife(const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr) override;

      /**
       *  @brief measure the jackknife resampling of the two-point correlation
       *  function, &xi;(r)         
       *
       *  @param dd vector of data-data pairs, divider per regions
       *
       *  @param rr vector of random-random pairs, divider per regions
       *
       *  @param dr vector of random-random pairs, divider per regions   *
       *
       *  @return none
       */
      vector<shared_ptr<Data> > XiJackknife(const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr) override;

      /**
       *  @brief measure the bootstrap resampling of the two-point correlation
       *  function, &xi;(r)  
       *
       *  @param nMocks number of bootstrap resampling
       *
       *  @param dd vector of data-data pairs, divider per regions
       *
       *  @param rr vector of random-random pairs, divider per regions      
       *
       *  @return none
       */
      vector<shared_ptr<Data> > XiBootstrap(const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr) override;

      /**
       *  @brief measure the bootstrap resampling of the two-point correlation
       *  function, &xi;(r)  
       *
       *  @param nMocks number of bootstrap resampling
       *
       *  @param dd vector of data-data pairs, divider per regions
       *
       *  @param rr vector of random-random pairs, divider per regions 
       *
       *  @param dr vector of random-random pairs, divider per regions  
       *
       *  @return none
       */
      vector<shared_ptr<Data> > XiBootstrap(const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr) override;
      
      ///@}


    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class TwoPointCorrelation2D
       */
      TwoPointCorrelation2D () { m_dataset = Data::Create(dataType::_2D_data_); }

      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @return object of class TwoPointCorrelation2D
       */
      TwoPointCorrelation2D (const catalogue::Catalogue data, const catalogue::Catalogue random) 
	: TwoPointCorrelation(data, random) { m_dataset = Data::Create(dataType::_2D_data_); }

      /**
       *  @brief default destructor
       *  @return none
       */
      ~TwoPointCorrelation2D () = default;

      ///@}

        /**
       *  @name Member functions to get the private/protected members
       */
      ///@{

      /**
       *  @brief get the protected member Data::m_x
       *  @return the x coordinates
       */
      vector<double> xx () const override { return m_dataset->xx(); }

      /**
       *  @brief get the protected member m_y
       *  @return the y coordinates
       */
      vector<double> yy () const override { return m_dataset->yy(); }

      /**
       *  @brief get the protected member m_fxy
       *  @return the binned correlation function 
       */
      vector<vector<double> > xi2D () const override { return m_dataset->fxy(); }

      /**
       *  @brief get the protected member m_error_fxy
       *  @return the error on the binned correlation function
       *  function
       */
      vector<vector<double> > error2D () const override { return m_dataset->error_fxy(); }
      
      ///@}


      /**
       *  @name Member functions to count measure the two-point correlation function
       */
      ///@{

      /**
       *  @brief measure the two-point correlation function
       *
       *  @param errType type of &xi;(r) error
       *  
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to
       *  store the number of pairs (if the pairs are read from files)
       *
       *  @param dir_output_ResampleXi output directory of the
       *  resampled &xi;(r)
       *
       *  @param nMocks number of resampling for bootstrap
       *
       *  @param count_dd 1 &rarr; count the number of data-data
       *  opairs; 0 &rarr; read the number of data-data pairs from
       *  file
       *
       *  @param count_rr 1 &rarr; count the number of random-random
       *  opairs; 0 &rarr; read the number of random-random pairs from
       *  file
       *
       *  @param count_dd 1 &rarr; count the number of data-random
       *  opairs; 0 &rarr; read the number of data-random pairs from
       *  file
       *
       *  @param count_rr 1 &rarr; count the number of random-random
       *  pairs; 0 &rarr; read the number of random-random pairs
       *
       *  @param count_dr 1 &rarr; count the number of data-random
       *  pairs; 0 &rarr; read the number of data-random pairs
       *
       *  @param tcount 1 &rarr; activate the time counter; 0 &rarr;
       *  don't activate the time counter; 
       *
       *  @return none
       */
      virtual void measure (const ErrorType errType=ErrorType::_Poisson_, const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_ResampleXi=par::defaultString, const int nMocks=0, const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1) = 0;
      
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
      virtual void read (const string dir, const string file)
      { ErrorMsg("Error in read() of TwoPointCorrelation2D.h"); }	

      /**
       *  @brief write the measured two-point correlation
       *  @param dir output directory
       *  @param file output file
       *  @param rank cpu index (for MPI usage)
       *  @return none
       */
      virtual void write (const string dir, const string file, const int rank=0) const
      { ErrorMsg("Error in write() of TwoPointCorrelation2D.h"); }	

      ///@}

    };
  }
}

#endif
