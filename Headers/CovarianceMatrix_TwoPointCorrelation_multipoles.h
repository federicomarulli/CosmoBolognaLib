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
 *  @file Headers/CovarianceMatrix_TwoPointCorrelation_multipoles.h
 *
 *  @brief The class TwoPointCorrelation
 *
 *  This file defines the interface of the class CovarianceMatrix_TwoPointCorrelation_multipoles,
 *  used to measure the covariance of two-point correlation function
 *
 *  @author Alfonso Veropalumbo
 *
 *  @author alfonso.veropalumbo@uniroma3.it
 */

#ifndef __COVMAT_TWOPOINT_MULTI__
#define __COVMAT_TWOPOINT_MULTI__


#include "Data.h"
#include "CovarianceMatrix.h"
#include "ChainMesh_Catalogue.h"


// ===================================================================================================


namespace cbl {

  namespace measure {

    /**
     *  @brief The namespace of the <B> covariance matrix
     *  </B> measure
     *
     *  The \e measure::cov namespace contains all the functions and
     *  classes to measure cosmological covariance matrices
     */
    namespace covmat {

      /**
       *  @class CovarianceMatrix_TwoPointCorrelation_multipoles CovarianceMatrix_TwoPointCorrelation_multipoles.h
       *  "Headers/CovarianceMatrix_TwoPointCorrelation_multipoles.h"
       *
       *  @brief The class CovarianceMatrix_TwoPointCorrelation_multipoles
       *
       *  This is the class used to measure the covariance of two-point
       *  correlation function monopole
       *
       */
      class CovarianceMatrix_TwoPointCorrelation_multipoles {

        protected :

        /**
         *  @name Random catalogue
         */
        ///@{

        /// random catalogue
        std::shared_ptr<catalogue::Catalogue> m_random;

        ///@}

        /**
         *  @name Chain mesh
         */
        ///@{

        /// First Chain-mesh
        chainmesh::ChainMesh_Catalogue m_chainMesh_rMax;

        /// second Chain-mesh
        chainmesh::ChainMesh_Catalogue m_chainMesh_rCut;

        ///@}

        /**
        * @name twopcf parameters
        */

        ///@{

        /// Minimum scale
        double m_rMin;

        /// Maximum scale
        double m_rMax;

        /// Bin size
        double m_binSize;

        /// Number of bins
        size_t m_nBins;

        /// Bin type
        BinType m_binType;

        ///@}

        ///@{

        std::shared_ptr<cbl::glob::FuncGrid> m_interpXi_0;

        std::shared_ptr<cbl::glob::FuncGrid> m_interpXi_2;

        std::shared_ptr<cbl::glob::FuncGrid> m_interpXi_4;

        /// Minimum separation, below this scale two point
        /// correlation function will be set to 0
        double m_minSep;

        /// Maximum separation, above this scale two point
        /// correlation function will be set to 0
        double m_maxSep;

        ///@}

        /**
         *  @name Random catalogue
         */
        ///@{

        /// covariance matrix
        Eigen::MatrixXd m_random_random;

        /// First term of the 2pcf covariance model
        Eigen::MatrixXd m_C4;

        /// Second term of the 2pcf covariance model
        Eigen::MatrixXd m_C3;

        /// Third term of the 2pcf covariance model
        Eigen::MatrixXd m_C2;

        ///@}

        /**
         * @brief compute random-random pairs
         * and fill the m_random_random Matrix
         * @param tcount	 true &rarr; activate the time counter; false
         *  &rarr; no time counter
         * @return None
         */
        void m_compute_RR (const bool tcount);

        /**
         * @brief compute C4 term
         * and fill the m_C4 Matrix
         * @param tcount	 true &rarr; activate the time counter; false
         *  &rarr; no time counter
         * @return None
         */
        void m_compute_C4 (const bool tcount);

        /**
         * @brief compute C3 term
         * and fill the m_C3 Matrix
         * @param tcount	 true &rarr; activate the time counter; false
         *  &rarr; no time counter
         * @return None
         */
        void m_compute_C3 (const bool tcount);

        /**
         * @brief compute C2 term
         * and fill the m_C2 Matrix
         * @param tcount	 true &rarr; activate the time counter; false
         *  &rarr; no time counter
         * @return None
         */
        void m_compute_C2 (const bool tcount);

        /**
         * @brief write the Matrix given in input
         * @param the matrix to be wrote
         * @param dir the output directory
         * @param file the output file
         * @param the number of digits
         * @return None
         */
        void write_matrix (Eigen::MatrixXd matrix, const std::string dir, const std::string file, const int nDigits=20);

        /**
         * @brief set catalogue and chain mesh
         * @param random random catalogue
         * @return None
         */
        void m_set_catalogues(const catalogue::Catalogue random);

      public:

        /**
         * @brief default constructor
         * @return object of type CovarianceMatrix_TwoPointCorrelation_multipoles
         */
        CovarianceMatrix_TwoPointCorrelation_multipoles () {}

       /**
        * @brief Constructor of CovarianceMatrix_TwoPointCorrelation_multipoles
        * @param random random catalogue
        * @param binType the bin binType
        * @param rMin minimum separation
        * @param rMax maximum separation
        * @param nBins number of bins
        * @param interpXi_0 function to interpolate the xi monopole
        * @param interpXi_2 function to interpolate the xi quadrupole
        * @param interpXi_4 function to interpolate the xi hexadecapole
        * @param minSeparation  Minimum separation, below this scale two point
        * correlation function will be set to 0
        * @param minSeparation  Maximum separation, above this scale two point
        * correlation function will be set to 0
        * @return object of type CovarianceMatrix_TwoPointCorrelation_multipoles
        */
        CovarianceMatrix_TwoPointCorrelation_multipoles (const catalogue::Catalogue random,
          const BinType binType,
          const double rMin,
          const double rMax,
          const int nBins,
          const glob::FuncGrid interpXi_0,
          const glob::FuncGrid interpXi_2,
          const glob::FuncGrid interpXi_4,
          const double minSeparation=0.05,
          const double maxSeparation=200);

       /**
        * @brief Constructor of CovarianceMatrix_TwoPointCorrelation_multipoles
        * @param random random catalogue
        * @param binType the bin binType
        * @param rMin minimum separation
        * @param rMax maximum separation
        * @param binSize the bin size
        * @param interpXi_0 function to interpolate the xi monopole
        * @param interpXi_2 function to interpolate the xi quadrupole
        * @param interpXi_4 function to interpolate the xi hexadecapole
	* @param minSeparation  Minimum separation, below this scale two point
        * correlation function will be set to 0
        * @param minSeparation  Maximum separation, above this scale two point
        * correlation function will be set to 0
        * @return object of type CovarianceMatrix_TwoPointCorrelation_multipoles
        */
        CovarianceMatrix_TwoPointCorrelation_multipoles (const catalogue::Catalogue random,
          const BinType binType,
          const double rMin,
          const double rMax,
          const double binSize,
          const glob::FuncGrid interpXi_0,
          const glob::FuncGrid interpXi_2,
          const glob::FuncGrid interpXi_4,
          const double minSeparation=0.05,
          const double maxSeparation=200);

       /**
        * @brief default destructor
        * @return None
        */
        virtual ~CovarianceMatrix_TwoPointCorrelation_multipoles ()  = default;

       /**
        * @brief compute the terms
        * @param output_dir directory to store the output
        * @param compute_RR true &rarr; count the number of random-random
     	  *  pairs; false &rarr; read the number of random-random pairs
        * @param compute_C4  true &rarr; compute the C4 term; false &rarr; read the C4 term
        * @param compute_C3  true &rarr; compute the C3 term; false &rarr; read the C3 term
        * @param compute_C2  true &rarr; compute the C2 term; false &rarr; read the C2 term
        * @param tcount	 true &rarr; activate the time counter; false
        *  &rarr; no time counter
        * @return object of type CovarianceMatrix_TwoPointCorrelation_multipoles
        */
        void compute_terms (const std::string output_dir,
                            const bool compute_RR=true,
                            const bool compute_C4=true,
                            const bool compute_C3=true,
                            const bool compute_C2=true,
                            const bool tcount=false);

        /**
         * @brief return the covariance matrix
         * @brief alpha parameter to contorl non gaussianity
         * @return object of type cbl::data::CovarianceMatrix
         */
        data::CovarianceMatrix operator() (const double alpha);
      };
    }
  }
}

#endif
