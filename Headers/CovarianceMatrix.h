/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/CovarianceMatrix.h
 *
 *  @brief The class CovarianceMatrix
 *
 *  This file defines the interface of the class CovarianceMatrix
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __COVMAT__
#define __COVMAT__

#include "Func.h"

namespace cbl {

  /**
   *  @brief The namespace of the functions and classes used to handle
   *  <B> data </B> of any kind
   *  
   *  The \e data namespace contains all the main functions and
   *  classes to handle data of any kind
   */
  namespace data {

     
    /**
     *  @class CovarianceMatrix CovarianceMatrix.h
     *  "Headers/CovarianceMatrix.h"
     *
     *  @brief The class CovarianceMatrix
     *
     *  This is the base class used to manage 
     *  covariance matrices
     */
    class CovarianceMatrix
    {
      
    protected:

      /// number of data
      size_t m_order;

      /// covariance matrix
      Eigen::MatrixXd m_matrix;
      
      /// precision matrix
      Eigen::MatrixXd m_precision;
      
      /// correlation matrix
      Eigen::MatrixXd m_correlation;

      /// diagonal of the covariance matrix
      Eigen::VectorXd m_variance;

      /// standard deviation
      Eigen::VectorXd m_std;

      /// The hartlap factor, only set when
      // the covariance is measured from multiple dataset
      double m_hartlap_factor;

      /**
       * @brief set internal attributes to
       * default values
       */
      void m_set_default ();

      /**
       * @brief set internal attributes
       *
       * @param matrix the covariance matrix
       *
       * @param nmeasures number of measures
       *
       * @param prec the precision required in the inversion of the
       * covariance matrix
       *
       * @return none, or an error message if the derived object does
       * not have this member
       */
      virtual void m_set (const std::vector<double> matrix, const double nmeasures=-1, const double prec=1.e-10);

      /**
       * @brief compute the hartlap 
       * factor. This is used to de-bias
       * precision matrix measured from covariance
       * measured from limited number of datasets
       *
       *  \f[
       *  \hat{\Psi}=\left(1-\frac{N_{\mathrm{b}}+1}{N_{\mathrm{s}}-1}\right) \hat{\mathrm{c}}^{-1}
       *  \f]
       *
       * @param order the matrix order (or the number of data points)
       *
       * @param nmeasures number of measures used to compute the covariance
       *
       * @return the hartlap factor
       */
      inline double hartlap_factor (const double order, const double nmeasures=-1) { return (nmeasures>1) ? 1-(order+1)/(nmeasures-1) : 1; }

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *
       *  
       */
      CovarianceMatrix () { m_set_default(); }
      
      /**
       *  @brief constructor which sets the
       *  covariance matrix
       *
       *  @param covariance_matrix array containing the covariance matrix
       *
       *  @param nmeasures number of measures used to compute the covariance
       *
       *  @param prec the precision required in the inversion of the
       *  covariance matrix
       *
       *  
       */
      CovarianceMatrix (std::vector<std::vector<double>> covariance_matrix, const double nmeasures=-1, const double prec=1.e-10) 
      { set_from_matrix(covariance_matrix, nmeasures, prec); }

      /**
       *  @brief constructor which gets the data from an input vector
       *
       *  @param standard_deviation the standard deviation
       *
       *  @param nmeasures number of measures used to compute the covariance
       *
       *  
       */
      CovarianceMatrix (std::vector<double> standard_deviation, const double nmeasures=-1)
      {set_from_standard_deviation(standard_deviation, nmeasures);}

      /**
       *  @brief constructs with sets the covariance
       *  matrix reading from an input file
       *
       *  @param filename file containing the covariance matrix in the
       *  format: column 0 \f$ \rightarrow \f$ x<SUB>i</SUB>, column 1
       *  \f$ \rightarrow \f$ x<SUB>j</SUB>, column cov_col &rarr;
       *  cov(x<SUB>i</SUB>,x<SUB>j</SUB>)
       *
       *  @param cov_col covariance matrix column, starting from 0
       *
       *  @param skipped_lines comment lines to be skipped
       *
       *  @param nmeasures number of measures used to compute the covariance
       *
       *  @param prec the precision required in the inversion of the
       *  	 covariance matrix
       *
       *  
       */
      CovarianceMatrix (const std::string filename, const int cov_col=2, const int skipped_lines=0, const double nmeasures=-1, const double prec=1.e-10) 
      {read(filename, cov_col, skipped_lines, prec, nmeasures); }


      /**
       *  @brief default destructor
       */
      virtual ~CovarianceMatrix () = default;

      ///@}

      /**
       *  @name Member functions to get the private/protected members
       */
      ///@{

      /**
       *  @brief get the value of the covariance matrix at index i,j
       *
       *  @param i index
       *
       *  @param j index
       *
       *  @return the value of the covariance matrix
       */
      double operator() (const int i, const int j) const {return m_matrix(i, j);}

      /**
       *  @brief get the covariance matrix
       *  @return the covariance matrix
       */
      std::vector<std::vector<double>> operator() () const;
      
      /**
       *  @brief get the value of the correlation matrix at index i,j
       *
       *  @param i index
       *
       *  @param j index
       *
       *  @return the value of the correlation
       * \f$ Corr_{i,j} = \frac{Cov_{i,j}}{\sqrt{Cov_{i,i} \cdot Cov{j,j}}} \f$
       */
      double correlation (const int i, const int j) const {return m_correlation(i, j);}

      /**
       *  @brief get the value of the correlation matrix at index i,j
       *
       *  @return the value of the correlation
       * \f$ Corr_{i,j} = \frac{Cov_{i,j}}{\sqrt{Cov_{i,i} \cdot Cov{j,j}}} \f$
       */
      std::vector<std::vector<double>> correlation () const;
      
      /**
       *  @brief get the value of the precision matrix at index i,j
       *
       *  @param i index
       *
       *  @param j index
       *
       *  @return the value of the precision matrix at position i,j
       */
      double precision (const int i, const int j) const {return m_precision(i, j);} 

      /**
       *  @brief get the precision matrix
       *
       *  @return the std::vector containing the precision matrix
       */
      std::vector<std::vector<double>> precision () const; 

      /**
       *  @brief get the value of the precision matrix at index i,j
       *  times the Hartlap factor
       *
       *  This function returns the precision matrix times the Hartlap
       *  factor at index i,j:
       *
       *  \f[
       *  \hat{\Psi}=\left(1-\frac{N_{\mathrm{b}}+1}{N_{\mathrm{s}}-1}\right) \hat{\mathrm{c}}^{-1}
       *  \f]
       *
       *  where \f$N_b\f$ are the number of bins (the order of the matrix) and
       *  \f$N_s\f$ the number of data used to measure the covariance
       *
       *  @param i index
       *
       *  @param j index
       *
       *  @return the value of the precision matrix at position i,j
       */
      double precision_hartlap (const int i, const int j) const {return m_hartlap_factor*m_precision(i, j);} 

      /**
       *  @brief get the value of the precision matrix at index i,j
       *  times the Hartlap factor
       *
       *  This function returns the precision matrix times the Hartlap
       *  factor at index i,j:
       *
       *  \f[
       *  \hat{\Psi}=\left(1-\frac{N_{\mathrm{b}}+1}{N_{\mathrm{s}}-1}\right) \hat{\mathrm{c}}^{-1}
       *  \f]
       *
       *  where \f$N_b\f$ are the number of bins (the order of the matrix) and
       *  \f$N_s\f$ the number of data used to measure the covariance
       *
       *  @return the std::vector containing the precision matrix
       */
      std::vector<std::vector<double>> precision_hartlap () const; 

      /**
       *  @brief get value of the standard deviation at index i
       *
       *  @param i index
       *
       *  @return the value of the standard deviation vector at position i
       */
       double standard_deviation (const int i) const { return m_std(i); }

      /**
       *  @brief get the standard deviation
       *
       *  @return the standard deviation
       */
      std::vector<double> standard_deviation () const;

      /**
       *  @brief get value of the variance at index i
       *
       *  @param i index
       *
       *  @return the value of the variance vector at position i
       */
      double variance (const int i) const { return m_variance(i); }

      /**
       *  @brief get the variance
       *
       *  @return the variance
       */
      std::vector<double> variance () const;

      /**
       *   @brief return the covariance matrix order
       *
       *   @return the covariance matrix order 
       */
      size_t order () const { return m_order; }

      ///@}
      

      /**
       *  @name Member functions to set the private/protected members
       */
      ///@{
      
      /**
       *  @brief set the covariance matrix by passing a 
       *  std::vector<std::vector<double>> object;
       *
       *  @param covariance std::vector<std::vector<double>> containing the covariance matrix
       *
       *  @param prec the precision required in the inversion of the
       *  	 covariance matrix
       *
       *  @param nmeasures number of measures used to compute the covariance
       */
      void set_from_matrix (const std::vector<std::vector<double>> covariance, const double nmeasures=-1, const double prec=1.e-10) 
      { m_set(cbl::flatten(covariance), nmeasures, prec); }

      /**
       *  @brief set the covariance by passing the diagonal square root
       *
       *  @param standard_deviation std::vector containing the standard deviation
       *
       *  @param nmeasures number of measures used to compute the covariance
       */
      void set_from_standard_deviation (const std::vector<double> standard_deviation, const double nmeasures=-1);

      /**
       * @brief measure the covariance from a collection of dataset
       *
       * This use the standard estimator of the covariance
       *
       * \f[
       * \hat{C}_{i j}=\frac{f}{N_{\mathrm{s}}-1} \sum_{k=1}^{N_{\mathrm{s}}}
       * \left(D_{i}^{k}-\bar{D}_{i}\right)\left(D_{j}^{k}-\bar{D}_{j}\right)
       * \f]
       *
       * where \f$D_{i}\f$ is the i-th measure \f$ \bar{D}_{i}\f$ the mean.
       * \f$f\f$ is a normalization factor the user can provide in input.
       *
       * @param dataset vector of pointers of object of type Data
       *
       * @param normalization the normalization factor
       *
       * @param prec the precision required in the inversion of the
       * covariance matrix
       */
      void measure (const std::vector<std::shared_ptr<Data>> dataset, const double normalization=1, const double prec=1.e-10);

      /**
       * @brief measure the covariance from a collection of dataset
       *
       * This use the standard estimator of the covariance
       *
       * \f[
       * \hat{C}_{i j}=\frac{f}{N_{\mathrm{s}}-1} \sum_{k=1}^{N_{\mathrm{s}}}
       * \left(D_{i}^{k}-\bar{D}_{i}\right)\left(D_{j}^{k}-\bar{D}_{j}\right)
       * \f]
       *
       * where \f$D_{i}\f$ is the i-th measure \f$ \bar{D}_{i}\f$ the mean.
       * \f$f\f$ is a normalization factor the user can provide in input.
       *
       * This function takes a vector of measurement. Each measurement is a vector
       * of dataset. This should be used when multiple measurements are made
       *
       * @param dataset vector of pointers of object of type Data
       *
       * @param normalization the normalization factor
       *
       * @param prec the precision required in the inversion of the
       * covariance matrix
       */
      void measure (const std::vector<std::vector<std::shared_ptr<Data>>> dataset, const double normalization=1, const double prec=1.e-10);
      
      ///@}

      /**
       *  @name Member functions used for Input/Output
       */
      ///@{

      /**
       *  @brief set the covariance matrix, reading from
       *  an input file
       *
       *  @param filename file containing the covariance matrix in the
       *  format: column 0 \f$ \rightarrow \f$ x<SUB>i</SUB>, column 1
       *  \f$ \rightarrow \f$ x<SUB>j</SUB>, column cov_col &rarr;
       *  cov(x<SUB>i</SUB>,x<SUB>j</SUB>)
       *
       *  @param cov_col covariance matrix column, starting from 0
       *
       *  @param skipped_lines comment lines to be skipped
       *
       *  @param nmeasures number of measures used to compute the covariance
       *
       *  @param prec the precision required in the inversion of the
       *  	 covariance matrix
       */
      void read (const std::string filename, const int cov_col=2, const int skipped_lines=0, const double nmeasures=-1, const double prec=1.e-10);
      
      /**
       *  @brief write the covariance matrix
       *
       *  @param dir output directory
       *
       *  @param file output file
       *
       *  the first line of the output file
       *
       *  @param precision the float precision
       *
       *  @param rank cpu index (for MPI usage)
       */
      void write (const std::string dir, const std::string file, const int precision=4, const int rank=0) const;

      ///@}
      
      /**
       *  @name Member functions for covariance matrix cut
       */

      ///@{

      /**
       * @brief cut the data, for Data1D
       *
       * @param [in] mask std::vector containing values to be masked
       *
       * @return pointer to an object of type Data1D
       */
      CovarianceMatrix cut (const std::vector<bool> mask) const;

      ///@}

      /**
       *  @name Member functions to add two (or more) covariance matrices
       */

      ///@{
      
      /**
       * @brief overloading of the += operator, to sum two catalogues
       *
       * This function add two covariance matrices of order N and M.
       * It creates a block matrix of order N+M
       *
       * @param covariance object of class CovarianceMatrix
       *
       * @return object of class CovarianceMatrix
       */
      CovarianceMatrix operator += (const CovarianceMatrix covariance) const;

      /**
       * @brief overloading of the += operator, to sum two covariance
       * matrices
       * 
       * This function add two covariance matrices of order N and M.
       * It creates a block matrix of order N+M
       *
       * @param covariance object of class CovarianceMatrix
       *
       * @return object of class CovarianceMatrix
       */
      CovarianceMatrix operator += (const std::shared_ptr<CovarianceMatrix> covariance) const;

      /**
       * @brief overloading of the += operator, to sum several covariance
       * matrices, passed in a vector
       *
       * This function add the current covariance matrix of order N
       * with several covariance matrices of order M_i.
       * It creates a block matrix of order \f$ N+\sum_i^n M_i\f$.
       *
       * @param covariance object of class CovarianceMatrix
       *
       * @return object of class CovarianceMatrix
       */
      CovarianceMatrix operator += (const std::vector<CovarianceMatrix> covariance) const;

      /**
       * @brief overloading of the += operator, to sum several covariance
       * matrices, passed in a vector
       *
       * This function add the current covariance matrix of order N
       * with several covariance matrices of order M_i.
       * It creates a block matrix of order \f$ N+\sum_i^n M_i\f$.
       *
       * @param covariance object of class CovarianceMatrix
       *
       * @return object of class CovarianceMatrix
       */
      CovarianceMatrix operator += (const std::vector<std::shared_ptr<CovarianceMatrix>> covariance) const;

      ///@}
      
      
    };


  }
}

#endif
