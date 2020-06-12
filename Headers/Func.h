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
 *  @file Headers/Func.h
 *
 *  @brief Useful generic functions
 *
 *  This file contains the prototypes of a large set of useful
 *  functions of wide used
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __FUNC__
#define __FUNC__

#include "RandomNumbers.h"


// ============================================================================================


namespace cbl {

  
  /**
   *  @name Functions of generic use  
   */
  ///@{
  
  /**
   *  @brief 1D interpolation
   *
   *  @param [in] _xx the point where the input function will be
   *  interpolated or extrapolated
   *
   *  @param [in] xx std::vector containing the binned values of x
   *
   *  @param [in] yy std::vector containing the binned values of the
   *  function, y(x), to be interpolated or extrapolated
   *
   *  @param [in] type the method used to interpolate or extrapolate:
   *  "Linear" &rarr; linear interpolation; "Poly" &rarr; polynomial
   *  interpolation; "Spline" &rarr; cubic spline interpolation; "Rat"
   *  &rarr; diagonal rational function interpolation; "BaryRat"
   *  &rarr; barycentric rational interpolation
   *
   *  @return the interpolated value of the input function
   *
   *  @warning if _xx is outside the range of the input std::vector
   *  xx, the returned value is the extrapolation
   */
  double interpolated (const double _xx, const std::vector<double> xx, const std::vector<double> yy, const std::string type);
  
  /**
   *  @brief 2D interpolation
   *
   *  @param [in] _x1 the point in the first dimension where the input
   *  function will be interpolated
   *
   *  @param [in] _x2 the point in the second dimension where the
   *  input function will be interpolated
   *
   *  @param [in] x1 std::vector containing the binned values of x in the
   *  first dimension
   *
   *  @param [in] x2 std::vector containing the binned values of x in the
   *  second dimension
   *
   *  @param [in] yy std::vector containing the binned values of the
   *  function, y(x), to be interpolated or extrapolated
   *
   *  @param [in] type the method used to interpolate or extrapolate:
   *  "Linear" &rarr; linear interpolation; "Poly" &rarr; polynomial
   *  interpolation; "Spline" &rarr; cubic spline interpolation; "Rat"
   *  &rarr; diagonal rational function interpolation; "BaryRat"
   *  &rarr; barycentric rational interpolation
   *
   *  @return the interpolated or extrapolated value of the
   *  input function
   *
   *  @warning if _x1 and/or _x2 are outside the range of the input
   *  std::vectors x1 and/or x2, the returned value is the extrapolatation
   *
   */
  double interpolated_2D (const double _x1, const double _x2, const std::vector<double> x1, const std::vector<double> x2, const std::vector<std::vector<double>> yy, const std::string type);
  
  /**
   *  @brief filter W(r/r<SUB>c</SUB>), used e.g. for filtering the
   *  correlation function
   *  @param r the scale at which calculate the filter value
   *  @param rc the characteristic filter scale
   *  @return the filter
   */
  double Filter (const double r, const double rc);

  /**
   *  @brief the order l Legendre polynomial 
   *  @param mu the variable mu
   *  @param l the order l Legendre polynomial
   *  @return the order l Legendre polynomial
   */
  double legendre_polynomial (const double mu, const int l);

  /**
   *  @brief the order l Legendre polynomial integrand
   *  @param mu the variable mu
   *  @param params the parameters for the function
   *  @return the order l Legendre polynomial integrand
   */
  double legendre_polynomial_integral (double mu, void *params);

  /**
   *  @brief the average of the Legendre polynomial
   *  of the l-th order over the mu range
   *  @param ll the order of the Legendre polynomial
   *  @param mu the order of the Legendre polynomial
   *  @param delta_mu the order of the Legendre polynomial
   *  @return the average of the Legendre polynomial
   *  of the l-th order over the mu range
   */
  double Legendre_polynomial_mu_average (const int ll, const double mu, const double delta_mu);

  /**
   *  @brief the average of the Legendre polynomial
   *  of the l-th order over the \f$\mu=\cos(\theta)\f$ range
   *  @param mu_min the lower limit of integration of the Legendre polynomial
   *  @param mu_max the upper limit of integration of the Legendre polynomial
   *  @param ll the order of the Legendre polynomial
   *  @return the average of the Legendre polynomial
   *  of the l-th order over the mu range
   */
  double Legendre_polynomial_mu_average (const double mu_min, const double mu_max, const int ll);

  /**
   *  @brief the average of the Legendre polynomial
   *  of the l-th order over the \f$\theta\f$ range
   *  @param theta_min the lower limit of integration of the Legendre polynomial
   *  @param theta_max the upper limit of integration of the Legendre polynomial
   *  @param ll the order of the Legendre polynomial
   *  @return the average of the Legendre polynomial
   *  of the l-th order over the mu range
   */
  double Legendre_polynomial_theta_average (const double theta_min, const double theta_max, const int ll);

  /**
   *  @brief the average of the Legendre polynomial of the l-th order
   *  over the \f$r_{12}, r_{13}, r_{23} \f$
   *
   *  @param r12_min the lower limit of integration for \f$r_{12}\f$
   *  @param r12_max the upper limit of integration for \f$r_{12}\f$
   *  @param r13_min the lower limit of integration for \f$r_{13}\f$
   *  @param r13_max the upper limit of integration for \f$r_{13}\f$
   *  @param r23_min the lower limit of integration for \f$r_{23}\f$
   *  @param r23_max the upper limit of integration for \f$r_{23}\f$
   *  @param ll the order of the Legendre polynomial
   *  @param rel_err the relative error
   *  @param abs_err the absolute error
   *  @param nevals the maximum number of integrals evaluation
   *
   *  @return the average of the Legendre polynomial of the l-th order
   *  over the mu range
   */
  double Legendre_polynomial_triangles_average (const double r12_min, const double r12_max, const double r13_min, const double r13_max, const double r23_min, const double r23_max, const int ll, const double rel_err=1.e-5, const double abs_err=1.e-8, const int nevals=100);

  /**
   *  @brief the average of the Legendre polynomial up to a maximum
   *  order lMax of all triangles with sides r12, r13
   *
   *  @param r12 first triangle side
   *  @param r13 second triangle side
   *  @param deltaR the bin width
   *  @param lMax the maximum Legedre polynomial order
   *  @param rel_err the relative error
   *  @param abs_err the absolute error
   *  @param nevals the maximum number of integrals evaluation
   *
   *  @return the average of the Legendre polynomial of the l-th order
   *  over the mu range
   */
  std::vector<std::vector<double>> Legendre_polynomial_triangles_average (const double r12, const double r13, const double deltaR, const int lMax, const double rel_err=1.e-5, const double abs_err=1.e-8, const int nevals=100);

  /**
   *  @brief the order l, degree m spherical harmonics
   * 
   *  @param l the degree l 
   *
   *  @param m the order m
   *
   *  @param xx the variable x
   *
   *  @param yy the variable y
   *
   *  @param zz the variable z
   *
   *  @return the order l, degree m spherical harmonics
   */
  std::complex<double> spherical_harmonics (const int l, const int m, const double xx, const double yy, const double zz);

  /**
   *  @brief the spherical harmonics up to \f$l_{max}\f$
   * 
   *  @param lmax the maximum degree l
   *
   *  @param xx the variable x
   *
   *  @param yy the variable y
   *
   *  @param zz the variable z
   *
   *  @return the spherical harmonics up to \f$l_{max}\f$
   */
  std::vector<std::vector<std::complex<double>>> spherical_harmonics (const int lmax, const double xx, const double yy, const double zz);

  /**
   *  @brief the spherical harmonics up to \f$l_{max}\f$
   * 
   *  @param lmax the maximum degree l 
   *
   *  @param xx the variable x
   *
   *  @param yy the variable y
   *
   *  @param zz the variable z
   *
   *  @return  the spherical harmonics up to \f$l_{max}\f$
   */
  std::vector<std::complex<double>> spherical_harmonics_array (const int lmax, const double xx, const double yy, const double zz);

  /**
   *  @brief the l=0 spherical Bessel function 
   *  @param xx the variable x
   *  @return the l=0 spherical Bessel function
   */
  double j0 (const double xx);

  /**
   *  @brief the l=2 spherical Bessel function 
   *  @param xx the variable x
   *  @return the l=2 spherical Bessel function
   */
  double j2 (const double xx);

  /**
   *  @brief the l=4 spherical Bessel function 
   *  @param xx the variable x
   *  @return the l=4 spherical Bessel function
   */
  double j4 (const double xx);

  /**
   *  @brief the order l spherical Bessel function 
   *  @param xx the variable x
   *  @param order the order l of spherical Bessel function
   *  @return the order l spherical Bessel function
   */
  double jl (const double xx, const int order);

  /**
   *  @brief the distance average l=0 spherical Bessel function this
   *  function is used to obtain the analytic twop monopole covariance
   *  @param kk the variable k
   *  @param r_down the lower limit of the twopcf bin
   *  @param r_up the upper limit of the twopcf bin
   *  @return the distance average l=0 spherical Bessel function
   */
  double j0_distance_average (const double kk, const double r_down, const double r_up);

  /**
   *  @brief the distance average l=2 spherical Bessel function this
   *  function is used to obtain the analytic twop quadrupole
   *  covariance
   *  @param kk the variable k
   *  @param r_down the lower limit of the twopcf bin
   *  @param r_up the upper limit of the twopcf bin
   *  @return the distance average l=2 spherical Bessel function
   */
  double j2_distance_average (const double kk, const double r_down, const double r_up);
     
  /**
   *  @brief the generic integrand to obtain the distance average 
   *   spherical Bessel function of order l
   *  @param rr the variable r
   *  @param params the parameters for the function
   *  @return the distance average l=2 spherical Bessel function
   */
  double jl_spherical_integrand (double rr, void *params);
  
  /**
   *  @brief the distance average for the order l-th spherical Bessel function 
   *  @param kk the variable k
   *  @param order the shperical Bessel function order
   *  @param r_down the lower limit of the twopcf bin
   *  @param r_up the upper limit of the twopcf bin
   *  @return the distance average l spherical Bessel function
   */
  double jl_distance_average (const double kk, const int order, const double r_down, const double r_up);

  /**
   *  @brief integral, computed with the trapezoid rule, using ordered
   *  data
   *
   *  @param xx the point in which function is defined
   *  @param yy values of the function 
   *  @return the definite integral of the function
   */
  double trapezoid_integration (const std::vector<double> xx, const std::vector<double> yy);

  /**
   *  @brief compute the Wigner 3-j symbol
   *  
   * compute the Wigner 3-j symbol of type
   * \f[
   *	\left(\begin{array}{ccc}{l}} & {l^'} & {l_2}} 
   *	\\ {m_{1}} & {m_{2}} & {m_{3}}\end{array}\right)
   * \f]
   *
   * @param l the first index
   * @param l_prime the second index
   * @param l2 the third index
   * @return the Wigner 3-j symbol
   */
  double coupling_3j (const int l, const int l_prime, const int l2);

  /// @cond glob
  void gauleg (const double, const double, double *, double *, const int);
  /// @endcond


  // ============================================================================================
  
  /**
   *  @brief read a vector from file
   *
   *  @param [in] file_vector input file where the vector is
   *  stored
   *
   *  @param [out] xx the variable
   *
   *  @param [out] vector the vector 
   *
   *  @param [in] col the columns to be read
   *
   *  @return none
   *
   *  @author Alfonso Veropalumbo
   *  @author alfonso.veropalumbo@unibo.it
   */
  void read_vector (const std::string file_vector, std::vector<double> &xx, std::vector<double> &vector, const std::vector<int> col={});
  
  /**
   *  @brief read a matrix from file
   *
   *  @param [in] file_matrix input file where the matrix is
   *  stored
   *
   *  @param [out] xx the variable in row
   *
   *  @param [out] yy the variable in column
   *
   *  @param [out] matrix the matrix 
   *
   *  @param [in] col the columns to be read
   *
   *  @return none
   *
   *  @author Alfonso Veropalumbo
   *  @author alfonso.veropalumbo@unibo.it
   */
  void read_matrix (const std::string file_matrix, std::vector<double> &xx, std::vector<double> &yy, std::vector<std::vector<double>> &matrix, const std::vector<int> col={});

  /**
   *  @brief read a data from a file ASCII
   *
   *  @param file_name the name of the file to read
   *
   *  @param path_name the path where the file is stored
   *
   *  @param column_data vector containing the indices of the columns
   *  to read, starting the counting from 1
   *
   *  @param skip_nlines the number of lines to skip 
   *
   *  @return a vector of vectors containing the columns (first index)
   *  and the lines (second index) read from the file
   *
   *  @author Sofia Contarini
   *  @author sofia.contarini3@unibo.it
   */
  std::vector<std::vector<double>> read_file (const std::string file_name, const std::string path_name, const std::vector<int> column_data, const int skip_nlines=0);

  /**
   *  @brief compute the determinant of a matrix
   *  @param mat the matrix
   *  @return the determinant
   */
  double determinant_matrix (const std::vector<std::vector<double>> mat); 

  /**
   *  @brief function to invert a matrix
   *
   *  this function implements the inversion of a given matrix using
   *  the GSL
   *
   *  if the input matrix is a covariance estimated with a finite
   *  number (Nres>0) of resamples (e.g. via jackknife or bootstrap),
   *  the inverted matrix is corrected as follows (Hartlap, Simon and
   *  Schneider 2006):
   *
   *  \f[ \hat{C}^{-1}=\left(1-\frac{n_b+1}{N_{res}-1}\right)C^{-1} \f]
   *
   *  where \f$n_b\f$ is the number of bins and \f$N_{res}\f$ is the
   *  number of resamplings
   *
   *  @param [in] mat the matrix to be inverted
   *
   *  @param [out] mat_inv the inverted matrix
   *
   *  @param [in] prec the precision required 
   *
   *  @param [in] Nres \f$N_{res}\f$, the number of catalogue
   *  resamplings used to estimate the covariance matrix;
   *  \f$N_{res}=-1\f$ if the covariance matrix has not been estimated
   *  with resampling methods
   *
   *  @return none
   */
  void invert_matrix (const std::vector<std::vector<double>> mat, std::vector<std::vector<double>> &mat_inv, const double prec=1.e-10, const int Nres=-1); 

  /**
   *  @brief function to invert a sub-matrix, extracted from a given
   *  matrix
   *
   *  this function implements the inversion of a sub-matrix,
   *  extracted from a given matrix, using the GSL
   *
   *  if the input matrix is a covariance estimated with a finite
   *  number (Nres>0) of resamples (e.g. via jackknife or bootstrap),
   *  the inverted matrix is corrected as follows (Hartlap, Simon and
   *  Schneider 2006):
   *
   *  \f[ \hat{C}^{-1}=\left(1-\frac{n_b+1}{N_{res}-1}\right)C^{-1} \f]
   *
   *  where \f$n_b\f$ is the number of bins and \f$N_{res}\f$ is the
   *  number of resamplings
   *
   *  @param [in] mat the matrix to be inverted
   *
   *  @param [out] mat_inv the inverted matrix
   *
   *  @param [in] i1 minimum index considered
   *
   *  @param [in] i2 maximum index considered
   *
   *  @param [in] prec the precision required
   *
   *  @param [in] Nres \f$N_{res}\f$, the number of catalogue
   *  resamplings used to estimate the covariance matrix;
   *  \f$N_{res}=-1\f$ if the covariance matrix has not been estimated
   *  with resampling methods
   *
   *  @return none
   */
  void invert_matrix (const std::vector<std::vector<double>> mat, std::vector<std::vector<double>> &mat_inv, const int i1, const int i2, const double prec=1.e-10, const int Nres=-1); 

  /**
   *  @brief compute the covariance matrix from an input dataset
   *
   *  @param [in] mat the data input matrix
   *
   *  @param [out] cov the output covariance matrix
   *
   *  @param [in] JK false &rarr; normalize to 1/(n-1); true &rarr;
   *  normalize to n-1/n (for Jackknife)
   *
   *  @return none
   */
  void covariance_matrix (const std::vector<std::vector<double>> mat, std::vector<std::vector<double>> &cov, const bool JK=false);

  /**
   *  @brief compute the covariance matrix, reading the dataset from files
   *
   *  @param [in] file the std::vector containing the input files
   *   
   *  @param [out] rad the std::vector containing the binned radii
   *
   *  @param [out] mean the std::vector containing the mean values
   *
   *  @param [out] cov the output covariance matrix
   *
   *  @param [in] JK false &rarr; normalize to 1/(n-1); true &rarr;
   *  normalize to n-1/n (for Jackknife)
   *
   *  @return none
   */
  void covariance_matrix (const std::vector<std::string> file, std::vector<double> &rad, std::vector<double> &mean, std::vector<std::vector<double>> &cov, const bool JK=false);

  /**
   *  @brief compute the covariance matrix, reading the dataset from
   *  files, and store it in a file
   *
   *  @param [in] file the std::vector containing the input files
   *
   *  @param [out] covariance_matrix_file the output covariance matrix
   *  file
   *
   *  @param [in] JK false &rarr; normalize to 1/(n-1); true &rarr;
   *  normalize to n-1/n (for Jackknife)
   *
   *  @return none
   */
  void covariance_matrix (const std::vector<std::string> file, const std::string covariance_matrix_file, const bool JK=0);

  /**
   *  @brief read and invert the covariance matrix
   *
   *  @param [in] filecov input file where the covariance matrix is
   *  stored
   *
   *  @param [out] cov the covariance matrix 
   *  
   *  @param [out] cov_inv the inverse of the covariance matrix
   * 
   *  @param [in] i1 mininum index
   *
   *  @param [in] i2 maximum index
   *
   *  @param [in] prec the precision required 
   *
   *  @param [in] Nres \f$N_{res}\f$, the number of catalogue
   *  resamplings used to estimate the covariance matrix;
   *  \f$N_{res}=-1\f$ if the covariance matrix has not been estimated
   *  with resampling methods
   *
   *  @return none
   *
   *  @author Alfonso Veropalumbo
   *  @author alfonso.veropalumbo@unibo.it
   */
  void read_invert_covariance (const std::string filecov, std::vector<std::vector<double>> &cov, std::vector<std::vector<double>> &cov_inv, const size_t i1, const size_t i2, const double prec=1.e-10, const int Nres=-1);

  /**
   *  @brief return a number sampled from a given distribution
   *
   *  @param xx std::vector containing the x variables
   *
   *  @param fx std::vector containing the f(x) variables
   *
   *  @param xmin minimum value of the variable
   *
   *  @param xmax maximum value of the variable
   *
   *  @param seed random seed
   *
   *  @return a std::vector containing the numbers extracted from the
   *  given distribution
   *
   *  @author Alfonso Veropalumbo
   *  @author alfonso.veropalumbo@unibo.it
   */
  double number_from_distribution (const std::vector<double> xx, const std::vector<double> fx, const double xmin, const double xmax, const int seed);

  /**
   *  @brief return a std::vector of numbers sampled from a given
   *  distribution
   *
   *  @param nRan number of elements to be extracted 
   *
   *  @param xx std::vector containing the x variables
   *
   *  @param fx std::vector containing the f(x) variables
   *
   *  @param xmin minimum value of the variable
   *
   *  @param xmax maximum value of the variable
   *
   *  @param seed random seed
   *
   *  @return a std::vector containing the numbers extracted from the
   *  given distribution
   *
   *  @author Alfonso Veropalumbo
   *  @author alfonso.veropalumbo@unibo.it
   */
  std::vector<double> vector_from_distribution (const int nRan, const std::vector<double> xx, const std::vector<double> fx, const double xmin, const double xmax, const int seed);

  /**
   *  @brief return the std::vector indexes corresponding to a given
   *  interval
   *
   *  @param xx input std::vector
   *
   *  @param x_min minimum x value
   *
   *  @param x_max maximum x value
   *
   *  @return the 2D std::vector containing the indexes
   *
   *  @author Alfonso Veropalumbo
   *  @author alfonso.veropalumbo@unibo.it
   */
  std::vector<size_t> minimum_maximum_indexes (const std::vector<double> xx, const double x_min, const double x_max);

  /**
   * @brief get the cosine of the angle between two sides of a
   * triangle
   *
   * \f$ mu = \frac{r_1^2+r_2^2-r_3^3}\frac{2*r_1*r_2}\f$
   *
   * @param r1 the first side of the triangle
   * @param r2 the second side of the triangle
   * @param r3 the third side of the triangle
   *
   * @return the cosine of the angle between two sides of a triangle
   */
  double get_mu (const double r1, const double r2, const double r3);

  /**
   * @brief the unnormalized window function
   *
   * @param x the function variable
   * @param min window function minimum
   * @param max windof function maximum
   *
   * @return  unnormalized window function
   */
  double window_function (const double x, const double min=-1, const double max=1);

  /**
   * @brief get the binomial coefficient
   *
   * @param n first integer
   * @param m second integer
   *
   * @return the binomial coefficient
   */
  double binomial_coefficient(const int n, const int m);

  /**
   * @brief get the Clebsh-Gordan coefficient in the notation
   * \f$ \sum_{l}\left\langle l_1 l_2 m_1 m_2  | l_3 m_3 \right\rangle \f$
   *
   * @param j1 index 
   * @param j2 index
   * @param j3 index
   * @param m1 index
   * @param m2 index
   * @param m3 index
   *
   * @return the Clebsh-Gordan coefficient
   */ 
  double clebsh_gordan(const int j1, const int j2, const int m1, const int m2, const int j3, const int m3);

  /**
   * @brief Wigner \f$3-j\f$ symbol
   *
   * @param j1 index 
   * @param j2 index
   * @param j3 index
   * @param m1 index
   * @param m2 index
   * @param m3 index
   *
   * @return the Clebsh-Gordan coefficient
   */ 
  double wigner_3j(const int j1, const int j2, const int j3, const int m1, const int m2, const int m3);

  /**
   * @brief Wigner \f$6-j\f$ symbol
   *
   * @param j1 index 
   * @param j2 index
   * @param j3 index
   * @param j4 index
   * @param j5 index
   * @param j6 index
   *
   * @return the Clebsh-Gordan coefficient
   */ 
  double wigner_6j(const int j1, const int j2, const int j3, const int j4, const int j5, const int j6);

  /**
   * @brief compute the integral of three spherical bessel function, 
   * from Mehrem, 2011
   *
   * \f[ \int_{0}^{\infty} k^{2} j_{L_{1}}\left(k r_{1}\right)
   *  j_{L_{2}}\left(k r_{2}\right) j_{L_{3}}\left(k r_{3}\right) d k
   *  = \frac{\pi \beta(\Delta)}{8 \pi^2 r_{1} r_{2} r_{3} left\langle
   *  L_{1} L_{2} 00 | L_{3} 0\right\rangle}(i)^{L_{1}+L_{2}+L_{3}}
   *  \left(2 L_{3}+1\right)\left(\frac{r_{1}}{r_{3}}\right)^{L_{3}}
   *  \sum_{L=0}^{L_{3}}\left(\begin{array}{c}{2 L_{3}} \\ {2
   *  L}\end{array}\right)^{1 / 2}\left(\frac{r_{2}}{r_{1}}\right)^{L}
   *  \times \sum_{l}\left\langle L_{1}\left(L_{3}-L\right) 00 |
   *  0\right\rangle\left\langle L_{2} L 00 | l 0\right\rangle\left\
   *  {\begin{array}{lll}{L_{1}} & {L_{2}} & {L_{3}} \\ {L} &
   *  {L_{3}-L} & {l}\end{array}\right\} P_{l}(\Delta) \f]
   *
   * where \f$\sum_{l}\left\langle l_1 l_2 m_1 m_2 | l_3 m_3
   * \right\rangle\f$ is the Clebsh-Gordan coefficient, computed by
   * cbl::clebsh_gordan, and \f${\begin{array}{lll}{L_{1}} & {L_{2}} &
   * {L_{3}} \\ {L} & {L_{3}-L} & {l}\end{array}\right\}\f$ is the
   * \f$6-j\f$ Wigner symbol.
   *
   * @param r1
   * @param r2
   * @param r3
   * @param L1 the order of the first spherical bessel function
   * @param L2 the order of the second spherical bessel function
   * @param L3 the order of the third spherical bessel function
   *
   * @return the integral of three spherical bessel function
   */
  double three_spherical_bessel_integral (const double r1, const double r2, const double r3, const int L1, const int L2, const int L3);

  /**
   *  @brief generate a covariant sample of n points using a
   *  covariance matrix
   *
   *  @param mean the mean values for the sample
   *
   *  @param covariance the covariance matrix of the sample
   *
   *  @param idum seed for random number generator
   *
   *  @return std::vector containing a correlated sample of given mean and
   *  covariance
   *
   *  @author Alfonso Veropalumbo
   *  @author alfonso.veropalumbo@unibo.it
   */
  std::vector<double> generate_correlated_data (const std::vector<double> mean, const std::vector<std::vector<double>> covariance, const int idum =213123);

  /**
   *  @brief generate a covariant sample of n points using a
   *  covariance matrix
   *
   *  @param nExtractions the number of correlated samples to extract
   *
   *  @param mean the mean values for the sample
   *
   *  @param covariance the covariance matrix of the sample
   *
   *  @param idum seed for random number generator
   *
   *  @return std::vector containing a correlated samples of given
   *  mean and covariance
   *
   *  @author Alfonso Veropalumbo
   *  @author alfonso.veropalumbo@unibo.it
   */
  std::vector<std::vector<double>> generate_correlated_data (const int nExtractions, const std::vector<double> mean, const std::vector<std::vector<double>> covariance, const int idum=12312);

  /**
   *  @brief reads a vector from a binary file
   *
   *  @param fin input stream
   *
   *  @param vec the vector container where data will be stored
   *
   *  @param NN the number of elements to be read 
   *
   *  @return none
   *
   *  @author Tommaso Ronconi
   *  @author tronconi@sissa.it
   */
  template <typename T> 
    void vectorReadFromBinary (std::ifstream &fin, std::vector< T > &vec, size_t NN)
    {
      for (size_t ii = 0; ii<NN; ii++) {
	T var;
	fin.read((char *)&var, sizeof(T));
	vec.push_back(var);
      }
    }

  ///@}


  // ============================================================================================


  /**
   *  @name Functions for statistical analyses
   */
  ///@{
  
  /**
   *  @brief the average of a std::vector
   *
   *  for the derivation of the formulae used here for numerically
   *  stable calculation see Chan et al. 1979, Finch 2009 and
   *  reference therein
   *
   *  @param vect the input std::vector
   *  @return the average of vect
   */
  double Average (const std::vector<double> vect);
  
  /**
   *  @brief the weighted average of a std::vector
   *
   *  for the derivation of the formulae used here for numerically
   *  stable calculation see Chan et al. 1979, Finch 2009 and
   *  reference therein
   *
   *  @param vect the input std::vector
   *  @param weight the weight
   *  @return the weighted average of vect
   */
  double Average (const std::vector<double> vect, const std::vector<double> weight);
  
  /**
   *  @brief the standard deviation of a std::vector
   *
   *  for the derivation of the formulae used here for numerically
   *  stable calculation see Chan et al. 1979, Finch 2009 and
   *  reference therein
   *
   *  @param vect the input std::vector
   *  @return &sigma;
   */
  double Sigma (const std::vector<double> vect);
  
  /**
   *  @brief the weighted standard deviation of a std::vector
   *
   *  for the derivation of the formulae used here for numerically
   *  stable calculation see Chan et al. 1979, Finch 2009 and
   *  reference therein
   *
   *  @param vect the input std::vector
   *  @param weight the weight
   *  @return &sigma;
   */
  double Sigma (const std::vector<double> vect, const std::vector<double> weight);
  
  /**
   *  @brief the first, second and third quartiles of a std::vector
   *  @param Vect the input std::vector
   *  @return a std::vector containing the first, second and third
   *  quartiles
   */
  std::vector<double> Quartile (const std::vector<double> Vect);

  /**
   *  @brief compute the moments of a set of data
   *  @param [in] data the std::vector containing the input data
   *  @param [out] ave the mean
   *  @param [out] adev the average deviation
   *  @param [out] sdev the standard deviation
   *  @param [out] var the variance
   *  @param [out] skew the skewness
   *  @param [out] curt the kurtosis
   *  @return none
   */
  void Moment (const std::vector<double> data, double &ave, double &adev, double &sdev, double &var, double &skew, double &curt);

  ///@}


  // ============================================================================================
  
  
  /**
   *  @name Special functions
   */
  ///@{

 
  /**
   *  @brief the quadratic function 
   *  @param xx the variable x
   *  @param pp a void pointer 
   *  @param par a std::vector containing the coefficients
   *  @return the quadratic function: par[0]*x<SUP>2</SUP>+par[1]*x+par[2]
   *  @warning pp is not used, it is necessary only for GSL operations
   */
  template <typename T> 
    T Pol2 (T xx, std::shared_ptr<void> pp, std::vector<double> par)
    {
      (void)pp;
      return par[0]*pow(xx,2)+par[1]*xx+par[2];
    }

  /**
   *  @brief the cubic function 
   *  @param xx the variable x
   *  @param pp a void pointer 
   *  @param par a std::vector containing the coefficients
   *  @return the cubic function:
   *  par[0]*x<SUP>3</SUP>+par[1]*x<SUP>2</SUP>+par[2]*x+par[3]
   *  @warning pp is not used, it is necessary only for GSL operations
   */
  template <typename T> 
    T Pol3 (T xx, void *pp, std::vector<double> par) 
    {
      return par[0]*pow(xx,3)+par[1]*pow(xx,2)+par[2]*xx+par[3];
    }

  /**
   *  @brief linear function
   *  @param xx the coordinate x
   *  @return a std::vector containing [1, x]
   */
  inline std::vector<double> linearfit (const double xx) 
  {
    std::vector<double> vect(2);
    vect[0] = 1.;
    for (int i=1; i<2; i++) vect[i] = xx*vect[i-1];
    return vect;
  }

  /**
   *  @brief quadratic function
   *  @param xx the coordinate x
   *  @return a std::vector containing [1, x, x<SUP>2</SUP>]
   */
  inline std::vector<double> quadratic (const double xx) 
  {
    std::vector<double> vect(3);
    vect[0] = 1.;
    for (int i=1; i<3; i++) vect[i] = xx*vect[i-1];
    return vect;
  }

  /**
   *  @brief cubic function
   *  @param xx the coordinate x
   *  @return a std::vector of containing [1, x, x<SUP>2</SUP>,
   *  x<SUP>3</SUP>]
   */
  inline std::vector<double> cubicfit (const double xx) 
  {
    std::vector<double> vect(4);
    vect[0] = 1.;
    for (int i=1; i<4; i++) vect[i] = xx*vect[i-1];
    return vect;
  }

  /**
   *  @brief the Identity function  
   *  @param xx the variable x
   *  @param pp a void pointer 
   *  @param par a std::vector 
   *  @return 1.
   *
   *  @warning pp and par are not used, but they are necessary in the
   *  function template
   */
  template <typename T> 
    T identity (T xx, std::shared_ptr<void> pp, std::vector<double> par)
    {
      (void)xx; (void)pp; (void)par;
      return 1.;
    }

  /**
   *  @brief the rectangular distribution 
   *  @param xx the variable x
   *  @param pp a void pointer 
   *  @param par a std::vector containing the coefficients: par[0]=lower limix,
   *  par[1]=upper limit;
   *  @return the probability of x
   *
   *  @warning pp is not used, but they are necessary in the
   *  function template
   */
  template <typename T> 
    T rectangular (T xx, std::shared_ptr<void> pp, std::vector<double> par)
    {
      if (xx>par[0] && par[1]>xx)
	return 1./(par[1]-par[0]);
      else return 0.;
    }


  /**
   *  @brief probability of the closest element to x from a list with weights 
   *  @param xx the variable x
   *  @param pp a void pointer 
   *  @param par a std::vector containing the coefficients: 
   *  @return the weight of closest element from a discrete list to x
   *  @warning par is not used, it is necessary only for GSL operations
   */
  double closest_probability (double xx, std::shared_ptr<void> pp, std::vector<double> par);

  /**
   *  @brief probability of an interpolated distribution 
   *  @param xx the variable x
   *  @param pp a void pointer 
   *  @param par a std::vector containing the coefficients: 
   *  @return the weight of closest element from a discrete list to x
   *  @warning par is not used, it is necessary only for GSL operations
   */
  double distribution_probability (double xx, std::shared_ptr<void> pp, std::vector<double> par);

  /**
   *  @brief the Gaussian function 
   *  @param xx the variable x
   *  @param pp a void pointer 
   *  @param par a std::vector containing the coefficients: par[0]=mean,
   *  par[1]=&sigma;
   *  @return the Gaussian function
   *  @warning pp is not used, it is necessary only for GSL operations
   */
  template <typename T> 
    T gaussian (T xx, std::shared_ptr<void> pp, std::vector<double> par)
    {
      (void)pp;
      T gauss = 1./(par[1]*sqrt(2.*par::pi))*exp(-pow(xx-par[0],2)/(2.*par[1]*par[1]));
      return (par.size()==2) ? gauss : gauss*par[2];
    }

  /**
   *  @brief the poisson distribution 
   *  @param xx the variable x
   *  @param pp a void pointer 
   *  @param par a std::vector containing the coefficients: par[0]=mean,
   *  @return the poisson distribution
   *  @warning pp is not used, it is necessary only for GSL operations
   */
  template <typename T>
    T poisson (T xx, std::shared_ptr<void> pp, std::vector<double> par)
    {
      (void)pp;
      T pois = exp(int(xx) * log(par[0]) - lgamma(int(xx) + 1.0) - par[0]);
      return pois;
    }

  /**
   *  @brief the Maxwellian distribution
   *  @param vel velocity 
   *  @param sigma &sigma;
   *  @return the Maxwellian distribution, P(vel)
   */
  template <typename T> 
    T maxwellian_distr (const T vel, const T sigma) 
    {
      return sqrt(54./par::pi)*pow(vel/sigma,2)*exp(-1.5*pow(vel/sigma,2))/sigma;
    }

  /**
   *  @brief the power-law function
   *  @param xx the variable x
   *  @param x0 the normalization, x<SUB>0</SUB>
   *  @param gamma the slope, &gamma;
   *  @return the power-law function: (x/x<SUB>0</SUB>)<SUP>&gamma;</SUP>
   */
  template <typename T> 
    T powerlaw (const T xx, const T x0, const T gamma)
    {
      return pow(xx/x0,-gamma);
    }

  /**
   *  @brief the double power-law function
   *  @param xx the variable x
   *  @param x0 the normalization, x<SUB>0</SUB>
   *  @param alpha the slope, &alpha;
   *  @param beta the slope, &beta;
   *  @return the double power-law function
   */
  template <typename T> 
    T double_powerlaw (const T xx, const T x0, const T alpha, const T beta) 
    {
      return pow(2.,beta-alpha)/(pow(xx/x0,alpha)*pow(1.+xx/x0,beta-alpha));
    }

  /**
   *  @brief the top-hat window function
   *  @param kR the variable k*R
   *  @return the top-hat window function
   */
  template <typename T> 
    T TopHat_WF (const T kR) 
    {
      return 3.*(sin(kR)-kR*cos(kR))/pow(kR,3);
    }
  
  /**
   *  @brief the derivative of the top-hat 
   *  window function
   *  @param kR the variable k*R
   *  @return the derivative of the top-hat 
   *  window function
   */
  template <typename T> 
    T TopHat_WF_D1 (const T kR) 
    {
      return (3.*(kR*kR-3.)*sin(kR)+9.*kR*cos(kR))*pow(kR, -4);
    }

  /**
   *  @brief the radius of a sphere of a given mass and density
   *  @param Mass the mass
   *  @param Rho the density
   *  @return the radius
   */
  template <typename T> 
    T Radius (const T Mass, const T Rho) 
    {
      return pow(3.*Mass/(4.*par::pi*Rho),1./3.); 
    }

  /**
   *  @brief the volume of a sphere of a given radius
   *  @param RR the radius
   *  @return the volume
   */
  template <typename T> 
    T volume_sphere (const T RR) 
    {
      return 4./3.*par::pi*pow(RR, 3);
    }

  /**
   *  @brief the mass of a sphere of a given radius and density
   *  @param RR the radius
   *  @param Rho the density
   *  @return the mass
   */
  template <typename T> 
    T Mass (const T RR, const T Rho) 
    {
      return 4./3.*par::pi*Rho*pow(RR,3);
    }
  
  /**
   *  @brief the radial velocity
   *  @param vx the velocity along the x direction
   *  @param vy the velocity along the y direction
   *  @param vz the velocity along the z direction
   *  @param ra the Right Ascension
   *  @param dec the Declination
   *  @return the radial velocity
   */
  template <typename T> 
    T radial_velocity (const T vx, const T vy, const T vz, const T ra, const T dec)
    {
      return vx*cos(dec)*sin(ra)+vy*cos(dec)*cos(ra)+vz*sin(dec);
    }

  /**
   *  @brief the Legendre polynomial P<SUB>2</SUB>
   *  @param x the variable x
   *  @return P<SUB>2</SUB>
   */
  template <typename T> 
    T P_2 (const T x) 
    {
      return (3.*x*x-1.)*0.5;
    }
      
  /**
   *  @brief the Legendre polynomial P<SUB>4</SUB>
   *  @param x the variable x
   *  @return P<SUB>4</SUB>
   */
  template <typename T> 
    T P_4 (const T x) 
    {
      return (35.*x*x*x*x-30.*x*x+3.)*0.125;
    }
  
  /**
   *  @brief the Legendre polynomial P<SUB>6</SUB>
   *  @param x the variable x
   *  @return P<SUB>6</SUB>
   */
  template <typename T> 
    T P_6 (const T x) 
    {
      return (231.*x*x*x*x*x*x-315.*x*x*x*x+105.*x*x-5.)*0.0625;
    }

  ///@}
  

  // ============================================================================================


  /**
   *  @name Generic operations on functions 
   */
  ///@{

  /**
   *  @brief measure the var function (where "var" could be the mass,
   *  the luminosity, the radius, etc.)
   *  @param [in] var std::vector containing the set of data
   *  @param [in] bin the number of bin used
   *  @param [in] V_min the minimum value of the range
   *  @param [in] V_max the maximum value of the range
   *  @param [in] Volume the volume
   *  @param [out] Var std::vector containing the binned values of "var" 
   *  @param [out] Phi std::vector containing the binned values of the var
   *  function
   *  @param [out] err std::vector containing the Poisson errors
   *  @return none
   */
  void measure_var_function (const std::vector<double> var, const int bin, const double V_min, const double V_max, const double Volume, std::vector<double> &Var, std::vector<double> &Phi, std::vector<double> &err);

  /**
   *  @brief derive and store the number distribution of a given
   *  std::vector 
   *  @param [out] xx std::vector containing the binned values of the
   *  variable 
   *  @param [out] fx std::vector containing the binned values of the
   *  distribution
   *  @param [out] err std::vector containing the binned Poisson errors
   *  @param [in] FF std::vector containing the given set of data
   *  @param [in] WW std::vector containing the weights
   *  @param [in] nbin the number of bins
   *  @param [in] linear true &rarr; linear binning; false &rarr; logarithmic
   *  binning
   *  @param [in] file_out the output file where the distribution is
   *  stored
   *  @param [in] fact factor used to normalized the distribution
   *  @param [in] V1 the minimum limit of the distribution
   *  @param [in] V2 the maximum limit of the distribution
   *  @param [in] bin_type "Linear" &rarr; dn/dvar; "Log10" &rarr; dn/dlog(var); "Log" &rarr; dn/dln(var)
   *  @param [in] conv true &rarr; compute the Gaussian convolvolution of
   *  the distribution; false &rarr; do not convolve
   *  @param [in] sigma &sigma; of the Gaussian kernel
   *  @return none
   */
  void distribution (std::vector<double> &xx, std::vector<double> &fx, std::vector<double> &err, const std::vector<double> FF, const std::vector<double> WW, const int nbin, const bool linear=true, const std::string file_out=par::defaultString, const double fact=1., const double V1=par::defaultDouble, const double V2=par::defaultDouble, const std::string bin_type="Linear", const bool conv=false, const double sigma=0.);

  /**
   *  @brief simple Monte Carlo integration of f(x)
   *  @param func the function f(x)
   *  @param x1 minimum limit of the integral
   *  @param x2 maximum limit of the integral
   *  @param seed the seed for random number generation
   *  @return \f$\int_{x1}^{x2} f(x)dx\f$
   */
  double MC_Int (double func(const double), const double x1, const double x2, const int seed=3213); 

  /**
   *  @brief simple Monte Carlo integration of f(x,A)
   *  @param func the function f(x,A)
   *  @param AA parameter of the function, A
   *  @param x1 minimum limit of the integral
   *  @param x2 maximum limit of the integral
   *  @param seed the seed for random number generation
   *  @return \f$\int_{x1}^{x2} f(x)dx\f$
   */
  double MC_Int (double func(const double, const double AA), const double AA, const double x1, double x2, const int seed=3213); 
  
  /**
   *  @brief simple Monte Carlo integration of f(x,A,B,C,D,E)
   *  @param func the function f(x,A,B,C,D,E)
   *  @param AA parameter of the function, A
   *  @param BB parameter of the function, B
   *  @param CC parameter of the function, C
   *  @param DD parameter of the function, D
   *  @param EE parameter of the function, E
   *  @param x1 minimum limit of the integral
   *  @param x2 maximum limit of the integral
   *  @param seed the seed for random number generation
   *  @return \f$\int_{x1}^{x2} f(x,A,B,C,D,E)dx\f$
   */
  double MC_Int (double func(const double, const double AA, const double BB, const double CC, const double DD, const double EE), const double AA, const double BB, const double CC, const double DD, const double EE, const double x1, const double x2, const int seed=3213); 

  
  // ============================================================================================

  
  /**
   *  @brief create a 1D grid given an input function
   *  @param file_grid the file with the input function
   *  @param [in] func the input function
   *  @param [in] par the function parameters
   *  @param [in] bin the number of bins in the grid
   *  @param [in] x_min the minimum limit of the grid
   *  @param [in] x_max the maximum limit of the grid
   *  @param [in] binning the binning type: it can be "lin", "loglin" or
   *  "log"
   *  @param [in,out] xx std::vector containing the grid points
   *  @param [in,out] yy std::vector containing the values of the function at the
   *  grid points
   *  @return none
   */
  void bin_function (const std::string file_grid, double func(double, void*), void *par, const int bin, const double x_min, const double x_max, const std::string binning, std::vector<double> &xx, std::vector<double> &yy);

  /**
   *  @brief create a 2D grid given an input function
   *  @param [in] file_grid the file with the input function
   *  @param [in] func the input function
   *  @param [in] par the function parameters
   *  @param [in] bin the number of bins in the grid
   *  @param [in] x1_min the minimum limit of the grid in one direction
   *  @param [in] x1_max the maximum limit of the grid in one direction
   *  @param [in] x2_min the minimum limit of the grid in one direction
   *  @param [in] x2_max the maximum limit of the grid in one direction
   *  @param [in] binning the binning type: it can be "lin", "loglin" or
   *  "log"
   *  @param [in,out] xx1 std::vector containing the grid points in one direction
   *  @param [in,out] xx2 std::vector containing the grid points in one direction
   *  @param [in,out] yy std::vector containing the values of the function at the
   *  grid points
   *  @return none
   */
  void bin_function_2D (const std::string file_grid, double func(double *, size_t, void *), void * par, const int bin, const double x1_min, const double x1_max, const double x2_min, const double x2_max, const std::string binning, std::vector<double> &xx1, std::vector<double> &xx2, std::vector<std::vector<double>> &yy);

  /// @cond glob
  double func_grid_lin (double, void *);
  double func_grid_loglin (double, void *);
  double func_grid_log (double, void *);
  double func_grid_lin_2D (double *, size_t, void *);
  double func_grid_loglin_2D (double *, size_t, void *);
  double func_grid_log_2D (double *, size_t, void *);
  /// @endcond

  
  /**
   *  @brief convolution of two functions
   *
   *  Get the convolution of the two functions f<SUB>1</SUB>(x) and
   *  f<SUB>2</SUB>(x), and store it in the output std::vector res. The two
   *  functions have to be defined on the same x-axis range, with
   *  equal number of points, &Delta;x =
   *  (x<SUB>max</SUB>-x<SUB>min</SUB>)/n<SUB>x</SUB>.
   *
   *  @param [in] f1 first function, f<SUB>1</SUB>(x)
   *  @param [in] f2 second function, f<SUB>2</SUB>(x)
   *  @param [out] res convolution function
   *  @param [in] deltaX &Delta;x =
   *  (x<SUB>max</SUB>-x<SUB>min</SUB>)/n<SUB>x</SUB>
   *  @return none
   *  
   *  @author Alfonso Veropalumbo
   *  @author alfonso.veropalumbo@unibo.it
   */
  // 
  void convolution (const std::vector<double> f1, const std::vector<double> f2, std::vector<double> &res, const double deltaX);

  ///@}
  

  // ============================================================================================


  /**
   *  @name Functions to calculate distances 
   */
  ///@{

  /**
   *  @brief conversion to degrees
   *  @param angle the input angle 
   *  @param inputUnits the units of the input angle
   *  @return the angle in degrees 
   */
  double degrees (const double angle, const CoordinateUnits inputUnits=CoordinateUnits::_radians_);
  
  /**
   *  @brief conversion to radians
   *  @param angle the input angle
   *  @param inputUnits the units of the input angle
   *  @return the angle in radians
   */
  double radians (const double angle, const CoordinateUnits inputUnits=CoordinateUnits::_degrees_);
  
  /**
   *  @brief conversion to arcseconds
   *  @param angle the input angle 
   *  @param inputUnits the units of the input angle
   *  @return the angle in arcseconds
   */
  double arcseconds (const double angle, const CoordinateUnits inputUnits=CoordinateUnits::_radians_);
  
  /**
   *  @brief conversion to arcminutes
   *  @param angle the input angle
   *  @param inputUnits the units of the input angle
   *  @return the angle in arcminutes
   */
  double arcminutes (const double angle, const CoordinateUnits inputUnits=CoordinateUnits::_radians_);

  /**
   *  @brief conversion to angle units
   *  @param angle the input angle
   *  @param inputUnits the units of the input angle
   *  @param outputUnits the units of the output angle
   *  @return the angle in the converted units
   */
  double converted_angle (const double angle, const CoordinateUnits inputUnits=CoordinateUnits::_radians_, const CoordinateUnits outputUnits=CoordinateUnits::_degrees_);
    
  /**
   *  @brief conversion from Cartesian coordinates to polar
   *  coordinates
   *
   *  @param [in] XX the Cartesian coordinate x
   *  @param [in] YY the Cartesian coordinate y
   *  @param [in] ZZ the Cartesian coordinate z
   *  @param [out] ra the Right Ascension [radians]
   *  @param [out] dec the Declination [radians]
   *  @param [out] dd the comoving distance
   *  @return none
   */
  void polar_coord (const double XX, const double YY, const double ZZ, double &ra, double &dec, double &dd); 

  /**
   *  @brief conversion from polar coordinates to Cartesian
   *  coordinates
   *
   *  @param [in] ra the Right Ascension [radians]
   *  @param [in] dec the Declination [radians]
   *  @param [in] dd the comoving distance 
   *  @param [out] XX the Cartesian coordinate x
   *  @param [out] YY the Cartesian coordinate y
   *  @param [out] ZZ the Cartesian coordinate z
   *  @return none
   */
  void cartesian_coord (const double ra, const double dec, const double dd, double &XX, double &YY, double &ZZ);

  /**
   *  @brief conversion from Cartesian coordinates to polar
   *  coordinates used for a set of objects
   *
   *  @param [in] XX std::vector containing the Cartesian coordinates x
   *  @param [in] YY std::vector containing the Cartesian coordinates y
   *  @param [in] ZZ std::vector containing the Cartesian coordinates z
   *  @param [out] ra std::vector containing the Right Ascension values [radians]
   *  @param [out] dec std::vector containing the Declination values [radians]
   *  @param [out] dd std::vector containing the comoving distances
   *  @return none
   */
  void polar_coord (const std::vector<double> XX, const std::vector<double> YY, const std::vector<double> ZZ, std::vector<double> &ra, std::vector<double> &dec, std::vector<double> &dd); 

  /**
   *  @brief conversion from polar coordinates to Cartesian
   *  coordinates used for a set of objects
   *
   *  @param [in] ra std::vector containing the Right Ascension values [radians]
   *  @param [in] dec std::vector containing the Declination values [radians]
   *  @param [in] dd std::vector containing the comoving distances
   *  @param [out] XX std::vector containing the Cartesian coordinates x
   *  @param [out] YY std::vector containing the Cartesian coordinates y
   *  @param [out] ZZ std::vector containing the Cartesian coordinates z
   *  @return none
   */
  void cartesian_coord (const std::vector<double> ra, const std::vector<double> dec, const std::vector<double> dd, std::vector<double> &XX, std::vector<double> &YY, std::vector<double> &ZZ);

  /**
   *  @brief the Euclidean distance in 3D relative to the origin
   *  (0,0,0), i.e. the Euclidean norm
   *
   *  @param x1 x coordinate of the first object
   *  @param x2 x coordinate of the second object
   *  @param y1 y coordinate of the first object
   *  @param y2 y coordinate of the second object
   *  @param z1 z coordinate of the first object
   *  @param z2 z coordinate of the second object
   *  @return the Euclidean distance
   */
  double Euclidean_distance (const double x1, const double x2, const double y1, const double y2, const double z1, const double z2);

  /**
   *  @brief the perpendicular separation, r<SUB>p</SUB>
   *
   *  @param ra1 the Right Ascension of the first object [radians]
   *  @param ra2 the Right Ascension of the second object [radians]
   *  @param dec1 the Declination of the first object [radians]
   *  @param dec2 the Declination of the second object [radians]
   *  @param d1 the comoving distance of the first object
   *  @param d2 the comoving distance of the second object
   *  @return the perpendicular separation, r<SUB>p</SUB>
   */
  double perpendicular_distance (const double ra1, const double ra2, const double dec1, const double dec2, const double d1, const double d2);
  
  /**
   *  @brief the angular separation in 3D 
   *  @author Alfonso Veropalumbo
   *  @author alfonso.veropalumbo@unibo.it
   *
   *  @param x1 x coordinate of the first object
   *  @param x2 x coordinate of the second object
   *  @param y1 y coordinate of the first object
   *  @param y2 y coordinate of the second object
   *  @param z1 z coordinate of the first object
   *  @param z2 z coordinate of the second object
   *  @return the angular separation [radians]
   */
  double angular_distance (const double x1, const double x2, const double y1, const double y2, const double z1, const double z2);

  /**
   *  @brief the haversine angular separation
   *  @author Alfonso Veropalumbo
   *  @author alfonso.veropalumbo@unibo.it
   *
   *  @param ra1 the Right Ascension of the first object [radians]
   *  @param ra2 the Right Ascension of the second object [radians]
   *  @param dec1 the Declination of the first object [radians]
   *  @param dec2 the Declination of the second object [radians]
   *  @return the haversine angular separation [radians]
   */
  double haversine_distance (const double ra1, const double ra2, const double dec1, const double dec2);
  

  /* ======== Alfonso Veropalumbo ======== */

  /**
   * @brief check if ra coordinate
   * is inside the boundaries
   *
   * @param angle the angle
   *
   * @param minval angle lower limit
   *
   * @param maxval angle upper limit
   *
   * @return none
   */
  void sdss_atbound (double &angle, const double minval, const double maxval);

  /**
   * @brief set the angular coordinates in the
   * SDSS boundaries
   *
   * @param theta the first angular coordinate
   *
   * @param phi the second angular coordinate
   *
   * @return none
   */
  void sdss_atbound2 (double &theta, double &phi);

  /**
   * @brief convert from equatorial
   * coordinates R.A., Dec to sdss coordinates
   * \f$ \lambda, \eta \f$
   *
   * @param ra vector containing R.A. values
   *
   * @param dec vector containing Dec. values
   *
   * @param lambda vector containing the \f$ \lambda \f$ values
   *
   * @param eta vector containing the \f$ \eta \f$ values
   *
   * @return none
   */
  void eq2sdss (const std::vector<double> ra, const std::vector<double> dec, std::vector<double> &lambda, std::vector<double> &eta); 

  /**
   * @brief convert sdss coordinates
   * \f$ \lambda, \eta \f$ to R.A., Dec.
   *
   * @param lambda vector containing the \f$ \lambda \f$ values
   *
   * @param eta vector containing the \f$ \eta \f$ values
   *
   * @param ra vector containing R.A. values
   *
   * @param dec vector containing Dec. values
   *
   * @return none
   */
  void sdss2eq (const std::vector<double> lambda, const std::vector<double> eta, std::vector<double> &ra, std::vector<double> &dec);
  
  /**
   * @brief compute the SDSS stripe given 
   * SDSS coordinates \f$ \lambda, \eta \f$
   *
   * @param eta vector containing the \f$ \eta \f$ values
   *
   * @param lambda vector containing the \f$ \lambda \f$ values
   *
   * @param stripe vector containing the stripe value
   *
   * @param str_u vector containing the list of stripes
   *
   * @return none
   */
  void sdss_stripe (const std::vector<double> eta, const std::vector<double> lambda, std::vector<int> &stripe, std::vector<int> &str_u);


  ///@}


  // ============================================================================================
  
  
  /**
   *  @name Functions to model the correlation function
   */
  ///@{
  
 
  /**
   *  @brief the two-point correlation function computed from the
   *  Fourier transform of the power spectrum
   * 
   *  @param rr the comoving separation, r
   *
   *  @param lgkk std::vector containing the logarithm of the wave std::vectors,
   *  log<SUB>10</SUB>k
   *
   *  @param lgPk std::vector containing the logarithm of the power
   *  spectrum, log<SUB>10</SUB>P(k)
   *
   *  @param k_min the minimum value of the wave std::vector used in the
   *  integral of the Fourier transform
   *
   *  @param k_max the maximum value of the wave std::vector used in the
   *  integral of the Fourier transform
   *
   *  @param aa parameter used to smooth the integrand, given by the
   *  eq. 24 of Anderson et al. 2012
   *
   *  @param prec accuracy of the GSL integration 
   *
   *  @return the two-point correlation function, &xi;(r)
   */
  double xi_from_Pk (const double rr, const std::vector<double> lgkk, const std::vector<double> lgPk, const double k_min=0., const double k_max=100., const double aa=0., const double prec=1.e-2);

  /**
   *  @brief the two-point correlation function computed from the
   *  Fourier transform of the power spectrum read from a file
   * 
   *  @param rr the comoving separation, r
   *
   *  @param file name of the file where the power spectrum is stored
   *
   *  @param c1 the column of the file corresponding to the wave
   *  std::vector, k
   *
   *  @param c2 the column of the file corresponding to the power
   *  spectrum, P(k)
   *
   *  @param k_min the minimum value of the wave std::vector used in the
   *  integral of the Fourier transform
   *
   *  @param k_max the maximum value of the wave std::vector used in the
   *  integral of the Fourier transform
   *
   *  @param aa parameter used to smooth the integrand, given by the
   *  eq. 24 of Anderson et al. 2012
   *
   *  @param prec accuracy of the GSL integration 
   *
   *  @return the two-point correlation function, &xi;(r)
   */
  double xi_from_Pk (const double rr, const std::string file, const int c1=1, const int c2=2, const double k_min=0., const double k_max=100., const double aa=0., const double prec=1.e-2);

  /**
   *  @brief the power spectrum computed from the Fourier transform of
   *  the two-point correlation function
   * 
   *  @param kk the wave std::vector, k
   *
   *  @param lgrr std::vector containing the logarithm of the comoving
   *  separations, log<SUB>10</SUB>r
   *
   *  @param lgxi std::vector containing the logarithm of the two-point
   *  correlation function, log<SUB>10</SUB>&xi;(r)
   *
   *  @param r_min the minimum value of the comoving separation used
   *  in the integral of the Fourier transform
   *
   *  @param r_max the maximum value of the comoving separation used
   *  in the integral of the Fourier transform
   *
   *  @return the power spectrum, P(k)
   */
  double Pk_from_xi (const double kk, const std::vector<double> lgrr, const std::vector<double> lgxi, const double r_min=0.03, const double r_max=100.); 

  /**
   *  @brief the power spectrum computed from the Fourier transform of
   *  the two-point correlation function read from a file
   * 
   *  @param kk the wave std::vector, k
   *
   *  @param file name of the file where the two-point correlation
   *  function is stored
   *
   *  @param c1 the column of the file corresponding to the comoving
   *  separation, r
   *
   *  @param c2 the column of the file corresponding to the two-point
   *  correlation function, &xi;(r)
   *
   *  @param r_min the minimum value of the comoving separation used
   *  in the integral of the Fourier transform
   *
   *  @param r_max the maximum value of the comoving separation used
   *  in the integral of the Fourier transform
   *
   *  @return the power spectrum, P(k)
   */
  double Pk_from_xi (const double kk, const std::string file, const int c1=1, const int c2=2, const double r_min=0.03, const double r_max=100.); 

  namespace glob {
    /// @cond glob
    double func_xi_GSL (double, void *);
    double func_SSM_GSL (double, void *);
    /// @endcond
  }

  /**
   *  @brief the projected two-point correlation function
   * 
   *  this function estimates the projected correlation function by
   *  integrating a given two-point correlation function as follows:
   *
   *  \f[
   *  w_p(r_p)=2\int_{r_p}^{r_{max}}\frac{\xi(r)}{\sqrt{r^2-r_p^2}}r{\rm
   *  d}r \f]
   *
   *  @param rp r<SUB>p</SUB>: comoving separation perpendicular to
   *  the line-of-sight
   *
   *  @param rr std::vector containing the central values of the binned
   *  comoving separations, r
   *
   *  @param xi std::vector containing the central values of the binned
   *  two-point correlation function, &xi;(r)
   *
   *  @param r_max the maximum value of the comoving separation used
   *  in the integral 
   *
   *  @return the projected correlation function, w(r<SUB>p</SUB>)
   */
  double wp (const double rp, const std::vector<double> rr, const std::vector<double> xi, const double r_max=100.); 

  /**
   *  @brief the projected two-point correlation function
   * 
   *  this function estimates the projected correlation function by
   *  integrating a given two-point correlation function as follows:
   *  \f[
   *  w_p(r_p)=\int_{r_p}^{r_{max}}\frac{\xi(r)}{\sqrt{r^2-r_p^2}}rdr
   *  \f]
   *
   *  the two-point correlation function is read from a file
   *
   *  @param rp r<SUB>p</SUB>: comoving separation perpendicular to
   *  the line-of-sight
   *
   *  @param file name of the file where the two-point correlation
   *  function is stored
   *
   *  @param r_max the maximum value of the comoving separation used
   *  in the integral 
   *
   *  @return the projected correlation function,
   *  w<SUB>p</SUB>(r<SUB>p</SUB>)
   */
  double wp (const double rp, const std::string file, const double r_max=100.); 

  /**
   *  @brief the rms mass fluctuation within a given radius
   * 
   *  computed with the spherically averaged correlation function as
   *  follows:
   *
   *  \f[
   *  \sigma_R^2=\int_0^2dy\,y^2\xi(yR)\left(3-\frac{9y}{4}+\frac{3y^3}{16}\right)
   *  \f]
   *
   *  or with the projected correlation function as follows:
   *
   *  \f[
   *  \sigma_R^2=\frac{1}{R^3}\int_0^\infty dr_p\,r_pw_p(r_p)g(r_p/R)
   *  \f]
   *
   *  where g(x) is:
   *
   *  \f[
   *  \left\{
   *  \begin{array}{ll}
   *  \frac{1}{2\pi}\left(3\pi-9x+x^3\right) & \mbox{for}\,x\le2 \\
   *  \frac{1}{2\pi}\left(\frac{-x^4+11x^2-28}{\sqrt{x^2-4}}+x^3-9x+6\sin^{-1}(2/x)\right) & \mbox{for}\,x>2 \\
   *  \end{array}
   *  \right.
   *  \f] 
   *  see e.g. Zehavi et al. 2005, ApJ, 621, 22
   *
   *  @param RR the radius R [Mpc/h]
   *
   *  @param corrType 1 &rarr; the spherically averaged correlation
   *  function is used; 2 &rarr; the projected correlation function is
   *  used
   *
   *  @param rr std::vector containing comoving separations
   *
   *  @param corr std::vector containing the two-point correlation function
   *
   *  @return &sigma;<SUB>R</SUB>: the rms mass fluctuation within a radius R [Mpc/h] 
   */
  double sigmaR (const double RR, const int corrType, const std::vector<double> rr, const std::vector<double> corr);

  /**
   *  @brief the projected correlation function,
   *  w<SUB>p</SUB>(r<SUB>p</SUB>), computed by modelling the
   *  two-point correlation function, &xi;(r), as a power-law
   * 
   *  @param rp r<SUB>p</SUB>: comoving separation perpendicular to
   *  the line-of-sight
   *
   *  @param r0 r<SUB>0</SUB>: the clustering normalization
   *
   *  @param gamma &gamma;: the clustering slope
   *
   *  @return the projected correlation function
   *  w<SUB>p</SUB>(r<SUB>p</SUB>)
   */
  double xi_projected_powerlaw (const double rp, const double r0, const double gamma); 

  /**
   *  @brief the ratio between the redshift-space and real-space
   *  correlation functions
   *
   *  as predicted by the large-scale limit of the Kaiser/Hamilton
   *  model:
   * 
   *  \f[ \frac{\xi(s)}{\xi(r)} = 1 + \frac{2\beta}{3} +
   *  \frac{\beta^2}{5} \f]
   *
   *  @param beta &beta;=f/b, where f is the linear growth rate and b
   *  is the bias
   *
   *  @return &xi;(s)/&xi;(r)
   */
  double xi_ratio (const double beta);                         

  /**
   *  @brief the ratio between the redshift-space and real-space
   *  correlation functions
   *
   *  as predicted by the large-scale limit of the Kaiser/Hamilton
   *  model:
   * 
   *  \f[ \frac{\xi(s)}{\xi(r)} = 1 +
   *  \frac{2}{3}\frac{f\sigma_8}{b\sigma_8} +
   *  \frac{1}{5}\left(\frac{f\sigma_8}{b\sigma_8}\right)^2 \f]
   *
   *  @param f_sigma8 f*&sigma;<SUB>8</SUB>
   *  
   *  @param bias_sigma8 b*&sigma;<SUB>8</SUB>
   *
   *  @return &xi;(s)/&xi;(r)
   */
  double xi_ratio (const double f_sigma8, const double bias_sigma8);                

  /// @cond glob
  /**
   *  @brief the ratio between the redshift-space and real-space
   *  correlation functions
   *
   *  as predicted by the large-scale limit of the Kaiser/Hamilton
   *  model:

   *  \f[ \frac{\xi(s)}{\xi(r)} = 1 + \frac{2\beta}{3} +
   *  \frac{\beta^2}{5} \f]
   *
   *  @param xx coordinate x
   *  @param pp void pointer
   *  @param par std::vector containing one or two parameters
   *  @return &xi;(s)/&xi;(r)
   */
  double xi_ratio (double xx, std::shared_ptr<void> pp, std::vector<double> par); // for the chi2 function
  /// @endcond

  /**
   *  @brief error on the ratio between the redshift-space and
   *  real-space correlation functions
   *
   *  as predicted by the large-scale limit of the Kaiser/Hamilton
   *  model:
   * 
   *  \f[ \delta\left[\frac{\xi(s)}{\xi(r)}\right] = \left(\frac{2}{3}
   *  + \frac{2\beta}{5}\right)\cdot\delta(\beta) \f]
   *
   *  @param beta &beta;=f/b
   *  
   *  @param error_beta error on &beta;
   *
   *  @return &delta;[&xi;(s)/&xi;(r)]
   */
  double error_xi_ratio (const double beta, const double error_beta); 

  /**
   *  @brief the barred correlation function
   *
   *  computed with a simple rectangular integration 
   *
   *  (see e.g. Hamilton 1992)
   *
   *  \f[ \overline{\xi}(r) = \frac{3}{r^3}\int^r_0dr'\xi(r')r'{^2}
   *  \f]
   *
   *  @param RR the comoving separation, r, where the barred
   *  correlation function is computed
   *
   *  @param rr std::vector containing the input comoving separations
   *
   *  @param xi std::vector containing the input two-point correlation
   *  function
   *
   *  @param rAPP comoving scale below which a power-law model for the
   *  two-point correlation function is assumed in the integral
   *
   *  @param r0 the power-law normalization, r<SUB>0</SUB>  
   *
   *  @param gamma the power-law slope, &gamma;
   *
   *  @return \f$ \overline{\xi}(r) \f$
   */
  double barred_xi_direct (const double RR, const std::vector<double> rr, const std::vector<double> xi, const double rAPP=0., const double r0=-1., const double gamma=1.); 

  /**
   *  @brief the double barred correlation function
   *
   *  computed with a simple rectangular integration 
   *
   *  (see e.g. Hamilton 1992)
   *
   *  \f[ \overline{\overline{\xi}}(r) = \frac{5}{r^5}\int^r_0dr'\xi(r')r'{^4}
   *  \f]
   *
   *  @param RR the comoving separation, r, where the barred
   *  correlation function is computed
   *
   *  @param rr std::vector containing the input comoving separations
   *
   *  @param xi std::vector containing the input two-point correlation
   *  function
   *
   *  @param rAPP comoving scale below which a power-law model for the
   *  two-point correlation function is assumed in the integral
   *
   *  @param r0 the power-law normalization, r<SUB>0</SUB>  
   *
   *  @param gamma the power-law slope, &gamma;
   *
   *  @return \f$ \overline{\overline{\xi}}(r) \f$
   */
  double barred_xi__direct (const double RR, const std::vector<double> rr, const std::vector<double> xi, const double rAPP=0., const double r0=-1., const double gamma=1.); 

  /**
   *  @brief the barred correlation function
   *
   *  (see e.g. Hamilton 1992)
   *
   *  \f[ \overline{\xi}(r) = \frac{3}{r^3}\int^r_0dr'\xi(r')r'{^2}
   *  \f]
   *
   *  @param RR the comoving separation, r, where the barred
   *  correlation function is computed
   *
   *  @param rr std::vector containing the input comoving separations
   *
   *  @param xi std::vector containing the input two-point correlation
   *  function
   *
   *  @param rAPP comoving scale below which a power-law model for the
   *  two-point correlation function is assumed in the integral
   *
   *  @param r0 the power-law normalization, r<SUB>0</SUB>  
   *
   *  @param gamma the power-law slope, &gamma;
   *
   *  @return \f$ \overline{\xi}(r) \f$
   */
  double barred_xi_ (const double RR, const std::vector<double> rr, const std::vector<double> xi, const double rAPP=0., const double r0=-1., const double gamma=1.); 

  /**
   *  @brief the double barred correlation function
   *
   *  (see e.g. Hamilton 1992)
   *
   *  \f[ \overline{\overline{\xi}}(r) = \frac{5}{r^5}\int^r_0dr'\xi(r')r'{^4}
   *  \f]
   *
   *  @param RR the comoving separation, r, where the barred
   *  correlation function is computed
   *
   *  @param rr std::vector containing the input comoving separations
   *
   *  @param xi std::vector containing the input two-point correlation
   *  function
   *
   *  @param rAPP comoving scale below which a power-law model for the
   *  two-point correlation function is assumed in the integral
   *
   *  @param r0 the power-law normalization, r<SUB>0</SUB>  
   *
   *  @param gamma the power-law slope, &gamma;
   *
   *  @return \f$ \overline{\overline{\xi}}(r) \f$
   */
  double barred_xi__ (const double RR, const std::vector<double> rr, const std::vector<double> xi, const double rAPP=0., const double r0=-1., const double gamma=1.); 

  /**
   *  @brief xi<SUB>0</SUB>(s) from &xi;(r,&mu;)
   *
   *  \f[ \xi_0(s) = \frac{1}{2}\int_{-1}^1\xi(s,\mu)d\mu \f]
   *
   *  @param indexR index correspondent to the comoving separation
   *  where the multipole is computed
   *
   *  @param mu std::vector containing the angle between the separation
   *  std::vector and the line of sight
   *
   *  @param xi matrix containing the values of &xi;(r,&mu;)
   *
   *  @return xi<SUB>0</SUB>(s)
   */
  double multipole_xi0 (const int indexR, const std::vector<double> mu, const std::vector<std::vector<double>> xi);
  
  /**
   *  @brief xi<SUB>2</SUB>(s) from &xi;(r,&mu;)
   *
   *  \f[ \xi_2(s) = \frac{5}{2}\int_{-1}^1\xi(s,\mu)L_2(\mu)d\mu \f]
   *
   *  @param indexR index correspondent to the comoving separation
   *  where the multipole is computed
   *
   *  @param mu std::vector containing the angle between the
   *  separation std::vector and the line of sight
   *
   *  @param xi matrix containing the values of &xi;(r,&mu;)
   *
   *  @return xi<SUB>2</SUB>(s)
   */
  double multipole_xi2 (const int indexR, const std::vector<double> mu, const std::vector<std::vector<double>> xi);
 
  /**
   *  @brief xi<SUB>4</SUB>(s) from &xi;(r,&mu;)
   *
   *  \f[ \xi_4(s) = \frac{9}{2}\int_{-1}^1\xi(s,\mu)L_4(\mu)d\mu \f]
   *
   *  @param indexR index correspondent to the comoving separation
   *  where the multipole is computed
   *
   *  @param mu std::vector containing the angle between the separation
   *  std::vector and the line of sight
   *
   *  @param xi matrix containing the values of &xi;(r,&mu;)
   *
   *  @return xi<SUB>4</SUB>(s)
   */
  double multipole_xi4 (const int indexR, const std::vector<double> mu, const std::vector<std::vector<double>> xi);
  
  /**
   *  @brief error on xi<SUB>0</SUB>(s) from &xi;(r,&mu;)
   *
   *  @param indexR index correspondent to the comoving separation
   *  where the multipole is computed
   *
   *  @param mu std::vector containing the angle between the separation
   *  std::vector and the line of sight
   *
   *  @param error matrix containing the errors of &xi;(r,&mu;)
   *
   *  @return error on xi<SUB>0</SUB>(s)
   */
  double error_multipole_xi0 (const int indexR, const std::vector<double> mu, const std::vector<std::vector<double>> error);

  /**
   *  @brief error on xi<SUB>2</SUB>(s) from &xi;(r,&mu;)
   *
   *  @param indexR index correspondent to the comoving separation
   *  where the multipole is computed
   *
   *  @param mu std::vector containing the angle between the separation
   *  std::vector and the line of sight
   *
   *  @param error matrix containing the errors of &xi;(r,&mu;)
   *
   *  @return error on xi<SUB>2</SUB>(s)
   */
  double error_multipole_xi2 (const int indexR, const std::vector<double> mu, const std::vector<std::vector<double>> error);

  /**
   *  @brief error on xi<SUB>4</SUB>(s) from &xi;(r,&mu;)
   *
   *  @param indexR index correspondent to the comoving separation
   *  where the multipole is computed
   *
   *  @param mu std::vector containing the angle between the separation
   *  std::vector and the line of sight
   *
   *  @param error matrix containing the errors of &xi;(r,&mu;)
   *
   *  @return error on xi<SUB>4</SUB>(s)
   */
  double error_multipole_xi4 (const int indexR, const std::vector<double> mu, const std::vector<std::vector<double>> error);

  /**
   *  @brief xi<SUB>0</SUB>(s) from &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  \f[ \xi_0(s) = \frac{1}{2}\int_{-1}^1\xi(s,\mu)d\mu \f]
   *
   *  @param ss comoving scale where the multipole is computed
   *  @param rp std::vector containing the values of r<SUB>p</SUB>
   *  @param pi std::vector containing the values of &pi;
   *  @param xi matrix containing the values of &xi;(r,&mu;)
   *  @param delta_s bin size 
   *  @return xi<SUB>0</SUB>(s)
   */
  double multipole_xi0 (const double ss, const std::vector<double> rp, const std::vector<double> pi, const std::vector<std::vector<double>> xi, const double delta_s);

  /**
   *  @brief xi<SUB>2</SUB>(s) from &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  \f[ \xi_2(s) = \frac{5}{2}\int_{-1}^1\xi(s,\mu)L_2(\mu)d\mu \f]
   *
   *  @param ss comoving scale where the multipole is computed
   *  @param rp std::vector containing the values of r<SUB>p</SUB>
   *  @param pi std::vector containing the values of &pi;
   *  @param xi matrix containing the values of &xi;(r,&mu;)
   *  @param delta_s bin size 
   *  @return xi<SUB>2</SUB>(s)
   */
  double multipole_xi2 (const double ss, const std::vector<double> rp, const std::vector<double> pi, const std::vector<std::vector<double>> xi, const double delta_s);

  /**
   *  @brief xi<SUB>4</SUB>(s) from &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  \f[ \xi_4(s) = \frac{9}{2}\int_{-1}^1\xi(s,\mu)L_4(\mu)d\mu \f]
   *
   *  @param ss comoving scale where the multipole is computed
   *  @param rp std::vector containing the values of r<SUB>p</SUB>
   *  @param pi std::vector containing the values of &pi;
   *  @param xi matrix containing the values of &xi;(r,&mu;)
   *  @param delta_s bin size 
   *  @return xi<SUB>4</SUB>(s)
   */
  double multipole_xi4 (const double ss, const std::vector<double> rp, const std::vector<double> pi, const std::vector<std::vector<double>> xi, const double delta_s);

  /**
   *  @brief error on xi<SUB>0</SUB>(s) from &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  @param ss comoving scale where the multipole is computed
   *  @param rp std::vector containing the values of r<SUB>p</SUB>
   *  @param pi std::vector containing the values of &pi;
   *  @param error matrix containing the errors of &xi;(r,&mu;)
   *  @param delta_s bin size 
   *  @return error on xi<SUB>0</SUB>(s)
   */
  double error_multipole_xi0 (const double ss, const std::vector<double> rp, const std::vector<double> pi, const std::vector<std::vector<double>> error, const double delta_s);

  /**
   *  @brief error on xi<SUB>2</SUB>(s) from &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  @param ss comoving scale where the multipole is computed
   *  @param rp std::vector containing the values of r<SUB>p</SUB>
   *  @param pi std::vector containing the values of &pi;
   *  @param error matrix containing the errors of &xi;(r,&mu;)
   *  @param delta_s bin size 
   *  @return error on xi<SUB>2</SUB>(s)
   */
  double error_multipole_xi2 (const double ss, const std::vector<double> rp, const std::vector<double> pi, const std::vector<std::vector<double>> error, const double delta_s);

  /**
   *  @brief error on xi<SUB>4</SUB>(s) from &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  @param ss comoving scale where the multipole is computed
   *  @param rp std::vector containing the values of r<SUB>p</SUB>
   *  @param pi std::vector containing the values of &pi;
   *  @param error matrix containing the errors of &xi;(r,&mu;)
   *  @param delta_s bin size 
   *  @return error on xi<SUB>4</SUB>(s)
   */
  double error_multipole_xi4 (const double ss, const std::vector<double> rp, const std::vector<double> pi, const std::vector<std::vector<double>> error, const double delta_s);

  /// @cond glob
  /**
   *  @brief multipoles (&xi;<SUB>0</SUB> + &xi;<SUB>2</SUB>) of the
   *  two-point correlation function used in the &chi;<SUP>2</SUP>
   *
   *  @param rr the comoving separation
   *  @param pp a void pointer 
   *  @param par a std::vector containing the coefficients
   *
   *  @return multipoles (&xi;<SUB>0</SUB> + &xi;<SUB>2</SUB>) of the
   *  two-point correlation
   */
  double multipoles (double rr, std::shared_ptr<void> pp, std::vector<double> par);
  /// @endcond

  /**
   *  @brief the model multipole &xi;<SUB>0</SUB> of the two-point
   *  correlation function
   *
   *  @param beta &beta;=f/b, where f is the linear growth rate and b
   *  is the bias 
   *
   *  @param xi_real the real-space two-point correlation function
   *
   *  @return the multipole &xi;<SUB>0</SUB>
   */
  double multipole_xi0_model (const double beta, const double xi_real);

  /**
   *  @brief the model multipole &xi;<SUB>0</SUB> of the two-point
   *  correlation function
   *
   *  @param f_sigma8 f*&sigma;<SUB>8</SUB>
   *  
   *  @param bias_sigma8 b*&sigma;<SUB>8</SUB>
   *
   *  @param sigma8z &sigma;<SUB>8</SUB>
   *  
   *  @param xi_DM the real-space two-point correlation function of
   *  the dark matter
   *
   *  @return the multipole &xi;<SUB>0</SUB>
   */
  double multipole_xi0_model (const double f_sigma8, const double bias_sigma8, const double sigma8z, const double xi_DM);

  /// @cond glob
  /**
   *  @brief the model multipole &xi;<SUB>0</SUB> of the two-point
   *  correlation function
   *
   *  @param xx coordinate x
   *  @param pp void pointer
   *  @param par std::vector containing one or two parameters
   *
   *  @return the multipole &xi;<SUB>0</SUB>
   */
  double multipole_xi0_model (double xx, std::shared_ptr<void> pp, std::vector<double> par);
  /// @endcond

  /**
   *  @brief the model multipole &xi;<SUB>2</SUB> of the two-point
   *  correlation function
   *
   *  @param beta &beta;=f/b, where f is the linear growth rate and b
   *  is the bias 
   *
   *  @param xi_real the real-space two-point correlation function
   *
   *  @param xi_ \f$ \overline{\xi}(r) \f$
   *
   *  @return the multipole &xi;<SUB>2</SUB>
   */
  double multipole_xi2_model (const double beta, const double xi_real, const double xi_); 

  /**
   *  @brief the model multipole &xi;<SUB>4</SUB> of the two-point
   *  correlation function
   *
   *  @param beta &beta;=f/b, where f is the linear growth rate and b
   *  is the bias 
   *
   *  @param xi_real the real-space two-point correlation function
   *
   *  @param xi_ \f$ \overline{\xi}(r) \f$
   *
   *  @param xi__ \f$ \overline{\overline{\xi}}(r) \f$
   *
   *  @return the multipole &xi;<SUB>4</SUB>
   */
  double multipole_xi4_model (const double beta, const double xi_real, const double xi_, const double xi__);

  /// @cond glob
  // theoretical model for the linear xi(rp,pi)
  double xi2D_lin_model (double, double, std::shared_ptr<void> , std::vector<double>);
  /// @endcond

  /**
   *  @brief the linear dispersion model for
   *  &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  @param beta &beta;=f/b, where f is the linear growth rate and b
   *  is the bias
   *
   *  @param bias the bias
   *
   *  @param xi_real the real-space correlation function
   *
   *  @param xi_ \f$ \overline{\xi}(r) \f$
   *
   *  @param xi__ \f$ \overline{\overline{\xi}}(r) \f$
   *
   *  @param P_2 the Legendre polynomial P<SUB>2</SUB>
   *
   *  @param P_4 the Legendre polynomial P<SUB>4</SUB>
   *
   *  @return &xi;(r<SUB>p</SUB>,&pi;)
   */
  double xi2D_lin_model (const double beta, const double bias, const double xi_real, const double xi_, const double xi__, const double P_2, const double P_4);

  /**
   *  @brief the linear dispersion model for
   *  &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  @param rp r<SUB>p</SUB>: comoving separation perpendicular to
   *  the line-of-sight
   *
   *  @param pi &pi;: comoving separation parallel to the
   *  line-of-sight
   *
   *  @param beta &beta;=f/b, where f is the linear growth rate and b
   *  is the bias
   *
   *  @param bias the bias  
   *
   *  @param rad_real std::vector containing the binnend values of the
   *  comoving separations
   *
   *  @param xi_real std::vector containing the binnend values of the
   *  real-space correlation function 
   *
   *  @param xi_ std::vector containing the binnend values of \f$
   *  \overline{\xi}(r) \f$
   *
   *  @param xi__ std::vector containing the binnend values of \f$
   *  \overline{\overline{\xi}}(r) \f$
   *
   *  @param index index for internal use
   *
   *  @param bias_nl 0 \f$ \rightarrow \f$ linear bias; \f$ \rightarrow \f$ 1 non-linear bias 
   *
   *  @param bA the parameter b<SUB>A</SUB> used to model the bias
   *
   *  @return &xi;(r<SUB>p</SUB>,&pi;)
   */
  double xi2D_lin_model (const double rp, const double pi, const double beta, const double bias, const std::vector<double> rad_real, const std::vector<double> xi_real, const std::vector<double> xi_, const std::vector<double> xi__, const int index=-1, const bool bias_nl=0, const double bA=0.);

  /// @cond glob
  // dispersion model for xi(rp,pi)
  double xi2D_model (double, double, std::shared_ptr<void>, std::vector<double>);
  /// @endcond

  /**
   *  @brief the non-linear dispersion model for
   *  &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  @param rp r<SUB>p</SUB>: comoving separation perpendicular to
   *  the line-of-sight
   *
   *  @param pi &pi;: comoving separation parallel to the
   *  line-of-sight
   *
   *  @param beta &beta;=f/b, where f is the linear growth rate and b
   *  is the bias
   *
   *  @param bias the bias
   *
   *  @param sigma12 &sigma;<SUB>12</SUB>
   *
   *  @param rad_real std::vector containing the binnend values of the
   *  comoving separations
   *
   *  @param xi_real std::vector containing the binnend values of the
   *  real-space correlation function 
   *
   *  @param xi_ std::vector containing the binnend values of \f$
   *  \overline{\xi}(r) \f$
   *
   *  @param xi__ std::vector containing the binnend values of \f$
   *  \overline{\overline{\xi}}(r) \f$
   *
   *  @param var 1/[H(z)a(z)]
   *
   *  @param FV 0 \f$ \rightarrow \f$ exponential; \f$ \rightarrow \f$
   *  Gaussian
   *
   *  @param index index for internal use
   *
   *  @param bias_nl 0 \f$ \rightarrow \f$ linear bias; \f$
   *  \rightarrow \f$ non-linear bias
   *
   *  @param bA the parameter b<SUB>A</SUB> used to model the bias
   *
   *  @param v_min the minimum value of the velocity used in the
   *  convolution
   *
   *  @param v_max the maximum value of the velocity used in the
   *  convolution
   *
   *  @param step_v the step of the convolution integral
   *
   *  @return &xi;(r<SUB>p</SUB>,&pi;)
   */
  double xi2D_model (const double rp, const double pi, const double beta, const double bias, const double sigma12, const std::vector<double> rad_real, const std::vector<double> xi_real, const std::vector<double> xi_, const std::vector<double> xi__, const double var, const int FV, int index=-1, const bool bias_nl=0, const double bA=0., const double v_min=-3000., const double v_max=3000., const int step_v=500);
  
  /**
   *  @brief the linear dispersion model for
   *  &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  @param rp r<SUB>p</SUB>: comoving separation perpendicular to
   *  the line-of-sight
   *
   *  @param pi &pi;: comoving separation parallel to the
   *  line-of-sight
   *
   *  @param beta &beta;=f/b, where f is the linear growth rate and b
   *  is the bias
   *
   *  @param bias the bias  
   *
   *  @param funcXiR pointer to an object of type FuncGrid, to
   *  interpolate on \f$ \xi(r) \f$
   *
   *  @param funcXiR_ pointer to an object of type FuncGrid, to
   *  interpolate on \f$ \overline{\xi}(r) \f$
   *
   *  @param funcXiR__ pointer to an object of type FuncGrid, to
   *  interpolate on \f$ \overline{\overline{\xi}} (r) \f$
   *
   *  @param bias_nl 0 \f$ \rightarrow \f$ linear bias; \f$ \rightarrow \f$ 1 non-linear bias 
   *
   *  @param bA the parameter b<SUB>A</SUB> used to model the bias
   *
   *  @return &xi;(r<SUB>p</SUB>,&pi;)
   */
  double xi2D_lin_model (const double rp, const double pi, const double beta, const double bias,  const std::shared_ptr<void> funcXiR, const std::shared_ptr<void> funcXiR_, const std::shared_ptr<void> funcXiR__, const bool bias_nl=0, const double bA=0.);

  /**
   *  @brief the non-linear dispersion model for
   *  &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  @param rp r<SUB>p</SUB>: comoving separation perpendicular to
   *  the line-of-sight
   *
   *  @param pi &pi;: comoving separation parallel to the
   *  line-of-sight
   *
   *  @param beta &beta;=f/b, where f is the linear growth rate and b
   *  is the bias
   *
   *  @param bias the bias
   *
   *  @param sigma12 &sigma;<SUB>12</SUB>
   *
   *  @param funcXiR pointer to an object of type FuncGrid, to
   *  interpolate on \f$ \xi(r) \f$
   *
   *  @param funcXiR_ pointer to an object of type FuncGrid, to
   *  interpolate on \f$ \overline{\xi}(r) \f$
   *
   *  @param funcXiR__ pointer to an object of type FuncGrid, to
   *  interpolate on \f$ \overline{\overline{\xi}} (r) \f$
   *
   *  @param var 1/[H(z)a(z)]
   *
   *  @param FV 0 \f$ \rightarrow \f$ exponential; \f$ \rightarrow \f$ 1 gaussian 
   *
   *  @param bias_nl 0 \f$ \rightarrow \f$ linear bias; \f$ \rightarrow \f$ 1 non-linear bias 
   *
   *  @param bA the parameter b<SUB>A</SUB> used to model the bias
   *
   *  @param v_min the minimum value of the velocity used in the
   *  convolution
   *
   *  @param v_max the maximum value of the velocity used in the
   *  convolution
   *
   *  @param step_v the step of the convolution integral
   *
   *  @return &xi;(r<SUB>p</SUB>,&pi;)
   */
  double xi2D_model (const double rp, const double pi, const double beta, const double bias, const double sigma12, const std::shared_ptr<void> funcXiR, const std::shared_ptr<void> funcXiR_, const std::shared_ptr<void> funcXiR__, const double var, const int FV, const bool bias_nl=0, const double bA=0., const double v_min=-3000., const double v_max=3000., const int step_v=500);

  /**
   *  @brief pairwise velocity distribution
   *
   *  @param vel comoving velocity
   *
   *  @param sigma12 &sigma;<SUB>12</SUB>
   *
   *  @param FV 0 \f$ \rightarrow \f$ exponential; \f$ \rightarrow \f$
   *  1 gaussian
   *
   *  @return f(v)
   */
  double f_v (const double vel, const double sigma12, const int FV);

  /**
   *  @brief pairwise velocity distribution
   *
   *  (see Chuang&Wang 2012 (1209.0210), eq. 25)
   *
   *  @param vel comoving velocity
   *
   *  @param rp r<SUB>p</SUB>: comoving separation perpendicular to
   *  the line-of-sight
   *
   *  @param pi &pi;: comoving separation parallel to the
   *  line-of-sight
   *
   *  @param var 1/[H(z)a(z)]
   *  @param sigmav0 &sigma;<SUB>v,0</SUB>
   *  @param cmu C<SUB>&mu;</SUB>
   *  @param cs1 C<SUB>&sigma;1</SUB>
   *  @param cs2 C<SUB>&sigma;2</SUB>
   *
   *  @return f(v)
   */
  double f_v (const double vel, const double rp, const double pi, const double var, const double sigmav0, const double cmu, const double cs1, const double cs2);

  /**
   *  @brief velocity distribution used to model BAO
   *  
   *  (see Chuang&Wang 2012, 1209.0210)
   *
   *  @param xx s
   *  @param f_g f<SUB>g</SUB>
   *  @param k_star k<SUB>*</SUB>
   *
   *  @return f<SUB>*</SUB>
   */
  double f_star (const double xx, const double f_g, const double k_star);

  /**
   *  @brief a possible parameterization of the non-linear bias
   *  
   *  the non-linear bias (see Chuang&Wang 2012 (1209.0210),
   *  eqs. 20-21): \f$ b(r) = r^{\frac{b_A}{1+(\frac{r}{b_B})^{b_C}}}
   *  \f$
   *
   *  @param rr the comoving scale
   *  @param bA b<SUB>A</SUB>
   *  @param bB b<SUB>B</SUB>
   *  @param bC b<SUB>C</SUB>
   *
   *  @return b(r)
   */
  double b_nl (const double rr, const double bA, const double bB=10., const double bC=4.);

  /**
   *  @brief estimated relative error on \f$\beta=f/b\f$
   *
   *  \f$
   *  \delta\beta/\beta\sim Cb^{0.7}V^{-0.5}
   *  \exp\left[n_0/(b^2n)\right]
   *  \f$
   *
   *  where \f$n_0=1.7\times10^{-4} h^3 Mpc^{-3}\f$ and
   *  \f$C=4.9\times10^2 h^{-1.5} Mpc^{1.5}\f$ 
   *
   *  see Eq. 20 of Bianchi et al. 2012, http://arxiv.org/abs/1203.1545
   *
   *  @param bias the linear galaxy bias, b
   *  @param Volume the survey volume, V
   *  @param density the galaxy density, n
   *  @return \f$\delta\beta/\beta\f$
   */
  double relative_error_beta (const double bias, const double Volume, const double density); 

  /**
   * @brief integrand of the 2d power spectrum to obtain power
   * spectrum multipole
   *
   * @param mu the value mu
   * @param parameters the parameters for the integration
   * @return the power spectrum multipoles integrand
   */
  double Pkl_Kaiser_integrand(const double mu, void *parameters);

  /**
   * @brief integrand of the 2d power spectrum to obtain sigma^2(k)
   * @param mu the value mu
   * @param parameters the parameters for the integration
   * @return the sigma2 integrand
   */
  double sigma2_integrand (const double mu, void *parameters);

  /**
   * @brief integrand to obtain the 2PCF multipoles
   * @param kk the value kk
   * @param parameters the parameters for the integration
   * @return the 2pcf multipoles integrand
   */
  double covariance_XiMultipoles_integrand (const double kk, void *parameters);

  /**
   * @brief integrand to obtain covariance for the 2PCF multipoles
   * @param kk the value kk
   * @param parameters the parameters for the integration
   * @return the covariance of 2pcf multipoles integrand
   */
  double XiMultipoles_integrand (const double kk, void *parameters);

  /**
   * @brief integrand to obtain the 2PCF multipoles from 2D 2pcf
   * in polar coordinates
   * @param mu the value mu
   * @param parameters the parameters for the integration
   * @return the 2pcf multipoles integrand
   */
  double XiMultipoles_from_Xi2D_integrand (const double mu, void *parameters);

  /**
   * @brief function to obtain the Kaiser factor 
   * @param order the expansion order
   * @param bias the bias factor
   * @param f the linear growth factor
   * @return the Kaiser integral
   */
  double Pkl_Kaiser_integral(const int order, const double bias, const double f);

  /**
   * @brief function to obtain the linear RSD power spectrum
   * monopole
   * @param kk the scales k
   * @param Pk the power spectrum
   * @param bias the bias factor
   * @param f the linear growth factor
   * @return the linear RSD power spectrum monopole
   */
  std::vector<double> Pk0_Kaiser(const std::vector<double> kk, const std::vector<double> Pk, const double bias, const double f);

  /**
   * @brief function to obtain the linear RSD 
   * power spectrum quadrupole
   * @param kk the scales k
   * @param Pk the power spectrum
   * @param bias the bias factor
   * @param f the linear growth factor
   * @return the linear RSD power spectrum quadrupole
   */
  std::vector<double> Pk2_Kaiser(const std::vector<double> kk, const std::vector<double> Pk, const double bias, const double f);

  /**
   * @brief function to obtain the linear RSD 
   * power spectrum hexadecapole
   * @param kk the scales k
   * @param Pk the power spectrum
   * @param bias the bias factor
   * @param f the linear growth factor
   * @return the linear RSD power spectrum hexadecapole
   */
  std::vector<double> Pk4_Kaiser(const std::vector<double> kk, const std::vector<double> Pk, const double bias, const double f);

  /**
   * @brief function to obtain Pk multipoles from linear RSD (Kaiser)
   * @param orders the l-th multipole desired
   * @param kk the scales k
   * @param Pk the power spectrum
   * @param bias the bias factor
   * @param f the linear growth factor
   * @return the power spectrum multipoles
   */
  std::vector<std::vector<double>> Pkl_Kaiser(const std::vector<int> orders, const std::vector<double> kk, const std::vector<double> Pk, const double bias, const double f);

  /**
   * @brief function to obtain the two point correlation
   * funciton monopole
   * @param r the scales r
   * @param kk the scales k
   * @param Pk0 the power spectrum monopole
   * @param k_cut the k scale to cut the integrand
   * @param cut_pow the power of the integrand cut
   * @param IntegrationMethod method for integration
   *
   * @return the linear RSD power spectrum monopole
   */
  std::vector<double> Xi0(const std::vector<double> r, const std::vector<double> kk, const std::vector<double> Pk0, const double k_cut=0.7, const double cut_pow=2, const int IntegrationMethod = 1);

  /**
   * @brief function to obtain the two point correlation
   * function quadrupole
   * @param rr the scales r
   * @param kk the scales k
   * @param Pk2 the power spectrum quadrupole
   * @param k_cut the k scale to cut the integrand
   * @param cut_pow the power of the integrand cut
   * @param IntegrationMethod method for integration
   *
   * @return the linear RSD power spectrum quadrupole
   */
  std::vector<double> Xi2 (const std::vector<double> rr, const std::vector<double> kk, const std::vector<double> Pk2, const double k_cut=0.58, const double cut_pow=4, const int IntegrationMethod = 1);

  /**
   * @brief function to obtain the two point correlation 
   * function hexadecapole
   * @param rr the scales r
   * @param kk the scales k
   * @param Pk4 the power spectrum hexadecapole
   * @param k_cut the k scale to cut the integrand
   * @param cut_pow the power of the integrand cut
   * @param IntegrationMethod method for integration
   *
   * @return the linear RSD power spectrum hexadecapole
   */
  std::vector<double> Xi4 (const std::vector<double> rr, const std::vector<double> kk, const std::vector<double> Pk4, const double k_cut=0.6, const double cut_pow=2, const int IntegrationMethod = 1);

  /**
   * @brief function to obtain the monopole and
   * quadrupole of the two point correlation function 
   *
   * @param alpha_perpendicular the shift along the line of sight
   * @param alpha_parallel the shift parallel to the line of sight
   * @param rr the scales r
   * @param rl the scales at which the multipoles are defined
   * @param Xi0 the 2pfc monopole
   * @param Xi2 the 2pfc quadrupole
   *
   * @return the monopole and quadrupole 
   * of the two point correlation function 
   */
  std::vector<std::vector<double>> Xi02_AP (const double alpha_perpendicular, const double alpha_parallel, const std::vector<double> rr, const std::vector<double> rl, const std::vector<double> Xi0, const std::vector<double> Xi2);

  /**
   * @brief function to obtain the monopole, quadrupole
   * and hexadecapole of the two-point correlation function 
   * 
   * @param alpha_perpendicular the shift along the line of sight
   * @param alpha_parallel the shift parallel to the line of sight
   * @param rr the scales r
   * @param rl the scales at which the multipoles are defined
   * @param Xi0 the 2pfc monopole
   * @param Xi2 the 2pfc quadrupole
   * @param Xi4 the 2pfc hexadecapole
   *
   * @return the monopole, quadrupole and hexadecapole of the two
   * point correlation function
   */
  std::vector<std::vector<double>> Xi024_AP (const double alpha_perpendicular, const double alpha_parallel, const std::vector<double> rr, const std::vector<double> rl, const std::vector<double> Xi0, const std::vector<double> Xi2, const std::vector<double> Xi4);

  /**
   * @brief function to obtain the 2pcf wedges
   *
   * @param mu_min the lower limit of integration for wedges
   * @param delta_mu the mu width for wedges
   * @param alpha_perpendicular the shift along the line of sight
   * @param alpha_parallel the shift parallel to the line of sight
   * @param rr the scales r
   * @param rl the scales at which the multipoles are defined
   * @param Xi0 the 2pfc monopole
   * @param Xi2 the 2pfc quadrupole
   * @param Xi4 the 2pfc hecadecapole
   *
   * @return the 2pcf wedges 
   */
  std::vector<std::vector<double>> XiWedges_AP (const std::vector<double> mu_min, const std::vector<double> delta_mu, const double alpha_perpendicular, const double alpha_parallel, const std::vector<double> rr, const std::vector<double> rl, const std::vector<double> Xi0, const std::vector<double> Xi2, const std::vector<double> Xi4);

  /**
   * @brief multipole expansion of the per-mode covariance sigma2_k
   * (see Grieb et al. 2016, eq. 15 https://arxiv.org/pdf/1509.04293)
   *
   * \f[
   *
   *   \sigma_{\ell_{1} \ell_{2}}^{2}(k) \equiv \frac{\left(2
   *   \ell_{1}+1\right)\left(2 \ell_{2}+1\right)}{V_{\mathrm{s}}}
   *   \times \int_{-1}^{1}\left[P(k,
   *   \mu)+\frac{1}{\bar{n}}\right]^{2} \mathcal{L}_{\ell_{1}}(\mu)
   *   \mathcal{L}_{\ell_{2}}(\mu) \mathrm{d} \mu
   *
   * \f]
   *
   * where \f$P(k, \mu)\f$ is the polar power spectrum computed from
   * input power spectrum multipoles and \f$\mathcal{L}\f$ are the
   * Legendre polynomials.
   *
   * @param nObjects number of objects in the sample
   * @param Volume the sample volume
   * @param kk the scales kk
   * @param Pk_multipoles the power spectrum multipoles 
   * @param orders the power spectrum multipoles orders
   * @return the sigma2_k (see i.e. Grieb et al. 2016, eq. 15)
   */
  std::vector< std::vector<double>> sigma2_k (const double nObjects, const double Volume, const std::vector<double> kk, const std::vector<std::vector<double>> Pk_multipoles, const std::vector<int> orders);

  /**
   * @brief Covariance matrix for two-point correlation multipoles
   * (see Grieb et al. 2016, Eq. 18 https://arxiv.org/pdf/1509.04293)
   *
   * \f[ C_{\ell_{1} \ell_{2}}^{\epsilon}\left(s_{i},
   *  s_{j}\right)=\frac{\mathrm{i}^{\ell_{1}+\ell_{2}}}{2 \pi^{2}}
   *  \int_{0}^{\infty} k^{2} \sigma_{\ell_{1} \ell_{2}}^{2}(k)
   *  \bar{\jmath}_{\ell_{1}}\left(k s_{i}\right)
   *  \bar{\jmath}_{\ell_{2}}\left(k s_{j}\right) \mathrm{d} k \f]
   *
   * where \f$\sigma_{\ell_{1} \ell_{2}}^{2}(k)\f$ is the multipole
   * expansion of the per-mode covariance sigma2_k, computed by
   * cbl::sigma2_k and \f$\bar{\jmath}\f$ is the shell-averaged Bessel
   * function, computed by cbl::jl_distance_average.
   *
   * @param rr output scales
   * @param covariance analytic covariance matrix
   * @param nbins number of configuration space bins
   * @param rMin minimum configuration space scale
   * @param rMax maximum configuration space scale
   * @param nObjects number of objects in the sample
   * @param Volume the sample volume
   * @param kk the scales kk
   * @param Pk_multipoles the power spectrum multipoles 
   * @param orders the power spectrum multipoles orders
   * @param bin_type the bin type
   * @return none
   */
  void Covariance_XiMultipoles (std::vector<double> &rr, std::vector<std::vector<double>> &covariance, const int nbins, const double rMin, const double rMax, const double nObjects, const double Volume, const std::vector<double> kk, const std::vector<std::vector<double>> Pk_multipoles, const std::vector<int> orders, const cbl::BinType bin_type=cbl::BinType::_linear_);

  /**
   * @brief Covariance matrix for two-point correlation wedges (see
   * Grieb et al. 2016, Eq. 19 (https://arxiv.org/pdf/1509.04293)
   *
   * \f[ C_{\mu \mu^{\prime}}^{\xi}\left(s_{i}, s_{j}\right)=
   *    \sum_{\ell_{1}, \ell_{2}}
   *    \frac{\mathrm{i}^{\ell_{1}+\ell_{2}}}{2 \pi^{2}}
   *    \overline{\mathcal{L}}_{\ell_{1}, \mu}
   *    \overline{\mathcal{L}}_{\ell_{2}, \mu^{\prime}} \times
   *    \int_{0}^{\infty} k^{2} \sigma_{\ell_{1} \ell_{2}}^{2}(k)
   *    \bar{\jmath}_{\ell_{1}}\left(k s_{i}\right)
   *    \bar{\jmath}_{\ell_{2}}\left(k s_{j}\right) \mathrm{d} k \f]
   *
   * where \f$\sigma_{\ell_{1} \ell_{2}}^{2}(k)\f$ is the multipole
   * expansion of the per-mode covariance computed by cbl::sigma2_k,
   * \f$\bar{\jmath}\f$ is the shell-averaged Bessel function,
   * computed by cbl::jl_distance_average an
   * \f$\overline{\mathcal{L}}\f$ is the average of the Legendre
   * polynomial over the \f$\mu\f$ bin, computed by
   * cbl::Legendre_polynomial_mu_average.
   *
   * @param rr output scales
   * @param covariance analytic covariance matrix
   * @param mu the lower wedge integratin limit
   * @param delta_mu the wedge mu bin size
   * @param nbins number of configuration space bins
   * @param rMin minimum configuration space scale
   * @param rMax maximum configuration space scale
   * @param nObjects number of objects in the sample
   * @param Volume the sample volume
   * @param kk the scales kk
   * @param Pk_multipoles the power spectrum multipoles 
   * @param orders the power spectrum multipoles orders
   * @param bin_type the bin type
   * @return none
   */
  void Covariance_XiWedges (std::vector<double> &rr, std::vector<std::vector<double>> &covariance, const std::vector<double> mu, const std::vector<double> delta_mu, const int nbins, const double rMin, const double rMax, const double nObjects, const double Volume, const std::vector<double> kk, const std::vector<std::vector<double>> Pk_multipoles, const std::vector<int> orders, const cbl::BinType bin_type=cbl::BinType::_linear_);

  ///@}


// ============================================================================================


  /**
   *  @brief The namespace of the functions and classes of <B>
   *  internal auxiliary use </B>
   *  
   *  The \e glob namespace contains all classes, structures and
   *  functions of internal auxiliary use
   */
  namespace glob {

    struct STR_grid
    {
      std::vector<double> _xx, _yy;
    };

    struct STR_grid_2D
    {
      std::vector<double> _xx1, _xx2;
      std::vector<std::vector<double>> _yy;
    };

    struct STR_xi0_model
    {
      double bias_sigma8;
      double sigma8z;
      std::vector<double> xi_DM;
    };

    struct STR_xi2D_model
    {
      std::vector<double> rp, pi, xi_real, xi_, xi__, P2, P4, vel;
      std::vector<int> lim_index_fit, type;
      int step_v, FV, dim;
      double delta_v;
      bool bias_nl;
    };

    struct STR_xi 
    {
      double rr, aa;
      std::vector<double> lgkk, lgPk;
    };

    struct STR_SSM 
    {
      int unit;
      double hh, mass, rho, n_spec;
      std::vector<double> lgkk, lgPk;
    };
  
    struct STR_jl_distance_average
    {
      int order;
      double k;
    };
  
    struct STR_Pl_mu_integral
    {
      int order;
    };

    struct STR_Pkl_Kaiser_integrand
    {
      int l;
      double Pk, bias, f;
    };

    struct STR_closest_probability
    {
      std::vector<double> values;
      std::vector<double> weights;
    };

    struct STR_distribution_probability
    {
      std::shared_ptr<cbl::glob::FuncGrid> func;
    };

    struct STR_sigma2_integrand
    {
      int l1, l2;
      double density_inv,kk;
      std::vector<int> orders;
      std::vector<FuncGrid> Pk_multipoles_interp;
    };

    struct STR_XiMultipoles_integrand
    {
      double r;
      int l;
      FuncGrid *Pkl;
      double k_cut;
      double cut_pow;
    };

    struct STR_xi2D_smu_integrand
    {
      FuncGrid *func;
      int order;
    };

    struct STR_covariance_XiMultipoles_integrand
    {
      FuncGrid *s2, *jl1r1, *jl2r2;
    };
  }
  

  /**
   *  @name Functions to model the correlation function
   */
  ///@{
  
  /**
   * @brief function to obtain the monopole and
   * quadrupole of the two point correlation function 
   *
   * @param alpha_perpendicular the shift along the line of sight
   * @param alpha_parallel the shift parallel to the line of sight
   * @param rr the scales r
   * @param xi0_interp the xi0 interpolator
   * @param xi2_interp the xi2 interpolator
   *
   * @return the monopole and quadrupole 
   * of the two point correlation function 
   */
  std::vector<std::vector<double>> Xi02_AP (const double alpha_perpendicular, const double alpha_parallel, const std::vector<double> rr, const std::shared_ptr<glob::FuncGrid> xi0_interp, const std::shared_ptr<glob::FuncGrid> xi2_interp);

  /**
   * @brief function to obtain the monopole, quadrupole
   * and hexadecapole of the two-point correlation function 
   * 
   * @param alpha_perpendicular the shift along the line of sight
   * @param alpha_parallel the shift parallel to the line of sight
   * @param rr the scales r
   * @param xi0_interp the xi0 interpolator
   * @param xi2_interp the xi2 interpolator
   * @param xi4_interp the xi4 interpolator
   *
   * @return the monopole, quadrupole and hexadecapole of the two
   * point correlation function
   */
  std::vector< std::vector<double>> Xi024_AP (const double alpha_perpendicular, const double alpha_parallel, const std::vector<double> rr, const std::shared_ptr<glob::FuncGrid> xi0_interp, const std::shared_ptr<glob::FuncGrid> xi2_interp, const std::shared_ptr<glob::FuncGrid> xi4_interp);

  /**
   * @brief function to obtain the 2pcf wedges
   *
   * @param mu_min the lower limit of integration for wedges
   * @param delta_mu the mu width for wedges
   * @param alpha_perpendicular the shift along the line of sight
   * @param alpha_parallel the shift parallel to the line of sight
   * @param rr the scales r
   * @param xi0_interp the xi0 interpolator
   * @param xi2_interp the xi2 interpolator
   * @param xi4_interp the xi4 interpolator
   *
   * @return the 2pcf wedges 
   */
  std::vector<std::vector<double>> XiWedges_AP (const std::vector<double> mu_min, const std::vector<double> delta_mu, const double alpha_perpendicular, const double alpha_parallel, const std::vector<double> rr, const std::shared_ptr<glob::FuncGrid> xi0_interp, const std::shared_ptr<glob::FuncGrid> xi2_interp, const std::shared_ptr<glob::FuncGrid> xi4_interp);

  ///@}

}

#endif
