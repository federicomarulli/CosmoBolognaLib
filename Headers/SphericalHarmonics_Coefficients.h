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
 *  @file Headers/SphericalHarmonics_Coefficients.h
 *
 *  @brief Generic functions that use one or more classes of the
 *  CosmoBolognaLib
 *
 *  This file contains the prototypes of the class SphericalHarmonics_Coefficients
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __SHCOEFF__
#define __SHCOEFF__

#include "Func.h"

namespace cbl {

 namespace glob {

   /// Eigen 4D array
   typedef Eigen::Array<double, 4, 1> dVector;

   /**
    *  @enum SpHarMethod
    *  @brief method to compute the spherical harmonics
    */
   enum class SpHarMethod {

     /// use standard method to compute
     _STANDARD_,

     /// Use GSL to compute spherical harmonics
     _GSL_,

     /// Use eigen to vectorize the spherical harmonics computation
     _EIGEN_

   };

   /**
    * @brief return a vector containing the
    * SpHarMethod names
    * @return a vector containing the
    * SpHarMethod names
    */
   inline std::vector<std::string> SpHarMethodNames () { return {"Standard", "GSL", "Eigen"}; }

   /**
    * @brief struct to manage the spherical
    * harmonics computation
    */
   struct Alm {

     /// \f$l\f$ order
     int ell = 0;

     /// \f$m\f$ order
     int emm = 0;

     /// normalization factor
     double Norm = 0;

     /// value of associated Legendre polynomial \f$P_{l,m}(\theta)\f$
     double Plm = 0;

     /// \f$ \cos(m*\phi) \f$
     double realPart = 0;

     /// \f$ \sin(m*\phi) \f$
     double imagPart = 0;

     /*
      * @brief reset the internal variables
      * @return None
      */
     void reset()
     {
       Plm = 0;
       realPart = 0;
       imagPart = 0;
     }

     /*
      * @brief set the internal variables
      * @param _Plm <double>, legendre polynomial value
      * @param real <double>, the real part
      * @param imag <double>, the imaginary part
      * @return None
      */
     void set(const double _Plm, const double real=1., const double imag=0.)
     {
       Plm = _Plm;
       realPart = Plm*real;
       imagPart = Plm*imag;
     }

     /*
      * @brief add input to the internal variables
      * @param _Plm legendre polynomial value
      * @param real <double>, the real part
      * @param imag <double>,the imaginary part
      * @return None
      */
     void add(const double _Plm, const double real=1., const double imag=0.)
     {
       Plm = _Plm;
       realPart += Plm*real;
       imagPart += Plm*imag;
     }

     /*
      * @brief return the spherical harmonics real part
      * @return the spherical harmonics real part
      */
     double real() { return Norm*realPart;}

     /*
      * @brief return the spherical harmonics imaginary part
      * @return the spherical harmonics imaginary part
      */
     double imag() { return Norm*imagPart;}

     /*
      * @brief return the spherical harmonics unnormalized real part
      * @return the spherical harmonics unnormalized real part
      */
     double realUnnorm() { return realPart;}

     /*
      * @brief return the spherical harmonics unnormalized imaginary part
      * @return the spherical harmonics unnormalized imaginary part
      */
     double imagUnnorm() { return imagPart;}

     /**
      * @brief function to print alm
      * @return none
      */
     void print() { std::cout << ell<< " " << emm << " " << Norm*realPart << " " << Norm*imagPart << std::endl; }
   };

   /**
    * @brief struct to manage the spherical
    * harmonics computation using Eigen
    */
   struct AlmEigen {
     EIGEN_MAKE_ALIGNED_OPERATOR_NEW

       /// \f$l\f$ order
       int ell = 0;

     /// \f$m\f$ order
     int emm = 0;

     /// normalization factor
     dVector Norm;

     /// value of associated Legendre polynomial \f$P_{l,m}(\theta)\f$
     dVector Plm;

     /// \f$ \cos(m*\phi) \f$
     dVector realPart;

     /// \f$ \sin(m*\phi) \f$
     dVector imagPart;

     /*
      * @brief reset the internal variables
      * @return None
      */
     void reset()
     {
       Plm = {0, 0, 0, 0};
       realPart = {0, 0, 0, 0};
       imagPart = {0, 0, 0, 0};
     }

     /*
      * @brief set the internal variables
      * @param _Plm <double>, legendre polynomial value
      * @return None
      */
     void set(const dVector &_Plm)
     {
       Plm = _Plm;
       realPart = Plm;
       imagPart = {0, 0, 0, 0};
     }

     /*
      * @brief set the internal variables
      * @param _Plm <double>, legendre polynomial value
      * @param real <double>, the real part
      * @param imag <double>, the imaginary part
      * @return None
      */
     void set(const dVector &_Plm, const dVector &real, const dVector &imag)
     {
       Plm = _Plm;
       realPart = _Plm*real;
       imagPart = _Plm*imag;
     }

     /*
      * @brief add input to the internal variables
      * @param _Plm legendre polynomial value
      * @param real <double>, the real part
      * @param imag <double>,the imaginary part
      * @return None
      */
     void add(const dVector &_Plm, const dVector &real, const dVector &imag)
     {
       Plm = _Plm;
       realPart = realPart+_Plm*real;
       imagPart = imagPart+_Plm*imag; 
     }

     /*
      * @brief return the spherical harmonics real part
      * @return the spherical harmonics real part
      */
     dVector real() { return Norm*realPart;}

     /*
      * @brief return the spherical harmonics imaginary part
      * @return the spherical harmonics imaginary part
      */
     dVector imag() { return Norm*imagPart;}

     /*
      * @brief return the spherical harmonics unnormalized real part
      * @return the spherical harmonics unnormalized real part
      */
     dVector realUnnorm() { return realPart;}

     /*
      * @brief return the spherical harmonics unnormalized imaginary part
      * @return the spherical harmonics unnormalized imaginary part
      */
     dVector imagUnnorm() { return imagPart;}

     /**
      * @brief function to print alm
      * @return none
      */
     void print() { 
       for (size_t i=0; i<4; i++) 
	 std::cout << ell << " " << emm << " " << i << " "  << Norm[i]*realPart[i] << " " << Norm[i]*imagPart[i] << std::endl;
     }
   };

   /**
    * @class SphericalHarmonicsArray_Standard
    * @brief
    *
    */
   class SphericalHarmonicsArray {

     public:
       /**
	* @brief default construtor
	*/
       SphericalHarmonicsArray () {}

       /**
	* @brief Destructor
	*/
       virtual ~SphericalHarmonicsArray() = default;

       /**
	* @brief compute the a_lm
	*
	* @param x the x coordinate
	* @param y the y coordinate
	* @param z the z coordinate
	* @param weight the weight
	*/ 
       virtual void compute (const double &x, const double &y, const double &z, const double weight=1)
       { (void)x; (void)y; (void)z; (void)weight; exit(1);}

       /**
	* @brief compute the a_lm
	*
	* @param x the x coordinate
	* @param y the y coordinate
	* @param z the z coordinate
	* @param weight the weight
	*/ 
       virtual void compute () { exit(1);}

       virtual void print() {exit(1);}

       static std::shared_ptr<SphericalHarmonicsArray> factory(const SpHarMethod method, const int lMax);

       virtual size_t nOrders () const {return m_nOrders; }

     protected:

       /// numbero of Legendre polynomial orders
       size_t m_nOrders;

       /// number of spherical harmonics = \f$(l_{max}+1)*(l_{max}+2)/2\f$
       size_t m_nSpH;

       /**
	* @param set internal variables m_nOrders, m_nSph, m_alm
	* from the maximum order of expansion
	* @param lMax maximum order of the  expansion
	* @return None
	*/
       virtual void m_setAlm(const int lMax) = 0;

       /**
	* @brief reset internal values of m_alm
	* @return None
	*/
       virtual void m_resetAlm() = 0;
   };

   /**
    * @class SphericalHarmonicsArray_Standard
    * @brief
    *
    */
   class SphericalHarmonicsArray_Standard : public SphericalHarmonicsArray {

     public:

       /**
	* @brief default construtor
	*/
       SphericalHarmonicsArray_Standard() {}

       /**
	* @brief default construtor
	*
	* @param lMax maximum order of the  expansion
	*/
       SphericalHarmonicsArray_Standard(const int lMax) : SphericalHarmonicsArray() { m_setAlm(lMax);}

       /**
	* @brief Destructor
	*/
       virtual ~SphericalHarmonicsArray_Standard() = default;

       /**
	* @brief compute the a_lm
	*
	* @param x the x coordinate
	* @param y the y coordinate
	* @param z the z coordinate
	* @param weight the weight
	*/ 
       void compute (const double &x, const double &y, const double &z, const double weight=1);

       Alm operator() (const int i) {return m_alm[i];}

       std::vector<Alm> operator() () {return m_alm;}

       void print();

     protected:

       /// \f$a_{l, m}\f$ vector
       std::vector<Alm> m_alm;

       /// numbero of Legendre polynomial orders
       size_t m_nOrders;

       /// number of spherical harmonics = \f$(l_{max}+1)*(l_{max}+2)/2\f$
       size_t m_nSpH;

       /**
	* @param set internal variables m_nOrders, m_nSph, m_alm
	* from the maximum order of expansion
	* @param lMax maximum order of the  expansion
	* @return None
	*/
       void m_setAlm(const int lMax);

       /**
	* @brief reset internal values of m_alm
	* @return None
	*/
       void m_resetAlm();


   }; /* End of SphericalHarmonicsArray_Standard class */

   /**
    * @class SphericalHarmonicsArray
    * @brief compute spherical harmonics using
    * GSL
    *
    */
   class SphericalHarmonicsArray_GSL : public SphericalHarmonicsArray_Standard {

     public:

       /**
	* @brief default construtor
	*/
       SphericalHarmonicsArray_GSL() {}

       /**
	* @brief default construtor
	*
	* @param lMax maximum order of the  expansion
	*/
       SphericalHarmonicsArray_GSL (const int lMax) : SphericalHarmonicsArray_Standard (lMax) {}

       /**
	* @brief Destructor[i]
	*/
       virtual ~SphericalHarmonicsArray_GSL() = default;

       /**
	* @brief compute the a_lm
	*
	* @param x the x coordinate
	* @param y the y coordinate
	* @param z the z coordinate
	* @param weight the weight
	*/ 
       void compute (const double &x, const double &y, const double &z, const double weight=1);
   };

   /**
    * @class SphericalHarmonicsArray_Standard
    * @brief
    *
    */
   class SphericalHarmonicsArray_Eigen : public SphericalHarmonicsArray {

     public:
       EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	 /**
	  * @brief default construtor
	  */
	 SphericalHarmonicsArray_Eigen() {}

       /**
	* @brief default construtor
	*
	* @param lMax maximum order of the  expansion
	*/
       SphericalHarmonicsArray_Eigen(const int lMax) : SphericalHarmonicsArray() { m_setAlm(lMax);}

       /**
	* @brief Destructor
	*/
       virtual ~SphericalHarmonicsArray_Eigen() = default;

       /**
	* @brief compute the a_lm
	*
	* @param x the x coordinate
	* @param y the y coordinate
	* @param z the z coordinate
	* @param weight the weight
	*/ 
       void compute (const double &x, const double &y, const double &z, const double weight=1);

       /**
	* @brief compute the a_lm whatever the buffer size
	*
	* @return none
	*/ 
       void compute ();

       void print();

     protected:

       dVector m_zero;

       dVector m_one;

       dVector m_xx;

       dVector m_yy;

       dVector m_zz;

       dVector m_weight;

       int m_buffer_size;

       /// \f$a_{l, m}\f$ vector
       std::vector<AlmEigen, Eigen::aligned_allocator<AlmEigen>> m_alm;

       /// numbero of Legendre polynomial orders
       size_t m_nOrders;

       /// number of spherical harmonics = \f$(l_{max}+1)*(l_{max}+2)/2\f$
       size_t m_nSpH;


       void m_buffer_reset ();

       /**
	* @param set internal variables m_nOrders, m_nSph, m_alm
	* from the maximum order of expansion
	* @param lMax maximum order of the  expansion
	* @return None
	*/
       void m_setAlm(const int lMax);

       /**
	* @brief reset internal values of m_alm
	* @return None
	*/
       void m_resetAlm();

   };

   /**
    * @brief compute legendre coefficients
    * using the spherical harmonics addition theorem
    *
    * @param Ylm_1 first spherical harmonics expansion
    * @param Ylm_2 second spherical harmonics expansion
    * @return a vector containing the legendre coefficients
    */
   std::vector<double> legendre_polynomials (const SphericalHarmonicsArray& Ylm_1, const SphericalHarmonicsArray& Ylm_2);

   /**
    *  @class SphericalHarmonics_Coefficients SphericalHarmonics_Coefficients.h "Headers/SphericalHarmonics_Coefficients.h"
    *
    *  @brief The class SphericalHarmonics_Coefficients
    *
    *  This class is used to handle objects of type 
    *  <EM> SphericalHarmonics_Coefficients </EM>. 
    *  It contains all methods to compute coefficients of
    *  spherical harmonics expansion \f$a_{lm}\f$ in any position of the
    *  unit sphere.
    *
    *  Coefficients can be binned according to the magnitude of the separation
    *  vector and accumulated. 
    */
   class SphericalHarmonics_Coefficients {

     protected:

       /// the number of separation bins
       int m_nbins;

       /// the number of multipoles, \f$ l_{max}+1 \f$
       int m_norder;

       /// the maximum multipole \f$ l_{max} \f$
       int m_lmax;

       /// the total number of spherical harmonics
       int m_n_sph;

       /// the number of spherical harmonics for a given choice of \f$l\f$
       std::vector<int> m_n_sph_l;

       /// the spherical harmonics expansion coefficients in separation bins
       std::vector<std::vector<std::complex<double>>> m_alm;

       /// the normalization
       std::vector<double> m_normalization;

       /// vector for temporary computation of associated legendre polynomials
       std::vector<double> m_Plm;

       /// vector for temporary computation of spherical harmonics
       std::vector<std::complex<double>> m_sph;

     public:
       /**
	*  @name Constructors/destructors
	*/
       ///@{

       /**
	* @brief default constructor
	*
	* @return object of type SphericalHarmonics_Coefficients
	*/
       SphericalHarmonics_Coefficients () {}

       /**
	* @brief default constructor
	*
	* @param norder the number of multipoles, \f$ l_{max}+1 \f$
	*
	* @param nbins the number of separation bins
	*
	* @return object of type SphericalHarmonics_Coefficients
	*/
       SphericalHarmonics_Coefficients (const int norder, const int nbins=1) { initialize(norder, nbins); }

       /**
	* @brief default descructor
	*
	* @return none
	*/
       ~SphericalHarmonics_Coefficients () {}

       ///@}

       /**
	* @brief return the real part of the n-th coefficient
	* of the expansion for a given separation bin
	*
	* @param n the n-th coefficient of the spherical harmonics
	* expansion
	*
	* @param bin the separation bin
	*
	* @return  the real part of the n-th coefficient
	* of the expansion for a given separation bin
	*/
       double real (const int n, const int bin=0) { return m_alm[bin][n].real();} 

       /**
	* @brief return the imaginary part of the n-th coefficient
	* of the expansion for a given separation bin
	*
	* @param n the n-th coefficient of the spherical harmonics
	* expansion
	*
	* @param bin the separation bin
	*
	* @return  the imaginary part of the n-th coefficient
	* of the expansion for a given separation bin
	*/
       double imag (const int n, const int bin=0) { return m_alm[bin][n].imag();} 

       /**
	* @brief initialize the internal quantities
	*
	* @param norder the number of multipoles, \f$ l_{max}+1 \f$
	*
	* @param nbins the number of separation bins 
	*
	* @return none
	*/
       void initialize (const int norder, const int nbins=1);

       /**
	* @brief reset the internal quantities
	*
	* @return none
	*/
       void reset ();

       /**
	* @brief compute the \f$ a_{lm}\f$ for the normalized
	* coordinates \f$ \lbrace x, y, z \rbrace \f$
	*
	* @param xx the x coordinate
	*
	* @param yy the y coordinate
	*
	* @param zz the z coordinate
	*
	* @return vector containing the \f$ a_{lm}\f$ for the 
	* normalized coordinates \f$ \lbrace x, y, z \rbrace \f$
	*/
       std::vector<std::complex<double>> alm(const double xx, const double yy, const double zz);

       /**
	* @brief add the \f$ a_{lm}\f$ to a specific separation
	* bin with a weight
	*
	* @param alm the spherical harmonics expansion coefficients
	*
	* @param ww the weight
	*
	* @param bin the separation bin
	*
	* @return none
	*/
       void add (const std::vector<std::complex<double>> alm, const double ww, const int bin=0);

       /**
	* @brief compute add the \f$ a_{lm}\f$ for the normalized
	* coordinates \f$ \lbrace x, y, z \rbrace \f$ to a specific 
	* separation bin with a weight
	*
	* @param xx the x coordinate
	*
	* @param yy the y coordinate
	*
	* @param zz the z coordinate
	*
	* @param ww the weight
	*
	* @param bin the separation bin
	*
	* @return none
	*/
       void add (const double xx, const double yy, const double zz, const double ww, const int bin=0);

       /**
	* @brief compute the product of the 
	* \f$ a_{lm} \f$ in two separation bin.
	*
	* This function computes the product of the 
	* \f$ a_{lm} \f$ in two separation bin:
	*
	* \f[
	*   \zeta_l(r_1, r_2) = \sum_{m=-l}^l a_{lm}(r_1) a^*_{lm} (r_2)
	* \f]
	*
	* @param l the coefficient order
	*
	* @param bin1 the first separation bin
	*
	* @param bin2 the second separation bin
	*
	* @return none
	*/
       double power (const int l, const int bin1, const int bin2);
   };

 }
}

#endif
