/*******************************************************************
 *  Copyright (C) 2016 by Alfonso Veropalumbo                      *
 *  alfonso.veropalumbo@unibo.it                                   *
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
 *  @file Headers/FFTlog.h
 *
 *  @brief Wrapper for fftlog wripper
 *
 *  This file contains the prototypes of the class FFTlog,
 *  wrapping functions for fast Hankel transform 
 *  implemented in the FFTlog library (Hamilton 2000)
 *
 *  @author Alfonso Veropalumbo 
 *
 *  @author alfonso.veropalumbo@unbo.it
 */

#ifndef __FFTlog__
#define __FFTlog__

#include "FuncGrid.h"

namespace cbl {
  
  namespace wrapper {
  
    /**
     *  @brief The namespace of the <B> FFTlog wrappers </B>
     *  
     *  The \e fftlog namespace contains all the wrapper functions of
     *  the FFTlog routines, by Hamilton (see Hamilton 2000, Appendix B)
     *
     *  For a complete description of the routines please refer to
     *  http://casa.colorado.edu/~ajsh/FFTLog/
     */
    namespace fftlog {

      extern"C"
      {
	/**
	 *  @brief wrapper of the fhti subroutine contained in
	 *  External/fftlog-f90-master/fftlog.f This is an
	 *  initialization routine
	 *
	 *  @param [in] _n number of points in the array to be transformed
	 *  @param [in] _mu index of \f$J_\mu\f$ in Hankel transform
	 *  @param [in] _q exponent of power law bias
	 *  @param [in] _dlnr separation between natural log of points
	 *
	 *  @param [in] _kr \f$k_c\cdot r_c\f$ where c is the central
	 *  point of the array \f$kr = k_j r_{n+1-j} = k_{n+1-j} r_j\f$
	 *
	 *  @param [in] _kropt 0 &rarr; use input kr as is; 1 &rarr;
	 *  change kr to nearest low-ringing kr, quietly; 2 &rarr;
	 *  change kr to nearest low-ringing kr, verbosely; 3 &rarr;
	 *  change kr interactively
	 *
	 *  @param [out] _wsave working array
	 *  @param [out] _ok 1 &rarr; all went ok; 0 &rarr; error in initializations
	 *
	 *  @return none
	 */
	void fhti_ (int *_n, double *_mu, double *_q, double *_dlnr, double *_kr, int *_kropt, double *_wsave, int *_ok);

	/**
	 *  @brief wrapper of the fftl subroutine contained in
	 *  External/fftlog-f90-master/fftlog.f This subroutine computes
	 *  the discrete Fourier sine or cosine transform of a
	 *  logarithmically spaced periodic sequence
	 *
	 *  @param _n number of points in the array to be transformed
	 *
	 *  @param _a input &rarr; the array A to transform, output
	 *  &rarr; the transformed array
	 *
	 *  @param _rk \f$r_c/k_c\f$
	 *
	 *  @param _dir 1 &rarr; forward transform, -1 &rarr; backward
	 *  transform
	 *
	 *  @param _wsave working array
	 *
	 *  @return none
	 */
	void fftl_ (int *_n, double *_a, double *_rk, int *_dir, double *_wsave);
      }

      /**
       *  @brief wrapper of the FFTlog to compute  
       *  the discrete Fourier sine or cosine transform of a 
       *  logarithmically spaced periodic sequence
       *
       *  @param yy array containing the output positions
       *
       *  @param dir &rarr; forward transform, -1 &rarr; backward
       *  transform
       *
       *  @param xx the input position of the array to transform
       *  @param fx the array to transform
       *  @param mu index of \f$J_\mu\f$ in Hankel transform
       *  @param q exponent of power law bias
       *
       *  @param kr \f$k_c\cdot r_c\f$ where c is the central point of
       *  the array \f$kr = k_j r_{n+1-j} = k_{n+1-j} r_j\f$
       *
       *  @param kropt 0 &rarr; use input kr as is; 1 &rarr; change kr
       *  to nearest low-ringing kr, quietly; 2 &rarr; change kr to
       *  nearest low-ringing kr, verbosely; 3 &rarr; change kr
       *  interactively.
       *
       *  @return the transformed array
       */
      std::vector<double> transform_FFTlog (const std::vector<double> yy, const int dir, const std::vector<double> xx, const std::vector<double> fx, const double mu=0, const double q=0, const double kr=1, const int kropt=0);

      /**
       *  @brief wrapper of the FFTlog to compute the discrete Fourier
       *  sine or cosine transform of a logarithmically spaced periodic
       *  sequence
       *
       *  @param yy array containing the output positions
       *  @param fy the transformed array 
       *
       *  @param dir &rarr; forward transform, -1 &rarr; backward
       *  transform
       *
       *  @param xx the input position of the array to transform
       *  @param fx the array to transform
       *  @param mu index of \f$J_\mu\f$ in Hankel transform
       *  @param q exponent of power law bias
       *
       *  @param kr \f$k_c\cdot r_c\f$ where c is the central point of
       *  the array \f$kr = k_j r_{n+1-j} = k_{n+1-j} r_j\f$ .
       *
       *  @param kropt 0 &rarr; use input kr as is; 1 &rarr; change kr
       *  to nearest low-ringing kr, quietly; 2 &rarr; change kr to
       *  nearest low-ringing kr, verbosely; 3 &rarr; change kr
       *  interactively.
       *
       *  @return the transformed array
       */
      void transform_FFTlog (std::vector<double> &yy, std::vector<double> &fy, const int dir, const std::vector<double> xx, const std::vector<double> fx, const double mu=0, const double q=0, const double kr=1, const int kropt=0);

    }
  }

}

#endif
