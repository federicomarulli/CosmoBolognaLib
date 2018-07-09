/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/ModelFunction_TwoPointCorrelation_projected.h
 *
 *  @brief Functions to model the projected two-point correlation
 *  function
 *
 *  This file contains all the prototypes of the functions used to
 *  model the projected two-point correlation function
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODFUNCTWOPPRO__
#define __MODFUNCTWOPPRO__

#include "ModelFunction_TwoPointCorrelation1D_monopole.h"


// ============================================================================


namespace cbl {

  namespace modelling {

    namespace twopt {

      /**
       *  @name Functions of the Halo Occupation Distribution (HOD) models 
       */
      ///@{
      
      /**
       *  @brief function used to compute the projected two-point
       *  correlation function
       *
       *  this function is used to compute the projected correlation
       *  function by integrating the real-space spherically averaged
       *  two-point correlation function, as follows:
       *
       *  \f[w_{p}(r_p, z) = 2\int_{r_p}^{r_{out}}\xi(r,
       *  z)\frac{r\,{\rm d}r} {\sqrt{r^2-r_p^2}} \f]
       *
       *  where \f$r_{out}=\sqrt{r_p^2+r_{max}^2}\f$ and \f$\xi(r,
       *  z)\f$ is either the 1-halo, the 2-halo, or the full-shape
       *  real-space correlation function, computed by either
       *  cbl::modelling::twopt::xi_1halo,
       *  cbl::modelling::twopt::xi_2halo, or
       *  cbl::modelling::twopt::xi_HOD, respectively
       *
       *  in the limit \f$r_{max}\rightarrow\infty\f$, this function
       *  is completely independent of peculiar velocities
       *
       *  @param func the two-point real-space correlation function
       *  that will be integrated (it can be either the 1-halo, the
       *  2-halo, or the full-shape correlation function)
       *
       *  @param rp \f$r_p\f$: the scale perpendicular to the line of
       *  sight at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the projected two-point correlation function
       *
       *  @warning This function does not account for residual
       *  redshift-space distortions (RRSD) caused by the finite
       *  integration range used (i.e. \f$r_{max}<\infty\f$). If it is
       *  used to model the projected correlation function measured
       *  from real data, it introduces a bias casued by neglecting
       *  the RRSD (see e.g. sec. 2.3 of van den Bosch et
       *  al. 2012). It is thus useful only for testing these errors,
       *  or for general theoretical investigations not dealing with
       *  real data.
       */
      std::vector<double> wp_from_xi_approx (FunctionVectorVectorPtrVectorRef func, const std::vector<double> rp, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       *  @brief model for the 1-halo term of the projected two-point
       *  correlation function
       *
       *  this function estimates the 1-halo term of the projected
       *  correlation function by integrating the 1-halo term of the
       *  real-space spherically averaged two-point correlation
       *  function, as follows:
       *
       *  \f[w_{p, 1halo}(r_p, z) =
       *  2\int_{r_p}^{r_{out}}\xi_{1halo}(r, z)\frac{r\,{\rm d}r}
       *  {\sqrt{r^2-r_p^2}} \f]
       *
       *  where \f$r_{out}=\sqrt{r_p^2+r_{max}^2}\f$ and
       *  \f$\xi_{1halo}(r, z)\f$ is computed by
       *  cbl::modelling::twopt::xi_1halo
       *
       *  in the limit \f$r_{max}\rightarrow\infty\f$, this function
       *  is completely independent of peculiar velocities
       *
       *  @param rp \f$r_p\f$: the scale perpendicular to the line of
       *  sight at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the 1-halo term of the projected two-point
       *  correlation function
       *
       *  @warning This function does not account for residual
       *  redshift-space distortions (RRSD) caused by the finite
       *  integration range used (i.e. \f$r_{max}<\infty\f$). If it is
       *  used to model the projected correlation function measured
       *  from real data, it introduces a bias casued by neglecting
       *  the RRSD (see e.g. sec. 2.3 of van den Bosch et
       *  al. 2012). It is thus useful only for testing these errors,
       *  or for general theoretical investigations not dealing with
       *  real data.
       */
      std::vector<double> wp_1halo_approx (const std::vector<double> rp, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       *  @brief model for the 2-halo term of the projected two-point
       *  correlation function
       *
       *  this function estimates the 2-halo term of the projected
       *  correlation function by integrating the 2-halo term of the
       *  real-space spherically averaged two-point correlation
       *  function, as follows:
       *
       *  \f[w_{p, 2halo}(r_p, z) =
       *  2\int_{r_p}^{r_{out}}\xi_{2halo}(r, z)\frac{r\,{\rm d}r}
       *  {\sqrt{r^2-r_p^2}} \f]
       *
       *  where \f$r_{out}=\sqrt{r_p^2+r_{max}^2}\f$ and
       *  \f$\xi_{2halo}(r, z)\f$ is computed by
       *  cbl::modelling::twopt::xi_2halo
       *
       *  in the limit \f$r_{max}\rightarrow\infty\f$, this function
       *  is completely independent of peculiar velocities
       *
       *  @param rp \f$r_p\f$: the scale perpendicular to the line of
       *  sight at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the 2-halo term of the projected two-point
       *  correlation function
       *
       *  @warning This function does not account for residual
       *  redshift-space distortions (RRSD) caused by the finite
       *  integration range used (i.e. \f$r_{max}<\infty\f$). If it is
       *  used to model the projected correlation function measured
       *  from real data, it introduces a bias casued by neglecting
       *  the RRSD (see e.g. sec. 2.3 of van den Bosch et
       *  al. 2012). It is thus useful only for testing these errors,
       *  or for general theoretical investigations not dealing with
       *  real data.
       */
      std::vector<double> wp_2halo_approx (const std::vector<double> rp, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       *  @brief HOD model of the projected two-point correlation
       *  function
       *
       *  the function computes:
       *
       *  \f[w_p(r_p, z) = w_{p,1halo}(r_p, z)+w_{p,2halo}(r_p, z)\f]
       *
       *  where \f$w_{p,1halo}(r_p)\f$ and \f$w_{p,2halo}(r_p)\f$ are
       *  computed by cbl::modelling::twopt::wp_1halo_approx and
       *  cbl::modelling::twopt::wp_2halo_approx, respectively
       *
       *  in the limit \f$r_{max}\rightarrow\infty\f$, this function
       *  is completely independent of peculiar velocities
       *
       *  @param rp \f$r_p\f$: the scale perpendicular to the line of
       *  sight at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the HOD projected two-point correlation function
       *
       *  @warning This function does not account for residual
       *  redshift-space distortions (RRSD) caused by the finite
       *  integration range used (i.e. \f$r_{max}<\infty\f$). If it is
       *  used to model the projected correlation function measured
       *  from real data, it introduces a bias casued by neglecting
       *  the RRSD (see e.g. sec. 2.3 of van den Bosch et
       *  al. 2012). It is thus useful only for testing these errors,
       *  or for general theoretical investigations not dealing with
       *  real data.
       */
      std::vector<double> wp_HOD_approx (const std::vector<double> rp, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       *  @brief function used to compute the projected two-point
       *  correlation function
       *
       *  this function is used to compute the projected
       *  correlation function by integrating the redshift-space 2D
       *  correlation function (e.g. van den Bosch et al. 2012):
       *
       *  \f[w_{p}(r_p, z) = 2\int_{0}^{\pi_{max}}\xi(r_p, \pi,
       *  z)\,{\rm d}\pi\f]
       *
       *  where the redshift-space galaxy correlation function,
       *  \f$\xi(r_p, \pi, z)\f$, is either the 1-halo, the 2-halo, or
       *  the full-shape redshift-space correlation function, computed
       *  by either cbl::modelling::twopt::xi_1halo_zspace,
       *  cbl::modelling::twopt::xi_2halo_zspace, or
       *  cbl::modelling::twopt::xi_HOD_space, respectively
       *
       *  @param func the redshift-space two-point correlation
       *  function that will be integrated (it can be either the
       *  1-halo, the 2-halo, or the full-shape correlation function)
       *
       *  @param rp \f$r_p\f$: the scale perpendicular to the line of
       *  sight at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the projected two-point correlation function
       *
       *  @warning This function accounts for residual redshift-space
       *  distortions (RRSD) caused by the finite integration range
       *  used with real data (see e.g. sec. 2.3 of van den Bosch et
       *  al. 2012). Its accuracy depends on the accuracy of the
       *  redshift-space distortion model used to describe
       *  \f$\xi_g(r_p, \pi, z)\f$.
       */
      std::vector<double> wp_from_xi (FunctionDoubleDoubleDoublePtrVectorRef func, const std::vector<double> rp, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       *  @brief model for the 1-halo term of the projected two-point
       *  correlation function
       *
       *  this function computes the 1-halo term of the projected
       *  correlation function by integrating the redshift-space 2D
       *  correlation function (e.g. van den Bosch et al. 2012):
       *
       *  \f[w_{p, 1halo}(r_p, z) =
       *  2\int_{0}^{\pi_{max}}\xi_{1halo}(r_p, \pi, z)\,{\rm d}\pi\f]
       *
       *  where the 1-halo term of the redshift-space galaxy
       *  correlation function, \f$\xi_{1halo}(r_p, \pi, z)\f$, is
       *  computed by cbl::modelling::twopt::xi_1halo_zspace
       *
       *  @param rp \f$r_p\f$: the scale perpendicular to the line of
       *  sight at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the 1-halo term of the projected two-point
       *  correlation function
       *
       *  @warning This function accounts for residual redshift-space
       *  distortions (RRSD) caused by the finite integration range
       *  used with real data (see e.g. sec. 2.3 of van den Bosch et
       *  al. 2012). Its accuracy depends on the accuracy of the
       *  redshift-space distortion model used to describe
       *  \f$\xi_g(r_p, \pi, z)\f$.
       */
      std::vector<double> wp_1halo (const std::vector<double> rp, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       *  @brief model for the 2-halo term of the projected two-point
       *  correlation function
       *
       *  this function computes the 2-halo term of the projected
       *  correlation function by integrating the redshift-space 2D
       *  correlation function (e.g. van den Bosch et al. 2012):
       *
       *  \f[w_{p, 2halo}(r_p, z) =
       *  2\int_{0}^{\pi_{max}}\xi_{2halo}(r_p, \pi, z)\,{\rm d}\pi\f]
       *
       *  where the 2-halo term of the redshift-space galaxy
       *  correlation function, \f$\xi_{2halo}(r_p, \pi, z)\f$, is
       *  computed by cbl::modelling::twopt::xi_2halo_zspace
       *
       *  @param rp \f$r_p\f$: the scale perpendicular to the line of
       *  sight at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the 2-halo term of the projected two-point
       *  correlation function
       *
       *  @warning This function accounts for residual redshift-space
       *  distortions (RRSD) caused by the finite integration range
       *  used with real data (see e.g. sec. 2.3 of van den Bosch et
       *  al. 2012). Its accuracy depends on the accuracy of the
       *  redshift-space distortion model used to describe
       *  \f$\xi_g(r_p, \pi, z)\f$.
       */
      std::vector<double> wp_2halo (const std::vector<double> rp, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       *  @brief HOD model of the projected two-point
       *  correlation function
       *
       *  this function computes:
       *
       *  \f[w_p(r_p, z) = w_{p,1halo}(r_p, z)+w_{p,2halo}(r_p, z)\f]
       *
       *  where \f$w_{p,1halo}(r_p)\f$ and \f$w_{p,2halo}(r_p)\f$ are
       *  computed by cbl::modelling::twopt::wp_1halo and
       *  cbl::modelling::twopt::wp_2halo, respectively
       *
       *  @param rp \f$r_p\f$: the scale perpendicular to the line of
       *  sight at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  fixed input data used to construct the model
       *
       *  @param parameter vector containing the model parameters
       *
       *  @return the HOD projected two-point correlation function
       *
       *  @warning This function accounts for residual redshift-space
       *  distortions (RRSD) caused by the finite integration range
       *  used with real data (see e.g. sec. 2.3 of van den Bosch et
       *  al. 2012). Its accuracy depends on the accuracy of the
       *  redshift-space distortion model used to describe
       *  \f$\xi_g(r_p, \pi, z)\f$.
       */
      std::vector<double> wp_HOD (const std::vector<double> rp, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      ///@}
      
    }
  }
}

#endif
