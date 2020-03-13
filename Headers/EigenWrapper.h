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
 *  @file Headers/EigenWrapper.h
 *
 *  @brief functions that wrap Eigen routines
 *
 *  This file contains the wrappers of Eigen routines for 
 *  vector and matrix manipulation
 *
 *  @author Alfonso Veropalumbo 
 *
 *  @author alfonso.veropalumbo@uniroma3.it
 */

#ifndef __EIGwrap__
#define __EIGwrap__

#include "Kernel.h"

namespace cbl {

  /**
   *  @brief The namespace of the <B> wrappers </B>
   *  
   *  The \e gsl namespace contains all the wrapper functions of the
   *  CBL
   */
  namespace wrapper {

    /**
     *  @brief The namespace of the <B> Eigen wrappers </B>
     *  
     *  The \e eigen namespace contains all the wrapper functions of the
     *  Eigen routines
     */
    namespace eigen {

      /**
       * @brief convert an Eigen::MatrixXd to
       * a std::vector<double>
       *
       * @param vec object of type Eigen::MatrixXd
       *
       * @return object of type std::vector<double>
       */
      std::vector<double> EigenToVector (const Eigen::MatrixXd vec);

      /**
       * @brief convert an Eigen::MatrixXd to
       * a std::vector<std::vector<double>>
       *
       * @param mat object of type Eigen::MatrixXd
       *
       * @return object of type std::vector<std::vector<double>>
       */
      std::vector<std::vector<double>> EigenToMatrix (const Eigen::MatrixXd mat);

      /**
       * @brief convert a std::vector<double> to
       * an Eigen::MatrixXd object
       *
       * @param vec object of type vector<double>
       *
       * @return object of type Eigen::MatrixXd
       */
      Eigen::MatrixXd VectorToEigen (const std::vector<double> vec);

      /**
       * @brief convert a std::vector<std::vector<double>> to
       * an Eigen::MatrixXd object
       *
       * @param mat object of type vector<double>
       *
       * @return object of type Eigen::MatrixXd
       */
      Eigen::MatrixXd MatrixToEigen (const std::vector<std::vector<double>> mat);

      /**
       * @brief convert a std::vector<double> to
       * an Eigen::MatrixXd object
       *
       * This function converts a std::vector<double> to an
       * Eigen::MatrixXd object. It assumes that the input vector is a
       * flat square matrix of order equal to the square root of the
       * vector size.
       *
       * @param mat object of type vector<double>
       *
       * @return object of type Eigen::MatrixXd
       */
      Eigen::MatrixXd SquareMatrixToEigen (const std::vector<double> mat);
    }
  }
}

#endif
