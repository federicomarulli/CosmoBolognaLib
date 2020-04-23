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
 *  @file Headers/Modelling2D.h
 *
 *  @brief The class Modelling2D
 *
 *  This file defines the interface of the class Modelling2D, used for
 *  modelling any kind of measurements
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLING2D__
#define __MODELLING2D__


#include "Modelling.h"


// ===================================================================================================


namespace cbl {

  /**
   *  @brief The namespace of the functions and classes used for <B>
   *  modelling </B>
   *  
   * The \e modelling namespace contains all the functions and classes
   * used to model any kind of measurements
   */
  namespace modelling {

    /**
     *  @class Modelling2D Modelling2D.h
     *  "Headers/Modelling2D.h"
     *
     *  @brief The class Modelling2D
     *
     *  This file defines the interface of the base class Modelling2D,
     *  used for modelling any kind of measurements
     *
     */
    class Modelling2D : public Modelling {

      public:

        /**
         *  @name Constructors/destructors
         */
        ///@{

        /**
         *  @brief default constuctor
         *  @return object of class Modelling
         */
        Modelling2D () {}

        /**
         *  @brief default destructor
         *  @return none
         */
        ~Modelling2D () = default;

        ///@}

        /**
         *  @name Member functions used to set internal parameters
         */
        ///@{

        /**
         *  @brief set the fit range 
         *
         *  @param xmin minimum x value used for the fit
         *
         *  @param xmax maximum x value used for the fit
         *
         *  @param ymin minimum y value used for the fit
         *
         *  @param ymax maximum y value used for the fit
         *
         *  @return none
         */
        void set_fit_range (const double xmin, const double xmax, const double ymin, const double ymax);

        ///@}

        /**
         *  @brief write the model at xx, yy
         *  for given parameters
         *
         *  @param output_dir the output directory
         *  @param output_file the output file
         *  @param xx vector of points at which the model is computed,
         *  first axis
         *  @param yy vector of points at which the model is computed,
         *  second axis
         *  @param parameters vector containing the input parameters
         *  used to compute the model; if this vector is not provided,
         *  the model will be computed using the best-fit parameters
         *
         *  @return none
         */
        void write_model (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> yy, const std::vector<double> parameters);

        /**
         *  @brief write the model at xx, yy with best-fit parameters
         *  obtained from posterior maximization
         *  
         *  @param output_dir the output directory
         *  @param output_file the output file
         *  @param xx vector of points at which the model is computed,
         *  first axis
         *  @param yy vector of points at which the model is computed,
         *  second axis
         *
         *  @return none
         */
        void write_model_at_bestfit (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> yy);

        /**
         *  @brief write the model at xx, yy
         *  computing 16th, 50th and 84th percentiles
         *  from the chains.
         *
         *  @param output_dir the output directory
         *  @param output_file the output file
         *  @param xx vector of points at which the model is computed,
         *  first axis
         *  @param yy vector of points at which the model is computed,
         *  second axis
         *  @param start the starting position for each chain
         *  @param thin the position step
         *
         *  @return none
         */
        void write_model_from_chains (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> yy, const int start=0, const int thin=1);

        ///@}

  };
}
}

#endif
