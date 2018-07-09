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
 *  @file Headers/Modelling.h
 *
 *  @brief The class Modelling
 *
 *  This file defines the interface of the class Modelling, used for
 *  modelling any kind of measurements
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLING1D__
#define __MODELLING1D__


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
     *  @class Modelling1D Modelling1D.h
     *  "Headers/Modelling1D.h"
     *
     *  @brief The class Modelling1D
     *
     *  This file defines the interface of the base class Modelling1D,
     *  used for modelling any kind of measurements
     *
     */
    class Modelling1D : public Modelling {

      public:

        /**
         *  @name Constructors/destructors
         */
        ///@{

        /**
         *  @brief default constuctor
         *  @return object of class Modelling
         */
        Modelling1D () {}

        /**
         *  @brief default destructor
         *  @return none
         */
        ~Modelling1D () = default;

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
         *  @return none
         */
        void set_fit_range (const double xmin, const double xmax);

        ///@}


        /**
         *  @brief write the model at xx
         *  for given parameters
         *
         *  @param output_dir the output directory
         *  @param output_file the output file
         *  @param xx vector of points at which the model is computed,
         *  @param parameters vector containing the input parameters
         *  used to compute the model; if this vector is not provided,
         *  the model will be computed using the best-fit parameters
         *
         *  @return none
         */
        void write_model (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> parameters);

        /**
         *  @brief write the model at xx 
         *  with best-fit parameters obtained from likelihood
         *  maximization
         *
         *  @param output_dir the output directory
         *  @param output_file the output file
         *  @param xx vector of points at which the model is computed,
         *
         *  @return none
         */
        void write_model_at_bestfit (const std::string output_dir, const std::string output_file, const std::vector<double> xx);

        /**
         *  @brief write the model at xx 
         *  computing 16th, 50th and 84th percentiles
         *  from the chains.
         *
         *  @param output_dir the output directory
         *  @param output_file the output file
         *  @param xx vector of points at which the model is computed,
         *  @param start the starting position for each chain
         *  @param thin the position step
         *
         *  @return none
         */
        void write_model_from_chains (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const int start=0, const int thin=1);

    };
  }
}

#endif
