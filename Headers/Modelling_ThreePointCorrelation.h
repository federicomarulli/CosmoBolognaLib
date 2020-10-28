/********************************************************************
 *  Copyright (C) 2017 by Federico Marulli                          *
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
 *  @file Headers/Modelling_ThreePointCorrelation.h
 *
 *  @brief The class Modelling_ThreePointCorrelation
 *
 *  This file defines the interface of the class
 *  Modelling_ThreePointCorrelation, used to model three-point
 *  correlation functions of any kind
 *
 *  @author Federico Marulli, Michele Moresco
 *
 *  @author federico.marulli3@unibo.it, michele.moresco@unibo.it
 */

#ifndef __MODELLINGTHREEPOINT__
#define __MODELLINGTHREEPOINT__


#include "ThreePointCorrelation.h"
#include "Modelling1D.h"
#include "ModelFunction_ThreePointCorrelation.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    /**
     *  @brief The namespace of the <B> three-point correlation
     *  function modelling </B>
     *  
     *  The \e modelling::threept namespace contains all the functions
     *  and classes to model the three-point correlation function
     */
    namespace threept {
    
      /**
       *  @class Modelling_ThreePointCorrelation
       *  Modelling_ThreePointCorrelation.h
       *  "Headers/Modelling_ThreePointCorrelation.h"
       *
       *  @brief The class Modelling_ThreePointCorrelation
       *
       *  This file defines the interface of the base class
       *  Modelling_ThreePointCorrelation, used for modelling any kind
       *  of three-point correlation function measurements
       *
       */
      class Modelling_ThreePointCorrelation : public Modelling1D
      {
      
      protected:
	
	/// the three-point correlation function type
	measure::threept::ThreePType m_threePType;

	/// the container of parameters for three-point correlation function model computation
	modelling::threept::STR_data_model_threept m_data_model;

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  ThreePointCorrelation
	 */
	Modelling_ThreePointCorrelation () = default;
	
	/**
	 *  @brief constuctor
	 *  @param threep the three-point correlation function to model
	 *  _ThreePointCorrelation
	 */
	Modelling_ThreePointCorrelation (const std::shared_ptr<cbl::measure::threept::ThreePointCorrelation> threep)
	  { m_data = threep->dataset(); }
	
	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~Modelling_ThreePointCorrelation () = default;

	/**
	 *  @brief static factory used to construct modelling of 
	 *  three-point correlation functions of any type
	 *
	 *  @param threep the three-point correlation function to
	 *  model
	 *
	 *  @return a pointer to an object of class 
	 *  Modelling_ThreePointCorrelation of a given type
	 */
	static std::shared_ptr<Modelling_ThreePointCorrelation> Create (const std::shared_ptr<measure::threept::ThreePointCorrelation> threep);

	/**
	 *  @brief static factory used to construct modelling of 
	 *  three-point correlation functions of any type
	 *
	 *  @param threePType type of the three-point correlation
	 *  function
	 *
	 *  @param threept_dataset the dataset containing the
	 *  three-point correlation function to model
	 *
	 *  @return a pointer to an object of class 
	 *  Modelling_ThreePointCorrelation of a given type
	 */
	static std::shared_ptr<Modelling_ThreePointCorrelation> Create (const measure::threept::ThreePType threePType, const std::shared_ptr<data::Data> threept_dataset);

	///@}

	
	/**
	 *  @brief return the type of correlation function
	 *  @return the type of correlation function
	 */
	measure::threept::ThreePType threePType () { return m_threePType; }


	// ============================================================================================

	
	/**
	 *  @brief set the data model for the three-point correlation
	 *  function 
	 *
	 *  @param Q_DM vector contaning the DM reduced three-point
	 *  correlation function
	 *
	 *  
	 */
	void set_data_model (const std::vector<double> Q_DM);

	/**
	 *  @brief set the data model for the three-point correlation
	 *  function (see Slepian, Eisenstein 2017)
         *
	 *  @param r1 the first triangle side
	 *  
	 *  @param r2 the second triangle side
	 *
	 *  @param cosmology the fiducial cosmology
	 *
	 *  @param redshift the redshift
	 *
	 *  @param method_Pk method used to compute the power
	 *  spectrum; valid choices for method_Pk are: CAMB
	 *  [http://camb.info/], CLASS [http://class-code.net/],
	 *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465],
	 *  EisensteinHu [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *
	 *  @param NL 0 &rarr; linear power spectrum; 1 &rarr;
	 *  non-linear power spectrum
	 *
	 *  @param max_ll the maximum legendre multipole in 3pcf model
         *    
	 *  @param k_min minimum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *  
	 *  @param k_max maximum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *
	 *  @param step_k number of steps used to compute the binned
	 *  dark matter correlation function
         *
	 *  @param r_min minimum separation up to which the integrals of the
         *  tree-level 3pcf are computed
	 *
	 *  @param r_max maximum separation up to which the integrals of the
         *  tree-level 3pcf are computed
         *
         *  @param step_r number of steps used to compute the integrals of the
         *  tree-level 3pcf
         *
         *  @param force_realSpace  \f$ \rightarrow \f$ redshift-space model; 1
	 *  \f$ \rightarrow \f$ real-space model
	 *
	 *  @param use_k false \f$ \rightarrow \f$ do not compute k-integrals;
	 *  true \f$ \rightarrow \f$ compute k-integrals
         *
         *  @param use_k  &rarr; do not compute the k-integrals; 1
	 *  &rarr; compute the k-integrals
	 *
	 *  @param output_dir the output_directory
	 *
	 *  @param store_output if true the output files created by
	 *  the Boltmann solver are stored; if false the output files
	 *  are removed
	 *
	 *  @param output_root the output file root
	 *
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *
	 *  @param prec accuracy of the GSL integration
         *
	 *  
	 */
	void set_data_model_zeta_RSD (const double r1, const double r2, const cbl::cosmology::Cosmology cosmology, const double redshift, const std::string method_Pk="CAMB", const bool NL=false, const int max_ll=5, const double k_min=1.e-4, const double k_max=100, const int step_k=500, const double r_min=1.e-4, const double r_max=200, const int step_r=200, const bool force_realSpace=false, const bool use_k=false, const std::string output_dir=cbl::par::defaultString, const bool store_output=true, const std::string output_root=cbl::par::defaultString, const int norm=-1, const double prec=1.e-4);

	/**
	 *  @brief set the data model for the three-point correlation
	 *  function with non-local contributions
	 *
	 *  @param cosmology the fiducial cosmology
         *
	 *  @param r1 the first triangle side
	 *
	 *  @param r2 the second triangle side
	 *
	 *  @param theta vector of theta
	 *
	 *  @param model method used to compute the 3pcf
	 *
	 *  @param kk vector containing wavevector moduls
	 *
	 *  @param Pk_DM vector containing the theoretical dark matter 
         *  power-spectrum
	 *
	 *  
	 */
	void set_data_Q_nonlocal (const cosmology::Cosmology cosmology, const double r1, const double r2, const std::vector<double> theta, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_DM);

      };
    }
  }
}

#endif
