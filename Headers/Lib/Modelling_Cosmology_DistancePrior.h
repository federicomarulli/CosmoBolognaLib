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
 *  @file Headers/Lib/Modelling_Cosmology.h
 *
 *  @brief The class Modelling_Cosmology
 *
 *  This file defines the interface of the class Modelling, used for
 *  modelling cosmological measurements
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGCOSMDP__
#define __MODELLINGCOSMDP__


#include "Cosmology.h"
#include "Data1D.h"

// ===================================================================================================


namespace cosmobl {

  namespace modelling {

    namespace cosmology{

      class CMB_DistancePrior
      {

	protected:

	  /// dataset
	  shared_ptr<data::Data> m_dataset;

	public:

	  CMB_DistancePrior() {}

	  ~CMB_DistancePrior() {}

	  static shared_ptr<CMB_DistancePrior> Create (const string distance_prior_name);

	  virtual shared_ptr<data::Data> dataset()
	  {ErrorCBL("Error in dataset() of CMB_DistancePrior. No dataset for base class!"); return NULL;}
	  
	  virtual vector<double> model (const cosmobl::cosmology::Cosmology cosmology)
	  {(void)cosmology; ErrorCBL("Error in model() of CMB_DistancePrior. No model for base class!"); vector<double> vv; return vv;}

      };

      class Aubourg15_Planck15 : public CMB_DistancePrior
      {

	public:

	  Aubourg15_Planck15 () {
	    vector<double> redshift = {0., 0., 1090.};
	    vector<double> data = {0.02245, 0.1386, 94.33};
	    vector<vector<double>> covariance(3, vector<double>(3, 0));
	    covariance[0][0] =  1.286e-7;
	    covariance[0][1] =  -6.033e-7;
	    covariance[0][2] =  -1.443e-5;

	    covariance[1][0] =  -6.033e-7;
	    covariance[1][1] =  7.542e-6;
	    covariance[1][2] =  -3.605e-5;

	    covariance[2][0] =  1.443e-5;
	    covariance[2][1] =  -3.605e-5;
	    covariance[2][2] =  0.004264;

	    m_dataset = make_shared<data::Data1D>(data::Data1D(redshift, data, covariance));
	  }

	  shared_ptr<data::Data> dataset() {return m_dataset;}

	  vector<double> model (const cosmobl::cosmology::Cosmology cosmology)
	  {
	    vector<double> mm(3, 0);

	    mm[0] = cosmology.Omega_baryon()*cosmology.hh()*cosmology.hh();
	    mm[1] = cosmology.Omega_matter()*cosmology.hh()*cosmology.hh();
	    mm[2] = cosmology.D_M(m_dataset->xx(2))/cosmology.rs_CAMB();

	    return mm;
	  }

      };

      class Aubourg15_WMAP09 : public CMB_DistancePrior
      {

	public:

	  Aubourg15_WMAP09 () {

	    vector<double> redshift = {0., 0., 1090.};
	    vector<double> data = {0.02259, 0.1354, 94.51};
	    vector<vector<double>> covariance(3, vector<double>(3, 0));
	    covariance[0][0] =  2.864-7;
	    covariance[0][1] =  -4.809e-7;
	    covariance[0][2] =  -1.111-5;

	    covariance[1][0] = -4.809e-7;
	    covariance[1][1] = 1.908e-5;
	    covariance[1][2] = -7.495e-6;

	    covariance[2][0] = -1.111e-5;
	    covariance[2][1] = -7.495e-6;
	    covariance[2][2] = 0.02542;

	    m_dataset = make_shared<data::Data1D>(data::Data1D(redshift, data, covariance));
	  }

	  shared_ptr<data::Data> dataset() {return m_dataset;}

	  vector<double> model (const cosmobl::cosmology::Cosmology cosmology)
	  {
	    vector<double> mm(3, 0);

	    mm[0] = cosmology.Omega_baryon()*cosmology.hh()*cosmology.hh();
	    mm[1] = cosmology.Omega_matter()*cosmology.hh()*cosmology.hh();
	    mm[2] = cosmology.D_M(m_dataset->xx(2))/cosmology.rs_CAMB();

	    return mm;
	  }
      };
    }
  }
}

#endif
