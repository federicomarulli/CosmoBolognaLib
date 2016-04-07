/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo  *
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
 *******************************************************************/

/**
 *  @file Headers/Lib/LogNormal.h
 *
 *  @brief Implementation of the lognormal data structure
 *
 *  This file defines the interface of the class LogNormal
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __LOGNORMAL__
#define __LOGNORMAL__

#include "Catalogue.h"
#include <fftw3.h>
#include <erf.h>


namespace cosmobl {

  /**
   *  @class LogNormal LogNormal.h "Headers/Lib/LogNormal.h"
   *
   *  @brief The class LogNormal
   *
   *  This class is used to handle objects of type <EM> lognormal
   *  </EM>
   */
  class LogNormal {
    
  protected:

    /// the number of lognormal realization to produce/read
    int m_nLN;
    
    /// pointer to the input data 
    shared_ptr<catalogue::Catalogue> m_data;
    
    /// pointer to the random catalogues
    shared_ptr<catalogue::Catalogue> m_random;
    
    /// cector containing pointers to the LogNormal realizations
    vector<shared_ptr<catalogue::Catalogue> > m_LNCat;  

    /// 0 &rarr; the input is a cosmology; 1 &rarr; the input is &xi;(r)
    bool m_withxi;

    /// vector containing the binnend separations of the two-point corrleation function used to create the density field
    vector<double> m_rmodel;
    
    /// vector containing the starting two-point correlation function used to create the density field
    vector<double> m_ximodel;

    /// approximate cell size of the density field
    double m_rmin;
    
    /// bias of the lognormal density field to be realized
    double m_bias;
    
    /// pointer to the fiducial cosmology
    shared_ptr<Cosmology> m_cosmology;
    
    /// 0 &rarr; redshift-space (only monopole distortion); 1 &rarr; real-space
    bool m_Real;
    
    /// string containing the author for the model power spectrum
    string m_author;
    
    /// 0 &rarr; compute the linear power spectrum; 1 &rarr; compute the non-linear power spectrum
    bool m_NL;

    /// the cosmological model used to compute distances
    string m_model;

    
  public:

    /**
     *  @brief default constructor
     *  @return object of class LogNormal
     */
    LogNormal() {};

    /**
     *  @brief default destructor
     *  @return none
     */
    ~LogNormal() {};

    /**
     *  @brief constructor 
     *  @param data input data catalogue
     *
     *  @param random input random catalogue (should be much larger
     *  than the random catalogue used to measure &xi;(r))
     *
     *  @param nLN number of lognormal realizations
     *  @return object of class LogNormal
     */
    LogNormal (const shared_ptr<catalogue::Catalogue> data, const shared_ptr<catalogue::Catalogue> random, const int nLN) : m_nLN(nLN), m_data(data), m_random(random)
    {
      m_LNCat.resize(m_nLN);
    }

    /**
     *  @brief set data and random catalogues
     *  @param data input data catalogue
     *  @param random input random catalogue
     *  @return none
     */
    void setCatalogues (const shared_ptr<catalogue::Catalogue> data, const shared_ptr<catalogue::Catalogue> random);
    
    /**
     *  @brief set the starting two-point correlation function
     *  @param rr binned comoving separation 
     *  @param xi input two-point correlation function
     *  @return none
     */
    void setParameters_from_xi (const vector<double>, const vector<double>);

    /**
     *  @brief set the parameters to compute a prediction of &xi;(r)
     *  @param cosmology the input cosmology
     *  @param bias the bias parameter
     *  @param Real 0 &rarr; redshift-space (only monopole distortion); 1 &rarr; real-space
     *  @param author the method used to compute dark matter power spectrum
     *  @param NL 0 &rarr; compute the linear power spectrum; 1 &rarr; compute the non-linear power spectrum
     *  @param model the cosmological model used to compute distances
     *  @return none
     */
    void setParameters_from_model (const shared_ptr<Cosmology>, const double, const bool Real=1, const string author="CAMB", const bool NL=0, const string model="LCDM");
    
    /**
     *  @brief set the total number of realizations
     *  @param nLN the number of realizations
     *  @return none
     */
    void set_nLN (const int);
    
    /**
     *  @brief get the private member LogNormal::m_nLN
     *  @return the number of LogNormal realizations
     */
    int nLN () { return m_nLN; }

    /**
     *  @brief get the private member LogNormal::m_LNCat[i]
     *  @param i index of the LogNormal realization
     *  @return the i-th LogNormal realization
     */
    shared_ptr<catalogue::Catalogue> LNCat (const int i) { return m_LNCat[i]; }
    
    /**
     *  @brief generate the LogNormal mock catalogues
     *  @param rmin the cell size in comoving coordinates
     *  @param dir the output directory
     *  @param start the starting index of the mock to be created
     *  @param filename the prefix of the ouput file containing the
     *  LogNormal realizations
     *  @return none
     */
    void generate_LogNormal_mock (const double, const string, const int start=0, const string filename="lognormal_");
    
  };
}

#endif
