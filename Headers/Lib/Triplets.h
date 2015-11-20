/*******************************************************************
 *  Copyright (C) 2015 by Michele Moresco, Federico Marulli        *
 *  and Alfonso Veropalumbo                                        *
 *                                                                 *
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
 *  @file Headers/Lib/Triplets.h
 *
 *  @brief The class Triplets
 *
 *  This file defines the interface of the class Triplets, used to handle
 *  triplets of objects to compute the two-point correlation function
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#ifndef __TRIPLETS__
#define __TRIPLETS__


#include "TwoPointCorrelation.h"


// ===================================================================================================


namespace cosmobl {
  
  /**
   *  @class Triplets Triplets.h "Headers/Lib/Triplets.h"
   *
   *  @brief The class Triplets
   *
   *  This class is used to handle objects of type <EM> Triplets
   *  </EM>. It contains all virtual methods implemented in the
   *  derived classes Triplets2D and Triplets3D
   */
  class Triplets {

  public:

    /**
     *  @brief default destructor
     *  @return none
     */
    virtual ~Triplets () {}
    
    /**
     *  @brief get the private member \e m_binsize
     *  @return the bin size, or an error message if the
     *  derived object does not have this member
     */
    virtual double binsize () const { cosmobl::ErrorMsg("Error in Triplets::binsize() of Triplets.h!"); return 0; }

    /**
     *  @brief get the private member \e m_side_s
     *  @return the size of r12, or an error message if the
     *  derived object does not have this member
     */
    virtual double side_s () const { cosmobl::ErrorMsg("Error in Triplets::side_s() of Triplets.h!"); return 0; }

    /**
     *  @brief get the private member \e m_side_u
     *  @return the ratio r13/r12, or an error message if the
     *  derived object does not have this member
     */
    virtual double side_u () const { cosmobl::ErrorMsg("Error in Triplets::side_u() of Triplets.h!"); return 0; }

    /**
     *  @brief get the private member \e m_perc_increase
     *  @return the ratio &Delta;r12/r12=&Delta;r13/r13, or an error
     *  message if the derived object does not have this member
     */
    virtual double perc_increase () const { cosmobl::ErrorMsg("Error in Triplets::perc_increase() of Triplets.h!"); return 0; }

    /**
     *  @brief get the private member \e m_TT[i]
     *  @param i the bin index
     *  @return the number of triplets in the i-th angular bin, or an
     *  error message if the derived object does not have this member
     */
    virtual double TT (const int i) const { cosmobl::ErrorMsg("Error in Triplets::TT() of Triplets.h!"); return 0; }

    /**
     *  @brief get the private member \e m_TT
     *  @return the vector containing the number of triplets in linear
     *  bins, or an error message if the derived object does not have
     *  this member
     */
    virtual vector<double> TT () const { cosmobl::ErrorMsg("Error in Triplets::TT() of Triplets.h!"); vector<double> TT; return TT;}
  
    /**
     *  @brief sum the number of triplets
     *  @param tt pointer to an object of class Triplets
     *  @param ww the weight
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void sum (const shared_ptr<Triplets> tt, const double ww=1.) { cosmobl::ErrorMsg("Error in Triplets::sum() of Triplets.h!"); }

    /**
     *  @brief estimate the distance between two objects and update
     *	the triplet vectors accordingly
     *  @param r12 distance between object1 and object2
     *  @param r13 distance between object1 and object3
     *  @param r23 distance between object2 and object3
     *  @param ww the weight
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void put (const double r12, const double r13, const double r23, const double ww=1.)
    { cosmobl::ErrorMsg("Error in Triplets::put() of Triplets.h!"); }

    /**
     *  @brief estimate the distance between three objects and update
     *	the triplet vectors accordingly
     *  @param obj1 pointer to an object of class Object
     *  @param obj2 pointer to an object of class Object
     *  @param obj3 pointer to an object of class Object
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2, const shared_ptr<Object> obj3)
    { cosmobl::ErrorMsg("Error in Triplets::put() of Triplets.h!"); }
    
  };


  // ============================================================================

  /**
   *  @class Triplets2D Triplets.h "Headers/Lib/Triplets.h"
   *
   *  @brief The class Triplets2D
   *
   *  This class is used to handle objects of type <EM> Triplets2D
   *  </EM>, used to measure the angular two-point correlation
   *  function
   */
  class Triplets2D : public Triplets
  {

  protected:

    /// the bin type
    string m_type_binning;

    /// the bin size
    double m_binsize;

    /// the size of r12
    double m_side_s;

    /// the ratio r13/r12
    double m_side_u; 

    /// the ratio &Delta;r12/r12=&Delta;r13/r13
    double m_perc_increase;

    /// the number of triplets;
    vector<double> m_TT;

  
  public:
  
    /**
     *  @brief default constructor
     *  @return object of class Triplets2D, with protected members set to -1
     */
    Triplets2D () 
      : m_type_binning(""), m_binsize(-1), m_side_s(-1), m_side_u(-1.), m_perc_increase(-1.) {}
    
    /**
     *  @brief constructor
     *  @param type_binning the type of binning
     *  @param binsize the bin size
     *  @param side_s the size of r12
     *  @param side_u the ratio r13/r12
     *  @param perc_increase the ratio &Delta;r12/r12=&Delta;r13/r13
     *  @return object of class Triplets2D
     */
    Triplets2D (string const type_binning, const double binsize, const double side_s, const double side_u, const double perc_increase) 
      : m_type_binning(type_binning), m_binsize(binsize), m_side_s(side_s), m_side_u(side_u), m_perc_increase(perc_increase)
    {
      if (type_binning=="ang")
	m_TT.resize(nint(par::pi/binsize)+1, 0.);
      else if (type_binning=="lin")
	m_TT.resize(nint((side_s*(1+side_u)-side_s*(side_u-1))/binsize)+1, 0.);
    }
  
    /**
     *  @brief get the protected member Triplets2D::m_binsize
     *  @return the bin size
     */
    double binsize () const override { return m_binsize; }

    /**
     *  @brief get the protected member Triplets2D::m_side_s
     *  @return the size of r12
     */
    double side_s () const override { return m_side_s; }

    /**
     *  @brief get the protected member Triplets2D::m_side_u
     *  @return the ratio r13/r12
     */
    double side_u () const override { return m_side_u; }

    /**
     *  @brief get the protected member Triplets2D::m_perc_increase
     *  @return the ratio &Delta;r12/r12=&Delta;r13/r13
     */    
    double perc_increase () const override { return m_perc_increase; }

    /**
     *  @brief get the private member Triplets2D::m_TT[i]
     *  @param i the bin index
     *  @return the number of triplets in the i-th bin
     */
    double TT (const int i) const override { return m_TT[i]; }

    /**
     *  @brief get the private member Triplets2D::m_TT
     *  @return the vector containing the number of triplets
     */
    vector<double> TT () const override { return m_TT; }
  
    /**
     *  @brief sum the number of triplets
     *  @param tt pointer to an object of class Triplets
     *  @param ww the weight
     *  @return none
     */
    void sum (const shared_ptr<Triplets>, const double ww=1) override;

    /**
     *  @brief estimate the distance between three objects and update
     *	the triplet vectors accordingly
     *  @param r12 distance between object1 and object2
     *  @param r13 distance between object1 and object3
     *  @param r23 distance between object2 and object3
     *  @param ww total weight
     *  @return none
     */
    void put (const double, const double, const double, const double ww=1.) override;
    
    /**
     *  @brief estimate the distance between three objects and update
     *	the triplet vectors accordingly
     *  @param obj1 pointer to an object of class Object
     *  @param obj2 pointer to an object of class Object
     *  @param obj3 pointer to an object of class Object
     *  @return none
     */
    void put (const shared_ptr<Object>, const shared_ptr<Object>, const shared_ptr<Object>) override;
  };


  // ============================================================================
  
  /**
   *  @class Triplets3D Triplets.h "Headers/Lib/Triplets.h"
   *
   *  @brief The class Triplets3D
   *
   *  This class is used to handle objects of type <EM> Triplets3D
   *  </EM>, used to measure the 3D two-point correlation function
   */
  class Triplets3D : public Triplets
  {
 
  protected:

    /// the bin type
    string m_type_binning;
    
    /// the bin size
    double m_binsize;

    /// the size of r12
    double m_side_s;

    /// the ratio r13/r12
    double m_side_u; 

    /// the ratio &Delta;r12/r12=&Delta;r13/r13
    double m_perc_increase;

    /// the number of triplets;
    vector<double> m_TT;

  
  public:
  
    /**
     *  @brief default constructor
     *  @return object of class Triplets3D, with protected members set to -1
     */
    Triplets3D () 
      : m_type_binning(""), m_binsize(-1), m_side_s(-1), m_side_u(-1.), m_perc_increase(-1.) {}
    
    /**
     *  @brief constructor
     *  @param type_binning the type of binning
     *  @param binsize the bin size
     *  @param side_s the size of r12
     *  @param side_u the ratio r13/r12
     *  @param perc_increase the ratio &Delta;r12/r12=&Delta;r13/r13
     *  @return object of class Triplets3D
     */
    Triplets3D (const string type_binning, const double binsize, const double side_s, const double side_u, const double perc_increase) 
      : m_type_binning(type_binning), m_binsize(binsize), m_side_s(side_s), m_side_u(side_u), m_perc_increase(perc_increase)
    {
      if (type_binning=="ang")
	m_TT.resize(nint(par::pi/binsize)+1, 0.);
      else if (type_binning=="lin")
	m_TT.resize(nint((side_s*(1+side_u)-side_s*(side_u-1))/binsize)+1, 0.);
    }

    /**
     *  @brief get the protected member Triplets3D::m_binsize
     *  @return the bin size
     */ 
    double binsize () const override { return m_binsize; }

    /**
     *  @brief get the protected member Triplets3D::m_side_s
     *  @return the size of r12
     */
    double side_s () const override { return m_side_s; }

    /**
     *  @brief get the protected member Triplets3D::m_side_u
     *  @return the ratio r13/r12
     */
    double side_u () const override { return m_side_u; }

    /**
     *  @brief get the protected member Triplets3D::m_perc_increase
     *  @return the ratio &Delta;r12/r12=&Delta;r13/r13
     */    
    double perc_increase () const override { return m_perc_increase; }

    /**
     *  @brief get the private member Triplets3D::m_TT[i]
     *  @param i the bin index
     *  @return the number of triplets in the i-th bin
     */
    double TT (const int i) const override { return m_TT[i]; }

    /**
     *  @brief get the private member Triplets3D::m_TT
     *  @return the vector containing the number of triplets
     */
    vector<double> TT () const override { return m_TT; }
  
    /**
     *  @brief sum the number of triplets
     *  @param tt pointer to an object of class Triplets
     *  @param ww the weight
     *  @return none
     */
    void sum (const shared_ptr<Triplets>, const double ww=1) override;

    /**
     *  @brief estimate the distance between three objects and update
     *	the triplet vectors accordingly
     *  @param r12 distance between object1 and object2
     *  @param r13 distance between object1 and object3
     *  @param r23 distance between object2 and object3
     *  @param ww total weight
     *  @return none
     */
    void put (const double, const double, const double, const double ww=1.) override;
    
    /**
     *  @brief estimate the distance between three objects and update
     *	the triplet vectors accordingly
     *  @param obj1 pointer to an object of class Object
     *  @param obj2 pointer to an object of class Object
     *  @param obj3 pointer to an object of class Object
     *  @return none
     */
    void put (const shared_ptr<Object>, const shared_ptr<Object>, const shared_ptr<Object>) override;
  };
}

#endif
