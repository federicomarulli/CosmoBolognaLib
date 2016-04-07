/********************************************************************
 *  Copyright (C) 2015 by Michele Moresco, Federico Marulli         *
 *  and Alfonso Veropalumbo                                         *
 *                                                                  *
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
 *  @file Headers/Lib/Triplet.h
 *
 *  @brief The class Triplet
 *
 *  This file defines the interface of the class Triplet, used to handle
 *  triplets of objects to compute the three-point correlation function
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#ifndef __TRIPLET__
#define __TRIPLET__


#include "Catalogue.h"


// ===================================================================================================


namespace cosmobl {
  
  /**
   *  @brief The namespace of the triplets 
   *  
   *  The \e triplets namespace contains all the functions and classes
   *  to handle triplets of objects
   */
  namespace triplets {
    
    /**
     * @enum TripletType
     * @brief the triplet type
     */
    enum TripletType {

      /// 1D triplet in comoving coordinates and angular bins 
      _comoving_theta_,

      /// 1D triplet in comoving coordinates and linear bins 
      _comoving_side_

    };
  
  
    /**
     *  @class Triplet Triplet.h "Headers/Lib/Triplet.h"
     *
     *  @brief The class Triplet
     *
     *  This class is used to handle objects of type <EM> Triplet
     *  </EM>. It contains all virtual methods implemented in the
     *  derived classes Triplet2D and Triplet3D
     */
    class Triplet {

    private:
      
      /**
       *  @name Member functions used to set the binning parameters (customized in all the derived classes) 
       */
      ///@{
  
      /**
       *  @brief set the binning parameters
       *  @return none
       */
      virtual void set_parameters () = 0;
  
      ///@}

      
    protected:
      
      /// the dimension of the triplet vectors 
      Dim m_tripletDim;
     
      /// triplet type
      TripletType m_tripletType;
      
      
    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Triplet () {}

      /**
       *  @brief static factory used to construct triplets of any type
       *  
       *  @param type the triplet type; it can be: _comoving_theta,
       *  _comoving_side
       *  @param side_s the size of r<SUB>12</SUB>
       *  @param side_u the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
       *  @param perc_increase the ratio
       *  &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
       *  @param nbins the number of bins
       *  @return a pointer to an object of class Triplet of a given type
       */
      static shared_ptr<Triplet> Create (const TripletType type, const double side_s, const int side_u, const double perc_increase, const int nbins);

      ///@}
    

      /**
       *  @name Member functions used to get private/protected parameters
       */
      ///@{
      
      /**
       *  @brief get the dimension of the triplet vectors 
       *  @return the dimension of the triplet vectors 
       */
      Dim tripletDim () const { return m_tripletDim; }

      /**
       *  @brief get the triplet type
       *  @return the triplet type
       */
      TripletType tripletType () const { return m_tripletType; }

      /**
       *  @brief get the protected member \e m_scale[i]
       *  @param i the bin index
       *  @return the i-th binned scale
       */
      virtual double scale (const int i) const { cosmobl::ErrorMsg("Error in Triplet::scale() of Triplet.h!"); return 0; }

      /**
       *  @brief get the protected member \e m_scale
       *  @return the vector containing the binned scales
       */
      virtual vector<double> scale () const { cosmobl::ErrorMsg("Error in Triplet::scale() of Triplet.h!"); vector<double> vv; return vv; }
      
      /**
       *  @brief get the private member \e m_TT1D[i]
       *  @param i the bin index
       *  @return the number of triplets in the i-th angular bin, or an
       *  error message if the derived object does not have this member
       */
      virtual double TT1D (const int i) const { cosmobl::ErrorMsg("Error in Triplet::TT1D() of Triplet.h!"); return 0; }
      
      /**
       *  @brief get the private member \e m_TT1D
       *  @return the vector containing the number of triplets in linear
       *  bins, or an error message if the derived object does not have
       *  this member
       */
      virtual vector<double> TT1D () const { cosmobl::ErrorMsg("Error in Triplet::TT1D() of Triplet.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the protected member \e m_scale_D1[i]
       *  @param i the bin index in the first dimension
       *  @return the i-th binned scale
       */
      virtual double scale_D1 (const int i) const { cosmobl::ErrorMsg("Error in Triplet::scale_D1() of Triplet.h!"); return 0; }

      /**
       *  @brief get the protected member \e m_scale_D1
       *  @return the vector containing the binned scales
       */
      virtual vector<double> scale_D1 () const { cosmobl::ErrorMsg("Error in Triplet::scale_D1() of Triplet.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the protected member \e m_scale_D2[i]
       *  @param i the bin index in the first dimension
       *  @return the i-th binned scale
       */
      virtual double scale_D2 (const int i) const { cosmobl::ErrorMsg("Error in Triplet::scale_D2() of Triplet.h!"); return 0; }

      /**
       *  @brief get the protected member \e m_scale_D2
       *  @return the vector containing the binned scales
       */
      virtual vector<double> scale_D2 () const { cosmobl::ErrorMsg("Error in Triplet::scale_D2() of Triplet.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the protected member \e m_TT2D[i]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the number of pairs in the i-th bin
       */
      virtual double TT2D (const int i, const int j) const { cosmobl::ErrorMsg("Error in Triplet::TT2D() of Triplet.h!"); return 0; }

      /**
       *  @brief get the protected member \e m_TT2D
       *  @return the vector containing the binned number of pairs
       */
      virtual vector<vector<double> > TT2D () const { cosmobl::ErrorMsg("Error in Triplet::TT2D() of Triplet.h!"); vector<vector<double> > vv; return vv; }
      
      /**
       *  @brief get the private member \e m_side_s
       *  @return the size of r<SUB>12</SUB>, or an error message if the
       *  derived object does not have this member
       */
      virtual double side_s () const { cosmobl::ErrorMsg("Error in Triplet::side_s() of Triplet.h!"); return 0; }

      /**
       *  @brief get the private member \e m_side_u
       *  @return the ratio r<SUB>13</SUB>/r<SUB>12</SUB>, or an error
       *  message if the derived object does not have this member
       */
      virtual double side_u () const { cosmobl::ErrorMsg("Error in Triplet::side_u() of Triplet.h!"); return 0; }

      /**
       *  @brief get the private member \e m_perc_increase
       *  @return the ratio
       *  &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual double perc_increase () const { cosmobl::ErrorMsg("Error in Triplet::perc_increase() of Triplet.h!"); return 0; }

      /**
       *  @brief get the private member \e m_nbins
       *  @return the number of bins, or an error message if the
       *  derived object does not have this member
       */
      virtual int nbins () const { cosmobl::ErrorMsg("Error in Triplet::nbins() of Triplet.h!"); return 0; }

      /**
       *  @brief get the private member \e m_binSize
       *  @return the bin size, or an error message if the
       *  derived object does not have this member
       */
      virtual double binSize () const { cosmobl::ErrorMsg("Error in Triplet::binSize() of Triplet.h!"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_side_s_D1
       *  @return the size of r<SUB>12</SUB> in the first dimension,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual double side_s_D1 () const { cosmobl::ErrorMsg("Error in Triplet::m_side_s_D1() of Triplet.h!"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_side_u_D1
       *  @return the ratio r<SUB>13</SUB>/r<SUB>12</SUB> in the first
       *  dimension, or an error message if the derived object does
       *  not have this member
       */
      virtual double side_u_D1 () const { cosmobl::ErrorMsg("Error in Triplet::side_u_D1() of Triplet.h!"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_perc_increase_D1
       *  @return the ratio
       *  &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
       *  in the first dimension, or an error message if the derived
       *  object does not have this member
       */    
      virtual double perc_increase_D1 () const { cosmobl::ErrorMsg("Error in Triplet::perc_increase_D1() of Triplet.h!"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_nbins_D1
       *  @return the number of bins the first dimension, or an
       *  error message if the derived object does not have this
       *  member
       */    
      virtual int nbins_D1 () const { cosmobl::ErrorMsg("Error in Triplet::nbins_D1() of Triplet.h!"); return 0; }
      
      /**
       *  @brief get the protected member Triplet1D::m_binSize_D1
       *  @return the bin size in the first dimension, or an
       *  error message if the derived object does not have this
       *  member
       */    
      virtual double binSize_D1 () const { cosmobl::ErrorMsg("Error in Triplet::binSize_D1() of Triplet.h!"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_side_s_D2
       *  @return the size of r<SUB>12</SUB> in the second dimension,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual double side_s_D2 () const { cosmobl::ErrorMsg("Error in Triplet::side_s_D2() of Triplet.h!"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_side_u_D2
       *  @return the ratio r<SUB>13</SUB>/r<SUB>12</SUB> in the
       *  second dimension, or an error message if the derived object
       *  does not have this member
       */
      virtual double side_u_D2 () const { cosmobl::ErrorMsg("Error in Triplet::side_u_D2() of Triplet.h!"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_perc_increase_D2
       *  @return the ratio
       *  &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
       *  in the second dimension, or an error message if the derived
       *  object does not have this member
       */    
      virtual double perc_increase_D2 () const { cosmobl::ErrorMsg("Error in Triplet::perc_increase_D2() of Triplet.h!"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_nbins_D2
       *  @return the number of bins the second dimension, or an
       *  error message if the derived object does not have this
       *  member
       */    
      virtual int nbins_D2 () const { cosmobl::ErrorMsg("Error in Triplet::nbins_D2() of Triplet.h!"); return 0; }
      
      /**
       *  @brief get the protected member Triplet1D::m_binSize_D2
       *  @return the bin size in the second dimension, or an
       *  error message if the derived object does not have this
       *  member
       */    
      virtual double binSize_D2 () const { cosmobl::ErrorMsg("Error in Triplet::binSize_D2() of Triplet.h!"); return 0; }

      ///@}
      
      
      /**
       *  @name Member functions used to set private/protected members
       */
      ///@{

      /**
       *  @brief set the member m_TT1D[i]
       *  @param i the bin index
       *  @param tt the number of triplets in the bin
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_TT1D (const int i, const double tt) { cosmobl::ErrorMsg("Error in Triplet::set_TT() of Triplet.h!"); }

      /**
       *  @brief set the protected member Triplet1D::m_TT1D[i] adding
       *  the number of triplets
       *  @param i the bin index
       *  @param tt the number of triplets in the bin
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void add_TT1D (const int i, const double tt) { cosmobl::ErrorMsg("Error in Triplet::add_TT() of Triplet.h!"); }

      ///@}

      
      /**
       *  @name Member functions used to handle object triplets (customized in all the derived classes) 
       */
      ///@{
    
      /**
       *  @brief estimate the distance between two objects and update
       *  the triplet vectors accordingly
       *  @param r12 distance between object1 and object2
       *  @param r13 distance between object1 and object3
       *  @param r23 distance between object2 and object3
       *  @param ww the weight
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void put (const double r12, const double r13, const double r23, const double ww=1.) = 0;

      /**
       *  @brief estimate the distance between three objects and
       *  update the triplet vectors accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @param obj3 pointer to an object of class Object
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2, const shared_ptr<catalogue::Object> obj3) = 0;
    
      /**
       *  @brief sum the number of triplets
       *  @param tt pointer to an object of class Triplet
       *  @param ww the weight
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void Sum (const shared_ptr<Triplet> tt, const double ww=1.) = 0;

      ///@}

    }; 
  }
}

#endif
