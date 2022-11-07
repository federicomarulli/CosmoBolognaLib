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
 *  @file Headers/Triplet.h
 *
 *  @brief The class Triplet
 *
 *  This file defines the interface of the class Triplet, used to handle
 *  triplets of objects to compute the three-point correlation function
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#ifndef __TRIPLET__
#define __TRIPLET__


#include "Catalogue.h"


// ===================================================================================================


namespace cbl {
  
  /**
   *  @brief The namespace of the functions and classes used to handle
   *  <B> triplets of objects </B>
   *  
   *  The \e triplets namespace contains all the functions and classes
   *  to handle triplets of objects
   */
  namespace triplets {
    
    /**
     * @enum TripletType
     * @brief the triplet type
     */
    enum class TripletType {

      /// 1D triplet in comoving coordinates and angular bins 
      _comoving_theta_,

      /// 1D triplet in comoving coordinates and linear bins 
      _comoving_side_,

      /// 1D triplet in comoving coordinates and linear bins of the cosine of theta
      _comoving_costheta_,

      /// multipoles of the triplets 
      _multipoles_direct_,

    };

    /**
     * @brief return a vector containing the
     * TripletType names
     * @return a vector containing the
     * TripletType names
     */
    inline std::vector<std::string> TripletTypeNames ()
    { return {"comoving_theta", "comoving_side", "comoving_costheta", "multipoles_direct"}; }

    /**
     * @brief cast an enum of type TripletType
     * from its index
     * @param tripletPTypeIndex the tripletPType index
     * @return object of class TripletType
     */
    inline TripletType TripletTypeCast (const int tripletPTypeIndex)
    { return castFromValue<TripletType>(tripletPTypeIndex); }

    /**
     * @brief cast an enum of type TripletType
     * from its name
     * @param tripletPTypeName the tripletPType name
     * @return object of class TripletType
     */
    inline TripletType TripletTypeCast (const std::string tripletPTypeName)
    { return castFromName<TripletType>(tripletPTypeName, TripletTypeNames()); }

    /**
     * @brief cast an enum of type TripletType
     * from indeces
     * @param tripletPTypeIndeces the tripletPType indeces
     * @return vector of objects of class TripletType
     */
    inline std::vector<TripletType> TripletTypeCast (const std::vector<int> tripletPTypeIndeces)
    { return castFromValues<TripletType>(tripletPTypeIndeces); } 

    /**
     * @brief cast an enum of type TripletType
     * from thier names
     * @param tripletPTypeNames the tripletPType names
     * @return vector of TripletType enums
     */
    inline std::vector<TripletType> TripletTypeCast (const std::vector<std::string> tripletPTypeNames)
    { return castFromNames<TripletType>(tripletPTypeNames, TripletTypeNames()); }


    /**
     *  @class Triplet Triplet.h "Headers/Triplet.h"
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
       *  @name Member functions used to set the binning parameters
       *  (customized in all the derived classes)
       */
      ///@{
  
      /**
       *  @brief set the binning parameters
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
       *  
       */
      virtual ~Triplet () {}

      /**
       *  @brief static factory used to construct triplets of any type
       *  
       *  @param type the triplet type; it can be: _comoving_theta,
       *  _comoving_side
       *  @param r12 the size of r<SUB>12</SUB>
       *  @param r12_binSize the r<SUB>12</SUB> bin size
       *  @param r13 the ratio r<SUB>13</SUB>
       *  @param r13_binSize the r<SUB>12</SUB> bin size
       *  @param nbins the number of bins
       *  @return a pointer to an object of class Triplet of a given type
       */
      static std::shared_ptr<Triplet> Create (const TripletType type, const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int nbins);

      ///@}
    

      /**
       *  @name Member functions used to get private/protected parameters
       */
      ///@{
      
      /**
       *  @brief get the dimension of the triplet vectors 
       *  @return the dimension of the triplet vectors 
       */
      Dim tripletDim () const
      { return m_tripletDim; }

      /**
       *  @brief get the triplet type
       *  @return the triplet type
       */
      TripletType tripletType () const
      { return m_tripletType; }

      /**
       *  @brief get the protected member \e m_scale[i]
       *  @param i the bin index
       *  @return the i-th binned scale
       */
      virtual double scale (const int i) const
      { (void)i; cbl::ErrorCBL("", "scale", "Triplet.h"); return 0; }

      /**
       *  @brief get the protected member \e m_scale
       *  @return the vector containing the binned scales
       */
      virtual std::vector<double> scale () const
      { cbl::ErrorCBL("", "scale", "Triplet.h"); std::vector<double> vv; return vv; }
      
      /**
       *  @brief get the private member \e m_TT1D[i]
       *  @param i the bin index
       *  @return the number of triplets in the i-th angular bin, or an
       *  error message if the derived object does not have this member
       */
      virtual double TT1D (const int i) const { (void)i; cbl::ErrorCBL("", "TT1D", "Triplet.h"); return 0; }
      
      /**
       *  @brief get the private member \e m_TT1D
       *  @return the vector containing the number of triplets in linear
       *  bins, or an error message if the derived object does not have
       *  this member
       */
      virtual std::vector<double> TT1D () const { cbl::ErrorCBL("", "TT1D", "Triplet.h"); std::vector<double> vv; return vv; }

      /**
       *  @brief get the protected member \e m_scale_D1[i]
       *  @param i the bin index in the first dimension
       *  @return the i-th binned scale
       */
      virtual double scale_D1 (const int i) const { (void)i; cbl::ErrorCBL("", "scale_D1", "Triplet.h"); return 0; }

      /**
       *  @brief get the protected member \e m_scale_D1
       *  @return the vector containing the binned scales
       */
      virtual std::vector<double> scale_D1 () const { cbl::ErrorCBL("", "scale_D1", "Triplet.h"); std::vector<double> vv; return vv; }

      /**
       *  @brief get the protected member \e m_scale_D2[i]
       *  @param i the bin index in the first dimension
       *  @return the i-th binned scale
       */
      virtual double scale_D2 (const int i) const { (void)i; cbl::ErrorCBL("", "scale_D2", "Triplet.h"); return 0; }

      /**
       *  @brief get the protected member \e m_scale_D2
       *  @return the vector containing the binned scales
       */
      virtual std::vector<double> scale_D2 () const { cbl::ErrorCBL("", "scale_D2", "Triplet.h"); std::vector<double> vv; return vv; }

      /**
       *  @brief get the protected member \e m_TT2D[i]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the number of pairs in the i-th bin
       */
      virtual double TT2D (const int i, const int j) const { (void)i; (void)j; cbl::ErrorCBL("", "TT2D", "Triplet.h"); return 0; }

      /**
       *  @brief get the protected member \e m_TT2D
       *  @return the vector containing the binned number of pairs
       */
      virtual std::vector<std::vector<double> > TT2D () const { cbl::ErrorCBL("", "TT2D", "Triplet.h"); std::vector<std::vector<double> > vv; return vv; }
      
      /**
       *  @brief get the private member \e m_r12
       *  @return the size of r<SUB>12</SUB>, or an error message if the
       *  derived object does not have this member
       */
      virtual double r12 () const { cbl::ErrorCBL("", "r12", "Triplet.h"); return 0; }
      
      /**
       *  @brief get the private member \e m_r12
       *  @return the size of the r<SUB>12</SUB> bin, or an error message if the
       *  derived object does not have this member
       */
      virtual double r12_binSize () const { cbl::ErrorCBL("", "r12_binSize", "Triplet.h"); return 0; }

      /**
       *  @brief get the private member \e m_side_u
       *  @return the size of r<SUB>13</SUB>, or an error message if the
       *  derived object does not have this member
       */
      virtual double r13 () const { cbl::ErrorCBL("", "r13", "Triplet.h"); return 0; }

      /**
       *  @brief get the private member \e m_side_u
       *  @return the size of r<SUB>13</SUB> bin, or an error message if the
       *  derived object does not have this member
       */
      virtual double r13_binSize () const { cbl::ErrorCBL("", "r13_binSize", "Triplet.h"); return 0; }

      /**
       *  @brief get the private member \e m_nbins
       *  @return the number of bins, or an error message if the
       *  derived object does not have this member
       */
      virtual int nbins () const { cbl::ErrorCBL("", "nbins", "Triplet.h"); return 0; }

      /**
       *  @brief get the private member \e m_binSize
       *  @return the bin size, or an error message if the
       *  derived object does not have this member
       */
      virtual double binSize () const { cbl::ErrorCBL("", "binSize", "Triplet.h"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_r12_D1
       *  @return the size of r<SUB>12</SUB> in the first dimension,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual double r12_D1 () const { cbl::ErrorCBL("", "r12_D1", "Triplet.h"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_r12_binSize_D1
       *  @return the size of r<SUB>12</SUB> bin in the first dimension,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual double r12_binSize_D1 () const { cbl::ErrorCBL("", "r12_binSize_D1", "Triplet.h"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_r13_D1
       *  @return the size of r<SUB>13</SUB> in the first dimension,
       *  or an error message if the derived object does
       *  not have this member
       */
      virtual double r13_D1 () const { cbl::ErrorCBL("", "r13_D1", "Triplet.h"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_r13_binSize_D1
       *  @return the size of r<SUB>13</SUB> bin in the first dimension,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual double r13_binSize_D1 () const { cbl::ErrorCBL("", "r13_binSize_D1", "Triplet.h"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_nbins_D1
       *  @return the number of bins the first dimension, or an
       *  error message if the derived object does not have this
       *  member
       */    
      virtual int nbins_D1 () const { cbl::ErrorCBL("", "nbins_D1", "Triplet.h"); return 0; }
      
      /**
       *  @brief get the protected member Triplet1D::m_binSize_D1
       *  @return the bin size in the first dimension, or an
       *  error message if the derived object does not have this
       *  member
       */    
      virtual double binSize_D1 () const { cbl::ErrorCBL("", "binSize_D1", "Triplet.h"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_r12_D1
       *  @return the size of r<SUB>12</SUB> in the second dimension,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual double r12_D2 () const { cbl::ErrorCBL("", "r12_D2", "Triplet.h"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_r12_binSize_D1
       *  @return the size of r<SUB>12</SUB> bin in the second dimension,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual double r12_binSize_D2 () const { cbl::ErrorCBL("", "r12_binSize_D2", "Triplet.h"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_r13_D2
       *  @return the size of r<SUB>13</SUB> in the second dimension,
       *  or an error message if the derived object does
       *  not have this member
       */
      virtual double r13_D2 () const { cbl::ErrorCBL("", "r13_D2", "Triplet.h"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_r13_binSize_D2
       *  @return the size of r<SUB>13</SUB> bin in the second dimension,
       *  or an error message if the derived object does not have this
       *  member
       */
      virtual double r13_binSize_D2 () const { cbl::ErrorCBL("", "r13_binSize_D2", "Triplet.h"); return 0; }

      /**
       *  @brief get the protected member Triplet1D::m_nbins_D2
       *  @return the number of bins the second dimension, or an
       *  error message if the derived object does not have this
       *  member
       */    
      virtual int nbins_D2 () const { cbl::ErrorCBL("", "nbins_D2", "Triplet.h"); return 0; }
      
      /**
       *  @brief get the protected member Triplet1D::m_binSize_D2
       *  @return the bin size in the second dimension, or an
       *  error message if the derived object does not have this
       *  member
       */    
      virtual double binSize_D2 () const { cbl::ErrorCBL("", "binSize_D2", "Triplet.h"); return 0; }

      ///@}
      
      
      /**
       *  @name Member functions used to set private/protected members
       */
      ///@{

      /**
       *  @brief set the member m_TT1D[i]
       *  @param i the bin index
       *  @param tt the number of triplets in the bin
       */
      virtual void set_TT1D (const int i, const double tt)
      { (void)i; (void)tt; cbl::ErrorCBL("", "set_TT1D", "Triplet.h"); }

      /**
       *  @brief set the protected member Triplet1D::m_TT1D[i] adding
       *  the number of triplets
       *  @param i the bin index
       *  @param tt the number of triplets in the bin
       */
      virtual void add_TT1D (const int i, const double tt)
      { (void)i; (void)tt; cbl::ErrorCBL("", "add_TT1D", "Triplet.h"); }

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
       *  @param klin triplet bin
       */
      virtual void get_triplet (const double r12, const double r13, const double r23, int &klin) = 0;   

      /**
       *  @brief update the triplet
       *  @param klin triplet bin
       *  @param ww the weight
       */
      virtual void set_triplet (const int klin, const double ww=1.) = 0;   

      /**
       *  @brief estimate the distance between two objects and update
       *  the triplet vectors accordingly
       *  @param r12 distance between object1 and object2
       *  @param r13 distance between object1 and object3
       *  @param r23 distance between object2 and object3
       *  @param ww the weight
       */
      virtual void put (const double r12, const double r13, const double r23, const double ww=1.) = 0;

      /**
       *  @brief estimate the distance between three objects and
       *  update the triplet vectors accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @param obj3 pointer to an object of class Object
       */
      virtual void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2, const std::shared_ptr<catalogue::Object> obj3) = 0;
    
      /**
       *  @brief sum the number of triplets
       *  @param tt pointer to an object of class Triplet
       *  @param ww the weight
       */
      virtual void Sum (const std::shared_ptr<Triplet> tt, const double ww=1.) = 0;

      ///@}

    }; 
  }
}

#endif
