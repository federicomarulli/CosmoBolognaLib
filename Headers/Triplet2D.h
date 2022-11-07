/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli, Michele Moresco,        *
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
 *  @file Headers/Triplet2D.h
 *
 *  @brief The class Triplet2D
 *
 *  This file defines the interface of the class Triplet2D*,
 *  used to handle 2D triplets of objects
 *  
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#ifndef __TRIPLET2D__
#define __TRIPLET2D__


#include "Triplet.h"


// ===================================================================================================


namespace cbl {
  
  namespace triplets {
    
    /**
     *  @class Triplet2D Triplet2D.h
     *  "Headers/Triplet2D.h"
     *
     *  @brief The class Triplet2D
     *
     *  This class is used to handle objects of type <EM>
     *  Triplet2D </EM>, used to handle 2D triplets of
     *  objects 
     *
     *  @warning This class has not been fully implemented yet
     */
    class Triplet2D : public Triplet {

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
      
      /// the binned scales in the first dimension
      std::vector<double> m_scale_D1;
      
      /// the binned scales in the second dimension
      std::vector<double> m_scale_D2;
      
      /// the number of binned triplets
      std::vector<std::vector<double> > m_TT2D;


       /**
       *  @name Binning parameters
       */
      ///@{

      /// the size of r<SUB>12</SUB> in the first dimension
      double m_r12_D1;

      /// the size of r<SUB>12</SUB> bin in the first dimension
      double m_r12_binSize_D1; 

      /// the size of r<SUB>13</SUB> in the first dimension
      double m_r13_D1;

      /// the size of r<SUB>13</SUB> bin in the first dimension
      double m_r13_binSize_D1; 

      /// the number of bins in the first dimension
      int m_nbins_D1;
      
      /// the bin size in the first dimension
      double m_binSize_D1;

      /// the size of r<SUB>12</SUB> in the second dimension
      double m_r12_D2;

      /// the size of r<SUB>12</SUB> bin in the second dimension
      double m_r12_binSize_D2; 

      /// the size of r<SUB>13</SUB> in the second dimension
      double m_r13_D2;

      /// the size of r<SUB>13</SUB> bin in the second dimension
      double m_r13_binSize_D2; 

      /// the number of bins in the second dimension
      int m_nbins_D2;
      
      /// the bin size in the second dimension
      double m_binSize_D2;
      
      ///@}

  
    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      Triplet2D () 
	: m_r12_D1(-1), m_r12_binSize_D1(-1), m_r13_D1(-1), m_r13_binSize_D1(-1), m_nbins_D1(-1), m_binSize_D1(-1), m_r12_D2(-1), m_r12_binSize_D2(-1), m_r13_D2(-1), m_r13_binSize_D2(-1), m_nbins_D2(-1), m_binSize_D2(-1)
	{ m_tripletDim = Dim::_2D_; }
                                                                
      /**
       *  @brief constructor
       *  @param r12_D1 the size of r<SUB>12</SUB> in the first
       *  dimension
       *  @param r12_binSize_D1 the size of r<SUB>12</SUB> bin in
       *  the first dimension
       *  @param r13_D1 the size of r<SUB>13</SUB> in the first
       *  dimension
       *  @param r13_binSize_D1 the size of r<SUB>13</SUB> bin in
       *  the first dimension
       *  @param nbins_D1 the number of bins in the first dimension
       *  @param r12_D2 the size of r<SUB>12</SUB> in the second
       *  dimension
       *  @param r12_binSize_D2 the size of r<SUB>12</SUB> bin in
       *  the second dimension
       *  @param r13_D2 the size of r<SUB>13</SUB> in the second
       *  dimension
       *  @param r13_binSize_D2 the size of r<SUB>13</SUB> bin in
       *  the second dimension
       *  @param nbins_D2 the number of bins in the second dimension
       */
      Triplet2D (const double r12_D1, const double r12_binSize_D1, const double r13_D1, const double r13_binSize_D1, const int nbins_D1, const double r12_D2, const double r12_binSize_D2, const double r13_D2, const double r13_binSize_D2, const int nbins_D2) 
	: m_r12_D1(r12_D1), m_r12_binSize_D1(r12_binSize_D1), m_r13_D1(r13_D1), m_r13_binSize_D1(r13_binSize_D1), m_nbins_D1(nbins_D1), m_r12_D2(r12_D2), m_r12_binSize_D2(r12_binSize_D2), m_r13_D2(r13_D2), m_r13_binSize_D2(r13_binSize_D2), m_nbins_D2(nbins_D2)
      { m_tripletDim = Dim::_2D_; m_TT2D.resize(m_nbins_D1+1, std::vector<double>(m_nbins_D2+1, 0.)); }

      /**
       *  @brief default destructor
       *  
       */
      ~Triplet2D () = default;
      
      ///@}

      
      /**
       *  @name Member functions used to get the protected members
       */
      ///@{

      /**
       *  @brief get the protected member \e m_scale_D1[i]
       *  @param i the bin index in the first dimension
       *  @return the i-th binned scale
       */
      double scale_D1 (const int i) const override { (void)i; return m_scale_D1[i]; }

      /**
       *  @brief get the protected member \e m_scale_D1
       *  @return the vector containing the binned scales
       */
      std::vector<double> scale_D1 () const override { return m_scale_D1; }

      /**
       *  @brief get the protected member \e m_scale_D2[i]
       *  @param i the bin index in the first dimension
       *  @return the i-th binned scale
       */
      double scale_D2 (const int i) const override { (void)i; return m_scale_D2[i]; }

      /**
       *  @brief get the protected member \e m_scale_D2
       *  @return the vector containing the binned scales
       */
      std::vector<double> scale_D2 () const override { return m_scale_D2; }

      /**
       *  @brief get the protected member \e m_PP2D[i]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the number of pairs in the i-th bin
       */
      double TT2D (const int i, const int j) const override { (void)i; (void)j; return m_TT2D[i][j]; }

      /**
       *  @brief get the protected member \e m_PP2D
       *  @return the vector containing the binned number of pairs
       */
      std::vector<std::vector<double> > TT2D () const override { return m_TT2D; }
      
      /**
       *  @brief get the protected member Triplet1D::m_r12_D1
       *  @return the size of r<SUB>12</SUB> in the first dimension
       */
      double r12_D1 () const override { return m_r12_D1; }

      /**
       *  @brief get the protected member Triplet1D::m_r12_binSize_D1
       *  @return the size of r<SUB>12</SUB> bin in the first dimension
       */
      double r12_binSize_D1 () const override { return m_r12_binSize_D1; }

      /**
       *  @brief get the protected member Triplet1D::m_r13_D1
       *  @return the size of r<SUB>13</SUB> in the first dimension
       */
      double r13_D1 () const override { return m_r13_D1; }

      /**
       *  @brief get the protected member Triplet1D::m_r13_binSize_D1
       *  @return the size of r<SUB>13</SUB> bin in the first dimension
       */
      double r13_binSize_D1 () const override { return m_r13_binSize_D1; }

      /**
       *  @brief get the protected member \e m_binSize_D1
       *  @return the bin size in the first dimension
       */    
      double binSize_D1 () const override { return m_binSize_D1; }

      /**
       *  @brief get the protected member Triplet1D::m_r12_D2
       *  @return the size of r<SUB>12</SUB> in the second dimension
       */
      double r12_D2 () const override { return m_r12_D2; }

      /**
       *  @brief get the protected member Triplet1D::m_r12_binSize_D2
       *  @return the size of r<SUB>12</SUB> bin in the second dimension
       */
      double r12_binSize_D2 () const override { return m_r12_binSize_D2; }

      /**
       *  @brief get the protected member Triplet1D::m_r13_D2
       *  @return the size of r<SUB>13</SUB> in the second dimension
       */
      double r13_D2 () const override { return m_r13_D2; }

      /**
       *  @brief get the protected member Triplet1D::m_r13_binSize_D2
       *  @return the size of r<SUB>13</SUB> bin in the second dimension
       */
      double r13_binSize_D2 () const override { return m_r13_binSize_D2; }

      /**
       *  @brief get the protected member \e m_binSize_D2
       *  @return the bin size in the second dimension
       */    
      double binSize_D2 () const override { return m_binSize_D2; }
      
      ///@}
      
    
      /**
       *  @name Member functions used to handle triplets
       */
      ///@{
    
      /**
       *  @brief sum the number of triplets
       *  @param tt pointer to an object of class Triplet
       *  @param ww the weight
       *  
       *
       *  @warning This method has not been implemented yet
       */
      void Sum (const std::shared_ptr<Triplet> tt, const double ww=1) override
      { (void)tt; (void)ww; ErrorCBL("", "Sum", "Triple2D.h", cbl::glob::ExitCode::_workInProgress_); }
     
      /**
       *  @brief estimate the distance between two objects and update
       *  the triplet vectors accordingly
       *  @param r12 distance between object1 and object2
       *  @param r13 distance between object1 and object3
       *  @param r23 distance between object2 and object3
       *  @param klin triplet bin
       */
      virtual void get_triplet (const double r12, const double r13, const double r23, int &klin) override
      { (void)r12; (void)r13; (void)r23; (void)klin; ErrorCBL("", "get_triplet", "Triple2D.h", cbl::glob::ExitCode::_workInProgress_); }

      /**
       *  @brief update the triplet
       *  @param klin triplet bin
       *  @param ww the weight
       */
      virtual void set_triplet (const int klin, const double ww=1.) override   
      { (void)klin; (void)ww; ErrorCBL("", "set_triplet", "Triple2D.h", cbl::glob::ExitCode::_workInProgress_); }

      /**
       *  @brief estimate the distance between three objects and
       *  update the triplet vectors accordingly
       *  @param r12 distance between object1 and object2
       *  @param r13 distance between object1 and object3
       *  @param r23 distance between object2 and object3
       *  @param ww total weight
       *
       *  @warning this method has not been implemented yet
       */
      void put (const double r12, const double r13, const double r23, const double ww=1.) override
      { (void)r12; (void)r13; (void)r23; (void)ww; ErrorCBL("", "put", "Triple2D.h", cbl::glob::ExitCode::_workInProgress_); }
    
      /**
       *  @brief estimate the distance between three objects and
       *  update the triplet vectors accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @param obj3 pointer to an object of class Object
       *  
       *  @warning This method has not been implemented yet
       */
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2, const std::shared_ptr<catalogue::Object> obj3) override
      { (void)obj1; (void)obj2; (void)obj3; ErrorCBL("", "put", "Triple2D.h", cbl::glob::ExitCode::_workInProgress_); }
      
      ///@}
    
    };
  }
}

#endif
