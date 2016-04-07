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
 *  @file Headers/Lib/Triplet2D.h
 *
 *  @brief The class Triplet2D
 *
 *  This file defines the interface of the class Triplet2D*,
 *  used to handle 2D triplets of objects
 *  
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#ifndef __TRIPLET2D__
#define __TRIPLET2D__


#include "Triplet.h"


// ===================================================================================================


namespace cosmobl {
  
  namespace triplets {
    
    /**
     *  @class Triplet2D Triplet2D.h
     *  "Headers/Lib/Triplet2D.h"
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
      
      /// the binned scales in the first dimension
      vector<double> m_scale_D1;
      
      /// the binned scales in the second dimension
      vector<double> m_scale_D2;
      
      /// the number of binned triplets
      vector<vector<double> > m_TT2D;


       /**
       *  @name Binning parameters
       */
      ///@{

      /// the size of r<SUB>12</SUB> in the first dimension
      double m_side_s_D1;

      /// the ratio r<SUB>13</SUB>/<SUB>r12</SUB> in the first dimension
      double m_side_u_D1; 

      /// the ratio &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB> in the first dimension
      double m_perc_increase_D1;

      /// the number of bins in the first dimension
      int m_nbins_D1;
      
      /// the bin size in the first dimension
      double m_binSize_D1;

       /// the size of r<SUB>12</SUB> in the second dimension
      double m_side_s_D2;

      /// the ratio r<SUB>13</SUB>/<SUB>r12</SUB> in the second dimension
      double m_side_u_D2; 

      /// the ratio &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB> in the second dimension
      double m_perc_increase_D2;

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
       *  @return object of class Triplet2D, with protected
       *  members set to -1
       */
      Triplet2D () 
	: m_side_s_D1(-1), m_side_u_D1(-1.), m_perc_increase_D1(-1.), m_nbins_D1(-1), m_binSize_D1(-1), m_side_s_D2(-1), m_side_u_D2(-1.), m_perc_increase_D2(-1.), m_nbins_D2(-1), m_binSize_D2(-1)
	{ m_tripletDim = _2D_; }
                                                                
      /**
       *  @brief constructor
       *  @param side_s_D1 the size of r<SUB>12</SUB> in the first
       *  dimension
       *  @param side_u_D1 the ratio r<SUB>13</SUB>/r<SUB>12</SUB> in
       *  the first dimension
       *  @param perc_increase_D1 the ratio
       *  &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
       *  in the first dimension
       *  @param nbins_D1 the number of bins in the first dimension
       *  @param side_s_D2 the size of r<SUB>12</SUB> in the second
       *  dimension
       *  @param side_u_D2 the ratio r<SUB>13</SUB>/r<SUB>12</SUB> in
       *  the second dimension
       *  @param perc_increase_D2 the ratio
       *  &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
       *  in the second dimension
       *  @param nbins_D2 the number of bins in the second dimension
       *  @return object of class Triplet2D
       */
      Triplet2D (const double side_s_D1, const double side_u_D1, const double perc_increase_D1, const int nbins_D1, const double side_s_D2, const double side_u_D2, const double perc_increase_D2, const int nbins_D2) 
	: m_side_s_D1(side_s_D1), m_side_u_D1(side_u_D1), m_perc_increase_D1(perc_increase_D1), m_nbins_D1(nbins_D1), m_side_s_D2(side_s_D2), m_side_u_D2(side_u_D2), m_perc_increase_D2(perc_increase_D2), m_nbins_D2(nbins_D2)
      { m_tripletDim = _2D_; m_TT2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); }

      /**
       *  @brief default destructor
       *  @return none
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
      double scale_D1 (const int i) const override { return m_scale_D1[i]; }

      /**
       *  @brief get the protected member \e m_scale_D1
       *  @return the vector containing the binned scales
       */
      vector<double> scale_D1 () const override { return m_scale_D1; }

      /**
       *  @brief get the protected member \e m_scale_D2[i]
       *  @param i the bin index in the first dimension
       *  @return the i-th binned scale
       */
      double scale_D2 (const int i) const override { return m_scale_D2[i]; }

      /**
       *  @brief get the protected member \e m_scale_D2
       *  @return the vector containing the binned scales
       */
      vector<double> scale_D2 () const override { return m_scale_D2; }

      /**
       *  @brief get the protected member \e m_PP2D[i]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @return the number of pairs in the i-th bin
       */
      double TT2D (const int i, const int j) const override { return m_TT2D[i][j]; }

      /**
       *  @brief get the protected member \e m_PP2D
       *  @return the vector containing the binned number of pairs
       */
      vector<vector<double> > TT2D () const override { return m_TT2D; }
      
      /**
       *  @brief get the protected member \e m_side_s_D1
       *  @return the size of r<SUB>12</SUB> in the first dimension
       */
      double side_s_D1 () const override { return m_side_s_D1; }
      
      /**
       *  @brief get the protected member \e m_side_u_D1
       *  @return the ratio r<SUB>13</SUB>/r<SUB>12</SUB> in the first
       *  dimension
       */
      double side_u_D1 () const override { return m_side_u_D1; }

      /**
       *  @brief get the protected member \e m_perc_increase_D1
       *  @return the ratio
       *  &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
       *  in the first dimension
       */    
      double perc_increase_D1 () const override { return m_perc_increase_D1; }

      /**
       *  @brief get the protected member \e m_binSize_D1
       *  @return the bin size in the first dimension
       */    
      double binSize_D1 () const override { return m_binSize_D1; }

      /**
       *  @brief get the protected member \e m_side_s_D2
       *  @return the size of r<SUB>12</SUB> in the second dimension
       */
      double side_s_D2 () const override { return m_side_s_D2; }

      /**
       *  @brief get the protected member \e m_side_u_D2
       *  @return the ratio r<SUB>13</SUB>/r<SUB>12</SUB> in the second
       *  dimension
       */
      double side_u_D2 () const override { return m_side_u_D2; }

      /**
       *  @brief get the protected member \e m_perc_increase_D2
       *  @return the ratio
       *  &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
       *  in the second dimension
       */    
      double perc_increase_D2 () const override { return m_perc_increase_D2; }

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
       *  @return none
       *
       *  @warning This method has not been implemented yet
       */
      void Sum (const shared_ptr<Triplet> tt, const double ww=1) override { ErrorMsg("Work in progress..."); }

      /**
       *  @brief estimate the distance between three objects and
       *  update the triplet vectors accordingly
       *  @param r12 distance between object1 and object2
       *  @param r13 distance between object1 and object3
       *  @param r23 distance between object2 and object3
       *  @param ww total weight
       *  @return none
       *
       *  @warning This method has not been implemented yet
       */
      void put (const double r12, const double r13, const double r23, const double ww=1.) override { ErrorMsg("Work in progress..."); }
    
      /**
       *  @brief estimate the distance between three objects and
       *  update the triplet vectors accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @param obj3 pointer to an object of class Object
       *  @return none
       *
       *  @warning This method has not been implemented yet
       */
      void put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2, const shared_ptr<catalogue::Object> obj3) override { ErrorMsg("Work in progress..."); }
      
      ///@}
    
    };
  }
}

#endif
