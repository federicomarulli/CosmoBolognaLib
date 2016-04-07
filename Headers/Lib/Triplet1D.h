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
 *  @file Headers/Lib/Triplet1D.h
 *
 *  @brief The class Triplet1D
 *
 *  This file defines the interface of the class Triplet1D*,
 *  used to handle 1D triplets of objects
 *  
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#ifndef __TRIPLET1D__
#define __TRIPLET1D__


#include "Triplet.h"


// ===================================================================================================


namespace cosmobl {
  
  namespace triplets {
    
    /**
     *  @class Triplet1D Triplet1D.h
     *  "Headers/Lib/Triplet1D.h"
     *
     *  @brief The class Triplet1D
     *
     *  This class is used to handle objects of type <EM>
     *  Triplet1D </EM>, used to handle 1D triplets of
     *  objects 
     */
    class Triplet1D : public Triplet {

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

      /// the binned scales
      vector<double> m_scale;
      
      /// the number of binned triplets
      vector<double> m_TT1D;

    
      /**
       *  @name Binning parameters
       */
      ///@{

      /// the size of r<SUB>12</SUB>
      double m_side_s;

      /// the ratio r<SUB>13</SUB>/<SUB>r12</SUB>
      double m_side_u; 

      /// the ratio &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
      double m_perc_increase;
      
      /// the number of bins
      int m_nbins;

      /// the bin size
      double m_binSize;
      
      ///@}

  
    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Triplet1D, with protected
       *  members set to -1
       */
      Triplet1D () 
	: m_side_s(-1), m_side_u(-1.), m_perc_increase(-1.), m_nbins(-1), m_binSize(-1)
	{ m_tripletDim = _1D_; }
                                                                
      /**
       *  @brief constructor
       *  @param side_s the size of r<SUB>12</SUB>
       *  @param side_u the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
       *  @param perc_increase the ratio &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
       *  @param nbins the number of bins
       *  @return object of class Triplet1D
       */
      Triplet1D (const double side_s, const double side_u, const double perc_increase, const int nbins) 
	: m_side_s(side_s), m_side_u(side_u), m_perc_increase(perc_increase), m_nbins(nbins)
	{ m_tripletDim = _1D_; }

      /**
       *  @brief default destructor
       *  @return none
       */
      ~Triplet1D () = default;
      
      ///@}

    
      /**
       *  @name Member functions used to get the protected members
       */
      ///@{

      /**
       *  @brief get the protected member \e m_scale[i]
       *  @param i the bin index
       *  @return the i-th binned scale
       */
      double scale (const int i) const override { return m_scale[i]; }

      /**
       *  @brief get the protected member \e m_scale
       *  @return the vector containing the binned scales
       */
      vector<double> scale () const override { return m_scale; }
      
      /**
       *  @brief get the private member \e m_TT1D[i]
       *  @param i the bin index
       *  @return the number of triplets in the i-th bin
       */
      double TT1D (const int i) const override { return m_TT1D[i]; }

      /**
       *  @brief get the private member \e m_TT1D
       *  @return the vector containing the number of triplets
       */
      vector<double> TT1D () const override { return m_TT1D; }
      
      /**
       *  @brief get the protected member \e m_side_s
       *  @return the size of r<SUB>12</SUB>
       */
      double side_s () const override { return m_side_s; }

      /**
       *  @brief get the protected member \e m_side_u
       *  @return the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
       */
      double side_u () const override { return m_side_u; }

      /**
       *  @brief get the protected member \e m_perc_increase
       *  @return the ratio &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
       */    
      double perc_increase () const override { return m_perc_increase; }

      /**
       *  @brief get the protected member \e m_nbins
       *  @return the number of bins
       */    
      int nbins () const override { return m_nbins; }
      
      /**
       *  @brief get the protected member \e m_binSize
       *  @return the bin size
       */    
      double binSize () const override { return m_binSize; }

      ///@}

      
      /**
       *  @name Member functions used to set private/protected members
       */
      ///@{

      /**
       *  @brief set the member m_TT1D[i]
       *  @param i the bin index
       *  @param tt the number of triplets in the bin
       *  @return none
       */
      void set_TT1D (const int i, const double tt) { checkDim(m_TT1D, i, "m_TT1D"); m_TT1D[i] = tt; }

      /**
       *  @brief set the protected member \e m_TT1D[i] adding
       *  the number of triplets
       *  @param i the bin index
       *  @param tt the number of triplets in the bin
       *  @return none
       */
      void add_TT1D (const int i, const double tt) { checkDim(m_TT1D, i, "m_TT1D"); m_TT1D[i] += tt; }
      
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
       */
      void Sum (const shared_ptr<Triplet> tt, const double ww=1) override;

      ///@}
    
    };

    
    // ============================================================================================
    // ============================================================================================

    
    /**
     *  @class Triplet1D_comoving Triplet1D.h
     *  "Headers/Lib/Triplet1D.h"
     *
     *  @brief The class Triplet1D_comoving
     *
     *  This class is used to handle objects of type <EM>
     *  Triplet1D_comoving </EM>, used to handle 1D triplets of
     *  objects in comoving coordinates
     */
    class Triplet1D_comoving : public Triplet1D {

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

      
    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Triplet1D_comoving
       */
      Triplet1D_comoving () = default;
                                                                
      /**
       *  @brief constructor
       *  @param side_s the size of r<SUB>12</SUB>
       *  @param side_u the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
       *  @param perc_increase the ratio &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
       *  @param nbins the number of bins
       *  @return object of class Triplet1D_comoving
       */
      Triplet1D_comoving (const double side_s, const double side_u, const double perc_increase, const int nbins) 
	: Triplet1D(side_s, side_u, perc_increase, nbins) {}
      
      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Triplet1D_comoving () = default;

      ///@}
    
    };

    
    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Triplet1D_comoving_theta Triplet1D.h
     *  "Headers/Lib/Triplet1D.h"
     *
     *  @brief The class Triplet1D_comoving_theta
     *
     *  This class is used to handle objects of type <EM>
     *  Triplet1D_comoving_theta </EM>, used to handle 1D triplets of
     *  objects in comoving coordinates, in angular bins
     */
    class Triplet1D_comoving_theta : public Triplet1D_comoving {
      
    private:

      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters
       *  @return none
       *
       *  @warning This method has not been implemented yet
       */
      void set_parameters () override;
      
      ///@}

      
    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Triplet1D_comoving_theta, with protected
       *  members set to -1
       */
      Triplet1D_comoving_theta () 
	{ m_tripletType = _comoving_theta_; }
                                                                
      /**
       *  @brief constructor
       *  @param side_s the size of r<SUB>12</SUB>
       *  @param side_u the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
       *  @param perc_increase the ratio &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
       *  @param nbins the number of bins
       *  @return object of class Triplet1D_comoving_theta
       */
      Triplet1D_comoving_theta (const double side_s, const double side_u, const double perc_increase, const int nbins) 
	: Triplet1D_comoving(side_s, side_u, perc_increase, nbins)
	{ m_tripletType = _comoving_theta_; set_parameters(); m_scale.resize(m_nbins+1, 0.); m_TT1D.resize(nbins+1, 0.); }

      ///@}

    
      /**
       *  @name Member functions used to handle triplets
       */
      ///@{

      /**
       *  @brief estimate the distance between three objects and
       *  update the triplet vectors accordingly
       *  @param r12 distance between object1 and object2
       *  @param r13 distance between object1 and object3
       *  @param r23 distance between object2 and object3
       *  @param ww total weight
       *  @return none
       */
      void put (const double r12, const double r13, const double r23, const double ww=1.) override;
    
      /**
       *  @brief estimate the distance between three objects and
       *  update the triplet vectors accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @param obj3 pointer to an object of class Object
       *  @return none
       */
      void put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2, const shared_ptr<catalogue::Object> obj3) override;

      ///@}
    
    };


    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Triplet1D_comoving_side Triplet1D.h
     *  "Headers/Lib/Triplet1D.h"
     *
     *  @brief The class Triplet1D_comoving_side
     *
     *  This class is used to handle objects of type <EM>
     *  Triplet1D_comoving_side </EM>, used to handle 1D triplets of
     *  objects in comoving coordinates linearly binned
     */
    class Triplet1D_comoving_side : public Triplet1D_comoving {

    private:

      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters
       *  @return none
       *
       *  @warning This method has not been implemented yet
       */
      void set_parameters () override;

      ///@}

      
    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Triplet1D_comoving_side, with
       *  protected members set to -1
       */
      Triplet1D_comoving_side () 
	{ m_tripletType = _comoving_side_; }
      
      /**
       *  @brief constructor
       *  @param side_s the size of r<SUB>12</SUB>
       *  @param side_u the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
       *  @param perc_increase the ratio &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
       *  @param nbins the number of bins
       *  @return object of class Triplet1D_comoving_side
       */
      Triplet1D_comoving_side (const double side_s, const double side_u, const double perc_increase, const int nbins) 
	: Triplet1D_comoving(side_s, side_u, perc_increase, nbins)
      { m_tripletType = _comoving_side_; set_parameters(); m_TT1D.resize(nbins+1, 0.); }

      ///@}

      
      /**
       *  @name Member functions used to handle triplets
       */
      ///@{

      /**
       *  @brief estimate the distance between three objects and
       *  update the triplet vectors accordingly
       *  @param r12 distance between object1 and object2
       *  @param r13 distance between object1 and object3
       *  @param r23 distance between object2 and object3
       *  @param ww total weight
       *  @return none
       */
      void put (const double r12, const double r13, const double r23, const double ww=1.) override;
    
      /**
       *  @brief estimate the distance between three objects and
       *  update the triplet vectors accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @param obj3 pointer to an object of class Object
       *  @return none
       */
      void put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2, const shared_ptr<catalogue::Object> obj3) override;

      ///@}
      
    };

    
    // ============================================================================================
    // ============================================================================================


     /**
     *  @class Triplet1D_angular Triplet1D.h
     *  "Headers/Lib/Triplet1D.h"
     *
     *  @brief The class Triplet1D_angular
     *
     *  This class is used to handle objects of type <EM>
     *  Triplet1D_angular </EM>, used to handle 1D triplets of
     *  objects in angular coordinates
     *
     *  @warning This class has not been fully implemented yet
     */
    class Triplet1D_angular : public Triplet1D {

    private:

      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters
       *  @return none
       *
       *  @warning This method has not been implemented yet
       */
      virtual void set_parameters () = 0;

      ///@}

      
    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Triplet1D_angular
       */
      Triplet1D_angular () = default;
                                                                
      /**
       *  @brief constructor
       *  @param side_s the size of r<SUB>12</SUB>
       *  @param side_u the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
       *  @param perc_increase the ratio &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
       *  @param nbins the number of bins
       *  @return object of class Triplet1D_angular
       */
      Triplet1D_angular (const double side_s, const double side_u, const double perc_increase, const int nbins) 
	: Triplet1D(side_s, side_u, perc_increase, nbins) {}

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Triplet1D_angular () = default;

      ///@}
    
    };

  }
}

#endif
