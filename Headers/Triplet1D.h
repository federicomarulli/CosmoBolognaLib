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
 *  @file Headers/Triplet1D.h
 *
 *  @brief The class Triplet1D
 *
 *  This file defines the interface of the class Triplet1D*,
 *  used to handle 1D triplets of objects
 *  
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#ifndef __TRIPLET1D__
#define __TRIPLET1D__


#include "Triplet.h"


// ===================================================================================================


namespace cbl {
  
  namespace triplets {
    
    /**
     *  @class Triplet1D Triplet1D.h
     *  "Headers/Triplet1D.h"
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
      std::vector<double> m_scale;
      
      /// the number of binned triplets
      std::vector<double> m_TT1D;

    
      /**
       *  @name Binning parameters
       */
      ///@{

      /// the size of r<SUB>12</SUB>
      double m_r12;

      /// the size of r<SUB>12</SUB> bin
      double m_r12_binSize;

      /// the size of r<SUB>13</SUB>
      double m_r13;
   
      /// the size of r<SUB>13</SUB> bin
      double m_r13_binSize;   

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
	: m_r12(-1), m_r12_binSize(-1), m_r13(-1), m_r13_binSize(-1), m_nbins(-1), m_binSize(-1)
	{ m_tripletDim = Dim::_1D_; }
                                                                
      /**
       *  @brief constructor
       *  @param r12 the size of r<SUB>12</SUB>
       *  @param r12_binSize the size of r<SUB>12</SUB> bin
       *  @param r13 the size of r<SUB>13</SUB>
       *  @param r13_binSize the size of r<SUB>13</SUB> bin
       *  @param nbins the number of bins
       *  @return object of class Triplet1D
       */
      Triplet1D (const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int nbins) 
	: m_r12(r12), m_r12_binSize(r12_binSize), m_r13(r13), m_r13_binSize(r13_binSize) , m_nbins(nbins)
	{ m_tripletDim = Dim::_1D_; }

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
      std::vector<double> scale () const override { return m_scale; }
      
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
      std::vector<double> TT1D () const override { return m_TT1D; }
      
      /**
       *  @brief get the private member \e m_r12
       *  @return the size of r<SUB>12</SUB>
       */
      double r12 () const override { return m_r12; }
      
      /**
       *  @brief get the private member \e m_r12_binSize
       *  @return the size of the r<SUB>12</SUB> bin
       */
      double r12_binSize () const override { return m_r12_binSize; }

      /**
       *  @brief get the private member \e m_r13
       *  @return the size of r<SUB>13</SUB>
       */
      double r13 () const override { return m_r13;  }

      /**
       *  @brief get the private member \e m_r13_binSize
       *  @return the size of r<SUB>13</SUB> bin
       */
      double r13_binSize () const override { return m_r13_binSize; }

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
      void Sum (const std::shared_ptr<Triplet> tt, const double ww=1) override;

      ///@}
    
    };

    
    // ============================================================================================
    // ============================================================================================

    
    /**
     *  @class Triplet1D_comoving Triplet1D.h
     *  "Headers/Triplet1D.h"
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
       *  @param r12 the size of r<SUB>12</SUB>
       *  @param r12_binSize the size of r<SUB>12</SUB> bin
       *  @param r13 the size of r<SUB>13</SUB>
       *  @param r13_binSize the size of r<SUB>13</SUB> bin
       *  @param nbins the number of bins
       *  @return object of class Triplet1D_comoving
       */
      Triplet1D_comoving (const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int nbins) 
	: Triplet1D(r12, r12_binSize, r13, r13_binSize, nbins) {}
      
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
     *  "Headers/Triplet1D.h"
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
	{ m_tripletType = TripletType::_comoving_theta_; }
                                                                
      /**
       *  @brief constructor
       *  @param r12 the size of r<SUB>12</SUB>
       *  @param r12_binSize the size of r<SUB>12</SUB> bin
       *  @param r13 the size of r<SUB>13</SUB>
       *  @param r13_binSize the size of r<SUB>13</SUB> bin
       *  @param nbins the number of bins
       *  @return object of class Triplet1D_comoving_theta
       */
      Triplet1D_comoving_theta (const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int nbins) 
	: Triplet1D_comoving(r12, r12_binSize, r13, r13_binSize, nbins)
	{ m_tripletType = TripletType::_comoving_theta_; set_parameters(); m_scale.resize(m_nbins+1, 0.); m_TT1D.resize(nbins+1, 0.); }

      ///@}

    
      /**
       *  @name Member functions used to handle triplets
       */
      ///@{
     
      /**
       *  @brief estimate the distance between two objects and update
       *  the triplet vectors accordingly
       *  @param r12 distance between object1 and object2
       *  @param r13 distance between object1 and object3
       *  @param r23 distance between object2 and object3
       *  @param klin triplet bin
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      void get_triplet (const double r12, const double r13, const double r23, int &klin) override;   

      /**
       *  @brief update the triplet
       *  @param klin triplet bin
       *  @param ww the weight
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      void set_triplet (const int klin, const double ww=1.) override;

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
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2, const std::shared_ptr<catalogue::Object> obj3) override;

      ///@}
    
    };

    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Triplet1D_multipoles_direct Triplet1D.h
     *  "Headers/Triplet1D.h"
     *
     *  @brief The class Triplet1D_multipoles_direct
     *
     *  This class is used to handle objects of type <EM>
     *  Triplet1D_comoving_theta </EM>, used to handle 1D triplets of
     *  objects in comoving coordinates, in angular bins
     */
    class Triplet1D_multipoles_direct : public Triplet1D_comoving {
      
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
       *  @return object of class Triplet1D_multipoles_direct, with protected
       *  members set to -1
       */
      Triplet1D_multipoles_direct () 
	{ m_tripletType = TripletType::_multipoles_direct_; }
                                                                
      /**
       *  @brief constructor
       *  @param r12 the size of r<SUB>12</SUB>
       *  @param r12_binSize the size of r<SUB>12</SUB> bin
       *  @param r13 the size of r<SUB>13</SUB>
       *  @param r13_binSize the size of r<SUB>13</SUB> bin
       *  @param norders the number of Legendre multipoles
       *  @return object of class Triplet1D_multipoles_direct
       */
      Triplet1D_multipoles_direct (const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int norders) 
	: Triplet1D_comoving(r12, r12_binSize, r13, r13_binSize, norders)
	{ m_tripletType = TripletType::_multipoles_direct_; set_parameters(); m_scale.resize(m_nbins+1, 0.); m_TT1D.resize(m_nbins+1, 0.); }

      ///@}

    
      /**
       *  @name Member functions used to handle triplets
       */
      ///@{

      /**
       *  @brief estimate the distance between two objects and update
       *  the triplet vectors accordingly
       *  @param r12 distance between object1 and object2
       *  @param r13 distance between object1 and object3
       *  @param r23 distance between object2 and object3
       *  @param klin triplet bin
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      void get_triplet (const double r12, const double r13, const double r23, int &klin)    
      { (void)r12; (void)r13; (void)r23; (void)klin; ErrorCBL("", "get_triplet", "Triplet1D.h!"); }

      /**
       *  @brief update the triplet
       *  @param klin triplet bin
       *  @param ww the weight
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      void set_triplet (const int klin, const double ww=1.)   
      { (void)klin; (void)ww; cbl::ErrorCBL("", "get_triplet", "Triplet1D.h!"); }
   
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
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2, const std::shared_ptr<catalogue::Object> obj3) override;

      ///@}
    
    };


    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Triplet1D_comoving_side Triplet1D.h
     *  "Headers/Triplet1D.h"
     *
     *  @brief The class Triplet1D_comoving_side
     *
     *  This class is used to handle objects of type <EM>
     *  Triplet1D_comoving_side </EM>, used to handle 1D triplets of
     *  objects in comoving coordinates linearly binned
     */
    class Triplet1D_comoving_side : public Triplet1D_comoving {

    private:

      /// minimum scale
      double m_min;

      /// maximum scale
      double m_max;

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
      { m_tripletType = TripletType::_comoving_side_; }

      /**
       *  @brief constructor
       *  @param r12 the size of r<SUB>12</SUB>
       *  @param r12_binSize the size of r<SUB>12</SUB> bin
       *  @param r13 the size of r<SUB>13</SUB>
       *  @param r13_binSize the size of r<SUB>13</SUB> bin
       *  @param nbins the number of bins
       *  @return object of class Triplet1D_comoving_side
       */
      Triplet1D_comoving_side (const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int nbins) 
	: Triplet1D_comoving(r12, r12_binSize, r13, r13_binSize, nbins)
      { m_tripletType = TripletType::_comoving_side_; set_parameters(); m_scale.resize(m_nbins+1, 0.); m_TT1D.resize(nbins+1, 0.); }

      ///@}

      
      /**
       *  @name Member functions used to handle triplets
       */
      ///@{
     
      /**
       *  @brief estimate the distance between two objects and update
       *  the triplet vectors accordingly
       *  @param r12 distance between object1 and object2
       *  @param r13 distance between object1 and object3
       *  @param r23 distance between object2 and object3
       *  @param klin triplet bin
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      void get_triplet (const double r12, const double r13, const double r23, int &klin) override;   

      /**
       *  @brief update the triplet
       *  @param klin triplet bin
       *  @param ww the weight
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      void set_triplet (const int klin, const double ww=1.) override;

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
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2, const std::shared_ptr<catalogue::Object> obj3) override;

      ///@}
      
    };

    
    // ============================================================================================
    // ============================================================================================


   /**
     *  @class Triplet1D_comoving_costheta Triplet1D.h
     *  "Headers/Triplet1D.h"
     *
     *  @brief The class Triplet1D_comoving_costheta
     *
     *  This class is used to handle objects of type <EM>
     *  Triplet1D_comoving_costheta </EM>, used to handle 1D triplets of
     *  objects in comoving coordinates, in bins of the cosine of the angle
     */
    class Triplet1D_comoving_costheta : public Triplet1D_comoving {
      
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
       *  @return object of class Triplet1D_comoving_costheta, with protected
       *  members set to -1
       */
      Triplet1D_comoving_costheta () 
	{ m_tripletType = TripletType::_comoving_costheta_; }
                                                                
      /**
       *  @brief constructor
       *  @param r12 the size of r<SUB>12</SUB>
       *  @param r12_binSize the size of r<SUB>12</SUB> bin
       *  @param r13 the size of r<SUB>13</SUB>
       *  @param r13_binSize the size of r<SUB>13</SUB> bin
       *  @param nbins the number of bins
       *  @return object of class Triplet1D_comoving_costheta
       */
      Triplet1D_comoving_costheta (const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int nbins) 
	: Triplet1D_comoving(r12, r12_binSize, r13, r13_binSize, nbins)
	{ m_tripletType = TripletType::_comoving_costheta_; set_parameters(); m_scale.resize(m_nbins+1, 0.); m_TT1D.resize(m_nbins+1, 0.); }

      ///@}

    
      /**
       *  @name Member functions used to handle triplets
       */
      ///@{
     
      /**
       *  @brief estimate the distance between two objects and update
       *  the triplet vectors accordingly
       *  @param r12 distance between object1 and object2
       *  @param r13 distance between object1 and object3
       *  @param r23 distance between object2 and object3
       *  @param klin triplet bin
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      void get_triplet (const double r12, const double r13, const double r23, int &klin) override;   

      /**
       *  @brief update the triplet
       *  @param klin triplet bin
       *  @param ww the weight
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      void set_triplet (const int klin, const double ww=1.) override;

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
      void put (const std::shared_ptr<catalogue::Object> obj1, const std::shared_ptr<catalogue::Object> obj2, const std::shared_ptr<catalogue::Object> obj3) override;

      ///@}
    
    };

    
    // ============================================================================================
    // ============================================================================================


     /**
     *  @class Triplet1D_angular Triplet1D.h
     *  "Headers/Triplet1D.h"
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
       *  @param r12 the size of r<SUB>12</SUB>
       *  @param r12_binSize the size of r<SUB>12</SUB> bin
       *  @param r13 the size of r<SUB>13</SUB>
       *  @param r13_binSize the size of r<SUB>13</SUB> bin
       *  @param nbins the number of bins
       *  @return object of class Triplet1D_angular
       */
      Triplet1D_angular (const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int nbins) 
	: Triplet1D(r12, r12_binSize, r13, r13_binSize, nbins) {}

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
