/*******************************************************************
 *  Copyright (C) 2015 by Federico Marulli                         *
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
 *  @file Headers/Lib/Pair2D.h
 *
 *  @brief The classes Pair2D*
 *
 *  This file defines the interface of the base class Pair2D, used to
 *  handle 2D pairs of objects of any kind
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unbo.it
 */

#ifndef __PAIR2D__
#define __PAIR2D__


#include "Pair.h"


// ===================================================================================================


namespace cosmobl {
  
  namespace pairs {

    /**
     *  @class Pair2D Pair2D.h "Headers/Lib/Pair2D.h"
     *
     *  @brief The class Pair2D
     *
     *  This class is used to handle objects of type <EM> Pair2D
     *  </EM>.
     */
    class Pair2D : public Pair {

    private:
      
      /**
       *  @name Member functions used to set the binning parameters (customized in all the derived classes) 
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  @return none
       */
      virtual void m_set_parameters_nbins () = 0;
  
      /**
       *  @brief set the binning parameters given the bin size
       *  @return none
       */
      virtual void m_set_parameters_binSize () = 0;
  
      ///@}

      
    protected:

      /// the binned scales in the first dimension
      vector<double> m_scale_D1;
      
      /// the binned scales in the second dimension
      vector<double> m_scale_D2;
      
      /// the number of binned pairs 
      vector<vector<double> > m_PP2D;

      /**
       *  @name Binning parameters
       */
      ///@{
  
      /// the inverse of the bin size in the first dimension
      double m_binSize_inv_D1;
  
      /// number of bins in the first dimension
      int m_nbins_D1;
  
      /// radial shift used to centre the output bins in the first dimension
      double m_shift_D1;
  
      /// the inverse of the bin size in the second dimension
      double m_binSize_inv_D2;
  
      /// number of bins in the second dimension
      int m_nbins_D2;
  
      /// radial shift used to centre the output bins in the second dimension
      double m_shift_D2;
  
      ///@}
  
    public:
  
      /**
       *  @name Constructors/destructors
       */
      ///@{
  
      /**
       *  @brief default constuctor
       *  @return object of class Pair2D
       */
      Pair2D ()
	: m_binSize_inv_D1(1.), m_nbins_D1(50), m_shift_D1(0.5), m_binSize_inv_D2(1.), m_nbins_D2(50), m_shift_D2(0.5)
	{ m_pairDim = _2D_; m_angularUnits = _radians_; m_angularWeight = nullptr; }

      /**
       *  @brief constuctor
       *  @param binSize_D1 the bin size in the first dimension
       *  @param nbins_D1 number of bins in the first dimension
       *  @param shift_D1 radial shift used to centre the output bins
       *  in the first dimension
       *  @param binSize_D2 the bin size in the second dimension
       *  @param nbins_D2 number of bins in the second dimension
       *  @param shift_D2 radial shift used to centre the output bins
       *  in the second dimension
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D
       */
      Pair2D (const double binSize_D1, const int nbins_D1, const double shift_D1, const double binSize_D2, const int nbins_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: m_binSize_inv_D1(1./binSize_D1), m_nbins_D1(nbins_D1), m_shift_D1(shift_D1), m_binSize_inv_D2(1./binSize_D2), m_nbins_D2(nbins_D2), m_shift_D2(shift_D2)
      { m_pairDim = _2D_; m_angularUnits = angularUnits; m_angularWeight = angularWeight; }
  
      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Pair2D () = default;
  
      ///@}
  
  
      /**
       *  @name Member functions to get the protected members 
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
      double PP2D (const int i, const int j) const override { return m_PP2D[i][j]; }

      /**
       *  @brief get the protected member \e m_PP2D
       *  @return the vector containing the binned number of pairs
       */
      vector<vector<double> > PP2D () const override { return m_PP2D; }
    
      /**
       *  @brief get the protected member Pair2D::m_binSize_inv_D1
       *  @return the inverse of the bin size in the first dimension
       */
      double binSize_inv_D1 () const override { return m_binSize_inv_D1; }

      /**
       *  @brief get the protected member Pair2D::m_nbins_D1
       *  @return the number of bins in the first dimension
       */
      int nbins_D1 () const override { return m_nbins_D1; }
    
      /**
       *  @brief get the protected member Pair2D::m_shift_D1
       *  @return the radial shift in the first dimension used to centre
       *  the output bins
       */
      double shift_D1 () const override { return m_shift_D1; }
    
      /**
       *  @brief get the protected member Pair2D::m_binSize_inv_D2
       *  @return the inverse of the bin size in the second dimension
       */
      double binSize_inv_D2 () const override { return m_binSize_inv_D2; }

      /**
       *  @brief get the protected member Pair2D::m_nbins_D2
       *  @return the number of bins in the second dimension
       */
      int nbins_D2 () const override { return m_nbins_D2; }
    
      /**
       *  @brief get the protected member Pair2D::m_shift_D2
       *  @return the radial shift in the second dimension used to centre
       *  the output bins
       */
      double shift_D2 () const override { return m_shift_D2; }
  
      ///@}

    
      /**
       *  @name Member functions used to set the protected members
       */
      ///@{

      /**
       *  @brief set the protected member Pair2D::m_scale_D1[i]
       *  @param i the bin index in the first dimension
       *  @param pp the binned scales
       *  @return none
       */
      void set_scale_D1 (const int i, const double pp) { checkDim(m_scale_D1, i, "m_scale_D1"); m_scale_D1[i] = pp; }

      /**
       *  @brief set the protected member Pair2D::m_scale_D2[i]
       *  @param i the bin index in the first dimension
       *  @param pp the binned scales
       *  @return none
       */
      void set_scale_D2 (const int i, const double pp) { checkDim(m_scale_D2, i, "m_scale_D2"); m_scale_D2[i] = pp; }
      
      /**
       *  @brief set the protected member Pair2D::m_PP2D[i][j]
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @param pp the number of pairs in the bin
       *  @return none
       */
      void set_PP2D (const int i, const int j, const double pp) { checkDim(m_PP2D, i, j, "m_PP2D"); m_PP2D[i][j] = pp; }

      /**
       *  @brief set the protected member PairD1::m_PP2D[i][j] adding
       *  the number of pairs
       *  @param i the bin index in the first dimension
       *  @param j the bin index in the second dimension
       *  @param pp the number of pairs in the bin
       *  @return none
       */
      void add_PP2D (const int i, const int j, const double pp) { checkDim(m_PP2D, i, j, "m_PP2D"); m_PP2D[i][j] += pp; }

      ///@}


      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
    
      /**
       *  @brief sum the number of binned pairs
       *  @param pp an object of class Pair
       *  @param ww the weight
       *  @return none
       */
      void Sum (const shared_ptr<Pair> pp, const double ww=1) override;

      ///@}
    
    };


    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingCartesian Pair2D.h
     *  "Headers/Lib/Pair2D.h"
     *
     *  @brief The class Pair2D_comovingCartesian
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingCartesian </EM>.
     */
    class Pair2D_comovingCartesian : public Pair2D {
      
    private:
      
      /**
       *  @name Member functions used to set the binning parameters (customized in all the derived classes) 
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  @return none
       */
      virtual void m_set_parameters_nbins () = 0;
  
      /**
       *  @brief set the binning parameters given the bin size
       *  @return none
       */
      virtual void m_set_parameters_binSize () = 0;
  
      ///@}
      
      
    protected:
  
      /**
       *  @name Binning parameters
       */
      ///@{

      /// minimum perpendicular separation used to count the pairs
      double m_rpMin;
  
      /// maximum perpendicular separation used to count the pairs
      double m_rpMax;

      /// minimum parallel separation used to count the pairs
      double m_piMin;
  
      /// maximum parallel separation used to count the pairs
      double m_piMax;

      ///@}


    public:
  
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Pair2D_comovingCartesian
       */
      Pair2D_comovingCartesian () : m_rpMin(0.1), m_rpMax(50.), m_piMin(0.1), m_piMax(50.) {}

      /**
       *  @brief constuctor   
       *  @param rpMin minimum perpendicular separation used to count
       *  the pairs
       *  @param rpMax maximum perpendicular separation used to count
       *  the pairs
       *  @param nbins_rp number of bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular
       *  separation, i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param nbins_pi number of bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian
       */
      Pair2D_comovingCartesian (const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D(1., nbins_rp, shift_rp, 1., nbins_pi, shift_pi, angularUnits, angularWeight), m_rpMin(rpMin), m_rpMax(rpMax), m_piMin(piMin), m_piMax(piMax) {}
  
      /**
       *  @brief constuctor
       *  @param rpMin minimum perpendicular separation used to count
       *  the pairs
       *  @param rpMax maximum perpendicular separation used to count
       *  the pairs
       *  @param binSize_rp size of the bins in the perpendicular
       *  separation
       *  @param shift_rp shift parameter in the perpendicular
       *  separation, i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param binSize_pi size of the bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian
       */
      Pair2D_comovingCartesian (const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D(binSize_rp, 50, shift_rp, binSize_pi, 50, shift_pi, angularUnits, angularWeight), m_rpMin(rpMin), m_rpMax(rpMax), m_piMin(piMin), m_piMax(piMax) {}
  
      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Pair2D_comovingCartesian () = default;

      ///@}

  
      /**
       *  @name Member functions used to get the protected parameters
       */
      ///@{

      /**
       *  @brief get the protected member Pair2D::m_rpMin
       *  @return the minimum perpendicular separation used to count the
       *  pairs
       */
      double sMin_D1 () const override { return m_rpMin; }
      
      /**
       *  @brief get the protected member Pair2D::m_rpMax
       *  @return the maximum perpendicular separation used to count the
       *  pairs
       */
      double sMax_D1 () const override { return m_rpMax; }

      /**
       *  @brief get the protected member Pair2D::m_piMin
       *  @return the minimum parallel separation used to count the
       *  pairs
       */
      double sMin_D2 () const override { return m_piMin; }

      /**
       *  @brief get the protected member Pair2D::m_piMax
       *  @return the maximum parallel separation used to count the
       *  pairs
       */
      double sMax_D2 () const override { return m_piMax; }

      ///@}
    
    };


    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingCartesian_linlin Pair2D.h
     *  "Headers/Lib/Pair2D.h"
     *
     *  @brief The class Pair2D_comovingCartesian_linlin
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingCartesian_linlin </EM>.
     */
    class Pair2D_comovingCartesian_linlin : public Pair2D_comovingCartesian {

    private:
    
      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  @return none
       */
      void m_set_parameters_nbins () override;
    
      /**
       *  @brief set the binning parameters given the bin size
       *  @return none
       */
      void m_set_parameters_binSize () override;

      ///@}
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Pair2D_comovingCartesian_linlin
       */
      Pair2D_comovingCartesian_linlin () { m_pairType = _comovingCartesian_linlin_; } 

      /**
       *  @brief constuctor   
       *  @param rpMin minimum perpendicular separation used to count the
       *  pairs
       *  @param rpMax maximum perpendicular separation used to count the
       *  pairs
       *  @param nbins_rp number of bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular separation,
       *  i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param nbins_pi number of bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian
       */
      Pair2D_comovingCartesian_linlin (const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingCartesian(rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight)
	{ m_pairType = _comovingCartesian_linlin_; m_set_parameters_nbins(); m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); }
  
      /**
       *  @brief constuctor
       *  @param rpMin minimum perpendicular separation used to count the
       *  pairs
       *  @param rpMax maximum perpendicular separation used to count the
       *  pairs
       *  @param binSize_rp size of the bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular separation,
       *  i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param binSize_pi size of the bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian
       */
      Pair2D_comovingCartesian_linlin (const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingCartesian(rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight)
	{ m_pairType = _comovingCartesian_linlin_; m_set_parameters_binSize(); m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); } 
  
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Pair2D_comovingCartesian_linlin () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      void put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2) override;
  
      ///@}

    };


    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingCartesian_loglin Pair2D.h
     *  "Headers/Lib/Pair2D.h"
     *
     *  @brief The class Pair2D_comovingCartesian_loglin
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingCartesian_loglin </EM>.
     */
    class Pair2D_comovingCartesian_loglin : public Pair2D_comovingCartesian {

    private:
    
      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  @return none
       */
      void m_set_parameters_nbins () override;
    
      /**
       *  @brief set the binning parameters given the bin size
       *  @return none
       */
      void m_set_parameters_binSize () override;

      ///@}
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Pair2D_comovingCartesian_loglin
       */
      Pair2D_comovingCartesian_loglin () { m_pairType = _comovingCartesian_loglin_; m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); }

      /**
       *  @brief constuctor   
       *  @param rpMin minimum perpendicular separation used to count the
       *  pairs
       *  @param rpMax maximum perpendicular separation used to count the
       *  pairs
       *  @param nbins_rp number of bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular separation,
       *  i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param nbins_pi number of bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian
       */
      Pair2D_comovingCartesian_loglin (const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingCartesian(rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight)
	{ m_pairType = _comovingCartesian_loglin_; m_set_parameters_nbins(); m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); }
  
      /**
       *  @brief constuctor
       *  @param rpMin minimum perpendicular separation used to count the
       *  pairs
       *  @param rpMax maximum perpendicular separation used to count the
       *  pairs
       *  @param binSize_rp size of the bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular separation,
       *  i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param binSize_pi size of the bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian
       */
      Pair2D_comovingCartesian_loglin (const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingCartesian(rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight)
	{ m_pairType = _comovingCartesian_loglin_; m_set_parameters_binSize(); m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); } 
  
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Pair2D_comovingCartesian_loglin () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      void put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2) override;
  
      ///@}

    };

    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingCartesian_linlog Pair2D.h
     *  "Headers/Lib/Pair2D.h"
     *
     *  @brief The class Pair2D_comovingCartesian_linlog
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingCartesian_linlog </EM>.
     */
    class Pair2D_comovingCartesian_linlog : public Pair2D_comovingCartesian {

    private:
    
      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  @return none
       */
      void m_set_parameters_nbins () override;
    
      /**
       *  @brief set the binning parameters given the bin size
       *  @return none
       */
      void m_set_parameters_binSize () override;

      ///@}
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Pair2D_comovingCartesian_linlog
       */
      Pair2D_comovingCartesian_linlog () { m_pairType = _comovingCartesian_linlog_; m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); }

      /**
       *  @brief constuctor   
       *  @param rpMin minimum perpendicular separation used to count the
       *  pairs
       *  @param rpMax maximum perpendicular separation used to count the
       *  pairs
       *  @param nbins_rp number of bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular separation,
       *  i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param nbins_pi number of bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian
       */
      Pair2D_comovingCartesian_linlog (const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingCartesian(rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight)
	{ m_pairType = _comovingCartesian_linlog_; m_set_parameters_nbins(); m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); }
  
      /**
       *  @brief constuctor
       *  @param rpMin minimum perpendicular separation used to count the
       *  pairs
       *  @param rpMax maximum perpendicular separation used to count the
       *  pairs
       *  @param binSize_rp size of the bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular separation,
       *  i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param binSize_pi size of the bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian
       */
      Pair2D_comovingCartesian_linlog (const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingCartesian(rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight)
	{ m_pairType = _comovingCartesian_linlog_; m_set_parameters_binSize(); m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); } 
  
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Pair2D_comovingCartesian_linlog () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      void put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2) override;
  
      ///@}

    };

    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingCartesian_loglog Pair2D.h
     *  "Headers/Lib/Pair2D.h"
     *
     *  @brief The class Pair2D_comovingCartesian_loglog
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingCartesian_loglog </EM>.
     */
    class Pair2D_comovingCartesian_loglog : public Pair2D_comovingCartesian {

    private:
    
      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  @return none
       */
      void m_set_parameters_nbins () override;
    
      /**
       *  @brief set the binning parameters given the bin size
       *  @return none
       */
      void m_set_parameters_binSize () override;

      ///@}
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Pair2D_comovingCartesian_loglog
       */
      Pair2D_comovingCartesian_loglog () { m_pairType = _comovingCartesian_loglog_; m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); }

      /**
       *  @brief constuctor   
       *  @param rpMin minimum perpendicular separation used to count the
       *  pairs
       *  @param rpMax maximum perpendicular separation used to count the
       *  pairs
       *  @param nbins_rp number of bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular separation,
       *  i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param nbins_pi number of bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian
       */
      Pair2D_comovingCartesian_loglog (const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingCartesian(rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight)
	{ m_pairType = _comovingCartesian_loglog_; m_set_parameters_nbins(); m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); }
  
      /**
       *  @brief constuctor
       *  @param rpMin minimum perpendicular separation used to count the
       *  pairs
       *  @param rpMax maximum perpendicular separation used to count the
       *  pairs
       *  @param binSize_rp size of the bins in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular separation,
       *  i.e. the shift is binSize*shift
       *  @param piMin minimum parallel separation used to count the pairs
       *  @param piMax maximum parallel separation used to count the pairs
       *  @param binSize_pi size of the bins in the parallel separation
       *  @param shift_pi shift parameter in the parallel separation,
       *  i.e. the shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingCartesian
       */
      Pair2D_comovingCartesian_loglog (const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingCartesian(rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight)
	{ m_pairType = _comovingCartesian_loglog_; m_set_parameters_binSize(); m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); } 
  
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Pair2D_comovingCartesian_loglog () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      void put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2) override;
  
      ///@}

    };

  
    // ============================================================================================
    // ============================================================================================

  
    /**
     *  @class Pair2D_comovingPolar Pair2D.h "Headers/Lib/Pair2D.h"
     *
     *  @brief The class Pair2D_comovingPolar
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingPolar </EM>.
     */
    class Pair2D_comovingPolar : public Pair2D {

    private:
      
      /**
       *  @name Member functions used to set the binning parameters (customized in all the derived classes) 
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  @return none
       */
      virtual void m_set_parameters_nbins () = 0;
  
      /**
       *  @brief set the binning parameters given the bin size
       *  @return none
       */
      virtual void m_set_parameters_binSize () = 0;
  
      ///@}

      
    protected:
  
      /**
       *  @name Binning parameters
       */
      ///@{

      /// minimum separation used to count the pairs
      double m_rMin;
  
      /// maximum separation used to count the pairs
      double m_rMax;

      /// minimum angle used to count the pairs
      double m_muMin;
  
      /// maximum angle used to count the pairs
      double m_muMax;

      ///@}


    public:
  
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Pair2D_comovingPolar
       */
      Pair2D_comovingPolar () : m_rMin(0.1), m_rMax(50.), m_muMin(0.1), m_muMax(50.) {}

      /**
       *  @brief constuctor   
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins_D1 number of bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param nbins_D2 number of bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar
       */
      Pair2D_comovingPolar (const double rMin, const double rMax, const int nbins_D1, const double shift_D1, const double muMin, const double muMax, const int nbins_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D(1., nbins_D1, shift_D1, 1., nbins_D2, shift_D2, angularUnits, angularWeight), m_rMin(rMin), m_rMax(rMax), m_muMin(muMin), m_muMax(muMax) {}
  
      /**
       *  @brief constuctor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize_D1 size of the bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param binSize_D2 size of the bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar
       */
      Pair2D_comovingPolar (const double rMin, const double rMax, const double binSize_D1, const double shift_D1, const double muMin, const double muMax, const double binSize_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D(binSize_D1, 50, shift_D1, binSize_D2, 50, shift_D2, angularUnits, angularWeight), m_rMin(rMin), m_rMax(rMax), m_muMin(muMin), m_muMax(muMax) {}
  
      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Pair2D_comovingPolar () = default;

      ///@}

  
      /**
       *  @name Member functions used to get the protected parameters
       */
      ///@{

      /**
       *  @brief get the protected member Pair2D::m_rMin
       *  @return the minimum separation used to count the pairs
       */
      double sMin_D1 () const override { return m_rMin; }

      /**
       *  @brief get the protected member Pair2D::m_rMax
       *  @return the maximum separation used to count the pairs
       */
      double sMax_D1 () const override { return m_rMax; }

      /**
       *  @brief get the protected member Pair2D::m_muMin
       *  @return the minimum angle used to count the pairs
       */
      double sMin_D2 () const override { return m_muMin; }

      /**
       *  @brief get the protected member Pair2D::m_muMax
       *  @return the maximum angle used to count the pairs
       */
      double sMax_D2 () const override { return m_muMax; }

      ///@}
    
    };

  
    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingPolar_linlin Pair2D.h
     *  "Headers/Lib/Pair2D.h"
     *
     *  @brief The class Pair2D_comovingPolar_linlin
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingPolar_linlin </EM>.
     */
    class Pair2D_comovingPolar_linlin : public Pair2D_comovingPolar {

    private:
    
      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  @return none
       */
      void m_set_parameters_nbins () override;
    
      /**
       *  @brief set the binning parameters given the bin size
       *  @return none
       */
      void m_set_parameters_binSize () override;

      ///@}
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Pair2D_comovingPolar_linlin
       */
      Pair2D_comovingPolar_linlin () { m_pairType = _comovingPolar_linlin_; } 

      /**
       *  @brief constuctor   
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins_D1 number of bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param nbins_D2 number of bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar
       */
      Pair2D_comovingPolar_linlin (const double rMin, const double rMax, const int nbins_D1, const double shift_D1, const double muMin, const double muMax, const int nbins_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingPolar(rMin, rMax, nbins_D1, shift_D1, muMin, muMax, nbins_D2, shift_D2, angularUnits, angularWeight)
	{ m_pairType = _comovingPolar_linlin_; m_set_parameters_nbins(); m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); }
  
      /**
       *  @brief constuctor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize_D1 size of the bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param binSize_D2 size of the bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar
       */
      Pair2D_comovingPolar_linlin (const double rMin, const double rMax, const double binSize_D1, const double shift_D1, const double muMin, const double muMax, const double binSize_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingPolar(rMin, rMax, binSize_D1, shift_D1, muMin, muMax, binSize_D2, shift_D2, angularUnits, angularWeight)
	{ m_pairType = _comovingPolar_linlin_; m_set_parameters_binSize(); m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); } 
  
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Pair2D_comovingPolar_linlin () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      void put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2) override;
  
      ///@}

    };

  
    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingPolar_loglin Pair2D.h
     *  "Headers/Lib/Pair2D.h"
     *
     *  @brief The class Pair2D_comovingPolar_loglin
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingPolar_loglin </EM>.
     */
    class Pair2D_comovingPolar_loglin : public Pair2D_comovingPolar {

    private:
    
      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  @return none
       */
      void m_set_parameters_nbins () override;
    
      /**
       *  @brief set the binning parameters given the bin size
       *  @return none
       */
      void m_set_parameters_binSize () override;

      ///@}
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Pair2D_comovingPolar_loglin
       */
      Pair2D_comovingPolar_loglin () { m_pairType = _comovingPolar_loglin_; m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); } 

      /**
       *  @brief constuctor   
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins_D1 number of bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param nbins_D2 number of bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar
       */
      Pair2D_comovingPolar_loglin (const double rMin, const double rMax, const int nbins_D1, const double shift_D1, const double muMin, const double muMax, const int nbins_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingPolar(rMin, rMax, nbins_D1, shift_D1, muMin, muMax, nbins_D2, shift_D2, angularUnits, angularWeight)
	{ m_pairType = _comovingPolar_loglin_; m_set_parameters_nbins(); m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); }
  
      /**
       *  @brief constuctor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize_D1 size of the bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param binSize_D2 size of the bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar
       */
      Pair2D_comovingPolar_loglin (const double rMin, const double rMax, const double binSize_D1, const double shift_D1, const double muMin, const double muMax, const double binSize_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingPolar(rMin, rMax, binSize_D1, shift_D1, muMin, muMax, binSize_D2, shift_D2, angularUnits, angularWeight)
	{ m_pairType = _comovingPolar_loglin_; m_set_parameters_binSize(); m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); } 
  
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Pair2D_comovingPolar_loglin () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      void put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2) override;
  
      ///@}

    };

  
    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingPolar_linlog Pair2D.h
     *  "Headers/Lib/Pair2D.h"
     *
     *  @brief The class Pair2D_comovingPolar_linlog
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingPolar_linlog </EM>.
     */
    class Pair2D_comovingPolar_linlog : public Pair2D_comovingPolar {

    private:
    
      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  @return none
       */
      void m_set_parameters_nbins () override;
    
      /**
       *  @brief set the binning parameters given the bin size
       *  @return none
       */
      void m_set_parameters_binSize () override;

      ///@}
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Pair2D_comovingPolar_linlog
       */
      Pair2D_comovingPolar_linlog () { m_pairType = _comovingPolar_linlog_; m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); } 

      /**
       *  @brief constuctor   
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins_D1 number of bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param nbins_D2 number of bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar
       */
      Pair2D_comovingPolar_linlog (const double rMin, const double rMax, const int nbins_D1, const double shift_D1, const double muMin, const double muMax, const int nbins_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingPolar(rMin, rMax, nbins_D1, shift_D1, muMin, muMax, nbins_D2, shift_D2, angularUnits, angularWeight)
	{ m_pairType = _comovingPolar_linlog_; m_set_parameters_nbins(); m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); }
  
      /**
       *  @brief constuctor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize_D1 size of the bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param binSize_D2 size of the bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar
       */
      Pair2D_comovingPolar_linlog (const double rMin, const double rMax, const double binSize_D1, const double shift_D1, const double muMin, const double muMax, const double binSize_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingPolar(rMin, rMax, binSize_D1, shift_D1, muMin, muMax, binSize_D2, shift_D2, angularUnits, angularWeight)
	{ m_pairType = _comovingPolar_linlog_; m_set_parameters_binSize(); m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); } 
  
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Pair2D_comovingPolar_linlog () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      void put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2) override;
  
      ///@}

    };

    // ============================================================================================
    // ============================================================================================


    /**
     *  @class Pair2D_comovingPolar_loglog Pair2D.h
     *  "Headers/Lib/Pair2D.h"
     *
     *  @brief The class Pair2D_comovingPolar_loglog
     *
     *  This class is used to handle objects of type <EM>
     *  Pair2D_comovingPolar_loglog </EM>.
     */
    class Pair2D_comovingPolar_loglog : public Pair2D_comovingPolar {

    private:
    
      /**
       *  @name Member functions used to set the binning parameters
       */
      ///@{
  
      /**
       *  @brief set the binning parameters given the number of bins
       *  @return none
       */
      void m_set_parameters_nbins () override;
    
      /**
       *  @brief set the binning parameters given the bin size
       *  @return none
       */
      void m_set_parameters_binSize () override;

      ///@}
  
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Pair2D_comovingPolar_loglog
       */
      Pair2D_comovingPolar_loglog () { m_pairType = _comovingPolar_loglog_; m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); } 

      /**
       *  @brief constuctor   
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins_D1 number of bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param nbins_D2 number of bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar
       */
      Pair2D_comovingPolar_loglog (const double rMin, const double rMax, const int nbins_D1, const double shift_D1, const double muMin, const double muMax, const int nbins_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingPolar(rMin, rMax, nbins_D1, shift_D1, muMin, muMax, nbins_D2, shift_D2, angularUnits, angularWeight)
	{ m_pairType = _comovingPolar_loglog_; m_set_parameters_nbins(); m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); }
  
      /**
       *  @brief constuctor
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize_D1 size of the bins in the first dimension
       *  @param shift_D1 shift parameter in the first dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param muMin minimum angle used to count the pairs
       *  @param muMax maximum angle used to count the pairs
       *  @param binSize_D2 size of the bins in the second dimension
       *  @param shift_D2 shift parameter in the second dimension,
       *  i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class Pair2D_comovingPolar
       */
      Pair2D_comovingPolar_loglog (const double rMin, const double rMax, const double binSize_D1, const double shift_D1, const double muMin, const double muMax, const double binSize_D2, const double shift_D2, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: Pair2D_comovingPolar(rMin, rMax, binSize_D1, shift_D1, muMin, muMax, binSize_D2, shift_D2, angularUnits, angularWeight)
	{ m_pairType = _comovingPolar_loglog_; m_set_parameters_binSize(); m_PP2D.resize(m_nbins_D1+1, vector<double>(m_nbins_D2+1, 0.)); } 
  
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Pair2D_comovingPolar_loglog () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to handle pairs
       */
      ///@{
  
      /**
       *  @brief estimate the distance between two objects and update the
       *  pair vector accordingly
       *  @param obj1 pointer to an object of class Object
       *  @param obj2 pointer to an object of class Object
       *  @return none
       */
      void put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2) override;
  
      ///@}

    };  
  }
}

#endif
