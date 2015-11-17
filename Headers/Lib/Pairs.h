/*******************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo *
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
 *  @file Headers/Lib/Pairs.h
 *
 *  @brief The class Pairs
 *
 *  This file defines the interface of the class Pairs, used to handle
 *  pairs of objects to compute the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __PAIRS__
#define __PAIRS__


#include "Catalogue.h"


// ===================================================================================================


namespace cosmobl {
  
  /**
   *  @class Pairs Pairs.h "Headers/Lib/Pairs.h"
   *
   *  @brief The class Pairs
   *
   *  This class is used to handle objects of type <EM> Pairs
   *  </EM>. It contains all virtual methods implemented in the
   *  derived classes Pairs2D and Pairs3D
   */
  class Pairs {

  public:
    
    /**
     *  @brief default destructor
     *  @return none
     */
    virtual ~Pairs () {}

    /**
     *  @brief get the private member \e m_nlog
     *  @return the number of logarithmic bins, or an error message if
     *  the derived object does not have this member
     */
    virtual int nlog () const { cosmobl::ErrorMsg("Error in Pairs::nlog() of Pairs.h!"); return 0; }

    /**
     *  @brief get the private member \e m_nlin
     *  @return the number of linear bins, or an error message if the
     *  derived object does not have this member
     */
    virtual int nlin () const { cosmobl::ErrorMsg("Error in Pairs::nlin() of Pairs.h!"); return 0; }

    /**
     *  @brief get the private member \e m_ncos
     *  @return the number of angular bins, or an error message if the
     *  derived object does not have this member
     */
    virtual int ncos () const { cosmobl::ErrorMsg("Error in Pairs::ncos() of Pairs.h!"); return 0; }

     /**
     *  @brief get the protected member \e m_thetaMin
     *  @return the minimum value of the angle &theta; used to count
     *  the number of pairs
     */
    virtual double thetaMin () const { cosmobl::ErrorMsg("Error in Pairs::thetaMin() of Pairs.h!"); return 0; }

    /**
     *  @brief get the protected member \e m_thetaMax
     *  @return the maximum value of the angle &theta; used to count
     *  the number of pairs
     */    
    virtual double thetaMax () const { cosmobl::ErrorMsg("Error in Pairs::thetaMax() of Pairs.h!"); return 0; }
    
    /**
     *  @brief get the private member \e m_rMin
     *  @return the minimum separation used to count pairs, or an
     *  error message if the derived object does not have this member
     */
    virtual double rMin () const { cosmobl::ErrorMsg("Error in Pairs::rMin() of Pairs.h!"); return 0; }

    /**
     *  @brief get the private member \e m_rMax
     *  @return the maximum separation used to count pairs, or an
     *  error message if the derived object does not have this member
     */
    virtual double rMax () const { cosmobl::ErrorMsg("Error in Pairs::rMax() of Pairs.h!"); return 0; }

    /**
     *  @brief get the private member \e m_PPlog[i]
     *  @param i the bin index
     *  @return the number of pairs in the i-th logarithmic bin, or an
     *  error message if the derived object does not have this member
     */
    virtual double PPlog (const int i) const { cosmobl::ErrorMsg("Error in Pairs::PPlog() of Pairs.h!"); return 0; }

    /**
     *  @brief get the private member \e m_PPlin[i]
     *  @param i the bin index
     *  @return the number of pairs in the i-th linear bin, or an
     *  error message if the derived object does not have this member
     */
    virtual double PPlin (const int i) const { cosmobl::ErrorMsg("Error in Pairs::PPlin() of Pairs.h!"); return 0; }

    /**
     *  @brief get the private member \e m_PPlog
     *  @return the vector containing the number of pairs in
     *  logarithmic bins, or an error message if the derived object
     *  does not have this member
     */
    virtual vector<double> PPlog () const { cosmobl::ErrorMsg("Error in Pairs::PPlog() of Pairs.h!"); vector<double> PP; return PP; }

    /**
     *  @brief get the private member \e m_PPlin
     *  @return the vector containing the number of pairs in linear
     *  bins, or an error message if the derived object does not have
     *  this member
     */
    virtual vector<double> PPlin () const { cosmobl::ErrorMsg("Error in Pairs::PPlin() of Pairs.h!"); vector<double> PP; return PP; }
    
    /**
     *  @brief get the private member \e m_PP2d[i][j]
     *  @param i the bin index in the direction perpendicular to the
     *  line-of-sight
     *  @param j the bin index in the direction parallel to the
     *  line-of-sight
     *  @return the number of pairs in the i-th linear bin
     *  perpendicular to the line-of-sight, and in the j-th linear bin
     *  parallel to the line-of-sight, or an error message if the
     *  derived object does not have this member
     */
    virtual double PP2d (const int i, const int j) const { cosmobl::ErrorMsg("Error in Pairs::PP2d() of Pairs.h!"); return 0; }

    /**
     *  @brief get the private member \e m_PPslog[i][j]
     *  @param i the logarithmic bin index in the direction
     *  perpendicular to the line-of-sight
     *  @param j the linear bin index in the direction parallel to the
     *  line-of-sight
     *  @return the number of pairs in the i-th linear bin
     *  perpendicular to the line-of-sight, and in the j-th
     *  logarithmic bin parallel to the line-of-sight, or an error
     *  message if the derived object does not have this member
     */
    virtual double PPslog (const int i, const int j) const { cosmobl::ErrorMsg("Error in Pairs::PPslog() of Pairs.h!"); return 0; }     

    /**
     *  @brief get the private member \e m_PPcoslog[i][j]
     *  @param i the logarithmic bin index 
     *  @param j the angular bin index
     *  @return the number of pairs in the i-th logarithmic bin in r,
     *  and in the j-th linear angular bin in cos(&theta;), or an
     *  error message if the derived object does not have this member
     */
    virtual double PPcoslog (const int i, const int j) const { cosmobl::ErrorMsg("Error in Pairs::PPcoslog() of Pairs.h!"); return 0; }

    /**
     *  @brief get the private member \e m_PPcoslin[i][j]
     *  @param i the linear bin index 
     *  @param j the angular bin index
     *  @return the number of pairs in the i-th linear bin in r, and
     *  in the j-th linear angular bin in cos(&theta;), or an error
     *  message if the derived object does not have this member
     */
    virtual double PPcoslin (const int i, const int j) const { cosmobl::ErrorMsg("Error in Pairs::PPcoslin() of Pairs.h!"); return 0; }
   
    /**
     *  @brief get the private member \e m_PP2d
     *  @return a matrix containing the number of pairs in linear bins
     *  perpendicular to the line-of-sight, and linear bins parallel
     *  to the line-of-sight, or an error message if the derived
     *  object does not have this member
     */
    virtual vector<vector<double> > PP2d () const { cosmobl::ErrorMsg("Error in Pairs::PP2d() of Pairs.h!"); vector<vector<double> > PP; return PP; }
   
    /**
     *  @brief get the private member \e m_PPslog
     *  @return a matrix containing the number of pairs in linear bins
     *  perpendicular to the line-of-sight, and in logarithmic bins
     *  parallel to the line-of-sight, or an error message if the
     *  derived object does not have this member
     */
    virtual vector<vector<double> > PPslog () const { cosmobl::ErrorMsg("Error in Pairs::PPslog() of Pairs.h!"); vector<vector<double> > PP; return PP; }

    /**
     *  @brief get the private member \e m_PPcoslog
     *  @return a matrix containing the number of pairs in logarithmic
     *  bins in r and linear angular bins in cos(&theta;), or an error
     *  message if the derived object does not have this member
     */
    virtual vector<vector<double> > PPcoslog () const { cosmobl::ErrorMsg("Error in Pairs::PPcoslog() of Pairs.h!"); vector<vector<double> > PP; return PP; }

    /**
     *  @brief get the private member \e m_PPcoslin
     *  @return a matrix containing the number of pairs in linear bins
     *  in r and linear angular bins in cos(&theta;), or an error
     *  message if the derived object does not have this member
     */
    virtual vector<vector<double> > PPcoslin () const { cosmobl::ErrorMsg("Error in Pairs::PPcoslin() of Pairs.h!"); vector<vector<double> > PP; return PP; }

    /**
     *  @brief set the private member \e m_PPlog[i]
     *  @param i the bin index
     *  @param pp the number of pairs in the bin
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_PPlog (const int i, const double pp) { cosmobl::ErrorMsg("Error in Pairs::PPlog() of Pairs.h!"); }

    /**
     *  @brief set the private member \e m_PPlin[i]
     *  @param i the bin index
     *  @param pp the number of pairs in the bin
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_PPlin (const int i, const double pp) { cosmobl::ErrorMsg("Error in Pairs::PPlin() of Pairs.h!"); }

    /**
     *  @brief set the private member \e m_PP2d[i][j]
     *  @param i the bin index in the direction perpendicular to the
     *  line-of-sight
     *  @param j the bin index in the direction parallel to the
     *  line-of-sight
     *  @param pp the number of pairs in the bin 
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_PP2d (const int i, const int j, const double pp) { cosmobl::ErrorMsg("Error in Pairs::PP2d() of Pairs.h!"); }

    /**
     *  @brief get the private member \e m_PPslog[i][j]
     *  @param i the logarithmic bin index in the direction
     *  perpendicular to the line-of-sight
     *  @param j the linear bin index in the direction parallel to the
     *  line-of-sight
     *  @param pp the number of pairs in the bin    
     *  @return none, or an error message if the derived object does not
     *  have this member
     */
    virtual void set_PPslog (const int i, const int j, const double pp) { cosmobl::ErrorMsg("Error in Pairs::PPslog() of Pairs.h!"); }     

    /**
     *  @brief get the private member \e m_PPcoslog[i][j]
     *  @param i the logarithmic bin index 
     *  @param j the angular bin index
     *  @param pp the number of pairs in the bin 
     *  @return none, or an error message if the derived object does not
     *  have this member
     */
    virtual void set_PPcoslog (const int i, const int j, const double pp) { cosmobl::ErrorMsg("Error in Pairs::PPcoslog() of Pairs.h!"); }

    /**
     *  @brief get the private member \e m_PPcoslin[i][j]
     *  @param i the linear bin index 
     *  @param j the angular bin index
     *  @param pp the number of pairs in the bin 
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_PPcoslin (const int i, const int j, const double pp) { cosmobl::ErrorMsg("Error in Pairs::PPcoslin() of Pairs.h!"); }
    
    /**
     *  @brief sum the number of pairs
     *  @param pp pointer to an object of class Pairs
     *  @param ww the weight
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void sum (const shared_ptr<Pairs> pp, const double ww) { cosmobl::ErrorMsg("Error in Pairs::sum() of Pairs.h!"); }

    /**
     *  @brief estimate the distance between two objects and update
     *	the pair vectors accordingly
     *  @param obj1 pointer to an object of class Object
     *  @param obj2 pointer to an object of class Object
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) { cosmobl::ErrorMsg("Error in Pairs::put() of Pairs.h!"); }

    /**
     *  @brief estimate the distance between two objects and update
     *	the logarithmic binned pair vector accordingly
     *  @param obj1 pointer to an object of class Object
     *  @param obj2 pointer to an object of class Object
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void put_log (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) { cosmobl::ErrorMsg("Error in Pairs::put_log() of Pairs.h!"); }

  };


  // ============================================================================

  /**
   *  @class Pairs2D Pairs.h "Headers/Lib/Pairs.h"
   *
   *  @brief The class Pairs2D
   *
   *  This class is used to handle objects of type <EM> Pairs2D
   *  </EM>, used to measure the angular two-point correlation
   *  function
   */
  class Pairs2D : public Pairs
  {

  protected:

    /// the number of logarithmic bins
    int m_nlog;

    /// the number of linear bins
    int m_nlin;

    /// the minimum value of the angle &theta; used to count the number of pairs
    double m_thetaMin; 

    /// the maximum value of the angle &theta; used to count the number of pairs
    double m_thetaMax;

    /// the inverse of the logarithmic bin size
    double m_logbinsz_inv;
    
    /// the inverse of the linear bin size
    double m_linbinsz_inv;

    /// the number of pairs in logarithmic bins of &theta;
    vector<double> m_PPlog;

    /// the number of pairs in linear bins of &theta;
    vector<double> m_PPlin;

  
  public:
  
    /**
     *  @brief default constructor
     *  @return object of class Pairs2D, with protected members set to -1
     */
    Pairs2D () 
      : m_nlog(-1), m_nlin(-1), m_thetaMin(-1.), m_thetaMax(-1.), m_logbinsz_inv(-1.), m_linbinsz_inv(-1.) {}
    
    /**
     *  @brief constructor
     *  @param nlog the number of logarithmic bins
     *  @param nlin the number of linear bins
     *  @param thetaMin the minimum value of the angle &theta; used to
     *  count the number of pairs
     *  @param thetaMax the maximum value of the angle &theta; used to
     *  count the number of pairs
     *  @param logbinsz the logarithmic bin size
     *  @param linbinsz the linear bin size
     *  @return object of class Pairs2D
     */
    Pairs2D (const int nlog, const int nlin, const double thetaMin, const double thetaMax, const double logbinsz, const double linbinsz) 
      : m_nlog(nlog), m_nlin(nlin), m_thetaMin(thetaMin), m_thetaMax(thetaMax), m_logbinsz_inv(1./logbinsz), m_linbinsz_inv(1./linbinsz)
      { 
	m_PPlog.resize(m_nlog+1, 0.);
	m_PPlin.resize(m_nlin+1, 0.);
      }
  
    /**
     *  @brief get the protected member Pairs2D::m_nlog
     *  @return the number of logarithmic bins
     */
    int nlog () const /*override*/ { return m_nlog; }

    /**
     *  @brief get the protected member Pairs2D::m_nlin
     *  @return the number of linear bins
     */
    int nlin () const /*override*/ { return m_nlin; }

    /**
     *  @brief get the protected member Pairs2D::m_thetaMin
     *  @return the minimum value of the angle &theta; used to count
     *  the number of pairs
     */
    double thetaMin () const /*override*/ { return m_thetaMin; }

    /**
     *  @brief get the protected member Pairs2D::m_thetaMax
     *  @return the maximum value of the angle &theta; used to count
     *  the number of pairs
     */    
    double thetaMax () const /*override*/ { return m_thetaMax; }

    /**
     *  @brief get the private member Pairs2D::m_PPlog[i]
     *  @param i the bin index
     *  @return the number of pairs in the i-th logarithmic bin
     */
    double PPlog (const int i) const /*override*/ { return m_PPlog[i]; }

    /**
     *  @brief get the private member Pairs2D::m_PPlin[i]
     *  @param i the bin index
     *  @return the number of pairs in the i-th linear bin
     */
    double PPlin (const int i) const /*override*/ { return m_PPlin[i]; }

    /**
     *  @brief get the private member Pairs2D::m_PPlog
     *  @return the vector containing the number of pairs in
     *  logarithmic bins
     */
    vector<double> PPlog () const /*override*/ { return m_PPlog; }

    /**
     *  @brief get the private member Pairs2D::m_PPlin
     *  @return the vector containing the number of pairs in linear
     *  bins
     */
    vector<double> PPlin () const /*override*/ { return m_PPlin; }

    /**
     *  @brief set the private member Pairs2D::m_PPlog[i]
     *  @param i the bin index
     *  @param pp the number of pairs in the bin
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    void set_PPlog (const int i, const double pp) /*override*/ { m_PPlog[i] = pp; }

    /**
     *  @brief set the private member Pairs2D::m_PPlin[i]
     *  @param i the bin index
     *  @param pp the number of pairs in the bin
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    void set_PPlin (const int i, const double pp) /*override*/ { m_PPlin[i] = pp; }
    
    /**
     *  @brief sum the number of pairs
     *  @param pp pointer to an object of class Pairs
     *  @param ww the weight
     *  @return none
     */
    void sum (const shared_ptr<Pairs>, const double ww=1);

    /**
     *  @brief estimate the distance between two objects and update
     *	the pair vectors accordingly
     *  @param obj1 pointer to an object of class Object
     *  @param obj2 pointer to an object of class Object
     *  @return none
     */
    void put (const shared_ptr<Object>, const shared_ptr<Object>);
  
  };


  // ============================================================================
  
  /**
   *  @class Pairs3D Pairs.h "Headers/Lib/Pairs.h"
   *
   *  @brief The class Pairs3D
   *
   *  This class is used to handle objects of type <EM> Pairs3D
   *  </EM>, used to measure the 3D two-point correlation function
   */
  class Pairs3D : public Pairs
  {
 
  protected:
    
    /// the number of logarithmic bins
    int m_nlog;

    /// the number of linear bins
    int m_nlin;

    /// the number of angular bins
    int m_ncos;

    /// the minimum separation used to count pairs
    double m_rMin;

    /// the maximum separation used to count pairs
    double m_rMax;

    /// the inverse of the logarithmic bin size
    double m_logbinsz_inv;

    /// the inverse of the linear bin size
    double m_linbinsz_inv;

    /// the inverse of the angular bin size
    double m_cosbinsz_inv;

    /// the number of pairs in logarithmic bins of &theta;
    vector<double> m_PPlog;

    /// the number of pairs in linear bins of &theta;
    vector<double> m_PPlin;

    /// the number of pairs in linear bins both perpendicular and parallel to the line-of-sight
    vector<vector<double> > m_PP2d;

    /// the number of pairs in linear bins perpendicular to the line-of-sight, and in logarithmic bins parallel to the line-of-sight
    vector<vector<double> > m_PPslog;

    /// the number of pairs in logarithmic bins in r and linear angular bins in cos(&theta;)
    vector<vector<double> > m_PPcoslog;

    /// the number of pairs in linear bins in r and linear angular bins in cos(&theta;)
    vector<vector<double> > m_PPcoslin;

  
  public:

    /**
     *  @brief default constructor
     *  @return object of class Pairs3D, with protected members set to -1
     */
    Pairs3D () 
      : m_nlog(-1), m_nlin(-1), m_ncos(-1), m_rMin(-1.), m_rMax(-1.), m_logbinsz_inv(-1.), m_linbinsz_inv(-1.), m_cosbinsz_inv(-1.) {}

    /**
     *  @brief constructor
     *  @param nlog the number of logarithmic bins
     *  @param nlin the number of linear bins
     *  @param ncos the number of angular bins
     *  @param rMin the minimum separation used to count pairs
     *  @param rMax the maximum separation used to count pairs
     *  @param logbinsz the logarithmic bin size
     *  @param linbinsz the linear bin size
     *  @param cosbinsz the angular bin size
     *  @return object of class Pairs3D
     */
    Pairs3D (const int nlog, const int nlin, const int ncos, const double rMin, const double rMax, const double logbinsz, const double linbinsz, const double cosbinsz) 
      : m_nlog(nlog), m_nlin(nlin), m_ncos(ncos), m_rMin(rMin), m_rMax(rMax), m_logbinsz_inv(1./logbinsz), m_linbinsz_inv(1./linbinsz), m_cosbinsz_inv(1./cosbinsz)
      {
	m_PPlog.resize(m_nlog+1, 0.);
	m_PPlin.resize(m_nlin+1, 0.);
	m_PP2d.resize(m_nlin+1, vector<double>(m_nlin+1, 0.));
	m_PPslog.resize(m_nlog+1, vector<double>(m_nlin+1, 0.));
	m_PPcoslog.resize(m_nlog+1, vector<double>(m_ncos+1, 0.));
	m_PPcoslin.resize(m_nlin+1, vector<double>(m_ncos+1, 0.));
      }

    /**
     *  @brief get the private member Pairs3D::m_nlog
     *  @return the number of logarithmic bins
     */
    int nlog () const /*override*/ { return m_nlog; }
    
    /**
     *  @brief get the private member Pairs3D::m_nlin
     *  @return the number of linear bins
     */
    int nlin () const /*override*/ { return m_nlin; }

    /**
     *  @brief get the private member Pairs3D::m_ncos
     *  @return the number of angular bins
     */
    int ncos () const /*override*/ { return m_ncos; }

    /**
     *  @brief get the private member Pairs3D::m_rMin
     *  @return the minimum separation used to count pairs
     */
    double rMin () const /*override*/ { return m_rMin; }

    /**
     *  @brief get the private member Pairs3D::m_rMax
     *  @return the maximum separation used to count pairs
     */
    double rMax () const /*override*/ { return m_rMax; }

    /**
     *  @brief get the private member Pairs3D::m_PPlog[i]
     *  @param i the bin index
     *  @return the number of pairs in the i-th logarithmic bin
     */
    double PPlog (const int i) const /*override*/ { return m_PPlog[i]; }

    /**
     *  @brief get the private member Pairs3D::m_PPlin[i]
     *  @param i the bin index
     *  @return the number of pairs in the i-th linear bin
     */
    double PPlin (const int i) const /*override*/ { return m_PPlin[i]; }

    /**
     *  @brief get the private member Pairs3D::m_PP2d[i][j]
     *  @param i the bin index in the direction perpendicular to the
     *  line-of-sight
     *  @param j the bin index in the direction parallel to the
     *  line-of-sight
     *  @return the number of pairs in the i-th linear bin
     *  perpendicular to the line-of-sight, and in the j-th linear bin
     *  parallel to the line-of-sight
     */
    double PP2d (const int i, const int j) const /*override*/ { return m_PP2d[i][j]; }

    /**
     *  @brief get the private member Pairs3D::m_PPslog[i][j]
     *  @param i the logarithmic bin index in the direction
     *  perpendicular to the line-of-sight
     *  @param j the linear bin index in the direction parallel to the
     *  line-of-sight
     *  @return the number of pairs in the i-th linear bin
     *  perpendicular to the line-of-sight, and in the j-th
     *  logarithmic bin parallel to the line-of-sight
     */
    double PPslog (const int i, const int j) const /*override*/ { return m_PPslog[i][j]; }     

    /**
     *  @brief get the private member Pairs3D::m_PPcoslog[i][j]
     *  @param i the logarithmic bin index 
     *  @param j the angular bin index
     *  @return the number of pairs in the i-th logarithmic bin in r,
     *  and in the j-th linear angular bin in cos(&theta;)
     */
    double PPcoslog (const int i, const int j) const /*override*/ { return m_PPcoslog[i][j]; }

    /**
     *  @brief get the private member Pairs3D::m_PPcoslin[i][j]
     *  @param i the linear bin index 
     *  @param j the angular bin index
     *  @return the number of pairs in the i-th linear bin in r, and
     *  in the j-th linear angular bin in cos(&theta;)
     */
    double PPcoslin (const int i, const int j) const /*override*/ { return m_PPcoslin[i][j]; }   
  
    /**
     *  @brief get the private member Pairs3D::m_PPlog
     *  @return the vector containing the number of pairs in
     *  logarithmic bins
     */
    vector<double> PPlog () const /*override*/ { return m_PPlog; }

    /**
     *  @brief get the private member Pairs3D::m_PPlin
     *  @return the vector containing the number of pairs in linear
     *  bins
     */
    vector<double> PPlin () const /*override*/ { return m_PPlin; }

    /**
     *  @brief get the private member Pairs3D::m_PP2d
     *  @return a matrix containing the number of pairs in linear bins
     *  perpendicular to the line-of-sight, and linear bins parallel
     *  to the line-of-sight
     */
    vector<vector<double> > PP2d () const /*override*/ { return m_PP2d; }

    /**
     *  @brief get the private member Pairs3D::m_PPslog
     *  @return a matrix containing the number of pairs in linear bins
     *  perpendicular to the line-of-sight, and in logarithmic bins
     *  parallel to the line-of-sight
     */
    vector<vector<double> > PPslog () const /*override*/ { return m_PPslog; }

    /**
     *  @brief get the private member Pairs3D::m_PPcoslog
     *  @return a matrix containing the number of pairs in logarithmic
     *  bins in r and linear angular bins in cos(&theta;)
     */
    vector<vector<double> > PPcoslog () const /*override*/ { return m_PPcoslog; }

    /**
     *  @brief get the private member Pairs3D::m_PPcoslin
     *  @return a matrix containing the number of pairs in linear bins
     *  in r and linear angular bins in cos(&theta;)
     */
    vector<vector<double> > PPcoslin () const /*override*/ { return m_PPcoslin; }

    /**
     *  @brief set the private member Pairs3D::m_PPlog[i]
     *  @param i the bin index
     *  @param pp the number of pairs in the bin
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    void set_PPlog (const int i, const double pp) /*override*/ { m_PPlog[i] = pp; }

    /**
     *  @brief set the private member Pairs3D::m_PPlin[i]
     *  @param i the bin index
     *  @param pp the number of pairs in the bin
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    void set_PPlin (const int i, const double pp) /*override*/ { m_PPlin[i] = pp; }

    /**
     *  @brief set the private member Pairs3D::m_PP2d[i][j]
     *  @param i the bin index in the direction perpendicular to the
     *  line-of-sight
     *  @param j the bin index in the direction parallel to the
     *  line-of-sight
     *  @param pp the number of pairs in the bin 
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    void set_PP2d (const int i, const int j, const double pp) /*override*/ { m_PP2d[i][j] = pp; }

    /**
     *  @brief get the private member Pairs3D::m_PPslog[i][j]
     *  @param i the logarithmic bin index in the direction
     *  perpendicular to the line-of-sight
     *  @param j the linear bin index in the direction parallel to the
     *  line-of-sight
     *  @param pp the number of pairs in the bin 
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    void set_PPslog (const int i, const int j, const double pp) /*override*/ {m_PPslog[i][j] = pp; }     

    /**
     *  @brief get the private member Pairs3D::m_PPcoslog[i][j]
     *  @param i the logarithmic bin index 
     *  @param j the angular bin index
     *  @param pp the number of pairs in the bin 
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    void set_PPcoslog (const int i, const int j, const double pp) /*override*/ { m_PPcoslog[i][j] = pp; }

    /**
     *  @brief get the private member Pairs3D::m_PPcoslin[i][j]
     *  @param i the linear bin index 
     *  @param j the angular bin index
     *  @param pp the number of pairs in the bin 
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    void set_PPcoslin (const int i, const int j, const double pp) /*override*/ { m_PPcoslin[i][j] = pp; }
    
    /**
     *  @brief sum the number of pairs
     *  @param pp pointer to an object of class Pairs
     *  @param ww the weight
     *  @return none
     */
    void sum (const shared_ptr<Pairs>, const double ww=1);

    /**
     *  @brief estimate the distance between two objects and update
     *	the pair vectors accordingly
     *  @param obj1 pointer to an object of class Object
     *  @param obj2 pointer to an object of class Object
     *  @return none
     */
    void put (const shared_ptr<Object>, const shared_ptr<Object>);
    
    /**
     *  @brief estimate the distance between two objects and update
     *	the pair vectors in logarithmic bins accordingly
     *  @param obj1 pointer to an object of class Object
     *  @param obj2 pointer to an object of class Object
     *  @return none
     */
    void put_log (const shared_ptr<Object>, const shared_ptr<Object>);
  };
}

#endif
