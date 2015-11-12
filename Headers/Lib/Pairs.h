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
    virtual ~Pairs() {}

    /**
     *  @brief get the private member \e m_nlog
     *  @return the number of logarithmic bins, or an error message if
     *  the derived object does not have this member
     */
    virtual int nlog () { cosmobl::ErrorMsg("Error in Pairs::nlog() of Pairs.h!"); return 0; }

    /**
     *  @brief get the private member \e m_nlin
     *  @return the number of linear bins, or an error message if the
     *  derived object does not have this member
     */
    virtual int nlin () { cosmobl::ErrorMsg("Error in Pairs::nlin() of Pairs.h!"); return 0; }

    /**
     *  @brief get the private member \e m_ncos
     *  @return the number of angular bins, or an error message if the
     *  derived object does not have this member
     */
    virtual int ncos () { cosmobl::ErrorMsg("Error in Pairs::ncos() of Pairs.h!"); return 0; }

    /**
     *  @brief get the private member \e m_rMin
     *  @return the minimum separation used to count pairs, or an
     *  error message if the derived object does not have this member
     */
    virtual double rMin () { cosmobl::ErrorMsg("Error in Pairs::rMin() of Pairs.h!"); return 0; }

    /**
     *  @brief get the private member \e m_rMax
     *  @return the maximum separation used to count pairs, or an
     *  error message if the derived object does not have this member
     */
    virtual double rMax () { cosmobl::ErrorMsg("Error in Pairs::rMax() of Pairs.h!"); return 0; }

    /**
     *  @brief get the private member \e m_PPlog[i]
     *  @param i the bin index
     *  @return the number of pairs in the i-th logarithmic bin, or an
     *  error message if the derived object does not have this member
     */
    virtual double PPlog (int i) { cosmobl::ErrorMsg("Error in Pairs::PPlog() of Pairs.h!"); return 0; }

    /**
     *  @brief get the private member \e m_PPlin[i]
     *  @param i the bin index
     *  @return the number of pairs in the i-th linear bin, or an
     *  error message if the derived object does not have this member
     */
    virtual double PPlin (int i) { cosmobl::ErrorMsg("Error in Pairs::PPlin() of Pairs.h!"); return 0; }

    /**
     *  @brief get the private member \e m_PPlog
     *  @return the vector containing the number of pairs in
     *  logarithmic bins, or an error message if the derived object
     *  does not have this member
     */
    virtual vector<double> PPlog() { cosmobl::ErrorMsg("Error in Pairs::PPlog() of Pairs.h!"); vector<double> PP; return PP; }

    /**
     *  @brief get the private member \e m_PPlin
     *  @return the vector containing the number of pairs in linear
     *  bins, or an error message if the derived object does not have
     *  this member
     */
    virtual vector<double> PPlin() { cosmobl::ErrorMsg("Error in Pairs::PPlin() of Pairs.h!"); vector<double> PP; return PP; }
    
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
    virtual double PP2d (int i, int j) { cosmobl::ErrorMsg("Error in Pairs::PP2d() of Pairs.h!"); return 0; }

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
    virtual double PPslog (int i, int j) { cosmobl::ErrorMsg("Error in Pairs::PPslog() of Pairs.h!"); return 0; }     

    /**
     *  @brief get the private member \e m_PPcoslog[i][j]
     *  @param i the logarithmic bin index 
     *  @param j the angular bin index
     *  @return the number of pairs in the i-th logarithmic bin in r,
     *  and in the j-th linear angular bin in cos(&theta;), or an
     *  error message if the derived object does not have this member
     */
    virtual double PPcoslog (int i, int j) { cosmobl::ErrorMsg("Error in Pairs::PPcoslog() of Pairs.h!"); return 0; }

    /**
     *  @brief get the private member \e m_PPcoslin[i][j]
     *  @param i the linear bin index 
     *  @param j the angular bin index
     *  @return the number of pairs in the i-th linear bin in r, and
     *  in the j-th linear angular bin in cos(&theta;), or an error
     *  message if the derived object does not have this member
     */
    virtual double PPcoslin (int i, int j) { cosmobl::ErrorMsg("Error in Pairs::PPcoslin() of Pairs.h!"); return 0; }
   
    /**
     *  @brief get the private member \e m_PP2d
     *  @return a matrix containing the number of pairs in linear bins
     *  perpendicular to the line-of-sight, and linear bins parallel
     *  to the line-of-sight, or an error message if the derived
     *  object does not have this member
     */
    virtual vector<vector<double> > PP2d() { cosmobl::ErrorMsg("Error in Pairs::PP2d() of Pairs.h!"); vector<vector<double> > PP; return PP; }
   
    /**
     *  @brief get the private member \e m_PPslog
     *  @return a matrix containing the number of pairs in linear bins
     *  perpendicular to the line-of-sight, and in logarithmic bins
     *  parallel to the line-of-sight, or an error message if the
     *  derived object does not have this member
     */
    virtual vector<vector<double> > PPslog() { cosmobl::ErrorMsg("Error in Pairs::PPslog() of Pairs.h!"); vector<vector<double> > PP; return PP; }

    /**
     *  @brief get the private member \e m_PPcoslog
     *  @return a matrix containing the number of pairs in logarithmic
     *  bins in r and linear angular bins in cos(&theta;), or an error
     *  message if the derived object does not have this member
     */
    virtual vector<vector<double> > PPcoslog() { cosmobl::ErrorMsg("Error in Pairs::PPcoslog() of Pairs.h!"); vector<vector<double> > PP; return PP; }

    /**
     *  @brief get the private member \e m_PPcoslin
     *  @return a matrix containing the number of pairs in linear bins
     *  in r and linear angular bins in cos(&theta;), or an error
     *  message if the derived object does not have this member
     */
    virtual vector<vector<double> > PPcoslin() { cosmobl::ErrorMsg("Error in Pairs::PPcoslin() of Pairs.h!"); vector<vector<double> > PP; return PP; }

    /**
     *  @brief set the private member \e m_PPlog[i]
     *  @param i the bin index
     *  @param pp the number of pairs in the bin
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_PPlog (int i, double pp) { cosmobl::ErrorMsg("Error in Pairs::PPlog() of Pairs.h!"); }

    /**
     *  @brief set the private member \e m_PPlin[i]
     *  @param i the bin index
     *  @param pp the number of pairs in the bin
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_PPlin (int i, double pp) { cosmobl::ErrorMsg("Error in Pairs::PPlin() of Pairs.h!"); }

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
    virtual void set_PP2d (int i, int j, double pp) { cosmobl::ErrorMsg("Error in Pairs::PP2d() of Pairs.h!"); }

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
    virtual void set_PPslog (int i, int j, double pp) { cosmobl::ErrorMsg("Error in Pairs::PPslog() of Pairs.h!"); }     

    /**
     *  @brief get the private member \e m_PPcoslog[i][j]
     *  @param i the logarithmic bin index 
     *  @param j the angular bin index
     *  @param pp the number of pairs in the bin 
     *  @return none, or an error message if the derived object does not
     *  have this member
     */
    virtual void set_PPcoslog (int i, int j, double pp) { cosmobl::ErrorMsg("Error in Pairs::PPcoslog() of Pairs.h!"); }

    /**
     *  @brief get the private member \e m_PPcoslin[i][j]
     *  @param i the linear bin index 
     *  @param j the angular bin index
     *  @param pp the number of pairs in the bin 
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void set_PPcoslin (int i, int j, double pp) { cosmobl::ErrorMsg("Error in Pairs::PPcoslin() of Pairs.h!"); }
    
    /**
     *  @brief sum the number of pairs
     *  @param pp pointer to an object of class Pairs
     *  @param ww the weight
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void sum (shared_ptr<Pairs> pp, double ww=1) { cosmobl::ErrorMsg("Error in Pairs::sum() of Pairs.h!"); }

    /**
     *  @brief estimate the distance between two objects and update
     *	the pair vectors accordingly
     *  @param obj1 pointer to an object of class Object
     *  @param obj2 pointer to an object of class Object
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void put (shared_ptr<Object> obj1, shared_ptr<Object> obj2) { cosmobl::ErrorMsg("Error in Pairs::put() of Pairs.h!"); }

    /**
     *  @brief estimate the distance between two objects and update
     *	the logarithmic binned pair vector accordingly
     *  @param obj1 pointer to an object of class Object
     *  @param obj2 pointer to an object of class Object
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    virtual void put_log (shared_ptr<Object> obj1, shared_ptr<Object> obj2) { cosmobl::ErrorMsg("Error in Pairs::put_log() of Pairs.h!"); }

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
    Pairs2D (int &nlog, int &nlin, double &thetaMin, double &thetaMax, double &logbinsz, double &linbinsz) 
      : m_nlog(nlog), m_nlin(nlin), m_thetaMin(thetaMin), m_thetaMax(thetaMax), m_logbinsz_inv(1./logbinsz), m_linbinsz_inv(1./linbinsz)
      { 
	m_PPlog.resize(m_nlog+1, 0.);
	m_PPlin.resize(m_nlin+1, 0.);
      }
  
    /**
     *  @brief get the protected member Pairs2D::m_nlog
     *  @return the number of logarithmic bins
     */
    int nlog () { return m_nlog; }

    /**
     *  @brief get the protected member Pairs2D::m_nlin
     *  @return the number of linear bins
     */
    int nlin () { return m_nlin; }

    /**
     *  @brief get the protected member Pairs2D::m_thetaMin
     *  @return the minimum value of the angle &theta; used to count
     *  the number of pairs
     */
    double thetaMin() { return m_thetaMin; }

    /**
     *  @brief get the protected member Pairs2D::m_thetaMax
     *  @return the maximum value of the angle &theta; used to count
     *  the number of pairs
     */    
    double thetaMax() { return m_thetaMax; }

    /**
     *  @brief get the private member Pairs2D::m_PPlog[i]
     *  @param i the bin index
     *  @return the number of pairs in the i-th logarithmic bin
     */
    double PPlog (int i) { return m_PPlog[i]; }

    /**
     *  @brief get the private member Pairs2D::m_PPlin[i]
     *  @param i the bin index
     *  @return the number of pairs in the i-th linear bin
     */
    double PPlin (int i) { return m_PPlin[i]; }

    /**
     *  @brief get the private member Pairs2D::m_PPlog
     *  @return the vector containing the number of pairs in
     *  logarithmic bins
     */
    vector<double> PPlog() { return m_PPlog; }

    /**
     *  @brief get the private member Pairs2D::m_PPlin
     *  @return the vector containing the number of pairs in linear
     *  bins
     */
    vector<double> PPlin() { return m_PPlin; }

    /**
     *  @brief set the private member Pairs2D::m_PPlog[i]
     *  @param i the bin index
     *  @param pp the number of pairs in the bin
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    void set_PPlog (int i, double pp) { m_PPlog[i] = pp; }

    /**
     *  @brief set the private member Pairs2D::m_PPlin[i]
     *  @param i the bin index
     *  @param pp the number of pairs in the bin
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    void set_PPlin (int i, double pp) { m_PPlin[i] = pp; }
    
    /**
     *  @brief sum the number of pairs
     *  @param pp pointer to an object of class Pairs
     *  @param ww the weight
     *  @return none
     */
    void sum (shared_ptr<Pairs>, double ww=1);

    /**
     *  @brief estimate the distance between two objects and update
     *	the pair vectors accordingly
     *  @param obj1 pointer to an object of class Object
     *  @param obj2 pointer to an object of class Object
     *  @return none
     */
    void put (shared_ptr<Object>, shared_ptr<Object>);
  
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
    Pairs3D (int &nlog, int &nlin, int &ncos, double &rMin, double &rMax, double &logbinsz, double &linbinsz, double &cosbinsz) 
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
    int nlog() { return m_nlog; }
    
    /**
     *  @brief get the private member Pairs3D::m_nlin
     *  @return the number of linear bins
     */
    int nlin() { return m_nlin; }

    /**
     *  @brief get the private member Pairs3D::m_ncos
     *  @return the number of angular bins
     */
    int ncos() { return m_ncos; }

    /**
     *  @brief get the private member Pairs3D::m_rMin
     *  @return the minimum separation used to count pairs
     */
    double rMin() { return m_rMin; }

    /**
     *  @brief get the private member Pairs3D::m_rMax
     *  @return the maximum separation used to count pairs
     */
    double rMax() { return m_rMax; }

    /**
     *  @brief get the private member Pairs3D::m_PPlog[i]
     *  @param i the bin index
     *  @return the number of pairs in the i-th logarithmic bin
     */
    double PPlog (int i) { return m_PPlog[i]; }

    /**
     *  @brief get the private member Pairs3D::m_PPlin[i]
     *  @param i the bin index
     *  @return the number of pairs in the i-th linear bin
     */
    double PPlin (int i) { return m_PPlin[i]; }

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
    double PP2d (int i, int j) { return m_PP2d[i][j]; }

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
    double PPslog (int i, int j) { return m_PPslog[i][j]; }     

    /**
     *  @brief get the private member Pairs3D::m_PPcoslog[i][j]
     *  @param i the logarithmic bin index 
     *  @param j the angular bin index
     *  @return the number of pairs in the i-th logarithmic bin in r,
     *  and in the j-th linear angular bin in cos(&theta;)
     */
    double PPcoslog (int i, int j) { return m_PPcoslog[i][j]; }

    /**
     *  @brief get the private member Pairs3D::m_PPcoslin[i][j]
     *  @param i the linear bin index 
     *  @param j the angular bin index
     *  @return the number of pairs in the i-th linear bin in r, and
     *  in the j-th linear angular bin in cos(&theta;)
     */
    double PPcoslin (int i, int j) { return m_PPcoslin[i][j]; }   
  
    /**
     *  @brief get the private member Pairs3D::m_PPlog
     *  @return the vector containing the number of pairs in
     *  logarithmic bins
     */
    vector<double> PPlog() { return m_PPlog; }

    /**
     *  @brief get the private member Pairs3D::m_PPlin
     *  @return the vector containing the number of pairs in linear
     *  bins
     */
    vector<double> PPlin() { return m_PPlin; }

    /**
     *  @brief get the private member Pairs3D::m_PP2d
     *  @return a matrix containing the number of pairs in linear bins
     *  perpendicular to the line-of-sight, and linear bins parallel
     *  to the line-of-sight
     */
    vector<vector<double> > PP2d() { return m_PP2d; }

    /**
     *  @brief get the private member Pairs3D::m_PPslog
     *  @return a matrix containing the number of pairs in linear bins
     *  perpendicular to the line-of-sight, and in logarithmic bins
     *  parallel to the line-of-sight
     */
    vector<vector<double> > PPslog() { return m_PPslog; }

    /**
     *  @brief get the private member Pairs3D::m_PPcoslog
     *  @return a matrix containing the number of pairs in logarithmic
     *  bins in r and linear angular bins in cos(&theta;)
     */
    vector<vector<double> > PPcoslog() { return m_PPcoslog; }

    /**
     *  @brief get the private member Pairs3D::m_PPcoslin
     *  @return a matrix containing the number of pairs in linear bins
     *  in r and linear angular bins in cos(&theta;)
     */
    vector<vector<double> > PPcoslin() { return m_PPcoslin; }

    /**
     *  @brief set the private member Pairs3D::m_PPlog[i]
     *  @param i the bin index
     *  @param pp the number of pairs in the bin
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    void set_PPlog (int i, double pp) { m_PPlog[i] = pp; }

    /**
     *  @brief set the private member Pairs3D::m_PPlin[i]
     *  @param i the bin index
     *  @param pp the number of pairs in the bin
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    void set_PPlin (int i, double pp) { m_PPlin[i] = pp; }

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
    void set_PP2d (int i, int j, double pp) { m_PP2d[i][j] = pp; }

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
    void set_PPslog (int i, int j, double pp) {m_PPslog[i][j] = pp; }     

    /**
     *  @brief get the private member Pairs3D::m_PPcoslog[i][j]
     *  @param i the logarithmic bin index 
     *  @param j the angular bin index
     *  @param pp the number of pairs in the bin 
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    void set_PPcoslog (int i, int j, double pp) { m_PPcoslog[i][j] = pp; }

    /**
     *  @brief get the private member Pairs3D::m_PPcoslin[i][j]
     *  @param i the linear bin index 
     *  @param j the angular bin index
     *  @param pp the number of pairs in the bin 
     *  @return none, or an error message if the derived object does
     *  not have this member
     */
    void set_PPcoslin (int i, int j, double pp) { m_PPcoslin[i][j] = pp; }
    
    /**
     *  @brief sum the number of pairs
     *  @param pp pointer to an object of class Pairs
     *  @param ww the weight
     *  @return none
     */
    void sum (shared_ptr<Pairs>, double ww=1);

    /**
     *  @brief estimate the distance between two objects and update
     *	the pair vectors accordingly
     *  @param obj1 pointer to an object of class Object
     *  @param obj2 pointer to an object of class Object
     *  @return none
     */
    void put (shared_ptr<Object>, shared_ptr<Object>);
    
    void put_log (shared_ptr<Object>, shared_ptr<Object>);
  };
}

#endif
