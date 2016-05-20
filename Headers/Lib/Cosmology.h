/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 ********************************************************************/

/**
 *  @file Headers/Lib/Cosmology.h
 *
 *  @brief The class Cosmology 
 *
 *  This file defines the interface of the class Cosmology, used for
 *  all kind of cosmological calculations
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __COSM__
#define __COSM__


#include "Likelihood.h"
#include "power_whu.h"


// ===================================================================================================


namespace cosmobl {

  /**
   * @enum CosmoPar
   * @brief the cosmological parameters
   */
  enum CosmoPar {

    /// &Omega;<SUB>M</SUB>: the density of baryons, cold dark matter and massive neutrinos (in units of the critical density) at z=0, LCDM case
    _Omega_matter_LCDM_,

    /// &Omega_,<SUB>M</SUB>: the density of baryons, cold dark matter and massive neutrinos (in units of the critical density) at z=0
    _Omega_matter_,
  
    /// &Omega_,<SUB>b</SUB>: the baryon density at z=0
    _Omega_baryon_,         

    /// &Omega_,<SUB>&nu_,</SUB>: the density of massive neutrinos at z=0
    _Omega_neutrinos_,      

    /// N<SUB>eff</SUB>: the effective number (for QED + non-instantaneous decoupling)
    _massless_neutrinos_,   

    /// the number of degenerate massive neutrino species 
    _massive_neutrinos_,    

    /// &Omega_,<SUB>DE</SUB>: the dark energy density at z=0
    _Omega_DE_,             

    /// &Omega_,<SUB>rad</SUB>: the radiation density at z=0
    _Omega_radiation_,     
             
    /// H<SUB>0</SUB>: the Hubble constant at z=0 [km/sec/Mpc] 
    _H0_,

    /// A<SUB>s</SUB>: the initial scalar amplitude of the power spectrum
    _scalar_amp_,           

    /// n<SUB>spec</SUB>: the primordial spectral index
    _n_spec_,               

    /// w<SUB>0</SUB>: the parameter of the dark energy equation of state (CPL parameterisation)
    _w0_,

    /// w<SUB>a</SUB>: the parameter of the dark energy equation of state (CPL parameterisation)
    _wa_,             

    /// f<SUB>NL</SUB>: the non-Gaussian amplitude
    _fNL_,                  

    /// &sigma_,<SUB>8</SUB>: the power spectrum normalization
    _sigma8_ 

  };


  /**
   *  @class Cosmology Cosmology.h "Headers/Lib/Cosmology.h"
   *
   *  @brief The class Cosmology
   *
   *  This class is used to handle objects of type <EM> cosmology
   *  </EM>. It can be used to estimate i) any kind of cosmic
   *  distances and volumes, ii) the halo mass function and bias, iii)
   *  the real-space and redshift-space dark matter power spectrum and
   *  correlation function (with both CAMB, CLASS, MPTbreeze and the
   *  analytic fitting forms by Eisenstein & Hu), iv) the accreted
   *  mass function, v) several BAO parameters, and vi) the halo mass
   *  function and bias in different non-Gaussian cosmological
   *  frameworks.
   */
  class Cosmology {

  private:
    
    /// &Omega;<SUB>M</SUB>: the density of baryons, cold dark matter and massive neutrinos (in units of the critical density) at z=0
    double m_Omega_matter;         
   
    /// &Omega;<SUB>b</SUB>: the baryon density at z=0
    double m_Omega_baryon;         

    /// &Omega;<SUB>&nu;</SUB>: the density of massive neutrinos at z=0
    double m_Omega_neutrinos;      

    /// N<SUB>eff</SUB>: the effective number (for QED + non-instantaneous decoupling)
    double m_massless_neutrinos;   

    /// the number of degenerate massive neutrino species 
    int m_massive_neutrinos;    

    /// &Omega;<SUB>DE</SUB>: the dark energy density at z=0
    double m_Omega_DE;             

    /// &Omega;<SUB>rad</SUB>: the radiation density at z=0
    double m_Omega_radiation;     

    /// &Omega;<SUB>k</SUB>: the density of curvature energy
    double m_Omega_k;              

    /// &Omega;<SUB>CDM</SUB>: the cold dark matter density at z=0
    double m_Omega_CDM;   
             
    /// H<SUB>0</SUB>: the Hubble constant at z=0 [km/sec/Mpc] 
    double m_H0;

    /// \e h: the Hubble parameter, H<SUB>0</SUB>/100
    double m_hh;

    /// t<SUB>H</SUB>: the Hubble time
    double m_t_H;                  

    /// D<SUB>H</SUB>: the Hubble distance
    double m_D_H;                  

    /// &sigma;<SUB>8</SUB>: the power spectrum normalization
    double m_sigma8;               

    /// A<SUB>s</SUB>: the initial scalar amplitude of the power spectrum
    double m_scalar_amp;           

    /// n<SUB>spec</SUB>: the primordial spectral index
    double m_n_spec;               

    /// w<SUB>0</SUB>: the parameter of the dark energy equation of state (CPL parameterisation)
    double m_w0;

    /// w<SUB>a</SUB>: the parameter of the dark energy equation of state (CPL parameterisation)
    double m_wa;             
    
    /// &rho;<SUB>0</SUB>: the mean density of the Universe at z=0 [Msun*Mpc^-3]
    double m_RhoZero;             

    /// f<SUB>NL</SUB>: the non-Gaussian amplitude
    double m_fNL;                  

    /// the non-Gaussian shape (type=1 local, type=2 equilateral, type=3 enfolded, type=4 orthogonal)
    int m_type_NG;    

    /// the normalization of the power spectrum for Eisenstein & Hu [http://background.uchicago.edu/~whu/transfer/transferpage.html]
    double m_Pk0_EH;

    /// the normalization of the power spectrum for CAMB [http://camb.info/]
    double m_Pk0_CAMB;

    /// the normalization of the power spectrum for MPTbreeze [http://arxiv.org/abs/1207.1465]
    double m_Pk0_MPTbreeze; 

    /// the normalization of the power spectrum for CLASS [http://class-code.net/]
    double m_Pk0_CLASS; 
        
    /// the cosmologial model used to compute distances
    string m_model;                

    /// 0 &rarr; phyical units; 1 &rarr; cosmological units (i.e. without \e h)
    bool m_unit;                   


    /**
     *  @brief auxiliary function to create a grid file with
     *  &sigma;(M) 
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *     
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return file_grid name of the file where the grid is stored
     */
    string create_grid_sigmaM (const string, const double, const string output_root="test", const string interpType="Linear", const int Num=-1, const double stepsize=100., const double k_max=100., const string file_par=par::defaultString) const;         

    /**
     *  @brief auxiliary function to compute the mass function
     *
     *  @param Mass mass
     *
     *  @param Sigma &sigma;(mass): the mass variance
     *
     *  @param Dln_Sigma dln&sigma;/dM: the derivative of the mass
     *  variance
     *
     *  @param redshift redshift
     *
     *  @param author_MF author(s) who proposed the mass function; valid
     *  authors are: PS (Press & Schechter), ST (Sheth & Tormen),
     *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
     *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
     *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
     *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
     *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
     *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
     *  al. 2012)
     *
     *  @param Delta &Delta;: the overdensity, defined as the mean
     *  interior density relative to the background
     *
     *  @return the mass function, d&Phi;/dM=dn(M)/dM
     *
     *  @warning the input parameter &Delta; is used only in the
     *  Tinker mass function; it is the background overdensity, not
     *  the critical overdensity &Delta;<SUB>crit</SUB>; the function
     *  cosmobl::Cosmology::DeltaR can be used to convert
     *  &Delta;<SUB>crit</SUB> into &Delta;
     */
    double MF_generator (const double Mass, const double Sigma, const double Dln_Sigma, const double redshift, const string author_MF, const double Delta=200.); 

    /**
     *  @brief auxiliary function to compute the halo bias
     *
     *  @param Sigma &sigma;(mass): the mass variance
     *
     *  @param redshift redshift
     *
     *  @param author author(s) who proposed the bias; valid authors
     *  are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen
     *  2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction
     *  of Warren 2004), Tinker (Tinker et al. 2010)
     *
     *  @param Delta &Delta;: the overdensity, defined as the mean
     *  interior density relative to the background
     *
     *  @return the halo bias
     *
     *  @warning the input parameter &Delta; is used only in the
     *  Tinker mass function; it is the background overdensity, not
     *  the critical overdensity &Delta;<SUB>crit</SUB>; the function
     *  cosmobl::Cosmology::DeltaR can be used to convert
     *  &Delta;<SUB>crit</SUB> into &Delta;
     */
    double bias_halo_generator (const double, const double, const string, const double Delta=200.) const;    

    /**
     *  @brief the incomplete elliptic integral 
     *
     *  this method computes the incomplete elliptic integral
     *  of the first kind F(&phi;|m), with shape parameter fixed to
     *  m=(2+3<SUP>0.5</SUP>)/4; it is used to compute the relation
     *  between comoving distance and redshift (see
     *  Numer. Math. (2010) 116:687–719, T. Fukushima "Fast
     *  computation of incomplete elliptic integral of first kind by
     *  half argument transformation")
     *
     *  @author Mauro Roncarelli
     *  @author mauro.roncarelli@unibo.it
     *
     *  @param phi &phi; scalar, argument
     *
     *  @return F(&phi;|m)
     *
     *  @warning this method works only for a fixed value of the shape
     *  parameter m=(2+3<SUP>0.5</SUP>)/4
     */
    double elf_dz (const double) const;

    /**
     *  @brief the inverse cosine amplitude of the Jacobian elliptic
     *  function
     *
     *  this method computes the inverse cosine amplitude of the
     *  Jacobian elliptic function cn<SUP>-1</SUP>(c|m), with shape
     *  parameter fixed to m=(2+3<SUP>0.5</SUP>)/4; namely
     *  cn<SUP>-1</SUP>(s|m) = F(arccos(c)|m); it is used to compute
     *  the relation between comoving distance and redshift (see
     *  Numer. Math. (2010) 116:687–719, T. Fukushima "Fast
     *  computation of incomplete elliptic integral of first kind by
     *  half argument transformation")
     *
     *  @author Mauro Roncarelli
     *  @author mauro.roncarelli@unibo.it
     *
     *  @param cc scalar, argument
     *
     *  @return cn<SUP>-1</SUP>(c|m)
     *
     *  @warning this method works only for a fixed value of the shape
     *  parameter m=(2+3<SUP>0.5</SUP>)/4; the argument c must be in
     *  the range 0<c<1
     */
    double acn_dz (const double) const;

    /**
     *  @brief the inverse sine amplitude of the Jacobian elliptic
     *  function
     *
     *  this method computes the inverse sine amplitude of the
     *  Jacobian elliptic function sn<SUP>-1</SUP>(s|m), with shape
     *  parameter fixed to m=(2+3<SUP>0.5</SUP>)/4, namely
     *  sn<SUP>-1</SUP>(s|m) = F(arcsin(s)|m); it is used to compute
     *  the relation between comoving distance and redshift (see
     *  Numer. Math. (2010) 116:687–719, T. Fukushima "Fast
     *  computation of incomplete elliptic integral of first kind by
     *  half argument transformation")
     *
     *  @author Mauro Roncarelli
     *  @author mauro.roncarelli@unibo.it
     *
     *  @param ss scalar, argument
     *
     *  @return sn<SUP>-1</SUP>(s|m)
     *
     *  @warning this method works only for a fixed value of the shape
     *  parameter m=(2+3<SUP>0.5</SUP>)/4; the argument s must be in
     *  the range 0<s<1
     */
    double asn_dz (const double) const;

    /**
     *  @brief the inverse truncated series necessary to compute
     *  sn<SUP>-1</SUP>(s|m) in ASN_DZ
     *
     *  this method computes the inverse truncated series necessary to
     *  compute sn<SUP>-1</SUP>(s|m) in ASN_DZ; the shape parameter is
     *  fixed to m=(2+3<SUP>0.5</SUP>)/4; it is used to compute the
     *  relation between comoving distance and redshift (see
     *  Numer. Math. (2010) 116:687–719, T. Fukushima "Fast
     *  computation of incomplete elliptic integral of first kind by
     *  half argument transformation"; the output of the function
     *  corresponds to the second factor of eq. (23))
     *
     *  @author Mauro Roncarelli
     *  @author mauro.roncarelli@unibo.it
     *
     *  @param yy scalar, argument
     *
     *  @return \f$\sum_{l=0}^L u_l(m)y^l\f$, with L fixed to 9
     *
     *  @warning this method works only for a fixed value of the shape
     *  parameter m=(2+3<SUP>0.5</SUP>)/4
     */
    double serf_dz (const double) const;


    // -----------------------------------------------------------------------


  public:
    
    /**
     *  @name Constructors/destructors
     */
    ///@{
    
    /**
     *  @brief constructor
     *
     *  @param Omega_matter &Omega;<SUB>M</SUB>: the density of
     *  baryons, cold dark matter and massive neutrinos (in units of
     *  the critical density) at z=0
     *
     *  @param Omega_baryon &Omega;<SUB>b</SUB>: the density of
     *  baryons at z=0
     *
     *  @param Omega_neutrinos &Omega;<SUB>&nu;</SUB>: the density of
     *  massive neutrinos at z=0
     *
     *  @param massless_neutrinos the effective number (for QED +
     *  non-instantaneous decoupling)
     *
     *  @param massive_neutrinos the number of degenerate massive
     *  neutrino species
     *
     *  @param Omega_DE &Omega;<SUB>DE</SUB>: the density of dark
     *  energy at z=0
     *
     *  @param Omega_radiation &Omega;<SUB>rad</SUB>: the density of
     *  radiation at z=0
     *
     *  @param hh \e h: the Hubble parameter, H<SUB>0</SUB>/100 
     *
     *  @param scalar_amp A<SUB>s</SUB>: the initial scalar amplitude
     *  of the power spectrum
     *
     *  @param n_spec n<SUB>spec</SUB>: the primordial spectral index
     *
     *  @param w0 w<SUB>0</SUB>: one of the two parameters of the dark energy
     *  equation of state (CPL parameterisation)
     *
     *  @param wa w<SUB>a</SUB>: one of the two parameters of the dark energy
     *  equation of state (CPL parameterisation)
     *
     *  @param fNL f<SUB>NL</SUB>: the non-Gaussian amplitude
     *
     *  @param type_NG the non-Gaussian shape (type=1 local, type=2
     *  equilateral, type=3 enfolded, type=4 orthogonal)
     *
     *  @param model the cosmologial model used to compute distances
     *
     *  @param unit 0 &rarr; phyical units; 1 &rarr; cosmological units
     *  (i.e. without \e h)
     *
     *  @return object of class Cosmology; by default:
     *  &Omega;<SUB>k</SUB>=1-&Omega;<SUB>M</SUB>-&Omega;<SUB>rad</SUB>-&Omega;<SUB>DE</SUB>,
     *  &Omega;<SUB>CDM</SUB>=&Omega;<SUB>M</SUB>-&Omega;<SUB>b</SUB>-&Omega;<SUB>&nu;</SUB>,
     *  t<SUB>H</SUB>=1/H<SUB>0</SUB>, D<SUB>H</SUB>=c/t<SUB>H</SUB>,
     *  &rho;<SUB>0</SUB>=&rho;<SUB>0</SUB>(&Omega;<SUB>M</SUB>,&Omega;<SUB>&nu;</SUB>),
     *  Pk0_*=1
     */
    Cosmology (const double Omega_matter=0.27, const double Omega_baryon=0.046, const double Omega_neutrinos=0., const double massless_neutrinos=3.04, const int massive_neutrinos=0, const double Omega_DE=0.73, const double Omega_radiation=0., const double hh=0.7, const double scalar_amp=2.46e-9, const double n_spec=0.96, const double w0=-1., const double wa=0., const double fNL=0., const int type_NG=1, const string model="LCDM", const bool unit=1);    

    /**
     *  @brief default destructor
     *  @return none
     */
    ~Cosmology () {};
    
    ///@}


    /**
     *  @name Functions to get the private members of the class
     */
    ///@{

    /**
     *  @brief get the private member specified by the enum CosmoPar
     *  
     *  @param parameter the cosmological parameter
     *
     *  @return value of the cosmological parameter
     */
    double value (const CosmoPar parameter) const;

    /**
     *  @brief get the private member Cosmology::m_Omega_matter
     *
     *  @return &Omega;<SUB>M</SUB>: the matter density, i.e. the
     *  density of baryons, cold dark matter, and massive neutrinos
     *  (in units of the critical density)
     */
    double Omega_matter () const { return m_Omega_matter; }; 

    /**
     *  @brief get the private member Cosmology::m_Omega_baryon
     *
     *  @return &Omega;<SUB>b</SUB>: the baryons density
     */
    double Omega_baryon () const { return m_Omega_baryon; };

    /**
     *  @brief get the private member Cosmology::m_Omega_neutrinos
     *
     *  @return &Omega;<SUB>&nu;</SUB>: the density of massive
     *  neutrinos
     */
    double Omega_neutrinos () const { return m_Omega_neutrinos; };

    /**
     *  @brief get the private member Cosmology::m_massless_neutrinos
     *
     *  @return N<SUB>eff</SUB>: the effective number (for QED + non-instantaneous decoupling)
     */
    double massless_neutrinos () const { return m_massless_neutrinos; };

    /**
     *  @brief get the private member Cosmology::m_massive_neutrinos
     *
     *  @return the number of degenerate massive neutrino species
     */
    int massive_neutrinos () const { return m_massive_neutrinos; };

    /**
     *  @brief get the private member Cosmology::m_Omega_DE
     *
     *  @return &Omega;<SUB>DE</SUB>: the dark energy density
     */
    double Omega_DE () const { return m_Omega_DE; };

    /**
     *  @brief get the private member Cosmology::m_Omega_radiation
     *
     *  @return &Omega;<SUB>rad</SUB>: the radiation density
     */
    double Omega_radiation () const { return m_Omega_radiation; };

    /**
     *  @brief get the private member Cosmology::m_Omega_k
     *
     *  @return &Omega;<SUB>k</SUB>: the density of curvature energy
     */
    double Omega_k () const { return m_Omega_k; };

    /**
     *  @brief get the private member Cosmology::m_Omega_CDM
     *
     *  @return &Omega;<SUB>CDM</SUB>: the cold dark matter density
     */
    double Omega_CDM () const { return m_Omega_CDM; };

    /**
     *  @brief get the private member Cosmology::m_H0
     *
     *  @return H<SUB>0</SUB>: the Hubble constant [km/sec/Mpc]
     */
    double H0 () const { return m_H0; };

    /**
     *  @brief get the private member Cosmology::m_hh
     *
     *  @return \e h: the Hubble parameter, H<SUB>0</SUB>/100
     */
    double hh () const { return m_hh; };

    /**
     *  @brief get the private member Cosmology::m_t_H
     *
     *  @return t<SUB>H</SUB>: the Hubble time
     */
    double t_H () const { return m_t_H; };

    /**
     *  @brief get the private member Cosmology::m_D_H
     *
     *  @return D<SUB>H</SUB>: the Hubble distance
     */
    double D_H () const { return m_D_H; };

    /**
     *  @brief get the private member Cosmology::m_sigma8
     *
     *  @return &sigma;<SUB>8</SUB>: the power spectrum normalization
     */
    double sigma8 () const { return m_sigma8; };

    /**
     *  @brief get the private member Cosmology::m_scalar_amp
     *
     *  @return A<SUB>s</SUB>: the initial scalar amplitude of the
     *  power spectrum
     */
    double scalar_amp () const { return m_scalar_amp; };

    /**
     *  @brief get the private member Cosmology::m_n_spec
     *
     *  @return n<SUB>spec</SUB>: the primordial spectral index
     */
    double n_spec () const { return m_n_spec; };

    /**
     *  @brief get the private member Cosmology::m_w0
     *
     *  @return w<SUB>0</SUB>: one of the parameters of the dark
     *  energy equation of state (CPL parameterisation)
     */
    double w0 () const { return m_w0; };

    /**
     *  @brief get the private member Cosmology::m_wa
     *
     *  @return w<SUB>a</SUB>: one of the parameters of the dark
     *  energy equation of state (CPL parameterisation)
     */
    double wa () const { return m_wa; };

    /**
     *  @brief get the private member Cosmology::m_RhoZero
     *
     *  @return &rho;<SUB>0</SUB>: the mean density of the Universe at
     *  z=0 [Msun*Mpc^-3]
     */
    double RhoZero () const { return m_RhoZero; };

    /**
     *  @brief get the private member Cosmology::m_fNL
     *
     *  @return f<SUB>NL</SUB>: the non-Gaussian amplitude
     */
    double fNL () const { return m_fNL; };

    /**
     *  @brief get the private member Cosmology::m_type_NG
     *
     *  @return the non-Gaussian shape (type=1 local, type=2
     *  equilateral, type=3 enfolded, type=4 orthogonal)
     */
    int type_NG () const { return m_type_NG; }; 

    /**
     *  @brief get the private member Cosmology::m_Pk0_EH
     *
     *  @return the normalization of the power spectrum for Eisenstein
     *  & Hu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     */
    double Pk0_EH () const { return m_Pk0_EH; };

    /**
     *  @brief get the private member Cosmology::m_Pk0_CAMB
     *
     *  @return the normalization of the power spectrum for CAMB [http://camb.info/]
     */
    double Pk0_CAMB () const { return m_Pk0_CAMB; };

    /**
     *  @brief get the private member Cosmology::m_Pk0_MPTbreeze 
     *
     *  @return the normalization of the power spectrum for MPTbreeze
     *  [http://arxiv.org/abs/1207.1465]
     */
    double Pk0_MPTbreeze () const { return m_Pk0_MPTbreeze; };

    /**
     *  @brief get the private member Cosmology::m_Pk0_CLASS
     *
     *  @return the normalization of the power spectrum for CLASS [http://class-code.net/]
     */
    double Pk0_CLASS () const { return m_Pk0_CLASS; };

    /**
     *  @brief get the private member Cosmology::m_model
     *
     *  @return the cosmologial model used to compute distances
     */
    string model () const { return m_model; }

    /**
     *  @brief get the private member Cosmology::m_unit
     *
     *  @return unit: 0 &rarr; phyical units; 1 &rarr; cosmological units
     *  (i.e. without \e h)
     */
    bool unit () const { return m_unit; };

    /**
     *  @brief print the values of the private members on the screen
     *
     *  @return none
     */
    void print_parameters () const {
      cout << "Omega_matter = " << m_Omega_matter << endl;
      cout << "Omega_baryon = " << m_Omega_baryon << endl;
      cout << "Omega_neutrinos = " << m_Omega_neutrinos << endl;
      cout << "massless_neutrinos = " << m_massless_neutrinos << endl;
      cout << "massive_neutrinos = " << m_massive_neutrinos << endl;
      cout << "Omega_DE = " << m_Omega_DE << endl;
      cout << "Omega_k = " << m_Omega_k << endl; 
      cout << "Omega_CDM = " << m_Omega_CDM << endl;
      cout << "h = " << m_hh << endl;
      cout << "scalar_amp " << m_scalar_amp << endl;
      cout << "n_spec = " << m_n_spec << endl;
      cout << "w0 = " << m_w0 << endl;
      cout << "wa = " << m_wa << endl;
      cout << "fNL = " << m_fNL << endl;
      cout << "type_NG = " << m_type_NG << endl;
    }
    ///@}


    /**
     *  @name Functions to set the private members of the class
     */
    ///@{

    /**
     *  @brief set the value of one cosmological paramter
     *  
     *  @param parameter cosmological parameter to set
     *  @param value the new value for the parameter 
     *
     *  @return none
     */
    void set_parameter (const CosmoPar parameter, const double value);

    /**
     *  @brief set the value of some cosmological paramters 
     *  
     *  @param parameter vector containing the cosmological parameters
     *  to set
     *  @param value vector containing the new values for the
     *  parameters
     *
     *  @return none
     */
    void set_parameter (const vector<CosmoPar> parameter, const vector<double> value);

    /**
     *  @brief set the value of &Omega;<SUB>M</SUB>, keeping
     * &Omega;<SUB>DE</SUB>=1-&Omega;<SUB>M</SUB>
     *
     *  @param Omega_matter &Omega;<SUB>M</SUB>: density of baryons,
     *  cold dark matter and massive neutrinos (in units of the
     *  critical density)
     *
     *  @return none
     */
    void set_Omega (const double Omega_matter=0.27) {
      m_Omega_matter = Omega_matter; 
      m_Omega_DE = 1.-m_Omega_matter; 
      m_Omega_k = 1.-m_Omega_matter-m_Omega_radiation-m_Omega_DE;
      m_Omega_CDM = m_Omega_matter-m_Omega_baryon-m_Omega_neutrinos;
    };

    /**
     * @brief set the value of &Omega;<SUB>b</SUB>, keeping
     * &Omega;<SUB>CDM</SUB>=&Omega;<SUB>M</SUB>-&Omega;<SUB>b</SUB>
     *
     *  @param Omega_baryon &Omega;<SUB>b</SUB>: density of baryons,
     *  (in units of the critical density)
     *
     *  @return none
     */
    void set_OmegaB (const double Omega_baryon=0.046) {
      m_Omega_baryon = Omega_baryon; 
      m_Omega_CDM = m_Omega_matter-m_Omega_baryon-m_Omega_neutrinos;
    };

    /**
     *  @brief set the value of &Omega;<SUB>M</SUB>
     *
     *  @param Omega_matter &Omega;<SUB>M</SUB>: density of baryons,
     *  cold dark matter and massive neutrinos (in units of the
     *  critical density)
     *
     *  @return none
     */
    void set_OmegaM (const double Omega_matter=0.27) {
      m_Omega_matter = Omega_matter; 
      m_Omega_k = 1.-m_Omega_matter-m_Omega_radiation-m_Omega_DE;
      m_Omega_CDM = m_Omega_matter-m_Omega_baryon-m_Omega_neutrinos;
    };

    /**
     *  @brief set the value of &Omega;<SUB>DE</SUB>
     *
     *  @param Omega_DE &Omega;<SUB>DE</SUB>: density of dark energy
     *
     *  @return none
     */
    void set_OmegaDE (const double Omega_DE=0.73) {
      m_Omega_DE = Omega_DE; 
      m_Omega_k = 1.-m_Omega_matter-m_Omega_radiation-m_Omega_DE;
      m_Omega_CDM = m_Omega_matter-m_Omega_baryon-m_Omega_neutrinos;
    };

    /**
     *  @brief set the value of &Omega;<SUB>&nu;</SUB>
     *  @param Omega_neutrinos &Omega;<SUB>&nu;</SUB>: density of
     *  massive neutrinos
     *  @param massless_neutrinos N<SUB>eff</SUB>: the effective
     *  number (for QED + non-instantaneous decoupling)
     *  @param massive_neutrinos the number of degenerate massive
     *  neutrino species
     *  @return none
     */
    void set_OmegaNu (const double Omega_neutrinos=0., const double massless_neutrinos=3.04, const int massive_neutrinos=0) {
      m_Omega_neutrinos = Omega_neutrinos; 
      m_massless_neutrinos = massless_neutrinos;
      m_massive_neutrinos = massive_neutrinos; 
      m_Omega_CDM = m_Omega_matter-m_Omega_baryon-m_Omega_neutrinos;
    };

    /**
     *  @brief set the private member Cosmology::m_Omega_radiation
     *
     *  @param Omega_radiation &Omega;<SUB>rad</SUB>: the radiation density
     *  @return none
     */
    void set_Omega_radiation (double Omega_radiation) { 
      m_Omega_radiation = Omega_radiation;
      m_Omega_k = 1.-m_Omega_matter-m_Omega_radiation-m_Omega_DE;
    };

    /**
     *  @brief set the value of H<SUB>0</SUB>
     *
     *  @param H0 H<SUB>0</SUB>: Hubble constant [km/sec/Mpc]
     *
     *  @return none
     */
    void set_H0 (const double H0=70.) {
      m_hh = H0/100.; 
      m_H0 = (m_unit) ? 100. : m_hh*100.;   
      m_t_H = 1./m_H0; 
      m_D_H = par::cc*m_t_H;
    };

    /**
     *  @brief set the value of &sigma;<SUB>8</SUB>
     *
     *  @param sigma8 &sigma;<SUB>8</SUB>: power spectrum normalization
     *
     *  @return none
     */
    void set_sigma8 (const double sigma8=-1.) { m_sigma8 = sigma8; };

    /**
     *  @brief set the value of A<SUB>s</SUB>
     *
     *  @param scalar_amp A<SUB>s</SUB>: initial scalar amplitude of the power spectrum
     *
     *  @return none
     */
    void set_scalar_amp (const double scalar_amp=2.46e-9) { m_scalar_amp = scalar_amp; }; 

    /**
     *  @brief set the value of n<SUB>spec</SUB>
     *
     *  @param n_spec n<SUB>spec</SUB>: the primordial spectral index
     *
     *  @return none
     */
    void set_n_spec (const double n_spec) { m_n_spec = n_spec; }; 

    /**
     *  @brief set the value of w<SUB>0</SUB>
     *
     *  @param w0 w<SUB>0</SUB>: parameter of the dark energy equation
     *  of state (CPL parameterisation)
     *
     *  @return none 
     */
    void set_w0 (const double w0=-1.) { m_w0 = w0; };  
    
    /**
     *  @brief set the value of w<SUB>a</SUB>
     *
     *  @param wa w<SUB>a</SUB>: parameter of the dark energy equation of state (CPL parameterisation)
     *
     *  @return none
     */
    void set_wa (const double wa=0.) { m_wa = wa; };  

    /**
     *  @brief set the value of &rho;<SUB>0</SUB>
     *
     *  @param RhoZero: the mean density of the Universe at z=0
     *  [Msun*Mpc^-3]
     *
     *  @return none
     */
    void set_RhoZero (const double RhoZero=7.5e10) { m_RhoZero = RhoZero; };  
    
    /**
     *  @brief set the value of f<SUB>NL</SUB>
     *
     *  @param fNL f<SUB>NL</SUB>: the non-Gaussian amplitude
     *
     *  @return none
     */
    void set_fNL (const double fNL=0.) { m_fNL = fNL; };  
    
    /**
     *  @brief set the value of the non-Gaussian shape
     *
     *  @param type_NG the non-Gaussian shape (type=1 local, type=2
     *  equilateral, type=3 enfolded, type=4 orthogonal)
     *
     *  @return none
     */
    void set_type_NG (const double type_NG=1) { m_type_NG = type_NG; };  
    
    /**
     *  @brief set the cosmologial model used to compute distances
     *
     *  @param model the cosmologial model used to compute distances
     *
     *  @return none
     */
    void set_model (const string model="LCDM") { m_model = model; };  

    /**
     *  @brief set the value of unit
     *
     *  @param unit: 0 &rarr; phyical units; 1 &rarr; cosmological units (i.e. without \e h)
     *
     *  @return none
     */
    void set_unit (const bool unit=1) { m_unit = unit; }

    ///@}

    
    /**
     *  @name Functions to estimate general cosmological parameters
     */
    ///@{

    /**
     *  @brief the matter density at a given redshift
     *
     *  @param redshift redshift
     *
     *  @return &Omega;<SUB>M</SUB>
     */
    double OmegaM (const double redshift=0.) const; 

    /**
     *  @brief the dark energy density at a given redshift
     *
     *  @param redshift redshift
     *
     *  @return &Omega;<SUB>DE</SUB>
     */
    double OmegaDE (const double redshift=0.) const;

    /**
     *  @brief the radiation density at a given redshift
     *
     *  @param redshift redshift
     *
     *  @return &Omega;<SUB>rad</SUB>
     */
    double OmegaR (const double redshift=0.) const;

    /**
     *  @brief the density of curvature energy at a given redshift
     *
     *  @param redshift redshift
     *
     *  @return &Omega;<SUB>k</SUB>
     */
    double OmegaK (const double redshift=0.) const;

    /**
     *  @brief the cosmic density at a given redshift
     *
     *  @param redshift redshift
     *
     *  @return &Omega;
     */
    double Omega (const double redshift=0.) const;

    /**
     *  @brief the mean cosmic density
     *
     *  @param Omega_matter &Omega;<SUB>M</SUB>: density of baryons,
     *  cold dark matter and massive neutrinos (in units of the
     *  critical density)
     *
     *  @param Omega_neutrinos &Omega;<SUB>&nu;</SUB>: density of
     *  massive neutrinos
     *
     *  @param unit1 0 &rarr; phyical units; 1 &rarr; cosmological units
     *  (i.e. without \e h)
     *
     *  @return &rho;<SUB>mean</SUB>: the mean cosmic density [Msun*Mpc^-3(*h^2)]
     */
    double Rho (const double Omega_matter=0.27, const double Omega_neutrinos=0., const bool unit1=1) const;  

    /**
     *  @brief the overdensity within a sphere of radius R
     *
     *	this function returns &Delta;<SUB>R</SUB>, the overdensity
     *  within a sphere of radius R, as a function of the critical
     *  overdensity
     *
     *  @param Delta_crit &Delta;<SUB>crit</SUB>: critical overdensity
     *
     *  @param redshift redshift
     *
     *  @return &Delta;<SUB>R</SUB>: the overdensity within a sphere
     *  of radius R
     */
    double DeltaR (const double Delta_crit=200., const double redshift=0.) const;  

    /**
     *  @brief the DE equation of state in the CPL parameterisation,
     *  as a function of redshift
     *
     *  @param redshift redshift
     *
     *  @return w: the DE equation of state in the CPL
     *  parameterisation
     */
    double w_CPL (const double redshift=0.) const;

    /**
     *  @brief auxiliary function used to compute the Hubble function
     *
     *  this function returns f<SUB>DE</SUB>(z), a parameter that
     *  multiplies &Omega;<SUB>DE</SUB> in the Hubble function (see
     *  e.g. Bassett & Hlozek 2010)
     *
     *  @param redshift redshift
     *
     *  @return f<SUB>DE</SUB>
     */
    double f_DE (const double redshift=0.) const;
   
    /**
     *  @brief auxiliary function used to compute the Hubble function
     *  @param redshift redshift
     *  @return E=H/H<SUB>0</SUB> 
     */
    double EE (const double redshift=0.) const;
   
    /**
     *  @brief the Hubble function
     *  @param redshift redshift
     *  @return H
     */
    double HH (const double redshift=0.) const;
 
    /**
     *  @brief auxiliary function used to estimate the linear growth
     factor
     *  @param redshift redshift
     *  @return g
     */
    double gg (const double redshift=0.) const;   

    /**
     *  @brief the linear growth factor at a given redshift
     *  @param redshift redshift
     *  @return D
     */
    double DD (const double redshift=0.) const;   


    /**
     *  @brief lookback time at a given redshift
     *  @param redshift redshift
     *  @return t<SUB>lookback</SUB> [Gyr]
     */
    double lookback_time (const double redshift=0.) const; 

    /**
     *  @brief cosmic time at a given redshift
     *  @param redshift redshift
     *  @return t<SUB>cosmic</SUB> [Gyr]
     */
    double cosmic_time (const double redshift=0.) const; 

    /**
     *  @brief auxiliary function used to compute the deceleration
     *  parameter
     * 
     *  see e.g. de Araujo 2005
     *  
     *  @param redshift redshift
     *  @return E<SUB>2</SUB>
     */
    double EE2 (const double redshift=0.) const;
    
    /**
     *  @brief the deceleration parameter at a given redshift
     *  @param redshift redshift
     *  @return q
     */
    double qq (const double redshift=0.) const;

    /**
     *  @brief derivative of the Hubble function at a given redshift
     *  @param redshift redshift
     *  @return dH/dz
     */
    double Hdot (const double redshift=0.) const;

    /**
     *  @brief redshift at which the Universe begins to accelerate 
     *
     *  see e.g. de Araujo 2005
     *
     *  @return z<SUB>acc</SUB>
     */
    double z_acc () const; 

    /**
     *  @brief redshift of matter-dark energy equality
     *
     *  see e.g. de Araujo 2005
     *
     *  @return z<SUB>acc</SUB>
     */
    double z_eq () const;
  
    // Maximum absolute magnitude to have a volume limited catalogue
    
    /**
     *  @brief maximum absolute magnitude to have a volume-limited
     *	catalogue
     *
     *  @param z_max maximum redshift
     *  @param mag_lim magnitude limit
     *  @return z<SUB>acc</SUB>
     */
    double Mag_Volume_limited (const double z_max=1., const double mag_lim=-20.) const;

    /**
     *  @brief bolometric luminosity
     *  @param redshift redshift
     *  @param flux flux
     *  @return L<SUB>bol</SUB>
     */
    double Lum_bol (const double redshift=0., const double flux=1.) const;  

    /**
     *  @brief redshift at a given comoving distance
     *  
     *  this method provides the redshift for a given comoving
     *  distance
     *  @param d_c line-of-sight comoving distance
     *
     *  @param z1_guess minimum redshift used to search the redshift
     *
     *  @param z2_guess maximum redshift used to search the redshift
     *
     *  @param prec precision of the computation ( prec =
     *  min(prec,1.e-5) )
     *
     *  @return redshift
     */
    double Redshift (const double d_c=1., const double z1_guess=0., const double z2_guess=10., const double prec=0.0001) const; 

    /**
     *  @brief redshift at a given comoving distance
     *
     *  this method provides the redshift for a given comoving
     *  distance; the iteration process exploits the properties of the
     *  f=d(z) function (i.e. f'=1/E(z), f"<0)
     *     
     *  @author Mauro Roncarelli
     *  @author mauro.roncarelli@unibo.it
     *
     *  @param d_c line-of-sight comoving distance
     *
     *  @param z1_guess minimum redshift of the region explored to
     *  search the redshift
     *
     *  @param z2_guess maximum redshift used to search the redshift
     *
     *  @param go_fast 0 &rarr; the method uses the function
     *  Cosmology::D_C; 1 &rarr; the method uses the function
     *  Cosmology::D_C_LCDM (much faster than D_C)
     *
     *  @param prec precision of the computation; ( prec =
     *  min(prec,1.e-5) )
     *
     *  @return redshift
     *
     *  @warning if go_fast=1, the method works correctly only for a
     *  flat &Lambda;CDM universe at low enough redshifts where
     *  &Omega;<SUB>r</SUB> is negligible;
     *
     *  @warning wrong redshift interval limits do not lead to an
     *  error, but just slow down the computation
     */
    double Redshift_LCDM (const double d_c=1., const double z1_guess=0., const double z2_guess=10., const bool go_fast=1, const double prec=0.0001) const; 

    /**
     *  @brief redshift at a given wf
     *
     *  this routine estimates the redshift from wf, given the parent
     *  halo mass at z=z', z', and its assembled fraction f
     *
     *  @author Carlo Giocoli
     *  @author cgiocoli@gmail.com
     *  @param mm mass
     *  @param redshift redshift
     *  @param ff assembled fraction
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass) const; valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param wwf rescaled variable w as in Lacey and Coles 1993
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @return redshift
     */
    double Redshift (const double, const double, const double, const string, const double, const string output_root="test") const; 

    /**
     *  @brief redshift at a given cosmic time
     *  @param time cosmic time
     *  @param z1_guess minimum redshift of the region explored to search the redshift
     *  @param z2_guess maximum redshift used to search the redshift
     *  @return redshift
     */
    double Redshift_time (const double, const double, const double) const;

    /**
     *  @brief spherical collapse density threshold at a given
     *  redshift
     *
     *  by Nakamura & Suto (1997)
     *
     *  @param redshift redshift
     *  @return &delta;<SUB>c</SUB>
     */
    double deltac (const double) const; 
 
    /**
     *  @brief virial overdensity
     *  @author Carlo Giocoli
     *  @author cgiocoli@gmail.com
     *  @param redshift redshift
     *  @return &Delta;<SUB>vir</SUB>
     */
    double Deltavir (const double redshift) const; 

    ///@}
    

    /**
     *  @name Functions to estimate cosmological distances and volumes
     */
    ///@{

    /**
     *  @brief the comoving line-of-sight distance at a given redshift
     *
     *  for demonstration, see \ref distances.cpp
     *
     *  @param redshift redshift
     *  @return D<SUB>C</SUB>
     */
    double D_C (const double) const;  

    /**
     *  @brief the comoving line-of-sight distance at a given redshift
     *
     *  this method provides the comoving distance for a given
     *  redshift using elliptic integrals (see Numer. Math. (2010)
     *  116:687–719, T. Fukushima "Fast computation of incomplete
     *  elliptic integral of first kind by half argument
     *  transformation")
     *
     *  d(z) is computed using the following formula, taken from
     *  A. Meszaros & J. Ripa 2013, arXiv:1306.4736,
     *
     *  d(z) = C*[F(&phi;<SUB>0</SUB>,m)-F(&phi;<SUB>1</SUB>,m)]
     *
     *  where
     *
     *  C = 3<SUP>-1/4</SUP> /
     *  (&Omega;<Sub>M</SUB><SUP>1/3</SUP>*&Omega;<SUB>&Lambda;</SUB><SUP>1/6</SUP>)
     *
     *  &phi;<SUB>0</SUB> = arccos(
     *  [1+(1-3<SUP>0.5</SUP>)*(&Omega;<SUB>&Lambda;</SUB>/&Omega;<Sub>M</SUB>)<SUP>1/3</SUP>]
     *  /
     *  [1+(1+3<SUP>0.5</SUP>)*(&Omega;<SUB>&Lambda;</SUB>/&Omega;<Sub>M</SUB>)<SUP>1/3</SUP>]
     *  )
     *
     *  &phi;<SUB>1</SUB> = arccos(
     *  [1+(1-3<SUP>0.5</SUP>)*(&Omega;<SUB>&Lambda;</SUB>/&Omega;<Sub>M</SUB>)<SUP>1/3</SUP>*(1+z)]
     *  /
     *  [1+(1+3<SUP>0.5</SUP>)*(&Omega;<SUB>&Lambda;</SUB>/&Omega;<Sub>M</SUB>)<SUP>1/3</SUP>*(1+z)]
     *  )
     *
     *  F(&phi;, m) is the incomplete elliptic integral of the first
     *  kind, with argument &phi;
     *
     *  m = (2+3<SUP>0.5</SUP>)/4 is the elliptic integral shape parameter
     *  
     *  for demonstration, see \ref distances.cpp
     *
     *  @author Mauro Roncarelli
     *  @author mauro.roncarelli@unibo.it
     *
     *  @param redshift redshift
     *  @return D<SUB>C</SUB>
     *
     *  @warning this method works only for a flat &Lambda;CDM
     *  universe at low enough redshifts where &Omega;<SUB>r</SUB> is
     *  negligible; it does not work for non-standard dark energy or
     *  non-flat models
     */
    double D_C_LCDM (const double) const;  


    // table of redshift -- comoving line-of-sight distance

    /**
     *  @brief create a table of [redshift, comoving line-of-sight
     * distance]
     *
     *  this function is used to create a table of [redshift, comoving
     *  line-of-sight distance], useful to speed up the analysis
     *
     *  @param [in] file_table name of the file where the table is stored
     *  @param [in] z_min minimum redshift of the table
     *  @param [in] z_max maximum redshift of the table
     *  @param [in] step redshift step
     *  @param [out] Redshift vector of redshifts
     *  @param [out] dc vector of comoving line-of-sight distances
     *  @return none
     */
    void D_C_table (const string, const double, const double, const int, vector<double> &, vector<double> &) const;

    /**
     *  @brief the comoving transverse distance at a given redshift
     *  @param redshift redshift
     *  @return D<SUB>M</SUB>
     */
    double D_M (const double) const;

    /**
     *  @brief the angular diameter distance at a given redshift
     *  @param redshift redshift
     *  @return D<SUB>A</SUB>
     */
    double D_A (const double) const; 
  
    /**
     *  @brief the luminosity distance at a given redshift
     *  @param redshift redshift
     *  @return D<SUB>L</SUB>
     */
    double D_L (const double) const; 
  
    /**
     *  @brief the average distance at a given redshift, used to
     rescale the correlation function
     *  @param redshift redshift
     *  @return D<SUB>V</SUB>
     */
    double D_V (const double) const;
  
    /**
     *  @brief the distance at a given redshift. Distance available are:
     D<SUB>C</SUB>,D<SUB>L</SUB>,D<SUB>A</SUB>,D<SUB>V</SUB>,
     D<SUB>V</SUB>/r<SUB>s</SUB>, r<SUB>s</SUB>/D<SUB>V</SUB>
     *  
     *  @author Alfonso Veropalumbo
     *  @author alfonso.veropalumbo@unibo.it
     *  @param redshift redshift
     *  @param distance_type the type of distance to return 
     *  @return Distance
     */
    double Distance (const double, const string) const;

    /**
     *  @brief comoving volume for a given redshift range and sky area
     *
     *  this function provides an approximated expression valid only
     *  for a LCDM model
     *
     *  @param z1 minimum redshift
     *  @param z2 maximum redshift
     *  @param Area sky area
     *  @return comoving volume
     */
    double Volume (const double, const double, const double) const; 

    /**
     *  @brief total comoving volume from z=0 to z
     *
     *  from Hogg 2000, Eq. 29
     *
     *  @author Alfonso Veropalumbo
     *  @author alfonso.veropalumbo@unibo.it
     *  @param zz redshift
     *  @return comoving volume
     */
    double Volume (const double) const;
  
    /**
     *  @brief maximum redshift for a given volume, sky area and minimum redshift
     *  @param Volume volume
     *  @param Area sky area
     *  @param z_min minimum redshift
     *  @return redshift
     */
    double max_redshift (const double, const double, const double) const; 

    /**
     *  @brief the derivative of the comoving volume,
     d<SUP>2</SUP>V/(dz*d&Omega;) at a given redshift
     *  @param redshift redshift
     *  @param angle_rad 0 &rarr; &Omega; in square degrees; 1 &rarr; &Omega;
     *  in steradians
     *
     *  @return d<SUP>2</SUP>V/(dz*d&Omega;)
     */
    double dV_dZdOmega (const double, const bool) const; 

    ///@}


    /**
     *  @name Functions to estimate the mass function
     */
    ///@{

    /**
     *  @brief the mass function of dark matter haloes (filaments and
     * sheets)
     *
     *  @param Mass mass
     *
     *  @param redshift redshift
     *
     *  @param author_MF author(s) who proposed the mass function; valid
     *  authors are: PS (Press & Schechter), ST (Sheth & Tormen),
     *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
     *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
     *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
     *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
     *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
     *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
     *  al. 2012)
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param Delta &Delta;: the overdensity, defined as the mean
     *  interior density relative to the background
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr; the
     *  GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return the mass function, d&Phi;/dM=dn(M)/dM
     *
     *  @warning the input parameter &Delta; is used only in the
     *  Tinker mass function; it is the background overdensity, not
     *  the critical overdensity &Delta;<SUB>crit</SUB>; the function
     *  cosmobl::Cosmology::DeltaR can be used to convert
     *  &Delta;<SUB>crit</SUB> into &Delta;
     */
    double mass_function (const double Mass, const double redshift, const string author_MF, const string method_SS, const string output_root="test", const double Delta=200., const string interpType="Linear", const int Num=-1, const double stepsize=100., const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    /**
     *  @brief the mass function of dark matter haloes (filaments and
     *  sheets) computed quickly using a grid
     *
     *  @param Mass mass
     *
     *  @param redshift redshift
     *
     *  @param author_MF author(s) who proposed the mass function; valid
     *  authors are: PS (Press & Schechter), ST (Sheth & Tormen),
     *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
     *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
     *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
     *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
     *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
     *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
     *  al. 2012)
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param Delta &Delta;: the overdensity, defined as the mean
     *  interior density relative to the background
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr; the
     *  GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return the mass function, d&Phi;/dM=dn(M)/dM
     *
     *  @warning the input parameter &Delta; is used only in the
     *  Tinker mass function; it is the background overdensity, not
     *  the critical overdensity &Delta;<SUB>crit</SUB>; the function
     *  cosmobl::Cosmology::DeltaR can be used to convert
     *  &Delta;<SUB>crit</SUB> into &Delta;
     */
    double mass_function_fast (const double Mass, const double redshift, const string author_MF, const string method_SS, const string output_root="test", const double Delta=200., const string interpType="Linear", const int Num=-1, const double stepsize=100., const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString); 

    /**
     *  @brief the mass function of dark matter haloes (filaments and
     *  sheets) computed quickly passing directly the mass variance
     *  and its derivative as inputs
     *
     *  @param Mass mass
     *
     *  @param Sigma &sigma;(mass): the mass variance
     *
     *  @param Dln_Sigma dln&sigma;/dM: the derivative of the mass
     *  variance
     *
     *  @param redshift redshift
     *
     *  @param author_MF author(s) who proposed the mass function;
     *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
     *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
     *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
     *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
     *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
     *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
     *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
     *  al. 2012)
     * 
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param Delta &Delta;: the overdensity, defined as the mean
     *  interior density relative to the background
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1
     *  &rarr; normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return the mass function, d&Phi;/dM=dn(M)/dM
     *
     *  @warning the input parameter &Delta; is used only in the
     *  Tinker mass function; it is the background overdensity, not
     *  the critical overdensity &Delta;<SUB>crit</SUB>; the function
     *  cosmobl::Cosmology::DeltaR can be used to convert
     *  &Delta;<SUB>crit</SUB> into &Delta;
     */
    double mass_function (const double, const double, const double, const double, const string, const string output_root="test", const double Delta=200., const string interpType="Linear", const int Num=-1, const double stepsize=100., const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string method_SS="CAMB", const string file_par=par::defaultString); 
    
    /**
     *  @brief number of dark matter haloes per steradian or square
     *  degree, for a given redshift range
     *
     *  @param Mass_min minimum mass
     *
     *  @param Mass_max maximum mass
     *
     *  @param z_min minimum redshift
     *
     *  @param z_max maximum redshift
     *
     *  @param angle_rad 0 &rarr; &Omega; in square degrees; 1 &rarr; &Omega;
     *  in steradians
     *
     *  @param author_MF author(s) who proposed the mass function;
     *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
     *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
     *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
     *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
     *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
     *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
     *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
     *  al. 2012)
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return n<SUB>haloes</SUB>: the number density of dark matter
     *  haloes (per steradian or square degree)
     */
    double n_haloes (const double, const double, const double, const double, const bool, const string, const string, const string output_root="test", const string interpType="Linear", const int Num=-1, const double stepsize=100., const double k_max=100., const string file_par=par::defaultString);

    /**
     *  @brief minimum halo mass, given the number of haloes in a
     *  given region of sky
     *
     *  @param n_halo number density of dark matter haloes
     *
     *  @param Area sky area
     *
     *  @param angle_rad 0 &rarr; &Omega; in square degrees; 1 &rarr; &Omega;
     *  in steradians
     *
     *  @param z_min minimum redshift
     * 
     *  @param z_max maximum redshift
     *
     *  @param Mmax maximum mass
     *  
     *  @param lgM1_guess logarithm of the minimum mass used by the
     *  root finder
     *
     *  @param lgM2_guess logarithm of the maximum mass used by the
     *  root finder
     *  
     *  @param author_MF author(s) who proposed the mass function;
     *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
     *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
     *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
     *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
     *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
     *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
     *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
     *  al. 2012)
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return minimum halo mass
     */
    double MhaloMin (const int, const double, const bool, const double, const double, const double, const double, const double, const string, const string, const string output_root="test", const string interpType="Linear", const int Num=-1, const double stepsize=100., const double k_max=100., const string file_par=par::defaultString) const;

    ///@}


    /**
     *  @name Functions to estimate the cosmic mass accretion history
     */
    ///@{

    /**
     *  @brief differential distribution 
     *
     *  this function provides the differential rescaled and
     *  generalized formation redshift distribution
     *
     *  @author Carlo Giocoli
     *  @author cgiocoli@gmail.com
     *  @param ww rescaled variable w as in Lacey and Coles 1993
     *  @param ff assembled fraction
     *  @param author valid authors are: NS (Nusser and Sheth), GTS
     *  (Giocoli et al. 2012)
     *  @return p(w): differential distribution 
     */
    double pw (const double, const double, const string) const; 

    /**
     *  @brief formation probability
     *
     *  this function provides the probability that a halo of a given
     *  mass m0 at redshift z0 make a mass fraction f at redshift z
     *
     *  @author Carlo Giocoli
     *  @author cgiocoli@gmail.com
     *  @param m0 halo mass
     *  @param z0 redshift when the halo has a mass m0
     *  @param frac mass fraction
     *  @param redshift redshift
     *
     *  @param author_model valid authors are: NS (Nusser and Sheth),
     *  GTS (Giocoli et al. 2012)
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *     
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @return p(z): formation probability
     */
    double pz (const double, const double, const double, const double, const string, const string, const string output_root="test") const; 

    /**
     *  @brief cumulative distribution 
     *
     *  this function provides the cumulative rescaled and generalized
     *  formation redshift distribution
     *
     *  @author Carlo Giocoli
     *  @author cgiocoli@gmail.com
     *  @param ww rescaled variable w as in Lacey and Coles 1993
     *  @param ff assembled fraction
     *  @param author valid authors are: NS (Nusser and Sheth), GTS
     *  (Giocoli et al. 2012)
     *  @return P(w): cumulative distribution
     */
    double cumPw (const double, const double, const string) const; 
    
    /**
     *  @brief median formation w     
     *  @author Carlo Giocoli
     *  @author cgiocoli@gmail.com
     *  @param [in] ff assembled fraction
     *  @param [in] author_model valid authors are: NS (Nusser and Sheth), GTS
     *  (Giocoli et al. 2012)
     *  @param [out] wf vector of w(f)
     *  @return none
     */
    void medianwf (const double, const string, vector<double> &) const; 

    /**
     *  @brief median formation z 
     *  @author Carlo Giocoli
     *  @author cgiocoli@gmail.com
     *  @param [in] ff assembled fraction
     *  @param [in] mass halo mass
     *  @param [in] z0 redshift when the halo has a mass mass
     *
     *  @param [in] author_model valid authors are: NS (Nusser and
     *  Sheth), GTS (Giocoli et al. 2012)
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param [out] zf vector of z(f)
     *      
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @return none
     */
    void medianzf (const double, const double, const double, const string, const string, vector<double> &, const string output_root="test") const; 
  
    /**
     *  @brief rescaled variable w as in Lacey and Coles 1993
     *
     *  this functions provides the conditional variable
     *  w=[&delta;<SUB>c</SUB>(zf) - &delta;<SUB>c</SUB>(z)] / &radic;
     *  [s(fm)-s(m)] where&delta;<SUB>c</SUB>(z) =
     *  &delta;<SUB>c0</SUB>(z)/D+(z)
     *
     *  @author Carlo Giocoli
     *  @author cgiocoli@gmail.com
     *  @param mm halo mass
     *  @param redshift redshift 
     *  @param ff assembled fraction
     *  @param zf redshift at which the mass is accreted
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *     
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @return the conditional variable w
     */
    double wf (const double, const double, const double, const double, const string, const string output_root="test") const; 

    /**
     *  @brief the unevolved mass function
     *  @author Carlo Giocoli
     *  @author cgiocoli@gmail.com
     *  @param mass_accr mass accreted
     *  @return the unevolved mass function
     */
    double unevolved_mass_function (const double) const; 

    /**
     *  @brief compute te halo concentration
     *  @author Carlo Giocoli
     *  @author cgiocoli@gmail.com
     *  @param Vmax V<SUB>max</SUB>
     *  @param Rmax R<SUB>max</SUB>
     *  @return the halo concentration
     */
    double concentration (const double, const double) const; 

    ///@}


    /**
     *  @name Functions to estimate the power spectrum and related quantities
     */
    ///@{

    /**
     *  @brief amplitude of the curvature perturbations
     *  
     *  this function provides an approximate value of A<SUB>s</SUB>,
     *  for a fiven value of &sigma;<SUB>8</SUB>, from Vikhlinin et
     *  al. 2009, Eq.3; it is valid only without massive neutrinos! 
     *  (see also Hu & Jain 2004)
     *
     *  @param sigma8 &sigma;<SUB>8</SUB>: the power spectrum normalization
     *  @return A<SUB>s</SUB>
     */
    double As (const double) const; 

    /**
     *  @brief unnormalized power spectrum
     *  @param kk wave vector module
     *  @param redshift redshift
     *
     *  @param method_Pk method used to compute the power spectrum; at
     *  the present state the only permitted value is "EisensteinHu"
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @return P<SUB>unnormalized</SUB>
     */
    double Pk_UnNorm (const double, const double, const string) const; // unnormalized power spectrum

    /**
     *  @brief return the path to the power spectrum output
     *
     *  @param code method used to compute the power spectrum
     * 
     *  @param NL 0 &rarr; linear power spectrum; 1 &rarr; non-linear
     *  power spectrum
     * 
     *  @param redshift the redshift
     *
     *  @param run true &rarr; write or read the table where the dark
     *  matter power spectrum is stored
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param k_max the maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will use be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return the path to the power spectrum output
     */
    string Pk_output_file (const string code, const bool NL, const double redshift, const bool run=0, const string output_root="test", const double k_max=100., const string file_par=par::defaultString);

    /**
     *  @brief run CAMB [http://camb.info/]
     *  
     *  this function runs CAMB [http://camb.info/], after editing the parameter file
     *  appropriately (if file_par=NULL)
     *
     *  @param NL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power
     *  spectrum
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will use be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return none
     */
    void run_CAMB (const bool, const double, const string output_root="test", const double k_max=100., const string file_par=par::defaultString) const; 

    /**
     *  @brief write or read the table where the dark matter power
     *  spectrum is stored
     *
     *  @param [in] code method used to compute the power spectrum;
     *  valid codes are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param [in] NL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power spectrum
     *
     *  @param [out] lgkk vector of log(k)
     *
     *  @param [out] lgPk vector of log(P(k))
     *
     *  @param [in] redshift redshift
     *
     *  @param [in] output_root output_root of the parameter file used
     *  to compute the power spectrum; it can be any name
     *
     *  @param [in] k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param [in] file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return none
     */
    void Table_PkCodes (const string code, const bool NL, vector<double> &lgkk, vector<double> &lgPk, const double redshift, const string output_root="test", const double k_max=100., string file_par=par::defaultString) const;

    /**
     *  @brief write or read the table where the dark matter two-point
     *  correlation function is stored
     *
     *  @param [in] code method used to compute the power spectrum;
     *  valid codes are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param [in] NL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power spectrum
     *
     *  @param [out] rr vector of comoving separations
     *
     *  @param [out] xi vector of the binned values of &xi;(r)
     *
     *  @param [in] redshift redshift
     *
     *  @param [in] output_root output_root of the parameter file used
     *  to compute the power spectrum; it can be any name
     *
     *  @param [in] k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param [in] file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return none
     */
    void Table_XiCodes (const string code, const bool NL, vector<double> &rr, vector<double> &xi, const double redshift, const string output_root, const double k_max, string file_par) const;

    /**
     *  @brief normalization of the power spectrum
     *
     *  this function sets the value of the private member m_Pk0_*,
     *  i.e. the normalization of the power spectrum
     *
     *  @param method_Pk method used to compute the power spectrum;
     *  valid choices for method_Pk are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr; the
     *  GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return none
     */
    void Pk_0 (const string, const double, const string output_root="test", const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString); 

    /**
     *  @brief normalized power spectrum
     *
     *  this function provides the power spectrum P(k); it can use
     *  either CAMB, CLASS, MPTbreeze or the analytic
     *  approximation by Eisenstein & Hu
     *
     *  @param kk the wave vector module
     *
     *  @param method_Pk method used to compute the power spectrum; valid
     *  choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param NL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power
     *  spectrum
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr; the
     *  GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return P(k)
     */
    double Pk (const double, const string, const bool, const double, const string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString); 

    /**
     *  @brief return power spectrum first three multipoles using linear kaiser model
     *
     *  this function provides the power spectrum P(k); it can use
     *  either CAMB, CLASS, MPTbreeze or the analytic
     *  approximation by Eisenstein & Hu
     *
     *  @param method_Pk method used to compute the power spectrum; valid
     *  choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param Pk0 the monopole of the power spectrum
     *
     *  @param Pk2 the quadrupole of the power spectrum
     *
     *  @param Pk4 the hexadecapole of the power spectrum
     *
     *  @param kk the scale at which compute the power spectrum
     *
     *  @param NL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power
     *  spectrum
     *
     *  @param redshift redshift
     * 
     *  @param bias the linear bias
     *
     *  @param sigma_NL the BAO damping parameter
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr; the
     *  GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return none
     */
    void Pk_Kaiser_multipoles (vector<double> &, vector<double> &, vector<double> &, const vector<double>, const string, const bool, const double, const double, const double sigma_NL = 0., const string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString); 

    /**
     *  @brief the dark matter power spectrum, de-wiggled (see e.g. Anderson et al 2014)
     *
     *  this function provides the De-Wiggled dark matter power spectrum,
     *
     *  @author Alfonso Veropalumbo
     *  @author alfonso.veropalumbo@unibo.it
     *  
     *  @param kk the wave vector module
     *
     *  @param redshift redshift
     *
     *  @param sigma_NL the non linear BAO damping
     *
     *  @param output_root output_root of the parameter file used to compute
     *  the power spectrum and &sigma;(mass); it can be any name
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1
     *  &rarr; normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param aa parameter \e a of Eq. 24 of Anderson et al. 2012
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @return P;<SUB>DW</SUB>(k): the De-Wiggled power
     *  spectrum of dark matter
     */
    double Pk_DeWiggle (const double , const double , const double , const string output_root = "test", const bool norm=1, const double k_min=0., const double k_max=100., const double aa=1., const double prec=1.e-2);

    /**
     *  @brief the unnormalized mass variance, &sigma;<SUP>2</SUP>(R)
     *
     *  @param RR the radius R
     *
     *  @param method_Pk method used to compute the power spectrum; valid
     *  choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return unnormalized &sigma;<SUP>2</SUP>(R)
     */
    double SSR (const double, const string, const double, const string output_root="test", const double k_max=100., const string file_par=par::defaultString) const; 
 
    /**
     *  @brief the mass variance, &sigma;<SUP>2</SUP>(R)
     *
     *  @param RR the radius R
     *
     *  @param method_Pk method used to compute the power spectrum;
     *  valid choices for method_Pk are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return unnormalized &sigma;<SUP>2</SUP>(R)
     */
    double SSR_norm (const double, const string, const double, const string output_root="test", const double k_max=100., const string file_par=par::defaultString) const; 

    /**
     *  @brief the unnormalized mass variance, &sigma;<SUP>2</SUP>(M)
     *
     *  @param MM the mass M
     *
     *  @param method_Pk method used to compute the power spectrum;
     *  valid choices for method_Pk are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return unnormalized &sigma;<SUP>2</SUP>(M)
     */
    double SSM (const double, const string, const double, const string output_root="test", const double k_max=100., const string file_par=par::defaultString) const; 

    /**
     *  @brief the mass variance, &sigma;<SUP>2</SUP>(M)
     *
     *  @param MM the mass M
     *
     *  @param method_Pk method used to compute the power spectrum; valid
     *  choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return &sigma;<SUP>2</SUP>(M)
     */
    double SSM_norm (const double, const string, const double, const string output_root="test", const double k_max=100., const string file_par=par::defaultString) const; 
 
    /**
     *  @brief the nth-order derivative of the mass variance,
     * d<SUP>n</SUP>&sigma;<SUP>2</SUP>(R)/dR<SUP>n</SUP>
     *
     *  @param nd the derivative order n
     *
     *  @param RR the radius R
     *
     *  @param method_Pk method used to compute the power spectrum; valid
     *  choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return d<SUP>n</SUP>&sigma;<SUP>2</SUP>(R)/dR<SUP>n</SUP>
     */
    double dnSR (const int nd, const double RR, const string method_Pk, const double redshift, const string output_root="test", const string interpType="Linear", const int Num=-1, const double stepsize=100., const double k_max=100., const string file_par=par::defaultString) const; 

    /**
     *  @brief the derivative of the mass variance,
     *  d&sigma;<SUP>nd</SUP>(M)/dM
     *
     *  @param nd order of the derivative (d&sigma;<SUP>nd</SUP>(M)/dM)
     * 
     *  @param MM the mass M
     *
     *  @param method_Pk method used to compute the power spectrum; valid
     *  choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return d&sigma;<SUP>2</SUP>(M)/dM
     */
    double dnSM (const int nd, const double MM, const string method_Pk, const double redshift, const string output_root="test", const string interpType="Linear", const int Num=-1, const double stepsize=100., const double k_max=100., const string file_par=par::defaultString) const; 

    ///@}


    /**
     *  @name Functions to estimate the two-point correlation function, bias and related quantities
     */
    ///@{

    /**
     *  @brief the dark matter two-point correlation function
     *
     *  this function provides the dark matter correlation function,
     *  obtained by Fourier transforming the matter power spectrum
     *
     *  @param rr the module of the comoving separation
     *
     *  @param method_Pk method used to compute the power spectrum; valid
     *  choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param NL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power
     *  spectrum
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1
     *  &rarr; normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param aa parameter \e a of Eq. 24 of Anderson et al. 2012
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return &xi;<SUB>DM</SUB>(r): the spherically
     *  averaged (monopole) of the two-point correlation function of
     *  dark matter
     */
    double xi_DM (const double, const string, const double, const string output_root="test", const bool NL=1, const int norm=-1, const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=0, const double prec=1.e-2, const string file_par=par::defaultString);

    /**
     *  @brief the dark matter two-point correlation function, de-wiggled (see e.g. Anderson et al 2014)
     *
     *  this function provides the dark matter correlation function,
     *  obtained by Fourier transforming the De-Wiggled matter power spectrum
     *
     *  @author Alfonso Veropalumbo
     *  @author alfonso.veropalumbo@unibo.it
     *  
     *  @param rr the module of the comoving separation
     *
     *  @param redshift redshift
     *
     *  @param sigma_NL the non linear BAO damping
     *
     *  @param output_root output_root of the parameter file used to compute
     *  the power spectrum and &sigma;(mass); it can be any name
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1
     *  &rarr; normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param aa parameter \e a of Eq. 24 of Anderson et al. 2012
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @return &xi;<SUB>DW</SUB>(r): the De-Wiggled spherically
     *  averaged (monopole) of the two-point correlation function of
     *  dark matter
     */
    double xi_DM_DeWiggle (const double , const double , const double , const string output_root = "test", const bool norm=1, const double k_min=0., const double k_max=100., const double aa=1., const double prec=1.e-2);

    /**
     *  @brief get the dark matter two-point correlation function
     *
     *  this function provides the dark matter correlation function,
     *  obtained by Fourier transforming the matter power spectrum
     *
     *  @param [out] rr vector of r, the module of the comoving
     *  separation
     *
     *  @param [out] Xi vector of &xi;(r), the two-point correlation
     *  function of dark matter
     *
     *  @param [in] method_Pk method used to compute the power spectrum;
     *  valid choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param [in] redshift redshift
     *
     *  @param [in] output_root output_root of the parameter file used
     *  to compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param [in] xiType 0 &rarr; standard; 1 &rarr; Chuang & Wang model
     *
     *  @param [in] k_star k<SUB>*</SUB> of the Chuang & Wang model
     *
     *  @param [in] xiNL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power spectrum
     *
     *  @param [in] norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param [in] r_min minimum separation up to which the
     *  correlation function is computed
     *
     *  @param [in] r_max maximum separation up to which the
     *  correlation function is computed
     *
     *  @param [in] k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param [in] k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param [in] aa parameter \e a of Eq. 24 of Anderson et
     *  al. 2012
     *
     *  @param [in] GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param [in] prec accuracy of the GSL integration
     *
     *  @param [in] file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, const ignoring the cosmological parameters of the object
     *
     *  @return none
     */
    void get_xi (vector<double> &, vector<double> &, const string, const double, const string output_root="test", const bool xiType=0, const double k_star=-1., const bool xiNL=0, const int norm=-1, const double r_min=0.1, const double r_max=150., const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);
  
    /**
     *  @brief get the barred dark matter correlation functions
     *
     *  this function provides the dark matter \e barred correlation
     *  functions, used to model the two-point correlation function in redshift-space
     *
     *  @param [in] rr vector of r, the module of the comoving
     *  separation
     *
     *  @param [in] Xi vector of &xi;(r), the two-point correlation
     *  function of dark matter
     *
     *  @param [out] Xi_ vector of barred &xi;(r),
     *
     *  @param [out] Xi__ vector of double-barred &xi;(r)
     *
     *  @param [in] method_Pk method used to compute the power spectrum;
     *  valid choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param [in] redshift redshift
     *
     *  @param [in] xiType 0 &rarr; standard; 1 &rarr; Chuang & Wang model
     *
     *  @param [in] k_star k<SUB>*</SUB> of the Chuang & Wang model
     *
     *  @param [in] xiNL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power spectrum
     *
     *  @param [in] norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param [in] r_min minimum separation up to which the
     *  correlation function is computed
     *
     *  @param [in] r_max maximum separation up to which the
     *  correlation function is computed
     *
     *  @param [in] k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param [in] k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param [in] aa parameter \e a of Eq. 24 of Anderson et
     *  al. 2012
     *
     *  @param [in] GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param [in] prec accuracy of the GSL integration
     *
     *  @param [in] file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return none
     */
    void get_barred_xi (vector<double> rr, vector<double> Xi, vector<double> &Xi_, vector<double> &Xi__, const string method_Pk, const double redshift, const bool xiType=0, const double k_star=-1., const bool xiNL=0, const int norm=-1, const double r_min=0.1, const double r_max=150., const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString) const;

    /**
     *  @brief the dark matter projected correlation function
     *
     *  this function provides the dark matter projected correlation
     *  functions, obtained by Fourier transforming the matter power
     *  spectrum
     *
     *  @param rp r<SUB>p</SUB>: projected separation
     *
     *  @param method_Pk method used to compute the power spectrum; valid
     *  choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param NL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power
     *  spectrum
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param r_min minimum separation up to which the
     *  correlation function is computed
     *
     *  @param r_max maximum separation up to which the
     *  correlation function is computed
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param aa parameter \e a of Eq. 24 of Anderson et
     *  al. 2012
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return w<SUB>p,DM</SUB>(&theta;): the projected correlation
     *  function of dark matter
     */
    double wp_DM (const double, const string, const double, const string output_root="test", const bool NL=1, const int norm=-1, const double r_min=1.e-3, const double r_max=350., const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    /**
     *  @brief the k<SUB>*</SUB> parameter 
     *
     *  this function provides the k<SUB>*</SUB> parameter used to
     *  model the BAO (see e.g. Chuang & Wang 2012, Crocce &
     *  Scoccimarro2006, Matsubara 2008)
     *
     *  @param method_Pk method used to compute the power spectrum; valid
     *  choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return k<SUB>*</SUB>
     */
    double k_star (const string, const double, const string output_root="test", const double k_max=100., const string file_par=par::defaultString) const; 


    /**
     *  @brief the dark matter rms mass fluctuation
     *
     *  @param RR radius inside which the dark matter rms mass
     *  fluctuation is computed
     *
     *  @param corrType 0 &rarr; the projected correlation function,
     *  w(&theta;), is used; 1 &rarr; the spherically averaged correlation
     *  function, &xi;(r), is used
     *
     *  @param method_Pk method used to compute the power spectrum; valid
     *  choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *  
     *  @param NL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power
     *  spectrum
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param r_min minimum separation up to which the
     *  correlation function is computed
     *
     *  @param r_max maximum separation up to which the
     *  correlation function is computed
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param aa parameter \e a of Eq. 24 of Anderson et
     *  al. 2012
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return &sigma;<SUB>R</SUB>: the dark matter rms mass fluctuation
     */
    double sigmaR_DM (const double, const int, const string, const double, const string output_root="test", const bool NL=1, const int norm=-1, const double r_min=1.e-3, const double r_max=350., const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString); 

    /**
     *  @brief the dark matter rms mass fluctuation within 8 Mpc/h
     *
     *  this function provides the rms mass fluctuation within 8
     *  Mpc/h, estimated directly from the power spectrum
     *
     *  @param method_Pk method used to compute the power spectrum; valid
     *  choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *  
     *  @param NL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power
     *  spectrum
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return &sigma;<SUB>8</SUB>: the dark matter rms mass
     *  fluctuation within 8 Mpc/h
     */
    double sigma8_Pk (const string, const double, const string output_root="test", const bool NL=0, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString) const; 

    /**
     *  @brief bias of dark matter haloes
     *
     *  @param Mass halo mass
     *
     *  @param redshift redshift
     *
     *  @param author author(s) who proposed the bias; valid authors
     *  are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen
     *  2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction
     *  of Warren 2004), Tinker (Tinker et al. 2010)
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param Delta &Delta;: the overdensity, defined as the mean
     *  interior density relative to the background
     *  
     *  @param kk wave vector module
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return b<SUB>halo</SUB>: the dark matter bias
     *
     *  @warning the input parameter &Delta; is used only in the
     *  Tinker mass function; it is the background overdensity, not
     *  the critical overdensity &Delta;<SUB>crit</SUB>; the function
     *  cosmobl::Cosmology::DeltaR can be used to convert
     *  &Delta;<SUB>crit</SUB> into &Delta;
     */
    double bias_halo (const double, const double, const string, const string, const string output_root="test", const double Delta=200., const double kk=-1., const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString); 

    /**
     *  @brief bias of dark matter haloes
     *
     *  @param Mass halo mass
     *
     *  @param Sigma &sigma;(mass): the mass variance
     *
     *  @param redshift redshift
     *
     *  @param author_bias author(s) who proposed the bias; valid authors
     *  are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen
     *  2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction
     *  of Warren 2004), Tinker (Tinker et al. 2010)
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param Delta &Delta;: the overdensity, defined as the mean
     *  interior density relative to the background
     *  
     *  @param kk wave vector module
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return b<SUB>halo</SUB>: the dark matter bias
     *
     *  @warning the input parameter &Delta; is used only in the
     *  Tinker mass function; it is the background overdensity, not
     *  the critical overdensity &Delta;<SUB>crit</SUB>; the function
     *  cosmobl::Cosmology::DeltaR can be used to convert
     *  &Delta;<SUB>crit</SUB> into &Delta;
     */
    double bias_halo (const double, const double, const double, const string, const string output_root="test", const double Delta=200., const double kk=-1., const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string method_SS="CAMB", const string file_par=par::defaultString);
  
    /**
     *  @brief effective bias of dark matter haloes
     *
     *  @param Mass_min minimum halo mass
     *
     *  @param Mass_max maximum halo mass
     *
     *  @param redshift redshift
     *
     *  @param author_bias author(s) who proposed the bias; valid
     *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
     *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
     *  correction of Warren 2004), Tinker (Tinker et al. 2010)
     *
     *  @param author_MF author(s) who proposed the mass function;
     *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
     *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
     *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
     *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
     *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
     *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
     *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
     *  al. 2012)
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param Delta &Delta;: the overdensity, defined as the mean
     *  interior density relative to the background
     *  
     *  @param kk wave vector module
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return b<SUB>eff</SUB>: the effective dark matter bias
     *
     *  @warning the input parameter &Delta; is used only in the
     *  Tinker mass function; it is the background overdensity, not
     *  the critical overdensity &Delta;<SUB>crit</SUB>; the function
     *  cosmobl::Cosmology::DeltaR can be used to convert
     *  &Delta;<SUB>crit</SUB> into &Delta;
     */
    double bias_eff (const double, const double, const double, const string, const string, const string, const string output_root="test", const double Delta=200., const double kk=-1., const string interpType="Linear", const int Num=-1, const double stepsize=100., const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);
 
    /**
     *  @brief effective bias of dark matter haloes
     *
     *  @param MM vector of halo masses
     *
     *  @param MF vector of mass function values, d&Phi;/dM=dn(M)/dM
     *
     *  @param redshift redshift
     *
     *  @param author_bias author(s) who proposed the bias; valid authors
     *  are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen
     *  2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction
     *  of Warren 2004), Tinker (Tinker et al. 2010)
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param Delta &Delta;: the overdensity, defined as the mean
     *  interior density relative to the background
     *  
     *  @param kk wave vector module
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return b<SUB>eff</SUB>: the effective dark matter bias
     *
     *  @warning the input parameter &Delta; is used only in the
     *  Tinker mass function; it is the background overdensity, not
     *  the critical overdensity &Delta;<SUB>crit</SUB>; the function
     *  cosmobl::Cosmology::DeltaR can be used to convert
     *  &Delta;<SUB>crit</SUB> into &Delta;
     */
    double bias_eff (const vector<double>, const vector<double>, const double, const string, const string, const string output_root="test", const double Delta=200., const double kk=-1., const string interpType="Linear", const int Num=-1, const double stepsize=100., const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    ///@}


    /**
     *  @name Functions to model redshift-space distortions
     */
    ///@{

    /**
     *  @brief the linear growth rate
     *  @param redshift redshift
     *  @param kk wave vector module
     *  @return f: the linear growth rate
     */
    double linear_growth_rate (const double, const double kk=-1.) const;

    /**
     *  @brief f*&sigma;<SUB>8</SUB>: the linear growth rate times
     *  the dark matter rms mass fluctuation within 8 Mpc/h
     *
     *  @param redshift redshift
     *
     *  @param method_Pk method used to compute the power spectrum and
     *  &sigma;(mass); valid choices for method_Pk are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *  
     *  @param output_root output_root of the parameter file used to compute
     *  the power spectrum and &sigma;(mass); it can be any name
     *
     *  @param kk wave vector module 
     *
     *  @param NL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power
     *  spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *   
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr; the GSL
     *  libraries are used
     *
     *  @param prec accuracy of the GSL integration
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return f*&sigma;<SUB>8</SUB>
     */
    double fsigma8 (const double, const string, const string output_root="test", const double kk=-1., const bool NL=0, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString) const;

    /**
     *  @brief the specific growth rate &beta;
     *  @param redshift redshift
     *  @param bias bias
     *  @param kk wave vector module
     *  @return &beta;=f/b, where f is the linear growth rate and b is
     *  the bias
     */
    double beta (const double, const double, const double kk=-1.) const;

    /**
     *  @brief the error on the specific growth rate &beta;
     *  @param redshift redshift
     *  @param bias bias
     *  @param err_bias error on the bias
     *  @param kk wave vector module
     *  @return error on &beta;=f/b, where f is the linear growth rate
     *  and b is the bias
     */
    double error_beta (const double, const double, const double, const double kk=-1.) const;

    /**
     *  @brief the error on the specific growth rate &beta;
     *
     *  @param Mass_min minimum halo mass
     * 
     *  @param Mass_max maximum halo mass
     *
     *  @param redshift redshift
     *
     *  @param author_bias author(s) who proposed the bias; valid
     *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
     *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
     *  correction of Warren 2004), Tinker (Tinker et al. 2010)
     *
     *  @param author_MF author(s) who proposed the mass function;
     *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
     *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
     *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
     *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
     *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
     *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
     *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
     *  al. 2012)
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param Delta &Delta;: the overdensity, defined as the mean
     *  interior density relative to the background
     *  
     *  @param kk wave vector module
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return &beta;=f/b, where f is the linear growth rate and b is
     *  the bias
     *
     *  @warning the input parameter &Delta; is used only in the
     *  Tinker mass function; it is the background overdensity, not
     *  the critical overdensity &Delta;<SUB>crit</SUB>; the function
     *  cosmobl::Cosmology::DeltaR can be used to convert
     *  &Delta;<SUB>crit</SUB> into &Delta;
     */
    double beta (const double, const double, const double, const string, const string, const string, const string output_root="test", const double Delta=200., const double kk=-1., const string interpType="Linear", const int Num=-1, const double stepsize=100., const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    /**
     *  @brief the specific growth rate &beta;
     *
     *  @param Mass_min minimum halo mass
     * 
     *  @param Mass_max maximum halo mass
     *
     *  @param redshift redshift
     *
     *  @param author_bias author(s) who proposed the bias; valid
     *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
     *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
     *  correction of Warren 2004), Tinker (Tinker et al. 2010)
     *
     *  @param author_MF author(s) who proposed the mass function;
     *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
     *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
     *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
     *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
     *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
     *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
     *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
     *  al. 2012)
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param err_bias error on the bias
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param Delta &Delta;: the overdensity, defined as the mean
     *  interior density relative to the background
     *  
     *  @param kk wave vector module
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return error on &beta;=f/b, where f is the linear growth rate and b is
     *  the bias
     *
     *  @warning the input parameter &Delta; is used only in the
     *  Tinker mass function; it is the background overdensity, not
     *  the critical overdensity &Delta;<SUB>crit</SUB>; the function
     *  cosmobl::Cosmology::DeltaR can be used to convert
     *  &Delta;<SUB>crit</SUB> into &Delta;
     */
    double error_beta (const double, const double, const double, const string, const string, const string, const double, const string output_root="test", const double Delta=200., const double kk=-1., const string interpType="Linear", const int Num=-1, const double stepsize=100., const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString); 
  
    /**
     *  @brief the specific growth rate &beta;
     *
     *  @param MM vector of halo masses
     *
     *  @param MF vector of mass function values, d&Phi;/dM=dn(M)/dM
     *
     *  @param redshift redshift
     *
     *  @param author_bias author(s) who proposed the bias; valid authors
     *  are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen
     *  2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction
     *  of Warren 2004), Tinker (Tinker et al. 2010)
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param Delta &Delta;: the overdensity, defined as the mean
     *  interior density relative to the background
     *  
     *  @param kk wave vector module
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return &beta;=f/b, where f is the linear growth rate and b is
     *  the bias
     *
     *  @warning the input parameter &Delta; is used only in the
     *  Tinker mass function; it is the background overdensity, not
     *  the critical overdensity &Delta;<SUB>crit</SUB>; the function
     *  cosmobl::Cosmology::DeltaR can be used to convert
     *  &Delta;<SUB>crit</SUB> into &Delta;
     */
    double beta (const vector<double>, const vector<double>, const double, const string, const string, const string output_root="test", const double Delta=200., const double kk=-1., const string interpType="Linear", const int Num=-1, const double stepsize=100., const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    /**
     *  @brief the error on the specific growth rate &beta;
     *
     *  @param MM vector of halo masses
     *
     *  @param MF vector of mass function values, d&Phi;/dM=dn(M)/dM
     *
     *  @param redshift redshift
     *
     *  @param author_bias author(s) who proposed the bias; valid authors
     *  are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen
     *  2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction
     *  of Warren 2004), Tinker (Tinker et al. 2010)
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param err_bias error on the bias
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param Delta &Delta;: the overdensity, defined as the mean
     *  interior density relative to the background
     *  
     *  @param kk wave vector module
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return error on &beta;=f/b, where f is the linear growth rate and b is
     *  the bias
     *
     *  @warning the input parameter &Delta; is used only in the
     *  Tinker mass function; it is the background overdensity, not
     *  the critical overdensity &Delta;<SUB>crit</SUB>; the function
     *  cosmobl::Cosmology::DeltaR can be used to convert
     *  &Delta;<SUB>crit</SUB> into &Delta;
     */
    double error_beta (const vector<double>, const vector<double>, const double, const string, const string, const double, const string output_root="test", const double Delta=200., const double kk=-1., const string interpType="Linear", const int Num=-1, const double stepsize=100., const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);
 
    /**
     *  @brief the error on the specific growth rate &beta; from
     * Bianchi et al. 2012
     *
     *  @param Volume comoving volume 
     *
     *  @param density comoving density
     *
     *  @param Mass_min minimum halo mass
     * 
     *  @param Mass_max maximum halo mass
     *
     *  @param redshift redshift
     *
     *  @param author_bias author(s) who proposed the bias; valid
     *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
     *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
     *  correction of Warren 2004), Tinker (Tinker et al. 2010)
     *
     *  @param author_MF author(s) who proposed the mass function;
     *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
     *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
     *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
     *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
     *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
     *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
     *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
     *  al. 2012)
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param Delta &Delta;: the overdensity, defined as the mean
     *  interior density relative to the background
     *  
     *  @param kk wave vector module
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return error on &beta;=f/b, where f is the linear growth rate
     *  and b is the bias
     *
     *  @warning the input parameter &Delta; is used only in the
     *  Tinker mass function; it is the background overdensity, not
     *  the critical overdensity &Delta;<SUB>crit</SUB>; the function
     *  cosmobl::Cosmology::DeltaR can be used to convert
     *  &Delta;<SUB>crit</SUB> into &Delta;
     */
    double error_beta_measured (const double, const double, const double, const double, const double, const string, const string, const string, const string output_root="test", const double Delta=200., const double kk=-1., const string interpType="Linear", const int Num=-1, const double stepsize=100., const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString); 

    /**
     *  @brief the normalized quadrupole Q
     *
     *  @param Mass_min minimum halo mass
     * 
     *  @param Mass_max maximum halo mass
     *
     *  @param redshift redshift
     *
     *  @param author_bias author(s) who proposed the bias; valid
     *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
     *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
     *  correction of Warren 2004), Tinker (Tinker et al. 2010)
     *
     *  @param author_MF author(s) who proposed the mass function;
     *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
     *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
     *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
     *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
     *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
     *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
     *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
     *  al. 2012)
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param Delta &Delta;: the overdensity, defined as the mean
     *  interior density relative to the background
     *  
     *  @param kk wave vector module
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return Q: the normalized quadrupole
     *
     *  @warning the input parameter &Delta; is used only in the
     *  Tinker mass function; it is the background overdensity, not
     *  the critical overdensity &Delta;<SUB>crit</SUB>; the function
     *  cosmobl::Cosmology::DeltaR can be used to convert
     *  &Delta;<SUB>crit</SUB> into &Delta;
     */
    double quadrupole (const double, const double, const double, const string, const string, const string, const string output_root="test", const double Delta=200., const double kk=-1., const string interpType="Linear", const int Num=-1, const double stepsize=100., const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString); 

    /**
     *  @brief the normalized quadrupole Q
     *
     *  @param MM vector of halo masses
     *
     *  @param MF vector of mass function values, d&Phi;/dM=dn(M)/dM
     *
     *  @param redshift redshift
     *
     *  @param author_bias author(s) who proposed the bias; valid authors
     *  are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen
     *  2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction
     *  of Warren 2004), Tinker (Tinker et al. 2010)
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param Delta &Delta;: the overdensity, defined as the mean
     *  interior density relative to the background
     *  
     *  @param kk wave vector module
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return Q: the normalized quadrupole
     *
     *  @warning the input parameter &Delta; is used only in the
     *  Tinker mass function; it is the background overdensity, not
     *  the critical overdensity &Delta;<SUB>crit</SUB>; the function
     *  cosmobl::Cosmology::DeltaR can be used to convert
     *  &Delta;<SUB>crit</SUB> into &Delta;
     */
    double quadrupole (const vector<double>, const vector<double>, const double, const string, const string, const string output_root="test", const double Delta=200., const double kk=-1., const string interpType="Linear", const int Num=-1, const double stepsize=100., const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    /**
     *  @brief the mean square bulk flow
     *
     *  @param rr comoving radius 
     *
     *  @param k_int_min minimum wave vector module up to which the
     *  integral is computed
     *
     *  @param method_Pk method used to compute the power spectrum and
     *  &sigma;(mass); valid choices for method_Pk are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift the redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return the mean square bulk flow
     */
    double square_bulk_flow (const double, const double, const string, const double, const string output_root="test", const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    /**
     *  @brief the mean square bulk flow
     *
     *  @param rr comoving radius 
     *
     *  @param k_int_min minimum wave vector module up to which the
     *  integral is computed
     *
     *  @param lgkk vector of log(k)
     *
     *  @param lgPk vector of log(P(k))
     *
     *  @param redshift redshift
     *
     *  @return the mean square bulk flow
     */
    double square_bulk_flow_Table (const double, const double, const vector<double>, const vector<double>, const double) const; 

    /**
     *  @brief the mean square velocity dispersion
     *
     *  @param rr comoving radius 
     *
     *  @param k_int_min minimum wave vector module up to which the
     *  integral is computed
     *
     *  @param method_Pk method used to compute the power spectrum and
     *  &sigma;(mass); valid choices for method_Pk are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift the redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return the mean square velocity dispersion
     */
    double square_velocity_dispersion (const double, const double, const string, const double, const string output_root="test", const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);
    
    /**
     *  @brief the Cosmic Mach Number
     *
     *  @param rr comoving radius 
     *
     *  @param k_int_min minimum wave vector module up to which the
     *  integral is computed
     *
     *  @param method_Pk method used to compute the power spectrum and
     *  &sigma;(mass); valid choices for method_Pk are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum and &sigma;(mass); it can be any
     *  name
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return the Cosmic Mach Number
     */
    double CMN (const double, const double, const string, const double, const string output_root="test", const double k_max=100., const string file_par=par::defaultString) const;

    /**
     *  @brief the hierarchical moments S<SUB>n</SUB>
     *
     *  this function provides the hierarchical moments S<SUB>n</SUB>
     *  given by the perturbation theory (see e.g. Juszkiewicz et
     *  al. 1993, Bernardeau 1994, Wolk 2013)
     *
     *  @param nn order of the moment
     *
     *  @param RR comoving separation
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param output_root output_root of the parameter file used to compute
     *  the power spectrum and &sigma;(mass); it can be any name
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return the hierarchical moments, S<SUB>n</SUB>, given by the
     *  perturbation theory
     */
    double Sn_PT (const int, const double, const string, const string output_root="test", const string interpType="Linear", const int Num=-1, const double stepsize=100., const double k_max=100., const string file_par=par::defaultString) const;
  
    /**
     *  @brief the deprojected hierarchical moments &Sigma;<SUB>n</SUB>
     *
     *  this function provides the deprojected hierarchical moments
     *  &Sigma;<SUB>n</SUB> given by the perturbation theory (see
     *  e.g. Juszkiewicz et al. 1993, Bernardeau 1994, Wolk 2013)
     *
     *  @param nn order of the moment
     *
     *  @param RR comoving separation
     *
     *  @param method_SS method used to compute the power spectrum and
     *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param output_root output_root of the parameter file used to compute
     *  the power spectrum and &sigma;(mass); it can be any name
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return the deprojected hierarchical moments,
     *  &Sigma;<SUB>n</SUB>, given by the perturbation theory
     */
    double Sigman_PT (const int, const double, const string, const string output_root="test", const string interpType="Linear", const int Num=-1, const double stepsize=100, const double k_max=100., const string file_par=par::defaultString) const;
    
    /**
     *  @brief 1D monopole in the Kaiser limit
     *
     *  this function provides the monopole &xi;<SUB>0</SUB>(r) predicted
     *  at large scales, in the Kaiser limit
     *
     *  @param rad comoving separation
     *
     *  @param f_sigma8 f*&sigma;<SUB>8</SUB>
     *
     *  @param bias_sigma8 b*&sigma;<SUB>8</SUB>
     *
     *  @param method_Pk method used to compute the power spectrum and
     *  &sigma;(mass); valid choices for method_Pk are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to compute
     *  the power spectrum and &sigma;(mass); it can be any name
     *
     *  @param xiType 0 &rarr; standard; 1 &rarr; Chuang & Wang model
     *
     *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
     *
     *  @param xiNL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power
     *  spectrum
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param r_min minimum separation up to which the
     *  correlation function is computed
     *
     *  @param r_max maximum separation up to which the
     *  correlation function is computed
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param aa parameter \e a of Eq. 24 of Anderson et
     *  al. 2012
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return &xi;<SUB>0</SUB>
     *
     */
    double xi0_Kaiser (const double, const double, const double, const string, const double, const string output_root="test", const bool xiType=0, const double k_star=-1., const bool xiNL=0, const int norm=-1, const double r_min=0.1, const double r_max=150., const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    /**
     *  @brief 2D correlation function, &xi;(r<SUB>p</SUB>,&pi;),
     *  predicted by the dispersion model
     *
     *  @param rp r<SUB>p</SUB>: the comoving separation perpendicular
     *  to the line-of-sight
     *
     *  @param pi &pi;: the comoving separation parallel to the
     *  line-of-sight
     *
     *  @param f_sigma8 f*&sigma;<SUB>8</SUB>
     *
     *  @param bias_sigma8 b*&sigma;<SUB>8</SUB>
     *
     *  @param sigma12 &sigma;<SUB>12</SUB>: pairwise peculiar velocity
     *  dispersion
     *
     *  @param method_Pk method used to compute the power spectrum and
     *  &sigma;(mass); valid choices for method_Pk are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param redshift redshift
     *
     *  @param FV 0 &rarr; exponential form for f(v); 1 &rarr; Gaussian form
     *  for f(v); where f(v) is the velocity distribution function
     *
     *  @param NL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power
     *  spectrum
     *
     *  @param rr vector of r, the module of the comoving
     *  separation
     *
     *  @param Xi vector of &xi;(r), the two-point correlation
     *  function of dark matter
     *
     *  @param Xi_ vector of barred &xi;(r),
     *
     *  @param Xi__ vector of double-barred &xi;(r)
     *
     *  @param output_root output_root of the parameter file used to compute
     *  the power spectrum and &sigma;(mass); it can be any name
     *
     *  @param index internal parameter used when minimizing the
     *  &chi;<SUB>2</SUB>
     *
     *  @param bias_nl 0 &rarr; linear bias; 1 &rarr; non-linear bias 
     *
     *  @param bA b<SUB>a</SUB> non-linear bias parameter
     *
     *  @param xiType 0 &rarr; standard; 1 &rarr; Chuang & Wang model
     *
     *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
     *
     *  @param xiNL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power
     *  spectrum
     *
     *  @param v_min minimum velocity used in the convolution of the
     *  correlation function
     *
     *  @param v_max maximum velocity used in the convolution of the
     *  correlation function
     *
     *  @param step_v number of steps used in the convolution of the
     *  correlation function
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1 &rarr;
     *  normalize the power spectrum
     *
     *  @param r_min minimum separation up to which the
     *  correlation function is computed
     *
     *  @param r_max maximum separation up to which the
     *  correlation function is computed
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param aa parameter \e a of Eq. 24 of Anderson et
     *  al. 2012
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return &xi;(r<SUB>p</SUB>,&pi;)
     */
    double xi2D_DispersionModel (const double, const double, const double, const double, const double, const string, const double, const int, const bool, vector<double>, vector<double> &, vector<double> &, vector<double> &, const string output_root="test", const int index=-1, const bool bias_nl=0, const double bA=-1., const bool xiType=0, const double k_star=-1., const bool xiNL=0, const double v_min=-3000., const double v_max=3000., const int step_v=500, const int norm=-1, const double r_min=0.1, const double r_max=150., const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    /**
     *  @brief the function &xi;<SUB>*</SUB> of the Chuang & Wang 2012
     *  model
     *
     *  see Chuang & Wang 2012, 1209.0210
     *
     *  @param rr comoving separation
     *
     *  @param redshift redshift
     *
     *  @param output_root output_root of the parameter file used to compute
     *  the power spectrum and &sigma;(mass); it can be any name
     *
     *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return &xi;<SUB>*</SUB>
     */
    double xi_star (const double, const double, const string output_root="test", const double k_star=-1., const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);
  
    /**
     *  @brief the function &xi;<SUB>g,nw</SUB>(s) of the Chuang &
     *  Wang 2012 model
     *
     *  see Chuang & Wang 2012, 1209.0210
     *
     *  @param rp r<SUB>p</SUB>: the comoving separation perpendicular
     *  to the line-of-sight
     *
     *  @param pi &pi;: the comoving separation parallel to the
     *  line-of-sight
     *      
     *  @param f_sigma8 f*&sigma;<SUB>8</SUB>
     *
     *  @param bias_sigma8 b*&sigma;<SUB>8</SUB>
     *
     *  @param bA b<SUB>a</SUB> non-linear bias parameter
     *
     *  @param redshift redshift
     *    
     *  @param rr vector of r, the module of the comoving
     *  separation
     *
     *  @param Xi vector of &xi;(r), the two-point correlation
     *  function of dark matter
     *
     *  @param Xi_ vector of barred &xi;(r),
     *
     *  @param Xi__ vector of double-barred &xi;(r)
     *
     *  @param output_root output_root of the parameter file used to compute
     *  the power spectrum and &sigma;(mass); it can be any name
     *
     *  @return &xi;<SUB>g,nw</SUB>(s)
     */
    double xisnl_gnw (const double, const double, const double, const double, const double, const double, const vector<double>, const vector<double>, vector<double> &, vector<double> &, const string output_root="test");
 
    /**
     *  @brief the function &xi;<SUB>g,BAO</SUB>(s) of the Chuang &
     *  Wang 2012 model
     *
     *  see Chuang & Wang 2012, 1209.0210
     *
     *  @param rp r<SUB>p</SUB>: the comoving separation perpendicular
     *  to the line-of-sight
     *
     *  @param pi &pi;: the comoving separation parallel to the
     *  line-of-sight
     *      
     *  @param f_sigma8 f*&sigma;<SUB>8</SUB>
     *
     *  @param bias_sigma8 b*&sigma;<SUB>8</SUB>
     *
     *  @param redshift redshift
     *
     *  @param rr vector of r, the module of the comoving
     *  separation
     *
     *  @param Xi vector of &xi;(r), the two-point correlation
     *  function of dark matter
     *
     *  @param Xi_ vector of barred &xi;(r),
     *
     *  @param Xi__ vector of double-barred &xi;(r)
     *
     *  @param output_root output_root of the parameter file used to compute
     *  the power spectrum and &sigma;(mass); it can be any name
     *
     *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
     *
     *  @param x_min minimum velocity used in the integral of the
     *  Chuang & Wang model
     *
     *  @param x_max maximum velocity used in the integral of the
     *  Chuang & Wang model
     *
     *  @param step_x number of steps in the integral of the Chuang &
     *  Wang model
     *
     *  @return &xi;<SUB>g,BAO</SUB>(s)
     */
    double xis_gBAO (const double, const double, const double, const double, const double, vector<double>, vector<double>, vector<double> &, vector<double> &, const string output_root="test", const double k_star=-1., const double x_min=-3000., const double x_max=3000., const int step_x=500);
 
    /**
     *  @brief 2D correlation function, &xi;(r<SUB>p</SUB>,&pi;),
     *  predicted by the Chuang & Wang model
     *
     *  @param rp r<SUB>p</SUB>: the comoving separation perpendicular
     *  to the line-of-sight
     *
     *  @param pi &pi;: the comoving separation parallel to the
     *  line-of-sight
     *
     *  @param beta &beta;=f/b, where f is the linear growth rate
     *  and b is the bias
     *
     *  @param bias_lin linear bias
     *
     *  @param bA b<SUB>a</SUB> non-linear bias parameter
     *
     *  @param sigmav0 &sigma;<SUB>0</SUB>(v): parameter of the
     *  velocity distribution function, f(v)
     *
     *  @param cmu parameter of the velocity distribution function,
     *  f(v)
     *
     *  @param cs1 parameter of the velocity distribution function,
     *  f(v)
     *
     *  @param cs2 parameter of the velocity distribution function,
     *  f(v)
     *
     *  @param redshift redshift
     *
     *  @param rr1 vector of r, the module of the comoving separation
     *
     *  @param Xi1 vector of &xi;(r), the two-point correlation
     *  function of dark matter
     *
     *  @param rr2 vector of r, the module of the comoving separation
     *
     *  @param Xi2 vector of &xi;(r), the two-point correlation
     *  function of dark matter
     *
     *  @param Xi1_ vector of barred &xi;(r),
     *
     *  @param Xi1__ vector of double-barred &xi;(r)
     *   
     *  @param Xi2_ vector of barred &xi;(r),
     *
     *  @param Xi2__ vector of double-barred &xi;(r)
     *
     *  @param output_root output_root of the parameter file used to compute
     *  the power spectrum and &sigma;(mass); it can be any name
     *
     *  @param BAO 0 &rarr; no BAO convolution; 1 &rarr; BAO convolution
     *
     *  @param xiType 0 &rarr; standard; 1 &rarr; Chuang & Wang model
     *
     *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
     *
     *  @param xiNL 0 &rarr; linear power spectrum; 1 &rarr; non-linear power
     *  spectrum
     *
     *  @param r_min minimum separation up to which the
     *  correlation function is computed
     *
     *  @param r_max maximum separation up to which the
     *  correlation function is computed
     *
     *  @param v_min minimum velocity used in the convolution of the
     *  correlation function
     *
     *  @param v_max maximum velocity used in the convolution of the
     *  correlation function
     *
     *  @param step_v number of steps used in the convolution of the
     *  correlation function
     *
     *  @param k_min minimum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the
     *  power spectrum is computed
     *
     *  @param x_min minimum velocity used in the integral of the
     *  Chuang & Wang model
     *
     *  @param x_max maximum velocity used in the integral of the
     *  Chuang & Wang model
     *
     *  @param step_x number of steps in the integral of the Chuang &
     *  Wang model
     *
     *  @param aa parameter \e a of Eq. 24 of Anderson et
     *  al. 2012
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration
     *
     *  @param file_par name of the parameter file; if a
     *  parameter file is provided (i.e. file_par!=NULL), it will be
     *  used, ignoring the cosmological parameters of the object
     *
     *  @return &xi;(r<SUB>p</SUB>,&pi;)
     */
    double xi2D_CW (const double, const double, const double, const double, const double, const double, const double, const double, const double, const double, vector<double>, vector<double>, vector<double>, vector<double>, vector<double> &, vector<double> &, vector<double> &, vector<double> &, const string output_root="test", const bool BAO=1, const bool xiType=0, const double k_star=-1, const bool xiNL=0, const double r_min=0.1, const double r_max=150., const double v_min=-3000., const double v_max=3000., const int step_v=500, const double k_min=0., const double k_max=100., const double x_min=-3000., const double x_max=3000., const int step_x=500, const double aa=0., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    ///@}


    /**
     *  @name Functions to model baryon acoustic oscillations
     */
    ///@{
    
    /**
     *  @brief the sound horizon at the drag epoch
     *  r<SUB>s</SUB>(z<SUB>d</SUB>), valid choices for method_Pk are:
     *  EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html],
     *  CAMB [http://camb.info/]
     *
     *  @author Alfonso Veropalumbo
     *  @author alfonso.veropalumbo@unibo.it
     *
     *  @param method_Pk the method to compute the sound horizon
     *  
     *  @param T_CMB T<SUB>CMB</SUB>: the present day CMB temperature [K]
     *     
     *  @return r<SUB>s</SUB>
     */
    double rs (const string, const double T_CMB=par::TCMB) const;

    /**
     *  @brief the sound horizon at the drag epoch predicted by
     *  Eisentein & Hu 1998 
     *
     *  see Eisentein & Hu 1998, Section 2.1
     *
     *  @author Alfonso Veropalumbo 
     *  @author alfonso.veropalumbo@unibo.it 
     *  @param T_CMB CMB temperature
     *  @return r<SUB>s</SUB>
     */
    double rs_EH (const double T_CMB=par::TCMB) const;

    /**
     *  @brief the sound horizon at the drag epoch estimated with CAMB [http://camb.info/],
     *  analytical formula by Aubourg et al. 2014
     *
     *  see Anderson et al 2014, Eq. 16
     *  
     *  @author Alfonso Veropalumbo
     *  @author alfonso.veropalumbo@unibo.it
     *  @return r<SUB>s</SUB>
     */
    double rs_CAMB () const;
  
    /**
     *  @brief the fiducial cosmology independent ratio
     *  r<SUB>s</SUB>/D<SUB>V</SUB>,  valid choices for method_Pk are: EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html],
     *  CAMB [http://camb.info/]
     *
     *  both r<SUB>s</SUB> and D<SUB>V</SUB> are in Mpc
     *
     *  @author Alfonso Veropalumbo
     *  @author alfonso.veropalumbo@unibo.it
     *
     *  @param redshift redshift
     *  @param method_Pk method used to compute the sound horizon;
     *  @param T_CMB CMB temperature
     *
     *  @return y<SUB>s</SUB>
     */
    double ys (const double, const string, const double T_CMB=par::TCMB) const;
  
    /**
     *  @brief the acoustic parameter 
     *
     *  see Eisenstein 2005 
     *
     *  @author Alfonso Veropalumbo
     *  @author alfonso.veropalumbo@unibo.it
     *  @param redshift redshift
     *  @return the acoustic parameter 
     */
    double Az (const double) const;

    ///@}


    /**
     *  @name Functions to model cosmological quantities in non-Gaussian cosmologies
     */
    ///@{

    /**
     *  @brief the amplitude of the matter power spectrum
     *
     *  @author Cosimo Fedeli
     *  @author cosimo.fedeli@oabo.inaf.it
     *
     *  @param method_Pk method used to compute the power spectrum;
     *  valid choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     * 
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1
     *  &rarr; normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return A<SUB>m</SUB>
     */
    double Am (const string, const string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString); 

    /**
     *  @brief the potential spectral amplitude 
     *
     *  @author Cosimo Fedeli
     *  @author cosimo.fedeli@oabo.inaf.it
     *
     *  @param method_Pk method used to compute the power spectrum;
     *  valid choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     * 
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1
     *  &rarr; normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return the potential spectral amplitude
     */
    double potential_spectral_amplitude (const string, const string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    /**
     *  @brief the bispectrum
     *
     *  @author Cosimo Fedeli
     *  @author cosimo.fedeli@oabo.inaf.it
     *
     *  @param kk wave vector module
     *  
     *  @param method_Pk method used to compute the power spectrum; 
     *  valid choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     * 
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1
     *  &rarr; normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return the potential spectral amplitude
     */
    double bispectrum (const vector<double>, const string, const string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);
    
    /**
     *  @brief auxiliary function to estimate cosmological quantities
     *  in non-Gaussian cosmologies
     *
     *  @author Cosimo Fedeli
     *  @author cosimo.fedeli@oabo.inaf.it
     *
     *  @param kk wave vector module
     *
     *  @param mass halo mass
     *  
     *  @param method_Pk method used to compute the power spectrum; 
     *  valid choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     * 
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1
     *  &rarr; normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return mrk
     */
    double mrk (const double, const double, const string, const string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    /**
     *  @brief auxiliary function to estimate cosmological quantities
     *  in non-Gaussian cosmologies
     *
     *  @author Cosimo Fedeli
     *  @author cosimo.fedeli@oabo.inaf.it
     *
     *  @param kk wave vector module
     *
     *  @param mass halo mass
     *  
     *  @param method_Pk method used to compute the power spectrum;
     *  valid choices for method_Pk are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     * 
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1
     *  &rarr; normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return frk
     */
    double frk (const double, const double, const string, const string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    /// @cond TEST_NG
    double bias_kernel (const double, void *); 

    double frk_test (const double, const double, const string, const string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);
    /// @endcond


    /**
     *  @brief correction to the halo bias in non-Gaussian cosmologies
     *
     *  @author Cosimo Fedeli
     *  @author cosimo.fedeli@oabo.inaf.it
     *
     *  @param kk wave vector module
     *
     *  @param mass halo mass
     *  
     *  @param method_Pk method used to compute the power spectrum; valid
     *  choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     * 
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1
     *  &rarr; normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return bias correction
     */
    double bias_correction (const double, const double, const string, const string  output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    /**
     *  @brief the skewness
     *
     *  @author Cosimo Fedeli
     *  @author cosimo.fedeli@oabo.inaf.it
     *
     *  @param mass halo mass
     *  
     *  @param method_Pk method used to compute the power spectrum; valid
     *  choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     * 
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1
     *  &rarr; normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return skewness
     */
    double skewness (const double, const string, const string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    /**
     *  @brief the derivative of the skewness, ds/dM
     *
     *  @author Cosimo Fedeli
     *  @author cosimo.fedeli@oabo.inaf.it
     *
     *  @param mass halo mass
     *  
     *  @param method_Pk method used to compute the power spectrum; valid
     *  choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
     *  [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     * 
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1
     *  &rarr; normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return derivative of the skewness
     */
    double dskewnessdM (const double, const string, const string  output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    /**
     *  @brief correction to the halo mass in non-Gaussian cosmologies
     *
     *  @author Cosimo Fedeli
     *  @author cosimo.fedeli@oabo.inaf.it
     *
     *  @param mass the halo mass
     *  
     *  @param redshift the redshift
     *
     *  @param method_Pk method used to compute the power spectrum;
     *  valid choices for method_Pk are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     * 
     *  @param output_root output_root of the parameter file used to
     *  compute the power spectrum; it can be any name
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *   
     *  @param norm 0 &rarr; don't normalize the power spectrum; 1
     *  &rarr; normalize the power spectrum
     *
     *  @param k_min minimum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
     *  the GSL libraries are used
     *
     *  @param prec accuracy of the GSL integration 
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *
     *  @return bias correction
     */
    double MF_correction (const double, const double, const string, const string output_root="test", const string interpType="Linear", const int Num=-1, const double stepsize=100., const int norm=-1, const double k_min=0., const double k_max=100., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

    ///@}


    /**
     *  @name Functions to estimate the void size function
     */
    ///@{

    /**
     *  @brief Volume of the sphere of corresponding radius
     *
     *  @author Tommaso Ronconi
     *  @author tommaso.ronconi@studio.unibo.it
     *
     *  @param R the radius of the sphere
     *  
     *  @return volume of the sphere
     */
    double VolS (const double) const;
    
    /**
     *  @brief Linear (under)density contrast
     *
     *  @author Tommaso Ronconi
     *  @author tommaso.ronconi@studio.unibo.it
     *
     *  @param rho_vm the non linear density contrast: \f$\rho_v/\rho_m\f$ (default value set to \f$0.2\f$)
     *  
     *  @return The linear density contrast used as second barrier in the excursion set formalism
     *  for voids, as given by Bernardeu (1994):
     *  \f$\delta_v \equiv \frac{\rho_v - \rho_m}{\rho_m} \approx C [1 - (\rho_v/\rho_m)^{- 1/C}]\f$
     *  where \f$\rho_v =\ \f$ average void density, 
     *  \f$\rho_m =\ \f$ average density of the surrounding Universe and \f$C = 1.594\f$, a costant.
     */
    double deltav (const double rho_vm = 0.2) const;

    /**
     *  @brief expansion factor
     *
     *  @author Tommaso Ronconi
     *  @author tommaso.ronconi@studio.unibo.it
     *
     *  @param rho_vm the non linear density contrast: \f$\rho_v/\rho_m\f$ (default value set to \f$0.2\f$)
     *  
     *  @return the expansion factor: \f$\frac{r}{r_L} = \biggl(\frac{\rho_v}{\rho_m}\biggr)^{-1/3}\f$
     */
    double r_rL (const double rho_vm = 0.2) const;

    /**
     *  @brief \f$f_{\ln \sigma}(\sigma)\f$ (approximation)
     *
     *  @author Tommaso Ronconi
     *  @author tommaso.ronconi@studio.unibo.it
     *
     *  @param SS variance of the linear density field (\f$\sigma^2(R)\f$)
     *
     *  @param del_v linear density contrast defining a void
     *
     *  @param del_c critical value of the linear density field
     *  
     *  @return the fraction of trajectories that evolve into voids,
     *  as given in equation (8) of Jennings et al. (2013)
     */
    double f_nu (const double, const double, const double) const;

    /**
     *  @brief the void size function
     *
     *  @author Tommaso Ronconi
     *  @author tommaso.ronconi@studio.unibo.it
     *
     *  @param R radius
     *
     *  @param redshift redshift
     *
     *  @param rho_vm the non linear density contrast: \f$\rho_v/\rho_m\f$ (default value set to \f$0.2\f$)
     *
     *  @param del_v linear density contrast defining a void
     *
     *  @param del_c critical value of the linear density field
     *
     *  @param method_Pk method used to compute the power spectrum;
     *  valid choices for method_Pk are: CAMB [http://camb.info/],
     *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
     *  [http://arxiv.org/abs/1207.1465], EisensteinHu
     *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
     *
     *  @param output_root output_root of the parameter file used to compute
     *  the power spectrum and &sigma;(mass); it can be any name
     *
     *  @param interpType method to interpolate the power spectrum
     *
     *  @param Num number of near points used in the interpolation
     *
     *  @param stepsize width of the steps used in the derivative
     *  method
     *
     *  @param k_max maximum wave vector module up to which the power
     *  spectrum is computed
     *
     *  @param file_par name of the parameter file; if a parameter
     *  file is provided (i.e. file_par!=NULL), it will be used,
     *  ignoring the cosmological parameters of the object
     *  
     *  @return the number density of voids as a function of radius.
     *  Volume Conserving Model, equation (17) from Jennings et al.(2013) 
     */
    double size_function (const double, const double, const double, const double, const double, const string, const string, const string, const int, const double, const double, const string) const;

    ///@}

    /**
     *  @name Functions to estimate the multipoles/wedges covariance matrix
     */
    ///@{

    vector<vector<double>> get_XiMonopole_covariance(const int nbins, const double rMin, const double rMax, const double nn, const double Volume, const vector<double> kk, const vector<double> Pk0, const int IntegrationMethod=1);

    vector<vector<double>> get_XiMultipoles(const int nbins, const double rMin, const double rMax, const vector<double> kk, const vector<double> Pk0, const vector<double> Pk2, const vector<double> Pk4, const int IntegrationMethod=1);

    vector<vector<double>> get_XiMultipoles_covariance(const int nbins, const double rMin, const double rMax, const double nn, const double Volume, const vector<double> kk, const vector<double> Pk0, const vector<double> Pk2, const vector<double> Pk4, const int IntegrationMethod=1);

    ///@}
  };


  // =====================================================================================
  

  namespace glob {

    /// @cond glob
    double GSL_bias_kernel_wrapper (double, void *);
    
    double func_xi_EH_GSL (double, void *);
    double func_SSM_EH_GSL (double, void *);
    
    double bias_kernel2 (double, void *); // non-Gaussian cosmologies
    double skewness_kernel (double *, size_t, void *); // non-Gaussian cosmologies
    /// @endcond 
    
    struct GSL_f_pars {
      double kt;
      double mass;
      string method_Pk;
      string output_root;
      int norm;
      double k_min;
      double k_max;
      bool GSL;
      double prec;
      string file_par;
      Cosmology *pt_Cosmology;
    };

    struct STR_xi_EH
    {
      double Omega_matter;
      double Omega_baryon;
      double Omega_neutrinos;
      double massless_neutrinos;
      double massive_neutrinos;
      double Omega_DE;
      double Omega_radiation;
      double hh;
      double scalar_amp;
      double n_spec;
      double w0;
      double wa;
      double fNL;
      int type_NG;
      string model;
      bool unit;
      string method_Pk;
      double rr;
      double redshift;
      double aa;
    };

    struct STR_SSM_EH
    {
      double Omega_matter;
      double Omega_baryon;
      double Omega_neutrinos;
      double massless_neutrinos;
      double massive_neutrinos;
      double Omega_DE;
      double Omega_radiation;
      double hh;
      double scalar_amp;
      double n_spec;
      double w0;
      double wa;
      double fNL;
      int type_NG;
      string model;
      bool unit;
      string method_Pk;
      double redshift;
      double mass;
      double rho;
    };

    struct STR_NG /* Cosimo Fedeli */
    {
      double Omega_matter;
      double Omega_baryon;
      double Omega_neutrinos;
      double massless_neutrinos;
      double massive_neutrinos;
      double Omega_DE;
      double Omega_radiation;
      double hh;
      double scalar_amp;
      double n_spec;
      double w0;
      double wa;
      double fNL;
      int type_NG;
      string model;
      bool unit;
      double kt;
      double mass;
      string method_Pk;
      string output_root;
      int norm;
      double k_min;
      double k_max;
      bool GSL;
      double prec;
      string file_par;
    };
  }

}

#include "CosmClassFunc.h"

#endif
