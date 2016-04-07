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
 *  @file Headers/Lib/Constants.h
 *
 *  @brief Constants of general use
 *
 *  This file contains the global parameters and constants of the
 *  CosmoBolognaLib
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

namespace cosmobl {
  
  
  /**
   *  @brief The namespace of the global parameters and constants of
   *  the CosmoBolognaLib
   *
   *  The \e cosmoblpar namespace contains the global parameters and the
   *  mathematical, physical and astronomical constants of the
   *  CosmoBolognaLib
   */
  namespace par {

    /**
     *  @defgroup conv conversion factors
     *
     *  @brief The constant factors used to change units and to
     *  convert numbers into strings
     *
     *  @{
     */

    ///conversion factor 
    static const double yotta = 1.e24;
    
    ///conversion factor 
    static const double zetta = 1.e21;

    ///conversion factor 
    static const double exa = 1.e18;

    ///conversion factor 
    static const double peta = 1.e15;

    ///conversion factor 
    static const double tera = 1.e12;

    ///conversion factor 
    static const double giga = 1.e9;

    ///conversion factor 
    static const double mega = 1.e6;

    ///conversion factor 
    static const double kilo = 1.e3;

    ///conversion factor 
    static const double ecto = 1.e2;

    ///conversion factor 
    static const double deca = 10.;

    ///conversion factor 
    static const double deci = 1.e-1;

    ///conversion factor 
    static const double centi = 1.e-2;

    ///conversion factor 
    static const double milli = 1.e-3;

    ///conversion factor 
    static const double micro = 1.e-6;

    ///conversion factor 
    static const double nano = 1.e-9;

    ///conversion factor 
    static const double pico = 1.e-12;

    ///conversion factor 
    static const double femto = 1.e-15;

    ///conversion factor 
    static const double atto = 1.e-18;

    ///conversion factor 
    static const double zepto = 1.e-21;

    ///conversion factor 
    static const double yocto = 1.e-24;


    /// conversion symbol for: int -> string
    static const char fINT[] = "%i";

    /// conversion symbol for: long -> string
    static const char fLONG[] = "%lli";

    /// conversion symbol for: double -> string
    static const char fDP0[] = "%1.0f"; 

    /// conversion symbol for: double -> string
    static const char fDP1[] = "%2.1f"; 

    /// conversion symbol for: double -> string
    static const char fDP2[] = "%3.2f"; 

    /// conversion symbol for: double -> string
    static const char fDP3[] = "%4.3f"; 

    /// conversion symbol for: double -> string
    static const char fDP4[] = "%5.4f"; 
  
    /// conversion symbol for: double -> string
    static const char fDP5[] = "%6.5f"; 
  
    /// conversion symbol for: double -> string
    static const char fDP6[] = "%7.6f"; 

    /// conversion symbol for: double -> string
    static const char ee3[] = "%4.3e";

    /**
     *  @} */


    // ============================================================================================


    /**
     *  @defgroup math mathematical constants 
     *
     *  @brief Mathematical constants of wide use
     *
     *  @{
     */

    /// \f$\pi\f$: the ratio of a circle's circumference to its diameter 
    static const double pi = 3.1415926535897932;

    /// e: the Euler's constant, defined as \f$\lim_{n \to \infty}\left(1+\frac{1}{n}\right)^n\f$
    static const double ee = 2.7182818284590452;

    /**
     *  @} */
  

    // ============================================================================================


    /**
     *  @defgroup phys physical and astronomical constants 
     *  @brief Physical and astronomical constants of wide use
     *  @{
     */

    /// \f$\hbar\f$: the reduced Planck constant, the quantum of action [J s]
    static const double hbar = 1.054571726e-34;

    /// \f${\it c}\f$: the speed of light in vacuum (the value is exact) [km/sec] 
    static const double cc = 299792.458;      

    /// \f${\it k_B}\f$: the Boltzmann constant, the conversion factor relating temperature and energy [eV K<SUP>-1</SUP>]
    static const double kB = 1.3806488e23;
    
    /// \f$\sigma_{SB}=\frac{\pi^2k_B^4}{60\hbar^3c^2}\f$: the Stefan-Boltzmann constant, the factor relating the emissive power of a black body to the fourth power of its temperature [W m<SUP>-2</SUP> K<SUP>-4</SUP>]
    static const double sSB = 5.670373e-8;

    /// \f${e}\f$: the electrical charge of the electron [C]
    static const double el = 1.602176565e-19;
    
    /// \f${\alpha}=\frac{e^2}{\hbar c}\f$: the fine-structure constant
    static const double alpha = 7.2973525698e-3;

    /// \f${\epsilon_0}\f$: the electric constant (the value is exact) [F/m]
    static const double epsilon0 = 8.854187817e-12;

    /// \f${\mu_0}\f$: the magnetic constant (the value is exact) [N A<SUP>-2</SUP>]
    static const double mu0 = 12.566370614e-7; 

    /// \f${N_A}\f$: the Avogadro constant [mol<SUP>-1</SUP>]
    static const double NAv = 6.02214129e23; 

    /// \f${G_N}\f$: the Newtonian constant of gravitation [m<SUP>3</SUP> Kg<SUP>-1</SUP> s<SUP>-2</SUP>]
    static const double GN = 6.6738480e-11;

    /// \f${g_n}\f$: the standard gravitational acceleration (the value is exact) [m s<SUP>-2</SUP>]
    static const double gn = 9.80665;

    /// \f${l_P}=\sqrt{\hbar G_N/c^3}\f$: the Planck length [m]
    static const double lP = 1.616199e-35; 

    /// \f${M_P}=\sqrt{\hbar c/G_N}\f$: the Planck mass [kg]
    static const double MP = 2.17651e-8;

    /// \f${M_\odot}\f$: the solar mass [kg]
    static const double Msol = 1.98855e30;

    /// \f${m_e}\f$: the mass of the electron [kg]
    static const double me = 9.10938291e-31; 

    /// \f${m_n}\f$: the mass of the neutron [kg]
    static const double mn = 1.674927351e-27; 

    /// \f${m_p}\f$: the mass of the proton [kg]
    static const double mp = 1.672621777e-27;

    /// \f${au}\f$: the astronomical unit (the value is exact) [m]
    static const double au = 149597870700; 

    /// \f${pc}\f$: the parsec, defined as 1au/1arc sec [m] 
    static const double pc = 3.0856775814671916e16;

    /// \f${T_{CMB}}\f$: the present day CMB temperature [K]
    static const double TCMB = 2.72548;

    /// \f${yr}\f$: one year [s]
    static const double yr = 31557600;        
  
    /**
     *  @} */
    

    // ============================================================================================

    
    /**
     *  @defgroup col colours
     *
     *  @brief Some colour string factors, used when printing
     *  something on the screen
     *
     *  @{
     */

    /// default colour (used when printing something on the screen)
    static const string col_default = "\033[0m";

    /// red colour (used when printing something on the screen)
    static const string col_red = "\033[0;31m";

    /// green colour (used when printing something on the screen)
    static const string col_green = "\033[0;32m";
    
    /// blue colour (used when printing something on the screen)
    static const string col_blue = "\033[0;34m";

    /**
     *  @} */

    
    // ============================================================================================

  
    /**
     *  @defgroup default default values
     *
     *  @brief Default values of different types
     *
     *  @{
     */

    /// default string value
    static const string defaultString = "NULL";
    
    /// default integer value
    static const int defaultInt = numeric_limits<int>::min();

    /// default long value
    static const long defaultLong = numeric_limits<long>::min();

    /// default float value
    static const float defaultFloat = -numeric_limits<float>::max();
    
    /// default double value
    static const double defaultDouble = -numeric_limits<double>::max(); 
    
    /**
     *  @} */

    
    // ============================================================================================

  
    /**
     *  @defgroup dir useful directories
     *
     *  @brief Useful directories of general use
     *
     *  @{
     */
  
    /// directory where the CosmoBolognaLib are stored
    extern string DirCosmo;
    
    /// local directory of the main code
    extern string DirLoc;
    
    /**
     *  @} */
  
  }

}

