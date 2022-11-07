//========================================================================================
// Author: Jens Chluba
// Last modification: Oct 2010
// CITA, University of Toronto
// All rights reserved.
//========================================================================================

#ifndef COSMOLOGY_RECFAST_H
#define COSMOLOGY_RECFAST_H

//========================================================================================
struct Input 
{
    double YP;      
    double fHe;
    double To;      // Temperature of CMB at z=0
    double OmegaM;    // Omega matter 
    double OmegaB;    // Omega Baryons 
    double OmegaK;    // Omega Curvature
    double OmegaL;    // Omega Lambda 
    double zstart;
    double zend;
    double F;           // fudge-factor; normally F=1.14
    double H0;
    double h100;
    double Nnu;     // effective number of neutrinos
    //
    double fDM;
    int switch_on_recombination_corrs;
};

//========================================================================================
// Hubble-function in 1/sec
//========================================================================================
double H_z(double z);

//========================================================================================
// hydrogen number density in m^-3
//========================================================================================
double NH(double z);

//========================================================================================
// CMB temperature at z
//========================================================================================
double TCMB(double z);

//========================================================================================
// compute total contribution from relativistic particles (photon & neutrinos)
//========================================================================================
double calc_Orel(double TCMB0, double Nnu, double h100);

#endif 
