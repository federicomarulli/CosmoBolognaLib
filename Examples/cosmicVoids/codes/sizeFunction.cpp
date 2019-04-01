// ==========================================================================
// Example code: how to compute the theoretical size function of cosmic voids
// ==========================================================================

#include "Cosmology.h"
 
// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
std::string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;
 
 
int main () {

  try {
   
    // -----------------------------------------------
    // ----- use default cosmological parameters ----- 
    // -----------------------------------------------
    
    cbl::cosmology::Cosmology cosm;


    // -----------------------------------------------
    // ------------ logarithmic binning -------------- 
    // -----------------------------------------------
    
    std::vector<double> RR; 
    int n_val = 10;   // number of values at which the size function is computed
    double R_min = 1.,  R_max = 30.;   // minimum and maximum values of effective void radii
    
    double min_log = std::log(R_min);
    double max_log = std::log(R_max);
    
    double log_increment = (max_log-min_log)/(n_val-1);
    double log_value = min_log;

    for(int i=0; i<n_val; i++){
      RR.push_back(std::exp(log_value));
      log_value += log_increment;
    }

    // -------------------------------------------------------------
    // ----- compute the Sheth & van de Weygaert size function -----
    // -------------------------------------------------------------

    double zz  = 0.;  // redshift of the sample
    double b_eff = 1.;   // effective bias of the mass tracers

    std::cout << "The size function at z = " << zz << " is :" << std::endl;

    for(int i=0; i<n_val; i++)
      std::cout << std::scientific << std::setprecision(6) << cosm.size_function(RR[i], zz, "SvdW", b_eff) << " (h/Mpc)^3 at R = " << std::defaultfloat << std::setprecision(3) << RR[i] << " Mpc/h " << std::endl;
  }
  
  catch (cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
   
  return 0;
}
