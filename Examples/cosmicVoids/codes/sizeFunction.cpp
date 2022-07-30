// ==========================================================================
// Example code: how to compute the theoretical size function of cosmic voids
// ==========================================================================

#include "Cosmology.h"
 
int main () {

  try {
     
    // -------------------------------------------------------------
    // ----- compute the Sheth & van de Weygaert size function -----
    // -------------------------------------------------------------
 
    // set the cosmology using Planck18 parameters 
    const cbl::cosmology::Cosmology Planck18 {cbl::cosmology::CosmologicalModel::_Planck18_};
    
    // logarithmic binning in the void effectiv radii
    const std::vector<double> radius = cbl::logarithmic_bin_vector(10, 1., 30.); 

    // redshift
    const double redshift  = 0.;

    // effective bias of the tracers (e.g. galaxies or galaxy clusters) used to detect the voids 
    const double b_eff = 1.; 


    // print the size function
    
    std::cout << "The size function at z = " << redshift << " is :" << std::endl;
    
    for (auto && rr : radius)
      std::cout << std::scientific << std::setprecision(6) << Planck18.size_function(rr, redshift, "SvdW", b_eff) << " (h/Mpc)^3 at R = " << std::defaultfloat << std::setprecision(3) << rr << " Mpc/h" << std::endl;

  }
  
  catch (cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
   
  return 0;
}
