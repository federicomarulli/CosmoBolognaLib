// ========================================================================================================================
// Example code: how to model the 2D two-point correlation function in Cartesian coordinates, xi(rp, pi), in redshift-space
// ========================================================================================================================

#include "TwoPointCorrelation2D.h"

int main () {

  try {
  
    // -------------------------------------------------------------
    // ---------------- set the cosmological parameters ------------
    // -------------------------------------------------------------
  
    const cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck15_};
  
  
    // ----------------------------------------------------------
    // ---------------- read the input catalogue ----------------
    // ----------------------------------------------------------
  
    const std::string file_catalogue = "../input/cat.dat";

    const cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmology};
  
  
    // --------------------------------------------------------------------------------------
    // ---------------- construct the random catalogue (with cubic geometry) ---------------- 
    // -------------------------------------------------------------------------------------- 

    const double N_R = 3.; // random/data ratio
 
    const cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::RandomType::_createRandom_box_, catalogue, N_R};
  
  
    // -----------------------------------------------------------------------------------------------
    // ---------------- measure the 2D two-point correlation in Cartesian coordinates ----------------
    // -----------------------------------------------------------------------------------------------

    // ----- output data ----- 

    const std::string dir_pairs = "../output/";
    const std::string dir_output = "../output/";

  
    // ----- measure the 2D correlation function in Cartesian coordinates, xi(rp,pi), and store the outputs ----- 
  
    const double rpMin = 5.;     // minimum separation in the first dimension
    const double rpMax = 50.;    // maximum separation in the first dimension 
    const int nbins_D1 = 10;     // number of bins in the first dimension
    const double shift_D1 = 0.5; // spatial shift used to set the bin centre in the first dimension
    const double piMin = 5.;     // minimum separation in the second dimension
    const double piMax = 50.;    // maximum separation in the second dimension 
    const int nbins_D2 = 10;     // number of bins in the second dimension
    const double shift_D2 = 0.5; // spatial shift used to set the bin centre in the second dimension
    
    // construct the object using a static factory
    const auto xi2DCart = cbl::measure::twopt::TwoPointCorrelation::Create(cbl::measure::twopt::TwoPType::_2D_Cartesian_, catalogue, random_catalogue, cbl::BinType::_linear_, rpMin, rpMax, nbins_D1, shift_D1, cbl::BinType::_linear_, piMin, piMax, nbins_D2, shift_D2);
    
    // measure the 2D correlation function and compute Poisson errors
    xi2DCart->measure(cbl::measure::ErrorType::_Poisson_, dir_pairs);
    
    // store the output
    xi2DCart->write(dir_output, "xi_rp_pi_linlin.dat");

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

