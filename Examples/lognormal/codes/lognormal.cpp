// ===================================================================
// Example code: how to construct a log-normal density field catalogue
// ===================================================================

#include "LogNormal.h"

int main () {

  try {

    // --------------------------------------
    // ----- set the cosmological model -----
    // --------------------------------------
    
    const cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck18_};


    // -------------------------------------------------------------
    // ----- set the object to construct log-normal catalogues -----
    // -------------------------------------------------------------
    
    // the random catalogues used to construct the mask
    const cbl::catalogue::Catalogue random {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_comoving_, {"../input/random.dat"}};

    // mean number of objects in the log-normal catalogues
    const int nObjects = 1110;

    // mean redshit of the log-normal catalogues
    const double redshift = 1.;
    
    // the bias of the log-normal density field catalogues
    const double bias = 1.;
    
    // the cell size in comoving scale [Mpc/h]
    const double cell_size = 5.;
    
    cbl::lognormal::LogNormal logNormal(random, cosmology, nObjects, redshift, bias, cell_size);


    // ---------------------------------------------------------------------------------------------
    // ----- construct one log-normal calalogue and get the number of objects in the catalogue -----
    // ---------------------------------------------------------------------------------------------

    if (system("mkdir -p ../output/")) {}
    
    logNormal.generate(1, "../output/");

    std::cout << "The number of objects in the log-normal catalogue is " << logNormal.catalogue(0)->nObjects() << std::endl;
    
  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }

  return 0;
}
