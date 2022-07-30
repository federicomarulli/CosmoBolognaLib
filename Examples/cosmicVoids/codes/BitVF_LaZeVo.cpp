// ==================================================================================================================
// Example code: how to to clean a cosmic void catalogue, to extract cosmological constraints from void counting file
// ==================================================================================================================

#include "Catalogue.h"

using namespace std;
using namespace cbl;
using namespace catalogue;

int main () {

  try {

    // ---------------------------------------
    // ----- build the tracer catalogue ------
    // ---------------------------------------

    // define the path of the tracers catalogue
    const std::string file_tracers = "../input/halo_catalogue.txt";
    
    const double scaleFact = 1.;
    const double nSub = 1.1;

    // std::vector containing the variable name list to read from file
    std::vector<cbl::catalogue::Var> var_names_tracers = {cbl::catalogue::Var::_X_, cbl::catalogue::Var::_Y_, cbl::catalogue::Var::_Z_};
      
    // std::vector containing the column corresponding to each attribute
    std::vector<int> columns_tracers = {1, 2, 3};

    // tracers catalogue constructor
    cbl::catalogue::Catalogue tracer_catalogue {cbl::catalogue::ObjectType::_Halo_, cbl::CoordinateType::_comoving_, var_names_tracers, columns_tracers, {file_tracers}, 0, nSub, scaleFact};
    
    // random catalogue constructor (empty if it is not given)
    Catalogue random_catalogue = {};

    // store the mean particle separation of the simulation
    double mps = tracer_catalogue.mps();

    // output directory 
    const string dir_output = "../output/";

    // output suffix 
    const string output = "output.dat";

    // chainmesh cell size
    const double cellsize = 4*mps;

    // number of reconstruction of the displacement field
    const int n_rec = 5;

    // cell size for the estimation of the divergence field (mps units)
    const double step_size = 2.5/3;

    // threshold at which to stop the reconstruction of the displacement field
    // the value must be between 0 and 1 (th->0 = more accuracy)
    const double threshold = 0.;

    // conditions for print the displacement and the divergence field 
    vector<bool> print={true,true}

    // catalogue constructor
    Catalogue void_catalogue = Catalogue(VoidAlgorithm::_LaZeVo_, tracer_catalogue, random_catalogue, dir_output, output, cellsize, n_rec, step_size, threshold, print);

  }
  
  catch (cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
   
  return 0;
}