// ==================================================================================================================
// Example code: how to to clean a cosmic void catalogue, to extract cosmological constraints from void counting file
// ==================================================================================================================

#include "Catalogue.h"

using namespace cbl;
using namespace catalogue;


int main () {

  try {

    // ------------------------------------------
    // ----- load the input void catalogue ------
    // ------------------------------------------

    // ASCII void catalogue 
    std::string file_voids = "../input/void_catalogue.txt";

    // std::vector containing the variable name list to read from file
    std::vector<cbl::catalogue::Var> var_names_voids = {cbl::catalogue::Var::_X_, cbl::catalogue::Var::_Y_, cbl::catalogue::Var::_Z_, cbl::catalogue::Var::_Radius_};
    
    // std::vector containing the columns corresponding to each attribute
    std::vector<int> columns_voids = {1, 2, 3, 4};
    
    // catalogue constructor
    cbl::catalogue::Catalogue void_catalogue = Catalogue(cbl::catalogue::ObjectType::_Void_, cbl::CoordinateType::_comoving_, var_names_voids, columns_voids, {file_voids}, 0);
     
    // ---------------------------------------
    // ----- build the tracer catalogue ------
    // ---------------------------------------
      
    // binary halo gadget snapshot 
    std::string file_tracers = "../input/halo_catalogue.txt";

    // std::vector containing the variable name list to read from file
    std::vector<cbl::catalogue::Var> var_names_tracers = {cbl::catalogue::Var::_X_, cbl::catalogue::Var::_Y_, cbl::catalogue::Var::_Z_};
      
    // std::vector containing the column corresponding to each attribute
    std::vector<int> columns_tracers = {1, 2, 3};

    // catalogue constructor
    cbl::catalogue::Catalogue tracers_catalogue = Catalogue(cbl::catalogue::ObjectType::_Halo_, cbl::CoordinateType::_comoving_, var_names_tracers, columns_tracers, {file_tracers}, 0);

    // store the mean particle separation of the simulation
    double mps = tracers_catalogue.mps();      

    // generate the chain mesh of the inpute tracer catalogue
    cbl::chainmesh::ChainMesh3D ChM(2*mps, tracers_catalogue.var(cbl::catalogue::Var::_X_), tracers_catalogue.var(cbl::catalogue::Var::_Y_), tracers_catalogue.var(cbl::catalogue::Var::_Z_), void_catalogue.Max(cbl::catalogue::Var::_Radius_));

    // make a shared pointer to tracers_catalogue
    auto input_tracersCata = std::make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(std::move(tracers_catalogue)));

    
    // --------------------------------------------
    // ----- build the cleaned void catalogue -----
    // --------------------------------------------

     double ratio = 0.8; // variable used to compute the central density and the density contrast of a void

    // ---- use built in catalogue constructor, selecting which step of the cleaning procedure to perform ----

    // clean[0] = true -> erase voids with underdensities higher than a given threshold
    // clean[1] = true -> erase voids with effective radii outside a given range
    // clean[2] = true -> erase voids with density contrast lower than a given value
    std::vector<bool> clean = {false, false, false};

    std::vector<double> delta_r = {17., 150.}; // the interval of accepted radii
    double threshold = 0.3;                   // the density threshold
    double relevance = 1.57;                  // the minimum accepted density contrast

    // catalogue constructor
    void_catalogue.clean_void_catalogue(clean, delta_r, threshold, relevance, true, input_tracersCata, ChM, ratio, true, cbl::catalogue::Var::_CentralDensity_);

    // store the obtained catalogue in an ASCII file
    var_names_voids.emplace_back(cbl::catalogue::Var::_CentralDensity_);
    std::string mkdir = "mkdir -p ../output/";
    if (system(mkdir.c_str())) {}

    std::string cata_out = "../output/cleaned_void_catalogue.out";
    void_catalogue.write_data(cata_out, var_names_voids);
    
  }
  
  catch (cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
   
  return 0;
}

