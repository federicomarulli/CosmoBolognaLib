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

    // define the path of the tracers catalogue (NOT PROVIDED)
    const std::string file_tracers = "../input/your_tracer_catalogue.txt";
    
    const double scaleFact = 1.;
    const double nSub = 1.1;

    // std::vector containing the variable name list to read from file (IDs ARE NEEDED)
    std::vector<cbl::catalogue::Var> var_names_tracers = {cbl::catalogue::Var::_X_, cbl::catalogue::Var::_Y_, cbl::catalogue::Var::_Z_, cbl::catalogue::Var::_ID_};
      
    // std::vector containing the column corresponding to each attribute
    std::vector<int> columns_tracers = {1, 2, 3, 4};

    // tracers catalogue constructor
    cbl::catalogue::Catalogue tracer_catalogue {cbl::catalogue::ObjectType::_Halo_, cbl::CoordinateType::_comoving_, var_names_tracers, columns_tracers, {file_tracers}, 0, nSub, scaleFact};
    
    // define the path of the random catalogue (NOT PROVIDED)
    const std::string file_randoms = "../input/your_random_catalogue.txt";

    // random catalogue constructor (must provides tracer positions at high redshift, for example through the first snapshot of a simulation)
    cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::ObjectType::_Halo_, cbl::CoordinateType::_comoving_, var_names_tracers, columns_tracers, {file_randoms}, 0, nSub, scaleFact};

    // store the mean particle separation of the simulation
    double mps = tracer_catalogue.mps();

    // output directory 
    const string dir_output = "../output/";

    // output suffix 
    const string output = "output.dat";

    // chainmesh cell size
    const double cellsize = 4*mps;

    // number of reconstruction of the displacement field (IN THIS CASE, MUST BE = 1)
    const int n_rec = 1;

    // cell size for the estimation of the divergence field (mps units)
    const double step_size = 2.5/3;

    // threshold at which to stop the reconstruction of the displacement field
    // the value must be between 0 and 1 (NOT USED)
    const double threshold = 0.;

    // conditions for print the displacement and the divergence field 
    vector<bool> print={true,true}


    // > > > > >	WARNING! CATALOGUES MUST HAVE THE SAME OBJECTS AND BE ORDERED FOR THE IDs < < < < < 

    // v v v v v UNCOMMENT IF CATALOGUES ARE NOT ORDERED FOR THE ID OF THE OBJECTS v v v v v 
    /*
    vector<int> criteriumOrder = tracer_catalogue.var_int(Var::_ID_);
    vector<int> indices(tracer_catalogue.nObjects());
    iota(indices.begin(), indices.end(), 0);
    sort(indices.begin(), indices.end(), [&](int A, int B) -> bool {return criteriumOrder[A] < criteriumOrder[B];});
    tracer_catalogue.Order(indices);

    vector<int> criteriumOrderR = random_catalogue.var_int(Var::_ID_);
    vector<int> indicesR(random_catalogue.nObjects());
    iota(indicesR.begin(), indicesR.end(), 0);
    sort(indicesR.begin(), indicesR.end(), [&](int A, int B) -> bool {return criteriumOrderR[A] < criteriumOrderR[B];});
    random_catalogue.Order(indicesR);

    cout << endl;
    coutCBL << "Catalogues sorted" << endl << endl;
    */

    // v v v v v UNCOMMENT IF CATALOGUES HAVE DIFFERENT OBJECTS v v v v v 
    /*
    int i=0, j=0;
    while (i<(int)tracer_catalogue.nObjects() && j<(int)random_catalogue.nObjects()) {
      if (tracer_catalogue.ID(i) == random_catalogue.ID(j)) {
        i++;
        j++;
      }
      else if (tracer_catalogue.ID(i) > random_catalogue.ID(j)) random_catalogue.remove_object(j);
      else if (tracer_catalogue.ID(i) < random_catalogue.ID(j)) tracer_catalogue.remove_object(i);
    }
    
    if ((int)tracer_catalogue.nObjects()>i) for (size_t ii=i; ii<tracer_catalogue.nObjects(); ii++) tracer_catalogue.remove_object(ii);
    if ((int)random_catalogue.nObjects()>j) for (size_t jj=j; jj<random_catalogue.nObjects(); jj++) random_catalogue.remove_object(jj);
    
    cout << endl;
    coutCBL << "Unpair removed" << endl << endl;
    */

    // catalogue constructor
    Catalogue void_catalogue = Catalogue(VoidAlgorithm::_Exact_, tracer_catalogue, random_catalogue, dir_output, output, cellsize, n_rec, step_size, threshold, print);

  }
  
  catch (cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
   
  return 0;
}
