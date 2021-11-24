// =============================================================================================
// Example code: how to use HOD catalogue constructors and display catalogue summary statistics
// =============================================================================================
#include "Catalogue.h"
#include <iostream>
#include "Cosmology.h"
#include <string>
#include "GSLwrapper.h"
#include <gsl/gsl_math.h>
#include <omp.h>
#include "Distribution.h"
#include <random>
#include <chrono>

using namespace std;
using namespace cbl;

string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;
 
int main () {


    try {
      //set input and output paths:
      std::string input_path = "../input/";
      std::string output_path = "../output/";

      // Setting the cosmology used to read the halo catalogue
      cbl::cosmology::Cosmology cosm{};

      // String with input halo catalogue name:
      string input_file = "haloCat.dat";

      // vector with variables to use for the catalogue
      vector<cbl::catalogue::Var> attribute = {
        cbl::catalogue::Var::_X_,
        cbl::catalogue::Var::_Y_,
        cbl::catalogue::Var::_Z_,
        cbl::catalogue::Var::_Mass_
      };

        // vector with input file columns to be used according to variables vector
        vector<int> column = {1,2,3,4};

        // Number of row/s to ignore form the start of the file:
        const int comments = 0;

        // Fraction of objects that will be randomly selected 1.0 = 100%
        const double nSub = 1.0;

        // Factor for coordinates:
        const double fact = 1;

        // Catalogue Constructor:
        cbl::catalogue::Catalogue haloCat {
          cbl::catalogue::ObjectType::_HostHalo_, 
          cbl::CoordinateType::_comoving_, 
          attribute, column,
          {input_path + input_file},
          comments,
          nSub, 
          fact,
          cosm
        };

//--------------Creating a HOD catalogue from an HALO one based on Moster10 model:------------------

    // Setting the minimum stellar for the galaxy sample
    double threshold = pow(10., 8);
    
    //Setting the parameters vector. Empty vector takes default values
    vector<double> parameters = {};
    cbl::catalogue::Catalogue HOD {haloCat, cosm, cbl::catalogue::HODType::_Moster10_, threshold, false, parameters};


//-----------------------------Check MAX value of variables:----------------------------------------

    coutCBL << par::col_green << "HOD CATALOGUE MIN/MAX VALUES PER VARIABLE:" << par::col_default << endl;
    cout << "Min X: " << HOD.Min(cbl::catalogue::Var::_X_) << " - " << "Max X: " << HOD.Max(cbl::catalogue::Var::_X_) << endl;
    cout << "Min Y: " << HOD.Min(cbl::catalogue::Var::_Y_) << " - " << "Max Y: " << HOD.Max(cbl::catalogue::Var::_Y_) << endl;
    cout << "Min Z: " << HOD.Min(cbl::catalogue::Var::_Z_) << " - " << "Max Z: " << HOD.Max(cbl::catalogue::Var::_Z_) << endl;
    cout << "Min Mass: " << HOD.Min(cbl::catalogue::Var::_Mass_) << " - " << "Max Mass: " << HOD.Max(cbl::catalogue::Var::_Mass_) << endl;
    cout << "Total number of Galaxies: " << HOD.nObjects() << " - " << "Total number of Halos: " << haloCat.nObjects() << endl;

    // variable stats:
    vector<double> stats = {};
    HOD.stats_var(cbl::catalogue::Var::_Mass_, stats);

    cout << endl;
    coutCBL << par::col_green << "HOD CATALOGUE MASS STATS [Mpc/h]:" << par::col_default << endl;
    cout << "The mean: " << stats[0] << endl;
    cout << "The median: " << stats[1] << endl;
    cout << "The standard deviation: " << stats[2] << endl;
    cout << "Difference between the third and first quartiles: " << stats[3] << endl;
    cout << endl;

// --------------------------------Write HOD catalogue to file---------------------------------------

    // Vector with HOD catalogue parameters:
    std::vector<cbl::catalogue::Var> varsToPrint = {
      cbl::catalogue::Var::_GalaxyTag_,
      cbl::catalogue::Var::_X_,
      cbl::catalogue::Var::_Y_,
      cbl::catalogue::Var::_Z_,
      cbl::catalogue::Var::_Mass_,
      cbl::catalogue::Var::_MassInfall_,
      cbl::catalogue::Var::_Mstar_
    };

    // Setting file separator (optional, dafault is: "   " for compatibility with older versions)
    std::string sep = "\t";

    //Setting columns names (optional)
    std::string HODFileHeader = "(1)galaxyTag" + sep +
                                  "(2)X [Mpc/h]" + sep +
                                  "(3)Y [Mpc/h]" + sep +
                                  "(4)Z [Mpc/h]" + sep +
                                  "(5)M_vir (M_sub)[M_sun/h]" + sep +
                                  "(6)M_infall [M_sun/h]" + sep +
                                  "(7)M_star [M_sun/h]";

    // Write HOD catalogue to file
    std::string HodCatalogueFileName = "Galaxies_1e8_Moster_test.dat";
    HOD.write_data(output_path + HodCatalogueFileName, varsToPrint, sep, HODFileHeader);

  }
  catch(cbl::glob::Exception &exc) {std::cerr << exc.what() << endl; exit(1);}

  return 0;
}
