// ====================================================================================
// Example code: how to construct a Halo Occupation Distribution (HOD) galaxy catalogue
// ====================================================================================

#include "Catalogue.h"
 
int main () {

  try {

    // --------------------------------------------------------------
    // ---------------- set the cosmological parameters  ------------
    // --------------------------------------------------------------

    const cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck15_};

    
    // ---------------------------------------------------------------------
    // ---------------- read the dark matter halo catalogue ----------------
    // ---------------------------------------------------------------------
      
    // set input and output paths
    const std::string input_path = "../input/";
    const std::string output_path = "../output/";
    const std::string input_file = "haloCat.dat";

    // halo catalogue constructor
    cbl::catalogue::Catalogue haloCat { cbl::catalogue::ObjectType::_HostHalo_, cbl::CoordinateType::_comoving_, { cbl::catalogue::Var::_X_, cbl::catalogue::Var::_Y_, cbl::catalogue::Var::_Z_, cbl::catalogue::Var::_Mass_ }, {1, 2, 3, 4}, {input_path+input_file}, 0, 1., 1., cosmology };


    // --------------------------------------------------------------------------------------
    // ---------------- construct the HOD galaxy catalogue from the halo one ----------------
    // --------------------------------------------------------------------------------------

    cbl::catalogue::Catalogue HOD { haloCat, cosmology, cbl::catalogue::HODType::_Moster10_, 1.e8, false };

    std::cout << cbl::par::col_green << "HOD catalogue MIN/MAX values per variable:" << cbl::par::col_default << std::endl;
    std::cout << "Min X: " << HOD.Min(cbl::catalogue::Var::_X_) << " - " << "Max X: " << HOD.Max(cbl::catalogue::Var::_X_) << std::endl;
    std::cout << "Min Y: " << HOD.Min(cbl::catalogue::Var::_Y_) << " - " << "Max Y: " << HOD.Max(cbl::catalogue::Var::_Y_) << std::endl;
    std::cout << "Min Z: " << HOD.Min(cbl::catalogue::Var::_Z_) << " - " << "Max Z: " << HOD.Max(cbl::catalogue::Var::_Z_) << std::endl;
    std::cout << "Min Mass: " << HOD.Min(cbl::catalogue::Var::_Mass_) << " - " << "Max Mass: " << HOD.Max(cbl::catalogue::Var::_Mass_) << std::endl;
    std::cout << "Total number of galaxies: " << HOD.nObjects() << " - " << "Total number of halos: " << haloCat.nObjects() << std::endl;
    
    const std::vector<double> stats = HOD.stats_var(cbl::catalogue::Var::_Mass_);
    std::cout << std::endl << cbl::par::col_green << "HOD catalogue mass statitics [Mpc/h]:" << cbl::par::col_default << std::endl;
    std::cout << "The mean: " << stats[0] << std::endl;
    std::cout << "The median: " << stats[1] << std::endl;
    std::cout << "The standard deviation: " << stats[2] << std::endl;
    std::cout << "Difference between the third and first quartiles: " << stats[3] << std::endl << std::endl;

    
    // -----------------------------------------------------------------
    // ---------------- write the HOD catalogue to file ----------------
    // -----------------------------------------------------------------

    // output file name
    std::string HodCatalogueFileName = "Galaxies_1e8_Moster_test.dat";
 
    // vector with HOD catalogue parameters
    const std::vector<cbl::catalogue::Var> varsToPrint = { cbl::catalogue::Var::_GalaxyTag_, cbl::catalogue::Var::_X_, cbl::catalogue::Var::_Y_, cbl::catalogue::Var::_Z_, cbl::catalogue::Var::_Mass_, cbl::catalogue::Var::_MassInfall_, cbl::catalogue::Var::_Mstar_ };

    // file separator (optional, dafault is: "   " for compatibility with older versions)
    const std::string sep = "\t";

    // column names (optional)
    const std::string HODFileHeader = "(1)galaxyTag" + sep + "(2)X [Mpc/h]" + sep + "(3)Y [Mpc/h]" + sep + "(4)Z [Mpc/h]" + sep + "(5)M_vir (M_sub)[M_sun/h]" + sep + "(6)M_infall [M_sun/h]" + sep + "(7)M_star [M_sun/h]";
    
    // store the HOD catalogue
    HOD.write_data(output_path+HodCatalogueFileName, varsToPrint, sep, HODFileHeader);
    
  }
  
  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }

  return 0;
}
