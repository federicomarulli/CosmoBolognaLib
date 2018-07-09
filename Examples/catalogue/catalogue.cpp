// ===================================================================
// Example code: how to construct a catalogue of extragalactic objects
// ===================================================================

#include "Catalogue.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
std::string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;

int main () {

  try {
  
    std::string file_catalogue = "cat.dat";
  
  
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- method I : construct a galaxy catalogue directly by reading an input file ----------------
    // ----------------------------------------------------------------------------------------.........----------
  
    cbl::catalogue::Catalogue catalogue1 {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_comoving_, {file_catalogue}};

    std::cout << "The coordinates of the first galaxy in the catalogue are: " <<
      catalogue1[0]->xx() << ", " << catalogue1[0]->yy() << ", " << catalogue1[0]->zz() << std::endl;

    
    // -----------------------------------------------------------------------------------------------------------------------
    // ------- method II : construct a galaxy catalogue directly from a file, choosing which quantity has to be read ---------
    // -----------------------------------------------------------------------------------------------------------------------

    // std::vector containing the quantities to be read
    std::vector<cbl::catalogue::Var> attribute = {cbl::catalogue::Var::_X_, cbl::catalogue::Var::_Y_, cbl::catalogue::Var::_Z_};

    // std::vector containing the columns corresponding to each quantity
    std::vector<int> column = {1, 2, 3};
    
    cbl::catalogue::Catalogue catalogue2 {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_comoving_, attribute, column, {file_catalogue}};

    
    // -------------------------------------------------------------------------------------------------------------------------
    // ---------------- method III : construct a galaxy catalogue using std::vectors to specify the galaxy properties ---------------
    // -------------------------------------------------------------------------------------------------------------------------

    std::ifstream fin; fin.open(file_catalogue.c_str()); cbl::checkIO(fin, file_catalogue);
  
    double X, Y, Z;
    std::vector<double> x, y, z;
  
    while (fin >> X >> Y >> Z) {
      x.emplace_back(X);
      y.emplace_back(Y);
      z.emplace_back(Z);
    }
  
    fin.clear();
  
    cbl::catalogue::Catalogue catalogue3 {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_comoving_, x, y, z};

  
    // ------------------------------------------------------------------------------------------------------------------
    // ---------------- method IV : construct a std::vector of galaxies and add them into an empty catalogue -----------------
    // ------------------------------------------------------------------------------------------------------------------

    std::vector<std::shared_ptr<cbl::catalogue::Object>> object;
  
    fin.seekg(std::ios::beg);
  
    while (fin >> X >> Y >> Z) {
      cbl::comovingCoordinates coord = {X, Y, Z};
      auto galaxy = std::make_shared<cbl::catalogue::Galaxy>(coord);
      object.emplace_back(galaxy);
    }
  
    fin.clear(); fin.close();

    cbl::catalogue::Catalogue catalogue4; catalogue4.add_objects(object);

  
    std::cout << "The number of galaxy in catalogue1 is " << catalogue1.nObjects() << std::endl;
    std::cout << "The number of galaxy in catalogue2 is " << catalogue2.nObjects() << std::endl;
    std::cout << "The number of galaxy in catalogue3 is " << catalogue3.nObjects() << std::endl;
    std::cout << "The number of galaxy in catalogue4 is " << catalogue4.nObjects() << std::endl;

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

