// ===================================================================
// Example code: how to construct a catalogue of extragalactic objects
// ===================================================================

#include "Catalogue.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;

int main () {

  try {
  
    string file_catalogue = "cat.dat";
  
  
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- method I : construct a galaxy catalogue directly by reading an input file ----------------
    // ----------------------------------------------------------------------------------------.........----------
  
    cosmobl::catalogue::Catalogue catalogue1 {cosmobl::catalogue::_Galaxy_, cosmobl::_comovingCoordinates_, {file_catalogue}};

  
    // -------------------------------------------------------------------------------------------------------------------------
    // ---------------- method II : store galaxy coordinates in vectors and use them to construct the catalogue ----------------
    // -------------------------------------------------------------------------------------------------------------------------

    ifstream fin; fin.open(file_catalogue.c_str()); cosmobl::checkIO(fin, file_catalogue);
  
    double X, Y, Z;
    vector<double> x, y, z;
  
    while (fin >> X >> Y >> Z) {
      x.emplace_back(X);
      y.emplace_back(Y);
      z.emplace_back(Z);
    }
  
    fin.clear();
  
    cosmobl::catalogue::Catalogue catalogue2 {cosmobl::catalogue::_Galaxy_, cosmobl::_comovingCoordinates_, x, y, z};

  
    // ------------------------------------------------------------------------------------------------------------------
    // ---------------- method III : construct a vector of galaxies and add them into an empty catalogue ----------------
    // ------------------------------------------------------------------------------------------------------------------

    vector<shared_ptr<cosmobl::catalogue::Object>> object;
  
    fin.seekg(ios::beg);
  
    while (fin >> X >> Y >> Z) {
      cosmobl::comovingCoordinates coord = {X, Y, Z};
      auto galaxy = make_shared<cosmobl::catalogue::Galaxy>(coord);
      object.emplace_back(galaxy);
    }
  
    fin.clear(); fin.close();

    cosmobl::catalogue::Catalogue catalogue3; catalogue3.add_objects(object);

  
    cout << "The number of galaxy in catalogue1 is " << catalogue1.nObjects() << endl;
    cout << "The number of galaxy in catalogue2 is " << catalogue2.nObjects() << endl;
    cout << "The number of galaxy in catalogue3 is " << catalogue3.nObjects() << endl;

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}

