// ===================================================================
// Example code: how to construct a catalogue of extragalactic objects
// ===================================================================

#include "TwoPointCorrelation1D_monopole.h"

using namespace cosmobl;
using namespace catalogue;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;


int main () {

  string file_catalogue = "cat.dat";

  ifstream fin; fin.open(file_catalogue.c_str());  checkIO(file_catalogue, 1);
  
  
  // -----------------------------------------------------------------------------------------------------------
  // ---------------- method I : construct a galaxy catalogue directly by reading an input file ----------------
  // ----------------------------------------------------------------------------------------.........----------
  
  Catalogue catalogue1 {_Galaxy_, {file_catalogue}};

  
  // -------------------------------------------------------------------------------------------------------------------------
  // ---------------- method II : store galaxy coordinates in vectors and use them to construct the catalogue ----------------
  // -------------------------------------------------------------------------------------------------------------------------

  double X, Y, Z;
  vector<double> x, y, z;
  
  while (fin >> X >> Y >> Z) {
    x.emplace_back(X);
    y.emplace_back(Y);
    z.emplace_back(Z);
  }
  
  fin.clear();
  
  Catalogue catalogue2 {_Galaxy_, x, y, z};
  
  
  // ------------------------------------------------------------------------------------------------------------------
  // ---------------- method III : construct a vector of galaxies and add them into an empty catalogue ----------------
  // ------------------------------------------------------------------------------------------------------------------

  vector<shared_ptr<Object>> object;

  fin.seekg(ios::beg);
  
  while (fin >> X >> Y >> Z) {
    auto galaxy = make_shared<Galaxy>(X, Y, Z);
    object.emplace_back(galaxy);
  }
  
  fin.clear(); fin.close();

  Catalogue catalogue3; catalogue3.add_objects(object);

  
  cout << "The number of galaxy in catalogue1 is " << catalogue1.nObjects() << endl;
  cout << "The number of galaxy in catalogue2 is " << catalogue2.nObjects() << endl;
  cout << "The number of galaxy in catalogue3 is " << catalogue3.nObjects() << endl;
  
  return 0;
}

