#include "Catalogue.h"
#include "Cosmology.h"

using namespace cbl;
using namespace std;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;

 
int main () {

  try {
    
    const double OmM = 0.2711;
    const double Omb = 0.0451;
    const double Omn = 0.;
    const double massless = 3.04;
    const int massive = 0;
    const double Omr = 0.;
    const double OmL = 0.7289;
    const double hh = 0.703;
    const double As = 2.194e-9;
    const double pivot = 0.05;
    const double ns = 0.96;
    const double w0 = -1.;
    const double wa = 0.;

    cbl::cosmology::Cosmology cosm {OmM, Omb, Omn, massless, massive, OmL, Omr, hh, As, pivot, ns, w0, wa};
    cosm.set_sigma8(0.809);

    vector<catalogue::Var> attributes = {catalogue::Var::_Mass_, catalogue::Var::_X_, catalogue::Var::_Y_, catalogue::Var::_Z_};
    vector<int> column = {3, 5, 6, 7}; 

    const string file_catalogue = "./LCDM_CoDECS_Groups_092.txt";
    const int comments = 0;
    const double nSub = 1.1; 
    const double fact = 0.001;
    
    catalogue::Catalogue catalogue_input {catalogue::ObjectType::_HostHalo_, CoordinateType::_comoving_, attributes, column, {file_catalogue}, comments, nSub, fact};
    cout << "Catalogue read! " << endl;
    
    for (size_t i=0; i<catalogue_input.nObjects(); i++)
      catalogue_input.set_var(i, catalogue::Var::_Mass_, catalogue_input.var(i, catalogue::Var::_Mass_)*1.e10);    

    cout << "\n X MIN : " << catalogue_input.Min(catalogue::Var::_X_) << " -- X MAX " << catalogue_input.Max(catalogue::Var::_X_) << endl;
    cout << " Y MIN : " << catalogue_input.Min(catalogue::Var::_Y_) << " -- Y MAX " << catalogue_input.Max(catalogue::Var::_Y_) << endl;
    cout << " Z MIN : " << catalogue_input.Min(catalogue::Var::_Z_) << " -- Z MAX " << catalogue_input.Max(catalogue::Var::_Z_) << endl;
    cout << " Mass MIN : " << catalogue_input.Min(catalogue::Var::_Mass_) << " -- Mass MAX " << catalogue_input.Max(catalogue::Var::_Mass_) << endl << endl;

    catalogue::Catalogue catalogue_output {catalogue_input, cosm, catalogue::HODType::_Zehavi05_, -19., false};
    
  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}


