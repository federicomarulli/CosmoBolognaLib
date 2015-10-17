// ===================================================================================
// How to measure the two-point correlation function and estimate the jackknife errors
// ===================================================================================

#include "RandomCatalogue.h"
#include "GlobalFunc.h"

using namespace cosmobl;

const string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;


int main () {
  
  // -----------------------------------------------------------------------------------------------
  // ---------------- set the cosmological parameters and create the object 'cosmology' ------------
  // -----------------------------------------------------------------------------------------------

  double OmegaM = 0.25;
  double Omega_b = 0.045;
  double Omega_nu = 0.;
  double massless_neutrinos = 3.04;
  int    massive_neutrinos = 0; 
  double OmegaL = 1.-OmegaM;
  double Omega_radiation = 0.;
  double hh = 0.73;
  double scalar_amp = 2.742e-9;
  double n_s = 1;
  double wa = 0.;
  double w0 = -1.;   

  Cosmology cosmology {OmegaM, Omega_b, Omega_nu, massless_neutrinos, massive_neutrinos, OmegaL, Omega_radiation, hh, scalar_amp, n_s, w0, wa};

  
  // --------------------------------------------------------------------
  // ---------------- Input/Output files and directories ----------------
  // --------------------------------------------------------------------
  
  string file_catalogue = par::DirLoc+"../input/cat.dat";

  string dir_output = par::DirLoc+"../output/";
  string dir_pairs = dir_output+"pairs/";
  string dir_random_cat = dir_output;
  string dir_covariance = dir_output+"covariance/";
  
  string MK = "mkdir -p "+dir_output+" "+dir_pairs+" "+dir_covariance; if (system(MK.c_str())) {};

  
  // -------------------------------------------------------------------------------------------
  // ---------------- read the input catalogue and create the object 'catalogue'----------------
  // -------------------------------------------------------------------------------------------

  cout << "I'm reading the input catalogue..." << endl;
  
  ifstream fin(file_catalogue.c_str()); checkIO(file_catalogue, 1);
  
  double RA, DEC, RED;
  
  vector<shared_ptr<Object>> data_obj;

  while (fin >> RA >> DEC >> RED) {
    shared_ptr<Galaxy> gal(new Galaxy{RA, DEC, RED, cosmology});
    data_obj.push_back(gal);
  }
  fin.close();
  
  shared_ptr<Catalogue> catalogue(new Catalogue{data_obj});


  // ---------------------------------------------------------------------------------------------------------
  // ---------------- construct the random catalogue and create the object 'random_catalogue' ----------------
  // ---------------------------------------------------------------------------------------------------------

  double N_R = 1.;

  int nRandom = (int)catalogue->nObjects()*N_R;

  auto random_catalogue = random_catalogue_box(catalogue, nRandom, dir_random_cat);
  

  // -----------------------------------------------------------------------------
  // ---------------- subdivide the two catalogues in sub-regions ----------------
  // -----------------------------------------------------------------------------

  // number of sub-regions in which the catalogue is divided (in each dimension)
  int nx = 3, ny = 3, nz = 3;
  
  // assign an index to each object to identify the sub-region where the object is located
  set_ObjectRegion_SubBoxes(catalogue, random_catalogue, nx, ny, nz); 
  
  
  // ----------------------------------------------------------------------------
  // ---------------- measure the two-point correlation function ----------------
  // ----------------------------------------------------------------------------
  
  // parameters of the method
  double rMIN = 1.;
  double rMAX = 50.;
  double logbinSize = 0.05;
  double binSize = 0.5;
  double cosSize = 0.02;

  // create the object used to measure the two-point correlation function
  TwoPointCorrelation TwoP {catalogue, random_catalogue};

  // set the parameters
  TwoP.setParameters(rMIN, rMAX, logbinSize, binSize, cosSize);

  // measure the two-point correlation function
  TwoP.measure_xi(dir_pairs, dir_covariance);

  // store the output data
  TwoP.write_xi(dir_output);

  return 0;
}

