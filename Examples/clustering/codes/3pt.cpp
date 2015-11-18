// ===================================================
// How to measure the three-point correlation function
// ===================================================

#include "RandomCatalogue.h"
#include "ThreePointCorrelation.h"

using namespace cosmobl;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;


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
  string dir_triplets = dir_output+"triplets/";
  string dir_random_cat = dir_output;
  string dir_2pt = dir_output;
  
  string MK = "mkdir -p "+dir_output+" "+dir_pairs+" "+dir_triplets+" "+dir_random_cat+" "+dir_2pt; if (system(MK.c_str())) {};

  
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

  auto random_catalogue = random_catalogue_box(catalogue, nRandom);

  
  // ----------------------------------------------------------------------------
  // ---------------- measure the three-point correlation function ----------------
  // ----------------------------------------------------------------------------
  
  // parameters of the method
  double side_s = 20.;     // 1st side of the triangle
  double side_u = 2.;      // ratio of the 2nd side of the triangle (u*s)
  double perc = 0.0225;    // tolerance
  int nbin = 15;           // number of bins
  string type_bin = "lin"; // type of binning

  // create the object
  ThreePointCorrelation ThreeP {catalogue, random_catalogue};

  // set the parameters
  ThreeP.setParameters3p(side_s, side_u, perc, type_bin, nbin);
  
  // measure the 3PCF
  ThreeP.measure_Q(dir_triplets, dir_2pt);

  // write the output
  ThreeP.write_Q(dir_output);
  
  return 0;
}

