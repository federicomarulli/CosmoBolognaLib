// =====================================================================================
// Example code: how to model the bias from the projected two-point correlation function
// =====================================================================================

#include "TwoPointCorrelation1D_monopole.h"
#include "Modelling_TwoPointCorrelation.h"
#include "GlobalFunc.h"

using namespace cosmobl;
using namespace catalogue;
using namespace statistics;
using namespace twopt;
using namespace modelling;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;


int main () {
  
  // -----------------------------------------------------------------
  // ---------------- use default cosmological parameters ------------
  // -----------------------------------------------------------------
  
  Cosmology cosmology;

  
  // --------------------------------------------------------------------
  // ---------------- Input/Output files and directories ----------------
  // --------------------------------------------------------------------
  
  string file_catalogue = par::DirLoc+"../input/cat.dat";

  string dir_output = par::DirLoc+"../output/";
  string dir_pairs = dir_output+"pairs/";
  string dir_covariance = dir_output+"covariance/";
  
  string MK = "mkdir -p "+dir_output+" "+dir_pairs+" "+dir_covariance; if (system(MK.c_str())) {};

  
  // -----------------------------------------------------------------------------------------------------------
  // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
  // -----------------------------------------------------------------------------------------------------------

  cout << "I'm reading the input catalogue..." << endl;

  Catalogue catalogue {_Galaxy_, _observedCoordinates_, {file_catalogue}, cosmology};

  
  // ----------------------------------------------------------------
  // ---------------- construct the random catalogue ----------------
  // ----------------------------------------------------------------

  double N_R = 1.; // random/data ratio
  
  Catalogue random_catalogue {_createRandom_box_, catalogue, N_R};


  // --------------------------------------------------------------------------------------------
  // ---------------- measure the monopole of the two-point correlation function ----------------
  // --------------------------------------------------------------------------------------------

  // binning parameters

  double rMin = 1.;            // minimum separation 
  double rMax = 50.;           // maximum separation 
  int nbins = 20;              // number of bins
  double shift = 0.5;          // spatial shift used to set the bin centre 
  double piMax_integral = 30.; // upper limit of the integral

  
  // measure the projected two-point correlation function and estimate poissonian error

  auto TwoP = TwoPointCorrelation::Create(TwoPType::_1D_projected_, catalogue, random_catalogue, _logarithmic_, rMin, rMax, nbins, shift, rMin, rMax, nbins, shift, piMax_integral);

  TwoP->measure(ErrorType::_Poisson_, dir_output, {dir_output});
  TwoP->write(dir_output, "wp");

  
  // --------------------------------------------------------------------------------------------
  // -------------------------- model the bias from projected 2pcf ------------------------------
  // --------------------------------------------------------------------------------------------

  double redshift = 1.; // redshift of the sample

  double bias_value = 1.2; // guess value for the bias

  vector<double> bias_limits = {0.8, 3.}; 
  Prior bias_prior { PriorType::_UniformPrior_, bias_limits[0], bias_limits[1] }; // flat prior for the bias
  
  int nChains = 100;    // number of chains
  int chain_size = 1000; // size of the chains

  vector<double> fit_limits = {10., 40.}; // limits of the fit

  auto model_twop = Modelling_TwoPointCorrelation::Create(TwoP, redshift, cosmology);
  model_twop->fit_bias(statistics::_GaussianLikelihood_Error_, fit_limits, bias_value, bias_prior, nChains, chain_size, dir_output);

  return 0;
}

