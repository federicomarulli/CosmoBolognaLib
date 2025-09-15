// =======================================================================================
// Example code: how to compute the matter two-point correlation function with an emulator
// =======================================================================================

#include "Emulator.h"
#include <chrono>

int main () {

  try {

    // 100 triples of input parameters
    const std::vector<double> OmegaM = cbl::linear_bin_vector(1000, 0.25, 0.35);
    const double OmegaDE = 0.7;
    const double As = 2.e-9;

    // the spatial scales
    const std::vector<double> rr = cbl::linear_bin_vector(100, 10., 150.);

    // define the directory where the trained network is stored (in CosmoBolognaLib/External/Emulator_trained/)
    std::string dir_input = "xi_r2_CAMB";

    // define an object emulator, which load all the required parameters
    cbl::emulator::Emulator Emu{dir_input};

    std::vector<std::vector<double>> xi_r2_emu(OmegaM.size());

    // start the time
    auto begin_emu = std::chrono::high_resolution_clock::now();

    // compute the emulated xi
    for (size_t i=0; i<OmegaM.size(); i++)
      xi_r2_emu[i] = Emu.model({OmegaM[i], OmegaDE, As, 0.}, rr);

    // stop the time
    auto end_emu = std::chrono::high_resolution_clock::now();
    std::cout << "Elapsed time to compute 1000 emulated xi = " << std::chrono::duration_cast<std::chrono::milliseconds>(end_emu - begin_emu).count() << "ms" << std::endl;
  }

  catch (cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }

  return 0;
}
