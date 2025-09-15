#include "MapCalculator.h"
#include "Catalogue.h"
#include "Cosmologies.h"
#include "Kernel.h"

using namespace cbl; 

// defining a function to read the mask from a file
template<typename T>
std::vector<T> readColumnFromFile (const std::string& filename, const int columnIndex) {
  std::ifstream file(filename);
  if (!file.is_open()) 
    throw std::runtime_error("Error in opening " + filename);
  std::cout << "Reading mask file..." << std::endl;
  std::vector<T> column;
  std::string line;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    T value;
    int currentColumnIndex = 0;
    while (iss >> value) {
      if (currentColumnIndex==columnIndex) {
        column.push_back(value);
        break;
      }
      currentColumnIndex++;
    }
  }
  file.close();
  return column;
}

int main() {
    
  try {

    // reading data from file

    std::string file_catalogue = "input_example.dat";
    const auto cosmo = std::make_shared<cbl::cosmology::LCDM>("Planck18");
    cbl::catalogue::Catalogue cat (cbl::catalogue::ObjectType::_Mock_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmo, cbl::CoordinateUnits::_degrees_ );

    // setting up the class

    int ns_value = 64;
    Healpix_Ordering_Scheme scheme = RING; 
    int nsg_value = 3;
    std::vector<double> zobs = cbl::arange(0.3, 0.7, 0.05);
    std::vector<double> sgzar (1, 0.15);
    cbl::measure::MapCalculator map (ns_value, scheme, nsg_value, zobs, sgzar, cat);

    // reading the mask from file

    std::string mask = "masks.dat";
    std::vector<double> mskin = readColumnFromFile<double>(mask, 0);

    // generating some maps

    std::vector<cbl::Tensor3Dd> hpx_maps = {map.ADF(mskin), std::get<1>(map.ADF_ARF(mskin)), map.SN(mskin, 1)[0][0]};

    // printing all the maps

    std::vector<std::string> fout =  {"ADF.dat", "ARF.dat", "SN00.dat"};

    for (size_t i = 0; i < fout.size(); ++i) {
      std::cout << "Printing on "<< fout[i] << std::endl;
      map.print_map(hpx_maps[i], fout[i]);
    }

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); } 
  
  return 0;
}