// ===========================================
// Example code: how to read/write a fits file
// ===========================================

#include "FITSwrapper.h"
#include "Func.h" 

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;

int main () {

  try {

    // input file to read
    const string input_file = "catalogue.fits";

    // number of the table extension in the FITS file; a FITS can
    // contain more the one extension; the extension 0 is the primary
    // header
    const int next = 1;

    // the value used to fill non-existing columns, if
    // fill_value==par::defaultDouble, an error is raised when a
    // column is not found
    const double fill_value = 1.;

    // name of the columns in the Table
    const vector<string> column_names = {"XX", "YY", "ZZ", "Weight"};
   
    // read the columns from the table searching by names
    vector<vector<double>> table = cosmobl::ccfitswrapper::read_table_fits(input_file, column_names, next, fill_value);
    cout << table[0][0] << " " <<table[1][0] << " " << table[2][0] << " " << table[3][0] << endl;
    
    // the unit of the columns, optional
    const vector<string> column_units = {"Mpc/h", "Mpc/h", "Mpc/h", ""};

    // write the output on a new fits file, adding the weight column
    cosmobl::ccfitswrapper::write_table_fits(cosmobl::par::DirLoc, "catalogue_with_weights.fits", column_names, table, column_units);
  }
  
  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}
