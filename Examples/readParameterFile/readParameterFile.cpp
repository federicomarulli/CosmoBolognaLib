// ===============================================================
// Example code: how to read parameters from a standard *.ini file
// ===============================================================

#include "ReadParameters.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;

int main () {

  try {

    // define the input parameter file name
    string parameter_file = "parameters.ini";

    // construct the object ReadParameter
    cosmobl::glob::ReadParameters parameters {parameter_file};

    // name of a parameter saved in a string variable
    string par_name = "string_parameter";

    // use the find member function of class ReadParameters to find a string parameter named par_name
    cout << "First parameter in file, a generic name: " << parameters.find<string>(par_name) << endl;

    // there is no particular order to read the parameters from the file
    string par_common_name = "_parameter";
    cout << "Third parameter in file, a float: " << parameters.find<float>("float"+par_common_name) << endl;
    cout << "Second parameter in file, an int: " << parameters.find<int>("int"+par_common_name) << endl;
    if (parameters.find<bool>("bool"+par_common_name)) cout << "The bool parameter is true." << endl;
  
  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}

