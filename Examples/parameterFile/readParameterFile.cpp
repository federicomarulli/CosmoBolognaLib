// ===============================================================
// Example code: how to read parameters from a standard *.ini file
// ===============================================================

#include "ReadParameters.h"

int main () {

  try {

    // define the input parameter file name
    std::string parameter_file = "parameters.ini";

    // construct the object ReadParameter
    cbl::glob::ReadParameters parameters {parameter_file};

    // name of a parameter saved in a std::string variable
    std::string par_name = "string_parameter";

    // use the find member function of class ReadParameters to find a std::string parameter named par_name
    std::cout << "First parameter in file, a generic name: " << parameters.find<std::string>(par_name) << std::endl;

    // there is no particular order to read the parameters from the file
    std::string par_common_name = "_parameter";
    std::cout << "Third parameter in file, a float: " << parameters.find<float>("float"+par_common_name) << std::endl;
    std::cout << "Second parameter in file, an int: " << parameters.find<int>("int"+par_common_name) << std::endl;
    if (parameters.find<bool>("bool"+par_common_name)) std::cout << "The bool parameter is true." << std::endl;
    
    // also std::vectors can be read
    std::vector<int> vect = parameters.find_vector<int> ("vector"+par_common_name);
    std::cout << "The fifth parameter in file is a std::vector which contains " << vect.size() << " integers." << std::endl;
   
  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

