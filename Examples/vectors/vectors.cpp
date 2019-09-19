// ===========================================================
// Example code: how to print std::vectors and remove elements
// ===========================================================

#include "Kernel.h"

int main () {

  try {
    
    // define three std::vectors
    std::vector<double> v1 {1, 2, 3, 4, 5}, v2 {2, 3, 4}, v3 {3, 4, 5, 6};
  
    // print the std::vectors
    std::cout << std::endl << " std::vector v1: " << std::endl;
    cbl::Print(v1); 
    std::cout << std::endl << " std::vector v2: " << std::endl;
    cbl::Print(v2); 
    std::cout << std::endl << " std::vector v3: " << std::endl;
    cbl::Print(v3); 

    // erase three lines of std::vector v1 and show the result
    std::vector<int> lines0 {0, 3, 1};
    cbl::Erase(v1, lines0); 
    std::cout << std::endl << " std::vector v1 without the lines 0, 1, 3: " << std::endl;
    cbl::Print(v1); 
  
    // define a matrix
    std::vector<std::vector<double>> Mat {v1, v2, v3};

    // print the matrix
    std::cout << std::endl << " matrix Mat: " << std::endl;
    cbl::Print(Mat); 

    // erase two lines of the matrix and show the result
    std::vector<int> lines {2, 0};
    cbl::Erase_lines(Mat, lines);
    std::cout << std::endl << " matrix Mat without the lines 0 and 2: " << std::endl;
    cbl::Print(Mat);

    // erase two columns of the matrix and show the result
    std::vector<int> columns {0, 1};
    cbl::Erase_columns(Mat, columns); 
    std::cout << std::endl << " matrix Mat without the columns 0 and 1: " << std::endl;
    cbl::Print(Mat);
    
  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }

  return 0;
}
