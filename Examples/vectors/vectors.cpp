// ========================================
// How to print vectors and remove elements
// ========================================

#include "Func.h"

using namespace cosmobl;

// define the two global variables containing the names of the CosmoBolognaLib and current directories
string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;


int main () {

  // define three vectors
  vector<double> v1 {1, 2, 3, 4, 5}, v2 {2, 3, 4}, v3 {3, 4, 5, 6};
  
  // print the vectors
  cout << endl << " vector v1: " << endl; print(v1); 
  cout << endl << " vector v2: " << endl; print(v2); 
  cout << endl << " vector v3: " << endl; print(v3); 

  // erase three lines of vector v1 and show the result
  vector<int> lines0 {0, 3, 1};
  Erase(v1, lines0); 
  cout << endl << " vector v1 without the lines 0, 1, 3: " << endl; print(v1); 
  
  // define a matrix
  vector<vector<double>> Mat {v1, v2, v3};

  // print the matrix
  cout << endl << " matrix Mat: " << endl; print(Mat); 

  // erase two lines of the matrix and show the result
  vector<int> lines {2, 0};
  Erase_lines(Mat, lines);
  cout << endl << " matrix Mat without the lines 0 and 2: " << endl; print(Mat);

  // erase two columns of the matrix and show the result
  vector<int> columns {0, 1};
  Erase_columns(Mat, columns); 
  cout << endl << " matrix Mat without the columns 0 and 1: " << endl; print(Mat);


  return 0;
}
