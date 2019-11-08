// ================================================================================================
// Example code: how to construct an object of class Data1D, used to handle 1D datasets of any type
// ================================================================================================

#include "Data1D.h"

using namespace std;

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;


// =====================================================================


int main () {

  try {

    // the input file
    const string file = cbl::par::DirLoc+"data.dat";
    
    // header lines to skip
    const int n_lines_header = 2;

    
    // --- method I: providing just the name of the input file --- 
    
    // with default columns
    const cbl::data::Data1D data1(file, n_lines_header);
    data1.Print(); cout << endl;
    
    // specifying that the x column is the second one (the y and error
    // columns will be the third and fourth ones, by default)
    const cbl::data::Data1D data2(file, n_lines_header, 2);
    data2.Print(); cout << endl;
    
    // specifying all the columns
    const cbl::data::Data1D data3(file, n_lines_header, 2, {3}, {4});
    data3.Print(); cout << endl;
      
    // reading two datasets (that will be added one after the other in
    // a single dataset)
    const cbl::data::Data1D data4(file, n_lines_header, 2, {3, 4}, {5, 6});
    data4.Print(); cout << endl;

    
    // --- method II: providing the x, y and error vectors ---

    std::ifstream fin(file.c_str()); cbl::checkIO(fin, file);
  
    double X, Y, ERROR, N;
    std::vector<double> x, y, error;

    string line;
    for (int i=0; i<n_lines_header; ++i)
      getline(fin, line);
    
    while (fin >> X >> Y >> ERROR >> N >> N >> N) {
      x.emplace_back(X);
      y.emplace_back(Y);
      error.emplace_back(ERROR);
    }
  
    fin.clear();

    const cbl::data::Data1D data5(x, y, error);
    data5.Print(); cout << endl;
    
  }
  
  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }

  return 0;
}


