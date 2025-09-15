// =================================================================
// Example code: how to generate a mesh grid, with different methods
// =================================================================

#include <typeinfo>
#include <iostream>

#include "MeshGrid.h"
#include "Catalogue.h"

using namespace std;

int main () {
    
  try {

    // halo ASCII catalogue
    string file_tracers = "../input/halo_catalogue.txt";

    // vector containing the variable name list to read from file
    vector<cbl::catalogue::Var> var_names_tracers = {cbl::catalogue::Var::_X_, cbl::catalogue::Var::_Y_, cbl::catalogue::Var::_Z_};
        
    // vector containing the column corresponding to each attribute
    vector<int> columns_tracers = {1, 2, 3};

    // catalogue constructor
    cbl::catalogue::Catalogue tracers_catalogue = cbl::catalogue::Catalogue(cbl::catalogue::ObjectType::_Halo_, cbl::CoordinateType::_comoving_, var_names_tracers, columns_tracers, {file_tracers}, 0);  
    
    // store the mean particle separation of the simulation
    double mps = tracers_catalogue.mps();      

    
    // -------------------------------------------------------------
    // -------------------- Creation of the objects ----------------
    // -------------------------------------------------------------
    
    cout << "\nCHAINMESH AND MESHGRID CREATION:" << endl;

    auto start1 = chrono::high_resolution_clock::now();
    
    // generate the chain mesh of the input tracer catalogue
    cbl::chainmesh::ChainMesh3D ChM(2*mps, tracers_catalogue.var(cbl::catalogue::Var::_X_), tracers_catalogue.var(cbl::catalogue::Var::_Y_), tracers_catalogue.var(cbl::catalogue::Var::_Z_), 20*mps);
    auto end1 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end1 - start1;
    cout << "Chainmesh creation time: " << elapsed.count() << " seconds" << endl;

    auto start2 = chrono::high_resolution_clock::now();
    
    // generate an equivalent mesh grid of the input tracer catalogue
    cbl::meshgrid::MeshGrid MG = cbl::meshgrid::MeshGrid(tracers_catalogue.var(cbl::catalogue::Var::_X_), tracers_catalogue.var(cbl::catalogue::Var::_Y_), tracers_catalogue.var(cbl::catalogue::Var::_Z_), 2*mps, {});
    auto end2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed2 = end2 - start2;
    cout << "MeshGrid creation time: " << elapsed2.count() << " seconds\n" << endl;

    
    // -------------------------------------------------------------
    // ---------------------- Closer object test -------------------
    // -------------------------------------------------------------

    cout << "CLOSER OBJECT TEST:" << endl;

    start1 = chrono::high_resolution_clock::now();
    
    // find the closer object from a given point with the Chainmesh method
    ChM.get_searching_region(2*mps);
    double dist = 1e30;
    unsigned int finalInd = 0;
    vector<long> close = ChM.close_objects({350., 280., 500.});
    for (auto& part : close) {
      double dist_temp = sqrt((tracers_catalogue.xx(part)-350.)*(tracers_catalogue.xx(part)-350.)+(tracers_catalogue.yy(part)-280)*(tracers_catalogue.yy(part)-280)+(tracers_catalogue.zz(part)-500)*(tracers_catalogue.zz(part)-500));
      if (dist_temp < dist) {
	dist = dist_temp;
	finalInd = part;
      }
    }
    end1 = chrono::high_resolution_clock::now();
    elapsed = end1 - start1;
    cout << "Chainmesh --- object index: " << finalInd << ", time used: " << elapsed.count() << " seconds" << endl;

    start2 = chrono::high_resolution_clock::now();
    
    // find the closer object from a given point with the Mesh Grid method
    unsigned int co = MG.closerObject(350., 280., 500.);
    end2 = chrono::high_resolution_clock::now();
    elapsed2 = end2 - start2;
    cout << "MeshGrid  --- object index: " << co << ", time used: " << elapsed2.count() << " seconds\n" << endl;

    
    // -------------------------------------------------------------
    // ---------------------- Close objects test -------------------
    // -------------------------------------------------------------

    cout << "CLOSE OBJECTS TEST:" << endl;

    start1 = chrono::high_resolution_clock::now();
    
    // find the close objects from a given point with the Chainmesh method
    ChM.get_searching_region(15*mps);
    vector<unsigned int> finalInd2;
    close = ChM.close_objects({450., 450., 500.});
    for (auto& part : close) {
      double dist_temp = sqrt((tracers_catalogue.xx(part)-450.)*(tracers_catalogue.xx(part)-450.)+(tracers_catalogue.yy(part)-450)*(tracers_catalogue.yy(part)-450)+(tracers_catalogue.zz(part)-500)*(tracers_catalogue.zz(part)-500));
      if (dist_temp <= 15*mps && dist_temp>12.5*mps) {
	finalInd2.push_back(part);
      }
    }
    end1 = chrono::high_resolution_clock::now();
    elapsed = end1 - start1;
    cout << "Chainmesh --- number of objects between 12.5 and 15 MPS: " << finalInd2.size() << ", time used: " << elapsed.count() << " seconds" << endl;

    start2 = chrono::high_resolution_clock::now();
    
    // find the close objects from a given point with the Mesh Grid method
    vector<unsigned int> co2 = MG.closeObjects(450., 450., 500., 15*mps, 12.5*mps);
    end2 = chrono::high_resolution_clock::now();
    elapsed2 = end2 - start2;
    cout << "MeshGrid  --- number of objects between 12.5 and 15 MPS: " << co2.size() << ", time used: " << elapsed2.count() << " seconds\n" << endl;

    
    // -------------------------------------------------------------
    // -------------------- N closer objects test ------------------
    // -------------------------------------------------------------

    cout << "N CLOSER OBJECTS TEST:" << endl;

    start1 = chrono::high_resolution_clock::now();
    
    // find the closest N (5) objects from a given point with the Chainmesh method
    ChM.get_searching_region(4*mps);
    vector<pair<double, unsigned int>> objects;
    close = ChM.close_objects({450., 450., 500.});
    for (auto& part : close) 
      objects.emplace_back(sqrt((tracers_catalogue.xx(part)-450.)*(tracers_catalogue.xx(part)-450.)+(tracers_catalogue.yy(part)-450)*
				(tracers_catalogue.yy(part)-450)+(tracers_catalogue.zz(part)-500)*(tracers_catalogue.zz(part)-500)), part);
    
    sort(objects.begin(), objects.end(), [](const auto &a, const auto &b) { return a.first < b.first; });

    objects.resize(5);        

    vector<unsigned int> fIndexes(5);
    for (unsigned int i=0; i<5; i++) 
      fIndexes[i] = objects[i].second;
    end1 = chrono::high_resolution_clock::now();
    
    elapsed = end1 - start1;
    cout << "Chainmesh --- Closest 5 objects from the point: " << endl;
    cbl::Print(fIndexes);
    cout << "time used: " << elapsed.count() << " seconds\n" << endl;

    start2 = chrono::high_resolution_clock::now();
    
    // find the closest N (5) objects from a given point with the Mesh Grid method
    vector<unsigned int> fIndexes2 = MG.N_closerObjects(5, 450., 450., 500.);
    end2 = chrono::high_resolution_clock::now();
    elapsed2 = end2 - start2;
    cout << "MeshGrid --- Closest 5 objects from the point: " << endl;
    cbl::Print(fIndexes2);
    cout << "time used: " << elapsed2.count() << " seconds\n" << endl;

  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }
  
  return 0;
}
