/********************************************************************
 *  Copyright (C) 2022 by Federico Marulli, Simone Sartori          *
 *  Sofia Contarini, Carlo Cannarozzo and Tommaso Ronconi           *
 *  federico.marulli3@unibo.it simone.sartori5@studio.unibo.it      *
 *  sofia.contarini3@unibo.it carlo.cannarozzo@studio.unibo.it      *
 *  tommaso.ronconi@studio.unibo.it                                 *
 *                                                                  *
 *  This program is free software; you can redistribute it and/or   *
 *  modify it under the terms of the GNU General Public License as  *
 *  published by the Free Software Foundation; either version 2 of  *
 *  the License, or (at your option) any later version.             *
 *                                                                  *
 *  This program is distributed in the hope that it will be useful, *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   *
 *  GNU General Public License for more details.                    *
 *                                                                  *
 *  You should have received a copy of the GNU General Public       *
 *  License along with this program; if not, write to the Free      *
 *  Software Foundation, Inc.,                                      *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.       *
 ********************************************************************/

/**
 *  @file Catalogue/VoidCatalogue.cpp
 *
 *  @brief Methods of the class Catalogue to construct Void catalogues
 *
 *  This file contains the implementation of the methods of the class
 *  Catalogue used to create Void catalogues
 *
 *  @author Simone Sartori, Sofia Contarini, Federico Marulli, Carlo Cannarozzo, Tommaso Ronconi
 *
 *  @author simone.sartori5@studio.unibo.it, sofia.contarini3@unibo.it, 
 *  federico.marulli3@unibo.it, carlo.cannarozzo@studio.unibo.it, tommaso.ronconi@studio.unibo.it
 */


#include "Func.h"
#include "Catalogue.h"
#include "ChainMesh_Catalogue.h"
#include "Object.h"
#include "Catalogue.h"
#include "Data1D.h"
#include "Posterior.h"
#include "TwoPointCorrelation.h"
#include "TwoPointCorrelation1D_monopole.h"
#include "CatalogueChainMesh.h"
#include <omp.h>
#include <chrono>

using namespace std;
using namespace cbl;

using namespace catalogue; 

cbl::catalogue::Catalogue::Catalogue (const VoidAlgorithm algorithm, Catalogue tracer_catalogue, Catalogue random_catalogue, const std::string dir_output, const std::string output, const double cellsize, const int n_rec, const double step_size, const double threshold)
{
  // -------------------------------------------- //
  // ---------------- First Step ---------------- //
  // -------------------------------------------- //
    
  const double start_time = omp_get_wtime();
  shared_ptr<Catalogue> halo_cat, rand_cat;
  vector<Catalogue> displacement_catalogue(n_rec);
  vector<vector<unsigned int>> randoms;
  cosmology::Cosmology cosm;

  if (algorithm==VoidAlgorithm::_LaZeVo_) {

    Catalogue rand_catalogue; 

    long int time = (long int)chrono::steady_clock::now().time_since_epoch().count()%1000000;
    int rndd = time; 

    if (random_catalogue.nObjects()==0) 
    {
      rand_catalogue = Catalogue(RandomType::_createRandom_box_, tracer_catalogue, 1., 10, cosm, false, 10., {}, {}, {}, 10, rndd);
      coutCBL << "... random catalogue created!" << endl;
    }
    else if (random_catalogue.nObjects() == tracer_catalogue.nObjects())
    {
      coutCBL << "Reading random catalogue of " << random_catalogue.nObjects() << " objects provides by user..." << endl;
      rand_catalogue = random_catalogue;
      coutCBL << "...Done!" << endl;
    }
    else {
      coutCBL << "Reading random catalogue of " << random_catalogue.nObjects() << " objects provides by user..." << endl;
      cout << endl;
      int diff = random_catalogue.nObjects() - tracer_catalogue.nObjects();
      if (random_catalogue.nObjects() > tracer_catalogue.nObjects()) coutCBL << "WARNING! Random catalogue has " << 
          diff << " more objects than the tracers catalogue. The number will be equalized by the function catalogue::equalize_random_LC." << endl;
      if (random_catalogue.nObjects() < tracer_catalogue.nObjects()) coutCBL << "WARNING! Random catalogue has " << 
          -diff << " less objects than the tracers catalogue. The number will be equalized by the function catalogue::equalize_random_LC." << endl;      
      rand_catalogue.equalize_random_box(tracer_catalogue, rndd);
      coutCBL << "...Done!" << endl;
    }

    // shared pointers to the tracer and the random catalogues
    rand_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(rand_catalogue)));
    halo_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(tracer_catalogue)));  

    // maximum for the number of objects in the tracer catalogue (unsigned int) 
    unsigned int num_objects = ((halo_cat->nObjects())%2==0) ? halo_cat->nObjects() : (halo_cat->nObjects())-1;

    if ((halo_cat->nObjects())%2 != 0) {
      rand_cat -> remove_object(num_objects);
      halo_cat -> remove_object(num_objects);
    }

    // vector for the indices of haloes
    std::vector<unsigned int> index_halo_cat(num_objects);
    std::iota (std::begin(index_halo_cat), std::end(index_halo_cat), 0); 
    cout << "num obj: " << num_objects << endl;
    cout << endl;
    coutCBL << "Creating ChainMesh... " << flush;
    // chainmesh setting
    CatalogueChainMesh ChainMesh_haloes = CatalogueChainMesh(cbl::Euclidean_distance, cellsize, halo_cat);
    CatalogueChainMesh ChainMesh_randoms = CatalogueChainMesh(cbl::Euclidean_distance, cellsize, rand_cat);
    unsigned int N_near_obj = 3*3*4*par::pi;
    vector<vector<unsigned int>> near_part = ChainMesh_haloes.N_nearest_objects_cat(N_near_obj); //the selected halo (that will have distance = 0) and the nearest 3))
		
		randoms.resize(n_rec, vector<unsigned int>(num_objects)); //haloes fixed
    for(int rec=0; rec<n_rec; rec++) {
      std::iota (std::begin(randoms[rec]), std::end(randoms[rec]), 0); 
      auto time2 = (long int)chrono::steady_clock::now().time_since_epoch().count();
      time2 = (time2%1000000);
      std::shuffle(randoms[rec].begin(), randoms[rec].end(), default_random_engine(time2));
    }

    cout << " Done! " << endl;
    cout << endl;	
    coutCBL << "Setting starting configuration... " << endl;

    double dist = 4*halo_cat->mps();
    double minX = ChainMesh_haloes.lim()[0][0]+dist/10, maxX = ChainMesh_haloes.lim()[0][1]-dist/10;
    double minY = ChainMesh_haloes.lim()[1][0]+dist/10, maxY = ChainMesh_haloes.lim()[1][1]-dist/10;
    double minZ = ChainMesh_haloes.lim()[2][0]+dist/10, maxZ = ChainMesh_haloes.lim()[2][1]-dist/10;

    for (int rec=0; rec<n_rec; rec++) {    
        
      auto timeS = ((long int)chrono::steady_clock::now().time_since_epoch().count())%1000000;
      random::UniformRandomNumbers randX(minX, maxX, timeS);
      random::UniformRandomNumbers randY(minY, maxY, timeS+1);
      random::UniformRandomNumbers randZ(minZ, maxZ, timeS+2);

      CatalogueChainMesh ChM_halo_copy = ChainMesh_haloes;
      CatalogueChainMesh ChM_rand_copy = ChainMesh_randoms;
      unsigned int removed = 0;

      while (removed < num_objects) {
        vector<double> pos = {randX(), randY(), randZ()};
        vector<unsigned int> close=ChM_halo_copy.Close_objects(pos, dist);
        unsigned int tr_to_rmv = std::min(100, (int)close.size());
        if (tr_to_rmv > 0) {
          vector<unsigned int> close_random = ChM_rand_copy.N_nearest_objects(pos, close.size());
          std::random_device dev;
          std::mt19937 rng(dev());
          removed += tr_to_rmv;
          while (tr_to_rmv > 0) {
            std::uniform_int_distribution<std::mt19937::result_type> dist(0, tr_to_rmv-1);
            int rnd1 = dist(rng), rnd2 = dist(rng);
            randoms[rec][close[rnd1]] = close_random[rnd2];
            ChM_halo_copy.deletePart(close[rnd1]);
            ChM_rand_copy.deletePart(close_random[rnd2]);
            close.erase(close.begin()+rnd1);
            close_random.erase(close_random.begin()+rnd2);
            tr_to_rmv--;
          }
        }
      }
    }    

    coutCBL << "Performing the iterations... " << endl;
    cout << endl;

    for (int rec=0; rec<n_rec; rec++) {
      auto index_halo_cat_copy = index_halo_cat;
      double ratio = 1.;
      unsigned int n_iter = 0;

      while (ratio > threshold) {

        auto time2 = (long int)chrono::steady_clock::now().time_since_epoch().count()%1000000;
        std::shuffle(index_halo_cat_copy.begin(), index_halo_cat_copy.end(), default_random_engine(time2));
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> dist(1,N_near_obj-1);
        vector<bool> index_bool(num_objects, false), used(num_objects, false);
    
#pragma omp parallel num_threads(omp_get_max_threads())
        {     
#pragma omp for schedule(static)
          for (size_t i=0; i<num_objects; i++) {
            vector<unsigned int> H(4), R_def(4), R(4), R_copy(4);
            unsigned int rand1 = dist(rng), rand2 = dist(rng), rand3 = dist(rng);
            while (rand1 == rand2 || rand1==rand3 || rand2 == rand3) {
              rand2 = dist(rng);
              rand3 = dist(rng);
            }
            if (used[near_part[index_halo_cat_copy[i]][0]] == false && used[near_part[index_halo_cat_copy[i]][rand1]] == false && 
            used[near_part[index_halo_cat_copy[i]][rand2]] == false && used[near_part[index_halo_cat_copy[i]][rand3]] == false) { 
            // selezione particelle coinvolte
              used[near_part[index_halo_cat_copy[i]][0]] = true;
              used[near_part[index_halo_cat_copy[i]][rand1]] = true;
              used[near_part[index_halo_cat_copy[i]][rand2]] = true;
              used[near_part[index_halo_cat_copy[i]][rand3]] = true;
              H[0]=near_part[index_halo_cat_copy[i]][0];
              H[1]=near_part[index_halo_cat_copy[i]][rand1];
              H[2]=near_part[index_halo_cat_copy[i]][rand2];
              H[3]=near_part[index_halo_cat_copy[i]][rand3];
              R[0]=randoms[rec][H[0]];
              R[1]=randoms[rec][H[1]];
              R[2]=randoms[rec][H[2]];
              R[3]=randoms[rec][H[3]];

              R_copy = R;
        // calcolo configurazione minima
              double dist_min = 1000000.;
              double dist;
              do {
                dist = 0.;
                for (size_t j=0; j<R.size(); j++) 
                  dist += cbl::Euclidean_distance(halo_cat->xx(H[j]), rand_cat->xx(R[j]), halo_cat->yy(H[j]), rand_cat->yy(R[j]), halo_cat->zz(H[j]), rand_cat->zz(R[j]));
                
                if (dist < dist_min) {
                  dist_min = dist;
                  R_def = R;
                } 
              } while(std::next_permutation(R.begin(), R.end()));
        
              //aggiornamento vettori haloes e randoms
              if(R_def != R_copy) index_bool[i] = true; 

              for (size_t j=0; j<4; j++) randoms[rec][H[j]] = R_def[j];

              used[H[0]] = false;
              used[H[1]] = false;
              used[H[2]] = false;
              used[H[3]] = false;
            }
            else continue;
          }
        }

        unsigned int changed_couples = count(index_bool.begin(), index_bool.end(), true);

        ratio = (double)changed_couples/num_objects;
        if (n_iter%100==0) cout << "ratio: " << ratio << " " << n_iter << endl;
        n_iter++;
      }

      double distanza = 0.;
      for (unsigned int i=0; i<num_objects; i++) {
      distanza += cbl::Euclidean_distance(halo_cat->xx(i), rand_cat->xx(randoms[rec][i]),
          halo_cat->yy(i), rand_cat->yy(randoms[rec][i]), halo_cat->zz(i), rand_cat->zz(randoms[rec][i]));
      }
      cout << "distanza finale: " << distanza/num_objects << endl;

      vector<shared_ptr<cbl::catalogue::Object>> displ_catalogue_object;

      for (size_t ii = 0; ii<randoms[rec].size(); ii++)
      { 
        cbl::comovingCoordinates displacement_values = {rand_cat->xx(randoms[rec][ii])-halo_cat->xx(ii), rand_cat->yy(randoms[rec][ii])-halo_cat->yy(ii), rand_cat->zz(randoms[rec][ii])-halo_cat->zz(ii)};
        auto halo_displ = make_shared<cbl::catalogue::Halo>(displacement_values);
        displ_catalogue_object.push_back(halo_displ); 
      }
  
////////////////////////////////////////////////////////////////////////////////////////////////////////
      displacement_catalogue[rec].add_objects(displ_catalogue_object);

      if (rec == 0) {
        stringstream ss;
        ss << rec;
        string str = ss.str();
        //       string output = ss.str();
        string name =  "../output/Displacement";
        name.append("_");
        name.append(str);
        name.append("_");
        //       output = "coord_LCDM.dat";
        name.append(output);
        ofstream fout_displ(name.c_str());
        cout << endl;

        coutCBL << "Writing the displacement ... " << endl;
        for (unsigned int i=0; i<num_objects; i++)
        {
          fout_displ.precision(7);
          fout_displ << halo_cat->xx(i) << " " << halo_cat->yy(i) << " " << halo_cat->zz(i) << " "
            << rand_cat->xx(randoms[rec][i]) << " " << rand_cat->yy(randoms[rec][i]) << " " << rand_cat->zz(randoms[rec][i]) << " "
            << displacement_catalogue[rec].xx(i) << " " << displacement_catalogue[rec].yy(i) << " " << displacement_catalogue[rec].zz(i) << endl;   
        } // modified
          
        fout_displ.clear(); fout_displ.close();  
      }
    }


  } // End of LaZeVo method
  else if (algorithm==VoidAlgorithm::_RIVA_) // Begin of RIVA method
  {
      
    cout << endl; coutCBL << par::col_green << "RIVA algorithm" << par::col_default << endl << endl;
    ErrorCBL ("Still not avaiable!", "Catalogue", "VoidCatalogue.cpp");
      /*
	(void)delta_movement;
	// * * * Two Point Correlation Function * * * //

	const double N_R = 1.; // random/data ratio

	const cbl::catalogue::Catalogue test_random_catalogue {cbl::catalogue::_box_, tracer_catalogue, N_R};

	const double rMin  = 1e3; // minimum separation 
	const double rMax  = 5e4; // maximum separation 
	const int nbins    = 50;  // number of bins
	const double shift = 0.5; // spatial shift used to set the bin centre

	string w_file = "xi_RIVA_pre.dat";
	cbl::twopt::TwoPointCorrelation1D_monopole TPCF_catalogue {tracer_catalogue, test_random_catalogue, cbl::_logarithmic_, rMin, rMax, nbins, shift};
	// // //     TPCF_catalogue.measure_for_RIVA(cbl::twopt::_Poisson_, dir_output);
	TPCF_catalogue.measure(cbl::twopt::_Poisson_, dir_output);

	TPCF_catalogue.write(dir_output, w_file);
    
    
	// // // // // WHILE CYCLE - - - while (2pcf_peak/max(xi[i])<1e-2) {...}
    
	auto halo_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(tracer_catalogue)));

	vector<shared_ptr<cbl::catalogue::Object>> final_catalogue_object; // vector shared pointer for the constructor 'final_catalogue'
    
	vector<double> x_position, y_position, z_position;

	for (size_t i = 0; i<num_objects; i++)
	{
	x_position.push_back(halo_cat->xx(i));
	y_position.push_back(halo_cat->yy(i));
	z_position.push_back(halo_cat->zz(i));
	}    

	srand(time(NULL));
	int randomDirection = 0;
    
	for (size_t i = 0; i<x_position.size(); i++)
	{
	randomDirection = (rand() % 3) + (-1);
	x_position[i]=x_position[i]+delta_movement*randomDirection;
	randomDirection = (rand() % 3) + (-1);
	y_position[i]=y_position[i]+delta_movement*randomDirection;
	randomDirection = (rand() % 3) + (-1);  
	z_position[i]=z_position[i]+delta_movement*randomDirection;
	}
    
	vector<shared_ptr<cbl::catalogue::Object>> temp_initial_catalogue_object;

	for (size_t i = 0; i<x_position.size(); i++)
	{
	cbl::comovingCoordinates initial_coordinates = {x_position[i], y_position[i], z_position[i]};
	auto inital_halo_coord = make_shared<cbl::catalogue::Halo>(initial_coordinates);
	temp_initial_catalogue_object.push_back(inital_halo_coord);
	}

	const cbl::catalogue::Catalogue test_random_catalogue_post {cbl::catalogue::_box_, tracer_catalogue, N_R};


	string w_file_post = "xi_RIVA_post.dat";
	cbl::twopt::TwoPointCorrelation1D_monopole TPCF_catalogue_post {temp_initial_catalogue_object, test_random_catalogue_post, cbl::_logarithmic_, rMin, rMax, nbins, shift};

	TPCF_catalogue_post.measure(cbl::twopt::_Poisson_, dir_output);

	TPCF_catalogue_post.write(dir_output, w_file_post);


	exit(1000);
      */
    } // End of RIVA method

  else ErrorCBL ("algorithm type is not correct!", "Catalogue", "VoidCatalogue.cpp");

  // --------------------------------------------- //
  // ---------------- Second Step ---------------- //
  // --------------------------------------------- //

  vector<shared_ptr<Catalogue>> displ_cat(n_rec);
  for (int i=0; i<n_rec; i++) displ_cat[i] = make_shared<Catalogue>(Catalogue(move(displacement_catalogue[i]))); 

  // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
  //   
  //   // * * * Two Point Correlation Function * * * //
  // 
  //   const double N_R = 1.; // random/data ratio
  //   
  //   const cbl::catalogue::Catalogue test_random_catalogue {cbl::catalogue::_box_, final_catalogue, N_R};
  // 
  //   const double rMin  = 1e3; // minimum separation 
  //   const double rMax  = 5e4; // maximum separation 
  //   const int nbins    = 20;  // number of bins
  //   const double shift = 0.5; // spatial shift used to set the bin centre
  //   
  //   string w_file = "xi.dat";
  //   
  //   cbl::twopt::TwoPointCorrelation1D_monopole TPCF_final_catalogue {final_catalogue, test_random_catalogue, cbl::_logarithmic_, rMin, rMax, nbins, shift};
  //       
  //   TPCF_final_catalogue.measure(cbl::twopt::_Poisson_, dir_output);
  //   
  //   TPCF_final_catalogue.write(dir_output, w_file);
  // 
  // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

  // ^ ^ ^ Uncomment the previous part to compute the Two Point Correlation Function of the final displacement with respect to a generic random field ^ ^ ^ //
    
  double volume = halo_cat->volume();
  double density = halo_cat->numdensity(); 
  double step = step_size*halo_cat->mps(); 

  coutCBL << "Volume = " << volume  << endl;
  coutCBL << "n      = " << density << endl;
  coutCBL << "step   = " << step    << endl;
  cout << endl;

  vector<double> min_X = {rand_cat->Min(Var::_X_)-0.5*step, halo_cat->Min(Var::_X_)-0.5*step};
  vector<double> min_Y = {rand_cat->Min(Var::_Y_)-0.5*step, halo_cat->Min(Var::_Y_)-0.5*step};
  vector<double> min_Z = {rand_cat->Min(Var::_Z_)-0.5*step, halo_cat->Min(Var::_Z_)-0.5*step};
  vector<double> min = {*min_element(min_X.begin(),min_X.end()), *min_element(min_Y.begin(),min_Y.end()), *min_element(min_Z.begin(),min_Z.end())};
  vector<double> max_X = {rand_cat->Max(Var::_X_)+0.5*step, halo_cat->Max(Var::_X_)+0.5*step};
  vector<double> max_Y = {rand_cat->Max(Var::_Y_)+0.5*step, halo_cat->Max(Var::_Y_)+0.5*step};
  vector<double> max_Z = {rand_cat->Max(Var::_Z_)+0.5*step, halo_cat->Max(Var::_Z_)+0.5*step};
  vector<double> max = {*max_element(max_X.begin(),max_X.end()), *max_element(max_Y.begin(),max_Y.end()), *max_element(max_Z.begin(),max_Z.end())};

  vector<double> cells(3);
  for(unsigned int&& i=0; i<3; i++) cells[i] = (max[i]-min[i])/step;
  unsigned int nCells = *max_element(cells.begin(), cells.end());
  auto min_copy = min;
  for(unsigned int&& i=0; i<3; i++) {
  	min[i] = min_copy[i] - (step*nCells-(max[i]-min_copy[i]))/2;
  	max[i] = max[i] + (step*nCells-(max[i]-min_copy[i]))/2;
  }
  
  vector<vector<vector<double>>> Divergence(nCells, vector<vector<double>>(nCells, vector<double>(nCells, 0.)));

  cout << endl;
  coutCBL << "Compute the Divergence field ..." << endl;
  cout << endl;

  for(int rec=0; rec<n_rec; rec++) {
    //cycling over the halo positions
    for (size_t&& i=0; i<randoms[rec].size(); i++) {
      vector<vector<double>> coord_part(2);
      coord_part[0] = {rand_cat->xx(randoms[rec][i]), rand_cat->yy(randoms[rec][i]), rand_cat->zz(randoms[rec][i])};
      coord_part[1] = {halo_cat->xx(i), halo_cat->yy(i), halo_cat->zz(i)};
      vector<double> displ = {displ_cat[rec]->xx(i), displ_cat[rec]->yy(i), displ_cat[rec]->zz(i)};
      vector<vector<unsigned int>> inds(2, vector<unsigned int>(3));
      // selezione cella
      for (unsigned int&& j=0; j<3; j++) {
        inds[0][j] = (coord_part[0][j]-min[j])/step;
        inds[1][j] = (coord_part[1][j]-min[j])/step;
      }
      if(inds[0]!=inds[1]) {
        //selezioni possibili facce per intersezione con il displacement   
        for (unsigned int part=0; part<2; part++) {
	  			for (unsigned int j=0; j<2; j++) { // 2 facce
	    			for (unsigned int&& k=0; k<3; k++) {
	      // trovo l'intersezione
	      			double t = (coord_part[0][k]-((inds[part][k]+j)*step+min[k]))/displ[k];
	      			if (t>=0. && t<=1.) {
								vector<double> inters_coord = {coord_part[0][0]-(coord_part[0][0]-coord_part[1][0])*t,
								coord_part[0][1]-(coord_part[0][1]-coord_part[1][1])*t, coord_part[0][2]-(coord_part[0][2]-coord_part[1][2])*t};
								if ((k==0) ? (inters_coord[1]>(inds[part][1]*step+min[1]) && inters_coord[1] <= ((inds[part][1]+1)*step+min[1]) &&
										inters_coord[2]>(inds[part][2]*step+min[2]) && inters_coord[2] <= ((inds[part][2]+1)*step+min[2])) :
										((k==1) ? (inters_coord[0]>(inds[part][0]*step+min[0]) && inters_coord[0] <= ((inds[part][0]+1)*step+min[0]) &&
										inters_coord[2]>(inds[part][2]*step+min[2]) && inters_coord[2] <= ((inds[part][2]+1)*step+min[2])) :
										(inters_coord[0]>(inds[part][0]*step+min[0]) && inters_coord[0] <= ((inds[part][0]+1)*step+min[0]) &&
										inters_coord[1]>(inds[part][1]*step+min[1]) && inters_coord[1] <= ((inds[part][1]+1)*step+min[1])))) {
									if(part==0) Divergence[inds[part][0]][inds[part][1]][inds[part][2]] -= abs(displ[k])/step;
									if(part==1) Divergence[inds[part][0]][inds[part][1]][inds[part][2]] += abs(displ[k])/step;
								}
							}
	     				else continue;
            }
          }
        }
      }
    }
  }

  coutCBL << "Smoothing..." << endl;

  vector<vector<vector<double>>> Divergence_copy = Divergence;

  for (unsigned int i=0; i<nCells; i++) {
    for (unsigned int j=0; j<nCells; j++) {
      for (unsigned int k=0; k<nCells; k++) {
        Divergence_copy[i][k][j] = Divergence_copy[i][k][j]/n_rec;
      }
    }
  }

  for (unsigned int i=0; i<nCells; i++) {
    for (unsigned int j=0; j<nCells; j++) {
      for (unsigned int k=0; k<nCells; k++) {
        double sum=0.;
        double sum_div=0.;
        for (int ii=std::max(0,(int)i-1); ii<std::min((int)nCells,(int)i+2); ii++) {
          for (int jj=std::max(0,(int)j-1); jj<std::min((int)nCells,(int)j+2); jj++) {
            for (int kk=std::max(0,(int)k-1); kk<std::min((int)nCells,(int)k+2); kk++) {
              double ex = exp(-pow(Euclidean_distance(min[0]+step*(i+0.5), min[0]+step*(ii+0.5), min[1]+step*(j+0.5), min[1]+step*(jj+0.5), min[2]+step*(k+0.5), min[2]+step*(kk+0.5)),2)/(2*step*step));
              sum += ex;
              sum_div += Divergence_copy[ii][jj][kk]*ex;
            }
          }
        }
        Divergence[i][j][k] = sum_div/sum;
      }
    }
  }  
  
  string name; 
  name = dir_output + "Divergence";
  name.append("_");
  name.append(output);
  ofstream fout_divergence(name.c_str());

  fout_divergence.precision(7);
  vector<vector<double>> coord_locmin(6);
  unsigned int cont = 0;

  for (unsigned int&& i=0; i<nCells; i++) {
    for (unsigned int&& j=0; j<nCells; j++) {
      for (unsigned int&& k=0; k<nCells; k++) {
        fout_divergence << min[0]+step*(i+0.5) << " " << min[1]+step*(j+0.5) << " " << min[2]+step*(k+0.5) << " " << Divergence[i][j][k] << endl;
        if (Divergence[i][j][k] < 0. && i>0 && j>0 && k>0 && i<nCells-1 && j<nCells-1 && k<nCells-1) {
          bool control = false; 
          double sum = 0.;
          vector<double> num(3, 0.), den(3, 0.);
          vector<int> delta_cells(3, 0);
          cont++;
          for (unsigned int&& ii=i-1; ii<i+2; ii++) {
            delta_cells[0] = i-ii; //inverse
            for (unsigned int&& jj=j-1; jj<j+2; jj++) {
              delta_cells[1] = j-jj; //inverse
              for (unsigned int&& kk=k-1; kk<k+2; kk++) {
                delta_cells[2] = k-kk; //inverse
                if (!(ii==i && jj==j && kk==k)) {
                  sum +=  Divergence[ii][jj][kk];
                  for (unsigned int&& t=0; t<3; t++) {
                    num[t] += Divergence[ii][jj][kk]*(step/Euclidean_distance((double)i, (double)ii, (double)j, (double)jj, (double)k, (double)kk))*delta_cells[t];
                    den[t] += abs(Divergence[ii][jj][kk]);
                  }
                  if (Divergence[ii][jj][kk]<=Divergence[i][j][k] || Divergence[ii][jj][kk] > 0.) {
                    control = true;
                    break;
                  }
                } 
              }
              if (control==true) break;
            }
            if (control==true) break;
          }
          if (control==false) {
          //if (control==false) {
						coord_locmin[0].push_back(min[0]+step*(i+0.5)+num[0]/den[0]);
						coord_locmin[1].push_back(min[1]+step*(j+0.5)+num[1]/den[1]);
						coord_locmin[2].push_back(min[2]+step*(k+0.5)+num[2]/den[2]);
						coord_locmin[3].push_back(i);
						coord_locmin[4].push_back(j);
						coord_locmin[5].push_back(k);
          }
        }
      }
    }
  }
  
  fout_divergence.clear(); fout_divergence.close();

  coutCBL << "Number of cells with negative divergence: " << cont << endl << endl;

  if (coord_locmin[0].size()==0) 
    ErrorCBL("No local minima found!", "Catalogue", "VoidCatalogue.cpp");

  coutCBL << "Number of local minima: " << coord_locmin[0].size() << endl << endl;
   
  // -------------------------------------------- //
  // ---------------- Third Step ---------------- //
  // -------------------------------------------- //
  
  coutCBL << "* * * Identification of voids * * *" << endl << endl;

  vector<double> radius_void(coord_locmin[0].size(), 0.);
  vector<double> temp_min = {halo_cat->Min(Var::_X_)-0.5*step, halo_cat->Min(Var::_Y_)-0.5*step, halo_cat->Min(Var::_Z_)-0.5*step};
  vector<double> temp_max = {halo_cat->Max(Var::_X_)+0.5*step, halo_cat->Max(Var::_Y_)+0.5*step, halo_cat->Max(Var::_Z_)+0.5*step};
  for (unsigned int&& i=0; i<coord_locmin[0].size(); i++) {
    double sum = 0., prev_sum=0.;
    int incr = 1, contCells = 0;
    while (sum<=0.) {
      int ii_min=coord_locmin[3][i]-incr;
      int jj_min=coord_locmin[4][i]-incr;
      int kk_min=coord_locmin[5][i]-incr;
      int ii_max=coord_locmin[3][i]+incr;
      int jj_max=coord_locmin[4][i]+incr;
      int kk_max=coord_locmin[5][i]+incr;
      if (ii_min>=0 && jj_min>=0 && kk_min>=0 && ii_max<(int)nCells && jj_max<(int)nCells && kk_max<(int)nCells) {
        int ii = ii_min;
        while (ii<=ii_max) {
          int jj = jj_min;
          while (jj<=jj_max) {
            int kk = kk_min;
            while (kk<=kk_max) {
              double dist = Euclidean_distance((int)coord_locmin[3][i], (int)ii, (int)coord_locmin[4][i], (int)jj, (int)coord_locmin[5][i], (int)kk);        
              if (dist<incr) sum+=Divergence[ii][jj][kk];
							contCells++;
              if (ii>ii_min+1 && ii<ii_max-1 && jj>jj_min+1 && jj<jj_max-1 && kk==kk_min+1) kk = kk_max-1;
              else kk++;
            }
            jj++;
          }
          ii++;
        }
        if (sum/contCells>=0.) {
          radius_void[i]=interpolated(0., {prev_sum, sum}, {(incr-1.5)*step, (incr-0.5)*step}, "Linear");
          break;
        }
        prev_sum = sum;
      }
      else 
      {
      	vector<double> t = {coord_locmin[0][i]-temp_min[0], coord_locmin[1][i]-temp_min[1], coord_locmin[2][i]-temp_min[2], 
        		temp_max[0]-coord_locmin[0][i], temp_max[1]-coord_locmin[1][i], temp_max[2]-coord_locmin[2][i]};
        radius_void[i]=*min_element(t.begin(), t.end());
      	break;
      }
      incr++;
    }
  }

  name = dir_output+"Voids_";
  name.append(output);
  ofstream fout_voids(name.c_str());
  for(int&& i=0; i<(int)radius_void.size(); ++i)  fout_voids << coord_locmin[0][i] << " " << coord_locmin[1][i] << " " << coord_locmin[2][i] << " " << radius_void[i] << endl;
  fout_voids.clear(); fout_voids.close();
  
  cout << endl;
  coutCBL << "Number of voids: " << (int)radius_void.size() << endl << endl;


  // ---------------------------------------------------------- //
  // ---------------- Computing execution time ---------------- //
  // ---------------------------------------------------------- //
  
  double seconds = omp_get_wtime() - start_time;
  
  cout << endl;
  coutCBL << "Time spent to compute: " << seconds       << " seconds " << endl ;
  coutCBL << "Time spent to compute: " << seconds/60.   << " minutes " << endl ;
  coutCBL << "Time spent to compute: " << seconds/3600. << " hours "   << endl ;
  cout << endl;
}

///========================================================================================

cbl::catalogue::Catalogue::Catalogue (const VoidAlgorithm algorithm, Catalogue tracer_catalogue, Catalogue random_catalogue, const std::string dir_output, const std::string output, const double cellsize, const cbl::cosmology::Cosmology cosm, const std::vector<double> RA_range, const std::vector<double> DEC_range, const int n_rec, const double step_size, const double threshold)
{
  // -------------------------------------------- //
  // ---------------- First Step ---------------- //
  // -------------------------------------------- //
    
  const double start_time = omp_get_wtime();
  shared_ptr<Catalogue> halo_cat, rand_cat;
  vector<Catalogue> displacement_catalogue(n_rec);
  vector<vector<unsigned int>> randoms;
  vector<double> cells_edges;
  std::setprecision(5);
  vector<vector<double>> prop = tracer_catalogue.compute_catalogueProperties_lightCone(cosm, RA_range, DEC_range, 25);


  auto poly = *[](const vector<double> x, vector<double> &parameter) -> const vector<double> {

    vector<double> model(x.size());
    for (size_t i=0; i<x.size(); i++)
      model[i] = parameter[0]*x[i]*x[i]*x[i] + parameter[1]*x[i]*x[i] + parameter[2]*x[i] + parameter[3];
		
    return model;
  };

  const function<vector<double>(vector<double>,vector<double>&, double)> mps_func = [&poly](const vector<double> x, std::vector<double> &parameter, double stepsize) {
    vector<double> model(x.size());
    for (size_t i=0; i<x.size(); i++) 
      model[i] = stepsize*poly({x[i]}, parameter)[0];
    
    return model;
  };

  if (algorithm==VoidAlgorithm::_LaZeVo_) {

    Catalogue rand_catalogue; 
      
    long int time = (long int)chrono::steady_clock::now().time_since_epoch().count();
    //time = (time%1000000);
    time = 2343;
    int rndd = time; 

    double mps_trcat = tracer_catalogue.mps();
      
    vector<vector<double>> data_trcat(3);
    data_trcat[0] = prop[0];
    for (size_t i=0; i<data_trcat[0].size(); i++) data_trcat[0][i] = cosm.D_C(data_trcat[0][i]);
    data_trcat[1]=prop[4];
    data_trcat[2]=prop[6]; 

    string name =  "../output/Mps_";
    name.append(output);
    ofstream fout_displ(name.c_str());
    cout << endl;
    for (size_t i=0; i<data_trcat[0].size(); i++)
    {
      fout_displ.precision(5);
      fout_displ << data_trcat[0][i] << " " << data_trcat[1][i] << " " << data_trcat[2][i] << endl;   
    } 
    fout_displ.clear(); fout_displ.close(); 
       
    const cbl::data::Data1D data(data_trcat[0], data_trcat[1], data_trcat[2]);
    shared_ptr<cbl::data::Data> ptr_data = make_shared<cbl::data::Data1D>(data);
      
    const int nparameters = 4;
    const vector<string> parNames = {"A", "B", "C", "D"};
    double valA = 0., valB = 0., valC = 0., valD = data_trcat[1][0];     
    
    vector<cbl::statistics::ParameterType> parType(nparameters, cbl::statistics::ParameterType::_Base_);
    auto ptr_modelInput = make_shared<cbl::cosmology::Cosmology>(cosm);
  
    const statistics::Model1D model(poly, nparameters, parType, parNames, ptr_modelInput);
    auto ptr_model = make_shared<cbl::statistics::Model1D>(model);
    
    cbl::statistics::Likelihood likelihood(ptr_data, ptr_model, cbl::statistics::LikelihoodType::_Gaussian_Error_);
  
    const vector<double> starting_p = {valA, valB, valC, valD};
    double minA = -1., maxA = 1.;
    double minB = -1., maxB = 1.;
    double minC = -1., maxC = 1.;
    double minD = 0., maxD = data_trcat[1][0];
    const vector<vector<double>> limits = { {minA, maxA}, {minB, maxB}, {minC, maxC}, {minD, maxD} };

    likelihood.maximize(starting_p, limits);

    const double start=cosm.D_C(cbl::Min(tracer_catalogue.var(Var::_Redshift_)))-data_trcat[1][0]/2;
    const double end=cosm.D_C(cbl::Max(tracer_catalogue.var(Var::_Redshift_)))+data_trcat[1][data_trcat[1].size()-1]/2;
    double front=start;
    double toll=data_trcat[1][0]*0.001;
    cells_edges.push_back(start);

    while (front < end) {
      vector<double> par = likelihood.parameters()->bestfit_value();
      par[2] = par[2]-1;
      par[3] = par[3]+front;
      while (mps_func({front-(data_trcat[0][1]-data_trcat[0][0])/10.}, par, step_size)[0]*
             mps_func({front+(data_trcat[0][1]-data_trcat[0][0])/10.}, par, step_size)[0] > 0) front+=(data_trcat[0][1]-data_trcat[0][0])/10.;    
      double left = front-(data_trcat[0][1]-data_trcat[0][0])/10., right = front+(data_trcat[0][1]-data_trcat[0][0])/10.;
      while ( !(mps_func({front}, par, step_size)[0] > -toll && mps_func({front}, par, step_size)[0] < toll) ) {
        if (left<0.) {
          if (mps_func({front}, par, step_size)[0]<0.) left = front;
          else right = front;
          front = (right+left)/2;
        }
        else {
          if (mps_func({front}, par, step_size)[0]<0.) right = front;
          else left = front;
          front = (right+left)/2;
        }
      }
      cells_edges.push_back(front);
    }
    if (random_catalogue.nObjects()==0) 
    {
      rand_catalogue = Catalogue(RandomType::_createRandom_homogeneous_LC_, tracer_catalogue, 1., cosm, RA_range, DEC_range, 100, rndd);
      coutCBL << "... random catalogue created!" << endl;
    }
    else if (random_catalogue.nObjects() == tracer_catalogue.nObjects())
    {
      coutCBL << "Reading random catalogue of " << random_catalogue.nObjects() << " objects provides by user..." << endl;
      rand_catalogue = random_catalogue;
      coutCBL << "...Done!" << endl;
    }
    else {
      coutCBL << "Reading random catalogue of " << random_catalogue.nObjects() << " objects provides by user..." << endl;
      cout << endl;
      int diff = random_catalogue.nObjects() - tracer_catalogue.nObjects();
      if (random_catalogue.nObjects() > tracer_catalogue.nObjects()) coutCBL << "WARNING! Random catalogue has " << 
          diff << " more objects than the tracers catalogue. The number will be equalized by the function catalogue::equalize_random_LC." << endl;
      if (random_catalogue.nObjects() < tracer_catalogue.nObjects()) coutCBL << "WARNING! Random catalogue has " << 
          -diff << " less objects than the tracers catalogue. The number will be equalized by the function catalogue::equalize_random_LC." << endl;      
      rand_catalogue.equalize_random_lightCone(tracer_catalogue, cosm, RA_range, DEC_range, rndd);
      coutCBL << "...Done!" << endl;
    }

    // shared pointers to the tracer and the random catalogues
    rand_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(rand_catalogue)));
    halo_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(tracer_catalogue)));  

    // maximum for the number of objects in the tracer catalogue (unsigned int) 
    unsigned int num_objects = ((halo_cat->nObjects())%2==0) ? halo_cat->nObjects() : (halo_cat->nObjects())-1;

    if ((halo_cat->nObjects())%2 != 0) {
      rand_cat -> remove_object(num_objects);
      halo_cat -> remove_object(num_objects);
    }

    // vector for the indices of haloes
    std::vector<unsigned int> index_halo_cat(num_objects);
    std::iota (std::begin(index_halo_cat), std::end(index_halo_cat), 0); 
    cout << "num obj: " << num_objects << endl;
    cout << endl;
    coutCBL << "Creating ChainMesh... " << flush;
    // chainmesh setting
    CatalogueChainMesh ChainMesh_haloes = CatalogueChainMesh(cbl::Euclidean_distance, cellsize, halo_cat);
    CatalogueChainMesh ChainMesh_randoms = CatalogueChainMesh(cbl::Euclidean_distance, cellsize, rand_cat);
    unsigned int N_near_obj = 3*3*4*par::pi;
    vector<vector<unsigned int>> near_part = ChainMesh_haloes.N_nearest_objects_cat(N_near_obj); //the selected halo (that will have distance = 0) and the nearest 3))
    randoms.resize(n_rec, vector<unsigned int>(num_objects)); //haloes fixed
    for(int rec=0; rec<n_rec; rec++) {
      std::iota (std::begin(randoms[rec]), std::end(randoms[rec]), 0); 
      auto time2 = (long int)chrono::steady_clock::now().time_since_epoch().count();
      time2 = (time2%1000000);
      std::shuffle(randoms[rec].begin(), randoms[rec].end(), default_random_engine(time2));
    }

    cout << " Done! " << endl;
    cout << endl;
    coutCBL << "Setting starting configuration... " << endl;

    double dist = 4*mps_trcat;
    double minX = ChainMesh_haloes.lim()[0][0]+dist/10, maxX = ChainMesh_haloes.lim()[0][1]-dist/10;
    double minY = ChainMesh_haloes.lim()[1][0]+dist/10, maxY = ChainMesh_haloes.lim()[1][1]-dist/10;
    double minZ = ChainMesh_haloes.lim()[2][0]+dist/10, maxZ = ChainMesh_haloes.lim()[2][1]-dist/10;

   // FACCIO STAMPARE DISPL1 E CONFRONTO CON DISPL0  
    // RISOLVERE QUESTIONE COPIA CHAIN MESH E POI TESTARE (SEMBRA VADA MALE E LENTA) E RISOLVERE PROBLEMA IDENTIFICAZIONE VUOTI
    for (int rec=0; rec<n_rec; rec++) {    
        
      auto timeS = ((long int)chrono::steady_clock::now().time_since_epoch().count())%1000000;
      random::UniformRandomNumbers randX(minX, maxX, timeS);
      random::UniformRandomNumbers randY(minY, maxY, timeS+1);
      random::UniformRandomNumbers randZ(minZ, maxZ, timeS+2);

      CatalogueChainMesh ChM_halo_copy = ChainMesh_haloes;
      CatalogueChainMesh ChM_rand_copy = ChainMesh_randoms;
      cout << ChM_halo_copy.part(7).size() << " " << ChainMesh_haloes.part(7).size() << endl;
      cout << ChM_halo_copy.cellsize() << endl;
      unsigned int removed = 0;

      while (removed < num_objects) {
        vector<double> pos = {randX(), randY(), randZ()};
        vector<unsigned int> close=ChM_halo_copy.Close_objects(pos, dist);
        unsigned int tr_to_rmv = std::min(100, (int)close.size());
        if (tr_to_rmv > 0) {
          vector<unsigned int> close_random = ChM_rand_copy.N_nearest_objects(pos, close.size());
          std::random_device dev;
          std::mt19937 rng(dev());
          removed += tr_to_rmv;
          while (tr_to_rmv > 0) {
            std::uniform_int_distribution<std::mt19937::result_type> dist(0, tr_to_rmv-1);
            int rnd1 = dist(rng), rnd2 = dist(rng);
            randoms[rec][close[rnd1]] = close_random[rnd2];
            ChM_halo_copy.deletePart(close[rnd1]);
            ChM_rand_copy.deletePart(close_random[rnd2]);
            close.erase(close.begin()+rnd1);
            close_random.erase(close_random.begin()+rnd2);
            tr_to_rmv--;
          }
        }
      }
    }    
  
    
    coutCBL << "Performing the iterations... " << endl;
    cout << endl;

    for (int rec=0; rec<n_rec; rec++) {
      auto index_halo_cat_copy = index_halo_cat;
      double ratio = 1;
      unsigned int n_iter = 0;

      while (ratio > threshold) {

        index_halo_cat_copy = index_halo_cat; 
        auto time2 = (long int)chrono::steady_clock::now().time_since_epoch().count();
        time2 = (time2%1000000);
        std::shuffle(index_halo_cat_copy.begin(), index_halo_cat_copy.end(), default_random_engine(time2));
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> dist(1,N_near_obj-1);
        vector<bool> index_bool(num_objects, false), used(num_objects, false);
        //vector<vector<unsigned int>> H(halo_cat_len_4,vector<unsigned int>(4)), R_def(halo_cat_len_4,vector<unsigned int>(4));
    
#pragma omp parallel num_threads(omp_get_max_threads())
        {     
#pragma omp for schedule(static)
          for (unsigned int i=0; i<num_objects; i++) {
            vector<unsigned int> H(4), R_def(4), R(4), R_copy(4);
            unsigned int rand1 = dist(rng);
            unsigned int rand2 = rand1;
            unsigned int rand3 = rand1;
            while (rand1 == rand2 || rand1==rand3 || rand2 == rand3) {
              rand2 = dist(rng);
              rand3 = dist(rng);
            }
            //vector<unsigned int> R(4);
            if (used[near_part[index_halo_cat_copy[i]][0]] == false && used[near_part[index_halo_cat_copy[i]][rand1]] == false && 
            used[near_part[index_halo_cat_copy[i]][rand2]] == false && used[near_part[index_halo_cat_copy[i]][rand3]] == false) { 
            // selezione particelle coinvolte
              used[near_part[index_halo_cat_copy[i]][0]] = true;
              used[near_part[index_halo_cat_copy[i]][rand1]] = true;
              used[near_part[index_halo_cat_copy[i]][rand2]] = true;
              used[near_part[index_halo_cat_copy[i]][rand3]] = true;
              H[0]=near_part[index_halo_cat_copy[i]][0];
              H[1]=near_part[index_halo_cat_copy[i]][rand1];
              H[2]=near_part[index_halo_cat_copy[i]][rand2];
              H[3]=near_part[index_halo_cat_copy[i]][rand3];
              R[0]=randoms[rec][H[0]];
              R[1]=randoms[rec][H[1]];
              R[2]=randoms[rec][H[2]];
              R[3]=randoms[rec][H[3]];

              R_copy = R;
        // calcolo configurazione minima
              double dist_min = 1000000.;
              double dist;
              do {
                dist = 0.;
                for (size_t j=0; j<R.size(); j++) 
                  dist += cbl::Euclidean_distance(halo_cat->xx(H[j]), rand_cat->xx(R[j]), halo_cat->yy(H[j]), rand_cat->yy(R[j]), halo_cat->zz(H[j]), rand_cat->zz(R[j]));
                
                if (dist < dist_min) {
                  dist_min = dist;
                  R_def = R;
                } 
              } while(std::next_permutation(R.begin(), R.end()));
        
              //aggiornamento vettori haloes e randoms
              if(R_def != R_copy) index_bool[i] = true; 

              for (size_t j=0; j<4; j++) randoms[rec][H[j]] = R_def[j];

              used[H[0]] = false;
              used[H[1]] = false;
              used[H[2]] = false;
              used[H[3]] = false;
            }
            else continue;
          }
        }

        unsigned int changed_couples = count(index_bool.begin(), index_bool.end(), true);

        ratio = (double)changed_couples/num_objects;
        if (n_iter%100==0) cout << setprecision(5) << "ratio: " << ratio << " " << n_iter << endl;
        n_iter++;
      }

      double distanza = 0.;
      for (unsigned int i=0; i<num_objects; i++) {
      distanza += cbl::Euclidean_distance(halo_cat->xx(i), rand_cat->xx(randoms[rec][i]),
          halo_cat->yy(i), rand_cat->yy(randoms[rec][i]), halo_cat->zz(i), rand_cat->zz(randoms[rec][i]));
      }
      cout << "distanza finale: " << distanza/num_objects << endl;

      vector<shared_ptr<cbl::catalogue::Object>> displ_catalogue_object;

      for (size_t ii = 0; ii<randoms[rec].size(); ii++)
      { 
        cbl::comovingCoordinates displacement_values = {rand_cat->xx(randoms[rec][ii])-halo_cat->xx(ii), rand_cat->yy(randoms[rec][ii])-halo_cat->yy(ii), rand_cat->zz(randoms[rec][ii])-halo_cat->zz(ii)};
        auto halo_displ = make_shared<cbl::catalogue::Halo>(displacement_values);
        displ_catalogue_object.push_back(halo_displ); 
      }

////////////////////////////////////////////////////////////////////////////////////////////////////////
      displacement_catalogue[rec].add_objects(displ_catalogue_object);

      if (rec == 0) {
        stringstream ss;
        ss << rec;
        string str = ss.str();
        //       string output = ss.str();
        string name =  "../output/Displacement";
        name.append("_");
        name.append(str);
        name.append("_");
        //       output = "coord_LCDM.dat";
        name.append(output);
        ofstream fout_displ(name.c_str());
        cout << endl;

        coutCBL << "Writing the displacement ... " << endl;
        for (unsigned int i=0; i<num_objects; i++)
        {
          fout_displ.precision(7);
          fout_displ << halo_cat->xx(i) << " " << halo_cat->yy(i) << " " << halo_cat->zz(i) << " "
            << rand_cat->xx(randoms[rec][i]) << " " << rand_cat->yy(randoms[rec][i]) << " " << rand_cat->zz(randoms[rec][i]) << " "
            << displacement_catalogue[rec].xx(i) << " " << displacement_catalogue[rec].yy(i) << " " << displacement_catalogue[rec].zz(i) << endl;   
        } // modified
          
        fout_displ.clear(); fout_displ.close();  
      }
    }


  } // End of LaZeVo method
  else if (algorithm==VoidAlgorithm::_RIVA_) // Begin of RIVA method
  {
      
    cout << endl; coutCBL << par::col_green << "RIVA algorithm" << par::col_default << endl << endl;
    ErrorCBL ("Still not avaiable!", "Catalogue", "VoidCatalogue.cpp");
      /*
	(void)delta_movement;
	// * * * Two Point Correlation Function * * * //

	const double N_R = 1.; // random/data ratio

	const cbl::catalogue::Catalogue test_random_catalogue {cbl::catalogue::_box_, tracer_catalogue, N_R};

	const double rMin  = 1e3; // minimum separation 
	const double rMax  = 5e4; // maximum separation 
	const int nbins    = 50;  // number of bins
	const double shift = 0.5; // spatial shift used to set the bin centre

	string w_file = "xi_RIVA_pre.dat";
	cbl::twopt::TwoPointCorrelation1D_monopole TPCF_catalogue {tracer_catalogue, test_random_catalogue, cbl::_logarithmic_, rMin, rMax, nbins, shift};
	// // //     TPCF_catalogue.measure_for_RIVA(cbl::twopt::_Poisson_, dir_output);
	TPCF_catalogue.measure(cbl::twopt::_Poisson_, dir_output);

	TPCF_catalogue.write(dir_output, w_file);
    
    
	// // // // // WHILE CYCLE - - - while (2pcf_peak/max(xi[i])<1e-2) {...}
    
	auto halo_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(tracer_catalogue)));

	vector<shared_ptr<cbl::catalogue::Object>> final_catalogue_object; // vector shared pointer for the constructor 'final_catalogue'
    
	vector<double> x_position, y_position, z_position;

	for (size_t i = 0; i<num_objects; i++)
	{
	x_position.push_back(halo_cat->xx(i));
	y_position.push_back(halo_cat->yy(i));
	z_position.push_back(halo_cat->zz(i));
	}    

	srand(time(NULL));
	int randomDirection = 0;
    
	for (size_t i = 0; i<x_position.size(); i++)
	{
	randomDirection = (rand() % 3) + (-1);
	x_position[i]=x_position[i]+delta_movement*randomDirection;
	randomDirection = (rand() % 3) + (-1);
	y_position[i]=y_position[i]+delta_movement*randomDirection;
	randomDirection = (rand() % 3) + (-1);  
	z_position[i]=z_position[i]+delta_movement*randomDirection;
	}
    
	vector<shared_ptr<cbl::catalogue::Object>> temp_initial_catalogue_object;

	for (size_t i = 0; i<x_position.size(); i++)
	{
	cbl::comovingCoordinates initial_coordinates = {x_position[i], y_position[i], z_position[i]};
	auto inital_halo_coord = make_shared<cbl::catalogue::Halo>(initial_coordinates);
	temp_initial_catalogue_object.push_back(inital_halo_coord);
	}

	const cbl::catalogue::Catalogue test_random_catalogue_post {cbl::catalogue::_box_, tracer_catalogue, N_R};


	string w_file_post = "xi_RIVA_post.dat";
	cbl::twopt::TwoPointCorrelation1D_monopole TPCF_catalogue_post {temp_initial_catalogue_object, test_random_catalogue_post, cbl::_logarithmic_, rMin, rMax, nbins, shift};

	TPCF_catalogue_post.measure(cbl::twopt::_Poisson_, dir_output);

	TPCF_catalogue_post.write(dir_output, w_file_post);


	exit(1000);
      */
    } // End of RIVA method

  else ErrorCBL ("algorithm type is not correct!", "Catalogue", "VoidCatalogue.cpp");

  // --------------------------------------------- //
  // ---------------- Second Step ---------------- //
  // --------------------------------------------- //

  vector<shared_ptr<Catalogue>> displ_cat(n_rec);
  for (int i=0; i<n_rec; i++) displ_cat[i] = make_shared<Catalogue>(Catalogue(move(displacement_catalogue[i]))); 

  // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
  //   
  //   // * * * Two Point Correlation Function * * * //
  // 
  //   const double N_R = 1.; // random/data ratio
  //   
  //   const cbl::catalogue::Catalogue test_random_catalogue {cbl::catalogue::_box_, final_catalogue, N_R};
  // 
  //   const double rMin  = 1e3; // minimum separation 
  //   const double rMax  = 5e4; // maximum separation 
  //   const int nbins    = 20;  // number of bins
  //   const double shift = 0.5; // spatial shift used to set the bin centre
  //   
  //   string w_file = "xi.dat";
  //   
  //   cbl::twopt::TwoPointCorrelation1D_monopole TPCF_final_catalogue {final_catalogue, test_random_catalogue, cbl::_logarithmic_, rMin, rMax, nbins, shift};
  //       
  //   TPCF_final_catalogue.measure(cbl::twopt::_Poisson_, dir_output);
  //   
  //   TPCF_final_catalogue.write(dir_output, w_file);
  // 
  // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

  // ^ ^ ^ Uncomment the previous part to compute the Two Point Correlation Function of the final displacement with respect to a generic random field ^ ^ ^ //
  
  double THETA_min = -DEC_range[1]+par::pi/2, THETA_max = -DEC_range[0]+par::pi/2;
  double PHI_min = RA_range[0], PHI_max = RA_range[1];

  cout << endl;
  coutCBL << "RA range: [" << PHI_min << ", " << PHI_max << "]" << endl; 
  coutCBL << "Dec range: [" << DEC_range[0] << ", " << DEC_range[1] << "]" << endl; 
  cout << endl;

  double delta_THETA = THETA_max - THETA_min, delta_PHI = RA_range[1]-RA_range[0];
  vector<double> volume = prop[2];
  vector<double> density = prop[3]; 
  vector<double> step(cells_edges.size()-1), z_cells_edges(cells_edges.size()), z_cells_centres(step.size()), cells_centre(step.size()), THETA_step(step.size()),
                 real_delta_THETA(step.size()), real_THETA_min(step.size());
  vector<unsigned int> THETA_ncells(step.size());
  
  for (size_t i=0; i<step.size(); i++) step[i]=cells_edges[i+1]-cells_edges[i];
  for (size_t i=0; i<cells_edges.size(); i++) z_cells_edges[i] = cosm.Redshift(cells_edges[i]);
  for (size_t i=0; i<step.size(); i++) cells_centre[i]=(cells_edges[i+1]+cells_edges[i])/2;
  for (size_t i=0; i<step.size(); i++) z_cells_centres[i]=cosm.Redshift(cells_centre[i]);
  for (size_t i=0; i<step.size(); i++) {
    THETA_step[i] = step[i]/cells_centre[i];
    if ((((int)(delta_THETA/(step[i]/cells_centre[i]))+1)*THETA_step[i]) > par::pi) {
      step[i]=par::pi*cells_centre[i]/((int)(delta_THETA/(step[i]/cells_centre[i]))+1);
      THETA_step[i] = step[i]/cells_centre[i];
    }
  }
  for (size_t i=0; i<step.size(); i++) THETA_ncells[i]=delta_THETA/THETA_step[i]+1;
  for (size_t i=0; i<step.size(); i++) real_delta_THETA[i]=THETA_ncells[i]*THETA_step[i];
  for (size_t i=0; i<step.size(); i++) real_THETA_min[i]= (THETA_min - (real_delta_THETA[i]-delta_THETA)/2 < 0.) ? 0. : THETA_min - (real_delta_THETA[i]-delta_THETA)/2;

  vector<vector<vector<double>>> Divergence(step.size());
  vector<vector<double>> real_delta_PHI(step.size()), real_PHI_min(step.size()), PHI_step(step.size());

  for (size_t i=0; i<step.size(); i++) {
    Divergence[i].resize(THETA_ncells[i]);
    PHI_step[i].resize(THETA_ncells[i]);
    real_PHI_min[i].resize(THETA_ncells[i]);
    real_delta_PHI[i].resize(THETA_ncells[i]);
    for (unsigned int j=0; j<THETA_ncells[i]; j++) {
      double angle = real_THETA_min[i] + (j+0.5)*THETA_step[i];
      double temp_PHI_step = THETA_step[i]/sin(angle);
      PHI_step[i][j] = temp_PHI_step;
      if ((((int)(delta_PHI/PHI_step[i][j])+1)*PHI_step[i][j]) > 2*par::pi) PHI_step[i][j] = 2*par::pi/((int)(delta_PHI/PHI_step[i][j])+1);  
      int temp_nCells = delta_PHI/PHI_step[i][j] + 1;
      real_delta_PHI[i][j] = temp_nCells*PHI_step[i][j]; 
      real_PHI_min[i][j] = PHI_min - (real_delta_PHI[i][j]-delta_PHI)/2;
      Divergence[i][j].resize(temp_nCells);
    }
  }

  vector<double> cell_vol(step.size()), sph_area(step.size()+1), udrl_area(step.size()); 
  for (size_t i=0; i<step.size(); i++) cell_vol[i] = (-cos(par::pi/2+THETA_step[i]/2)+cos(par::pi/2-THETA_step[i]/2))*THETA_step[i]/(4*par::pi)*
                                                      (volume_sphere(cells_edges[i+1])-volume_sphere(cells_edges[i]));
  for (size_t i=0; i<=step.size()+1; i++) sph_area[i] = (-cos(par::pi/2+THETA_step[i]/2)+cos(par::pi/2-THETA_step[i]/2))*THETA_step[i]*cells_edges[i]*cells_edges[i];
  for (size_t i=0; i<step.size(); i++) udrl_area[i] = THETA_step[i]*(cells_edges[i+1]*cells_edges[i+1]-cells_edges[i]*cells_edges[i])/2;                                                                                                      
  
  unsigned int nCells = 0;
  for (size_t i=0; i<step.size(); i++)
    for (size_t j=0; j<Divergence[i].size(); j++) nCells+=Divergence[i][j].size();

  cout << endl;
  coutCBL << "Compute the Divergence field on " << nCells << " cells ..." << endl;
  cout << endl;

  for(int rec=0; rec<n_rec; rec++) {
    for (size_t i=0; i<randoms[rec].size(); i++) {
      vector<vector<double>> coord_part(2);
      coord_part[0] = {rand_cat->redshift(randoms[rec][i]), -(rand_cat->dec(randoms[rec][i]))+par::pi/2, 
                      (DEC_range[0] < 0. && rand_cat->ra(randoms[rec][i]) > par::pi) ? rand_cat->ra(randoms[rec][i])-2*par::pi : rand_cat->ra(randoms[rec][i]),
                       rand_cat->xx(randoms[rec][i]), rand_cat->yy(randoms[rec][i]), rand_cat->zz(randoms[rec][i])};
      coord_part[1] = {halo_cat->redshift(i), -(halo_cat->dec(i))+par::pi/2, 
                      (DEC_range[0] < 0. && halo_cat->ra(i) > par::pi) ? halo_cat->ra(i)-2*par::pi : halo_cat->ra(i), 
                      halo_cat->xx(i), halo_cat->yy(i), halo_cat->zz(i)};
      vector<double> displ = {displ_cat[rec]->xx(randoms[rec][i]), displ_cat[rec]->yy(randoms[rec][i]), displ_cat[rec]->zz(randoms[rec][i])};
      vector<vector<unsigned int>> inds(2, vector<unsigned int>(3));
      // selezione cella
      for (int p=0; p<2; p++) {
        int j = 0;
        while (coord_part[p][0] > z_cells_edges[j+1]) j++;
        inds[p][0] = j;
        inds[p][1] = (coord_part[p][1]-real_THETA_min[j])/THETA_step[j];
        inds[p][2] = (coord_part[p][2]-real_PHI_min[j][inds[p][1]])/PHI_step[j][inds[p][1]];
      }

      if(inds[0]!=inds[1]) {
        for (int p=0; p<2; p++) {
          // possibili intersezioni 
          // intersezione facce sferiche 
          for (int sph=0; sph<2; sph++) {
            vector<double> coord_intersect;
            double a, b, c;
            a = (coord_part[1][3]-coord_part[0][3])*(coord_part[1][3]-coord_part[0][3]) + 
                (coord_part[1][4]-coord_part[0][4])*(coord_part[1][4]-coord_part[0][4]) + 
                (coord_part[1][5]-coord_part[0][5])*(coord_part[1][5]-coord_part[0][5]);
            b = 2*coord_part[0][3]*(coord_part[1][3]-coord_part[0][3]) + 2*coord_part[0][4]*(coord_part[1][4]-coord_part[0][4]) + 
                2*coord_part[0][5]*(coord_part[1][5]-coord_part[0][5]);
            c = coord_part[0][3]*coord_part[0][3] + coord_part[0][4]*coord_part[0][4] + coord_part[0][5]*coord_part[0][5] -  
                cells_edges[inds[p][0]+sph]*cells_edges[inds[p][0]+sph];
            double del = b*b - 4*a*c; 
            if (del>0)  {
              double u1, u2;
              u1 = (-b+sqrt(del))/(2*a);
              u2 = (-b-sqrt(del))/(2*a);
              if( u1<=1. && u1>=0. ) coord_intersect = {coord_part[0][3]+(coord_part[1][3]-coord_part[0][3])*u1, 
                                                       coord_part[0][4]+(coord_part[1][4]-coord_part[0][4])*u1, 
                                                       coord_part[0][5]+(coord_part[1][5]-coord_part[0][5])*u1};
              else if ( u2<=1. && u2>=0. ) coord_intersect = {coord_part[0][3]+(coord_part[1][3]-coord_part[0][3])*u2, 
                                      coord_part[0][4]+(coord_part[1][4]-coord_part[0][4])*u2, 
                                      coord_part[0][5]+(coord_part[1][5]-coord_part[0][5])*u2};
              else continue;
              //calcolo ora il versore alla superficie
              double phi, theta;
              phi = (RA_range[0] > 0. && atan(coord_intersect[1]/coord_intersect[0]) < 0.) ? 
                    atan(coord_intersect[1]/coord_intersect[0]) + 2*par::pi : atan(coord_intersect[1]/coord_intersect[0]);
              theta = acos(coord_intersect[2]/cells_edges[inds[p][0]+sph]);
              if (theta>inds[p][1]*THETA_step[inds[p][0]]+real_THETA_min[inds[p][0]] && theta <= (inds[p][1]+1)*THETA_step[inds[p][0]]+real_THETA_min[inds[p][0]] &&
                 phi>inds[p][2]*PHI_step[inds[p][0]][inds[p][1]]+real_PHI_min[inds[p][0]][inds[p][1]] && phi<=(inds[p][2]+1)*PHI_step[inds[p][0]][inds[p][1]]+real_PHI_min[inds[p][0]][inds[p][1]]) {
                double normaliz = sqrt(coord_intersect[0]*coord_intersect[0]+coord_intersect[1]*coord_intersect[1]+coord_intersect[2]*coord_intersect[2]);
                vector<double> norm = {coord_intersect[0]/normaliz, coord_intersect[1]/normaliz, coord_intersect[2]/normaliz};
                // calcolo il prodotto scalare tra normale e displ
                if (sph == 0) {
                  norm[0] = -norm[0];
                  norm[1] = -norm[1];
                  norm[2] = -norm[2];
                }
                double scal = norm[0]*displ[0] + norm[1]*displ[1] + norm[2]*displ[2];
                //computo DIVERGENZA
                Divergence[inds[p][0]][inds[p][1]][inds[p][2]] += scal*sph_area[inds[p][0]+sph]/cell_vol[inds[p][0]];
              }
            }
          }
          
          // calcolo divergenza sulle 2 facce della cella sopra e sotto (intersezione cono/retta)
          for (int k=0; k<2; k++) {
            double a, b, c, theta = THETA_step[inds[p][0]]*(inds[p][1]+k)+real_THETA_min[inds[p][0]];
            double l = coord_part[1][3]-coord_part[0][3], m = coord_part[1][4]-coord_part[0][4], n = coord_part[1][5]-coord_part[0][5];
            
            a = (l*l + m*m - n*n*tan(theta)*tan(theta));
            b = 2*(l*coord_part[0][3] + m*coord_part[0][4] - n*coord_part[0][5]*tan(theta)*tan(theta));
            c = coord_part[0][3]*coord_part[0][3] + coord_part[0][4]*coord_part[0][4] - coord_part[0][5]*coord_part[0][5]*tan(theta)*tan(theta);
          
            double t1, t2, delta;
            delta = b*b - 4*a*c;
          
            if (delta > 0) {
              t1 = (-b+sqrt(delta))/(2*a);
              t2 = (-b-sqrt(delta))/(2*a);
            }
            else continue;
            vector<double> coord_intersect;
            if( t1<=1. && t1>=0. ) coord_intersect = {coord_part[0][3]+(coord_part[1][3]-coord_part[0][3])*t1, 
                                                       coord_part[0][4]+(coord_part[1][4]-coord_part[0][4])*t1, 
                                                       coord_part[0][5]+(coord_part[1][5]-coord_part[0][5])*t1};
            else if ( t2<=1. && t2>=0. ) coord_intersect = {coord_part[0][3]+(coord_part[1][3]-coord_part[0][3])*t2, 
                                      coord_part[0][4]+(coord_part[1][4]-coord_part[0][4])*t2, 
                                      coord_part[0][5]+(coord_part[1][5]-coord_part[0][5])*t2};
            else continue;

            double r = sqrt(coord_intersect[0]*coord_intersect[0]+coord_intersect[1]*coord_intersect[1]+coord_intersect[2]*coord_intersect[2]);
            double phi = (RA_range[0] > 0. && atan(coord_intersect[1]/coord_intersect[0]) < 0.) ? atan(coord_intersect[1]/coord_intersect[0]) + 2*par::pi : atan(coord_intersect[1]/coord_intersect[0]);
            if (r>cells_edges[inds[p][0]] && r<cells_edges[inds[p][0]+1] && 
                phi>inds[p][2]*PHI_step[inds[p][0]][inds[p][1]]+real_PHI_min[inds[p][0]][inds[p][1]] && phi<=(inds[p][2]+1)*PHI_step[inds[p][0]][inds[p][1]]+real_PHI_min[inds[p][0]][inds[p][1]]) {
              vector<double> norm;
              if (k==0) norm = {2*coord_intersect[0]*cos(theta)*cos(theta), 2*coord_intersect[1]*cos(theta)*cos(theta), -2*coord_intersect[0]*sin(theta)*sin(theta)}; 
              if (k==1) norm = {-2*coord_intersect[0]*cos(theta)*cos(theta), -2*coord_intersect[1]*cos(theta)*cos(theta), 2*coord_intersect[0]*sin(theta)*sin(theta)};             
              double normaliz = sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
              norm = {norm[0]/normaliz, norm[1]/normaliz, norm[2]/normaliz};
              double scal = norm[0]*displ[0] + norm[1]*displ[1] + norm[2]*displ[2];
              Divergence[inds[p][0]][inds[p][1]][inds[p][2]] += scal*udrl_area[inds[p][0]+k]/cell_vol[inds[p][0]];
            }
          }
          vector<double> temp_PHI(2);
          vector<vector<double>> cent_coord(2);
          for(int fac=0; fac<2; fac++) temp_PHI[fac] = PHI_step[inds[p][0]][inds[p][1]]*(inds[p][2]+fac) + real_PHI_min[inds[p][0]][inds[p][1]];
          
          cent_coord[0] = {par::pi/2, temp_PHI[0]+3*par::pi/2};
          cent_coord[1] = {par::pi/2, temp_PHI[1]+par::pi/2};
         
          for (int n=0; n<2; n++) {
            //eq piano per tre punti (di cui uno (0,0,0)) e coeff versore normale
            double a,b,c;

            a = sin(cent_coord[n][0])*cos(cent_coord[n][1]);
            b = sin(cent_coord[n][0])*sin(cent_coord[n][1]);
            c = cos(cent_coord[n][0]);

            // intersezione
            double l = coord_part[1][3]-coord_part[0][3], m = coord_part[1][4]-coord_part[0][4], n_ = coord_part[1][5]-coord_part[0][5];

            vector<double> coord_intersect(3);
            double t = (a*coord_part[1][3] + b*coord_part[1][4] + c*coord_part[1][5])/(a*l + b*m + c*n_); 
            if (t>0 && t<1) coord_intersect = {coord_part[1][3]-(coord_part[1][3]-coord_part[0][3])*t, 
                                         coord_part[1][4]-(coord_part[1][4]-coord_part[0][4])*t, 
                                         coord_part[1][5]-(coord_part[1][5]-coord_part[0][5])*t};
            else continue;
              
            double r, theta;
            r = sqrt(coord_intersect[0]*coord_intersect[0] + coord_intersect[1]*coord_intersect[1] + coord_intersect[2]*coord_intersect[2]);
            theta = acos(coord_intersect[2]/r);

            //verifica intersezione in cella. ordine: dx, up, sx, down

            if (r >= cells_edges[inds[p][0]] && r < cells_edges[inds[p][0]+1] && theta > THETA_step[inds[p][0]]*inds[p][1]+real_THETA_min[inds[p][0]] 
                  && theta < THETA_step[inds[p][0]]*(inds[p][1]+1)+real_THETA_min[inds[p][0]]) {
              // calcolo il prodotto scalare tra normale e displ
              double scal = a*displ[0] + b*displ[1] + c*displ[2];
              //computo DIVERGENZA
              Divergence[inds[p][0]][inds[p][1]][inds[p][2]] += scal*udrl_area[inds[p][0]]/cell_vol[inds[p][0]];
            }
          }
        }
      }
    }
  }

  Catalogue Divergence_cell;

  int n_cell = 0;
  for (size_t i=0; i<Divergence.size(); i++) {
    for (size_t j=0; j<Divergence[i].size(); j++) {
      for (size_t k=0; k<Divergence[i][j].size(); k++) {
        double radius = ((cells_edges[i+1]+cells_edges[i])/2);
        double theta = (j+0.5)*THETA_step[i] + real_THETA_min[i];
        double phi = (k+0.5)*PHI_step[i][j] + real_PHI_min[i][j];
        vector<double> cell_values = {radius*sin(theta)*cos(phi), radius*sin(theta)*sin(phi), radius*cos(theta)};
        shared_ptr<cbl::catalogue::Object> cell_obj = make_shared<cbl::catalogue::Halo>(cell_values[0], cell_values[1], cell_values[2],
                                                      phi, theta, z_cells_centres[i], Divergence[i][j][k]/n_rec);
        Divergence_cell.add_object(cell_obj);
        n_cell++;
      }
    }
  }

  shared_ptr<Catalogue> div_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(Divergence_cell)));
  CatalogueChainMesh ChainMesh_div = CatalogueChainMesh(cbl::Euclidean_distance, cellsize, div_cat);

  coutCBL << "Smoothing..." << endl;

  n_cell = 0;

  for (size_t i=0; i<Divergence.size(); i++) {
    for (size_t j=0; j<Divergence[i].size(); j++) {
      for (size_t k=0; k<Divergence[i][j].size(); k++) {
        if ((i>0 && i<Divergence.size()-1) && (j>0 && j<Divergence[i].size()-1) && (k>0 && k<Divergence[i][j].size()-1)) {
          vector<unsigned int> close_cells = ChainMesh_div.Close_objects({div_cat->xx(n_cell), div_cat->yy(n_cell), div_cat->zz(n_cell)}, 2*step[i]);
          double sum=0.;
          double sum_div=0.;
          for(auto&& t : close_cells) {
            double ex = exp(-pow(Euclidean_distance(div_cat->xx(n_cell),div_cat->xx(t), div_cat->yy(n_cell), div_cat->yy(t),
                                                    div_cat->zz(n_cell), div_cat->zz(t)),2)/(0.5*step[i]*step[i]));
            sum += ex;
            sum_div += div_cat->var(t, Var::_Weight_)*ex;
          }
          div_cat->set_var(n_cell, Var::_Weight_, sum_div/sum);
        }
        n_cell++;
      }
    }
  }  
  
  string name;  
  name = dir_output + "Divergence";
  name.append("_");
  name.append(output);
  ofstream fout_divergence(name.c_str());

  fout_divergence.precision(7);
  vector<vector<double>> coord_locmin(9);
  n_cell = 0;
  int cont = 0;

  for (size_t i=0; i<Divergence.size(); i++) {
    for (size_t j=0; j<Divergence[i].size(); j++) {
      for (size_t k=0; k<Divergence[i][j].size(); k++) {
        fout_divergence << div_cat->xx(n_cell) << " " << div_cat->yy(n_cell) << " " << div_cat->zz(n_cell) << " " << 
                           div_cat->ra(n_cell) << " " << div_cat->dec(n_cell) << " " << div_cat->redshift(n_cell) << " " <<
                           div_cat->weight(n_cell) << endl;

        if (div_cat->weight(n_cell)<0. && (i>0 && i<Divergence.size()-1) && (j>0 && j<Divergence[i].size()-1) && (k>0 && k<Divergence[i][j].size()-1)) {
          vector<unsigned int> close_cells = ChainMesh_div.Close_objects({div_cat->xx(n_cell), div_cat->yy(n_cell), div_cat->zz(n_cell)}, 1.5*step[i]);
          bool control = false; 
          double sum = 0., den = 0.;
          vector<double> num(3, 0.);
          for(auto&& t : close_cells) {
            if ((int)t != n_cell) {
              if (div_cat->weight(t) < div_cat->weight(n_cell)) {
                control = true;
                break;
              }
              sum += div_cat->weight(t);
              num[0] += div_cat->weight(t)*(step[i]/Euclidean_distance(div_cat->xx(n_cell),div_cat->xx(t), div_cat->yy(n_cell), div_cat->yy(t),
                                                    div_cat->zz(n_cell), div_cat->zz(t))*(div_cat->xx(n_cell)-div_cat->xx(t)));
              num[1] += div_cat->weight(t)*(step[i]/Euclidean_distance(div_cat->xx(n_cell),div_cat->xx(t), div_cat->yy(n_cell), div_cat->yy(t),
                                                    div_cat->zz(n_cell), div_cat->zz(t))*(div_cat->yy(n_cell)-div_cat->yy(t)));
              num[2] += div_cat->weight(t)*(step[i]/Euclidean_distance(div_cat->xx(n_cell),div_cat->xx(t), div_cat->yy(n_cell), div_cat->yy(t),
                                                    div_cat->zz(n_cell), div_cat->zz(t))*(div_cat->zz(n_cell)-div_cat->zz(t)));
              den += abs(div_cat->weight(t));
            }
          }
          if (control==false) {
            coord_locmin[0].push_back(div_cat->xx(n_cell)+num[0]/den);
            coord_locmin[1].push_back(div_cat->yy(n_cell)+num[1]/den);
            coord_locmin[2].push_back(div_cat->zz(n_cell)+num[2]/den);
            double r, theta, phi;
            r = sqrt(coord_locmin[0][cont]*coord_locmin[0][cont]+coord_locmin[1][cont]*coord_locmin[1][cont]+coord_locmin[2][cont]*coord_locmin[2][cont]);
            phi = atan2(coord_locmin[1][cont], coord_locmin[0][cont]);
            theta = acos(coord_locmin[2][cont]/r); 
            coord_locmin[3].push_back(phi);
            coord_locmin[4].push_back(theta);
            coord_locmin[5].push_back(cosm.Redshift(r));
            coord_locmin[6].push_back(i);
            coord_locmin[7].push_back(j);
            coord_locmin[8].push_back(k);
            cont++;
          }
        }
        n_cell++;
      }
    }
  }
  
  fout_divergence.clear(); fout_divergence.close();

  unsigned int num_voids = coord_locmin[0].size();
  //int number = x_void.size();
  if (num_voids==0) 
    ErrorCBL("No local minima found!", "Catalogue", "VoidCatalogue.cpp");

  coutCBL << "Number of local minima: " << num_voids << endl << endl;

  // -------------------------------------------- //
  // ---------------- Third Step ---------------- //
  // -------------------------------------------- //
  
  // Identification of voids //
  

  coutCBL << "* * * Identification of voids * * *" << endl << endl;

  vector<double> radius_void(num_voids, 0.);
  for (unsigned int&& i=0; i<num_voids; i++) {
    double sum = -1., prev_sum;
    double incr = 0.;
    vector<double> temp_vect = {cosm.D_C(z_cells_edges[z_cells_edges.size()-1]-z_cells_edges[coord_locmin[6][num_voids]]), 
                      cosm.D_C(z_cells_edges[coord_locmin[6][num_voids]]-z_cells_edges[0]), 
                      (coord_locmin[4][num_voids]-THETA_min)*cells_centre[coord_locmin[6][num_voids]], (THETA_max-coord_locmin[4][num_voids])*cells_centre[coord_locmin[6][num_voids]],
                      (coord_locmin[3][num_voids]-PHI_min)*cells_centre[coord_locmin[6][num_voids]]*sin(coord_locmin[3][num_voids]),
                      (PHI_max-coord_locmin[3][num_voids])*cells_centre[coord_locmin[6][num_voids]]*sin(coord_locmin[3][num_voids])};
    double dist_lim = *min_element(temp_vect.begin(), temp_vect.end());
    while (sum<=0.) {
      prev_sum = sum;
      sum = 0.;
      incr += 0.5;
      vector<unsigned int> close_cells = ChainMesh_div.Close_objects({coord_locmin[0][num_voids], coord_locmin[1][num_voids], 
                                                                     coord_locmin[2][num_voids]}, incr*step[coord_locmin[6][num_voids]]);
      
      for (auto&& j : close_cells) sum += div_cat->weight(j);
      if (sum>0.) radius_void[i] = interpolated(0., {prev_sum, sum}, 
          {(incr-0.5)*step[coord_locmin[6][num_voids]], incr*step[coord_locmin[6][num_voids]]}, "Linear");  
      if (sum <0. && incr*step[coord_locmin[6][num_voids]]>=dist_lim) radius_void[i] = dist_lim;
    }
  }

  name = dir_output+"Voids_";
  name.append(output);
  ofstream fout_voids(name.c_str());

  for(int&& i=0; i<(int)radius_void.size(); ++i)  fout_voids << coord_locmin[0][i] << " " << coord_locmin[1][i] << " " << coord_locmin[2][i] << " " <<
                                coord_locmin[3][i] << " " << coord_locmin[4][i] << " " << coord_locmin[5][i] << " " << radius_void[i] << endl;
  
  fout_voids.clear(); fout_voids.close();
  
  cout << endl;
  coutCBL << "Number of voids: " << (int)radius_void.size() << endl << endl;


  // ---------------------------------------------------------- //
  // ---------------- Computing execution time ---------------- //
  // ---------------------------------------------------------- //
  
  //float seconds = float( clock () - begin_time ) / CLOCKS_PER_SEC;
  double seconds = omp_get_wtime() - start_time;
  
  cout << endl;
  coutCBL << "Time spent to compute: " << seconds       << " seconds " << endl ;
  coutCBL << "Time spent to compute: " << seconds/60.   << " minutes " << endl ;
  coutCBL << "Time spent to compute: " << seconds/3600. << " hours "   << endl ;
  cout << endl;
  cout << dir_output << endl;
}

// ===================================================================================

void cbl::catalogue::Catalogue::clean_void_catalogue (const std::vector<bool> clean, const std::vector<double> delta_r, const double threshold, const double statistical_relevance, bool rescale, const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const double ratio, const bool checkoverlap, const Var ol_criterion)
{
  const double start_time = omp_get_wtime();
    
  // ---------------------------------------------------- //
  // ---------------- Cleaning Procedure ---------------- //
  // ---------------------------------------------------- //

  if (clean.size() != 3) ErrorCBL("wrong vector size!", "Catalogue", "VoidCatalogue.cpp");  
  if (clean[0] || clean[1] || clean[2]) {
    vector<int> counter(clean.size(), 0);
    vector<bool> remove(nObjects(), false);
    cout << endl;
    coutCBL << par::col_blue << " *** CLEANING PROCEDURE STARTED *** \n" << par::col_default << endl;
    double cleaning_time = omp_get_wtime();
    coutCBL << "Voids in the Catalogue: " << nObjects() << endl;
    if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
    for (size_t i=0; i<nObjects(); i++) {
      //remove object i from catalogue if not belonging to a given interval of radii:
      if (clean[0] && (radius(i) < delta_r[0] || delta_r[1] < radius(i))) {
        remove[i] = true;
	      counter[0]++;
      }
      //remove object i from catalogue if central density is too high:
      if (clean[1] && centralDensity(i) > threshold && !remove[i]) {
        remove[i] = true;
        counter[1]++;
      }
      //remove object i from catalogue if not statistically significant:
      if (clean[2] && densityContrast(i) < statistical_relevance && !remove[i]) {
        remove[i] = true;
        counter[2]++;
      }
    }
    
    remove_objects(remove);

    cout << endl;
    coutCBL << par::col_green << " --- Removing spurious voids --- \n" << par::col_default << endl;
    coutCBL << "Removed voids: " << endl;
    cout << "\t r_min - r_max criterion : " << counter[0] << "\n" <<
      "\t central density too high: " << counter[1] << "\n" <<
      "\t statistically irrelevant: " << counter[2] << "\n" <<
      "\t total removed: " << counter[0]+counter[1]+counter[2] << endl;
    coutCBL << "Time spent by the cleaning procedure: " << omp_get_wtime()-cleaning_time << " seconds \n" << endl;
  }
  
  coutCBL << "Voids in the Catalogue: " << nObjects() << endl;
  if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
  
  // ---------------------------------------------------- //
  // ----------------- Radius Rescaling ----------------- //
  // ---------------------------------------------------- //
  
  if (rescale) {

    double density = tracers_catalogue->numdensity();
   
    cout << endl;
    coutCBL << par::col_green << " --- Rescaling radii --- \n" << par::col_default << endl;

    double rescaling_time = omp_get_wtime();
    
		vector<bool> remove(nObjects(), false), bad_rescale(nObjects(), false);

    double min_X = tracers_catalogue->Min(Var::_X_), max_X = tracers_catalogue->Max(Var::_X_);
    double min_Y = tracers_catalogue->Min(Var::_Y_), max_Y = tracers_catalogue->Max(Var::_Y_);
    double min_Z = tracers_catalogue->Min(Var::_Z_), max_Z = tracers_catalogue->Max(Var::_Z_);

#pragma omp parallel num_threads(omp_get_max_threads())
    {     
#pragma omp for ordered schedule(dynamic)
      for (size_t j=0; j<nObjects(); j++) {
				vector<double> temp_val = {2*radius(j), xx(j)-min_X, max_X-xx(j), yy(j)-min_Y, max_Y-yy(j), zz(j)-min_Z, max_Z-zz(j), delta_r[1]};
				double value = *min_element(temp_val.begin(), temp_val.end());
#pragma omp ordered
				ChM.get_searching_region(value);
				vector<long> close = ChM.close_objects(coordinate(j));
				vector<double> radii;
				for (auto&& k : close) {
	 			 	double distance = Catalogue::distance(j, tracers_catalogue->catalogue_object(k));
	  			if (distance < value) radii.emplace_back(distance);
				}
				if (radii.size() < 3) {
					remove[j] = true;
	 				bad_rescale[j]=true;
	 			 	continue;
	 			}
	
				std::sort (radii.begin(), radii.end());
 	
				vector<vector<double>> data(2, vector<double>(radii.size())); //distances, density contrast 

				for (size_t i=0; i<radii.size()-1; i++) {
					data[0][i]=(radii[i]+radii[i+1])/2;
	      	data[1][i]=(i+1)/(volume_sphere(data[0][i])*density);
				}
				data[0][radii.size()-1] = 2*data[0][radii.size()-2] - data[0][radii.size()-3];
				data[1][radii.size()-1] = radii.size()/(volume_sphere(data[0][radii.size()-1])*density);

				int N = radii.size()-1;
				while (!(data[1][N-1] < threshold && data[1][N] > threshold) && N > 0) N--;

				double new_radius;
				if (N==0) {
						remove[j] = true;
						bad_rescale[j]=true;
						continue;
				}
				else 
					new_radius = interpolated(threshold, {data[1][N-1], data[1][N]}, {data[0][N-1], data[0][N]}, "Linear"); 
		
				if (new_radius > 0) set_var(j, Var::_Radius_, new_radius);
				else {
					remove[j] = true;
					bad_rescale[j]=true;
					continue;
				}           
#pragma omp critical
	coutCBL << "..." << int(double(j)/double(nObjects())*100) << "% completed\r"; cout.flush();
				
			}
		}

    cout << endl;

    remove_objects(remove);
    coutCBL << "Removed voids (bad rescaled): " << count(bad_rescale.begin(), bad_rescale.end(), true) << endl;
    
    vector<bool> remove_outofrange(nObjects(), false);
    for (size_t ii = 0; ii<nObjects(); ii++) 
      if (radius(ii) < delta_r[0] || delta_r[1] < radius(ii)) 
				remove_outofrange[ii] = true;
      
    remove_objects(remove_outofrange);
    cout << "\t Out of range ["+conv(delta_r[0],par::fDP2)+","+conv(delta_r[1],par::fDP2)+"] : " << count(remove_outofrange.begin(), remove_outofrange.end(), true) << endl;
    
    compute_centralDensity(tracers_catalogue, ChM);
    compute_densityContrast(tracers_catalogue, ChM, ratio);
	
    coutCBL << "Time spent by the rescaling procedure: " << omp_get_wtime()-rescaling_time << " seconds \n" << endl;
    
  }
  
  coutCBL << "Voids in the Catalogue: " << nObjects() << endl;
  if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
    
  // ---------------------------------------------------- //
  // ------------------ Overlap Check ------------------- //
  // ---------------------------------------------------- //
  
  if (checkoverlap) {
    cout << endl;
    coutCBL << par::col_green << " --- Checking for overlapping voids --- \n" << par::col_default << endl;

		double ol_time = omp_get_wtime();
    vector<bool> remove(nObjects(), false);

    if (ol_criterion == Var::_CentralDensity_) sort(ol_criterion, false);
    else if (ol_criterion == Var::_DensityContrast_) sort(ol_criterion, true);
    else ErrorCBL("allowed overlap criteria are '_CentralDensity_' or '_DensityContrast_' .", "Catalogue", "VoidCatalogue.cpp");

    coutCBL << "Generating chain-mesh for void centres ..." << endl;

    chainmesh::ChainMesh3D ChM_voids(Min(Var::_Radius_), var(Var::_X_), var(Var::_Y_), var(Var::_Z_), 2*Max(Var::_Radius_));
    
    for (size_t i = 0; i<nObjects(); i++) {
      if (!remove[i]) {
        vector<long> close = ChM_voids.close_objects(coordinate(i));
        std::sort(close.begin(), close.end());
        for (size_t&& j : close) {
          if(!remove[j]) {
            double distance = Catalogue::distance(i, catalogue_object(j));
            if (distance < radius(i)+radius(j) && i!=j) {
              if (ol_criterion == Var::_CentralDensity_) {
                if (centralDensity(i) < centralDensity(j)) {
                  remove[j] = true;
                }
                else if (centralDensity(i) > centralDensity(j)) {
                  remove[i] = true;
                  break;
                }
                else if (centralDensity(i) == centralDensity(j)) {
                  if (densityContrast(i) < densityContrast(j)) {
                    remove[i] = true;
                    break;
                  }
                  else 
                    remove[j] = true;
                }
              }
              else if (ol_criterion == Var::_DensityContrast_) {
                if (densityContrast(i) < densityContrast(j)) {
                  remove[i] = true;
                  break;
                }
                else if (densityContrast(i) > densityContrast(j)) 
                  remove[j] = true;
                else if (densityContrast(i) == densityContrast(j)) {
                  if (centralDensity(i) < centralDensity(j)) 
                    remove[j] = true;
                  else {
                    remove[i] = true;
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }
 
    remove_objects(remove);
 
    coutCBL << "Voids removed to avoid overlap: " << count(remove.begin(), remove.end(), true) << endl;
    coutCBL << "Time spent by the overlap-checking procedure: " << omp_get_wtime()-ol_time << " seconds" << endl;
  }
  
  cout << endl;
  coutCBL << "Voids in the Catalogue: " << nObjects() << endl;
  if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp");
  
  cout << endl;
  coutCBL << "Total time spent: " << omp_get_wtime()-start_time << " seconds \n" << endl;

}


//********************************** VERSION 2.0 **********************************//

cbl::catalogue::Catalogue::Catalogue (const std::shared_ptr<Catalogue> input_voidCatalogue, const double Volume, const std::vector<bool> clean, const std::vector<double> delta_r, const double threshold, const double statistical_relevance, bool rescale, const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const double ratio, const bool checkoverlap, const Var ol_criterion)
{

  auto catalogue = input_voidCatalogue;
  //clock_t begin_time = clock();
  const double start_time = omp_get_wtime();
    
  // ---------------------------------------------------- //
  // ---------------- Cleaning Procedure ---------------- //
  // ---------------------------------------------------- //

  if (clean.size() != 3) ErrorCBL("wrong vector size!", "Catalogue", "VoidCatalogue.cpp");  
  if (clean[0] || clean[1] || clean[2]) {
    vector<int> counter(clean.size(), 0);
    vector<bool> remove(catalogue->nObjects(), false);
    cout << endl;
    coutCBL << par::col_blue << " *** CLEANING PROCEDURE STARTED *** \n" << par::col_default << endl;
    double cleaning_time = omp_get_wtime();
    coutCBL << "Voids in the Catalogue: " << catalogue->nObjects() << endl;
    if (catalogue->nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp");
    
    for (size_t i=0; i<catalogue->nObjects(); i++) {
      
      //remove object i from catalogue if not belonging to a given interval of radii:
      if (clean[0] && 
	  (catalogue->radius(i) < delta_r[0] || 
	   delta_r[1] < catalogue->radius(i))) {
	remove[i] = true;
	counter[0]++;
      }
      
      //remove object i from catalogue if central density is too high:
      if (clean[1] &&
	  catalogue->centralDensity(i) > threshold &&
	  !remove[i]) {
	remove[i] = true;
	counter[1]++;
      }
      
      //remove object i from catalogue if not statistically significant:
      if (clean[2] &&
	  catalogue->densityContrast(i) < statistical_relevance &&
	  !remove[i]) {
	remove[i] = true;
	counter[2]++;
      }
    }
    
    catalogue->remove_objects(remove);

    cout << endl;
    coutCBL << par::col_green << " --- Removing spurious voids --- \n" << par::col_default << endl;
    coutCBL << "Removed voids: " << endl;
    cout << "\t r_min - r_max criterion : " << counter[0] << "\n" <<
      "\t central density too high: " << counter[1] << "\n" <<
      "\t statistically irrelevant: " << counter[2] << "\n" <<
      "\t total removed: " << counter[0]+counter[1]+counter[2] << endl;

    //float cleaning_time = float( clock() - start_time ) / CLOCKS_PER_SEC;
    coutCBL << "Time spent by the cleaning procedure: " << omp_get_wtime()-cleaning_time << " seconds \n" << endl;
  }
  
  cout << endl;
  coutCBL << "Voids in the Catalogue: " << catalogue->nObjects() << endl;
  if (catalogue->nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
  
  // ---------------------------------------------------- //
  // ----------------- Radius Rescaling ----------------- //
  // ---------------------------------------------------- //
  
  if (rescale) {
    
    const double density = tracers_catalogue->nObjects()/Volume;
    const double mps = pow(density, -1./3.);
    coutCBL << "Sample volume = " << Volume << " (Mpc/h)^3" << endl;
    coutCBL << "Sample density = " << density << " (Mpc/h)^-3" << endl;
    coutCBL << "Sample mps = " << mps << " Mpc/h" << endl;
    
    cout << endl;
    coutCBL << par::col_green << " --- Rescaling radii --- \n" << par::col_default << endl;

    double rescaling_time = omp_get_wtime();
    
    //vector to memorize which element of the catalogue has to be removed at the end of the procedure:
    vector<bool> remove(catalogue->nObjects(), false);
    
    //counter for regions without any tracer:
    int void_voids = 0;

    //counter for voids that the procedure can't clean properly:
    int bad_rescaled = 0;

    int num_obj = catalogue->nObjects();

#pragma omp parallel num_threads(omp_get_max_threads())
    {
#pragma omp for ordered schedule(dynamic)
      for (int j = 0; j<num_obj; j++) {
	
#pragma omp critical
	coutCBL << "..." << int(double(j)/double(num_obj)*100) << "% completed\r"; cout.flush();
	
	double value = (3.*catalogue->radius(j) < delta_r[1]) ? 3.*catalogue->radius(j) : delta_r[1];
#pragma omp ordered
	ChM.get_searching_region(value);
	vector<long> close = ChM.close_objects(catalogue->coordinate(j));
	
	//compute distances between the void and the surrounding particles
	vector<double> distances;

	for (auto&& k : close) {
	  double distance = catalogue->distance(j, tracers_catalogue->catalogue_object(k));
	  if (distance < value) distances.emplace_back(distance);
	}
	
	// find radius at which the required density (threshold) is reached
	if (distances.size() > 0) {

	  // find starting radius
	  std::sort (distances.begin(), distances.end());
	  vector<double>::iterator up = std::upper_bound(distances.begin(), distances.end(), catalogue->radius(j));
	  auto kk = std::distance(distances.begin(), up);

	  // shrink or expand to match the required threshold
	  if (kk/(volume_sphere(distances[kk-1])*density) > threshold)
	    while (kk/(volume_sphere(distances[kk-1])*density) > threshold && kk > 1) kk--; // either you shrink
	  else while (kk/(volume_sphere(distances[kk-1])*density) < threshold && kk < (int) distances.size()-1) kk++; // or you expand (check -1)
	  
	  // linear interpolation:
	  double new_radius = interpolated(threshold,
					   {kk/(volume_sphere(distances[kk-1])*density), (kk+1)/(volume_sphere(distances[kk])*density)},
					   {distances[kk-1], distances[kk]}, "Linear"); // gsl function

	  if ((kk/(volume_sphere(new_radius)*density))>(threshold+0.15) || (kk/(volume_sphere(new_radius)*density))<(threshold-0.15)) {
	    remove[j] = true;
#pragma omp critical
	    bad_rescaled++;
	  }
	  
	  if (new_radius > 0) catalogue->set_var(j, Var::_Radius_, fabs(new_radius));
	  else
	    {
	      remove[j] = true;
#pragma omp critical
	      bad_rescaled++;
	    }  
	}
      
	else {
	  remove[j] = true;
#pragma omp critical
	  void_voids++;
	}
	 
      }//for
    }
    
    catalogue->remove_objects(remove);
    coutCBL << "Removed voids:   " << endl;
    cout << "\t Empty voids removed: " << void_voids << endl;
    cout << "\t Bad rescaled voids removed: " << bad_rescaled  << endl;
    //if (negative > 0) WarningMsg("Warning: there were "+conv(negative,par::fINT)+" negative radii.");
   
    
    //vector to memorize which element of the catalogue has to be removed at the end of the procedure:
    vector<bool> remove_outofrange(catalogue->nObjects(), false);
    int outofrange = 0;
    for (size_t ii = 0; ii<catalogue->nObjects(); ii++) {
      if (catalogue->radius(ii) < delta_r[0] ||
	  delta_r[1] < catalogue->radius(ii)) {
	remove_outofrange[ii] = true;
	outofrange++;
      }
    }
    
    catalogue->remove_objects(remove_outofrange);
    cout << "\t Out of range ["+conv(delta_r[0],par::fDP2)+","+conv(delta_r[1],par::fDP2)+"] : " << outofrange << endl;
    
    //compute new central density and density contrast:
    catalogue->compute_centralDensity(tracers_catalogue, ChM);
    catalogue->compute_densityContrast(tracers_catalogue, ChM, ratio);

    //float rescaling_time = float( clock() - begin_time ) / CLOCKS_PER_SEC;
    //rescaling_time = rescaling_time - cleaning_time;
    coutCBL << "Time spent by the rescaling procedure: " << omp_get_wtime()-rescaling_time << " seconds \n" << endl;
    
  }//rescale part
  
  coutCBL << "Voids in the Catalogue: " << catalogue->nObjects() << endl;
  if (catalogue->nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
    
  // ---------------------------------------------------- //
  // ------------------ Overlap Check ------------------- //
  // ---------------------------------------------------- //
  
  if (checkoverlap) {
    cout << endl;
    coutCBL << par::col_green << " --- Checking for overlapping voids --- \n" << par::col_default << endl;

    double ol_time = omp_get_wtime();

    if (ol_criterion == Var::_CentralDensity_) catalogue->sort(ol_criterion, false);
    else if (ol_criterion == Var::_DensityContrast_) catalogue->sort(ol_criterion, true);
    else ErrorCBL("allowed overlap criteria are '_CentralDensity_' or '_DensityContrast_' .", "Catalogue", "VoidCatalogue.cpp");
       
    vector<bool> remove(catalogue->nObjects(), false);
    
    coutCBL << "Generating chain-mesh for void centres ..." << endl;
    chainmesh::ChainMesh3D ChM_voids(catalogue->Min(Var::_Radius_),
				     catalogue->var(Var::_X_),
				     catalogue->var(Var::_Y_),
				     catalogue->var(Var::_Z_),
				     2.*catalogue->Max(Var::_Radius_));
    for (size_t i = 0; i<catalogue->nObjects(); i++) {
      vector<long> close = ChM_voids.close_objects(catalogue->coordinate(i));
      unsigned int h = 0;
      while (!remove[i] && h<close.size()) {
	double distance = catalogue->distance(i, catalogue->catalogue_object(close[h]));
	if ((distance < catalogue->radius(i)+catalogue->radius(close[h]) && (int) i < close[h]) ||
	    ((int) i > close[h] && remove[close[h]])) h++;
	else if (distance < catalogue->radius(i)+catalogue->radius(close[h]) && (int) i > close[h]) remove[i] = true;
	else h++;
      }//while
    }//for
    int overlap_removed = 0;
    for (size_t i = 0; i<remove.size(); i++) if (remove[i]) overlap_removed++;
    catalogue->remove_objects(remove);
    coutCBL << "Voids removed to avoid overlap: " << overlap_removed << endl;
    coutCBL << "Time spent by the overlap-checking procedure: " << omp_get_wtime()-ol_time << " seconds" << endl;
  }//overlap check
  cout << endl;
  coutCBL << "Voids in the Catalogue: " << catalogue->nObjects() << endl;
  if (catalogue->nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp");
  
  //float olchecking_time = float( clock() - begin_time ) / CLOCKS_PER_SEC;
  //olchecking_time = olchecking_time - rescaling_time - cleaning_time;
  cout << endl;
  coutCBL << "Total time spent: " << omp_get_wtime()-start_time << " seconds \n" << endl;

  m_object = catalogue->sample();
}


//********************************** VERSION 3.0 **********************************//

cbl::catalogue::Catalogue::Catalogue (const std::shared_ptr<Catalogue> input_voidCatalogue, const std::vector<double> par_numdensity, const std::vector<bool> clean, const std::vector<double> delta_r, const double threshold, const double statistical_relevance, bool rescale, const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const double ratio, const bool checkoverlap, const Var ol_criterion)
{

  auto catalogue = input_voidCatalogue;
  //clock_t begin_time = clock();
  const double start_time = omp_get_wtime();
    
  // ---------------------------------------------------- //
  // ---------------- Cleaning Procedure ---------------- //
  // ---------------------------------------------------- //

  if (clean.size() != 3) ErrorCBL("wrong vector size!", "Catalogue", "VoidCatalogue.cpp");  
  if (clean[0] || clean[1] || clean[2]) {
    vector<int> counter(clean.size(), 0);
    vector<bool> remove(catalogue->nObjects(), false);
    cout << endl;
    coutCBL << par::col_blue << " *** CLEANING PROCEDURE STARTED *** \n" << par::col_default << endl;
    double cleaning_time = omp_get_wtime();
    coutCBL << "Voids in the Catalogue: " << catalogue->nObjects() << endl;
    if (catalogue->nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
    for (size_t i=0; i<catalogue->nObjects(); i++) {
      
      //remove object i from catalogue if not belonging to a given interval of radii:
      if (clean[0] && 
	  (catalogue->radius(i) < delta_r[0] || 
	   delta_r[1] < catalogue->radius(i))) {
	remove[i] = true;
	counter[0]++;
      }
      
      //remove object i from catalogue if central density is too high:
      if (clean[1] &&
	  catalogue->centralDensity(i) > threshold &&
	  !remove[i]) {
	remove[i] = true;
	counter[1]++;
      }
      
      //remove object i from catalogue if not statistically significant:
      if (clean[2] &&
	  catalogue->densityContrast(i) < statistical_relevance &&
	  !remove[i]) {
	remove[i] = true;
	counter[2]++;
      }
    }
    
    catalogue->remove_objects(remove);

    cout << endl;
    coutCBL << par::col_green << " --- Removing spurious voids --- \n" << par::col_default << endl;
    coutCBL << "Removed voids: " << endl;
    cout << "\t r_min - r_max criterion : " << counter[0] << "\n" <<
      "\t central density too high: " << counter[1] << "\n" <<
      "\t statistically irrelevant: " << counter[2] << "\n" <<
      "\t total removed: " << counter[0]+counter[1]+counter[2] << endl;

    //float cleaning_time = float( clock() - start_time ) / CLOCKS_PER_SEC;
    coutCBL << "Time spent by the cleaning procedure: " << omp_get_wtime()-cleaning_time << " seconds \n" << endl;
  }
  
  coutCBL << "Voids in the Catalogue: " << catalogue->nObjects() << endl;
  if (catalogue->nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
  
  // ---------------------------------------------------- //
  // ----------------- Radius Rescaling ----------------- //
  // ---------------------------------------------------- //
  
  if (rescale) {
    
    cout << endl;
    coutCBL << par::col_green << " --- Rescaling radii --- \n" << par::col_default << endl;

    double rescaling_time = omp_get_wtime();
    
    //vector to memorize which element of the catalogue has to be removed at the end of the procedure:
    vector<bool> remove(catalogue->nObjects(), false);
    
    //counter for regions without any tracer:
    int void_voids = 0;

    //counter for voids that the procedure can't clean properly:
    int bad_rescaled = 0;

    int num_obj = catalogue->nObjects();

#pragma omp parallel num_threads(omp_get_max_threads())
    {
#pragma omp for ordered schedule(dynamic)
      for (int j = 0; j<num_obj; j++) {
	
#pragma omp critical
	coutCBL << "..." << int(double(j)/double(num_obj)*100) << "% completed\r"; cout.flush();
	
	double value = (3.*catalogue->radius(j) < delta_r[1]) ? 3.*catalogue->radius(j) : delta_r[1];
#pragma omp ordered
	ChM.get_searching_region(value);
	vector<long> close = ChM.close_objects(catalogue->coordinate(j));
	
	//compute distances between the void and the surrounding particles
	vector<double> distances;

	for (auto&& k : close) {
	  double distance = catalogue->distance(j, tracers_catalogue->catalogue_object(k));
	  if (distance < value) distances.emplace_back(distance);
	}	
	
	// find radius at which the required density (threshold) is reached
	if (distances.size() > 0) {

	  // redshift if the centre of the considered void
	  double zz = catalogue->redshift(j); 
	  
	  // compute the number density using a Nth-order polynomial
	  double density = 0.;
	  for (size_t N=par_numdensity.size(); N-->0;) density += par_numdensity[par_numdensity.size()-1-N]*pow(zz,N);
	  
	  // find starting radius
	  std::sort (distances.begin(), distances.end());
	  vector<double>::iterator up = std::upper_bound(distances.begin(), distances.end(), catalogue->radius(j));
	  auto kk = std::distance(distances.begin(), up);

	  // shrink or expand to match the required threshold
	  if (kk/(volume_sphere(distances[kk-1])*density) > threshold)
	    while (kk/(volume_sphere(distances[kk-1])*density) > threshold && kk > 1) kk--; // either you shrink
	  else while (kk/(volume_sphere(distances[kk-1])*density) < threshold && kk < (int) distances.size()-1) kk++; // or you expand (check -1)
	  
	
	  // linear interpolation:
	  double new_radius = interpolated(threshold,
					   {kk/(volume_sphere(distances[kk-1])*density), (kk+1)/(volume_sphere(distances[kk])*density)},
					   {distances[kk-1], distances[kk]}, "Linear"); // gsl function

	  if ((kk/(volume_sphere(new_radius)*density))>(threshold+0.15) || (kk/(volume_sphere(new_radius)*density))<(threshold-0.15)) {
	    remove[j] = true;
#pragma omp critical
	    bad_rescaled++;
	  }
	  
	  if (new_radius > 0) catalogue->set_var(j, Var::_Radius_, fabs(new_radius));
	  else
	    {
	      remove[j] = true;
#pragma omp critical
	      bad_rescaled++;
	    }  
	}
      
	else {
	  remove[j] = true;
#pragma omp critical
	  void_voids++;
	}
	 
      }//for
    }
    
    catalogue->remove_objects(remove);
    coutCBL << "Removed voids:   " << endl;
    cout << "\t Empty voids removed: " << void_voids << endl;
    cout << "\t Bad rescaled voids removed: " << bad_rescaled  << endl;
    //if (negative > 0) WarningMsg("Warning: there were "+conv(negative,par::fINT)+" negative radii.");
   
    
    //vector to memorize which element of the catalogue has to be removed at the end of the procedure:
    vector<bool> remove_outofrange(catalogue->nObjects(), false);
    int outofrange = 0;
    for (size_t ii = 0; ii<catalogue->nObjects(); ii++) {
      if (catalogue->radius(ii) < delta_r[0] ||
	  delta_r[1] < catalogue->radius(ii)) {
	remove_outofrange[ii] = true;
	outofrange++;
      }
    }
    
    catalogue->remove_objects(remove_outofrange);
    cout << "\t Out of range ["+conv(delta_r[0],par::fDP2)+","+conv(delta_r[1],par::fDP2)+"] : " << outofrange << endl;
    
    //compute new central density and density contrast:
    catalogue->compute_centralDensity(tracers_catalogue, ChM, par_numdensity, ratio);
    catalogue->compute_densityContrast(tracers_catalogue, ChM, ratio);

    //float rescaling_time = float( clock() - begin_time ) / CLOCKS_PER_SEC;
    //rescaling_time = rescaling_time - cleaning_time;
    coutCBL << "Time spent by the rescaling procedure: " << omp_get_wtime()-rescaling_time << " seconds \n" << endl;
    
  }//rescale part
  
  coutCBL << "Voids in the Catalogue: " << catalogue->nObjects() << endl;
  if (catalogue->nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp");  

    
  // ---------------------------------------------------- //
  // ------------------ Overlap Check ------------------- //
  // ---------------------------------------------------- //
  
  if (checkoverlap) {
    cout << endl;
    coutCBL << par::col_green << " --- Checking for overlapping voids --- \n" << par::col_default << endl;

    double ol_time = omp_get_wtime();

    if (ol_criterion == Var::_CentralDensity_) catalogue->sort(ol_criterion, false);
    else if (ol_criterion == Var::_DensityContrast_) catalogue->sort(ol_criterion, true);
    else ErrorCBL("allowed overlap criteria are '_CentralDensity_' or '_DensityContrast_' .", "Catalogue", "VoidCatalogue.cpp");
       
    vector<bool> remove(catalogue->nObjects(), false);
    
    coutCBL << "Generating chain-mesh for void centres ..." << endl;
    chainmesh::ChainMesh3D ChM_voids(catalogue->Min(Var::_Radius_),
				     catalogue->var(Var::_X_),
				     catalogue->var(Var::_Y_),
				     catalogue->var(Var::_Z_),
				     2.*catalogue->Max(Var::_Radius_));
    for (size_t i = 0; i<catalogue->nObjects(); i++) {
      vector<long> close = ChM_voids.close_objects(catalogue->coordinate(i));
      unsigned int h = 0;
      while (!remove[i] && h<close.size()) {
	double distance = catalogue->distance(i, catalogue->catalogue_object(close[h]));
	if ((distance < catalogue->radius(i)+catalogue->radius(close[h]) && (int) i < close[h]) ||
	    ((int) i > close[h] && remove[close[h]])) h++;
	else if (distance < catalogue->radius(i)+catalogue->radius(close[h]) && (int) i > close[h]) remove[i] = true;
	else h++;
      }//while
    }//for
    int overlap_removed = 0;
    for (size_t i = 0; i<remove.size(); i++) if (remove[i]) overlap_removed++;
    catalogue->remove_objects(remove);
    coutCBL << "Voids removed to avoid overlap: " << overlap_removed << endl;
    coutCBL << "Time spent by the overlap-checking procedure: " << omp_get_wtime()-ol_time << " seconds" << endl;
  }//overlap check
  cout << endl;
  coutCBL << "Voids in the Catalogue: " << catalogue->nObjects() << endl;
  if (catalogue->nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 

  //float olchecking_time = float( clock() - begin_time ) / CLOCKS_PER_SEC;
  //olchecking_time = olchecking_time - rescaling_time - cleaning_time;
  cout << endl;
  coutCBL << "Total time spent: " << omp_get_wtime()-start_time << " seconds \n" << endl;

  m_object = catalogue->sample();
}


//********************************** VERSION 4.0 **********************************//

cbl::catalogue::Catalogue::Catalogue (const std::shared_ptr<Catalogue> input_voidCatalogue, const std::vector<std::vector<double>> data_numdensity, const std::string method_interpolation, const std::vector<bool> clean, const std::vector<double> delta_r, const double threshold, const double statistical_relevance, bool rescale, const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const double ratio, const bool checkoverlap, const Var ol_criterion)
{

  auto catalogue = input_voidCatalogue;
  //clock_t begin_time = clock();
  const double start_time = omp_get_wtime();
    
  // ---------------------------------------------------- //
  // ---------------- Cleaning Procedure ---------------- //
  // ---------------------------------------------------- //

  if (clean.size() != 3) ErrorCBL("wrong vector size!", "Catalogue", "VoidCatalogue.cpp");
  if (data_numdensity.size() != 2) ErrorCBL("data_numdensity must have 2 rows: redshift values in data_numdensity[0] and number density values in data_numdensity[1]", "Catalogue", "VoidCatalogue.cpp");
  if (clean[0] || clean[1] || clean[2]) {
    vector<int> counter(clean.size(), 0);
    vector<bool> remove(catalogue->nObjects(), false);
    cout << endl;
    coutCBL << par::col_blue << " *** CLEANING PROCEDURE STARTED *** \n" << par::col_default << endl;
    double cleaning_time = omp_get_wtime();
    coutCBL << "Voids in the Catalogue: " << catalogue->nObjects() << endl;
    if (catalogue->nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
    for (size_t i=0; i<catalogue->nObjects(); i++) {
      
      //remove object i from catalogue if not belonging to a given interval of radii:
      if (clean[0] && 
	  (catalogue->radius(i) < delta_r[0] || 
	   delta_r[1] < catalogue->radius(i))) {
	remove[i] = true;
	counter[0]++;
      }
      
      //remove object i from catalogue if central density is too high:
      if (clean[1] &&
	  catalogue->centralDensity(i) > threshold &&
	  !remove[i]) {
	remove[i] = true;
	counter[1]++;
      }
      
      //remove object i from catalogue if not statistically significant:
      if (clean[2] &&
	  catalogue->densityContrast(i) < statistical_relevance &&
	  !remove[i]) {
	remove[i] = true;
	counter[2]++;
      }
    }
    
    catalogue->remove_objects(remove);

    cout << endl;
    coutCBL << par::col_green << " --- Removing spurious voids --- \n" << par::col_default << endl;
    coutCBL << "Removed voids: " << endl;
    cout << "\t r_min - r_max criterion : " << counter[0] << "\n" <<
      "\t central density too high: " << counter[1] << "\n" <<
      "\t statistically irrelevant: " << counter[2] << "\n" <<
      "\t total removed: " << counter[0]+counter[1]+counter[2] << endl;

    //float cleaning_time = float( clock() - start_time ) / CLOCKS_PER_SEC;
    coutCBL << "Time spent by the cleaning procedure: " << omp_get_wtime()-cleaning_time << " seconds \n" << endl;
  }
  
  coutCBL << "Voids in the Catalogue: " << catalogue->nObjects() << endl;
  if (catalogue->nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
  
  // ---------------------------------------------------- //
  // ----------------- Radius Rescaling ----------------- //
  // ---------------------------------------------------- //
  
  if (rescale) {
    
    cout << endl;
    coutCBL << par::col_green << " --- Rescaling radii --- \n" << par::col_default << endl;

    double rescaling_time = omp_get_wtime();
    
    //vector to memorize which element of the catalogue has to be removed at the end of the procedure:
    vector<bool> remove(catalogue->nObjects(), false);
    
    //counter for regions without any tracer:
    int void_voids = 0;

    //counter for voids that the procedure can't clean properly:
    int bad_rescaled = 0;

    int num_obj = catalogue->nObjects();

#pragma omp parallel num_threads(omp_get_max_threads())
    {
#pragma omp for ordered schedule(dynamic)
      for (int j = 0; j<num_obj; j++) {
	
#pragma omp critical
	coutCBL << "..." << int(double(j)/double(num_obj)*100) << "% completed\r"; cout.flush();
	
	double value = (3.*catalogue->radius(j) < delta_r[1]) ? 3.*catalogue->radius(j) : delta_r[1];
#pragma omp ordered
	ChM.get_searching_region(value);
	vector<long> close = ChM.close_objects(catalogue->coordinate(j));
	
	//compute distances between the void and the surrounding particles
	vector<double> distances;

	for (auto&& k : close) {
	  double distance = catalogue->distance(j, tracers_catalogue->catalogue_object(k));
	  if (distance < value) distances.emplace_back(distance);
	}	
	
	// find radius at which the required density (threshold) is reached
	if (distances.size() > 0) {

	  // redshift if the centre of the considered void
	  double zz = catalogue->redshift(j); 
	  
	  // compute the number density using a specified method
	  double density = interpolated(zz, data_numdensity[0], data_numdensity[1], method_interpolation);
	  
	  // find starting radius
	  std::sort (distances.begin(), distances.end());
	  vector<double>::iterator up = std::upper_bound(distances.begin(), distances.end(), catalogue->radius(j));
	  auto kk = std::distance(distances.begin(), up);

	  // shrink or expand to match the required threshold
	  if (kk/(volume_sphere(distances[kk-1])*density) > threshold)
	    while (kk/(volume_sphere(distances[kk-1])*density) > threshold && kk > 1) kk--; // either you shrink
	  else while (kk/(volume_sphere(distances[kk-1])*density) < threshold && kk < (int) distances.size()-1) kk++; // or you expand (check -1)
	  
	
	  // linear interpolation:
	  double new_radius = interpolated(threshold,
					   {kk/(volume_sphere(distances[kk-1])*density), (kk+1)/(volume_sphere(distances[kk])*density)},
					   {distances[kk-1], distances[kk]}, "Linear"); // gsl function

	  if ((kk/(volume_sphere(new_radius)*density))>(threshold+0.15) || (kk/(volume_sphere(new_radius)*density))<(threshold-0.15)) {
	    remove[j] = true;
#pragma omp critical
	    bad_rescaled++;
	  }
	  
	  if (new_radius > 0) catalogue->set_var(j, Var::_Radius_, fabs(new_radius));
	  else
	    {
	      remove[j] = true;
#pragma omp critical
	      bad_rescaled++;
	    }  
	}
      
	else {
	  remove[j] = true;
#pragma omp critical
	  void_voids++;
	}
	 
      }//for
    }
    
    catalogue->remove_objects(remove);
    coutCBL << "Removed voids:   " << endl;
    cout << "\t Empty voids removed: " << void_voids << endl;
    cout << "\t Bad rescaled voids removed: " << bad_rescaled  << endl;
    //if (negative > 0) WarningMsg("Warning: there were "+conv(negative,par::fINT)+" negative radii.");
   
    
    //vector to memorize which element of the catalogue has to be removed at the end of the procedure:
    vector<bool> remove_outofrange(catalogue->nObjects(), false);
    int outofrange = 0;
    for (size_t ii = 0; ii<catalogue->nObjects(); ii++) {
      if (catalogue->radius(ii) < delta_r[0] ||
	  delta_r[1] < catalogue->radius(ii)) {
	remove_outofrange[ii] = true;
	outofrange++;
      }
    }
    
    catalogue->remove_objects(remove_outofrange);
    cout << "\t Out of range ["+conv(delta_r[0],par::fDP2)+","+conv(delta_r[1],par::fDP2)+"] : " << outofrange << endl;
    
    //compute new central density and density contrast:
    catalogue->compute_centralDensity(tracers_catalogue, ChM, data_numdensity, method_interpolation, ratio);
    catalogue->compute_densityContrast(tracers_catalogue, ChM, ratio);

    //float rescaling_time = float( clock() - begin_time ) / CLOCKS_PER_SEC;
    //rescaling_time = rescaling_time - cleaning_time;
    coutCBL << "Time spent by the rescaling procedure: " << omp_get_wtime()-rescaling_time << " seconds \n" << endl;
    
  }//rescale part
  
  coutCBL << "Voids in the Catalogue: " << catalogue->nObjects() << endl;
  if (catalogue->nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp");  

    
  // ---------------------------------------------------- //
  // ------------------ Overlap Check ------------------- //
  // ---------------------------------------------------- //
  
  if (checkoverlap) {
    cout << endl;
    coutCBL << par::col_green << " --- Checking for overlapping voids --- \n" << par::col_default << endl;

    double ol_time = omp_get_wtime();

    if (ol_criterion == Var::_CentralDensity_) catalogue->sort(ol_criterion, false);
    else if (ol_criterion == Var::_DensityContrast_) catalogue->sort(ol_criterion, true);
    else ErrorCBL("allowed overlap criteria are '_CentralDensity_' or '_DensityContrast_' .", "Catalogue", "VoidCatalogue.cpp");
       
    vector<bool> remove(catalogue->nObjects(), false);
    
    coutCBL << "Generating chain-mesh for void centres ..." << endl;
    chainmesh::ChainMesh3D ChM_voids(catalogue->Min(Var::_Radius_),
				     catalogue->var(Var::_X_),
				     catalogue->var(Var::_Y_),
				     catalogue->var(Var::_Z_),
				     2.*catalogue->Max(Var::_Radius_));
    for (size_t i = 0; i<catalogue->nObjects(); i++) {
      vector<long> close = ChM_voids.close_objects(catalogue->coordinate(i));
      unsigned int h = 0;
      while (!remove[i] && h<close.size()) {
	double distance = catalogue->distance(i, catalogue->catalogue_object(close[h]));
	if ((distance < catalogue->radius(i)+catalogue->radius(close[h]) && (int) i < close[h]) ||
	    ((int) i > close[h] && remove[close[h]])) h++;
	else if (distance < catalogue->radius(i)+catalogue->radius(close[h]) && (int) i > close[h]) remove[i] = true;
	else h++;
      }//while
    }//for
    int overlap_removed = 0;
    for (size_t i = 0; i<remove.size(); i++) if (remove[i]) overlap_removed++;
    catalogue->remove_objects(remove);
    coutCBL << "Voids removed to avoid overlap: " << overlap_removed << endl;
    coutCBL << "Time spent by the overlap-checking procedure: " << omp_get_wtime()-ol_time << " seconds" << endl;
  }//overlap check
  cout << endl;
  coutCBL << "Voids in the Catalogue: " << catalogue->nObjects() << endl;
  if (catalogue->nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 

  //float olchecking_time = float( clock() - begin_time ) / CLOCKS_PER_SEC;
  //olchecking_time = olchecking_time - rescaling_time - cleaning_time;
  cout << endl;
  coutCBL << "Total time spent: " << omp_get_wtime()-start_time << " seconds \n" << endl;

  m_object = catalogue->sample();
}

// ============================================================================

void cbl::catalogue::Catalogue::compute_centralDensity (const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM)
{
  const double density = tracers_catalogue->numdensity();
	const double mps = tracers_catalogue->mps();

  for (size_t j=0; j<nObjects(); j++) {
    ChM.get_searching_region(2*mps);
		int cont=0;
		vector<long> close = ChM.close_objects(coordinate(j));
    for (auto&& k : close) {
      double distance = cbl::catalogue::Catalogue::distance(j, tracers_catalogue->catalogue_object(k));
      if (distance <= 2*mps) cont++;
    }
		set_var(j, Var::_CentralDensity_, cont/(volume_sphere(2*mps)*density));
	}
  
}


// ============================================================================


void cbl::catalogue::Catalogue::compute_centralDensity (const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const std::vector<double> par_numdensity, const double ratio)
{
  //vector to memorize which element of the catalogue has to be removed at the end of the procedure:
  vector<bool> remove(m_object.size(), false);

  //counter for regions without any tracer:
  int void_voids = 0;
#pragma omp parallel num_threads(omp_get_max_threads())
  {
#pragma omp for ordered schedule(static)
    for (size_t j = 0; j<m_object.size(); j++) {
#pragma omp ordered  
      ChM.get_searching_region(m_object[j]->radius());
      vector<long> close = ChM.close_objects(m_object[j]->coords());

      //compute distances between the void and the surrounding particles
      vector<double> distances;
      for (auto&& k : close) {
	double distance = cbl::catalogue::Catalogue::distance(j, tracers_catalogue->catalogue_object(k));
	if (distance < m_object[j]->radius()) distances.emplace_back(distance);
      }
    
      //FIND CENTRAL DENSITY
      if (distances.size() > 0) {
	
	// redshift if the centre of the considered void
	double zz = m_object[j]->redshift();   
	
	// compute the number density using a Nth-order polynomial
	double density = 0.;
	for (size_t N=par_numdensity.size(); N-->0;) density += par_numdensity[par_numdensity.size()-1-N]*pow(zz,N);

	std::sort (distances.begin(), distances.end());
	int NN = 0;
	while (distances[NN]<ratio*m_object[j]->radius() && NN<(int)distances.size()-1) NN++; // check -1
	
	if (NN > 0) {
	  double central_density;

	  if (NN/(volume_sphere(distances[NN-1])*density)<(NN+1)/(volume_sphere(distances[NN])*density))
	    central_density = interpolated(ratio*m_object[j]->radius(),
	    {distances[NN-1], distances[NN]},
	      {NN/(volume_sphere(distances[NN-1])*density), (NN+1)/(volume_sphere(distances[NN])*density)}, "Linear");
	  else
	    central_density = interpolated(ratio*m_object[j]->radius(),
	    {distances[NN], distances[NN-1]},
	      {(NN+1)/(volume_sphere(distances[NN])*density), NN/(volume_sphere(distances[NN-1])*density)}, "Linear");
	
	  m_object[j]->set_centralDensity(central_density);
	}
	else m_object[j]->set_centralDensity(0.);	
      }
      
      else {
	vector<double>().swap(distances);
#pragma omp critical
	void_voids ++;
	remove[j] = true;
      }
      
    }//for
  }

  if (void_voids > 0) {
    for (size_t j=m_object.size(); j --> 0;) {
      if (remove[j]) m_object.erase(m_object.begin()+j);
    }
  }

  coutCBL << "I removed " << void_voids << " voids in calculating the central density!" << endl;
  
}

// ============================================================================


void cbl::catalogue::Catalogue::compute_centralDensity (const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const std::vector<std::vector<double>> data_numdensity, const std::string method_interpolation, const double ratio)
{
  
  //vector to memorize which element of the catalogue has to be removed at the end of the procedure:
  vector<bool> remove(m_object.size(), false);

  //counter for regions without any tracer:
  int void_voids = 0;
#pragma omp parallel num_threads(omp_get_max_threads())
  {
#pragma omp for ordered schedule(static)
    for (size_t j = 0; j<m_object.size(); j++) {
#pragma omp ordered  
      ChM.get_searching_region(m_object[j]->radius());
      vector<long> close = ChM.close_objects(m_object[j]->coords());

      //compute distances between the void and the surrounding particles
      vector<double> distances;
      for (auto&& k : close) {
	double distance = cbl::catalogue::Catalogue::distance(j, tracers_catalogue->catalogue_object(k));
	if (distance < m_object[j]->radius()) distances.emplace_back(distance);
      }
    
      //FIND CENTRAL DENSITY
      if (distances.size() > 0) {
	
	// redshift if the centre of the considered void
	double zz = m_object[j]->redshift();

	// compute the number density using a specified method
	double density = interpolated(zz, data_numdensity[0], data_numdensity[1], method_interpolation);

	std::sort (distances.begin(), distances.end());
	int NN = 0;
	while (distances[NN]<ratio*m_object[j]->radius() && NN<(int)distances.size()-1) NN++; // check -1
	
	if (NN > 0) {
	  double central_density;

	  if (NN/(volume_sphere(distances[NN-1])*density)<(NN+1)/(volume_sphere(distances[NN])*density))
	    central_density = interpolated(ratio*m_object[j]->radius(),
	    {distances[NN-1], distances[NN]},
	      {NN/(volume_sphere(distances[NN-1])*density), (NN+1)/(volume_sphere(distances[NN])*density)}, "Linear");
	  else
	    central_density = interpolated(ratio*m_object[j]->radius(),
	    {distances[NN], distances[NN-1]},
	      {(NN+1)/(volume_sphere(distances[NN])*density), NN/(volume_sphere(distances[NN-1])*density)}, "Linear");
	
	  m_object[j]->set_centralDensity(central_density);
	}
	else m_object[j]->set_centralDensity(0.);	
      }
      
      else {
	vector<double>().swap(distances);
#pragma omp critical
	void_voids ++;
	remove[j] = true;
      }
      
    }//for
  }

  if (void_voids > 0) {
    for (size_t j=m_object.size(); j --> 0;) {
      if (remove[j]) m_object.erase(m_object.begin()+j);
    }
  }

  coutCBL << "I removed " << void_voids << " voids in calculating the central density!" << endl;
  
}


// ============================================================================


void cbl::catalogue::Catalogue::compute_densityContrast (const shared_ptr<Catalogue> tracers_catalogue, cbl::chainmesh::ChainMesh3D ChM, const double ratio) {
  
  vector<bool> remove(nObjects(), false), void_voids(nObjects(), false), cloud_in_void(nObjects(), false);

#pragma omp parallel num_threads(omp_get_max_threads())
  {
#pragma omp for ordered schedule(static)
    for (size_t j = 0; j<nObjects(); j++) {
#pragma omp ordered
      ChM.get_searching_region(radius(j));
      vector<long> close = ChM.close_objects(coordinate(j));
      
      //compute distances between the void and the surrounding particles
      vector<double> distances;
      for (auto&& k : close) {
				double distance = Catalogue::distance(j, tracers_catalogue->catalogue_object(k));
				if (distance <= radius(j)) distances.emplace_back(distance);
      }

      std::sort(distances.begin(), distances.end());
      
      if (distances.size() > 1) {
      
				//FIND CENTRAL DENSITY
				int NN = 0;
				while (distances[NN]<ratio*radius(j) && NN<(int)distances.size()-1) NN++; // check-1

				if (NN > 0 && NN<(int)distances.size()-1) {
					double delta_in = NN/cbl::volume_sphere(distances[NN-1]);
					double delta_out = distances.size()/cbl::volume_sphere(radius(j));
					if (delta_out/delta_in <= 1.) {
						remove[j] = true;
	    			cloud_in_void[j] = true;
	    			continue;
	 				}
	  			else set_var(j, Var::_DensityContrast_, delta_out/delta_in);
				}
				else if (NN == 0) {
					set_var(j, Var::_DensityContrast_, 9999.);
					//remove[j] = true;
					//void_voids[j] = true;
				}
				else if (NN == (int)distances.size()-1) {
					remove[j] = true;
					cloud_in_void[j] = true;
	    		continue;
  			}   
  		} 
  		else {
				remove[j] = true;
				void_voids[j] = true;
      } 
    }//for
  }
  
  remove_objects(remove);
  	
  coutCBL << "Cloud-in-void: " << std::count(cloud_in_void.begin(), cloud_in_void.end(), true) << endl;
  coutCBL << "I removed " << std::count(cloud_in_void.begin(), cloud_in_void.end(), true)+std::count(void_voids.begin(), void_voids.end(), true) << " voids in calculating the density contrast!" << endl;
  
}