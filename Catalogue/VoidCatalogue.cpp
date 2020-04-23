/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli, Carlo Cannarozzo        *
 *  and Tommaso Ronconi                                             *
 *  federico.marulli3@unibo.it carlo.cannarozzo@studio.unibo.it     *
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
 *  @author Federico Marulli, Carlo Cannarozzo, Tommaso Ronconi
 *
 *  @author federico.marulli3@unibo.it,
 *  carlo.cannarozzo@studio.unibo.it, tommaso.ronconi@studio.unibo.it
 */


#include "Func.h"
#include "Catalogue.h"
#include "ChainMesh_Catalogue.h"
#include "Object.h"
#include "TwoPointCorrelation.h"
#include "TwoPointCorrelation1D_monopole.h"
#include <omp.h>

using namespace std;
using namespace cbl;

//********************************** Carlo Cannarozzo **********************************//

cbl::catalogue::Catalogue::Catalogue (const VoidAlgorithm algorithm, const Catalogue tracer_catalogue, const double nSub, const std::vector<std::string> random_catalogue_vector, const std::string dir_output, const std::string output, const double r_max, const double cellsize, const int n_rec, const int n_iter, const bool swapping, const bool add_unpaired, const double convergence_fact, const double step_size, const double gaussian_smoothing, const double protovoid_distance)
{
  
  // -------------------------------------------- //
  // ---------------- First Step ---------------- //
  // -------------------------------------------- //
    
  // * * * Displacement field reconstruction * * * //     
  
  //const clock_t begin_time = clock();
  const double start_time = omp_get_wtime();
  
  cosmology::Cosmology cosm;
  
  Catalogue initial_catalogue;
  Catalogue final_catalogue;
  Catalogue displacement_catalogue;
  
  if (algorithm==VoidAlgorithm::_LaZeVo_) { // Begin of LaZeVo method
    cout << endl; coutCBL << par::col_green << "LaZeVo algorithm" << par::col_default << endl << endl;
    Catalogue rand_catalogue; // 'Catalogue' constructor for random halo objects
    
    // per ora leggo il primo elemento del vettore di string
    string rand_file = random_catalogue_vector[0];
    
    coutCBL << "* * * Looking for closest pairs particles-random particles * * *" << endl << endl;
    
    // LOOP FOR THE DIFFERENT REALISATIONS OF THE RANDOM CATALOGUE (different seeds)
      
    for (int rndd=1; rndd<=n_rec; ++rndd) { // 'for' cycle - Cycling to the number of realizations 'n_rec' set up in the 'input_param.ini' file 
      coutCBL << rndd << " of " << n_rec << " random realisations " << endl;

      if (rand_file=="not_provided") {
        Catalogue temp_catalogue {RandomType::_createRandom_box_, tracer_catalogue, 1., 10, cosm, false, 10., {}, {}, {}, 10, rndd};
        rand_catalogue = temp_catalogue;
      }

      // the random catalogue is created when the user does not provide it
      else if (rand_file!="not_provided" && n_rec<=int(random_catalogue_vector.size())) {
        coutCBL << "Reading random catalogue ... " << endl;
        Catalogue temp_catalogue {ObjectType::_Random_, CoordinateType::_comoving_, {random_catalogue_vector[rndd-1]}, 1, 2, 3, -1, -1, nSub};
        coutCBL << " Number of halos extracted from " << rand_file << " catalogue: " << rand_catalogue.nObjects() << endl << endl;
        rand_catalogue = temp_catalogue;
      }

      else ErrorCBL ("Check your inputs for the random catalogue!", "Catalogue", "VoidCatalogue.cpp");

      // shared pointers to the tracer and the random catalogues
      auto halo_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(tracer_catalogue)));
      auto rand_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(rand_catalogue)));

      // maximum for the number of objects in the tracer catalogue (unsigned int) 
      unsigned long int max_index = halo_cat->nObjects();
      if (max_index > 4294967295) coutCBL << "Wait! Number of objects exeeding maximum value for an unsigned int (4294967295)... subsample your catalogue or change the function cbl::Erase()!" << endl;

      unsigned int num_objects = halo_cat->nObjects();

      // vector for the indices of haloes
      std::vector<unsigned int> index_halo_cat(num_objects);
      std::iota (std::begin(index_halo_cat), std::end(index_halo_cat), 0);
      
      coutCBL << "Generating chain mesh ... " << flush;
      cbl::chainmesh::ChainMesh_Catalogue ChM (cellsize, rand_cat, r_max);
      cout << "done! " << endl;
      
      vector<shared_ptr<cbl::catalogue::Object>> initial_catalogue_object;
      vector<shared_ptr<cbl::catalogue::Object>> final_catalogue_object;
      vector<shared_ptr<cbl::catalogue::Object>> displ_catalogue_object;
      
      cout << endl;
      coutCBL << "Performing the iterations... " << endl;

      //  TIME
      double start_iter = omp_get_wtime();

      vector<unsigned int> vect_k, vect_i, vect_i_unpaired;
      unsigned int min_unpaired = 4294967295; // max value for unsigned int
      int counter = 0;

      // mpi for different iterations?
      unsigned int iteration_number = n_iter; 	
      for (unsigned int iter=0; iter<iteration_number; iter++) // <<<
	{
	  unsigned int paired = 0;
	  unsigned int unpaired = 0;	  
	  vector<unsigned int> temp_vect_k(num_objects);  // maximum reservation of memory that will possibily be used
	  vector<unsigned int> temp_vect_i(num_objects);

	  counter+=1;
	  coutCBL << counter << " / " << iteration_number << " \r"; cout.flush(); // <<<

	  vector<shared_ptr<cbl::catalogue::Object>> initial_catalogue_object;
	  vector<shared_ptr<cbl::catalogue::Object>> final_catalogue_object;
	  vector<shared_ptr<cbl::catalogue::Object>> displacement_catalogue_object;
	  vector<bool> used(num_objects, false);

	  std::vector<unsigned int> index_halo_cat_temp = index_halo_cat;
	  std::shuffle(index_halo_cat_temp.begin(), index_halo_cat_temp.end(), std::default_random_engine(iter));
  
#pragma omp parallel num_threads(omp_get_max_threads())
	  {	    
#pragma omp for ordered schedule(dynamic)
	    for (auto i=index_halo_cat_temp.begin(); i<index_halo_cat_temp.end(); i++) { // loop on each object of the catalogue
 
#pragma omp ordered
	      ChM.get_searching_region(r_max); // with parallelisation this has to be inside the loop...
	      vector<long int> close_objects = ChM.close_objects(halo_cat->coordinate(*i)); // searching for objects near the halo
	      
	      if (close_objects.size()>0) {

		double distance = r_max*10.;
		unsigned int k_dist = -1;
		
		for (auto&& k : close_objects) { // searching for the closer random object that is not already been paired
		  if (!used[k]) {
		    double new_distance = cbl::Euclidean_distance(halo_cat->xx(*i), rand_cat->xx(k),
								  halo_cat->yy(*i), rand_cat->yy(k),
								  halo_cat->zz(*i), rand_cat->zz(k));	
		    if (new_distance < distance) {
		      distance = new_distance;
		      k_dist = k;
		    }
		    
		  }
		}

#pragma omp ordered
		{
		  // if the minimum distance is "acceptable"
		  if (distance<r_max*10.) { // we pair the halo index with the index of the closer random ojbect (that we record it as "used")
		    temp_vect_i[paired] = *i;
		    temp_vect_k[paired] = k_dist;         
		    used[k_dist] = true;
#pragma omp critical
		    paired++; // counting the paired objects
		  }
		  
		  else {
#pragma omp critical
		    unpaired++; // counting the paired objects (no free partner in the neighborhood)
		  }
		  
		} // omp ordered
	      }
	      
	      else {;
#pragma omp critical
		unpaired++; // counting other unpaired objects --> this case should never happen
	      }
	      
	    } // for over halo indices
	  }// pragma parallel

	  
	  if (unpaired<min_unpaired) {
	    min_unpaired = unpaired;
	    vect_i = temp_vect_i;
	    vect_k = temp_vect_k;
	  }
	  
	} // for iterations

      // re-size of the vectors to make them occupy the exact memory they are using (new vector lenght)
      vect_i.resize(num_objects-min_unpaired);
      vect_k.resize(num_objects-min_unpaired);

      // TIME
      coutCBL << "Time spent during the iteration process: " << (omp_get_wtime()-start_iter)/60. << " minutes " << endl << endl;

      coutCBL << "Number of unpaired halos: " << min_unpaired << endl;
      coutCBL << "Percentage of unpaired halos: " << (float)min_unpaired/num_objects * 100 << "% " << endl;
      cout << endl;
      

      if (swapping==true) {

	      
	//**********************************//
	//**                              **//      
	//**   ANDRII -> swapping pairs   **//
	//**                              **//
	//**********************************//
      
	coutCBL << "* * * Swapping procedure started! * * * " << endl << endl;

	double start_swap = omp_get_wtime();
	int epsilon;
	vector<unsigned int> vect_k_tot = vect_k;
	vector<unsigned int> vect_i_tot = vect_i;
	  
	if (add_unpaired) {

	  vector<unsigned int> vect_k_unpaired = index_halo_cat; // vector with all the indices (we use the initial index vector created for haloes)
	  vector<unsigned int> vect_i_unpaired = index_halo_cat;
	  
	  // remove from this vector (containing all the indices) the indices of the paired random objects
	  std::sort(vect_k.begin(),vect_k.end());
	  for (size_t j = vect_k.size(); j --> 0;) {
	    vect_k_unpaired.erase(vect_k_unpaired.begin()+vect_k[j]);
	  }

	  // remove from this vector (containing all the indices) the indices of the paired objects
	  std::sort(vect_i.begin(),vect_i.end());
	  for (size_t j = vect_i.size(); j --> 0;) {
	    vect_i_unpaired.erase(vect_i_unpaired.begin()+vect_i[j]);
	  }
	
	  // now we minimise the action decreasing the distance between the (previously) unpaired objects trying to swap these couples //
	
	  coutCBL << "Coupling the unpaired tracers and minimising the distance between the pairs... " << flush;

	  // we create a vector of int to use as random numbers (shuffled)
	  std::vector<unsigned int> index_unpaired(min_unpaired);
	  std::iota (std::begin(index_unpaired), std::end(index_unpaired), 0);

	  // value that determines the level of convergence
	  epsilon = int(convergence_fact*min_unpaired);

	  for (int n_for=0; n_for<epsilon; n_for++) {
	    std::shuffle(index_unpaired.begin(), index_unpaired.end(), std::default_random_engine(n_for));
	    
#pragma omp parallel num_threads(omp_get_max_threads())
	    {	
#pragma omp for schedule(static)
	      for (size_t index=0; index<min_unpaired-1; index+=2) {
	  
		unsigned int index_halo_1 = vect_i_unpaired[index_unpaired[index]];
		unsigned int index_halo_2 = vect_i_unpaired[index_unpaired[index+1]];
		unsigned int index_random_1 = vect_k_unpaired[index_unpaired[index]];
		unsigned int index_random_2 = vect_k_unpaired[index_unpaired[index+1]];

		// distance between two pairs --> couples: halo[random_1]-rand[random_1] and halo[random_2]-rand[random_2]
		double distance_old =
		  cbl::Euclidean_distance(halo_cat->xx(index_halo_1), rand_cat->xx(index_random_1),
					  halo_cat->yy(index_halo_1), rand_cat->yy(index_random_1),
					  halo_cat->zz(index_halo_1), rand_cat->zz(index_random_1)) +
		  cbl::Euclidean_distance(halo_cat->xx(index_halo_2), rand_cat->xx(index_random_2),
					  halo_cat->yy(index_halo_2), rand_cat->yy(index_random_2),
					  halo_cat->zz(index_halo_2), rand_cat->zz(index_random_2));
	  
	  
	        // distance between two pairs after the swap --> couples: halo[random_1]-rand[random_2] and halo[random_2]-rand[random_1]
		double distance_new =
		  cbl::Euclidean_distance(halo_cat->xx(index_halo_1), rand_cat->xx(index_random_2),
					  halo_cat->yy(index_halo_1), rand_cat->yy(index_random_2),
					  halo_cat->zz(index_halo_1), rand_cat->zz(index_random_2)) +
		  cbl::Euclidean_distance(halo_cat->xx(index_halo_2), rand_cat->xx(index_random_1),
					  halo_cat->yy(index_halo_2), rand_cat->yy(index_random_1),
					  halo_cat->zz(index_halo_2), rand_cat->zz(index_random_1));


		if (distance_new < distance_old) { // if the new distance is smaller we exchange the 2 indices in the vector of random ojbects
		  vect_k_unpaired[index_unpaired[index]] = index_random_2;
		  vect_k_unpaired[index_unpaired[index+1]] = index_random_1;
	  
		}
	      }
	    }
	  }

	  cout << "done!" <<  endl;
 
	  // we put together the indices of the paired and previously unpaired objects

	  vect_i_tot.insert(vect_i_tot.end(), vect_i_unpaired.begin(), vect_i_unpaired.end());
	  vect_k_tot.insert(vect_k_tot.end(), vect_k_unpaired.begin(), vect_k_unpaired.end());

	}

      
	// swapping the couples once again (this time using vectors with all the indices) //

	if (add_unpaired) coutCBL << "Minimising the distance between all the pairs... " << flush;
	else coutCBL << "Minimising the distance between the pairs... " << flush;

	// new convergence parameter
	epsilon = int(convergence_fact*vect_k_tot.size());

	// we create a vector of idices from 0 to the maximum number of objects
	std::vector<unsigned int> index_tot(vect_k_tot.size());
	std::iota (std::begin(index_tot), std::end(index_tot), 0);
 
	// how many time we loop over all the haloes
	for (int n_for=0; n_for<epsilon; n_for++) {
	  std::shuffle(index_tot.begin(), index_tot.end(), std::default_random_engine(n_for));
	  
#pragma omp parallel num_threads(omp_get_max_threads())
	  {	
#pragma omp for schedule(static)
	    for (size_t index=0; index<vect_k_tot.size()-1; index+=2) {
	  
	      unsigned int index_halo_1 = vect_i_tot[index_tot[index]];
	      unsigned int index_halo_2 = vect_i_tot[index_tot[index+1]];
	      unsigned int index_random_1 = vect_k_tot[index_tot[index]];
	      unsigned int index_random_2 = vect_k_tot[index_tot[index+1]];

	      // distance between two pairs --> couples: halo[random_1]-rand[random_1] and halo[random_2]-rand[random_2]
	      double distance_old =
		cbl::Euclidean_distance(halo_cat->xx(index_halo_1), rand_cat->xx(index_random_1),
					halo_cat->yy(index_halo_1), rand_cat->yy(index_random_1),
					halo_cat->zz(index_halo_1), rand_cat->zz(index_random_1)) +
		cbl::Euclidean_distance(halo_cat->xx(index_halo_2), rand_cat->xx(index_random_2),
					halo_cat->yy(index_halo_2), rand_cat->yy(index_random_2),
					halo_cat->zz(index_halo_2), rand_cat->zz(index_random_2));
	  
	  
	      // distance between two pairs after the swap --> couples: halo[random_1]-rand[random_2] and halo[random_2]-rand[random_1]
	      double distance_new =
		cbl::Euclidean_distance(halo_cat->xx(index_halo_1), rand_cat->xx(index_random_2),
					halo_cat->yy(index_halo_1), rand_cat->yy(index_random_2),
					halo_cat->zz(index_halo_1), rand_cat->zz(index_random_2)) +
		cbl::Euclidean_distance(halo_cat->xx(index_halo_2), rand_cat->xx(index_random_1),
					halo_cat->yy(index_halo_2), rand_cat->yy(index_random_1),
					halo_cat->zz(index_halo_2), rand_cat->zz(index_random_1));


	      if (distance_new < distance_old) { // if the new distance is smaller we exchange the 2 indices in the vector of random ojbects
		vect_k_tot[index_tot[index]] = index_random_2;
		vect_k_tot[index_tot[index+1]] = index_random_1;		    
	      }
	    }
	  }
	}

	cout << "done!" << endl << endl;

	coutCBL << "* * * Swapping procedure ended! * * * " << endl << endl;


	// TIME
	coutCBL << "Time spent during the iteration process: " << (omp_get_wtime()-start_swap)/60. << " minutes " << endl;
      

	// creating catalogues with the initial and final positions and with the displacement values
	for (size_t ii=0; ii<vect_k_tot.size(); ii++) {
	
	  cbl::comovingCoordinates initial_coordinates = {halo_cat->xx(vect_i_tot[ii]), halo_cat->yy(vect_i_tot[ii]), halo_cat->zz(vect_i_tot[ii])};    
	  auto inital_halo_coord = make_shared<cbl::catalogue::Halo>(initial_coordinates);
	  initial_catalogue_object.push_back(inital_halo_coord);
            
	  cbl::comovingCoordinates final_coordinates = {rand_cat->xx(vect_k_tot[ii]), rand_cat->yy(vect_k_tot[ii]), rand_cat->zz(vect_k_tot[ii])};
	  auto halo_coord = make_shared<cbl::catalogue::Halo>(final_coordinates);
	  final_catalogue_object.push_back(halo_coord);
	    
	  cbl::comovingCoordinates displacement_values = {rand_cat->xx(vect_k_tot[ii])-halo_cat->xx(vect_i_tot[ii]), rand_cat->yy(vect_k_tot[ii])-halo_cat->yy(vect_i_tot[ii]), rand_cat->zz(vect_k_tot[ii])-halo_cat->zz(vect_i_tot[ii])};
	  auto halo_displ = make_shared<cbl::catalogue::Halo>(displacement_values);
	  displ_catalogue_object.push_back(halo_displ);
	}

      }

      else {

	for (size_t ii = 0; ii<vect_k.size(); ii++)
	  {
	    cbl::comovingCoordinates initial_coordinates = {halo_cat->xx(vect_i[ii]), halo_cat->yy(vect_i[ii]), halo_cat->zz(vect_i[ii])};
	    
	    auto inital_halo_coord = make_shared<cbl::catalogue::Halo>(initial_coordinates);
	    initial_catalogue_object.push_back(inital_halo_coord);
            
	    cbl::comovingCoordinates final_coordinates = {rand_cat->xx(vect_k[ii]), rand_cat->yy(vect_k[ii]), rand_cat->zz(vect_k[ii])};
	    auto final_halo_coord = make_shared<cbl::catalogue::Halo>(final_coordinates);
	    final_catalogue_object.push_back(final_halo_coord);
	    
	    cbl::comovingCoordinates displacement_values = {rand_cat->xx(vect_k[ii])-halo_cat->xx(vect_i[ii]), rand_cat->yy(vect_k[ii])-halo_cat->yy(vect_i[ii]), rand_cat->zz(vect_k[ii])-halo_cat->zz(vect_i[ii])};
	    auto halo_displ = make_shared<cbl::catalogue::Halo>(displacement_values);
	    displ_catalogue_object.push_back(halo_displ);
	  }
      

      }
      
      
      
      stringstream ss;
      ss << rndd;
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
      for (unsigned int i=0; i<final_catalogue_object.size(); i++)
	{
	  fout_displ.precision(7);
	  fout_displ << initial_catalogue_object[i]->xx() << " " << initial_catalogue_object[i]->yy() << " " << initial_catalogue_object[i]->zz() << " "
		     << final_catalogue_object[i]->xx() << " " << final_catalogue_object[i]->yy() << " " << final_catalogue_object[i]->zz() << " "
		     << displ_catalogue_object[i]->xx() << " " << displ_catalogue_object[i]->yy() << " " << displ_catalogue_object[i]->zz() << endl;   
	} // modified
      
      fout_displ.clear(); fout_displ.close();
      

      initial_catalogue.add_objects(initial_catalogue_object);
      final_catalogue.add_objects(final_catalogue_object);
      displacement_catalogue.add_objects(displ_catalogue_object);

    }// for random realisations
    
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
    
  auto start_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(initial_catalogue)));
  auto final_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(final_catalogue)));
  auto displ_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(displacement_catalogue)));

  coutCBL << "Number of all paired halos: " << final_catalogue.nObjects() << endl << endl;
  

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
  
  string name;
  //   string output;
  
  //   output = "coord_LCDM.dat";
  
  name = dir_output + "Divergence";
  name.append("_");
  name.append(output);
  ofstream fout_divergence(name.c_str());

  fout_divergence.precision(7);

  
  double min_x = final_cat->catalogue::Catalogue::Min(Var::_X_), max_x = final_cat->catalogue::Catalogue::Max(Var::_X_);
  double min_y = final_cat->catalogue::Catalogue::Min(Var::_Y_), max_y = final_cat->catalogue::Catalogue::Max(Var::_Y_);
  double min_z = final_cat->catalogue::Catalogue::Min(Var::_Z_), max_z = final_cat->catalogue::Catalogue::Max(Var::_Z_);
  
  double delta_x = max_x-min_x;
  double delta_y = max_y-min_y;
  double delta_z = max_z-min_z;
  
  double volume = tracer_catalogue.volume();
  double density = tracer_catalogue.nObjects()/volume; 
  
  double fraction = 1./3.;
  double step = step_size*tracer_catalogue.mps(); // 1.5*pow(density,-fraction) --> corresponding to m.i.s. (it can be decreased...)
  double substep = 1.e-3*step;
  double substep_1 = 1./substep;
  double substep10_1 = 0.1/substep;
  double fact_sub = 10*substep;

  double sigma = gaussian_smoothing*tracer_catalogue.mps(); // 1.5*step --> gaussian smoothing gaussiano to obtain a continuous vector field
  double sg = 1./(2.*sigma*sigma);

  coutCBL << "Δx     = " << delta_x << endl;
  coutCBL << "Δy     = " << delta_y << endl;
  coutCBL << "Δz     = " << delta_z << endl;
  coutCBL << "Volume = " << volume  << endl;
  coutCBL << "n      = " << density << endl;
  coutCBL << "step   = " << step    << endl;
  coutCBL << "σ      = " << sigma   << endl;
  cout << endl;
  
  int mem_alloc_x = int(delta_x/step); 
  int mem_alloc_y = int(delta_y/step);
  int mem_alloc_z = int(delta_z/step);

  Tensor3Di Pn         (mem_alloc_x, vector<vector<int>>          (mem_alloc_y, vector<int>          (mem_alloc_z, 0)));
  Tensor3Dd Divergence (mem_alloc_x, vector<vector<double>>       (mem_alloc_y, vector<double>       (mem_alloc_z, 0.)));
  Tensor3Dd Gradient_x (mem_alloc_x, vector<vector<double>>       (mem_alloc_y, vector<double>       (mem_alloc_z, 0.)));
  Tensor3Dd Gradient_y (mem_alloc_x, vector<vector<double>>       (mem_alloc_y, vector<double>       (mem_alloc_z, 0.)));
  Tensor3Dd Gradient_z (mem_alloc_x, vector<vector<double>>       (mem_alloc_y, vector<double>       (mem_alloc_z, 0.)));
  Tensor4Di Pp         (mem_alloc_x, vector<vector<vector<int>>>  (mem_alloc_y, vector<vector<int>>  (mem_alloc_z, vector<int>(0, 0))));

  
  // * * * Preparation for divergence field construction * * * //

  for (size_t o=0; o<final_cat->nObjects(); o++) // Looking in which grid halo locates
    {
      int io = round((final_cat->xx(o)-min_x)/step);  
      int jo = round((final_cat->yy(o)-min_y)/step);
      int ko = round((final_cat->zz(o)-min_z)/step);
    
      for (int i=-3; i<=3; ++i)
	{ 
	  for (int j=-3; j<=3; ++j)
	    {
	      for (int k=-3; k<=3; ++k)
		{            
		  if (i*i+j*j+k*k<=3.001*3.001 && io+i>0 && jo+j>0 && ko+k>0 && io+i<mem_alloc_x && jo+j<mem_alloc_y && ko+k<mem_alloc_z) // modified
		    {
		      Pn[io+i][jo+j][ko+k]++;
		      Pp[io+i][jo+j][ko+k].push_back(o); // modified
		    }
		}
	    }
	}
    }
  
  // * * * Divergence field construction * * * //
  
  cout << endl;
  coutCBL << "* * * Divergence field construction * * *" << endl << endl;
  coutCBL << "Wait...\r"; cout.flush();

  // searching for the cells with negative divergence
  
  int number = 0;
  
  double temp_x = 0., temp_y = 0., temp_z = 0.;
  vector<double> divergence(4, 0.0);

  vector<double> x_void, y_void, z_void, divergence_void; // modified
  vector<double> x_locmin, y_locmin, z_locmin;

  temp_x=min_x;
  while (temp_x<max_x)
    {
      double counter = 10./(max_x-step);	
      double percentage = temp_x-min_x;
      coutCBL << "..." << int(percentage*counter)*10 << "% completed \r"; cout.flush(); 
      temp_x+=step;
      temp_y=min_y;
    
      while (temp_y<max_y)
	{
	  temp_y+=step;
	  temp_z=min_z;
      
	  while (temp_z<max_z)
	    {
	      temp_z+=step;

	      vector<vector<double>> dist (4, vector<double>(4, 0.));
	      vector<vector<double>> v_0 (3, vector<double>(4, 0.));
	      vector<vector<double>> v_ (3, vector<double>(4, 0.));
	      vector<vector<double>> sum (4, vector<double>(4, 0.));
	      
	      int io = round((temp_x-min_x)/step);
	      int jo = round((temp_y-min_y)/step);
	      int ko = round((temp_z-min_z)/step);
        
	      if (io>0 && jo>0 && ko>0 && io<mem_alloc_x && jo<mem_alloc_y && ko<mem_alloc_z)
		{
		  for (int no=0; no<Pn[io][jo][ko]; ++no) // Cycle through all halos  // in origin this loop started from 1 and had the <= condition
		    {
		      int o = Pp[io][jo][ko][no]; 

		      vector<double> final_coord = {final_cat->xx(o), final_cat->yy(o), final_cat->zz(o)};
		      for (int LL=0; LL<=3; ++LL) {
			for (int KK=0; KK<=3; ++KK) {
			  vector<double> Temp = {temp_x, temp_y, temp_z};
			  if (KK!=0) Temp[KK-1] = Temp[KK-1]+substep; 
			  if (LL!=0) Temp[LL-1] = Temp[LL-1]+fact_sub;
			  dist[LL][KK] = pow((final_coord[0]-Temp[0]),2) + pow((final_coord[1]-Temp[1]),2) + pow((final_coord[2]-Temp[2]),2);
			}
		      }
		      
		      vector<double> displ_coord = {displ_cat->xx(o), displ_cat->yy(o), displ_cat->zz(o)};

		      for (int ii=0; ii<=3; ii++) {
			if (ii!=3) v_0[ii][0] += displ_coord[ii]*exp(-dist[0][0]*sg);
			sum[ii][0] += exp(-dist[0][ii]*sg);
			for (int nn=0; nn<=2; nn++) {
			  if (ii!=0) v_0[ii-1][nn+1] += displ_coord[nn]*exp(-dist[ii][0]*sg);
			  sum[ii][nn+1] += exp(-dist[nn+1][ii]*sg);
			  v_[nn][ii] += displ_coord[nn]*exp(-dist[ii][nn+1]*sg);
			}
		      }
		      
		    }
		  

		  for (int nn=0; nn<=2; nn++) {
		    v_0[nn][0] = v_0[nn][0]/sum[0][0];
		    v_[nn][0] = v_[nn][0]/sum[nn+1][0];
		    for (int ii=1; ii<=3; ii++) {
		      v_0[ii-1][nn+1] = v_0[ii-1][nn+1]/sum[0][ii];
		      v_[nn][ii] = v_[nn][ii]/sum[nn+1][ii];
		    }
		  }


		  divergence[0] = (v_[0][0]-v_0[0][0] + v_[1][0]-v_0[1][0] + v_[2][0]-v_0[2][0])*substep_1;	  
		  for (int ii=1; ii<=3; ii++) divergence[ii] = (v_[0][ii]-v_0[ii-1][1] + v_[1][ii]-v_0[ii-1][2] + v_[2][ii]-v_0[ii-1][3])*substep_1;
	  
		  Gradient_x[io][jo][ko] = (divergence[1]-divergence[0])*substep10_1;
		  Gradient_y[io][jo][ko] = (divergence[2]-divergence[0])*substep10_1;
		  Gradient_z[io][jo][ko] = (divergence[3]-divergence[0])*substep10_1;

		  Divergence[io][jo][ko] = divergence[0];

		  if (divergence[0]<0.)
		    {

		      // modified
		      
		      x_void.push_back(temp_x);
		      y_void.push_back(temp_y);
		      z_void.push_back(temp_z);
		      divergence_void.push_back(divergence[0]);

		      number++;
		    }
          
		  fout_divergence << temp_x << " " << temp_y << " " << temp_z << " " << divergence[0] << endl; // Divergence field for each grid point;
		}
	    }
	}
    }
  
  fout_divergence.clear(); fout_divergence.close();

  coutCBL << "Number of cells with negative divergence: " << number << endl << endl;

  // searching for the local minima

  temp_x = min_x;
  
  while (temp_x<max_x)
    {
      temp_x+=step;
      temp_y=min_y;
    
      while (temp_y<max_y)
	{
	  temp_y+=step;     
	  temp_z=min_z;
      
	  while (temp_z<max_z)
	    {
	      temp_z+=step;
        
	      int io = round((temp_x-min_x)/step);
	      int jo = round((temp_y-min_y)/step);
	      int ko = round((temp_z-min_z)/step);

	      int index = 0;
        
	      for (int i=-1; i<=1; ++i)
		{   
		  for (int j=-1; j<=1; ++j)
		    {
		      for (int k=-1; k<=1; ++k)
			{              
			  if (io+i>0 && jo+j>0 && ko+k>0 && io+i<mem_alloc_x && jo+j<mem_alloc_y && ko+k<mem_alloc_z && io>0 && jo>0 && ko>0 && io<mem_alloc_x && jo<mem_alloc_y && ko<mem_alloc_z)
			    {
			      if (Divergence[io+i][jo+j][ko+k]>Divergence[io][jo][ko] && Divergence[io][jo][ko]<0.) index++; // Condition on minimum value
			    }   // in origin Divergence[io][jo][ko]<-1e-10
			}
		    }
		}
        
	      if (index==26) // therefore if it completed all the if (3^3-1)
		{
		  x_locmin.push_back(io*step+min_x);
		  y_locmin.push_back(jo*step+min_y);
		  z_locmin.push_back(ko*step+min_z);

		}
	    }
	}
    }
  

  int num_voids = x_locmin.size();
  if (num_voids==0) 
    ErrorCBL("No local minima found!", "Catalogue", "VoidCatalogue.cpp");
  
   
  // -------------------------------------------- //
  // ---------------- Third Step ---------------- //
  // -------------------------------------------- //
  

  // Identification of voids //  --> this part has totally changed...

  coutCBL << "* * * Identification of voids * * *" << endl << endl;

  double range_cells = protovoid_distance*tracer_catalogue.mps();
  double cell_volume = step*step*step;
  vector<int> cells_in_subvoids(num_voids, 0);
  //double min_vol = pow(((3./(4.*M_PI))*sigma), fraction);
  vector<double> x_void_fin(num_voids);
  vector<double> y_void_fin(num_voids);
  vector<double> z_void_fin(num_voids);
  vector<double> radius_void(num_voids, 0.);

  for (int ID=0; ID<num_voids; ++ID) // loop over all the cells corresponding to the local minima
    {
      vector<int> peripheral_cells(number, -1), used_cells_id;
      vector<int> near_cells;
      
      for (int temp_number=0; temp_number<number; ++temp_number) // loop over all the cells with negative divergence
	{
	  // vectors containing all the positions of cells with div<0
	  temp_x = x_void[temp_number];
	  temp_y = y_void[temp_number];
	  temp_z = z_void[temp_number];

	  // for each void we record all the cells with div<0 and with a distance lower than r_max from the local minimum
	  if (cbl::Euclidean_distance(temp_x, x_locmin[ID], temp_y, y_locmin[ID], temp_z, z_locmin[ID])<r_max)  near_cells.push_back(temp_number);

	  // selecting all the cells inside a sphere with a radius lenght selected by the user (minimum size for a void)
	  if (cbl::Euclidean_distance(temp_x, x_locmin[ID], temp_y, y_locmin[ID], temp_z, z_locmin[ID])<range_cells)
	    {
	      cells_in_subvoids[ID]++;                     // counting the cells in each void
	      used_cells_id.push_back(temp_number);        // recording the indices of the used cells
	      peripheral_cells[temp_number]=temp_number;   // storage of the used cell indices (replacing in the vector the -1 with their index at the corresponding position)
		  
	    }
	}
     
      //if (pow(((3./(4.*M_PI))*cell_volume*double(cells_in_subvoids[ID])), fraction) < min_vol) break; // not good for openMP
      // if the volume of this prototype of void is too small we discard it

      //else { // otherwise, if the void is sufficiently large, we try to expand it adding to the enclosed cells new adjacent cells with div<0

      bool control = false;
	
      while (control==false) {
	  
	int new_counter = 0;

	for(auto&& temp_number : near_cells) { // loop over the indices of the cells thar are inside a radius of lenght r_max

	  temp_x = x_void[temp_number];
	  temp_y = y_void[temp_number];
	  temp_z = z_void[temp_number];
	    
	  if (peripheral_cells[temp_number]==-1) // check on all the cells with div<0 that have not already been added to this void
	    {
		
	      for(auto&& used_id : used_cells_id) // loop over the indices of the cells already enclosed in the (proto-)void
		// (does it a cell with div<0 adjacent to (or near) to any cell currently forming a void exist?) 
		{

		  if (cbl::Euclidean_distance(temp_x, x_void[used_id], temp_y, y_void[used_id], temp_z, z_void[used_id])<step*2.1)
		    {
		      used_cells_id.push_back(temp_number);
		      peripheral_cells[temp_number]=temp_number;
		      new_counter++;
		      break; // to avoid possible double countings
		    }
		}
	      
	    }
	    
	}// for near cells

	if (new_counter>0) { // if at least 1 new cell is found during this check the while loop continues (with the new enclosed cells to form the void)
	  cells_in_subvoids[ID]+=new_counter;
	}

	else {
	  control=true;  // when there are no more cells to enclose we exit from the while loop
	  radius_void[ID] = pow(((3./(4.*M_PI))*cell_volume*cells_in_subvoids[ID]), fraction); // compute the radius of the ID-th void
	  x_void_fin[ID] = x_locmin[ID];
	  y_void_fin[ID] = y_locmin[ID];
	  z_void_fin[ID] = z_locmin[ID];
	}
      }
	
      //} // else (adding adjacent cells)
      
    }// for ID
  
  name = dir_output+"Voids_";
  name.append(output);
  ofstream fout_voids(name.c_str());

  for(int ii=0; ii<(int)radius_void.size(); ++ii)  fout_voids << x_void_fin[ii] << " " << y_void_fin[ii] << " " << z_void_fin[ii] << " " << radius_void[ii] << endl;
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
}



//********************************** Tommaso Ronconi **********************************//

cbl::catalogue::Catalogue::Catalogue (const std::shared_ptr<Catalogue> input_voidCatalogue, const std::vector<bool> clean, const std::vector<double> delta_r, const double threshold, const double statistical_relevance, bool rescale, const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const double ratio, const bool checkoverlap, const Var ol_criterion)
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
    
    tracers_catalogue->compute_catalogueProperties();
    double density = tracers_catalogue->numdensity();
    double volume = tracers_catalogue->volume();
    
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
    catalogue->compute_centralDensity(tracers_catalogue, ChM, volume, ratio);
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
    catalogue->compute_centralDensity(tracers_catalogue, ChM, Volume, ratio);
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


// ============================================================================


void cbl::catalogue::Catalogue::compute_centralDensity (const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const double Volume, const double ratio)
{
  //vector to memorize which element of the catalogue has to be removed at the end of the procedure:
  vector<bool> remove(m_object.size(), false);

  //compute the number density of the tracer catalogue
  const double density = m_object.size()/Volume;
  
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


void cbl::catalogue::Catalogue::compute_densityContrast (const shared_ptr<Catalogue> tracers_catalogue, cbl::chainmesh::ChainMesh3D ChM, const double ratio) {
  
  //vector to memorize which element of the catalogue has to be removed at the end of the procedure:
  vector<bool> remove(m_object.size(), false);
    
  //counter for regions without any tracer:
  int void_voids = 0;
  int cloud_in_void = 0;

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
	std::sort(distances.begin(), distances.end());
	int NN = 0;
	while (distances[NN]<ratio*m_object[j]->radius() && NN<(int)distances.size()-1) NN++; // check -1

	if (NN > 0) {
	  double delta_in = NN/cbl::volume_sphere(distances[NN]);
	  double delta_out = distances.size()/cbl::volume_sphere(distances[distances.size()-1]);
	  if (delta_out/delta_in < 1.) {
	    remove[j] = true;
#pragma omp critical
	    cloud_in_void++;
	  }
	  else m_object[j]->set_densityContrast(delta_out/delta_in);
	}
	else m_object[j]->set_densityContrast(1.);
      }    
      else {
	vector<double>().swap(distances);
	remove[j] = true;
#pragma omp critical
	void_voids++;
      } 
    }//for
  }
  
  if (void_voids > 0 || cloud_in_void > 0) {
    for (size_t j=m_object.size(); j --> 0;) {
      if (remove[j]) m_object.erase(m_object.begin()+j);
    }
  }
  coutCBL << "Cloud-in-void: " << cloud_in_void << endl;
  coutCBL << "I removed " << void_voids+cloud_in_void << " voids in calculating the density contrast!" << endl;
  
}

