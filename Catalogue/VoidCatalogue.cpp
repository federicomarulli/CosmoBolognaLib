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
 *  @author federico.marulli3@unbo.it,
 *  carlo.cannarozzo@studio.unibo.it, tommaso.ronconi@studio.unibo.it
 */


// * * * * * * * * * * * * * * * * * * * * * //
// * * * LAST UPDATE -- 28th June 2017 * * * //
// * * * * * * * * * * * * * * * * * * * * * //


#include "Func.h"
#include "Catalogue.h"
#include "ChainMesh_Catalogue.h"
#include "Object.h"
#include "TwoPointCorrelation.h"
#include "TwoPointCorrelation1D_monopole.h"

using namespace std;

using namespace cbl;

/////////////////////////////////// Carlo Cannarozzo //////////////////////////////////////////

// ============================================================================

/// 3D Tensor of int
typedef vector< vector< vector <int> > >            Tensor3Di;

/// 3D Tensor of double
typedef vector< vector< vector <double> > >         Tensor3Dd;

/// 4D Tensor of int
typedef vector< vector< vector < vector <int> > > > Tensor4Di;  

/// @cond extvoid

cbl::catalogue::Catalogue::Catalogue (const VoidAlgorithm algorithm, const Catalogue halo_catalogue, const vector<string> file, const double nSub, const int n_rec, const string mode, const string dir_output, const string output, const double rmax, const int cellsize, const int n_iter, const double delta_movement)
{
  // -------------------------------------------- //
  // ---------------- First Step ---------------- //
  // -------------------------------------------- //
    
  // * * * Displacement field reconstruction * * * //     
  
  const clock_t begin_time = clock();
  
  cosmology::Cosmology cosm;
  
  Catalogue initial_catalogue;
  Catalogue total_catalogue;
  Catalogue displacement_catalogue;
  
  if (algorithm==VoidAlgorithm::_LaZeVo_) // Begin of LaZeVo method
  {
    cout << endl; coutCBL << par::col_green << "LaZeVo algorithm" << par::col_default << endl << endl;

    Catalogue rand_catalogue; // 'Catalogue' constructor for random halo objects

    string rand_file;
    
    for (size_t dd=0; dd<file.size(); ++dd) rand_file = file[dd];

    if (rand_file!="not_provided") 
    {
      Catalogue rand_tempcatalogue {ObjectType::_Random_, CoordinateType::_comoving_, {rand_file}, 1, 2, 3, -1, -1, nSub};
      coutCBL << " Number of halos extracted from " << rand_file << " catalogue: " << rand_tempcatalogue.nObjects() << endl << endl;
      rand_catalogue = rand_tempcatalogue;
    }
    
    coutCBL << "* * * Looking for closest pairs particles-random particles and calculating displacement field * * *" << endl << endl;
    
    for (int rndd=1; rndd<=n_rec; ++rndd) // 'for' cycle - Cycling to the number of realizations 'n_rec' set up in the 'input_param.ini' file 
    {
      int seed = rndd;
      
      coutCBL << rndd << " of " << n_rec << " random realizations " << endl;

      if (rand_file=="not_provided" && (mode=="non_periodic"||mode=="periodic"))
      {
        Catalogue rand_tempcatalogue {RandomType::_createRandom_box_, halo_catalogue, 1., 10, cosm, false, 10., {}, {}, {}, 10, seed};
        
        rand_catalogue = rand_tempcatalogue;      
      }
      
//       vector<bool> used(rand_catalogue.nObjects(), false);

//       stringstream ss;
//       ss << rndd;
//       string str = ss.str();
//       string output = ss.str();
//       string name = dir_output + "Displacement";
//       name.append("_");
//       name.append(str);
//       name.append("_");
//       output = "coord_LCDM";
//       name.append(output);
//       ofstream fout_displ(name.c_str());

      auto halo_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(halo_catalogue)));
      auto rand_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(rand_catalogue)));
      
      int unpaired=0; int old_unpaired = __INT_MAX__;
      
      vector<int> index_halo_cat;
      for (size_t i = 0; i<halo_cat->nObjects(); i++) index_halo_cat.push_back(i);
            
      cbl::chainmesh::ChainMesh_Catalogue ChM (cellsize, rand_cat, rmax);
      
      vector<shared_ptr<cbl::catalogue::Object>> initial_catalogue_object;
      vector<shared_ptr<cbl::catalogue::Object>> total_catalogue_object;
      vector<shared_ptr<cbl::catalogue::Object>> displacement_catalogue_object;
      
      cout << endl;
      unsigned int iteration_number= n_iter; /// <<<
      for (unsigned int iter=0; iter<iteration_number; iter++) /// <<<
      {        
        vector<shared_ptr<cbl::catalogue::Object>> temp_initial_catalogue_object;
        vector<shared_ptr<cbl::catalogue::Object>> temp_total_catalogue_object;
        vector<shared_ptr<cbl::catalogue::Object>> temp_displacement_catalogue_object;
        
        coutCBL << iter+1 << " / "<< iteration_number << " \r"; cout.flush(); /// <<<
        
        int iter_seed = iter;
        vector<bool> used(rand_catalogue.nObjects(), false);

        shuffle(index_halo_cat.begin(),index_halo_cat.end(), std::default_random_engine(iter_seed));

        unpaired = 0;

        for (auto&& i : index_halo_cat)
        {
          ChM.get_searching_region(rmax);
          vector<long> close_objects = ChM.close_objects(halo_cat->coordinate(i));
              
//           double counter = 100./halo_cat->nObjects();	
//           coutCBL << "..." << int(i*counter) << "% completed \r"; cout.flush(); 
          
          if (close_objects.size()>0)
          {
            vector<double> distance;
            vector<double> vect_i;
            vector<double> vect_k;
            vect_i.push_back(i);
                      
            for (auto&& k : close_objects)
            {
              if (!used[k])
              {
                distance.push_back(cbl::Euclidean_distance(halo_cat->xx(i), rand_cat->xx(k), halo_cat->yy(i), rand_cat->yy(k), halo_cat->zz(i), rand_cat->zz(k)));
                vect_k.push_back(k);
                vect_i.push_back(i);
              }
            }
            
            if (distance.size()>0)
            {
              vector<double>::iterator p1 = distance.begin(), p2 = vect_i.begin(), p3 = vect_k.begin();

              sort_3vectors(p1, p2, p3, distance.size());

  //             fout_displ.precision(7);
  //             fout_displ<<halo_cat->xx(vect_i[0])<<" "<<halo_cat->yy(vect_i[0])<<" "<<halo_cat->zz(vect_i[0])<<" "<<rand_cat->xx(vect_k[0])<<" "<<rand_cat->yy(vect_k[0])<<" "<<rand_cat->zz(vect_k[0])<<" "<<rand_cat->xx(vect_k[0])-halo_cat->xx(vect_i[0])<<" "<<rand_cat->yy(vect_k[0])-halo_cat->yy(vect_i[0])<<" "<<rand_cat->zz(vect_k[0])-halo_cat->zz(vect_i[0])<<" "<<distance[0]<<endl;

              cbl::comovingCoordinates initial_coordinates = {halo_cat->xx(vect_i[0]), halo_cat->yy(vect_i[0]), halo_cat->zz(vect_i[0])};
              auto inital_halo_coord = make_shared<cbl::catalogue::Halo>(initial_coordinates);
              temp_initial_catalogue_object.push_back(inital_halo_coord);
              
              cbl::comovingCoordinates final_coordinates = {rand_cat->xx(vect_k[0]), rand_cat->yy(vect_k[0]), rand_cat->zz(vect_k[0])};
              auto halo_coord = make_shared<cbl::catalogue::Halo>(final_coordinates);
              temp_total_catalogue_object.push_back(halo_coord);

              cbl::comovingCoordinates displacement_values = {rand_cat->xx(vect_k[0])-halo_cat->xx(vect_i[0]), rand_cat->yy(vect_k[0])-halo_cat->yy(vect_i[0]), rand_cat->zz(vect_k[0])-halo_cat->zz(vect_i[0])};
              auto halo_displ = make_shared<cbl::catalogue::Halo>(displacement_values);
              temp_displacement_catalogue_object.push_back(halo_displ);
              
              used[vect_k[0]] = true;            
            }
            
            else unpaired++;
          }
        }

        if(unpaired<old_unpaired)
        {
          initial_catalogue_object.erase(initial_catalogue_object.begin(), initial_catalogue_object.end());
          total_catalogue_object.erase(total_catalogue_object.begin(), total_catalogue_object.end());
          displacement_catalogue_object.erase(displacement_catalogue_object.begin(), displacement_catalogue_object.end());

          initial_catalogue_object=temp_initial_catalogue_object;
          total_catalogue_object=temp_total_catalogue_object;
          displacement_catalogue_object=temp_displacement_catalogue_object;

          old_unpaired = unpaired;

//           coutCBL << iter+1 << " / 2500 - lower: " << unpaired << endl;
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
      
      for (unsigned int i=0; i<total_catalogue_object.size(); i++)
      {
        fout_displ.precision(7);
        fout_displ<<initial_catalogue_object[i]->xx()<<" "<<initial_catalogue_object[i]->yy()<<" "<<initial_catalogue_object[i]->zz()<<" "<<total_catalogue_object[i]->xx()<<" "<<total_catalogue_object[i]->yy()<<" "<<total_catalogue_object[i]->zz()<<" "<<total_catalogue_object[i]->xx()-initial_catalogue_object[i]->xx()<<" "<<total_catalogue_object[i]->yy()-initial_catalogue_object[i]->yy()<<" "<<total_catalogue_object[i]->zz()-initial_catalogue_object[i]->zz()<<endl;//" "<<distance[0]<<endl;    
      }
      
      fout_displ.clear(); fout_displ.close();
      
      unpaired = old_unpaired;

      initial_catalogue.add_objects(initial_catalogue_object);
      total_catalogue.add_objects(total_catalogue_object);
      displacement_catalogue.add_objects(displacement_catalogue_object);

      coutCBL << "Number of unpaired halos: " << unpaired << endl;
      coutCBL << "Percentage of unpaired halos: " << (float)unpaired/halo_cat->nObjects() * 100 << "% " << endl;
      cout << endl;
    }
  } // End of LaZeVo method

  else if (algorithm==VoidAlgorithm::_RIVA_) // Begin of RIVA method
  {

    (void)delta_movement;
    /*
    cout << endl; coutCBL << par::col_green << "RIVA algorithm" << par::col_default << endl << endl;

    // * * * Two Point Correlation Function * * * //

    const double N_R = 1.; // random/data ratio

    const cbl::catalogue::Catalogue test_random_catalogue {cbl::catalogue::_box_, halo_catalogue, N_R};

    const double rMin  = 1e3; // minimum separation 
    const double rMax  = 5e4; // maximum separation 
    const int nbins    = 50;  // number of bins
    const double shift = 0.5; // spatial shift used to set the bin centre

    string w_file = "xi_RIVA_pre.dat";
    cbl::twopt::TwoPointCorrelation1D_monopole TPCF_catalogue {halo_catalogue, test_random_catalogue, cbl::_logarithmic_, rMin, rMax, nbins, shift};
// // //     TPCF_catalogue.measure_for_RIVA(cbl::twopt::_Poisson_, dir_output);
    TPCF_catalogue.measure(cbl::twopt::_Poisson_, dir_output);

    TPCF_catalogue.write(dir_output, w_file);
    
    
// // // // // WHILE CYCLE - - - while (2pcf_peak/max(xi[i])<1e-2) {...}
    
    auto halo_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(halo_catalogue)));

    vector<shared_ptr<cbl::catalogue::Object>> total_catalogue_object; // vector shared pointer for the constructor 'total_catalogue'
    
    vector<double> x_position, y_position, z_position;

    for (size_t i = 0; i<halo_cat->nObjects(); i++)
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

    const cbl::catalogue::Catalogue test_random_catalogue_post {cbl::catalogue::_box_, halo_catalogue, N_R};


    string w_file_post = "xi_RIVA_post.dat";
    cbl::twopt::TwoPointCorrelation1D_monopole TPCF_catalogue_post {temp_initial_catalogue_object, test_random_catalogue_post, cbl::_logarithmic_, rMin, rMax, nbins, shift};

    TPCF_catalogue_post.measure(cbl::twopt::_Poisson_, dir_output);

    TPCF_catalogue_post.write(dir_output, w_file_post);


    exit(1000);
  */
  } // End of RIVA method

  else ErrorCBL ("Error in cbl::catalogue::Catalogue::Catalogue() in VoidCatalogue.cpp: algorithm type is not correct!");

  // --------------------------------------------- //
  // ---------------- Second Step ---------------- //
  // --------------------------------------------- //
    
  auto start_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(initial_catalogue)));
  auto total_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(total_catalogue)));
  auto displ_cat = make_shared<cbl::catalogue::Catalogue> (cbl::catalogue::Catalogue(move(displacement_catalogue)));

  coutCBL << "Number of all paired halos: " << total_catalogue.nObjects() << endl << endl;
  

// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
//   
//   // * * * Two Point Correlation Function * * * //
// 
//   const double N_R = 1.; // random/data ratio
//   
//   const cbl::catalogue::Catalogue test_random_catalogue {cbl::catalogue::_box_, total_catalogue, N_R};
// 
//   const double rMin  = 1e3; // minimum separation 
//   const double rMax  = 5e4; // maximum separation 
//   const int nbins    = 20;  // number of bins
//   const double shift = 0.5; // spatial shift used to set the bin centre
//   
//   string w_file = "xi.dat";
//   
//   cbl::twopt::TwoPointCorrelation1D_monopole TPCF_total_catalogue {total_catalogue, test_random_catalogue, cbl::_logarithmic_, rMin, rMax, nbins, shift};
//       
//   TPCF_total_catalogue.measure(cbl::twopt::_Poisson_, dir_output);
//   
//   TPCF_total_catalogue.write(dir_output, w_file);
// 
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

// ^ ^ ^ Uncomment the previous part to compute the Two Point Correlation Function of the final displacement with respect to a generic random field ^ ^ ^ //
  
  string name;
//   string output;
  
//   output = "coord_LCDM.dat";
  
  name = dir_output + "Divergence";
  name.append("_");
  name.append(output);
  ofstream fout_divergence(name.c_str());

  name = dir_output + "Subvoids";
  name.append("_");
  name.append(output);
  ofstream fout_subvoids(name.c_str());

  name = dir_output + "Minima";
  name.append("_");
  name.append(output);
  ofstream fout_local_minima(name.c_str());
  
  fout_divergence.precision(7);
  fout_subvoids.precision(7);
  fout_local_minima.precision(7);

  
  double min_x       = total_cat->catalogue::Catalogue::Min(Var::_X_), max_x = total_cat->catalogue::Catalogue::Max(Var::_X_);
  double min_y       = total_cat->catalogue::Catalogue::Min(Var::_Y_), max_y = total_cat->catalogue::Catalogue::Max(Var::_Y_);
  double min_z       = total_cat->catalogue::Catalogue::Min(Var::_Z_), max_z = total_cat->catalogue::Catalogue::Max(Var::_Z_);
  
  double delta_x     = max_x-min_x;
  double delta_y     = max_y-min_y;
  double delta_z     = max_z-min_z;
  
  double volume      = delta_x*delta_y*delta_z;
  
  double density     = halo_catalogue.nObjects()/(volume); 
  
  double fraction    = 1./3.;
  double step        = 1.5*pow(density,-fraction); 
  double substep     = 1.e-3*step;
  double substep_1   = 1./substep;
  double substep10_1 = 0.1/substep;
// // //   double step_d      = step*1.75;

  double sigma       = 1.5*step; 
  double sg          = 1./(2.*sigma*sigma);
  int    grid_points = int(volume/(step*step*step))+100;

  coutCBL << "Δx     = " << delta_x << endl;
  coutCBL << "Δy     = " << delta_y << endl;
  coutCBL << "Δz     = " << delta_z << endl;
  coutCBL << "Volume = " << volume  << endl;
  coutCBL << "n      = " << density << endl;
  coutCBL << "step   = " << step    << endl;
  coutCBL << "σ      = " << sigma   << endl;
  cout << endl;
  
  int mem_alloc_x    = int(delta_x/step); 
  int mem_alloc_y    = int(delta_y/step);
  int mem_alloc_z    = int(delta_z/step);
  int mem_alloc_A    = round(400*n_rec*1.5*1.5*1.5); //1.5^3 test <<--<<--<<--

  Tensor3Di Pn         (mem_alloc_x, vector<vector<int> >          (mem_alloc_y, vector<int>          (mem_alloc_z, 0)));
  Tensor3Dd Divergence (mem_alloc_x, vector<vector<double> >       (mem_alloc_y, vector<double>       (mem_alloc_z, 0.)));
  Tensor3Dd Gradient_x (mem_alloc_x, vector<vector<double> >       (mem_alloc_y, vector<double>       (mem_alloc_z, 0.)));
  Tensor3Dd Gradient_y (mem_alloc_x, vector<vector<double> >       (mem_alloc_y, vector<double>       (mem_alloc_z, 0.)));
  Tensor3Dd Gradient_z (mem_alloc_x, vector<vector<double> >       (mem_alloc_y, vector<double>       (mem_alloc_z, 0.)));
  Tensor4Di Pp         (mem_alloc_x, vector<vector<vector<int> > > (mem_alloc_y, vector<vector<int> > (mem_alloc_z, vector<int>(mem_alloc_A, 0))));
  
  // * * * Preparation for divergence field construction * * * //

  for (size_t o = 0; o<total_cat->nObjects(); o++) // Looking in which grid halo locates
  {
    int io = round((total_cat->xx(o)-min_x)/step);  
    int jo = round((total_cat->yy(o)-min_y)/step);
    int ko = round((total_cat->zz(o)-min_z)/step);
    
    for (int i=-3; i<=3; ++i)
    { 
      for (int j=-3; j<=3; ++j)
      {
        for (int k=-3; k<=3; ++k)
        {            
          if (i*i+j*j+k*k<=3.001*3.001 && io+i>0 && jo+j>0 && ko+k>0 && io+i<mem_alloc_x && jo+j<mem_alloc_y && ko+k<mem_alloc_z && Pn[io+i][jo+j][ko+k]<mem_alloc_A)
          {
            Pn[io+i][jo+j][ko+k]++;
            Pp[io+i][jo+j][ko+k][Pn[io+i][jo+j][ko+k]]=o;
          }
        }
      }
    }
  }
  
  // * * * Divergence field construction * * * //
  
  cout << endl;
  coutCBL << "* * * Divergence field construction * * *" << endl << endl;
  coutCBL << "Wait...\r"; cout.flush();
  int number=0;
  
  double temp_x=0., temp_y=0., temp_z=0.;
  double divergence=0., x_divergence=0., y_divergence=0., z_divergence=0.;

  vector<double> x_void(grid_points), y_void(grid_points), z_void(grid_points), divergence_void(grid_points);
  vector<double> x_locmin, y_locmin, z_locmin;

  temp_x=min_x;
  while (temp_x<max_x)
  {
    double counter = 10./(max_x-step);	
    double percentage = temp_x-min_x;
    coutCBL << "..." << int(percentage*counter)*10 << "% completed\r"; cout.flush(); 
    temp_x+=step;

    temp_y=min_y;
    
    while (temp_y<max_y)
    {
      temp_y+=step;
      
      temp_z=min_z;
      
      while (temp_z<max_z)
      {
        temp_z+=step;

        double dist_0  = 0., dist_x  = 0., dist_y  = 0., dist_z  = 0.;
        double xdist_0 = 0., xdist_x = 0., xdist_y = 0., xdist_z = 0.;
        double ydist_0 = 0., ydist_x = 0., ydist_y = 0., ydist_z = 0.;
        double zdist_0 = 0., zdist_x = 0., zdist_y = 0., zdist_z = 0.;
        
        double vx0     = 0., vy0     = 0., vz0     = 0., sum0    = 0.;                                    
        double vx_x    = 0., sum_x   = 0., vy_y    = 0., sum_y   = 0., vz_z  = 0., sum_z  = 0.;                                    
        double xvx0    = 0., xvy0    = 0., xvz0    = 0., xsum0   = 0.;                                    
        double xvx_x   = 0., xsum_x  = 0., xvy_y   = 0., xsum_y  = 0., xvz_z = 0., xsum_z = 0.;                                    
        double yvx0    = 0., yvy0    = 0., yvz0    = 0., ysum0   = 0.;                                    
        double yvx_x   = 0., ysum_x  = 0., yvy_y   = 0., ysum_y  = 0., yvz_z = 0., ysum_z = 0.;                                    
        double zvx0    = 0., zvy0    = 0., zvz0    = 0., zsum0   = 0.;                                    
        double zvx_x   = 0., zsum_x  = 0., zvy_y   = 0., zsum_y  = 0., zvz_z = 0., zsum_z = 0.;

        int io = round((temp_x-min_x)/step);
        int jo = round((temp_y-min_y)/step);
        int ko = round((temp_z-min_z)/step);
        
        if (io>0 && jo>0 && ko>0 && io<mem_alloc_x && jo<mem_alloc_y && ko<mem_alloc_z)
        {
          for (int no=1; no<=Pn[io][jo][ko]; ++no) // Cycle through all halos
          {
            int o = Pp[io][jo][ko][no]; 
            
            dist_0  = pow((total_cat->xx(o)-temp_x),2)                    + pow((total_cat->yy(o)-temp_y),2)                    + pow((total_cat->zz(o)-temp_z),2);
            dist_x  = pow((total_cat->xx(o)-temp_x-substep),2)            + pow((total_cat->yy(o)-temp_y),2)                    + pow((total_cat->zz(o)-temp_z),2);
            dist_y  = pow((total_cat->xx(o)-temp_x),2)                    + pow((total_cat->yy(o)-temp_y-substep),2)            + pow((total_cat->zz(o)-temp_z),2);
            dist_z  = pow((total_cat->xx(o)-temp_x),2)                    + pow((total_cat->yy(o)-temp_y),2)                    + pow((total_cat->zz(o)-temp_z-substep),2);

            xdist_0 = pow((total_cat->xx(o)-temp_x-10*substep),2)         + pow((total_cat->yy(o)-temp_y),2)                    + pow((total_cat->zz(o)-temp_z),2);
            xdist_x = pow((total_cat->xx(o)-temp_x-substep-10*substep),2) + pow((total_cat->yy(o)-temp_y),2)                    + pow((total_cat->zz(o)-temp_z),2);
            xdist_y = pow((total_cat->xx(o)-temp_x-10*substep),2)         + pow((total_cat->yy(o)-temp_y-substep),2)            + pow((total_cat->zz(o)-temp_z),2);
            xdist_z = pow((total_cat->xx(o)-temp_x-10*substep),2)         + pow((total_cat->yy(o)-temp_y),2)                    + pow((total_cat->zz(o)-temp_z-substep),2);

            ydist_0 = pow((total_cat->xx(o)-temp_x),2)                    + pow((total_cat->yy(o)-temp_y-10*substep),2)         + pow((total_cat->zz(o)-temp_z),2);
            ydist_x = pow((total_cat->xx(o)-temp_x-substep),2)            + pow((total_cat->yy(o)-temp_y-10*substep),2)         + pow((total_cat->zz(o)-temp_z),2);
            ydist_y = pow((total_cat->xx(o)-temp_x),2)                    + pow((total_cat->yy(o)-temp_y-substep-10*substep),2) + pow((total_cat->zz(o)-temp_z),2);
            ydist_z = pow((total_cat->xx(o)-temp_x),2)                    + pow((total_cat->yy(o)-temp_y-10*substep),2)         + pow((total_cat->zz(o)-temp_z-substep),2);

            zdist_0 = pow((total_cat->xx(o)-temp_x),2)                    + pow((total_cat->yy(o)-temp_y),2)                    + pow((total_cat->zz(o)-temp_z-10*substep),2);
            zdist_x = pow((total_cat->xx(o)-temp_x-substep),2)            + pow((total_cat->yy(o)-temp_y),2)                    + pow((total_cat->zz(o)-temp_z-10*substep),2);
            zdist_y = pow((total_cat->xx(o)-temp_x),2)                    + pow((total_cat->yy(o)-temp_y-substep),2)            + pow((total_cat->zz(o)-temp_z-10*substep),2);
            zdist_z = pow((total_cat->xx(o)-temp_x),2)                    + pow((total_cat->yy(o)-temp_y),2)                    + pow((total_cat->zz(o)-temp_z-substep-10*substep),2);
                            
            vx0    += displ_cat->xx(o)*exp(-dist_0*sg);
            vy0    += displ_cat->yy(o)*exp(-dist_0*sg);            
            vz0    += displ_cat->zz(o)*exp(-dist_0*sg);
            sum0   += exp(-dist_0*sg);
            
            vx_x   += displ_cat->xx(o)*exp(-dist_x*sg);
            sum_x  += exp(-dist_x*sg);

            vy_y   += displ_cat->yy(o)*exp(-dist_y*sg);
            sum_y  += exp(-dist_y*sg);

            vz_z   += displ_cat->zz(o)*exp(-dist_z*sg);
            sum_z  += exp(-dist_z*sg);

            xvx0   += displ_cat->xx(o)*exp(-xdist_0*sg);
            xvy0   += displ_cat->yy(o)*exp(-xdist_0*sg);
            xvz0   += displ_cat->zz(o)*exp(-xdist_0*sg);
            xsum0  += exp(-xdist_0*sg);
            
            xvx_x  += displ_cat->xx(o)*exp(-xdist_x*sg);
            xsum_x += exp(-xdist_x*sg);

            xvy_y  += displ_cat->yy(o)*exp(-xdist_y*sg);
            xsum_y += exp(-xdist_y*sg);

            xvz_z  += displ_cat->zz(o)*exp(-xdist_z*sg);
            xsum_z += exp(-xdist_z*sg); 

            yvx0   += displ_cat->xx(o)*exp(-ydist_0*sg);
            yvy0   += displ_cat->yy(o)*exp(-ydist_0*sg);
            yvz0   += displ_cat->zz(o)*exp(-ydist_0*sg);
            ysum0  += exp(-ydist_0*sg);
            
            yvx_x  += displ_cat->xx(o)*exp(-ydist_x*sg);
            ysum_x += exp(-ydist_x*sg);

            yvy_y  += displ_cat->yy(o)*exp(-ydist_y*sg);
            ysum_y += exp(-ydist_y*sg);

            yvz_z  += displ_cat->zz(o)*exp(-ydist_z*sg);
            ysum_z += exp(-ydist_z*sg);

            zvx0   += displ_cat->xx(o)*exp(-zdist_0*sg);
            zvy0   += displ_cat->yy(o)*exp(-zdist_0*sg);
            zvz0   += displ_cat->zz(o)*exp(-zdist_0*sg);
            zsum0  += exp(-zdist_0*sg);
            
            zvx_x  += displ_cat->xx(o)*exp(-zdist_x*sg);
            zsum_x += exp(-zdist_x*sg);

            zvy_y  += displ_cat->yy(o)*exp(-zdist_y*sg);
            zsum_y += exp(-zdist_y*sg);

            zvz_z  += displ_cat->zz(o)*exp(-zdist_z*sg);
            zsum_z += exp(-zdist_z*sg);            
          }
          
          double sum0_1  = 1./sum0;
          
          vx0          = vx0*sum0_1;
          vy0          = vy0*sum0_1;
          vz0          = vz0*sum0_1;

          vx_x         = vx_x/sum_x;
          vy_y         = vy_y/sum_y;
          vz_z         = vz_z/sum_z;
          
          divergence   = (vx_x-vx0 + vy_y-vy0 + vz_z-vz0)*substep_1;
                                                                        
          xvx0         = xvx0/xsum0;
          xvy0         = xvy0/xsum0;
          xvz0         = xvz0/xsum0;

          xvx_x        = xvx_x/xsum_x;
          xvy_y        = xvy_y/xsum_y;
          xvz_z        = xvz_z/xsum_z;
          
          x_divergence = (xvx_x-xvx0 + xvy_y-xvy0 + xvz_z-xvz0)*substep_1;
          
          yvx0         = yvx0/ysum0;
          yvy0         = yvy0/ysum0;
          yvz0         = yvz0/ysum0;

          yvx_x        = yvx_x/ysum_x;
          yvy_y        = yvy_y/ysum_y;
          yvz_z        = yvz_z/ysum_z;

          y_divergence = (yvx_x-yvx0 + yvy_y-yvy0 + yvz_z-yvz0)*substep_1;                                        
          
          zvx0         = zvx0/zsum0;
          zvy0         = zvy0/zsum0;
          zvz0         = zvz0/zsum0;

          zvx_x        = zvx_x/zsum_x;
          zvy_y        = zvy_y/zsum_y;
          zvz_z        = zvz_z/zsum_z;
          
          z_divergence = (zvx_x-zvx0 + zvy_y-zvy0 + zvz_z-zvz0)*substep_1;

          Gradient_x[io][jo][ko] = (x_divergence-divergence)*substep10_1;
          Gradient_y[io][jo][ko] = (y_divergence-divergence)*substep10_1;
          Gradient_z[io][jo][ko] = (z_divergence-divergence)*substep10_1;

          Divergence[io][jo][ko] = divergence;

          if (divergence<0)
          {
            number++;
            
            x_void[number]          = temp_x;
            y_void[number]          = temp_y;
            z_void[number]          = temp_z;
            
            divergence_void[number] = divergence;
          }
          
          fout_divergence << temp_x << " " << temp_y << " " << temp_z << " " << divergence << endl; // Divergence field for each grid point;
        }
      }
    }
  }
  
  fout_divergence.clear(); fout_divergence.close();

  coutCBL << "Number of cells with negative divergence: " << number << endl << endl;

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
                if (Divergence[io+i][jo+j][ko+k]>Divergence[io][jo][ko] && Divergence[io][jo][ko]<-1e-10) index++; // Condition on minimum value
              }   
            }
          }
        }
        
        if (index==26)
        {
          x_locmin.push_back(io*step+min_x);
          y_locmin.push_back(jo*step+min_y);
          z_locmin.push_back(ko*step+min_z);
          
          fout_local_minima << io*step+min_x << " " << jo*step+min_y << " " << ko*step+min_z << endl; // local minima

        }
      }
    }
  }
  
  fout_local_minima.clear(); fout_local_minima.close();

  
//   coutCBL << "Number of local minima: " << x_locmin.size() << endl;
//   cout << endl;

//   coutCBL << "* * * Identification of subvoids * * *" << endl;
//   cout << endl;

  int control1 = 0;
  int control2 = 0;
  int control3 = 0;
  int control4 = 0;
  int control5 = 0;
  
  int subvoids_number = 0;
  vector <int> number_host_subvoids;
  vector <double> x_subvoid, y_subvoid, z_subvoid, div_subvoid;
  vector <int> id_subvoid;
  for (int temp_number=1; temp_number<=number; ++temp_number)
  {
    temp_x = x_void[temp_number];  // Initial coordinates
    temp_y = y_void[temp_number];
    temp_z = z_void[temp_number];
        
    int rejected = 0;
    int tumbler  = 0;
    
    while (tumbler==0)
    {
      int io = round((temp_x-min_x)/step); // Current coordinates
      int jo = round((temp_y-min_y)/step);
      int ko = round((temp_z-min_z)/step);

      if (Divergence[io][jo][ko]>0) // Border of subvoid
      {
        tumbler = 1;
        control1++;
      }

      if (rejected>1000) // Loop
      {
        tumbler=1;
        control2++;
      }

      for (int ID=0; ID<(int)x_locmin.size(); ++ID)
      {
        if (fabs(temp_x-x_locmin[ID])<1.1*step && fabs(temp_y-y_locmin[ID])<1.1*step && fabs(temp_z-z_locmin[ID])<1.1*step)
        {
          tumbler = 1;
          control3++; 
          fout_subvoids << x_locmin[ID]<< " " << y_locmin[ID]<< " " << z_locmin[ID]<< " " << temp_x << " " << temp_y << " " << temp_z << " " << x_void[temp_number] << " " << y_void[temp_number] << " " << z_void[temp_number] << " " << ID << " " << Divergence[io][jo][ko] << " " << rejected << endl;
          subvoids_number++;
          number_host_subvoids.push_back(ID);
          x_subvoid.push_back(x_locmin[ID]);
          y_subvoid.push_back(y_locmin[ID]);
          z_subvoid.push_back(z_locmin[ID]);
          id_subvoid.push_back(ID);
          div_subvoid.push_back(Divergence[io][jo][ko]);
        }
      }

      int measured = 0;
      if (temp_x>min_x && temp_x<max_x && temp_y>min_y && temp_y<max_y && temp_z>min_z && temp_z<max_z) measured = 1; 
      if (measured==0) // Border of box
      {
        tumbler = 1;
        control4++;
      }
                                 
      double v0_tot   = sqrt(Gradient_x[io][jo][ko]*Gradient_x[io][jo][ko]+Gradient_y[io][jo][ko]*Gradient_y[io][jo][ko]+Gradient_z[io][jo][ko]*Gradient_z[io][jo][ko]);
      double v0_tot_1 = 1/v0_tot;

      temp_x = temp_x-step*Gradient_x[io][jo][ko]*v0_tot_1;
      temp_y = temp_y-step*Gradient_y[io][jo][ko]*v0_tot_1;
      temp_z = temp_z-step*Gradient_z[io][jo][ko]*v0_tot_1;
      
      io = round((temp_x-min_x)/step);
      jo = round((temp_y-min_y)/step);
      ko = round((temp_z-min_z)/step);    
      
      if (io<=0 || jo<=0 || ko<=0 || io>=mem_alloc_x || jo>=mem_alloc_y || ko>=mem_alloc_z)
      {
        tumbler = 1;
        control5++;
      }

      rejected++;
    }
  }
 
  if (2*control3<control1+control2+control3+control4+control5) ErrorCBL("Error in cbl::catalogue::Catalogue::Catalogue() in VoidCatalogue.cpp: too sparse sample or too small/large range of values");
   
  fout_subvoids.clear(); fout_subvoids.close();

  // -------------------------------------------- //
  // ---------------- Third Step ---------------- //
  // -------------------------------------------- //
        
//   coutCBL << "Number of subvoids: " << subvoids_number << endl;
  cout << endl;
  coutCBL << "* * * Identification of voids * * *" << endl;
  
  if (subvoids_number==0) ErrorCBL ("Error in cbl::catalogue::Catalogue::Catalogue() in VoidCatalogue.cpp: no subvoids were found");
  
// // //   int max_ID = *max_element(number_host_subvoids.begin(), number_host_subvoids.end());
// // //   
// // //   double radius_lim = 3*step*pow(3./(4.*M_PI),fraction);
// // //   coutCBL << "Minimum effective radius: " << radius_lim << endl;
    
  name = dir_output+"Voids_";
  name.append(output);
  ofstream fout_voids(name.c_str());

  vector <int> id_unique_subvoid = id_subvoid;
  std::sort (id_unique_subvoid.begin(),id_unique_subvoid.end());  

  auto last = std::unique(id_unique_subvoid.begin(), id_unique_subvoid.end());

  id_unique_subvoid.erase(last, id_unique_subvoid.end()); 

  double cell_volume = step*step*step;
  vector <double> radius_void, counter_subvoid(id_unique_subvoid.size());
  vector <double> x_c_void(id_unique_subvoid.size()), y_c_void(id_unique_subvoid.size()), z_c_void(id_unique_subvoid.size()), div_c_void(id_unique_subvoid.size());
  
  for (int j=0; j<(int)id_subvoid.size();j++)
  {
    for (int i=0; i<(int)id_unique_subvoid.size(); i++)
    {
      if (id_subvoid[j]==id_unique_subvoid[i] && div_subvoid[j]<0.)
      {   
        counter_subvoid[i]=counter_subvoid[i]+1;    // counter of cells constituing a void
        
        x_c_void[i]=x_subvoid[j];                   // x-centre coordinate
        y_c_void[i]=y_subvoid[j];                   // y-centre coordinate
        z_c_void[i]=z_subvoid[j];                   // z-centre coordinate
        
        div_c_void[i]=div_subvoid[j];               // divergence value in the centre of void
      }   
    }
  }
  
  for (int i=0; i<(int)counter_subvoid.size(); i++) 
    
    radius_void.push_back(pow(((3/(4*M_PI))*cell_volume*counter_subvoid[i]), fraction)); // compute the radius of the i-th void

  for (int i=0; i<(int)x_c_void.size() ;i++)
  {
    fout_voids << x_c_void[i] << " " << y_c_void[i] << " " << z_c_void[i] << " " << radius_void[i] << " " << div_c_void[i] << endl;
  }
  
  fout_voids.clear(); fout_voids.close();

  cout << endl;
  coutCBL << "Number of voids: " << x_locmin.size() << endl << endl;

  
  // ---------------------------------------------------------- //
  // ---------------- Computing execution time ---------------- //
  // ---------------------------------------------------------- //
  
  float seconds = float( clock () - begin_time ) / CLOCKS_PER_SEC;
  
  cout << endl;
  coutCBL << "Time spent to compute: " << seconds      << " seconds " << endl ;
  coutCBL << "Time spent to compute: " << seconds/60   << " minutes " << endl ;
  coutCBL << "Time spent to compute: " << seconds/3600 << " hours "   << endl ;
  cout << endl;
}

/// @endcond

/////////////////////////////////// Tommaso Ronconi //////////////////////////////////////////

cbl::catalogue::Catalogue::Catalogue (const std::shared_ptr<Catalogue> input_voidCatalogue, const std::vector<bool> clean, const std::vector<double> delta_r, const double threshold, const double statistical_relevance, bool rescale, const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const double ratio, const bool checkoverlap, const Var ol_criterion)
{

  auto catalogue = input_voidCatalogue;
  clock_t begin_time = clock();
    
  // ---------------------------------------------------- //
  // ---------------- Cleaning Procedure ---------------- //
  // ---------------------------------------------------- //

  if (clean.size() != 3) ErrorCBL("Error in cbl::catalogue::Catalogue::Catalogue() in Catalogue.cpp: wrong vector size!");  
  if (clean[0] || clean[1] || clean[2]) {
    vector<int> counter(clean.size(), 0);
    vector<bool> remove(catalogue->nObjects(), false);
    coutCBL << "Input void catalogue cleaning procedure started ..." << endl;
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
    
    coutCBL << "############ Removed Voids ###########" << "\n" <<
      "\t r_min - r_max criterion : " << counter[0] << "\n" <<
      "\t central density too high: " << counter[1] << "\n" <<
      "\t statistically irrelevant: " << counter[2] << "\n" <<
      "\t total removed: " << counter[0]+counter[1]+counter[2] << endl;
  }
  coutCBL << "\n Voids in the Catalogue: " << catalogue->nObjects() << endl;
  
  float cleaning_time = float( clock() - begin_time ) / CLOCKS_PER_SEC;
  coutCBL << "Time spent by the cleaning procedure: " << cleaning_time << " seconds \n" << endl;
  
  // ---------------------------------------------------- //
  // ----------------- Radius Rescaling ----------------- //
  // ---------------------------------------------------- //
  
  if (rescale) {
    coutCBL << "Rescaling radii ... " << endl;
    
    tracers_catalogue->compute_catalogueProperties();
    double density = tracers_catalogue->numdensity();
    
    //vector to memorize which element of the catalogue has to be removed at the end of the procedure:
    vector<bool> remove(catalogue->nObjects(), false);
    
    //counter for regions without any tracer:
    int void_voids = 0;
    
    //counter for negative radius voids:
    int negative = 0;
    
    for (size_t j = 0; j<catalogue->nObjects(); j++) {
      double value = (3.*catalogue->radius(j) < delta_r[1]) ? 3.*catalogue->radius(j) : delta_r[1];
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
	else while (kk/(volume_sphere(distances[kk-1])*density) < threshold && kk < (int) distances.size()) kk++; // or you expand
	  
	// linear interpolation:
	double new_radius = interpolated(threshold,
					 {kk/(volume_sphere(distances[kk-1])*density), (kk+1)/(volume_sphere(distances[kk])*density)},
					 {distances[kk-1], distances[kk]}, "Linear"); // gsl function
	  
	  
	if (new_radius < 0) {
	  negative++;
	  remove[j] = true;
	}
	else catalogue->set_var(j, Var::_Radius_, abs(new_radius));
      }
      else {
	void_voids ++;
	remove[j] = true;
      }
    }//for
    
    catalogue->remove_objects(remove);
    coutCBL << "\t Empty voids removed: " << void_voids << endl;
    coutCBL << "\t Negative voids removed: " << negative << "\n" << endl;
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
    coutCBL << "\t Removed voids out of range ["+conv(delta_r[0],par::fDP2)+","+conv(delta_r[1],par::fDP2)+"] : " << outofrange << endl;
    
    //compute new central density and density contrast:
    catalogue->compute_centralDensity(tracers_catalogue, ChM, density, ratio);
    catalogue->compute_densityContrast(tracers_catalogue, ChM, ratio);
    
  }//rescale part
  
  coutCBL << "\n Voids in the Catalogue: " << catalogue->nObjects() << endl;

  float rescaling_time = float( clock() - begin_time ) / CLOCKS_PER_SEC;
  rescaling_time = rescaling_time - cleaning_time;
  coutCBL << "Time spent by the rescaling procedure: " << rescaling_time << " seconds \n" << endl;
    
  // ---------------------------------------------------- //
  // ------------------ Overlap Check ------------------- //
  // ---------------------------------------------------- //
  
  if (checkoverlap) {
    coutCBL << "Checking for overlapping voids ..." << endl;

    if (ol_criterion == Var::_CentralDensity_) catalogue->sort(ol_criterion, false);
    else if (ol_criterion == Var::_DensityContrast_) catalogue->sort(ol_criterion, true);
    else ErrorCBL("Error in cbl::catalogue::Catalogue::Catalogue() in VoidCatalogue.cpp: allowed overlap criteria are '_CentralDensity_' or '_DensityContrast_' .");
    
    vector<bool> remove(catalogue->nObjects(), false);
    
    coutCBL << "\t Generating chain-mesh for void centres ..." << endl;
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
    coutCBL << "\t Voids removed to avoid overlap: " << overlap_removed << endl;
  }//overlap check
  coutCBL << "\n Voids in the Catalogue: " << catalogue->nObjects() << endl;

  float olchecking_time = float( clock() - begin_time ) / CLOCKS_PER_SEC;
  olchecking_time = olchecking_time - rescaling_time - cleaning_time;
  coutCBL << "Time spent by the overlap-checking procedure: " << olchecking_time << " seconds" << endl;
  coutCBL << "\n Total time spent: " << float( clock() - begin_time ) / CLOCKS_PER_SEC << " seconds \n" << endl;
  
  m_object = catalogue->sample();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void cbl::catalogue::Catalogue::compute_centralDensity (const shared_ptr<Catalogue> tracers_catalogue, cbl::chainmesh::ChainMesh3D ChM, const double density, const double ratio) {
  ///
  //coutCBL << "check1" << endl;
  ///
  //vector to memorize which element of the catalogue has to be removed at the end of the procedure:
  vector<bool> remove(m_object.size(), false);
  ///
  //coutCBL << "check2" << endl;
  ///
  //counter for regions without any tracer:
  int void_voids = 0;
  
  for (size_t j = 0/*27971*/; j<m_object.size(); j++) {
    ///
    //vector<double> coords = m_object[j]->coords();
    //coutCBL << m_object[j]->radius() << "  " << coords[0] << "  " << coords[1] << "  " << coords[2] << endl;
    ///
    ChM.get_searching_region(m_object[j]->radius());
    vector<long> close = ChM.close_objects(m_object[j]->coords());
    ///
    //coutCBL << close.size() << endl;
    ///
    //compute distances between the void and the surrounding particles
      vector<double> distances;
      for (auto&& k : close) {
	double distance = cbl::catalogue::Catalogue::distance(j, tracers_catalogue->catalogue_object(k));
	if (distance < m_object[j]->radius()) distances.emplace_back(distance);
      }
      ///
      //coutCBL << distances.size() << endl;
      ///
      //FIND CENTRAL DENSITY
      if (distances.size() > 0) {
	std::sort (distances.begin(), distances.end());
	int NN = 0;
	while (distances[NN]<ratio*m_object[j]->radius()) NN++;
	///
	//coutCBL << (NN/cbl::volume_sphere(distances[NN]))/density << endl;
	///
	m_object[j]->set_centralDensity((NN/cbl::volume_sphere(distances[NN]))/density);
      }
      else {
	void_voids ++;
	remove[j] = true;
      }
      //coutCBL << j << "/" << m_object.size() << ": " << m_object << endl;
      /*
      if (distances.size() > 3) {
	std::sort (distances.begin(), distances.end());
	int NN = (int(ratio*distances.size())>3) ? int(ratio*distances.size())-1 : 3;
	m_object[j]->set_centralDensity((NN/cbl::volume_sphere(distances[NN]))/density);
      }   
      else {
	void_voids ++;
	remove[j] = true;
      } 
      */
      }//for
  ///
  //coutCBL << "check3" << endl;
  ///
  for (size_t j = 0; j<m_object.size(); j++)
    if (remove[j]) m_object.erase(m_object.begin()+j);
  ///
  //coutCBL << "check4" << endl;
  ///
  coutCBL << "I removed " << void_voids << " voids in calculating the central density!" << endl;
  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void cbl::catalogue::Catalogue::compute_densityContrast (const shared_ptr<Catalogue> tracers_catalogue, cbl::chainmesh::ChainMesh3D ChM, const double ratio) {
  
  //vector to memorize which element of the catalogue has to be removed at the end of the procedure:
  vector<bool> remove(m_object.size(), false);
    
  //counter for regions without any tracer:
  int void_voids = 0;
  int cloud_in_void = 0;
  
  for (size_t j = 0; j<m_object.size(); j++) {
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
	while (distances[NN]<ratio*m_object[j]->radius()) NN++;
	/*
	  if (distances.size() > 3) {
	  std::sort (distances.begin(), distances.end());
	  int NN = (int(ratio*distances.size())>3) ? int(ratio*distances.size())-1 : 3;*/
	if (NN > 0) {
	  double delta_in = NN/cbl::volume_sphere(distances[NN]);
	  double delta_out = distances.size()/cbl::volume_sphere(distances[distances.size()-1]);
	  if (delta_out/delta_in < 1.) {
	    cloud_in_void ++;
	    remove[j] = true;
	  }
	  else m_object[j]->set_densityContrast(delta_out/delta_in);
	}
	else m_object[j]->set_densityContrast(1.);
      }    
      else {
	void_voids ++;
	remove[j] = true;
      } 
  }//for
  
  for (size_t j = 0; j<m_object.size(); j++)
    if (remove[j]) m_object.erase(m_object.begin()+j);

  coutCBL << "Cloud-in-void: " << cloud_in_void << endl;
  coutCBL << "I removed " << void_voids+cloud_in_void << " voids in calculating the density contrast!" << endl;
  
}
