/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Carlo Cannarozzo     *
 *  federico.marulli3@unibo.it    carlo.cannarozzo@studio.unibo.it  *
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
 *  @author Federico Marulli, Carlo Cannarozzo
 *
 *  @author federico.marulli3@unbo.it, carlo.cannarozzo@studio.unibo.it
 */

#include "Func.h"
#include "Catalogue.h"
#include "ChainMesh_Catalogue.h"
#include "Object.h"
#include <time.h>
using namespace cosmobl;


// ============================================================================

typedef vector<vector <double> > Matrix;

/// @cond extvoid

cosmobl::catalogue::Catalogue::Catalogue (const VoidAlgorithm algorithm, const Catalogue halo_catalogue, const vector<string> file, const double nSub, const int n_rnd, const string mode, const double rmax, const int cellsize)
{
  // -------------------------------------------- //
  // ---------------- First Step ---------------- //
  // -------------------------------------------- //
  
  const clock_t begin_time = clock();
  
  cosmology::Cosmology cosm;
  
  Catalogue total_catalogue;
  vector<shared_ptr<cosmobl::catalogue::Object>> object;

  
  if (algorithm==_LaZeVo_)
  {
    cout << endl << par::col_green << "    > LaZeVo algorithm" << par::col_default << endl;

    Catalogue rand_catalogue; // 'Catalogue' constructor for random halo objects

    string rand_file;
    
    for (size_t dd=0; dd<file.size(); ++dd) rand_file = file[dd];

    if (rand_file!="not_provided") 
    {
      Catalogue rand_tempcatalogue {_RandomObject_, _comovingCoordinates_, {rand_file}, 1, 2, 3, -1, -1, nSub};
      cout << "    > The number of random halos in catalogue is " << rand_tempcatalogue.nObjects() << endl << endl;
      rand_catalogue = rand_tempcatalogue;
    }
    
    cout << endl << " * * * Looking for closest pairs particles-random particles and calculating displacement field * * * " << endl << endl;

    for (int rndd=1; rndd<=n_rnd; ++rndd) // 'for' cycle - Cycling to the number of realizations 'n_rnd' set up in the 'input_param.ini' file 
    {
      cout << "    > " << rndd << " of " << n_rnd << " random realizations " << endl;

      if (rand_file=="not_provided" && (mode=="non_periodic"||mode=="periodic"))
      {
        Catalogue rand_tempcatalogue {_createRandom_box_, halo_catalogue, 1., 10, cosm, false, 10., (int) time(NULL)}; // 2rand to pair all particles

        rand_catalogue = rand_tempcatalogue;      
      }

      stringstream ss;
      ss << rndd;
      string str = ss.str();
      string output = ss.str();
      string name = "../output/Displacement";
      name.append("_");
      name.append("_");
      name.append(output);
      ofstream fout_displ(name.c_str());

      auto halo_cat = make_shared<cosmobl::catalogue::Catalogue> (cosmobl::catalogue::Catalogue(move(halo_catalogue)));
      auto rand_cat = make_shared<cosmobl::catalogue::Catalogue> (cosmobl::catalogue::Catalogue(move(rand_catalogue)));
      
      for (size_t i = 0; i<halo_cat->nObjects(); i++)
      {
        cosmobl::chainmesh::ChainMesh_Catalogue ChM (cellsize, rand_cat, rmax);

        ChM.get_searching_region(rmax);
        vector<long> close_objects = ChM.close_objects(halo_cat->coordinate(i));
            
        double counter = 100./halo_cat->nObjects();	
        cout << " ..." << int(i*counter) << "% completed \r"; cout.flush(); 
        
        if (close_objects.size()>0)
        {
          vector<double> distance;
          vector<double> vect_i;
          vector<double> vect_k;
          vect_i.push_back(i);

          for (auto&& k : close_objects)
          {
            distance.push_back(cosmobl::Euclidean_distance(halo_cat->xx(i), rand_cat->xx(k), halo_cat->yy(i), rand_cat->yy(k), halo_cat->zz(i), rand_cat->zz(k)));
            vect_k.push_back(k);
            vect_i.push_back(i);
          }

          vector<double>::iterator p1 = distance.begin(), p2 = vect_i.begin(), p3 = vect_k.begin();

          sort_3vectors(p1, p2, p3, distance.size());

          fout_displ.precision(10);
          fout_displ<<halo_cat->xx(vect_i[0])<<" "<<halo_cat->yy(vect_i[0])<<" "<<halo_cat->zz(vect_i[0])<<" "<<rand_cat->xx(vect_k[0])<<" "<<rand_cat->yy(vect_k[0])<<" "<<rand_cat->zz(vect_k[0])<<" "<<distance[0]<<endl;
  
          cosmobl::comovingCoordinates coordzz = {rand_cat->xx(vect_k[0]), rand_cat->yy(vect_k[0]), rand_cat->zz(vect_k[0])};
          auto halo = make_shared<cosmobl::catalogue::Halo>(coordzz);
          object.emplace_back(halo);
                
          rand_cat->remove_object(vect_k[0]);
        }
      }
      cout << " Number of unpaired halos: " << rand_cat->nObjects() << endl;
      cout << " Percentage of unpaired halos: " << (float)rand_cat->nObjects()/halo_cat->nObjects() * 100 << "% " << endl;
      cout << endl;
    }
  }
  
  else if (algorithm==_RIVA_) {cout << endl << par::col_green << " RIVA algorithm is coming soon... " << par::col_default << endl; exit(1);}

  else ErrorCBL (" Algorithm type is not correct! ");

  // --------------------------------------------- //
  // ---------------- Second Step ---------------- //
  // --------------------------------------------- //
  
  total_catalogue.add_objects(object);
  cout << " Number of all paired halos " << total_catalogue.nObjects() << endl;

  auto total_cat = make_shared<cosmobl::catalogue::Catalogue> (cosmobl::catalogue::Catalogue(move(total_catalogue)));
  
  string name = /*dir_output+*/"../output/Divergence";
//   name.append(output);
  ofstream fout_divergence(name.c_str());

  name = /*dir_output+*/"../output/Subvoids";
//   name.append(output);
  ofstream fout_subvoids(name.c_str());

  fout_divergence.precision(7);
  fout_subvoids.precision(7);
  
  double min_x = total_cat->catalogue::Catalogue::Min(_X_), max_x = total_cat->catalogue::Catalogue::Max(_X_);
  double min_y = total_cat->catalogue::Catalogue::Min(_Y_), max_y = total_cat->catalogue::Catalogue::Max(_Y_);
  double min_z = total_cat->catalogue::Catalogue::Min(_Z_), max_z = total_cat->catalogue::Catalogue::Max(_Z_);
  
  double delta_x = max_x-min_x;
  double delta_y = max_y-min_y;
  double delta_z = max_z-min_z;
  
  double volume = delta_x*delta_y*delta_z;
  
  double density = total_cat->nObjects()/(n_rnd*volume);
   
  double step            = 1.5*pow(density,-0.333333); 
  double substep         = 1.e-3*step;
  double substep_1       = 1./substep;
  double substep10_1     = 0.1/substep;

  double sigma           = 1.5*step; 
  double sg              = 1./(2.*sigma*sigma);
  int grid_points        = int(volume/(step*step*step))+100;

  cout << endl;
  cout << " Δx           = " << delta_x     << endl;
  cout << " Δy           = " << delta_y     << endl;
  cout << " Δz           = " << delta_z     << endl;
  cout << endl;
  cout << " Volume       = " << volume      << endl;
  cout << endl;
  cout << " n            = " << density     << endl;
  cout << endl;
  cout << " Step         = " << step        << endl;
  cout << " Substep      = " << substep     << endl;
  cout << " Substep^-1   = " << substep_1   << endl;
  cout << " Substep10^-1 = " << substep10_1 << endl;
  cout << endl;
  cout << " σ            = " << sigma       << endl;
  cout << " Sg           = " << sg          << endl;
  cout << " Grid points  = " << grid_points << endl;
  cout << endl;

  int mem_alloc_x        = int(delta_x/step); 
  int mem_alloc_y        = int(delta_y/step);
  int mem_alloc_z        = int(delta_z/step);
  int memory_allocation  = round(400*n_rnd*1.5*1.5*1.5); //1.5^3 test <<--<<--<<-- by Andrii
  
  cout << " mem_alloc_x  = " << mem_alloc_x << endl;
  cout << " mem_alloc_y  = " << mem_alloc_y << endl;
  cout << " mem_alloc_z  = " << mem_alloc_z << endl;
  cout << " memory_all   = " << memory_allocation << endl;  
  cout << endl;
  
  Matrix Divergence(mem_alloc_x, vector<double> (mem_alloc_z,0.));
  
  cout << " --- " << endl;
  print(Divergence);
  cout << " --- " << endl;

  
  // -------------------------------------------- //
  // ---------------- Third Step ---------------- //
  // -------------------------------------------- //

  
  float diff = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
  
  cout << " Time spent to count the pairs: " << diff      << " seconds " << endl ;
  cout << " Time spent to count the pairs: " << diff/60   << " minutes " << endl ;
  cout << " Time spent to count the pairs: " << diff/3600 << " hours "   << endl ;
  cout << endl;

}

/// @endcond
