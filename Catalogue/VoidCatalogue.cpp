/********************************************************************
 *  Copyright (C) 2024 by Simone Sartori, Federico Marulli,         *
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
 *  @author Simone Sartori, Sofia Contarini, Federico Marulli, Carlo 
 *  Cannarozzo, Tommaso Ronconi
 *
 *  @author simone.sartori5@studio.unibo.it,
 *  sofia.contarini3@unibo.it, federico.marulli3@unibo.it,
 *  carlo.cannarozzo@studio.unibo.it, tommaso.ronconi@studio.unibo.it
 */

#include <omp.h>
#include <chrono>
#include <optional>
#include <thread>
#include <functional>

#include "Func.h"
#include "Catalogue.h"
#include "ChainMesh_Catalogue.h"
#include "Object.h"
#include "Catalogue.h"
#include "Data1D.h"
#include "Posterior.h"
#include "TwoPointCorrelation.h"
#include "TwoPointCorrelation1D_monopole.h"
#include "MeshGrid.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace glob;


// ============================================================================


cbl::catalogue::Catalogue::Catalogue (const VoidAlgorithm algorithm, cbl::catalogue::Catalogue tracer_catalogue, std::vector<cbl::catalogue::Catalogue> random_catalogue, const std::string dir_output, const std::string output, double sigma, const double cellsize, const int n_rec, const double step_size, const double threshold, const std::vector<bool> print, unsigned int seed)
{
  // -------------------------------------------- //
  // ---------------- First Step ---------------- //
  // -------------------------------------------- //

  cout << endl; coutCBL << "> > > > > >>>>>>>>>>   " << par::col_blue << "BACK-IN-TIME VOID FINDER" << par::col_default << "   <<<<<<<<<< < < < < <" << endl << endl;

  const double start_time = omp_get_wtime();
  const unsigned int num_objects = tracer_catalogue.nObjects();  
  const double mps = tracer_catalogue.mps();

  vector<vector<double>> coordinates(3, vector<double>(num_objects));

  for (unsigned int i=0; i<num_objects; i++) {
    coordinates[0][i] = tracer_catalogue.xx(i);
    coordinates[1][i] = tracer_catalogue.yy(i);
    coordinates[2][i] = tracer_catalogue.zz(i);
  }
  tracer_catalogue = {};
  
  double minTrCatX = *min_element(coordinates[0].begin(), coordinates[0].end()), maxTrCatX = *max_element(coordinates[0].begin(), coordinates[0].end());
  double minTrCatY = *min_element(coordinates[1].begin(), coordinates[1].end()), maxTrCatY = *max_element(coordinates[1].begin(), coordinates[1].end());
  double minTrCatZ = *min_element(coordinates[2].begin(), coordinates[2].end()), maxTrCatZ = *max_element(coordinates[2].begin(), coordinates[2].end());

  vector<vector<double>> random_final(3, vector<double>(num_objects));
  vector<vector<vector<double>>> displacement_temp(n_rec, vector<vector<double>>(num_objects, vector<double>(3, 0.)));
  vector<double> Displacement_finalStd(num_objects, 0.);
  vector<double> Displacement_finalAngStd(num_objects, 0.);
  vector<vector<double>> displacement(num_objects, vector<double>(3, 0.));
  vector<vector<bool>> displacement_mask(n_rec, vector<bool>(num_objects, true));
  vector<unsigned int> n_rec_eff(num_objects, n_rec);

  meshgrid::MeshGrid MeshGrid_tracers;
  coutCBL << "* * * Mesh Grid setting * * *" << endl;
  vector<vector<double>> lims(3, vector<double>(2));
  lims[0][0] = minTrCatX-0.1*mps; 
  lims[0][1] = maxTrCatX+0.1*mps; 
  lims[1][0] = minTrCatY-0.1*mps; 
  lims[1][1] = maxTrCatY+0.1*mps; 
  lims[2][0] = minTrCatZ-0.1*mps;  
  lims[2][1] = maxTrCatZ+0.1*mps; 
  MeshGrid_tracers = meshgrid::MeshGrid(coordinates[0], coordinates[1], coordinates[2], cellsize, lims);
  cout << endl;

  if (random_catalogue.size() == 0) {
    coutCBL << "* * * "<< n_rec << " random catalogue(s) will be generated with " << num_objects << " objects * * *" << endl;
  }  
  else {
    unsigned int randCounter = 0;
    for (int rec=0; rec<n_rec; rec++) {
      if (random_catalogue[rec].nObjects() == 0) {   
        randCounter++;
      }
    }     
    if (randCounter>0) {
      coutCBL << "* * * "<< n_rec-randCounter << " random catalogue(s) will be read from file(s)" << endl; 
      coutCBL << "* * * "<< randCounter << " random catalogue(s) will be generated with " << num_objects << " objects * * *" << endl;
    }
    else {
      coutCBL << "* * * "<< n_rec << " random catalogue(s) will be read from file(s)" << endl; 
      for (int rec=0; rec<n_rec; rec++) {
        if (random_catalogue[rec].nObjects() > 0 && (random_catalogue[rec].nObjects() > num_objects || random_catalogue[rec].nObjects() < num_objects)) {   
          if (random_catalogue[rec].nObjects() > num_objects) coutCBL << "WARNING! Random catalogue" << rec+1 << " has " << 
                        random_catalogue[rec].nObjects() - num_objects << " more objects than the tracers catalogue. The number will be equalized by the function catalogue::equalize_random_LC" << endl << endl;
          else coutCBL << "WARNING! Random catalogue has " << 
                        -(random_catalogue[rec].nObjects() - num_objects) << " less objects than the tracers catalogue. The number will be equalized by the function catalogue::equalize_random_LC" << endl << endl; 
        }
      }     
    }
  }

  cout << endl << endl;

  if (algorithm==VoidAlgorithm::_LaZeVo_) {
    
    coutCBL << par::col_green << "                        LaZeVo reconstruction" << par::col_default << endl << endl;
    coutCBL << "* * * Running "<< n_rec << " realizations on " << std::min(omp_get_max_threads(), n_rec) << " cores * * *" << endl << endl;
    
#pragma omp parallel num_threads(omp_get_max_threads())
    {
#pragma omp for    
      for (int rec=0; rec<n_rec; rec++) { 

        vector<vector<double>> coordinatesR(3, vector<double>(num_objects));   
        
        if (random_catalogue.size() == 0) {
          mt19937 gen;
          if (seed > 0) 
            gen.seed(seed + (unsigned int)1e6*rec);
          else {
            random_device rd;
            gen.seed(rd() + (unsigned int)1e6*rec);
          }
          uniform_real_distribution<> distX(minTrCatX, maxTrCatX);
          uniform_real_distribution<> distY(minTrCatY, maxTrCatY);
          uniform_real_distribution<> distZ(minTrCatZ, maxTrCatZ);
          for (size_t i = 0; i < num_objects; ++i) {
            coordinatesR[0][i] = distX(gen);
            coordinatesR[1][i] = distY(gen);
            coordinatesR[2][i] = distZ(gen);
          }
        }
        else if (random_catalogue[rec].nObjects() == 0) {        
          mt19937 gen;
          if (seed > 0) 
            gen.seed(seed + (unsigned int)1e6*rec);
          else {
            random_device rd;
            gen.seed(rd() + (unsigned int)1e6*rec);
          }
          uniform_real_distribution<> distX(minTrCatX, maxTrCatX);
          uniform_real_distribution<> distY(minTrCatY, maxTrCatY);
          uniform_real_distribution<> distZ(minTrCatZ, maxTrCatZ);
          for (size_t i = 0; i < num_objects; ++i) {
            coordinatesR[0][i] = distX(gen);
            coordinatesR[1][i] = distY(gen);
            coordinatesR[2][i] = distZ(gen);
          }
        }
        else if (random_catalogue[rec].nObjects() == num_objects) {
          for (unsigned int i=0; i<num_objects; i++) {
            coordinatesR[0][i] = random_catalogue[rec].xx(i);
            coordinatesR[1][i] = random_catalogue[rec].yy(i);
            coordinatesR[2][i] = random_catalogue[rec].zz(i);
          }
          random_catalogue[rec]={};
        }
        else {
          unsigned int rndd = (int)chrono::steady_clock::now().time_since_epoch().count();
          if (seed > 0) 
            rndd = seed;
          random_catalogue[rec].equalize_random_box(num_objects, {{minTrCatX,maxTrCatX},{minTrCatY,maxTrCatY},{minTrCatZ,maxTrCatZ}}, rndd);
          for (unsigned int i=0; i<num_objects; i++) {
            coordinatesR[0][i] = random_catalogue[rec].xx(i);
            coordinatesR[1][i] = random_catalogue[rec].yy(i);
            coordinatesR[2][i] = random_catalogue[rec].zz(i);
          }
          random_catalogue[rec]={};
        }

        vector<unsigned int> index_tracer_cat(num_objects);
        iota(begin(index_tracer_cat), end(index_tracer_cat), 0); 

        meshgrid::MeshGrid MeshGrid_randoms = meshgrid::MeshGrid(coordinatesR[0], coordinatesR[1], coordinatesR[2], cellsize, lims);

        unsigned int N_near_obj = 3*3*4*par::pi;
        vector<vector<unsigned int>> near_part(num_objects);

        for (size_t i=0; i<num_objects; i++) 
          near_part[i] = MeshGrid_tracers.N_closerObjects(N_near_obj, i); 
        
        vector<unsigned int> randoms(num_objects);
        iota(begin(randoms), end(randoms), 0); 
        mt19937 gen;
        if (seed > 0) 
          gen.seed(seed*(rec+1));
        else {
          random_device rd;
          gen.seed(rd() + (unsigned int)1e6*rec);
        }
        std::shuffle(randoms.begin(), randoms.end(), gen);

        double dist = 4*mps;
        unsigned int removed = 0;
        vector<unsigned int> startingIndexes(num_objects);
        meshgrid::MeshGrid MG_tracer_copy = MeshGrid_tracers;
        meshgrid::MeshGrid MG_random_copy = MeshGrid_randoms;
        iota(startingIndexes.begin(), startingIndexes.end(), 0);
        mt19937 rng;
        if (seed > 0) 
          rng.seed(seed*(rec+1)+1e4);
        else {
          random_device rd;
          rng.seed(rd() + (unsigned int)1e6*rec);
        }
        std::shuffle(startingIndexes.begin(), startingIndexes.end(), rng);
        vector<bool> maskStarters(num_objects, false);
        for (unsigned int i : startingIndexes) {
          if (maskStarters[i] == false) {
            vector<unsigned int> close = MG_tracer_copy.closeObjects(coordinates[0][i], coordinates[1][i], coordinates[2][i], dist);
            unsigned int tr_to_rmv = std::min(31, (int)close.size());
            if (tr_to_rmv > 0) {
              vector<unsigned int> close_random = MG_random_copy.N_closerObjects(close.size(), coordinates[0][i], coordinates[1][i], coordinates[2][i]);
              removed += tr_to_rmv;
              while (tr_to_rmv > 0) {
                uniform_int_distribution<std::mt19937::result_type> distrT(0, close.size()-1);
                uniform_int_distribution<std::mt19937::result_type> distrR(0, close_random.size()-1);
                unsigned int rnd1 = distrT(rng), rnd2 = distrR(rng); 
                randoms[close[rnd1]] = close_random[rnd2];
                MG_tracer_copy.deleteParticle(close[rnd1]);
                maskStarters[close[rnd1]] = true;
                MG_random_copy.deleteParticle(close_random[rnd2]);
                close.erase(close.begin()+rnd1);
                close_random.erase(close_random.begin()+rnd2);
                tr_to_rmv--;
              }
            }
          }
        }
        MG_tracer_copy = {};
        MG_random_copy = {};

        vector<double> distances(num_objects);
        vector<bool> isLonger(num_objects, false);
        double dist_sum = 0, dist_mean;
        for (unsigned int i=0; i<num_objects; i++) {
          distances[i] = sqrt((coordinates[0][i]-coordinatesR[0][randoms[i]])*(coordinates[0][i]-coordinatesR[0][randoms[i]]) + 
                        (coordinates[1][i]-coordinatesR[1][randoms[i]])*(coordinates[1][i]-coordinatesR[1][randoms[i]]) +
                        (coordinates[2][i]-coordinatesR[2][randoms[i]])*(coordinates[2][i]-coordinatesR[2][randoms[i]]));
          dist_sum += distances[i];
        } 
        dist_mean = dist_sum/num_objects;
        for (unsigned int i=0; i<num_objects; i++) {
          if (distances[i]>dist_mean) isLonger[i] = true;
        }

        double ratio = 1.;
        unsigned int n_iter = 0;

        while (ratio > threshold) {
          mt19937 rndd;
          if (seed > 0) 
            rndd.seed(seed*(rec+1)*num_objects);
          else {
            random_device rd;
            rndd.seed(rd()*(rec+1)*num_objects);
          }
          std::shuffle(index_tracer_cat.begin(), index_tracer_cat.end(), rndd);
          vector<bool> index_bool(num_objects, false);
          vector<unsigned int> to_use;

          for (size_t i : index_tracer_cat) 
            if (distances[i] > dist_mean) 
              to_use.push_back(i);
          
          mt19937 rng;
          if (seed > 0) 
            rng.seed(seed*(rec+1)*100);
          else {
            random_device rd;
            unsigned int thread = hash<thread::id>{}(this_thread::get_id()); 
            rng.seed(rd() * (unsigned int)(rec+1)*thread);
          }
          uniform_int_distribution<mt19937::result_type> dist(0, N_near_obj - 1);
          vector<int> permIndexes = {0, 1, 2, 3};
          vector<vector<int>> permutations(24,vector<int>(4));
          int cont=0;
          do {
            for (unsigned int i=0; i<permIndexes.size(); i++) {
              permutations[cont][i] = permIndexes[i];
            }
            cont++;
          } while (next_permutation(permIndexes.begin(), permIndexes.end()));

          for (unsigned int idx=0; idx<to_use.size(); idx++) {
            size_t i = to_use[idx];
            vector<unsigned int> H(4), R_def(4), R(4), R_copy(4);
            unsigned int rand1 = dist(rng), rand2 = dist(rng), rand3 = dist(rng);
            while (rand1 == rand2 || rand1==rand3 || rand2 == rand3) {
              rand2 = dist(rng);
              rand3 = dist(rng);
            }
            H[0]=near_part[i][0];
            H[1]=near_part[i][rand1];
            H[2]=near_part[i][rand2];
            H[3]=near_part[i][rand3];
            R[0]=randoms[H[0]];
            R[1]=randoms[H[1]];
            R[2]=randoms[H[2]];
            R[3]=randoms[H[3]];

            R_copy = R;
            double dist_min = 1.e30;
            for (unsigned int j=0; j<24; j++) {
              double dist=0;
              for (unsigned int k=0; k<4; k++) {
                for (unsigned int n=0; n<3; n++) {
                  dist += (coordinates[n][H[k]]-coordinatesR[n][R[permutations[j][k]]])*(coordinates[n][H[k]]-coordinatesR[n][R[permutations[j][k]]]);
                }
              }
              if (dist < dist_min) {
                dist_min = dist;
                R_def = vector<unsigned int>{R[permutations[j][0]], R[permutations[j][1]], R[permutations[j][2]], R[permutations[j][3]]};
              }          
            }
          
            if(R_def != R_copy) {
              double distSumPrev = 0, distSumPost = 0;
              index_bool[i] = true; 
              for (size_t j=0; j<4; j++) {
                distSumPrev += distances[H[j]];
                distances[H[j]] = sqrt((coordinates[0][H[j]]-coordinatesR[0][R_def[j]])*(coordinates[0][H[j]]-coordinatesR[0][R_def[j]]) + 
                    (coordinates[1][H[j]]-coordinatesR[1][R_def[j]])*(coordinates[1][H[j]]-coordinatesR[1][R_def[j]]) +
                    (coordinates[2][H[j]]-coordinatesR[2][R_def[j]])*(coordinates[2][H[j]]-coordinatesR[2][R_def[j]])); 
                distSumPost += distances[H[j]];
                randoms[H[j]] = R_def[j];
              }
              dist_mean = (dist_mean*num_objects - (distSumPrev-distSumPost))/num_objects;
            }
          }

          unsigned int changed_couples = count(index_bool.begin(), index_bool.end(), true);
          ratio = (double)changed_couples/num_objects;
          n_iter++;
        }

        for (size_t ii = 0; ii<randoms.size(); ii++) {
          displacement_temp[rec][ii][0] = coordinatesR[0][randoms[ii]]-coordinates[0][ii];
          displacement_temp[rec][ii][1] = coordinatesR[1][randoms[ii]]-coordinates[1][ii];
          displacement_temp[rec][ii][2] = coordinatesR[2][randoms[ii]]-coordinates[2][ii];
        }
      }
    } // end of parallelization

    coutCBL << "* * * Calculating the displacement field * * *" << endl;

    if (n_rec>1) {
      
      vector<double> Displacement_finalNorms(num_objects, 0.);
      vector<double> Displacement_varianceSum(num_objects, 0.);
      vector<double> Displacement_varianceSumAng(num_objects, 0.);
      vector<vector<double>> Displacement_angles(n_rec, vector<double>(num_objects, 0.));
      vector<vector<double>> displacementNorm(n_rec, vector<double>(num_objects, 0.));
      vector<vector<double>> displacementFinal_temp(num_objects, vector<double>(3, 0.));

      for (int rec=0; rec<n_rec; rec++ ) {
        for (size_t ii = 0; ii<num_objects; ii++) {
          displacementNorm[rec][ii] = sqrt(displacement_temp[rec][ii][0]*displacement_temp[rec][ii][0]+ 
                                displacement_temp[rec][ii][1]*displacement_temp[rec][ii][1]+
                                displacement_temp[rec][ii][2]*displacement_temp[rec][ii][2]);
          displacement[ii][0] += displacement_temp[rec][ii][0]; 
          displacement[ii][1] += displacement_temp[rec][ii][1];
          displacement[ii][2] += displacement_temp[rec][ii][2];
        }
      }
      
      for (size_t ii = 0; ii<num_objects; ii++) {
        displacementFinal_temp[ii][0] = displacement[ii][0]/n_rec;
        displacementFinal_temp[ii][1] = displacement[ii][1]/n_rec;
        displacementFinal_temp[ii][2] = displacement[ii][2]/n_rec;
        Displacement_finalNorms[ii] = sqrt(displacementFinal_temp[ii][0]*displacementFinal_temp[ii][0]+
                                displacementFinal_temp[ii][1]*displacementFinal_temp[ii][1]+
                                displacementFinal_temp[ii][2]*displacementFinal_temp[ii][2]);
      }

      for (int rec=0; rec<n_rec; rec++ ) {
        for (size_t ii = 0; ii<num_objects; ii++) {
          Displacement_angles[rec][ii] = acos((displacementFinal_temp[ii][0]*displacement_temp[rec][ii][0]+displacementFinal_temp[ii][1]*displacement_temp[rec][ii][1] +
                                    displacementFinal_temp[ii][2]*displacement_temp[rec][ii][2])/(Displacement_finalNorms[ii]*displacementNorm[rec][ii]));
        }
      }

      for (unsigned int i=0; i<num_objects; i++) {
        for (int rec=0; rec<n_rec; rec++) {
          Displacement_varianceSum[i] += (displacementNorm[rec][i] - Displacement_finalNorms[i])*(displacementNorm[rec][i] - Displacement_finalNorms[i]);
          Displacement_varianceSumAng[i] += Displacement_angles[rec][i]*Displacement_angles[rec][i];
        }
        Displacement_finalStd[i] = sqrt(Displacement_varianceSum[i]/(n_rec-1));
        Displacement_finalAngStd[i] = sqrt(Displacement_varianceSumAng[i]/(n_rec-1));
      }

      for (unsigned int i=0; i<num_objects; i++) {
        bool controlN = true;
        bool controlA = true;
        while (controlN || controlA) {
          for (int rec=0; rec<n_rec; rec++) {
            if (displacement_mask[rec][i]) {
              if (displacementNorm[rec][i] > Displacement_finalNorms[i]+Displacement_finalStd[i]*2 || displacementNorm[rec][i] < Displacement_finalNorms[i]-Displacement_finalStd[i]*2) {
                displacement_mask[rec][i] = false;
                displacement[i][0] -= displacement_temp[rec][i][0];
                displacement[i][1] -= displacement_temp[rec][i][1];
                displacement[i][2] -= displacement_temp[rec][i][2];
                Displacement_finalNorms[i] = Displacement_finalNorms[i] + (Displacement_finalNorms[i] - displacementNorm[rec][i])/(n_rec_eff[i]-1);
                n_rec_eff[i]--;
                continue;
              }
              else controlN = false;

              if (Displacement_angles[rec][i] > 2*Displacement_finalAngStd[i]) {
                displacement_mask[rec][i] = false;
                displacement[i][0] -= displacement_temp[rec][i][0];
                displacement[i][1] -= displacement_temp[rec][i][1];
                displacement[i][2] -= displacement_temp[rec][i][2];
                Displacement_finalNorms[i] = Displacement_finalNorms[i] + (Displacement_finalNorms[i] - displacementNorm[rec][i])/(n_rec_eff[i]-1);
                n_rec_eff[i]--;
              }
              else controlA = false;
            }
          }
          for (int rec=0; rec<n_rec; rec++) {
            if (displacement_mask[rec][i]) { 
              Displacement_varianceSum[i] += (displacementNorm[rec][i] - Displacement_finalNorms[i])*(displacementNorm[rec][i] - Displacement_finalNorms[i]);
              Displacement_varianceSumAng[i] += Displacement_angles[rec][i]*Displacement_angles[rec][i];
            }
          }
          Displacement_finalStd[i] = sqrt(Displacement_varianceSum[i]/(n_rec_eff[i]-1));
          Displacement_finalAngStd[i] = sqrt(Displacement_varianceSumAng[i]/(n_rec_eff[i]-1));
        }
      } 
    
      for (size_t ii = 0; ii<displacement.size(); ii++) {
        displacement[ii][0] = displacement[ii][0]/n_rec_eff[ii];
        random_final[0][ii] = coordinates[0][ii] + displacement[ii][0];
        displacement[ii][1] = displacement[ii][1]/n_rec_eff[ii];
        random_final[1][ii] = coordinates[1][ii] + displacement[ii][1];
        displacement[ii][2] = displacement[ii][2]/n_rec_eff[ii];
        random_final[2][ii] = coordinates[2][ii] + displacement[ii][2];
      }  
    }
    else {
      for (size_t ii = 0; ii<displacement.size(); ii++) {
        random_final[0][ii] = coordinates[0][ii] + displacement_temp[0][ii][0];
        random_final[1][ii] = coordinates[1][ii] + displacement_temp[0][ii][1];
        random_final[2][ii] = coordinates[2][ii] + displacement_temp[0][ii][2];
      }  
    }

  } // End of LaZeVo method
  else if (algorithm==VoidAlgorithm::_Exact_) { 
    
    coutCBL << par::col_green << "                         Exact reconstruction" << par::col_default << endl << endl;

    vector<unsigned int> randoms(num_objects);
    iota(begin(randoms), end(randoms), 0); 

    if (random_catalogue.size()!=1 || n_rec!=1) 
      ErrorCBL ("You should provide a vector with a sigle random catalogue. Moreover, n_rec must be = 1!", "Catalogue", "VoidCatalogue.cpp");

    if (random_catalogue[0].nObjects()!=num_objects) 
      ErrorCBL ("Random catalogue must have the same dimension of the tracers one!", "Catalogue", "VoidCatalogue.cpp");
    
    vector<vector<double>> coordinatesR(3, vector<double>(num_objects));   
    coutCBL << "* * * Reading random catalogue of " << random_catalogue[0].nObjects() << " objects provides by user * * *" << endl << endl;
#pragma omp parallel num_threads(omp_get_max_threads())
    {     
#pragma omp for schedule(static)
      for (unsigned int i=0; i<num_objects; i++) {
        coordinatesR[0][i] = random_catalogue[0].xx(i);
        coordinatesR[1][i] = random_catalogue[0].yy(i);
        coordinatesR[2][i] = random_catalogue[0].zz(i);
      }
    }
    random_catalogue[0] = {};

    for (size_t ii = 0; ii<randoms.size(); ii++) {
      displacement[ii][0] = coordinatesR[0][randoms[ii]]-coordinates[0][ii];
      random_final[0][0] = coordinatesR[0][randoms[ii]];
      displacement[ii][1] = coordinatesR[1][randoms[ii]]-coordinates[1][ii];
      random_final[1][1] = coordinatesR[1][randoms[ii]];
      displacement[ii][2] = coordinatesR[2][randoms[ii]]-coordinates[2][ii];
      random_final[2][2] = coordinatesR[2][randoms[ii]];
    }
  
    string name =  "../output/Displacement_";
    name.append(output);
    ofstream fout_displ(name.c_str());
    cout << endl << endl;
    coutCBL << "Writing the displacement field in " << name << " ..." << endl;
    
    for (unsigned int i=0; i<displacement.size(); i++) {
      fout_displ.precision(7);
          fout_displ << coordinates[0][i] << " " << coordinates[1][i] << " " << coordinates[2][i] << " "
                     << coordinatesR[0][i] << " " << coordinatesR[1][i] << " " << coordinatesR[2][i] << " "
                     << displacement[i][0] << " " << displacement[i][1] << " " << displacement[i][2] << endl;   
    } 
    fout_displ.clear(); fout_displ.close();  

  }
  else ErrorCBL ("algorithm type is not correct!", "Catalogue", "VoidCatalogue.cpp");

  // --------------------------------------------- //
  // ---------------- Second Step ---------------- //
  // --------------------------------------------- //
    
  double step = step_size*mps; 

  vector<double> min = {lims[0][0], lims[1][0], lims[2][0]};
  vector<double> max = {lims[0][1], lims[1][1], lims[2][1]};

  vector<double> cells(3);
  for(unsigned int&& i=0; i<3; i++) cells[i] = (max[i]-min[i])/step + 1;
  unsigned int nCells = *max_element(cells.begin(), cells.end());
  auto min_copy = min;
  for(unsigned int&& i=0; i<3; i++) 
    min[i] = min_copy[i] - (step*nCells-(max[i]-min_copy[i]))/2;
  
  vector<vector<vector<double>>> Divergence(nCells, vector<vector<double>>(nCells, vector<double>(nCells, 0.)));

  cout << endl;
  coutCBL << "* * * Estimation of the Divergence field * * *" << endl;
  cout << endl;
  for (int rec=0; rec<n_rec; rec++) {
    for (size_t i=0; i<random_final[0].size(); i++) {
      if (displacement_mask[rec][i]) { 
        vector<vector<double>> coord_part(2);
        coord_part[0] = {coordinates[0][i]+displacement_temp[rec][i][0], coordinates[1][i]+displacement_temp[rec][i][1], coordinates[2][i]+displacement_temp[rec][i][2]};
        coord_part[1] = {coordinates[0][i], coordinates[1][i], coordinates[2][i]};
        vector<double> displ = {displacement_temp[rec][i][0], displacement_temp[rec][i][1], displacement_temp[rec][i][2]};
        vector<vector<unsigned int>> inds(2, vector<unsigned int>(3));
        for (unsigned int&& j=0; j<3; j++) {
          inds[0][j] = (coord_part[0][j]-min[j])/step;
          inds[1][j] = (coord_part[1][j]-min[j])/step;
        }
        if(inds[0]!=inds[1]) {
          for (unsigned int part=0; part<2; part++) {
            for (unsigned int j=0; j<2; j++) { 
              for (unsigned int&& k=0; k<3; k++) {
                double t = (coord_part[0][k]-((inds[part][k]+j)*step+min[k]))/displ[k];
                if (t>=0. && t<=1.) {
                  vector<double> inters_coord = {coord_part[0][0]-(coord_part[0][0]-coord_part[1][0])*t,
                    coord_part[0][1]-(coord_part[0][1]-coord_part[1][1])*t, 
                    coord_part[0][2]-(coord_part[0][2]-coord_part[1][2])*t};
                  if ((k==0) ? (inters_coord[1]>(inds[part][1]*step+min[1]) && inters_coord[1] <= ((inds[part][1]+1)*step+min[1]) &&
                                inters_coord[2]>(inds[part][2]*step+min[2]) && inters_coord[2] <= ((inds[part][2]+1)*step+min[2])) :
                    ((k==1) ? (inters_coord[0]>(inds[part][0]*step+min[0]) && inters_coord[0] <= ((inds[part][0]+1)*step+min[0]) &&
                                inters_coord[2]>(inds[part][2]*step+min[2]) && inters_coord[2] <= ((inds[part][2]+1)*step+min[2])) :
                              (inters_coord[0]>(inds[part][0]*step+min[0]) && inters_coord[0] <= ((inds[part][0]+1)*step+min[0]) &&
                                inters_coord[1]>(inds[part][1]*step+min[1]) && inters_coord[1] <= ((inds[part][1]+1)*step+min[1])))) {
                    if(part==0) Divergence[inds[part][0]][inds[part][1]][inds[part][2]] -= abs(displ[k])/(n_rec_eff[i]*step);
                    if(part==1) Divergence[inds[part][0]][inds[part][1]][inds[part][2]] += abs(displ[k])/(n_rec_eff[i]*step);
                  }
                }
                else continue;
              }
            }
          }
        }
      }
    }
  }

  vector<unsigned int> tracCell(num_objects);        
  vector<unsigned int> randCell(num_objects);        
  
  for (unsigned int i=0; i<num_objects; i++) {
    tracCell[i] = ((coordinates[0][i]-min[0])/step)*nCells*nCells + ((coordinates[1][i]-min[1])/step)*nCells + ((coordinates[2][i]-min[2])/step);
    randCell[i] = ((random_final[0][i]-min[0])/step)*nCells*nCells + ((random_final[1][i]-min[1])/step)*nCells + ((random_final[2][i]-min[2])/step);
  }

  if (print[0]) {
    if (n_rec>1) {
      string name =  "../output/Displacement_";
      name.append(output);
      ofstream fout_displ(name.c_str());
      coutCBL << "Writing the displacement field in -> " << name << " ..." << endl << endl;
      fout_displ <<  "###   tracX   tracY   tracZ   randX   randY   randZ   displX   displY   displZ   stdDispl   stdAng   N_rec   divCellTrac   divCellRand  " << endl;
      for (unsigned int i=0; i<displacement.size(); i++) {
        fout_displ.precision(7);
              fout_displ << coordinates[0][i] << " " << coordinates[1][i] << " " << coordinates[2][i] << " "
                        << random_final[0][i] << " " << random_final[1][i] << " " << random_final[2][i] << " "
                        << displacement[i][0] << " " << displacement[i][1] << " " << displacement[i][2] << " " 
                        << Displacement_finalStd[i] << " " << Displacement_finalAngStd[i] << " " << n_rec_eff[i] << " " <<
                        tracCell[i] << " " << randCell[i] << endl;   
      } 
      fout_displ.clear(); fout_displ.close();  
    }
    else {
      string name =  "../output/Displacement_";
      name.append(output);
      ofstream fout_displ(name.c_str());
      cout << endl;
      coutCBL << "Writing the displacement field in -> " << name << " ..." << endl;
      fout_displ <<  "###   tracX   tracY   tracZ   randX   randY   randZ   displX   displY   displZ   divCellTrac   divCellRand  " << endl;
      for (unsigned int i=0; i<displacement.size(); i++) {
        fout_displ.precision(7);
            fout_displ << coordinates[0][i] << " " << coordinates[1][i] << " " << coordinates[2][i] << " "
                      << random_final[0][i] << " " << random_final[1][i] << " " << random_final[2][i] << " "
                      << displacement_temp[0][i][0] << " " << displacement_temp[0][i][1] << " " << displacement_temp[0][i][2] 
                      << tracCell[i] << " " << randCell[i] << endl;   
      } 
      fout_displ.clear(); fout_displ.close();  
    }
  }

  coutCBL << "* * * Interpolating * * *" << endl << endl;

  int emptyCellsCounter = 1;

  while (emptyCellsCounter > 0) {
    vector<vector<vector<double>>> Divergence_copy = Divergence;
    emptyCellsCounter = 0;
    for (unsigned int i=0; i<nCells; i++) {
      for (unsigned int j=0; j<nCells; j++) {
        for (unsigned int k=0; k<nCells; k++) {
          if (Divergence[i][j][k] == 0) {
            double num = 0;
            int den = 0;
            for (int ii=std::max(0,(int)i-1); ii<=std::min((int)nCells-1,(int)i+1); ii++) {
              for (int jj=std::max(0,(int)j-1); jj<=std::min((int)nCells-1,(int)j+1); jj++) {
                for (int kk=std::max(0,(int)k-1); kk<=std::min((int)nCells-1,(int)k+1); kk++) {
                  if (Divergence[ii][jj][kk] != 0) {
                    num += Divergence[ii][jj][kk];
                    den += 1;
                  }
                }
              }
            }
            if (den == 0) emptyCellsCounter++;
            else Divergence_copy[i][j][k] = num/den;
          }
        }
      }
    }
    Divergence = Divergence_copy;
  }

  if (sigma>0.) {
    coutCBL << "* * * Applying Gaussian smoothing with sigma = " << std::fixed << std::setprecision(3) << sigma*mps << " [Mpc/h] * * *" << endl << endl;

    vector<vector<vector<double>>> Divergence_copy(nCells, vector<vector<double>>(nCells, vector<double>(nCells, 0.)));

    sigma = sigma*mps;
    int sigmaCells = 3*sigma/step; 

    for (unsigned int i=0; i<nCells; i++) {
      for (unsigned int j=0; j<nCells; j++) {
        for (unsigned int k=0; k<nCells; k++) {
          for (int ii=std::max(0,(int)i-sigmaCells); ii<std::min((int)nCells,(int)i+sigmaCells+1); ii++) {
            for (int jj=std::max(0,(int)j-sigmaCells); jj<std::min((int)nCells,(int)j+sigmaCells+1); jj++) {
              for (int kk=std::max(0,(int)k-sigmaCells); kk<std::min((int)nCells,(int)k+sigmaCells+1); kk++) {
                double ex = step*step*step/pow(2*par::pi*sigma*sigma, 1.5)*exp(-((step*(i-ii))*(step*(i-ii))+(step*(j-jj))*(step*(j-jj))+(step*(k-kk))*(step*(k-kk)))/(2*sigma*sigma));
                Divergence_copy[ii][jj][kk] += Divergence[i][j][k]*ex;
              }
            }
          }
        }
      }
    }  
    Divergence = Divergence_copy;
  }

  coutCBL << "* * * Identification of voids * * *" << endl << endl;

  vector<vector<double>> coord_locmin(6);
  unsigned int cont = 0;

  for (unsigned int i=0; i<nCells; i++) {
    for (unsigned int j=0; j<nCells; j++) {
      for (unsigned int k=0; k<nCells; k++) {
        if (Divergence[i][j][k] < 0.) {
          bool control = false; 
          double sum = 0.;
          vector<double> num(3, 0.), den(3, 0.);
          vector<int> delta_cells(3, 0);
          cont++;
          for (int ii=i-1; ii<(int)i+2; ii++) {
            delta_cells[0] = i-ii; 
            for (int jj=j-1; jj<(int)j+2; jj++) {
              delta_cells[1] = j-jj; 
              for (int kk=k-1; kk<(int)k+2; kk++) {
                delta_cells[2] = k-kk; 
                if (!(ii==(int)i && jj==(int)j && kk==(int)k) && (ii>=0 && jj>=0 && kk>=0) && (ii<(int)nCells && jj<(int)nCells && kk<(int)nCells)) {
                  sum +=  Divergence[ii][jj][kk];
                  for (unsigned int t=0; t<3; t++) {
                    num[t] += Divergence[ii][jj][kk]*(step/Euclidean_distance((double)i, (double)ii, (double)j, (double)jj, (double)k, (double)kk))*delta_cells[t];
                    den[t] += abs(Divergence[ii][jj][kk]);
                  }
                  if (Divergence[ii][jj][kk]<=Divergence[i][j][k]) {
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
            coord_locmin[0].push_back(min[0]+step*(i+0.5)-num[0]/(den[0]*2));
            coord_locmin[1].push_back(min[1]+step*(j+0.5)-num[1]/(den[1]*2));
            coord_locmin[2].push_back(min[2]+step*(k+0.5)-num[2]/(den[2]*2));
            coord_locmin[3].push_back(i);
            coord_locmin[4].push_back(j);
            coord_locmin[5].push_back(k);
          }
          (void)sum;
        }
      }
    }
  }
  
  coutCBL << "Number of cells with negative divergence: " << cont << endl << endl;

  if (coord_locmin[0].size()==0) 
    ErrorCBL("No local minima found!", "Catalogue", "VoidCatalogue.cpp");

  //coutCBL << "Number of proto-voids: " << coord_locmin[0].size() << endl << endl;
   
  // -------------------------------------------- //
  // ---------------- Third Step ---------------- //
  // -------------------------------------------- //
  
  coutCBL << "* * * Running watershed algoritmh * * *" << endl << endl;
  vector<double> radius_void(coord_locmin[0].size());
  vector<vector<vector<int>>> voidIndexes(nCells, vector<vector<int>>(nCells, vector<int>(nCells, int(-1))));
  vector<vector<vector<int>>> cellIndexes(nCells, vector<vector<int>>(nCells, vector<int>(nCells, int(-1))));
  vector<vector<vector<vector<int>>>> flowCells(nCells, vector<vector<vector<int>>>(nCells, vector<vector<int>>(nCells, vector<int>(3, 0))));
  vector<vector<vector<std::optional<bool>>>> divergenceMask(nCells, vector<vector<std::optional<bool>>>(nCells, vector<std::optional<bool>>(nCells, std::optional<bool>(false))));

  for (unsigned int i=0; i<nCells; i++) {
    for (unsigned int j=0; j<nCells; j++) {
      for (unsigned int k=0; k<nCells; k++) {
        if (Divergence[i][j][k]<0.) divergenceMask[i][j][k] = true;
      }
    }
  }
  for (unsigned int i=0; i<coord_locmin[3].size(); i++) {
    divergenceMask[coord_locmin[3][i]][coord_locmin[4][i]][coord_locmin[5][i]] = nullopt;
    voidIndexes[coord_locmin[3][i]][coord_locmin[4][i]][coord_locmin[5][i]] = i;
  }

  for (unsigned int i=0; i<nCells; i++) {
    for (unsigned int j=0; j<nCells; j++) {
      for (unsigned int k=0; k<nCells; k++) {
        if (*divergenceMask[i][j][k]) {
          double min = 0.;
          vector<int> temp_coordMin = {0,0,0};
          for (int ii=(int)i-1; ii<(int)i+2; ii++) {
            for (int jj=(int)j-1; jj<(int)j+2; jj++) {
              for (int kk=(int)k-1; kk<(int)k+2; kk++) {
                if ((ii>=0 && jj>=0 && kk>=0) && (ii<(int)nCells && jj<(int)nCells && kk<(int)nCells)) {
                  if (Divergence[ii][jj][kk]<min) {
                    min = Divergence[ii][jj][kk];
                    temp_coordMin = {ii-(int)i, jj-(int)j, kk-(int)k};
                  }
                }
              }
            }
          }
          flowCells[i][j][k][0] = temp_coordMin[0];
          flowCells[i][j][k][1] = temp_coordMin[1];
          flowCells[i][j][k][2] = temp_coordMin[2];
        }
      }
    }
  }

  for (unsigned int i=0; i<nCells; i++) {
    for (unsigned int j=0; j<nCells; j++) {
      for (unsigned int k=0; k<nCells; k++) {
        if (*divergenceMask[i][j][k]) {
          bool stopCond = false;
          int di=0, dj=0, dk=0;
          while (!stopCond) {
            int di_temp = flowCells[i+di][j+dj][k+dk][0];
            int dj_temp = flowCells[i+di][j+dj][k+dk][1];
            int dk_temp = flowCells[i+di][j+dj][k+dk][2];
            di += di_temp; 
            dj += dj_temp;
            dk += dk_temp;
            if (!divergenceMask[i+di][j+dj][k+dk].has_value()) {
              stopCond = true;
              cellIndexes[i][j][k] = voidIndexes[(i+di)][(j+dj)][(k+dk)];
            }
          }
        }
      }
    }
  }

  if (print[1]) {
    string name; 
    name = dir_output + "Divergence";
    name.append("_");
    name.append(output);
    ofstream fout_divergence(name.c_str());

    fout_divergence.precision(7);
    
    coutCBL << "Writing the divergence field in -> " << name << " ..." << endl;
    cout << endl;

    for (unsigned int&& i=0; i<nCells; i++) 
      for (unsigned int&& j=0; j<nCells; j++) 
      	for (unsigned int&& k=0; k<nCells; k++) { 
	  fout_divergence << min[0]+step*(i+0.5) << " " << min[1]+step*(j+0.5) << " " << min[2]+step*(k+0.5) << " " << Divergence[i][j][k] << " " << cellIndexes[i][j][k] << endl;
	}
    fout_divergence.clear(); fout_divergence.close();
  }

  coutCBL << "* * * Voids rescaling * * *" << endl << endl;

  vector<unsigned int> cellsInVoids(coord_locmin[0].size(), 0);
  for (unsigned int i=0; i<nCells; i++) 
    for (unsigned int j=0; j<nCells; j++) 
      for (unsigned int k=0; k<nCells; k++) 
        if (*divergenceMask[i][j][k]) 
          cellsInVoids[cellIndexes[i][j][k]]++;
  
  //unsigned int voidThreshold = 4*par::pi*(1.5*1.5*1.5)/(3*step_size*step_size*step_size);
  vector<double> radii(coord_locmin[0].size(), 0.);
  for (unsigned int i=0; i<coord_locmin[0].size(); i++)
  //  if (cellsInVoids[i] > voidThreshold)
      radii[i] = step*pow(3*cellsInVoids[i]/(4*par::pi), 1./3);

  vector<shared_ptr<catalogue::Object>> void_obj_cat;
  for (size_t ii = 0; ii<coord_locmin[0].size(); ii++)
  { 
    comovingCoordinates coord_final_void = {coord_locmin[0][ii], coord_locmin[1][ii], coord_locmin[2][ii]};
    auto void_obj = make_shared<catalogue::Void>(coord_final_void);
    void_obj_cat.push_back(void_obj); 
  }
  add_objects(void_obj_cat);

/*
  string name;
  name = dir_output+"ProtoVoids_";
  name.append(output);
  coutCBL << "Saving proto-voids in -> " << name << " ..." << endl << endl;
  ofstream fout_protovoids(name.c_str());
  for(int i=0; i<(int)radii.size(); ++i) 
      fout_protovoids << coord_locmin[0][i] << " " << coord_locmin[1][i] << " " << coord_locmin[2][i] << " " << radii[i] << endl;
  fout_protovoids.clear(); fout_protovoids.close();
*/

  string name = dir_output+"Voids_";
  name.append(output);
  ofstream fout_voids(name.c_str());
  // int countVoids = 0;
  coutCBL << "Saving voids in -> " << name << " ..." << endl << endl;
  for(int i=0; i<(int)radii.size(); ++i) {
    //if(cellsInVoids[i] > voidThreshold) {
      set_var(i, Var::_Radius_, radii[i]);
      fout_voids << coord_locmin[0][i] << " " << coord_locmin[1][i] << " " << coord_locmin[2][i] << " " << radii[i] << endl;
    // countVoids++;
    //}
  }
  fout_voids.clear(); fout_voids.close();

  coutCBL << "Number of voids: " << coord_locmin[0].size() << endl << endl;

  // ---------------------------------------------------------- //
  // ---------------- Computing execution time ---------------- //
  // ---------------------------------------------------------- //
  
  double seconds = omp_get_wtime() - start_time;
  
  cout << endl;
  coutCBL << "Time spent to compute: " << seconds       << " seconds " << endl;
  coutCBL << "Time spent to compute: " << seconds/60.   << " minutes " << endl;
  coutCBL << "Time spent to compute: " << seconds/3600. << " hours "   << endl;
  cout << endl;
}

///========================================================================================

cbl::catalogue::Catalogue::Catalogue (const VoidAlgorithm algorithm, Catalogue tracer_catalogue, Catalogue random_catalogue, const std::string dir_output, const std::string output, const double cellsize, const std::shared_ptr<cosmology::Cosmology> cosmology, const std::vector<double> RA_range, const std::vector<double> DEC_range, const int n_rec, const double step_size, const double threshold)
{
  cout << endl; coutCBL << par::col_green << "BitVF for Light cones" << par::col_default << endl << endl;
  ErrorCBL ("Still not avaiable!", "Catalogue", "VoidCatalogue.cpp");
  cout << dir_output << output << cellsize << n_rec << step_size << threshold << cosmology->Omega_baryon() << endl;
  Print(RA_range);
  Print(DEC_range);
  if (algorithm == VoidAlgorithm::_LaZeVo_) cout << "aaa" << endl;
  if (tracer_catalogue.nObjects() == 0) cout << "aaa" << endl;
  if (random_catalogue.nObjects() == 0) cout << "aaa" << endl;

  /*    
  // -------------------------------------------- //
  // ---------------- First Step ---------------- //
  // -------------------------------------------- //
    
  const double start_time = omp_get_wtime();
  shared_ptr<Catalogue> tracer_cat, random_cat;
  vector<Catalogue> displacement_catalogue(n_rec);
  vector<vector<unsigned int>> randoms;
  vector<double> cells_edges;
  std::setprecision(5);
  vector<vector<double>> prop = tracer_catalogue.compute_catalogueProperties_lightCone(cosmology, RA_range, DEC_range, 25);


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

  Catalogue random_catalogue; 
      
  long int time = (long int)chrono::steady_clock::now().time_since_epoch().count();
  //time = (time%1000000);
  time = 2343;
  int rndd = time; 

  double mps_trcat = tracer_catalogue.mps();
      
  vector<vector<double>> data_trcat(3);
  data_trcat[0] = prop[0];
  for (size_t i=0; i<data_trcat[0].size(); i++) data_trcat[0][i] = cosmology->D_C(data_trcat[0][i]);
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
       
  const data::Data1D data(data_trcat[0], data_trcat[1], data_trcat[2]);
  shared_ptr<data::Data> ptr_data = make_shared<data::Data1D>(data);
      
  const int nparameters = 4;
  const vector<string> parNames = {"A", "B", "C", "D"};
  double valA = 0., valB = 0., valC = 0., valD = data_trcat[1][0];     
    
  vector<statistics::ParameterType> parType(nparameters, statistics::ParameterType::_Base_);
  auto ptr_modelInput = make_shared<cosmology::Cosmology>(cosm);
  
  const statistics::Model1D model(poly, nparameters, parType, parNames, ptr_modelInput);
  auto ptr_model = make_shared<statistics::Model1D>(model);
    
  statistics::Likelihood likelihood(ptr_data, ptr_model, statistics::LikelihoodType::_Gaussian_Error_);
  
  const vector<double> starting_p = {valA, valB, valC, valD};
  double minA = -1., maxA = 1.;
  double minB = -1., maxB = 1.;
  double minC = -1., maxC = 1.;
  double minD = 0., maxD = data_trcat[1][0];
  const vector<vector<double>> limits = { {minA, maxA}, {minB, maxB}, {minC, maxC}, {minD, maxD} };

  likelihood.maximize(starting_p, limits);

  const double start=cosmology->D_C(Min(tracer_catalogue.var(Var::_Redshift_)))-data_trcat[1][0]/2;
  const double end=cosmology->D_C(Max(tracer_catalogue.var(Var::_Redshift_)))+data_trcat[1][data_trcat[1].size()-1]/2;
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
  random_catalogue = Catalogue(RandomType::_createRandom_homogeneous_LC_, tracer_catalogue, 1., cosmology, RA_range, DEC_range, 100, rndd);
  coutCBL << "... random catalogue created!" << endl;
  }
  else if (random_catalogue.nObjects() == tracer_catalogue.nObjects())
  {
  coutCBL << "Reading random catalogue of " << random_catalogue.nObjects() << " objects provides by user..." << endl;
  random_catalogue = random_catalogue;
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
  random_catalogue.equalize_random_lightCone(tracer_catalogue, cosmology, RA_range, DEC_range, rndd);
  coutCBL << "...Done!" << endl;
  }

  // shared pointers to the tracer and the random catalogues
  random_cat = make_shared<catalogue::Catalogue> (catalogue::Catalogue(move(random_catalogue)));
  tracer_cat = make_shared<catalogue::Catalogue> (catalogue::Catalogue(move(tracer_catalogue)));  

  // maximum for the number of objects in the tracer catalogue (unsigned int) 
  unsigned int num_objects = ((tracer_cat->nObjects())%2==0) ? tracer_cat->nObjects() : (tracer_cat->nObjects())-1;

  if ((tracer_cat->nObjects())%2 != 0) {
  random_cat -> remove_object(num_objects);
  tracer_cat -> remove_object(num_objects);
  }

  // vector for the indices of tracers
  std::vector<unsigned int> index_tracer_cat(num_objects);
  std::iota (std::begin(index_tracer_cat), std::end(index_tracer_cat), 0); 
  cout << "num obj: " << num_objects << endl;
  cout << endl;
  coutCBL << "Creating ChainMesh... " << flush;
  // chainmesh setting
  vector<Var> variables = {Var::_X_, Var::_Y_, Var::_Z_};
  CatalogueChainMesh ChainMesh_tracers = CatalogueChainMesh(variables, cellsize, tracer_cat);
  CatalogueChainMesh ChainMesh_randoms = CatalogueChainMesh(variables, cellsize, random_cat);
  unsigned int N_near_obj = 3*3*4*par::pi;
  vector<vector<unsigned int>> near_part = ChainMesh_tracers.N_nearest_objects_cat(N_near_obj); //the selected halo (that will have distance = 0) and the nearest 3))
  randoms.resize(n_rec, vector<unsigned int>(num_objects)); //tracers fixed
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
  double minX = ChainMesh_tracers.lim()[0][0]+dist/10, maxX = ChainMesh_tracers.lim()[0][1]-dist/10;
  double minY = ChainMesh_tracers.lim()[1][0]+dist/10, maxY = ChainMesh_tracers.lim()[1][1]-dist/10;
  double minZ = ChainMesh_tracers.lim()[2][0]+dist/10, maxZ = ChainMesh_tracers.lim()[2][1]-dist/10;

  // FACCIO STAMPARE DISPL1 E CONFRONTO CON DISPL0  
  // RISOLVERE QUESTIONE COPIA CHAIN MESH E POI TESTARE (SEMBRA VADA MALE E LENTA) E RISOLVERE PROBLEMA IDENTIFICAZIONE VUOTI
  for (int rec=0; rec<n_rec; rec++) {    
        
  auto timeS = ((long int)chrono::steady_clock::now().time_since_epoch().count())%1000000;
  random::UniformRandomNumbers randX(minX, maxX, timeS);
  random::UniformRandomNumbers randY(minY, maxY, timeS+1);
  random::UniformRandomNumbers randZ(minZ, maxZ, timeS+2);

  CatalogueChainMesh ChM_tracer_copy = ChainMesh_tracers;
  CatalogueChainMesh ChM_random_copy = ChainMesh_randoms;
  cout << ChM_tracer_copy.part(7).size() << " " << ChainMesh_tracers.part(7).size() << endl;
  cout << ChM_tracer_copy.cellsize() << endl;
  unsigned int removed = 0;

  while (removed < num_objects) {
  vector<double> pos = {randX(), randY(), randZ()};
  vector<unsigned int> close=ChM_tracer_copy.Close_objects(pos, dist);
  unsigned int tr_to_rmv = std::min(100, (int)close.size());
  if (tr_to_rmv > 0) {
  vector<unsigned int> close_random = ChM_random_copy.N_nearest_objects(pos, close.size());
  std::random_device dev;
  std::mt19937 rng(dev());
  removed += tr_to_rmv;
  while (tr_to_rmv > 0) {
  std::uniform_int_distribution<std::mt19937::result_type> dist(0, tr_to_rmv-1);
  int rnd1 = dist(rng), rnd2 = dist(rng);
  randoms[rec][close[rnd1]] = close_random[rnd2];
  ChM_tracer_copy.deletePart(close[rnd1]);
  ChM_random_copy.deletePart(close_random[rnd2]);
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
  auto index_tracer_cat_copy = index_tracer_cat;
  double ratio = 1;
  unsigned int n_iter = 0;

  while (ratio > threshold) {

  index_tracer_cat_copy = index_tracer_cat; 
  auto time2 = (long int)chrono::steady_clock::now().time_since_epoch().count();
  time2 = (time2%1000000);
  std::shuffle(index_tracer_cat_copy.begin(), index_tracer_cat_copy.end(), default_random_engine(time2));
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<std::mt19937::result_type> dist(1,N_near_obj-1);
  vector<bool> index_bool(num_objects, false), used(num_objects, false);
  //vector<vector<unsigned int>> H(tracer_cat_len_4,vector<unsigned int>(4)), R_def(tracer_cat_len_4,vector<unsigned int>(4));
    
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
  if (used[near_part[index_tracer_cat_copy[i]][0]] == false && used[near_part[index_tracer_cat_copy[i]][rand1]] == false && 
  used[near_part[index_tracer_cat_copy[i]][rand2]] == false && used[near_part[index_tracer_cat_copy[i]][rand3]] == false) { 
  // selezione particelle coinvolte
  used[near_part[index_tracer_cat_copy[i]][0]] = true;
  used[near_part[index_tracer_cat_copy[i]][rand1]] = true;
  used[near_part[index_tracer_cat_copy[i]][rand2]] = true;
  used[near_part[index_tracer_cat_copy[i]][rand3]] = true;
  H[0]=near_part[index_tracer_cat_copy[i]][0];
  H[1]=near_part[index_tracer_cat_copy[i]][rand1];
  H[2]=near_part[index_tracer_cat_copy[i]][rand2];
  H[3]=near_part[index_tracer_cat_copy[i]][rand3];
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
  dist += Euclidean_distance(tracer_cat->xx(H[j]), random_cat->xx(R[j]), tracer_cat->yy(H[j]), random_cat->yy(R[j]), tracer_cat->zz(H[j]), random_cat->zz(R[j]));
                
  if (dist < dist_min) {
  dist_min = dist;
  R_def = R;
  } 
  } while(std::next_permutation(R.begin(), R.end()));
        
  //aggiornamento vettori tracers e randoms
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
  distanza += Euclidean_distance(tracer_cat->xx(i), random_cat->xx(randoms[rec][i]),
  tracer_cat->yy(i), random_cat->yy(randoms[rec][i]), tracer_cat->zz(i), random_cat->zz(randoms[rec][i]));
  }
  cout << "distanza finale: " << distanza/num_objects << endl;

  vector<shared_ptr<catalogue::Object>> displ_catalogue_object;

  for (size_t ii = 0; ii<randoms[rec].size(); ii++)
  { 
  comovingCoordinates displacement_values = {random_cat->xx(randoms[rec][ii])-tracer_cat->xx(ii), random_cat->yy(randoms[rec][ii])-tracer_cat->yy(ii), random_cat->zz(randoms[rec][ii])-tracer_cat->zz(ii)};
  auto tracer_displ = make_shared<catalogue::Halo>(displacement_values);
  displ_catalogue_object.push_back(tracer_displ); 
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
  fout_displ << tracer_cat->xx(i) << " " << tracer_cat->yy(i) << " " << tracer_cat->zz(i) << " "
  << random_cat->xx(randoms[rec][i]) << " " << random_cat->yy(randoms[rec][i]) << " " << random_cat->zz(randoms[rec][i]) << " "
  << displacement_catalogue[rec].xx(i) << " " << displacement_catalogue[rec].yy(i) << " " << displacement_catalogue[rec].zz(i) << endl;   
  } // modified
          
  fout_displ.clear(); fout_displ.close();  
  }
  }


  } // End of LaZeVo method
  else if (algorithm==VoidAlgorithm::_Exact_) 
  {
          
  } // End of Exact method
  else ErrorCBL ("algorithm type is not correct!", "Catalogue", "VoidCatalogue.cpp");

  // --------------------------------------------- //
  // ---------------- Second Step ---------------- //
  // --------------------------------------------- //

  vector<shared_ptr<Catalogue>> displ_cat(n_rec);
  for (int i=0; i<n_rec; i++) displ_cat[i] = make_shared<Catalogue>(Catalogue(move(displacement_catalogue[i]))); 
  
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
  for (size_t i=0; i<cells_edges.size(); i++) z_cells_edges[i] = cosmology->Redshift(cells_edges[i]);
  for (size_t i=0; i<step.size(); i++) cells_centre[i]=(cells_edges[i+1]+cells_edges[i])/2;
  for (size_t i=0; i<step.size(); i++) z_cells_centres[i]=cosmology->Redshift(cells_centre[i]);
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
  coord_part[0] = {random_cat->redshift(randoms[rec][i]), -(random_cat->dec(randoms[rec][i]))+par::pi/2, 
  (DEC_range[0] < 0. && random_cat->ra(randoms[rec][i]) > par::pi) ? random_cat->ra(randoms[rec][i])-2*par::pi : random_cat->ra(randoms[rec][i]),
  random_cat->xx(randoms[rec][i]), random_cat->yy(randoms[rec][i]), random_cat->zz(randoms[rec][i])};
  coord_part[1] = {tracer_cat->redshift(i), -(tracer_cat->dec(i))+par::pi/2, 
  (DEC_range[0] < 0. && tracer_cat->ra(i) > par::pi) ? tracer_cat->ra(i)-2*par::pi : tracer_cat->ra(i), 
  tracer_cat->xx(i), tracer_cat->yy(i), tracer_cat->zz(i)};
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
  shared_ptr<catalogue::Object> cell_obj = make_shared<catalogue::Halo>(cell_values[0], cell_values[1], cell_values[2],
  phi, theta, z_cells_centres[i], Divergence[i][j][k]/n_rec);
  Divergence_cell.add_object(cell_obj);
  n_cell++;
  }
  }
  }

  shared_ptr<Catalogue> div_cat = make_shared<catalogue::Catalogue> (catalogue::Catalogue(move(Divergence_cell)));
  vector<Var> variables = {Var::_X_, Var::_Y_, Var::_Z_};
  CatalogueChainMesh ChainMesh_div = CatalogueChainMesh(variables, cellsize, div_cat);

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
  coord_locmin[5].push_back(cosmology->Redshift(r));
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
  vector<double> temp_vect = {cosmology->D_C(z_cells_edges[z_cells_edges.size()-1]-z_cells_edges[coord_locmin[6][num_voids]]), 
  cosmology->D_C(z_cells_edges[coord_locmin[6][num_voids]]-z_cells_edges[0]), 
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
  */
}



//********************************** VERSION 1.0 **********************************//

void cbl::catalogue::Catalogue::clean_void_catalogue (const bool initial_radius, const std::vector<double> delta_r, const double threshold, bool rescale, const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const double ratio, const bool checkoverlap, const Var ol_criterion)
{
  const double start_time = omp_get_wtime();
    
  // ---------------------------------------------------- //
  // ---------------- Cleaning Procedure ---------------- //
  // ---------------------------------------------------- //

  cout << endl;
  coutCBL << par::col_blue << "> > > > > >>>>>>>>>> CLEANING PROCEDURE STARTED  <<<<<<<<< < < < < <" << par::col_default << endl << endl;
  coutCBL << "Voids in the initial Catalogue: " << nObjects() << endl << endl;
  if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 

  coutCBL << par::col_green << " --- Removing spurious voids --- " << par::col_default << endl << endl;
  coutCBL << "Removed voids: " << endl;

  double spurious_time = omp_get_wtime();
  if (initial_radius) {
    vector<bool> remove(nObjects(), false);
    for (size_t i=0; i<nObjects(); i++) 
      if (radius(i) < delta_r[0] || delta_r[1] < radius(i)) remove[i] = true;

    remove_objects(remove);
    cout << "\t r_min-r_max criterion: " << count(remove.begin(), remove.end(), true) << endl;
  }
  
  compute_centralDensity(tracers_catalogue, ChM, threshold);
  coutCBL << "Voids in the Catalogue: " << nObjects() << endl;
  if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
  coutCBL << "Time spent by the spurious voids-checking procedure: " << fixed << setprecision(3) << omp_get_wtime()-spurious_time << " seconds" << endl;


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
        vector<double> temp_val = {3.*radius(j), xx(j)-min_X, max_X-xx(j), yy(j)-min_Y, max_Y-yy(j), zz(j)-min_Z, max_Z-zz(j), delta_r[1]};
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
          bad_rescale[j] = true;
          continue;
        }

        std::sort(radii.begin(), radii.end());
        vector<vector<double>> data(2, vector<double>(radii.size())); //distances, density contrast 
        for (size_t i=0; i<radii.size()-1; i++) {
          data[0][i] = (radii[i]+radii[i+1])*0.5;
          data[1][i] = (i+1)/(volume_sphere(data[0][i])*density);
        }
        data[0][radii.size()-1] = 2*data[0][radii.size()-2] - data[0][radii.size()-3];
        data[1][radii.size()-1] = radii.size()/(volume_sphere(data[0][radii.size()-1])*density);
        
        unsigned int N = std::round(radii.size()*0.25);
        bool expand = false;

        if (N>1) {
          while (!(data[1][N-2] < threshold && data[1][N-1] > threshold)) 
            if (N>2) 
              N--;
            else {
              expand = true;
              break;
            }
          if (expand) {
            N = std::round(radii.size()*0.25);
            while (!(data[1][N-2] < threshold && data[1][N-1] > threshold)) 
              if (N<radii.size()-2) N++;
              else break;
          }
        }
      
        double new_radius;
        if (N<=2 or N>=radii.size()-2) {
          remove[j] = true;
          bad_rescale[j] = true;
          continue;
        }
        else
          new_radius = interpolated(threshold, {data[1][N-1], data[1][N]}, {data[0][N-1], data[0][N]}, "Linear");

        if (new_radius > 0) set_var(j, Var::_Radius_, new_radius);
        else {
          remove[j] = true;
          bad_rescale[j] = true;
          continue;
        }

	// -----------------------

        N = std::round(radii.size()*0.25);
        expand = false;

        if (N>1) {
          while (!(data[1][N-2] < 0.5 && data[1][N-1] > 0.5))
	    if (N>2)
	      N--;
	    else {
	      expand = true;
	      break;
	    }
          if (expand) {
            N = std::round(radii.size()*0.25);
            while (!(data[1][N-2] < 0.5 && data[1][N-1] > 0.5))
              if (N<radii.size()-2)
                N++;
              else
                break;
          }
        }
      
        double new_generic;
        if (N<=2 or N>=radii.size()-2) {
          remove[j] = true;
          bad_rescale[j] = true;
          continue;
        }
        else
	  new_generic = interpolated(threshold, {data[1][N-1], data[1][N]}, {data[0][N-1], data[0][N]}, "Linear");

        if (new_generic > 0) set_var(j, Var::_Generic_, new_generic);
        else {
          remove[j] = true;
          bad_rescale[j] = true;
          continue;
        }
      }
    }

    remove_objects(remove);
    coutCBL << "Removed voids:" << endl;
    cout << "\t Bad rescaled: " << count(bad_rescale.begin(), bad_rescale.end(), true) << endl;
    
    vector<bool> remove_outofrange(nObjects(), false);
    for (size_t ii = 0; ii<nObjects(); ii++) 
      if (radius(ii) < delta_r[0] || delta_r[1] < radius(ii)) 
	remove_outofrange[ii] = true;
      
    remove_objects(remove_outofrange);
    cout << "\t Out of range ["+conv(delta_r[0],par::fDP2)+","+conv(delta_r[1],par::fDP2)+"]: " << count(remove_outofrange.begin(), remove_outofrange.end(), true) << endl;

    compute_densityContrast(tracers_catalogue, ChM, ratio);

    coutCBL << "Time spent by the rescaling procedure: " << omp_get_wtime()-rescaling_time << " seconds" << endl;
    
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
    vector<double> criteriumOrder(nObjects());
    string criterium;

    if (ol_criterion == Var::_CentralDensity_) {
      criteriumOrder = var(Var::_CentralDensity_);
      criterium = "central density";
    }
    else if (ol_criterion == Var::_DensityContrast_) {
      criteriumOrder = var(Var::_DensityContrast_);
      criterium = "density contrast";
    }

    else ErrorCBL("allowed overlap criteria are '_CentralDensity_' or '_DensityContrast_' .", "Catalogue", "VoidCatalogue.cpp");

    std::vector<int> indices(nObjects(), 0);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](int A, int B) 
	      -> bool {return (ol_criterion == Var::_CentralDensity_) ? criteriumOrder[A] < criteriumOrder[B] : criteriumOrder[A] > criteriumOrder[B];});
    
    Order(indices);
    coutCBL << "Catalogue ordered according to the " << criterium << endl << endl;
    coutCBL << "* * * Generating ChainMesh for void centres * * *" << endl;

    chainmesh::ChainMesh3D ChM_voids(Min(Var::_Radius_), var(Var::_X_), var(Var::_Y_), var(Var::_Z_), 2*Max(Var::_Radius_));

    for (size_t i=0; i<nObjects(); i++) {
      if (!remove[i]) {
        vector<long> close = ChM_voids.close_objects(coordinate(i));
        std::sort(close.begin(), close.end());
        for (auto && j : close) {
          if (!remove[j]) {
            double distance = Catalogue::distance(i, catalogue_object(j));
            if (distance < radius(i)+radius(j) && (long)i!=j) {
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
 
    cout << endl;
    coutCBL << "Voids removed to avoid overlap: " << count(remove.begin(), remove.end(), true) << endl;
    coutCBL << "Time spent by the overlap-checking procedure: " << omp_get_wtime()-ol_time << " seconds" << endl;
  }
  
  cout << endl;
  coutCBL << "Voids in the final Catalogue: " << nObjects() << endl;
  if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp");
  
  cout << endl;
  coutCBL << "Total time spent: " << omp_get_wtime()-start_time << " seconds \n" << endl;

}

//********************************** VERSION 2.0 **********************************//

void cbl::catalogue::Catalogue::clean_void_catalogue (const std::vector<double> par_numdensity, const bool initial_radius, const std::vector<double> delta_r, const double threshold, bool rescale, const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const double ratio, const bool checkoverlap, const Var ol_criterion)
{
  const double start_time = omp_get_wtime();
    
  // ---------------------------------------------------- //
  // ---------------- Cleaning Procedure ---------------- //
  // ---------------------------------------------------- //

  cout << endl;
  coutCBL << par::col_blue << "> > > > > >>>>>>>>>> CLEANING PROCEDURE STARTED  <<<<<<<<< < < < < <" << par::col_default << endl << endl;
  coutCBL << "Voids in the initial Catalogue: " << nObjects() << endl << endl;
  if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 

  coutCBL << par::col_green << " --- Removing spurious voids --- " << par::col_default << endl << endl;
  coutCBL << "Removed voids: " << endl;
  
  if (initial_radius) {
    vector<bool> remove(nObjects(), false);
    for (size_t i=0; i<nObjects(); i++) 
      if (radius(i) < delta_r[0] || delta_r[1] < radius(i)) remove[i] = true;

    remove_objects(remove);
    cout << "\t r_min-r_max criterion: " << count(remove.begin(), remove.end(), true) << endl;
  }

  compute_centralDensity(tracers_catalogue, ChM, par_numdensity, threshold);
  coutCBL << "Voids in the Catalogue: " << nObjects() << endl;
  if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 

  // ---------------------------------------------------- //
  // ----------------- Radius Rescaling ----------------- //
  // ---------------------------------------------------- //
  
  if (rescale) {
   
    cout << endl;
    coutCBL << par::col_green << " --- Rescaling radii --- \n" << par::col_default << endl;

    double rescaling_time = omp_get_wtime();
    
    vector<bool> remove(nObjects(), false), bad_rescale(nObjects(), false);

#pragma omp parallel num_threads(omp_get_max_threads())
    {     
#pragma omp for ordered schedule(dynamic)
      for (size_t j=0; j<nObjects(); j++) {
        double value = 2*radius(j);
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
          bad_rescale[j] = true;
          continue;
        }

        double zz = redshift(j);
        double density = 0.;
        for (size_t N=par_numdensity.size(); N-->0;) density += par_numdensity[par_numdensity.size()-1-N]*pow(zz,N);

        std::sort(radii.begin(), radii.end());

        vector<vector<double>> data(2, vector<double>(radii.size())); //distances, density contrast 
        for (size_t i=0; i<radii.size()-1; i++) {
          data[0][i] = (radii[i]+radii[i+1])*0.5;
          data[1][i] = (i+1)/(volume_sphere(data[0][i])*density);
        }
        data[0][radii.size()-1] = 2*data[0][radii.size()-2] - data[0][radii.size()-3];
        data[1][radii.size()-1] = radii.size()/(volume_sphere(data[0][radii.size()-1])*density);
        
        unsigned int N = std::round(radii.size()*0.25);
        bool expand = false;

        if (N>1) {
          while (!(data[1][N-2] < threshold && data[1][N-1] > threshold))
            if (N>2) N--;
            else {
              expand = true;
              break;
            }
          if (expand) {
            N = std::round(radii.size()*0.25);
            while (!(data[1][N-2] < threshold && data[1][N-1] > threshold)) 
              if (N<radii.size()-2) N++;
              else break;
          }
        }
      
        double new_radius;
        if (N<=2 or N>=radii.size()-2) {
          remove[j] = true;
          bad_rescale[j] = true;
          continue;
        }
        else
          new_radius = interpolated(threshold, {data[1][N-1], data[1][N]}, {data[0][N-1], data[0][N]}, "Linear");

        if (new_radius > 0) set_var(j, Var::_Radius_, new_radius);
        else {
          remove[j] = true;
          bad_rescale[j] = true;
          continue;
        }

	// -----------------------

        N = std::round(radii.size()*0.25);
        expand = false;

        if (N>1) {
          while (!(data[1][N-2] < 0.5 && data[1][N-1] > 0.5))
	    if (N>2)
	      N--;
	    else {
	      expand = true;
	      break;
	    }
          if (expand) {
            N = std::round(radii.size()*0.25);
            while (!(data[1][N-2] < 0.5 && data[1][N-1] > 0.5))
              if (N<radii.size()-2)
                N++;
              else
                break;
          }
        }
      
        double new_generic;
        if (N<=2 or N>=radii.size()-2) {
          remove[j] = true;
          bad_rescale[j] = true;
          continue;
        }
        else
	  new_generic = interpolated(threshold, {data[1][N-1], data[1][N]}, {data[0][N-1], data[0][N]}, "Linear");

        if (new_generic > 0) set_var(j, Var::_Generic_, new_generic);
        else {
          remove[j] = true;
          bad_rescale[j] = true;
          continue;
        }
      }
    }

    remove_objects(remove);
    coutCBL << "Removed voids:" << endl;
    cout << "\t Bad rescaled: " << count(bad_rescale.begin(), bad_rescale.end(), true) << endl;

    vector<bool> remove_outofrange(nObjects(), false);
    for (size_t ii = 0; ii<nObjects(); ii++) 
      if (radius(ii) < delta_r[0] || delta_r[1] < radius(ii)) 
	remove_outofrange[ii] = true;
      
    remove_objects(remove_outofrange);
    cout << "\t Out of range ["+conv(delta_r[0],par::fDP2)+","+conv(delta_r[1],par::fDP2)+"]: " << count(remove_outofrange.begin(), remove_outofrange.end(), true) << endl;
    
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

    vector<double> criteriumOrder(nObjects());
    string criterium;

    if (ol_criterion == Var::_CentralDensity_) {
      criteriumOrder = var(Var::_CentralDensity_);
      criterium = "central density";
    }
    else if (ol_criterion == Var::_DensityContrast_) {
      criteriumOrder = var(Var::_DensityContrast_);
      criterium = "density contrast";
    }

    else ErrorCBL("allowed overlap criteria are '_CentralDensity_' or '_DensityContrast_' .", "Catalogue", "VoidCatalogue.cpp");

    std::vector<int> indices(nObjects(), 0);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](int A, int B) 
	      -> bool {return (ol_criterion == Var::_CentralDensity_) ? criteriumOrder[A] < criteriumOrder[B] : criteriumOrder[A] > criteriumOrder[B];});
    
    Order(indices);
    coutCBL << "Catalogue ordered according to the " << criterium << endl << endl;
    coutCBL << "* * * Generating ChainMesh for void centres * * *" << endl;

    chainmesh::ChainMesh3D ChM_voids(Min(Var::_Radius_), var(Var::_X_), var(Var::_Y_), var(Var::_Z_), 2*Max(Var::_Radius_));
    
    for (size_t i=0; i<nObjects(); i++) {
      if (!remove[i]) {
        vector<long> close = ChM_voids.close_objects(coordinate(i));
        std::sort(close.begin(), close.end());
        for (auto &&j : close) {
          if (!remove[j]) {
            double distance = Catalogue::distance(i, catalogue_object(j));
            if (distance<radius(i)+radius(j) && (long)i!=j) {
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
 
    cout << endl;
    coutCBL << "Voids removed to avoid overlap: " << count(remove.begin(), remove.end(), true) << endl;
    coutCBL << "Time spent by the overlap-checking procedure: " << omp_get_wtime()-ol_time << " seconds" << endl;
  }
  
  cout << endl;
  coutCBL << "Voids in the final Catalogue: " << nObjects() << endl;
  if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp");
  
  cout << endl;
  coutCBL << "Total time spent: " << omp_get_wtime()-start_time << " seconds \n" << endl;

}


//********************************** VERSION 3.0 **********************************//

void cbl::catalogue::Catalogue::clean_void_catalogue (const std::vector<std::vector<double>> data_numdensity, const std::string method_interpolation, const bool initial_radius, const std::vector<double> delta_r, const double threshold, bool rescale, const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const double ratio, const bool checkoverlap, const Var ol_criterion)
{

  const double start_time = omp_get_wtime();
    
  // ---------------------------------------------------- //
  // ---------------- Cleaning Procedure ---------------- //
  // ---------------------------------------------------- //

  cout << endl;
  coutCBL << par::col_blue << "> > > > > >>>>>>>>>> CLEANING PROCEDURE STARTED  <<<<<<<<< < < < < <" << par::col_default << endl << endl;
  coutCBL << "Voids in the initial Catalogue: " << nObjects() << endl << endl;
  if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 

  coutCBL << par::col_green << " --- Removing spurious voids --- " << par::col_default << endl << endl;
  coutCBL << "Removed voids: " << endl;
  
  if (initial_radius) {
    vector<bool> remove(nObjects(), false);
    for (size_t i=0; i<nObjects(); i++) 
      if (radius(i) < delta_r[0] || delta_r[1] < radius(i)) remove[i] = true;

    remove_objects(remove);
    cout << "\t r_min-r_max criterion: " << count(remove.begin(), remove.end(), true) << endl;
  }

  compute_centralDensity(tracers_catalogue, ChM, data_numdensity, method_interpolation, threshold);
  coutCBL << "Voids in the Catalogue: " << nObjects() << endl;
  if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
  
  // ---------------------------------------------------- //
  // ----------------- Radius Rescaling ----------------- //
  // ---------------------------------------------------- //
  
  if (rescale) {
   
    cout << endl;
    coutCBL << par::col_green << " --- Rescaling radii --- \n" << par::col_default << endl;

    double rescaling_time = omp_get_wtime();
    
    vector<bool> remove(nObjects(), false), bad_rescale(nObjects(), false);

#pragma omp parallel num_threads(omp_get_max_threads())
    {     
#pragma omp for ordered schedule(dynamic)
      for (size_t j=0; j<nObjects(); j++) {
        double value = 2*radius(j);
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
          bad_rescale[j] = true;
          continue;
        }

        double zz = redshift(j);
        double density = interpolated(zz, data_numdensity[0], data_numdensity[1], method_interpolation);

        std::sort(radii.begin(), radii.end());

        vector<vector<double>> data(2, vector<double>(radii.size())); //distances, density contrast 
        for (size_t i=0; i<radii.size()-1; i++) {
          data[0][i] = (radii[i]+radii[i+1])*0.5;
          data[1][i] = (i+1)/(volume_sphere(data[0][i])*density);
        }
        data[0][radii.size()-1] = 2*data[0][radii.size()-2] - data[0][radii.size()-3];
        data[1][radii.size()-1] = radii.size()/(volume_sphere(data[0][radii.size()-1])*density);
        
        unsigned int N = std::round(radii.size()*0.25);
        bool expand = false;

        if (N>1) {
          while (!(data[1][N-2] < threshold && data[1][N-1] > threshold))
            if (N>2) N--;
            else {
              expand = true;
              break;
            }
          if (expand) {
            N = std::round(radii.size()*0.25);
            while (!(data[1][N-2] < threshold && data[1][N-1] > threshold)) 
              if (N<radii.size()-2) N++;
              else break;
          }
        }
      
        double new_radius;
        if (N<=2 or N>=radii.size()-2) {
          remove[j] = true;
          bad_rescale[j] = true;
          continue;
        }
        else
          new_radius = interpolated(threshold, {data[1][N-1], data[1][N]}, {data[0][N-1], data[0][N]}, "Linear");

        if (new_radius > 0) set_var(j, Var::_Radius_, new_radius);
        else {
          remove[j] = true;
          bad_rescale[j] = true;
          continue;
        }

	// -----------------------

        N = std::round(radii.size()*0.25);
        expand = false;

        if (N>1) {
          while (!(data[1][N-2] < 0.5 && data[1][N-1] > 0.5))
	    if (N>2)
	      N--;
	    else {
	      expand = true;
	      break;
	    }
          if (expand) {
            N = std::round(radii.size()*0.25);
            while (!(data[1][N-2] < 0.5 && data[1][N-1] > 0.5))
              if (N<radii.size()-2)
                N++;
              else
                break;
          }
        }
      
        double new_generic;
        if (N<=2 or N>=radii.size()-2) {
          remove[j] = true;
          bad_rescale[j] = true;
          continue;
        }
        else
	  new_generic = interpolated(threshold, {data[1][N-1], data[1][N]}, {data[0][N-1], data[0][N]}, "Linear");

        if (new_generic > 0) set_var(j, Var::_Generic_, new_generic);
        else {
          remove[j] = true;
          bad_rescale[j] = true;
          continue;
        }
      }
    }

    remove_objects(remove);
    coutCBL << "Removed voids:" << endl;
    cout << "\t Bad rescaled: " << count(bad_rescale.begin(), bad_rescale.end(), true) << endl;
    
    vector<bool> remove_outofrange(nObjects(), false);
    for (size_t ii = 0; ii<nObjects(); ii++) 
      if (radius(ii) < delta_r[0] || delta_r[1] < radius(ii)) 
	remove_outofrange[ii] = true;
      
    remove_objects(remove_outofrange);
    cout << "\t Out of range ["+conv(delta_r[0],par::fDP2)+","+conv(delta_r[1],par::fDP2)+"]: " << count(remove_outofrange.begin(), remove_outofrange.end(), true) << endl;
    
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

    vector<double> criteriumOrder(nObjects());
    string criterium;

    if (ol_criterion == Var::_CentralDensity_) {
      criteriumOrder = var(Var::_CentralDensity_);
      criterium = "central density";
    }
    else if (ol_criterion == Var::_DensityContrast_) {
      criteriumOrder = var(Var::_DensityContrast_);
      criterium = "density contrast";
    }

    else ErrorCBL("allowed overlap criteria are '_CentralDensity_' or '_DensityContrast_' .", "Catalogue", "VoidCatalogue.cpp");

    std::vector<int> indices(nObjects(), 0);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](int A, int B) 
	      -> bool {return (ol_criterion == Var::_CentralDensity_) ? criteriumOrder[A] < criteriumOrder[B] : criteriumOrder[A] > criteriumOrder[B];});
    
    Order(indices);
    coutCBL << "Catalogue ordered according to the " << criterium << endl << endl;
    coutCBL << "* * * Generating ChainMesh for void centres * * *" << endl;

    chainmesh::ChainMesh3D ChM_voids(Min(Var::_Radius_), var(Var::_X_), var(Var::_Y_), var(Var::_Z_), 2*Max(Var::_Radius_));
    
    for (size_t i=0; i<nObjects(); i++) {
      if (!remove[i]) {
        vector<long> close = ChM_voids.close_objects(coordinate(i));
        std::sort(close.begin(), close.end());
        for (auto &&j : close) {
          if (!remove[j]) {
            double distance = Catalogue::distance(i, catalogue_object(j));
            if (distance < radius(i)+radius(j) && (long)i!=j) {
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

    cout << endl;
    coutCBL << "Voids removed to avoid overlap: " << count(remove.begin(), remove.end(), true) << endl;
    coutCBL << "Time spent by the overlap-checking procedure: " << omp_get_wtime()-ol_time << " seconds" << endl;
  }
  
  cout << endl;
  coutCBL << "Voids in the final Catalogue: " << nObjects() << endl;
  if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp");
  
  cout << endl;
  coutCBL << "Total time spent: " << omp_get_wtime()-start_time << " seconds \n" << endl;
  
}


// ============================================================================


void cbl::catalogue::Catalogue::compute_centralDensity (const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const double threshold)
{
  const double density = tracers_catalogue->numdensity();
  const double mps = tracers_catalogue->mps();
  vector<bool> remove(nObjects(), false);

  for (size_t j=0; j<nObjects(); j++) {
    ChM.get_searching_region(2*mps);
    int cont=0;
    vector<long> close = ChM.close_objects(coordinate(j));
    for (auto&& k : close) {
      double distance = catalogue::Catalogue::distance(j, tracers_catalogue->catalogue_object(k));
      if (distance <= 2*mps) cont++;
    }
    
    if (cont/(volume_sphere(2*mps)*density) > threshold) remove[j]=true;
    else
      set_var(j, Var::_CentralDensity_, cont/(volume_sphere(2*mps)*density));
  }
  remove_objects(remove);
  cout << "\t High central density: " << std::count(remove.begin(), remove.end(), true) << endl;
  
}


// ============================================================================


void cbl::catalogue::Catalogue::compute_centralDensity (const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const std::vector<double> par_numdensity, const double threshold)
{
  vector<bool> remove(nObjects(), false);
  for (size_t j=0; j<nObjects(); j++) {

    double zz = m_object[j]->redshift();
    double density = 0.;
    for (size_t N=par_numdensity.size(); N-->0;) density += par_numdensity[par_numdensity.size()-1-N]*pow(zz,N);
    double mps = pow(density, -1./3.);
    
    ChM.get_searching_region(2.*mps);
    int cont = 0;
    vector<long> close = ChM.close_objects(coordinate(j));
    for (auto&& k : close) {
      double distance = catalogue::Catalogue::distance(j, tracers_catalogue->catalogue_object(k));
      if (distance <= 2.*mps) cont++;
    }
    
    if (cont/(volume_sphere(2*mps)*density) > threshold) remove[j]=true;
    else
      set_var(j, Var::_CentralDensity_, cont/(volume_sphere(2*mps)*density));
  }
  remove_objects(remove);
  cout << "\t High central density: " << std::count(remove.begin(), remove.end(), true) << endl;
}


// ============================================================================


void cbl::catalogue::Catalogue::compute_centralDensity (const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const std::vector<std::vector<double>> data_numdensity, const std::string method_interpolation, const double threshold)
{
  vector<bool> remove(nObjects(), false);
  for (size_t j=0; j<nObjects(); j++) {

    double zz = m_object[j]->redshift();
    double density = interpolated(zz, data_numdensity[0], data_numdensity[1], method_interpolation);
    double mps = pow(density, -1./3.);
    
    ChM.get_searching_region(2.*mps);
    int cont = 0;
    vector<long> close = ChM.close_objects(coordinate(j));
    for (auto&& k : close) {
      double distance = catalogue::Catalogue::distance(j, tracers_catalogue->catalogue_object(k));
      if (distance <= 2.*mps) cont++;
    }
    
    if (cont/(volume_sphere(2*mps)*density) > threshold) remove[j]=true;
    else
      set_var(j, Var::_CentralDensity_, cont/(volume_sphere(2*mps)*density));
  }
  remove_objects(remove);
  cout << "\t High central density: " << std::count(remove.begin(), remove.end(), true) << endl;
  
}


// ============================================================================


void cbl::catalogue::Catalogue::compute_densityContrast (const shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const double ratio) {

  if (ratio==0. or ratio==1. or ratio>3.)
    ErrorCBL ("Please choose an appropriate value for the parameter ratio!", "compute_densityContrast", "VoidCatalogue.cpp");
  
  vector<bool> remove(nObjects(), false), void_voids(nObjects(), false), cloud_in_void(nObjects(), false);

#pragma omp parallel num_threads(omp_get_max_threads())
  {
#pragma omp for ordered schedule(static)
    for (size_t j = 0; j<nObjects(); j++) {

      double outer_radius = std::max(radius(j)*ratio, radius(j));
      double inner_radius = std::min(radius(j)*ratio, radius(j));

#pragma omp ordered
      ChM.get_searching_region(outer_radius);
      vector<long> close = ChM.close_objects(coordinate(j));
      
      //compute distances between the void and the surrounding particles
      vector<double> distances;
      for (auto&& k : close) {
        double distance = Catalogue::distance(j, tracers_catalogue->catalogue_object(k));
        if (distance <= outer_radius) distances.emplace_back(distance);
      }

      std::sort(distances.begin(), distances.end());
      
      if (distances.size() > 1) {
      
        //FIND CENTRAL DENSITY

        int NN_in = 0;
        int NN_out = 0;
        while (distances[NN_in]<inner_radius && NN_in<(int)distances.size()-1) NN_in++; // check-1
        while (distances[NN_out]<outer_radius && NN_out<(int)distances.size()-1) NN_out++;

        double delta_in = NN_in/volume_sphere(inner_radius);
        double delta_out = NN_out/volume_sphere(outer_radius);

        if (NN_in>0) {
          double compare = delta_out/delta_in;
          if (compare==0.) {
            remove[j] = true;
            void_voids[j] = true;
          }
          else if (compare>=1.)
            set_var(j, Var::_DensityContrast_, compare);

          else if (compare<1.) {
            remove[j] = true;
            cloud_in_void[j] = true;
          }
        }
        else set_var(j, Var::_DensityContrast_, 9999.);

      }
    }
  }
  remove_objects(remove);

  cout << "\t Cloud-in-void: " << std::count(cloud_in_void.begin(), cloud_in_void.end(), true) << endl;
  cout << "\t Empty voids: " << std::count(cloud_in_void.begin(), cloud_in_void.end(), true) << endl;

}
