/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo  *
 *  federico.marulli3@unibo.it                                      *
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
 *  @file
 *  Measure/CovarianceMatrix/CovarianceMatrix_TwoPointCorrelation_multipoles.cpp
 *
 *  @brief Methods of the class CovarianceMatrix_TwoPointCorrelation_multipoles used to
 *  measure the covariance of the monopole of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  CovarianceMatrix_TwoPointCorrelation_multipoles used to measure the covariance of the
 *  monopole of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "EigenWrapper.h"
#include "CovarianceMatrix_TwoPointCorrelation_multipoles.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace chainmesh;
using namespace measure;
using namespace covmat;


// ============================================================================================


void cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation_multipoles::m_compute_RR(const bool tcount)
{
  coutCBL << "I'm computing RR pairs" << endl;
  // timer
  time_t start; time (&start);
  int dp = cout.precision();
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);

  // number of objects in the first catalogue
  int nObj = m_random->nObjects();

  // pointer to the second catalogue (get from the chain mesh)
  shared_ptr<Catalogue> cat2 = m_chainMesh_rMax.catalogue();

  // factor used by the timer
  float fact_count = 100./nObj;
  (void)fact_count;

  // thread number
  int tid = 0;

  size_t nBins_mu = 30;
  double delta_mu_bin = double(nBins_mu);

  m_random_random.resize(m_nBins, nBins_mu);
  
#pragma omp parallel num_threads(omp_get_max_threads()) private(tid)
  {
    tid = omp_get_thread_num();

    std::vector<std::vector<double>> random_pairs_thread(m_nBins, std::vector<double>(nBins_mu, 0));
    // parallelized loop
#pragma omp for schedule(dynamic)
    for (int i=0; i<nObj; ++i) {

      Vector4D pos_i = m_random->eigen_coordinate(i);
      double ww_i = m_random->weight(i);

      vector<long int> close_objects = m_chainMesh_rMax.close_objects({pos_i[0], pos_i[1], pos_i[2]}, -1);

      for (auto &&j : close_objects) {
	Vector4D pos_j = cat2->eigen_coordinate(j);
	Vector4D sep = pos_i-pos_j;
	double rr = sep.norm();

	if (rr>m_rMax || rr==0)
	  continue;

	double ww_j = cat2->weight(j);
	double mu = fabs(m_random->dc(i)-cat2->dc(j))/rr;

	int rbin = static_cast<int>((rr-m_rMin)/m_binSize);
	int mubin = static_cast<int>(mu*delta_mu_bin);
	random_pairs_thread[rbin][mubin] += ww_i*ww_j;
      }

      // estimate the computational time and update the time count
      time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
      if (tcount && tid==0) { coutCBL << "\r" << float(i)*fact_count << "% completed (" << diff_temp << " seconds)\r"; cout.flush(); }
      if (i==int(nObj*0.25)) coutCBL << ".............25% completed" << endl;
      if (i==int(nObj*0.5)) coutCBL << ".............50% completed" << endl;
      if (i==int(nObj*0.75)) coutCBL << ".............75% completed"<< endl;
    }
#pragma omp critical
    {
      for (size_t i=0; i<m_nBins; ++i) 
	for (size_t j=0; j<nBins_mu; ++j) 
	  m_random_random(i, j) += random_pairs_thread[i][j]*delta_mu_bin;
    }
  }

  // show the time spent by the method
  time_t end; time (&end);
  double diff = difftime(end, start);
  if (tid==0) {
    if (diff<60) coutCBL << "   time spent to count the pairs: " << diff << " seconds" << endl ;
    else if (diff<3600) coutCBL << "   time spent to count the pairs: " << diff/60 << " minutes" << endl ;
    else coutCBL << "   time spent to count the pairs: " << diff/3600 << " hours" << endl ;
  }
  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(dp);
  coutCBL << "Done!" << endl;
}


// ============================================================================================


void cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation_multipoles::m_compute_C4(const bool tcount)
{
  coutCBL << "I'm computing the C4 term" << endl;

  // timer
  time_t start; time (&start);
  int dp = cout.precision();
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);

  vector<function<double(double)>> leg_pols(3);
  leg_pols[0] = [&] (const double mu) {(void)mu; return 1.;};
  leg_pols[1] = [&] (const double mu) {return 0.5 * (3*pow(mu,2)-1);};
  leg_pols[2] = [&] (const double mu) {return 0.125 * (35 * pow(mu, 4) - 30 * pow(mu, 2) +3) ;};

  // number of objects in the first catalogue
  int nObj = m_random->nObjects();

  // pointer to the second catalogue (get from the chain mesh)
  shared_ptr<Catalogue> cat2 = m_chainMesh_rMax.catalogue();

  // pointer to the second catalogue (get from the chain mesh)
  shared_ptr<Catalogue> cat3 = m_chainMesh_rCut.catalogue();

  // factor used by the timer
  float fact_count = 100./nObj;

  // thread number
  int tid = 0;
  
  std::vector<std::vector<double>> matrix(3*m_nBins, std::vector<double>(3*m_nBins, 0));

  double fraction = 0.1;
  cbl::random::UniformRandomNumbers rand(0., 1., 14231);

  #pragma omp parallel num_threads(omp_get_max_threads()) private(tid)
    {
      tid = omp_get_thread_num();

      std::vector<std::vector<double>> matrix_thread(3*m_nBins, std::vector<double>(3*m_nBins, 0));

      // parallelized loop
  #pragma omp for schedule(dynamic)
      for (int i=0; i<nObj; ++i) {

        Vector4D pos_i = m_random->eigen_coordinate(i);
        double ww_i = m_random->weight(i);

        vector<long int> close_objects_j = m_chainMesh_rMax.close_objects({pos_i[0], pos_i[1], pos_i[2]}, -1);
        size_t nEff_j = 0;
        vector<size_t> list_j (close_objects_j.size());
        vector<size_t> bin_r_j (close_objects_j.size());
        vector<size_t> bin_mu_j (close_objects_j.size());
	vector<double> mu_ij (close_objects_j.size());

        vector<long int> close_objects_k = m_chainMesh_rCut.close_objects({pos_i[0], pos_i[1], pos_i[2]}, -1);

        for (auto &&j : close_objects_j) {
          Vector4D pos_j = cat2->eigen_coordinate(j);
	  double rr_ij = (pos_i-pos_j).norm();

	  if ((rr_ij<m_rMin) || (rr_ij>m_rMax) || rr_ij==0)
	    continue;

	  if (rand() < fraction) {
	    mu_ij[nEff_j] = fabs(m_random->dc(i)-cat2->dc(j))/rr_ij;
	    list_j[nEff_j] = j;
	    bin_r_j[nEff_j] = static_cast<int>((rr_ij-m_rMin)/m_binSize);
	    bin_mu_j[nEff_j] = static_cast<int>(mu_ij[nEff_j]*30);
	    nEff_j +=1 ;
	  }
	}

	for (auto &&k : close_objects_k) {
	  Vector4D pos_k = cat3->eigen_coordinate(k);
	  double ww_k = cat3->weight(k);

	  double rr_ik = (pos_k-pos_i).norm();
	  double mu_ik = fabs(cat3->dc(k)-m_random->dc(i))/rr_ik;
	  double L2_ik = leg_pols[1](mu_ik);
	  double L4_ik = leg_pols[2](mu_ik);

	  if ((rr_ik<m_minSep) || (rr_ik>m_maxSep) || rr_ik==0) 
	    continue;

	  if (rand() < fraction) {
	    double xi0_ik = m_interpXi_0->operator()(rr_ik);
	    double xi2_ik = m_interpXi_2->operator()(rr_ik);
	    double xi4_ik = m_interpXi_4->operator()(rr_ik);
	    double xi_ik = xi0_ik + xi2_ik * L2_ik + xi4_ik * L4_ik; 

	    vector<long int> close_objects_l = m_chainMesh_rMax.close_objects({pos_k[0], pos_k[1], pos_k[2]}, -1);
	    for (auto &&l : close_objects_l) {
		Vector4D pos_l = cat2->eigen_coordinate(l);
		double ww_l = cat2->weight(l);

		double rr_kl = (pos_l-pos_k).norm();
		if ((rr_kl<m_rMin) || (rr_kl>m_rMax) || rr_kl==0)
		  continue;

		if (rand() < fraction) {
		  double mu_kl = fabs(cat3->dc(k)-cat2->dc(l))/rr_kl;
		  int bin2 = static_cast<int>((rr_kl-m_rMin)/m_binSize);
		  int mu2 = static_cast<int>(mu_kl*30);

		  for (size_t j=0; j<nEff_j; j++) {
		    double rr_jl = (pos_l - cat2->eigen_coordinate(list_j[j])).norm();
		    double mu_jl = fabs(cat2->dc(l)-cat2->dc(list_j[j]))/rr_jl;
		    double L2_jl = leg_pols[1](mu_jl);
		    double L4_jl = leg_pols[2](mu_jl);
		    double ww_j = cat2->weight(list_j[j]);

		    if ((rr_jl<m_minSep) || (rr_jl>m_maxSep) || rr_jl==0) 
		      continue;

		    double xi0_jl = m_interpXi_0->operator()(rr_jl);
		    double xi2_jl = m_interpXi_2->operator()(rr_jl);
		    double xi4_jl = m_interpXi_4->operator()(rr_jl);
		    double xi_jl = xi0_jl + xi2_jl * L2_jl + xi4_jl * L4_jl; 

		    double value = 2. * ww_i * ww_j * ww_k * ww_l * xi_ik * xi_jl;
		    value = value/(m_random_random(bin_r_j[j], bin_mu_j[j]) * m_random_random(bin2, mu2));
		    if (value!= value)
		      coutCBL << xi_ik << " " << xi_jl << endl;

		    for (size_t l1 = 0 ; l1<3; l1++)
		      for (size_t l2 = 0 ; l2<3; l2++) {
			double val = value * leg_pols[l1](mu_ij[j]) * leg_pols[l2](mu_kl) * (4*l1+1) * (4*l2+1);
			matrix_thread[m_nBins*l1+bin_r_j[j]][m_nBins*l2+bin2] += val;
			matrix_thread[m_nBins*l1+bin2][m_nBins*l2+bin_r_j[j]] += val;
		      }
		      
		  }
	      }
	    }
	  }
	}

        // estimate the computational time and update the time count
        time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
        if (tcount && tid==0) { coutCBL << "\r" << float(i)*fact_count << "% completed (" << diff_temp << " seconds)\r"; cout.flush(); }
        if (i==int(nObj*0.25)) coutCBL << ".............25% completed" << endl;
        if (i==int(nObj*0.5)) coutCBL << ".............50% completed" << endl;
        if (i==int(nObj*0.75)) coutCBL << ".............75% completed"<< endl;
      }
  #pragma omp critical
      {
        for (size_t i=0; i<3*m_nBins; ++i)
          for (size_t j=0; j<3*m_nBins; ++j)
            matrix[i][j] += matrix_thread[i][j]/pow(fraction, 3);
      }
    }

  m_C4 = cbl::wrapper::eigen::MatrixToEigen(matrix);

  // show the time spent by the method
  time_t end; time (&end);
  double diff = difftime(end, start);
  if (tid==0) {
    if (diff<60) coutCBL << "   time spent to count the pairs: " << diff << " seconds" << endl ;
    else if (diff<3600) coutCBL << "   time spent to count the pairs: " << diff/60 << " minutes" << endl ;
    else coutCBL << "   time spent to count the pairs: " << diff/3600 << " hours" << endl ;
  }
  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(dp);
  coutCBL << "Done!" << endl;
}


// ============================================================================================


void cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation_multipoles::m_compute_C3(const bool tcount)
{
  coutCBL << "I'm computing the C3 term" << endl;

  // timer
  time_t start; time (&start);
  int dp = cout.precision();
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);

  // number of objects in the first catalogue
  int nObj = m_random->nObjects();

  // pointer to the second catalogue (get from the chain mesh)
  shared_ptr<Catalogue> cat2 = m_chainMesh_rMax.catalogue();

  double fraction = 0.1;
  cbl::random::UniformRandomNumbers rand(0., 1., 14231);

  // factor used by the timer
  float fact_count = 100./nObj;

  // thread number
  int tid = 0;

  vector<function<double(double)>> leg_pols(3);
  leg_pols[0] = [&] (const double mu) {(void)mu; return 1.;};
  leg_pols[1] = [&] (const double mu) {return 0.5 * (3*pow(mu,2)-1);};
  leg_pols[2] = [&] (const double mu) {return 0.125 * (35 * pow(mu, 4) - 30 * pow(mu, 2) +3) ;};

  std::vector<std::vector<double>> matrix(3*m_nBins, std::vector<double>(3*m_nBins, 0));

  #pragma omp parallel num_threads(omp_get_max_threads()) private(tid)
    {
      tid = omp_get_thread_num();

      std::vector<std::vector<double>> matrix_thread(3*m_nBins, std::vector<double>(3*m_nBins, 0));

      // parallelized loop
#pragma omp for schedule(dynamic)
      for (int i=0; i<nObj; ++i) {
        Vector4D pos_i = m_random->eigen_coordinate(i);
        double ww_i = m_random->weight(i);

	vector<long int> close_objects = m_chainMesh_rMax.close_objects({pos_i[0], pos_i[1], pos_i[2]}, -1);

	for (auto &&j : close_objects) {
	  if (rand() < fraction) {
	    Vector4D pos_j = cat2->eigen_coordinate(j);
	    double ww_j = cat2->weight(j);

	    double rr_ij = (pos_i-pos_j).norm();
	    double mu_ij = fabs(m_random->dc(i)-cat2->dc(j))/rr_ij;

	    if ((rr_ij<m_rMin) || (rr_ij>m_rMax) || rr_ij==0)
	      continue;

	    int bin1 = static_cast<int>((rr_ij-m_rMin)/m_binSize);
	    int mu1 = static_cast<int>(mu_ij*30);

	    vector<long int> close_objects_j = m_chainMesh_rMax.close_objects({pos_j[0], pos_j[1], pos_j[2]}, -1);

	    for (auto &&k : close_objects_j) {
	      Vector4D pos_k = cat2->eigen_coordinate(k);
	      double ww_k = cat2->weight(k);

	      double rr_jk = (pos_k-pos_j).norm();

	      if ((rr_jk<m_rMin) || (rr_jk>m_rMax) || rr_jk==0)
		continue;

	      double mu_jk = fabs(cat2->dc(j)-cat2->dc(k))/rr_jk;
	      int bin2 = static_cast<int>((rr_jk-m_rMin)/m_binSize);
	      int mu2 = static_cast<int>(mu_jk*30);

	      double rr_ik = (pos_k-pos_i).norm();
	      double mu_ik = fabs(m_random->dc(i)-cat2->dc(k))/rr_ik;
	      double L2 = leg_pols[1](mu_ik);
	      double L4 = leg_pols[2](mu_ik);
	      double xi0 = m_interpXi_0->operator()(rr_ik);
	      double xi2 = m_interpXi_2->operator()(rr_ik);
	      double xi4 = m_interpXi_4->operator()(rr_ik);
	      double xi = (rr_ik>m_minSep && rr_ik<m_maxSep ) ? xi0 + L2*xi2 + L4*xi4 : 0.;
	      double value = 4 * ww_i * ww_j * ww_j * ww_k * xi ;

	      value /= (m_random_random(bin1, mu1) * m_random_random(bin2, mu2));

	      for (size_t l1 = 0 ; l1<3; l1++)
		for (size_t l2 = 0 ; l2<3; l2++) {
		  double val = value * leg_pols[l1](mu_ij) * leg_pols[l2](mu_jk) * (4*l1+1) * (4*l2+1);
		  matrix_thread[m_nBins*l1+bin1][m_nBins*l2+bin2] += val;
		  matrix_thread[m_nBins*l1+bin2][m_nBins*l2+bin1] += val;
		}

	    }
	  }
        }

        // estimate the computational time and update the time count
        time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
        if (tcount && tid==0) { coutCBL << "\r" << float(i)*fact_count << "% completed (" << diff_temp << " seconds)\r"; cout.flush(); }
        if (i==int(nObj*0.25)) coutCBL << ".............25% completed" << endl;
        if (i==int(nObj*0.5)) coutCBL << ".............50% completed" << endl;
        if (i==int(nObj*0.75)) coutCBL << ".............75% completed"<< endl;
      }
  #pragma omp critical
      {
        for (size_t i=0; i<3*m_nBins; ++i)
          for (size_t j=0; j<3*m_nBins; ++j)
            matrix[i][j] += matrix_thread[i][j]/fraction;
      }
    }

  m_C3 = cbl::wrapper::eigen::MatrixToEigen(matrix);

  // show the time spent by the method
  time_t end; time (&end);
  double diff = difftime(end, start);
  if (tid==0) {
    if (diff<60) coutCBL << "   time spent to count the pairs: " << diff << " seconds" << endl ;
    else if (diff<3600) coutCBL << "   time spent to count the pairs: " << diff/60 << " minutes" << endl ;
    else coutCBL << "   time spent to count the pairs: " << diff/3600 << " hours" << endl ;
  }
  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(dp);
  coutCBL << "Done!" << endl;
}


// ============================================================================================


void cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation_multipoles::m_compute_C2(const bool tcount)
{
  coutCBL << "I'm computing C2 term" << endl;
  // timer
  time_t start; time (&start);
  int dp = cout.precision();
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);

  // number of objects in the first catalogue
  int nObj = m_random->nObjects();

  // pointer to the second catalogue (get from the chain mesh)
  shared_ptr<Catalogue> cat2 = m_chainMesh_rMax.catalogue();

  // factor used by the timer
  float fact_count = 100./nObj;
  (void)fact_count;

  // thread number
  int tid = 0;

  double fraction = 1.1;
  cbl::random::UniformRandomNumbers rand(0., 1., 53324);

  std::vector<std::vector<double>> c2(3*m_nBins, std::vector<double>(3*m_nBins, 0));

#pragma omp parallel num_threads(omp_get_max_threads()) private(tid)
  {
    tid = omp_get_thread_num();

    std::vector<std::vector<double>> c2_thread(3*m_nBins, std::vector<double>(3*m_nBins, 0));

    // parallelized loop
#pragma omp for schedule(static, 2)
    for (int i=0; i<nObj; ++i) {
      Vector4D pos_i = m_random->eigen_coordinate(i);
      double ww_i = m_random->weight(i);

      vector<long int> close_objects = m_chainMesh_rMax.close_objects({pos_i[0], pos_i[1], pos_i[2]}, -1);

      for (auto &&j : close_objects) {
	if (rand() < fraction) {
	  Vector4D pos_j = cat2->eigen_coordinate(j);
	  double ww_j = cat2->weight(j);

	  double rr = (pos_i-pos_j).norm();
	  if ((rr<m_rMin) || (rr>m_rMax) || rr==0)
	    continue;

	  double mu = fabs(m_random->dc(i)-cat2->dc(j))/rr;
	  double L2 = 0.5 * (3 * mu*mu -1);
	  double L4 = 0.125 * ( 35 * pow(mu, 4) - 30 * pow(mu, 2) + 3);
	  vector<double> pl = {1, L2, L4};

	  double xi0 = m_interpXi_0->operator()(rr);
	  double xi2 = m_interpXi_2->operator()(rr);
	  double xi4 = m_interpXi_4->operator()(rr);
	  double xi = xi0 + L2*xi2 + L4*xi4;

	  int kk = static_cast<int>((rr-m_rMin)/m_binSize);
	  int mubin = static_cast<int>((mu)*30); 

	  double value =  ww_i*ww_j*ww_i*ww_j*(1+xi) / pow(m_random_random(kk, mubin), 2);

	  for (size_t l1=0; l1<3; l1++)
	    for (size_t l2=0; l2<3; l2++)
	      c2_thread[l1*m_nBins+kk][l2*m_nBins+kk] += 2 * (2*2*l1+1) * (2*2*l2+1) * value * pl[l1] * pl[l2];
	}
      }

      // estimate the computational time and update the time count
      time_t end_temp; time (&end_temp); double diff_temp = difftime(end_temp, start);
      if (tcount && tid==0) { coutCBL << "\r" << float(i)*fact_count << "% completed (" << diff_temp << " seconds)\r"; cout.flush(); }
      if (i==int(nObj*0.25)) coutCBL << ".............25% completed" << endl;
      if (i==int(nObj*0.5)) coutCBL << ".............50% completed" << endl;
      if (i==int(nObj*0.75)) coutCBL << ".............75% completed"<< endl;
    }
#pragma omp critical
    {
      for (size_t i=0; i<3*m_nBins; ++i)
	for (size_t j=0; j<3*m_nBins; ++j)
	  c2[i][j] += c2_thread[i][j]/fraction;
    }
  }
  m_C2 = cbl::wrapper::eigen::MatrixToEigen(c2);

  // show the time spent by the method
  time_t end; time (&end);
  double diff = difftime(end, start);
  if (tid==0) {
    if (diff<60) coutCBL << "   time spent to count the pairs: " << diff << " seconds" << endl ;
    else if (diff<3600) coutCBL << "   time spent to count the pairs: " << diff/60 << " minutes" << endl ;
    else coutCBL << "   time spent to count the pairs: " << diff/3600 << " hours" << endl ;
  }
  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(dp);
  coutCBL << "Done!" << endl;
}


// ============================================================================================


void cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation_multipoles::write_matrix (Eigen::MatrixXd matrix,
  const std::string dir,
  const std::string file,
  const int nDigits)
{
  makeDir(dir);
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  //fout << "### [1] bin1 # [2] bin2 # [3] matrix ### " << endl;

  for (int i=0; i<matrix.rows(); ++i){
    for (int j=0; j<matrix.cols(); ++j) {
      fout << setprecision(nDigits) << i << " " << j << " " << matrix(i, j) << endl;
    }
    fout << endl;
  }

  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl;
}


// ============================================================================================

void cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation_multipoles::m_set_catalogues(const catalogue::Catalogue random)
{
  m_random = std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(random));

  auto cat2 = std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(random));
  auto cat3 = std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(random));

  double rMAX = m_rMax*1.05;
  double cell_size = rMAX*0.1;
  m_chainMesh_rMax.set_par(cell_size, cat2, rMAX);

  rMAX = m_maxSep*1.05;
  cell_size = rMAX*0.1;
  m_chainMesh_rCut.set_par(cell_size, cat3, rMAX);
}

// ============================================================================================


cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation_multipoles::CovarianceMatrix_TwoPointCorrelation_multipoles (const catalogue::Catalogue random,
  const BinType binType,
  const double rMin,
  const double rMax,
  const int nBins,
  const glob::FuncGrid interpXi_0,
  const glob::FuncGrid interpXi_2,
  const glob::FuncGrid interpXi_4,
  const double minSeparation,
  const double maxSeparation)
{
  m_binType = binType;
  m_rMin = rMin;
  m_rMax = rMax;
  m_nBins = nBins;

  switch (m_binType) {
    case BinType::_linear_:
    {
      m_binSize = (m_rMax-m_rMin)/m_nBins;
      break;
    }

    case BinType::_logarithmic_:
    {
      double _rMin = log10(rMin);
      double _rMax = log10(rMax);
      m_binSize = (_rMax-_rMin)/m_nBins;
      break;
    }

    default:
      ErrorCBL("no such a bin type!", "constructor", "CovarianceMatrix_TwoPointCorrelation_multipoles.cpp");
  }

  m_interpXi_0 = std::make_shared<cbl::glob::FuncGrid>(interpXi_0);
  m_interpXi_2 = std::make_shared<cbl::glob::FuncGrid>(interpXi_2);
  m_interpXi_4 = std::make_shared<cbl::glob::FuncGrid>(interpXi_4);
  m_minSep = minSeparation;
  m_maxSep = maxSeparation;

  m_set_catalogues(random);
}


// ============================================================================================


cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation_multipoles::CovarianceMatrix_TwoPointCorrelation_multipoles (const catalogue::Catalogue random,
  const BinType binType,
  const double rMin,
  const double rMax,
  const double binSize,
  const glob::FuncGrid interpXi_0,
  const glob::FuncGrid interpXi_2,
  const glob::FuncGrid interpXi_4,
  const double minSeparation,
  const double maxSeparation)
{
  m_binType = binType;
  m_rMin = rMin;
  m_rMax = rMax;
  m_binSize = binSize;

  switch (m_binType) {

    case BinType::_linear_:
    {
      m_nBins = static_cast<int>((m_rMax-m_rMin)/m_binSize);
      break;
    }

    case BinType::_logarithmic_:
    {
      double _rMin = log10(rMin);
      double _rMax = log10(rMax);
      m_nBins = static_cast<int>((_rMax-_rMin)/m_binSize);
      break;
    }

    default:
      ErrorCBL("no such a bin type!", "constructor", "CovarianceMatrix_TwoPointCorrelation_multipoles.cpp");
  }

  m_interpXi_0 = std::make_shared<cbl::glob::FuncGrid>(interpXi_0);
  m_interpXi_2 = std::make_shared<cbl::glob::FuncGrid>(interpXi_2);
  m_interpXi_4 = std::make_shared<cbl::glob::FuncGrid>(interpXi_4);
  m_minSep = minSeparation;
  m_maxSep = maxSeparation;

  m_set_catalogues(random);
}


// ============================================================================================


void cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation_multipoles::compute_terms (const std::string output_dir,
  const bool compute_RR, const bool compute_C4, const bool compute_C3, const bool compute_C2, const bool tcount)
{
  if (compute_RR)
  {
    m_compute_RR(tcount);
    write_matrix(m_random_random, output_dir, "RR.dat");
  } else {
    std::vector<double> i1, i2;
    std::vector<std::vector<double>> matrix;
    read_matrix(output_dir+"RR.dat", i1, i2, matrix);
    m_random_random = wrapper::eigen::MatrixToEigen(matrix);
  }

  if (compute_C2)
  {
    m_compute_C2(tcount);
    write_matrix(m_C2, output_dir, "C2.dat");
  } else {
    std::vector<double> i1, i2;
    std::vector<std::vector<double>> matrix;
    read_matrix(output_dir+"C2.dat", i1, i2, matrix);
    m_C2 = wrapper::eigen::MatrixToEigen(matrix);
  }

  if (compute_C3)
  {
    m_compute_C3(tcount);
    write_matrix(m_C3, output_dir, "C3.dat");
  } else {
    std::vector<double> i1, i2;
    std::vector<std::vector<double>> matrix;
    read_matrix(output_dir+"C3.dat", i1, i2, matrix);
    m_C3 = wrapper::eigen::MatrixToEigen(matrix);
  }

  if (compute_C4)
  {
    m_compute_C4(tcount);
    write_matrix(m_C4, output_dir, "C4.dat");
  } else {
    std::vector<double> i1, i2;
    std::vector<std::vector<double>> matrix;
    read_matrix(output_dir+"C4.dat", i1, i2, matrix);
    m_C4 = wrapper::eigen::MatrixToEigen(matrix);
  }
}


// ============================================================================================


cbl::data::CovarianceMatrix cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation_multipoles::operator() (const double alpha)
{
  Eigen::MatrixXd matrix = (m_C4+alpha*m_C3+alpha*alpha*m_C2);

  return cbl::data::CovarianceMatrix(cbl::wrapper::eigen::EigenToMatrix(matrix));
}
