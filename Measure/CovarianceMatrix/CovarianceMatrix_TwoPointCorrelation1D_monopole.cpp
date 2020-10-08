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
 *  Measure/CovarianceMatrix/CovarianceMatrix_TwoPointCorrelation1D_monopole.cpp
 *
 *  @brief Methods of the class CovarianceMatrix_TwoPointCorrelation1D_monopole used to
 *  measure the covariance of the monopole of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  CovarianceMatrix_TwoPointCorrelation1D_monopole used to measure the covariance of the
 *  monopole of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "EigenWrapper.h"
#include "CovarianceMatrix_TwoPointCorrelation1D_monopole.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace chainmesh;
using namespace measure;
using namespace covmat;


// ============================================================================================


void cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation1D_monopole::m_compute_RR(const bool tcount)
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

  std::vector<double> random_pairs(m_nBins, 0);
  std::vector<std::vector<double>> rand_rand(m_nBins, std::vector<double>(m_nBins, 0));

#pragma omp parallel num_threads(omp_get_max_threads()) private(tid)
  {
    tid = omp_get_thread_num();

    std::vector<double> random_pairs_thread(m_nBins, 0);

    // parallelized loop
#pragma omp for schedule(static, 2)
    for (int i=0; i<nObj; ++i) {

      Vector4D pos_i = m_random->eigen_coordinate(i);
      double ww_i = m_random->weight(i);

      vector<long int> close_objects = m_chainMesh_rMax.close_objects({pos_i[0], pos_i[1], pos_i[2]}, -1);

      for (auto &&j : close_objects) {
        Vector4D pos_j = cat2->eigen_coordinate(j);
        double ww_j = cat2->weight(j);

        double rr = (pos_i-pos_j).norm();

        if ((rr<m_rMin) || (rr>m_rMax) || rr==0)
          continue;

        int kk = static_cast<int>((rr-m_rMin)/m_binSize);
        random_pairs_thread[kk] += ww_i*ww_j;
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
        random_pairs[i] += random_pairs_thread[i];
    }
  }

  for (size_t i=0; i<m_nBins; ++i)
    for (size_t j=0; j<m_nBins; ++j)
      rand_rand[i][j] = random_pairs[i]*random_pairs[j];

  m_random_random = cbl::wrapper::eigen::MatrixToEigen(rand_rand);

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


void cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation1D_monopole::m_compute_C4(const bool tcount)
{
  coutCBL << "I'm computing the C4 term" << endl;

  // timer
  time_t start; time (&start);
  int dp = cout.precision();
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(2);

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

  std::vector<std::vector<double>> matrix(m_nBins, std::vector<double>(m_nBins, 0));

  #pragma omp parallel num_threads(omp_get_max_threads()) private(tid)
    {
      tid = omp_get_thread_num();

      std::vector<std::vector<double>> matrix_thread(m_nBins, std::vector<double>(m_nBins, 0));

      // parallelized loop
  #pragma omp for schedule(static, 2)
      for (int i=0; i<nObj; ++i) {

        Vector4D pos_i = m_random->eigen_coordinate(i);
        double ww_i = m_random->weight(i);

        vector<long int> close_objects_j = m_chainMesh_rMax.close_objects({pos_i[0], pos_i[1], pos_i[2]}, -1);
        size_t nEff_j = 0;
        vector<size_t> list_j (close_objects_j.size());
        vector<size_t> bin_j (close_objects_j.size());

        vector<long int> close_objects_k = m_chainMesh_rCut.close_objects({pos_i[0], pos_i[1], pos_i[2]}, -1);

        for (auto &&j : close_objects_j) {
          Vector4D pos_j = cat2->eigen_coordinate(j);
          double rr_ij = (pos_i-pos_j).norm();

          if ((rr_ij<m_rMin) || (rr_ij>m_rMax) || rr_ij==0)
            continue;

          list_j[nEff_j] = j;
          bin_j[nEff_j] = static_cast<int>((rr_ij-m_rMin)/m_binSize);
          nEff_j +=1 ;
        }

        for (auto &&k : close_objects_k) {
          Vector4D pos_k = cat3->eigen_coordinate(k);
          double ww_k = cat3->weight(k);

          double rr_ik = (pos_k-pos_i).norm();
          double xi_ik = ((rr_ik<m_minSep) || (rr_ik>m_maxSep) ) ? 0. : m_interpXi->operator()(rr_ik);

          vector<long int> close_objects_l = m_chainMesh_rMax.close_objects({pos_k[0], pos_k[1], pos_k[2]}, -1);
          for (auto &&l : close_objects_l) {
            Vector4D pos_l = cat2->eigen_coordinate(l);
            double ww_l = cat2->weight(l);

            double rr_kl = (pos_l-pos_k).norm();
            if ((rr_kl<m_rMin) || (rr_kl>m_rMax) || rr_kl==0)
              continue;

            int bin2 = static_cast<int>((rr_kl-m_rMin)/m_binSize);

            for (size_t j=0; j<nEff_j; j++) {
              double rr_jl = (pos_l - cat2->eigen_coordinate(list_j[j])).norm();
              double ww_j = cat2->weight(list_j[j]);

              double xi_jl = ((rr_jl<m_minSep) || (rr_jl>m_maxSep)) ? 0. : m_interpXi->operator()(rr_jl);

              matrix_thread[bin_j[j]][bin2] += 2. * ww_i * ww_j * ww_k * ww_l * xi_ik * xi_jl;
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
        for (size_t i=0; i<m_nBins; ++i)
          for (size_t j=0; j<m_nBins; ++j)
            matrix[i][j] += matrix_thread[i][j];
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


void cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation1D_monopole::m_compute_C3(const bool tcount)
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

  // factor used by the timer
  float fact_count = 100./nObj;

  // thread number
  int tid = 0;

  std::vector<std::vector<double>> matrix(m_nBins, std::vector<double>(m_nBins, 0));

  #pragma omp parallel num_threads(omp_get_max_threads()) private(tid)
    {
      tid = omp_get_thread_num();

      std::vector<std::vector<double>> matrix_thread(m_nBins, std::vector<double>(m_nBins, 0));

      // parallelized loop
  #pragma omp for schedule(static, 2)
      for (int i=0; i<nObj; ++i) {
        Vector4D pos_i = m_random->eigen_coordinate(i);
        double ww_i = m_random->weight(i);

        vector<long int> close_objects = m_chainMesh_rMax.close_objects({pos_i[0], pos_i[1], pos_i[2]}, -1);

        for (auto &&j : close_objects) {
          Vector4D pos_j = cat2->eigen_coordinate(j);
          double ww_j = cat2->weight(j);

          double rr_ij = (pos_i-pos_j).norm();

          if ((rr_ij<m_rMin) || (rr_ij>m_rMax) || rr_ij==0)
            continue;

          int bin1 = static_cast<int>((rr_ij-m_rMin)/m_binSize);

          vector<long int> close_objects_j = m_chainMesh_rMax.close_objects({pos_j[0], pos_j[1], pos_j[2]}, -1);

          for (auto &&k : close_objects_j) {
            Vector4D pos_k = cat2->eigen_coordinate(k);
            double ww_k = cat2->weight(k);

            double rr_jk = (pos_k-pos_j).norm();

            if ((rr_jk<m_rMin) || (rr_jk>m_rMax) || rr_jk==0)
              continue;

            int bin2 = static_cast<int>((rr_jk-m_rMin)/m_binSize);
            double rr_ik = (pos_k-pos_i).norm();

            double xi = ((rr_ik<m_minSep) || (rr_ik>m_maxSep)) ? 0. : m_interpXi->operator()(rr_ik);
            matrix_thread[bin1][bin2] += ww_i * ww_j * ww_j * ww_k * xi ;
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
        for (size_t i=0; i<m_nBins; ++i)
          for (size_t j=0; j<m_nBins; ++j)
            matrix[i][j] += matrix_thread[i][j];
      }
    }

  m_C3 = 4*cbl::wrapper::eigen::MatrixToEigen(matrix);

  m_random->Order();

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


void cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation1D_monopole::m_compute_C2(const bool tcount)
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

  std::vector<double> c2(m_nBins, 0);
  std::vector<std::vector<double>> matrix(m_nBins, std::vector<double>(m_nBins, 0));

#pragma omp parallel num_threads(omp_get_max_threads()) private(tid)
  {
    tid = omp_get_thread_num();

    std::vector<double> c2_thread(m_nBins, 0);

    // parallelized loop
#pragma omp for schedule(static, 2)
    for (int i=0; i<nObj; ++i) {
      Vector4D pos_i = m_random->eigen_coordinate(i);
      double ww_i = m_random->weight(i);

      vector<long int> close_objects = m_chainMesh_rMax.close_objects({pos_i[0], pos_i[1], pos_i[2]}, -1);

      for (auto &&j : close_objects) {
        Vector4D pos_j = cat2->eigen_coordinate(j);
        double ww_j = cat2->weight(j);

        double rr = (pos_i-pos_j).norm();

        if ((rr<m_rMin) || (rr>m_rMax) || rr==0)
          continue;

        int kk = static_cast<int>((rr-m_rMin)/m_binSize);
        c2_thread[kk] += ww_i*ww_j*ww_i*ww_j*(1+m_interpXi->operator()(rr));
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
        c2[i] += c2_thread[i];
    }
  }

  for (size_t i=0; i<m_nBins; ++i)
      matrix[i][i] = 2 * c2[i];

  m_C2 = cbl::wrapper::eigen::MatrixToEigen(matrix);

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


void cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation1D_monopole::write_matrix (Eigen::MatrixXd matrix,
  const std::string dir,
  const std::string file,
  const int nDigits)
{
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  fout << "### [1] bin1 # [2] bin2 # [3] matrix ### " << endl;

  for (size_t i=0; i<m_nBins; ++i){
    for (size_t j=0; j<m_nBins; ++j) {
      fout << setprecision(nDigits) << i << " " << j << " " << matrix(i, j) << endl;
    }
    fout << endl;
  }

  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl;
}


// ============================================================================================

void cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation1D_monopole::m_set_catalogues(const catalogue::Catalogue random)
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


cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation1D_monopole::CovarianceMatrix_TwoPointCorrelation1D_monopole (const catalogue::Catalogue random,
  const BinType binType,
  const double rMin,
  const double rMax,
  const int nBins,
  const glob::FuncGrid interpXi,
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
      ErrorCBL("no such a bin type!", "constructor", "CovarianceMatrix_TwoPointCorrelation1D_monopole.cpp");
  }


  m_interpXi = std::make_shared<cbl::glob::FuncGrid>(interpXi);
  m_minSep = minSeparation;
  m_maxSep = maxSeparation;

  m_set_catalogues(random);
}


// ============================================================================================


cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation1D_monopole::CovarianceMatrix_TwoPointCorrelation1D_monopole (const catalogue::Catalogue random,
  const BinType binType,
  const double rMin,
  const double rMax,
  const double binSize,
  const glob::FuncGrid interpXi,
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
      ErrorCBL("no such a bin type!", "constructor", "CovarianceMatrix_TwoPointCorrelation1D_monopole.cpp");
  }

  m_interpXi = std::make_shared<cbl::glob::FuncGrid>(interpXi);
  m_minSep = minSeparation;
  m_maxSep = maxSeparation;

  m_set_catalogues(random);
}


// ============================================================================================


void cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation1D_monopole::compute_terms (const std::string output_dir,
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


cbl::data::CovarianceMatrix cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation1D_monopole::operator() (const double alpha)
{
  Eigen::MatrixXd matrix = (m_C4+alpha*m_C3+alpha*alpha*m_C2);
  matrix.cwiseQuotient(m_random_random);

  return cbl::data::CovarianceMatrix(cbl::wrapper::eigen::EigenToMatrix(matrix));
}
