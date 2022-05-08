/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Cosmology/Lib/PkXiNonLinear.cpp
 *
 *  @brief Methods of the class Cosmology used to model two-point
 *  statistics
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model the two-point correlation function and
 *  power spectrum
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unibo.it
 */

#include "Cosmology.h"

using namespace std;

using namespace cbl;
using namespace cosmology;


// =====================================================================================


double cbl::cosmology::Cosmology::f_k (const double kk, const shared_ptr<cbl::glob::FuncGrid> Pk, const double qmin, const double qmax, const double prec)
{
  auto integrand = [&] (const double qq) {
		     double _Pk = Pk->operator()(qq);
		     double fact = 6.*pow(kk, 7)*qq-79.*pow(kk, 5)*pow(qq, 3)+50.*pow(kk,3)*pow(qq, 5)-21.*kk*pow(qq, 7)+3./4*pow(kk*kk-qq*qq, 3)*(2.*kk*kk+7.*qq*qq)*log(pow(fabs(kk-qq), 2)/pow(fabs(kk+qq), 2));
		     return pow(504.*pow(kk, 3)*pow(qq, 5), -1)*fact*_Pk*qq*qq;
		   };
  
  return 4.*par::pi*wrapper::gsl::GSL_integrate_cquad(integrand, qmin, qmax, prec);
}


// =====================================================================================


double cbl::cosmology::Cosmology::g_k (const double kk, const shared_ptr<cbl::glob::FuncGrid> Pk, const double qmin, const double qmax, const double prec)
{
  auto integrand = [&] (const double qq) {
		     double _Pk = Pk->operator()(qq);
		     double fact = 6.*pow(kk, 7)*qq-41.*pow(kk, 5)*pow(qq, 3)+2.*pow(kk,3)*pow(qq, 5)-3.*kk*pow(qq, 7)+3./4*pow(kk*kk-qq*qq, 3)*(2.*kk*kk+qq*qq)*log(pow(fabs(kk-qq), 2)/pow(fabs(kk+qq), 2));
		     return pow(168.*pow(kk, 3)*pow(qq, 5), -1)*fact*_Pk*qq*qq;
		   };
  
  return 4.*par::pi*wrapper::gsl::GSL_integrate_cquad(integrand, qmin, qmax, prec);
}


// =====================================================================================


double cbl::cosmology::Cosmology::F2 (const double k, const double q, const double kq)
{
  return 5./7. + kq/2. *(k/q+q/k) + 2./7.*kq*kq;
}


// =====================================================================================


double cbl::cosmology::Cosmology::G2 (const double k, const double q, const double kq)
{
  return 3./7. + kq/2. *(k/q+q/k) + 4./7.*kq*kq;
}


// ============================================================================================


double cbl::cosmology::Cosmology::Pk_1loop (const double kk, const shared_ptr<cbl::glob::FuncGrid> Pk, const int corrtype, const double qmin, const double qmax, const double prec)
{
  function<double(double, double, double)> func1, func2;
 
  if (corrtype==0) {
    func1 = [&] (const double k, const double q, const double kq) {return F2(k, q, kq);};
    func2 = [&] (const double k, const double q, const double kq) {return F2(k, q, kq);};
  }
  else if (corrtype==1) {
    func1 = [&] (const double k, const double q, const double kq) {return F2(k, q, kq);};
    func2 = [&] (const double k, const double q, const double kq) {return G2(k, q, kq);};
  }
  else if (corrtype==2) {
    func1 = [&] (const double k, const double q, const double kq) {return G2(k, q, kq);};
    func2 = [&] (const double k, const double q, const double kq) {return G2(k, q, kq);};
  }
  else 
    ErrorCBL("the input value of corrtype is not allowed!", "Pk_1loop", "PkXiNonLinear.cpp");

  auto integrand = [&] (const double qq) 
		   {
		     auto integrand_intermediate = [&] (const double xx)
						   {
						     double kq = sqrt(kk*kk+qq*qq-2.*kk*qq*xx);
						     double akq = (kk*xx-qq)/kq;
						     return func1(kq, qq, akq)*func2(kq, qq, akq)*Pk->operator()(kq);
						   };
		     double res = wrapper::gsl::GSL_integrate_cquad(integrand_intermediate, -1, 1, prec);
		     return Pk->operator()(qq)*qq*qq*res;
		   };
  
  return 4.*par::pi*wrapper::gsl::GSL_integrate_cquad(integrand, qmin, qmax, prec);
}


// ============================================================================================


double cbl::cosmology::Cosmology::Pk_DeltaDelta (const double kk, const std::shared_ptr<cbl::glob::FuncGrid> Pk, const double qmin, const double qmax, const double prec)
{
  double GG = pow(exp(f_k(kk, Pk, qmin, qmax, prec)), 2);
  return pow(2*par::pi, 3)*GG*(Pk->operator()(kk)+Pk_1loop(kk, Pk, 0, qmin, qmax, prec));
}


// ============================================================================================


double cbl::cosmology::Cosmology::Pk_DeltaTheta (const double kk, const std::shared_ptr<cbl::glob::FuncGrid> Pk, const double qmin, const double qmax, const double prec)
{
  double GG = pow(exp(f_k(kk, Pk, qmin, qmax, prec)), 2); // to be checked: possible bug!!!
  return pow(2*par::pi, 3)*GG*(Pk->operator()(kk)+Pk_1loop(kk, Pk, 1, qmin, qmax, prec));
}


// ============================================================================================


double cbl::cosmology::Cosmology::Pk_ThetaTheta (const double kk, const std::shared_ptr<cbl::glob::FuncGrid> Pk, const double qmin, const double qmax, const double prec)
{
  double GG = exp(g_k(kk, Pk, qmin, qmax, prec)); // to be checked: possible bug!!!
  return pow(2*par::pi, 3)*GG*(Pk->operator()(kk)+Pk_1loop(kk, Pk, 2, qmin, qmax, prec));
}

// ============================================================================================
      

std::vector<double> cbl::cosmology::Cosmology::Pk_DeltaDelta (const std::vector<double> kk, const double redshift, const std::string method_Pk, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec, const std::string file_par, const bool unit1)
{
  vector<double> pkLin = Pk_matter(kk, method_Pk, false, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1);

  for (size_t i=0; i<kk.size(); i++)
    pkLin[i] /= pow(2*par::pi, 3);

  auto pkLinInterp = make_shared<glob::FuncGrid>(glob::FuncGrid(kk, pkLin, "Spline"));

  vector<double> pkDeltaDelta(kk.size(), 0);

  for (size_t i=0; i<kk.size(); i++)
    pkDeltaDelta[i] = Pk_DeltaDelta(kk[i], pkLinInterp, k_min, k_max, prec);

  return pkDeltaDelta;
}
      

// ============================================================================================


std::vector<double> cbl::cosmology::Cosmology::Pk_DeltaTheta (const std::vector<double> kk, const double redshift, const std::string method_Pk, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec, const std::string file_par, const bool unit1)
{
  vector<double> pkLin = Pk_matter(kk, method_Pk, false, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1);

  for (size_t i=0; i<kk.size(); i++)
    pkLin[i] /= pow(2*par::pi, 3);

  auto pkLinInterp = make_shared<glob::FuncGrid>(glob::FuncGrid(kk, pkLin, "Spline"));

  vector<double> pkDeltaTheta(kk.size(), 0);

  for (size_t i=0; i<kk.size(); i++)
    pkDeltaTheta[i] = Pk_DeltaTheta(kk[i], pkLinInterp, k_min, k_max, prec);

  return pkDeltaTheta;
}

// ============================================================================================


std::vector<double> cbl::cosmology::Cosmology::Pk_ThetaTheta (const std::vector<double> kk, const double redshift, const std::string method_Pk, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec, const std::string file_par, const bool unit1)
{
  vector<double> pkLin = Pk_matter(kk, method_Pk, false, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1);

  for (size_t i=0; i<kk.size(); i++)
    pkLin[i] /= pow(2*par::pi, 3);

  auto pkLinInterp = make_shared<glob::FuncGrid>(glob::FuncGrid(kk, pkLin, "Spline"));

  vector<double> pkThetaTheta(kk.size(), 0);

  for (size_t i=0; i<kk.size(); i++)
    pkThetaTheta[i] = Pk_ThetaTheta(kk[i], pkLinInterp, k_min, k_max, prec);

  return pkThetaTheta;
}


// =====================================================================================


std::vector<std::vector<double>> cbl::cosmology::Cosmology::Pk_TNS_AB_multipoles (std::vector<double> kk, const std::string method, const double redshift, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec)
{
  cbl::Path path;
  string dir = path.DirCosmo()+"/External/CPT_Library/";
  string output_tmpCPT = dir+"tmpCPT/";
  string MKout = "mkdir -p " + output_tmpCPT; if (system(MKout.c_str())) {}
  double sigma8_z0 = sigma8_Pk(method, 0., store_output, output_root);

  // Pklin_z0
  const vector<double> Pklin = Pk_matter(kk, method, false, 0., store_output, output_root, norm, k_min, k_max, prec);
  string file = "Pklin.dat";
  ofstream File_Pklin(output_tmpCPT + file);
  for (size_t nn=0; nn<kk.size(); ++nn)
    File_Pklin << kk[nn] << "\t" << Pklin[nn] << endl;
  File_Pklin.close();

  // setting parameters for PkA and PkB
  string File_par = "params_AB_termsTNS.ini";
  ofstream fsAB(output_tmpCPT + File_par);
  fsAB << redshift << "\n"
       << "1 100. \n"
       << "Pklin.dat \n"
       << "1 \n" << conv(m_hh, par::fDP6) <<"\n"
       << "2 \n" << par::TCMB << "\n"
       << "3 \n" << n_spec() << "\n"
       << "4 \n" << sigma8_z0 << "\n"
       << "5 \n" << conv(m_Omega_matter, par::fDP6) << "\n"
       << "6 \n" <<  conv(m_Omega_baryon, par::fDP6) << "\n"
       << "7 \n" << conv(m_w0, par::fDP6) << "\n"
       << "0 \n" << "1 \n" << "0. \n";
  fsAB.close();
  string calc_pk_correction = "cd " + output_tmpCPT + " && " + dir + "calc_pk_correction < " + File_par; if (system (calc_pk_correction.c_str())) {}
  string calc_pk_correction2 = "cd " + output_tmpCPT + " && " + dir + "calc_pk_correction2 < " + File_par; if (system (calc_pk_correction2.c_str())) {}

  double KK, PK, PK0EH, PK0corr, PK2EH, PK2corr, PK4EH, PK4corr;
  vector<double> kA, pkA0, pkA2, pkA4, kB, pkB0, pkB2, pkB4;

  const string filenameB = output_tmpCPT + "corr_pkred.dat"; // PkB terms
  ifstream finB(filenameB.c_str());
  while (finB >> KK >> PK >> PK0EH >> PK0corr >> PK2EH >> PK2corr >> PK4EH >> PK4corr){
    kB.emplace_back(KK);
    pkB0.emplace_back(PK0corr);
    pkB2.emplace_back(PK2corr);
    pkB4.emplace_back(PK4corr);
  }
  finB.clear();

  const string filenameA = output_tmpCPT + "corr_pkred2.dat"; // PkA terms
  ifstream finA(filenameA.c_str());
  while (finA >> KK >> PK >> PK0EH >> PK0corr >> PK2EH >> PK2corr >> PK4EH >> PK4corr){
    kA.emplace_back(KK);
    pkA0.emplace_back(PK0corr);
    pkA2.emplace_back(PK2corr);
    pkA4.emplace_back(PK4corr);
  }
  finA.clear();

  vector<double> pkA0_new(kk.size()), pkA2_new(kk.size()), pkA4_new(kk.size()), pkB0_new(kk.size()), pkB2_new(kk.size()), pkB4_new(kk.size());
  glob::FuncGrid interp_PkA0(kA, pkA0, "Spline");
  glob::FuncGrid interp_PkA2(kA, pkA2, "Spline");
  glob::FuncGrid interp_PkA4(kA, pkA4, "Spline");
  glob::FuncGrid interp_PkB0(kB, pkB0, "Spline");
  glob::FuncGrid interp_PkB2(kB, pkB2, "Spline");
  glob::FuncGrid interp_PkB4(kB, pkB4, "Spline");

  pkA0_new = interp_PkA0.eval_func(kk);
  pkA2_new = interp_PkA2.eval_func(kk);
  pkA4_new = interp_PkA4.eval_func(kk);
  pkB0_new = interp_PkB0.eval_func(kk);
  pkB2_new = interp_PkB2.eval_func(kk);
  pkB4_new = interp_PkB4.eval_func(kk);

  if (system (("rm -rf "+output_tmpCPT).c_str())) {}

  return {pkA0_new, pkA2_new, pkA4_new, pkB0_new, pkB2_new, pkB4_new};
}


// =====================================================================================


std::vector<std::vector<double>> cbl::cosmology::Cosmology::Pk_TNS_AB_1loop (std::vector<double> kk, const double mu, const std::string method, const double redshift, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec)
{
  vector<vector<double>> Pk_AB_multipoles = Pk_TNS_AB_multipoles(kk, method, redshift, store_output, output_root, norm, k_min, k_max, prec);
  vector<double> Pk_A, Pk_B;

  for (size_t ii=0; ii < Pk_AB_multipoles[0].size(); ++ii) {
    Pk_A.emplace_back(Pk_AB_multipoles[0][ii] + Pk_AB_multipoles[1][ii]*cbl::legendre_polynomial(mu, 2) + Pk_AB_multipoles[2][ii]*cbl::legendre_polynomial(mu, 4));
    Pk_B.emplace_back(Pk_AB_multipoles[3][ii] + Pk_AB_multipoles[4][ii]*cbl::legendre_polynomial(mu, 2) + Pk_AB_multipoles[5][ii]*cbl::legendre_polynomial(mu, 4));
  }
  
  return {Pk_A, Pk_B};
}


// =====================================================================================


std::vector<std::vector<double>> cbl::cosmology::Cosmology::Pk_TNS_AB_terms_1loop (std::vector<double> kk, const std::string method, const double redshift, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec)
{
  cbl::Path path;
  string dir = path.DirCosmo()+"/External/CPT_Library/";
  string output_tmpCPT = dir+"tmpCPT/";
  string MKout = "mkdir -p " + output_tmpCPT; if (system(MKout.c_str())) {}
  double sigma8_z0 = sigma8_Pk(method, 0., store_output, output_root);

  // input Pklin_z0
  const vector<double> Pklin = Pk_matter(kk, method, false, 0., store_output, output_root, norm, k_min, k_max, prec);
  string file = "Pklin.dat";
  ofstream File_Pklin(output_tmpCPT + file);
  for (size_t nn=0; nn<kk.size(); ++nn)
    File_Pklin << kk[nn] << "\t" << Pklin[nn] << endl;
  File_Pklin.close();

  // setting parameters for PkA and PkB
  string File_par = "params_AB_termsTNS.ini";
  ofstream fsAB(output_tmpCPT + File_par);
  fsAB << redshift << "\n"
       << "1 100. \n"
       << "Pklin.dat \n"
       << "1 \n" << conv(m_hh, par::fDP6) <<"\n"
       << "2 \n" << par::TCMB << "\n"
       << "3 \n" << n_spec() << "\n"
       << "4 \n" << sigma8_z0 << "\n"
       << "5 \n" << conv(m_Omega_matter, par::fDP6) << "\n"
       << "6 \n" <<  conv(m_Omega_baryon, par::fDP6) << "\n"
       << "7 \n" << conv(m_w0, par::fDP6) << "\n"
       << "0 \n" << "1 \n" << "0. \n";
  fsAB.close();
  string calc_pk_correction  = "cd " + output_tmpCPT + " && " + dir + "calc_pk_correction < "  + File_par; if (system (calc_pk_correction.c_str())) {}
  string calc_pk_correction2 = "cd " + output_tmpCPT + " && " + dir + "calc_pk_correction2 < " + File_par; if (system (calc_pk_correction2.c_str())) {}

  double KK, A11, A12, A22, A23, A33, B12, B13, B14, B22, B23, B24, B33, B34, B44;
  vector<double> kA, kB, pk_A11, pk_A12, pk_A22, pk_A23, pk_A33, pk_B12, pk_B13, pk_B14, pk_B22, pk_B23, pk_B24, pk_B33, pk_B34, pk_B44;

  const string filenameB = output_tmpCPT + "pkstd_corr_tree.dat"; // pk_B12, pk_B13, pk_B14, pk_B22, pk_B23, pk_B24, pk_B33, pk_B34, pk_B44

  ifstream finB(filenameB.c_str());
  while (finB >> KK >> B12 >> B13 >> B14 >> B22 >> B23 >> B24 >> B33 >> B34 >> B44){
    kB.emplace_back(KK);
    pk_B12.emplace_back(B12);
    pk_B13.emplace_back(B13);
    pk_B14.emplace_back(B14);
    pk_B22.emplace_back(B22);
    pk_B23.emplace_back(B23);
    pk_B24.emplace_back(B24);
    pk_B33.emplace_back(B33);
    pk_B34.emplace_back(B34);
    pk_B44.emplace_back(B44);
  }
  finB.clear();

  const string filenameA = output_tmpCPT + "pkstd_corr2_tree.dat"; // pk_A11, pk_A12, pk_A22, pk_A23, pk_A33
  ifstream finA(filenameA.c_str());
  while (finA >> KK >> A11 >> A12 >> A22 >> A23 >> A33){
    kA.emplace_back(KK);
    pk_A11.emplace_back(A11);
    pk_A12.emplace_back(A12);
    pk_A22.emplace_back(A22);
    pk_A23.emplace_back(A23);
    pk_A33.emplace_back(A33);
  }
  finA.clear();

  vector<double> pk_A11_new(kk.size()), pk_A12_new(kk.size()), pk_A22_new(kk.size()), pk_A23_new(kk.size()), pk_A33_new(kk.size()), pk_B12_new(kk.size()), pk_B13_new(kk.size()), pk_B14_new(kk.size()), pk_B22_new(kk.size()), pk_B23_new(kk.size()), pk_B24_new(kk.size()), pk_B33_new(kk.size()), pk_B34_new(kk.size()), pk_B44_new(kk.size());

  glob::FuncGrid interp_A11(kA, pk_A11, "Spline");
  glob::FuncGrid interp_A12(kA, pk_A12, "Spline");
  glob::FuncGrid interp_A22(kA, pk_A22, "Spline");
  glob::FuncGrid interp_A23(kA, pk_A23, "Spline");
  glob::FuncGrid interp_A33(kA, pk_A33, "Spline");

  glob::FuncGrid interp_B12(kB, pk_B12, "Spline");
  glob::FuncGrid interp_B13(kB, pk_B13, "Spline");
  glob::FuncGrid interp_B14(kB, pk_B14, "Spline");
  glob::FuncGrid interp_B22(kB, pk_B22, "Spline");
  glob::FuncGrid interp_B23(kB, pk_B23, "Spline");
  glob::FuncGrid interp_B24(kB, pk_B24, "Spline");
  glob::FuncGrid interp_B33(kB, pk_B33, "Spline");
  glob::FuncGrid interp_B34(kB, pk_B34, "Spline");
  glob::FuncGrid interp_B44(kB, pk_B44, "Spline");

  pk_A11_new = interp_A11.eval_func(kk);
  pk_A12_new = interp_A12.eval_func(kk);
  pk_A22_new = interp_A22.eval_func(kk);
  pk_A23_new = interp_A23.eval_func(kk);
  pk_A33_new = interp_A33.eval_func(kk);
  pk_B12_new = interp_B12.eval_func(kk);
  pk_B13_new = interp_B13.eval_func(kk);
  pk_B14_new = interp_B14.eval_func(kk);
  pk_B22_new = interp_B22.eval_func(kk);
  pk_B23_new = interp_B23.eval_func(kk);
  pk_B24_new = interp_B24.eval_func(kk);
  pk_B33_new = interp_B33.eval_func(kk);
  pk_B34_new = interp_B34.eval_func(kk);
  pk_B44_new = interp_B44.eval_func(kk);

  // define the normalization
  int Norm = norm;
  if (Norm==-1) Norm = (m_sigma8>0) ? 1 : 0;
  if (Norm==1) {
    double sigma8;
    const double RR = 8.;
    glob::FuncGrid interpPk(kk, Pklin, "Spline");
    auto func_sigma = [&] (double _k) { return pow(TopHat_WF(_k*RR)*_k, 2)*interpPk(_k); };
    sigma8 = sqrt(1./(2.*pow(par::pi, 2))*wrapper::gsl::GSL_integrate_qag (func_sigma, k_min, k_max, 1.e-5))/DN(redshift, 0.);
    m_Pk0_CAMB = pow(m_sigma8/sigma8,2);
  }
  else { m_Pk0_CAMB = 1.;}

  for (size_t i=0; i<kk.size(); i++) {
    pk_A11_new[i] *= m_Pk0_CAMB;
    pk_A12_new[i] *= m_Pk0_CAMB;
    pk_A22_new[i] *= m_Pk0_CAMB;
    pk_A23_new[i] *= m_Pk0_CAMB;
    pk_A33_new[i] *= m_Pk0_CAMB;
    pk_B12_new[i] *= m_Pk0_CAMB;
    pk_B13_new[i] *= m_Pk0_CAMB;
    pk_B14_new[i] *= m_Pk0_CAMB;
    pk_B22_new[i] *= m_Pk0_CAMB;
    pk_B23_new[i] *= m_Pk0_CAMB;
    pk_B24_new[i] *= m_Pk0_CAMB;
    pk_B33_new[i] *= m_Pk0_CAMB;
    pk_B34_new[i] *= m_Pk0_CAMB;
    pk_B44_new[i] *= m_Pk0_CAMB;
  }

  if (system (("rm -rf " + output_tmpCPT).c_str())) {}

  return {pk_A11_new, pk_A12_new, pk_A22_new, pk_A23_new, pk_A33_new, pk_B12_new, pk_B13_new, pk_B14_new, pk_B22_new, pk_B23_new, pk_B24_new, pk_B33_new, pk_B34_new, pk_B44_new};
}


// =====================================================================================


std::vector<std::vector<double>> cbl::cosmology::Cosmology::Pk_TNS_AB_1loop (std::vector<double> kk, const double mu, const double linear_growth_rate, const double bias, const std::string method, const double redshift, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec)
{
  double beta = linear_growth_rate/bias;

  vector<vector<double>> Pk_AB = Pk_TNS_AB_terms_1loop(kk, method, redshift, store_output, output_root, norm, k_min, k_max, prec);
  vector<double> Pk_A, Pk_B;

  for(size_t ii=0; ii < kk.size(); ++ii) {
    Pk_A.emplace_back(beta*pow(mu,2)*Pk_AB[0][ii] + pow(beta,2)*(pow(mu,2)*Pk_AB[1][ii] + pow(mu,4)*Pk_AB[2][ii]) + pow(beta,3)*(pow(mu,4)*Pk_AB[3][ii] + pow(mu,6)*Pk_AB[4][ii]));
    Pk_B.emplace_back(pow(mu,2)*(pow(beta,2)*Pk_AB[5][ii] + pow(beta,3)*Pk_AB[6][ii] + pow(beta,4)*Pk_AB[7][ii]) + pow(mu,4)*(pow(beta,2)*Pk_AB[8][ii] + pow(beta,3)*Pk_AB[9][ii] + pow(beta,4)*Pk_AB[10][ii]) + pow(mu,6)*(pow(beta,3)*Pk_AB[11][ii] + pow(beta,4)*Pk_AB[12][ii]) + pow(mu,8)*pow(beta,4)*Pk_AB[13][ii]);
  }

  return {Pk_A, Pk_B};
}


// =====================================================================================


std::vector<std::vector<double>> cbl::cosmology::Cosmology::Pk_TNS_dd_dt_tt (std::vector<double> kk, const std::string method, const double redshift, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec)
{
  cbl::Path path;
  string dir = path.DirCosmo()+"/External/CPT_Library/";
  string output_tmpCPT = dir+"tmpCPT/";
  string MKout = "mkdir -p " + output_tmpCPT; if (system(MKout.c_str())) {}
  double sigma8_z0 = sigma8_Pk(method, 0., store_output, output_root);

  // input Pklin_z0
  const vector<double> Pklin = Pk_matter(kk, method, false, 0., store_output, output_root, norm, k_min, k_max, prec);
  string file = "Pklin.dat";
  ofstream File_Pklin(output_tmpCPT + file);
  for (size_t nn=0; nn<kk.size(); ++nn)
    File_Pklin << kk[nn] << "\t" << Pklin[nn] << endl;
  File_Pklin.close();

  // setting input parameters
  string File_par = "params_stdPT2.ini";
  ofstream File_stdPT2(output_tmpCPT + File_par);
  File_stdPT2 << sigma8_z0 << "\n Pklin.dat";
  File_stdPT2.close();
  string stdPT2 = "cd " + output_tmpCPT + " && " + dir + "stdPT2 < " + File_par; if (system (stdPT2.c_str())) {}

  File_par = "params_read_pk2.ini";
  ofstream File_read(output_tmpCPT + File_par);
  File_read << "0 \n"
	    << "1 \n" << conv(m_Omega_matter, par::fDP6) <<"\n"
	    << "2 \n" << conv(m_Omega_baryon, par::fDP6) << "\n"
	    << "3 \n" << conv(m_hh, par::fDP6) << "\n"
	    << "4 \n" << conv(m_n_spec, par::fDP6) << "\n"
	    << "5 \n" << sigma8_z0 << "\n"
	    << "6 \n" << conv(m_w0, par::fDP6) << "\n"
	    << "0 \n" << "0 " << redshift << "\n"
	    << "1 \n" << "0 \n" << "0";
  File_read.close();
  
  string read_pk2 = "cd " + output_tmpCPT + " && " + dir + "read_pk2 < " + File_par; if (system (read_pk2.c_str())) {}

  double Kspt, PKlinNWspt, PKlinspt, PDDspt, PDTspt, PTTspt;
  vector<double> k_spt, PklinNW_spt, Pklin_spt, Pdd_spt, Pdt_spt, Ptt_spt;

  const string filenamePK = output_tmpCPT + "pkstd2.dat"; // Pkdd, Pkdt, Pktt
  ifstream finPK(filenamePK.c_str());

  while (finPK >> Kspt >> PKlinNWspt >> PKlinspt >> PDDspt >> PDTspt >> PTTspt)
    if (Kspt>0 && PKlinNWspt>0 && PKlinspt>0 && PDDspt>0 && PDTspt>0 && PTTspt>0) {
      k_spt.emplace_back(Kspt);
      PklinNW_spt.emplace_back(PKlinNWspt);
      Pklin_spt.emplace_back(PKlinspt);
      Pdd_spt.emplace_back(PDDspt);
      Pdt_spt.emplace_back(PDTspt);
      Ptt_spt.emplace_back(PTTspt);
    }
  finPK.clear();

  vector<double> Pkdd_new(kk.size()), Pkdt_new(kk.size()), Pktt_new(kk.size()), Pklin_new(kk.size());
  glob::FuncGrid interp_Pklin(k_spt, Pklin_spt, "Spline");
  glob::FuncGrid interp_Pkdd(k_spt, Pdd_spt, "Spline");
  glob::FuncGrid interp_Pkdt(k_spt, Pdt_spt, "Spline");
  glob::FuncGrid interp_Pktt(k_spt, Ptt_spt, "Spline");

  Pklin_new  = interp_Pklin.eval_func(kk);
  Pkdd_new = interp_Pkdd.eval_func(kk);
  Pkdt_new = interp_Pkdt.eval_func(kk);
  Pktt_new = interp_Pktt.eval_func(kk);

  // define the normalization
  int Norm = norm;
  if (Norm==-1) Norm = (m_sigma8>0) ? 1 : 0;
  if (Norm==1) {
    double sigma8;
    const double RR = 8.;
    glob::FuncGrid interpPk(kk, Pklin_new, "Spline");
    auto func_sigma = [&] (double _k) { return pow(TopHat_WF(_k*RR)*_k, 2)*interpPk(_k); };
    sigma8 = sqrt(1./(2.*pow(par::pi, 2))*wrapper::gsl::GSL_integrate_qag (func_sigma, k_min, k_max, 1.e-5))/DN(redshift, 0.);
    m_Pk0_CAMB = pow(m_sigma8/sigma8, 2);
  }
  else { m_Pk0_CAMB = 1.; }

  for (size_t i=0; i<kk.size(); i++) {
    Pkdd_new[i] *= m_Pk0_CAMB;
    Pkdt_new[i] *= m_Pk0_CAMB;
    Pktt_new[i] *= m_Pk0_CAMB;
  }

  if (system (("rm -rf "+output_tmpCPT).c_str())) {}

  return {Pkdd_new, Pkdt_new, Pktt_new};
}


// =====================================================================================


std::vector<std::vector<double>> cbl::cosmology::Cosmology::Pk_eTNS_terms_1loop (std::vector<double> kk, const std::string method, const double redshift, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec)
{
  cbl::Path path;
  string dir = path.DirCosmo()+"/External/CAMB_SPT_private/";
  string output_tmpCPT = dir+"tmpCPT_eTNS/";
  string MKout = "mkdir -p " + output_tmpCPT; if (system(MKout.c_str())) {}
  double HH0 = m_hh*100.;

  // input Pklin_z0
  const vector<double> Pklin = Pk_matter(kk, method, false, 0., store_output, output_root, norm, k_min, k_max, prec);
  string file = "Pklin.dat";
  ofstream File_Pklin(output_tmpCPT + file);
  for (size_t nn=0; nn<kk.size(); ++nn)
    File_Pklin << kk[nn] << "\t" << Pklin[nn] << endl;
  File_Pklin.close();

  // setting parameters for Pk_xy
  string File_par = "SPT_NLB_params.ini";
  ofstream fsPkt(output_tmpCPT + File_par);
  fsPkt << "output_root = SPT_NLB \n" // OutputRoot
	<< "get_scalar_cls = F \n"
    	<< "get_vector_cls = F \n"
    	<< "get_tensor_cls = F \n"
    	<< "COBE_normalize = F \n"
    	<< "CMB_outputscale = 7.4311e12\n"
    	<< "get_transfer = T \n"
    	<< "do_nonlinear = 3 \n"
    	<< "DoNonLocalBias = T \n"
    	<< "RSDmodel = 1 \n"
    	<< "w = " << conv(m_w0, par::fDP6) <<"\n"
    	<< "cs2_lam = 1 \n"
    	<< "hubble = " << conv(HH0, par::fDP6) <<"\n"
    	<< "use_physical = T \n"
    	<< "ombh2 = " << conv(m_Omega_baryon*m_hh*m_hh, par::fDP6) <<"\n"
    	<< "omch2 = " << conv(m_Omega_CDM*m_hh*m_hh, par::fDP6) <<"\n"
    	<< "omnuh2 = " << conv(m_Omega_neutrinos*m_hh*m_hh, par::fDP6) <<"\n"
    	<< "omk = " << conv(m_Omega_k, par::fDP6) <<"\n"
    	<< "temp_cmb = " << cbl::par::TCMB <<"\n"
    	<< "helium_fraction = 0.24 \n"
    	<< "massless_neutrinos = " << conv(m_massless_neutrinos, par::fDP6) <<"\n"
    	<< "massive_neutrinos = " << conv(m_massive_neutrinos, par::fINT) <<"\n"
    	<< "nu_mass_eigenstates = 1 \n"
    	<< "nu_mass_degeneracies = 0 \n"
    	<< "nu_mass_fractions = 1 \n"
    	<< "transfer_high_precision = T \n"
    	<< "transfer_kmax = " << k_max <<"\n"
    	<< "transfer_k_per_logint = 20 \n"
    	<< "transfer_num_redshifts = 1 \n"
    	<< "transfer_interp_matterpower = T \n"
    	<< "transfer_power_var = 7 \n"
    	<< "transfer_redshift(1) = " << cbl::conv(redshift, cbl::par::fDP1) <<"\n"
    	<< "transfer_filename(1) = Tk_z.dat \n"
    	<< "transfer_matterpower(1) = Pk_z.dat \n"
    	<< "reionization = T \n"
    	<< "re_use_optical_depth = T \n"
    	<< "re_optical_depth = " << conv(m_tau, par::fDP6) <<"\n"
    	<< "re_delta_redshift = 1.5 \n"
    	<< "re_ionization_frac = -1 \n"
    	<< "pivot_scalar = " << conv(m_scalar_pivot, par::fDP6) <<"\n"
    	<< "pivot_tensor = 0.002 \n"
    	<< "initial_power_num = 1 \n"
    	<< "scalar_spectral_index(1) = " << conv(m_n_spec, par::fDP6) <<"\n"
    	<< "scalar_nrun(1) = 0 \n"
    	<< "scalar_amp(1) = " << conv(m_scalar_amp, par::ee3) <<"\n"
    	<< "RECFAST_fudge_He = 0.86 \n"
    	<< "RECFAST_Heswitch = 6 \n"
    	<< "RECFAST_Hswitch = T \n"
    	<< "RECFAST_fudge = 1.14 \n"
    	<< "do_lensing_bispectrum = F \n"
    	<< "do_primordial_bispectrum = F \n"
    	<< "initial_condition = 1 \n"
    	<< "accurate_polarization = T \n"
    	<< "accurate_reionization = T \n"
    	<< "accurate_BB = F \n"
    	<< "do_late_rad_truncation = T \n"
    	<< "do_tensor_neutrinos = T \n"
    	<< "feedback_level = 1 \n"
    	<< "massive_nu_approx = 1 \n"
    	<< "number_of_threads = 0 \n"
    	<< "accuracy_boost = 2 \n"
    	<< "l_accuracy_boost = 1 \n"
    	<< "high_accuracy_default = F \n"
    	<< "l_sample_boost = 1 ";
  fsPkt.close();
  string calc_pk_eTNScorrection  = "cd " + output_tmpCPT + " && " + dir + "camb "  + File_par; if (system (calc_pk_eTNScorrection.c_str())) {}

  double Kspt, PKlin, PDD, PDV, PVV, PB2D, PB2V, PB22, PBS2D, PBS2V, PB2S2, PBS22, sigma32PKlin, BB1, BB2, BBS2;
  vector<double> k_spt, Pdd, Pdv, Pvv, Pb2d, Pb2v, Pb22, Pbs2d, Pbs2v, Pb2s2, Pbs22, sigma32Pklin, Bb1, Bb2, Bbs2;

  const string filenamePK = output_tmpCPT + "SPT_NLB_Pk_z.dat"; // k, Pklin, Pdd, Pdv, Pvv, Pb2d, Pb2v, Pb22, Pbs2d, Pbs2v, Pb2s2, Pbs22, sigma3^2Pklin, Bb1, Bb2, Bbs2

  ifstream finPK(filenamePK.c_str());
  while (finPK >> Kspt >> PKlin >> PDD >> PDV >> PVV >> PB2D >> PB2V >> PB22 >> PBS2D >> PBS2V >> PB2S2 >> PBS22 >> sigma32PKlin >> BB1 >> BB2 >> BBS2)
    if (Kspt>0 && PKlin>0 && PDD>0 && PDV>0 && PVV>0) {
      k_spt.emplace_back(Kspt);
      Pdd.emplace_back(PDD);
      Pdv.emplace_back(PDV);
      Pvv.emplace_back(PVV);
      Pb2d.emplace_back(PB2D);
      Pb2v.emplace_back(PB2V);
      Pb22.emplace_back(PB22);
      Pbs2d.emplace_back(PBS2D);
      Pbs2v.emplace_back(PBS2V);
      Pb2s2.emplace_back(PB2S2);
      Pbs22.emplace_back(PBS22);
      sigma32Pklin.emplace_back(sigma32PKlin);
      Bb1.emplace_back(BB1);
      Bb2.emplace_back(BB2);
      Bbs2.emplace_back(BBS2);
    }
  finPK.clear();

  vector<double> Pdd_new(kk.size()), Pdv_new(kk.size()), Pvv_new(kk.size()), Pb2d_new(kk.size()), Pb2v_new(kk.size()), Pb22_new(kk.size()), Pbs2d_new(kk.size()), Pbs2v_new(kk.size()), Pb2s2_new(kk.size()), Pbs22_new(kk.size()), sigma32Pklin_new(kk.size()), Bb1_new(kk.size()), Bb2_new(kk.size()), Bbs2_new(kk.size());

  glob::FuncGrid interp_Pdd(k_spt, Pdd, "Spline");
  glob::FuncGrid interp_Pdv(k_spt, Pdv, "Spline");
  glob::FuncGrid interp_Pvv(k_spt, Pvv, "Spline");
  glob::FuncGrid interp_Pb2d(k_spt, Pb2d, "Spline");
  glob::FuncGrid interp_Pb2v(k_spt, Pb2v, "Spline");
  glob::FuncGrid interp_Pb22(k_spt, Pb22, "Spline");
  glob::FuncGrid interp_Pbs2d(k_spt, Pbs2d, "Spline");
  glob::FuncGrid interp_Pbs2v(k_spt, Pbs2v, "Spline");
  glob::FuncGrid interp_Pb2s2(k_spt, Pb2s2, "Spline");
  glob::FuncGrid interp_Pbs22(k_spt, Pbs22, "Spline");
  glob::FuncGrid interp_sigma32Pklin(k_spt, sigma32Pklin, "Spline");
  glob::FuncGrid interp_Bb1(k_spt, Bb1, "Spline");
  glob::FuncGrid interp_Bb2(k_spt, Bb2, "Spline");
  glob::FuncGrid interp_Bbs2(k_spt, Bbs2, "Spline");

  Pdd_new = interp_Pdd.eval_func(kk);
  Pdv_new = interp_Pdv.eval_func(kk);
  Pvv_new = interp_Pvv.eval_func(kk);
  Pb2d_new = interp_Pb2d.eval_func(kk);
  Pb2v_new = interp_Pb2v.eval_func(kk);
  Pb22_new = interp_Pb22.eval_func(kk);
  Pbs2d_new = interp_Pbs2d.eval_func(kk);
  Pbs2v_new = interp_Pbs2v.eval_func(kk);
  Pb2s2_new = interp_Pb2s2.eval_func(kk);
  Pbs22_new = interp_Pbs22.eval_func(kk);
  sigma32Pklin_new = interp_sigma32Pklin.eval_func(kk);
  Bb1_new = interp_Bb1.eval_func(kk);
  Bb2_new = interp_Bb2.eval_func(kk);
  Bbs2_new = interp_Bbs2.eval_func(kk);

  if (system (("rm -rf "+output_tmpCPT).c_str())) {}

  return {Pdd_new, Pdv_new, Pvv_new, Pb2d_new, Pb2v_new, Pb22_new, Pbs2d_new, Pbs2v_new, Pb2s2_new, Pbs22_new, sigma32Pklin_new, Bb1_new, Bb2_new, Bbs2_new};
}
