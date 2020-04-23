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
 
  if(corrtype==0) {
    func1 = [&] (const double k, const double q, const double kq) {return F2(k, q, kq);};
    func2 = [&] (const double k, const double q, const double kq) {return F2(k, q, kq);};
  }
  else if(corrtype==1) {
    func1 = [&] (const double k, const double q, const double kq) {return F2(k, q, kq);};
    func2 = [&] (const double k, const double q, const double kq) {return G2(k, q, kq);};
  }
  else if(corrtype==2) {
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
  double GG = pow(exp(f_k(kk, Pk, qmin, qmax, prec)), 2);
  return pow(2*par::pi, 3)*GG*(Pk->operator()(kk)+Pk_1loop(kk, Pk, 1, qmin, qmax, prec));
}


// ============================================================================================


double cbl::cosmology::Cosmology::Pk_ThetaTheta (const double kk, const std::shared_ptr<cbl::glob::FuncGrid> Pk, const double qmin, const double qmax, const double prec)
{
  double GG = exp(g_k(kk, Pk, qmin, qmax, prec));
  return pow(2*par::pi, 3)*GG*(Pk->operator()(kk)+Pk_1loop(kk, Pk, 2, qmin, qmax, prec));
}

// ============================================================================================
      

std::vector<double> cbl::cosmology::Cosmology::Pk_DeltaDelta (const std::vector<double> kk, const double redshift, const std::string method_Pk, const std::string output_dir, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec, const std::string file_par, const bool unit1)
{
  vector<double> pkLin = Pk_DM(kk, method_Pk, false, redshift, output_dir, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1);

  for (size_t i=0; i<kk.size(); i++)
    pkLin[i] /= pow(2*par::pi, 3);

  auto pkLinInterp = make_shared<glob::FuncGrid>(glob::FuncGrid(kk, pkLin, "Spline"));

  vector<double> pkDeltaDelta(kk.size(), 0);

  for (size_t i=0; i<kk.size(); i++)
    pkDeltaDelta[i] = Pk_DeltaDelta(kk[i], pkLinInterp, k_min, k_max, prec);

  return pkDeltaDelta;
}
      

// ============================================================================================


std::vector<double> cbl::cosmology::Cosmology::Pk_DeltaTheta (const std::vector<double> kk, const double redshift, const std::string method_Pk, const std::string output_dir, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec, const std::string file_par, const bool unit1)
{
  vector<double> pkLin = Pk_DM(kk, method_Pk, false, redshift, output_dir, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1);

  for (size_t i=0; i<kk.size(); i++)
    pkLin[i] /= pow(2*par::pi, 3);

  auto pkLinInterp = make_shared<glob::FuncGrid>(glob::FuncGrid(kk, pkLin, "Spline"));

  vector<double> pkDeltaTheta(kk.size(), 0);

  for (size_t i=0; i<kk.size(); i++)
    pkDeltaTheta[i] = Pk_DeltaTheta(kk[i], pkLinInterp, k_min, k_max, prec);

  return pkDeltaTheta;
}

// ============================================================================================


std::vector<double> cbl::cosmology::Cosmology::Pk_ThetaTheta (const std::vector<double> kk, const double redshift, const std::string method_Pk, const std::string output_dir, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec, const std::string file_par, const bool unit1)
{
  vector<double> pkLin = Pk_DM(kk, method_Pk, false, redshift, output_dir, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1);

  for (size_t i=0; i<kk.size(); i++)
    pkLin[i] /= pow(2*par::pi, 3);

  auto pkLinInterp = make_shared<glob::FuncGrid>(glob::FuncGrid(kk, pkLin, "Spline"));

  vector<double> pkThetaTheta(kk.size(), 0);

  for (size_t i=0; i<kk.size(); i++)
    pkThetaTheta[i] = Pk_ThetaTheta(kk[i], pkLinInterp, k_min, k_max, prec);

  return pkThetaTheta;
}
