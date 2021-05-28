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
 *  \@file CosmoBolognaLib/Measure/ThreePointCorrelation/ThreePointCorrelation_comoving_multipoles.cpp
 *
 *  @brief Methods of the class
 *  ThreePointCorrelation_comoving_multipoles used to compute the 
 *  the multipoles of the three-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  ThreePointCorrelation_comoving_multipoles used to compute the 
 *  the multipoles of the three-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */


#include "ThreePointCorrelation_comoving_multipoles.h"

using namespace std;

using namespace cbl;

using namespace catalogue;
using namespace measure;
using namespace threept;
using namespace glob;

// ============================================================================


shared_ptr<cbl::catalogue::Catalogue> cbl::measure::threept::ThreePointCorrelation_comoving_multipoles::m_join_catalogues  (const cbl::catalogue::Catalogue& data, const cbl::catalogue::Catalogue& random, const size_t startPos, const size_t endPos) const
{
  // ---------------- check inputs ---------------------------------

  if (startPos > random.nObjects() || endPos > random.nObjects())
    ErrorCBL("input indeces are larger than random catalogue size!", "m_join_catalogues", "ThreePointCorrelation_comoving_multipoles.cpp");

  // ---------------------------------------------------------------
  // ---------------- construct the mixed catalogue ----------------
  // ---------------------------------------------------------------

  double fact = 0;
  for (size_t i=startPos; i<endPos; i++) 
    fact += random.weight(i);

  fact = -double(data.weightedN())/fact;

  catalogue::Catalogue joined_catalogue;

  for (size_t i=0; i<data.nObjects(); i++) {
    double xx = data.xx(i);
    double yy = data.yy(i);
    double zz = data.zz(i);
    double ww = data.weight(i);

    auto obj = make_shared<catalogue::Object>(catalogue::Object());
    obj->set_xx(xx);
    obj->set_yy(yy);
    obj->set_zz(zz);
    obj->set_weight(ww);
    joined_catalogue.add_object(obj);
  }

  for (size_t i=startPos; i<endPos; i++) {
    double xx = random.xx(i);
    double yy = random.yy(i);
    double zz = random.zz(i);
    double ww = random.weight(i)*fact;

    auto obj = make_shared<catalogue::Object>(catalogue::Object());
    obj->set_xx(xx);
    obj->set_yy(yy);
    obj->set_zz(zz);
    obj->set_weight(ww);
    joined_catalogue.add_object(obj);
  }

  return make_shared<Catalogue>(joined_catalogue);
}


// ============================================================================


vector<double> cbl::measure::threept::ThreePointCorrelation_comoving_multipoles::m_SzapudiSzalay_multipoles (const std::vector<double> nnn, const std::vector<double> rrr, const double normalization) const
{
  vector<double> _NNN = nnn;
  vector<double> _RRR = rrr;

  for (size_t i=0; i<m_nOrders; i++) {
    _NNN[i] *= 0.5*(2.*i+1.);
    _RRR[i] *= 0.5*(2.*i+1.);
  }

  vector<double> fl = _RRR;
  for (size_t i=1; i<m_nOrders; i++) {
    fl[i] = fl[i]/_RRR[0];
  }

  vector<vector<double>> A(m_nOrders, vector<double>(m_nOrders, 0)), A_inverse;

  for (size_t k=0; k<m_nOrders; k++)
    for (size_t l=0; l<m_nOrders; l++)
      for (size_t lp=1; lp<m_nOrders; lp++)
	A[k][l] += (2*k+1)*pow(gsl_sf_coupling_3j(2*l, 2*lp, 2*k, 0, 0, 0), 2)*fl[lp];

  for (size_t k=0; k<m_nOrders; k++)
    A[k][k] += 1.;

  invert_matrix(A, A_inverse, 1.e-4);

  vector<double> zeta_l(m_nOrders, 0.);
  for (size_t l=0; l<m_nOrders; l++)
    for (size_t k=0; k<m_nOrders; k++)
      zeta_l[l] += _NNN[k]*A_inverse[l][k]/_RRR[0]*normalization;

  return zeta_l;
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_comoving_multipoles::set_catalogues (cbl::catalogue::Catalogue catalogue, cbl::catalogue::Catalogue random_catalogue, const double split, const int seed)
{
  m_data = std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(catalogue)));
  m_random = std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(random_catalogue)));
  m_splitFactor = split;

  /// Some compatibility checks...
  if (m_splitFactor>0) {
    m_random->shuffle(seed);
    double nSplit = static_cast<double>(m_random->nObjects())/(m_data->nObjects()*m_splitFactor);

    double integerPart, fracPart;
    std::modf(nSplit, &integerPart);
    fracPart = nSplit-integerPart;
    m_nSplit = static_cast<size_t> (integerPart);

    if (integerPart<1) {
      ErrorCBL("random catalogue too small!", "set_catalogues", "ThreePointCorrelation_comoving_multipoles.cpp");
    } else if (fracPart > 0) {
      WarningMsgCBL("a fraction of random objects will not be used; the effective size of the random sample will be "+conv(m_nSplit*m_data->nObjects()*m_splitFactor, par::fINT)+"!="+conv(m_random->nObjects(), par::fINT), "set_catalogues", "ThreePointCorrelation_comoving_multipoles.cpp");
    }

    coutCBL << "Splitting activated! Three Point correlation function will be computed "<< m_nSplit << " times. Each time the random sample will have " << m_data->nObjects()*m_splitFactor << " objects." << endl;

  }
  else 
    m_nSplit = 1;

}
