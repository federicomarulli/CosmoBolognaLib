/********************************************************************
 *  Copyright (C) 2023 by Federico Marulli                          *
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

#include <iostream>
#include <memory>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Sigma
#include <boost/test/included/unit_test.hpp>

#include "Sigma.h"
#include "LCDM.h"

using namespace std;
using namespace cbl::cosmology;


// =========================================================


BOOST_AUTO_TEST_SUITE(Sigma)

BOOST_AUTO_TEST_CASE(sigma2M_function)
{
  auto cosmology = make_shared<LCDM>("Planck18");
  cbl::cosmology::Sigma SG(cosmology);
  
  const double mass = 1.e12;
  const double redshift = 0.5;

  double sigma2M_EH = SG.sigma2M({mass}, "EisensteinHu", redshift, true, "test", "Linear", 100., "NULL", true, false)[0];
  BOOST_REQUIRE_CLOSE(sigma2M_EH, 2.78454365, 1.e-3);
    
  double sigma2M_CAMB = SG.sigma2M({mass}, "CAMB", redshift, true, "test", "Linear", 100., "NULL", true, false)[0];
  BOOST_REQUIRE_CLOSE(sigma2M_CAMB, 2.80562257, 1.e-3);

  double sigma2M_CLASS = SG.sigma2M({mass}, "CLASS", redshift, true, "test", "Linear", 100., "NULL", true, false)[0];
  BOOST_REQUIRE_CLOSE(sigma2M_CLASS, 2.83418095, 1.e-3);

  cosmology->set_unit(false);
  cbl::cosmology::Sigma SG2(cosmology);

  sigma2M_EH = SG2.sigma2M({mass}, "EisensteinHu", redshift, true, "test", "Linear", 100., "NULL", true, false)[0];
  BOOST_REQUIRE_CLOSE(sigma2M_EH, 3.11871485, 1.e-3);

  sigma2M_CAMB = SG2.sigma2M({mass}, "CAMB", redshift, true, "test", "Linear", 100., "NULL", true, false)[0];
  BOOST_REQUIRE_CLOSE(sigma2M_CAMB, 3.14457084, 1.e-3);

  sigma2M_CLASS = SG2.sigma2M({mass}, "CLASS", redshift, true, "test", "Linear", 100., "NULL", true, false)[0];
  BOOST_REQUIRE_CLOSE(sigma2M_CLASS, 3.18256592, 1.e-3);
}


// =========================================================


BOOST_AUTO_TEST_CASE(sigma2R_function)
{
  auto cosmology = make_shared<LCDM>("Planck18");
  cbl::cosmology::Sigma SG(cosmology);
  
  const double radius = 10.;
  const double redshift = 0.5;

  double sigma2R_EH = SG.sigma2R({radius}, "EisensteinHu", redshift, true, "test", "Linear", 100., "NULL", true, false)[0];
  BOOST_REQUIRE_CLOSE(sigma2R_EH, 0.30290894, 1.e-3);

  double sigma2R_CAMB = SG.sigma2R({radius}, "CAMB", redshift, true, "test", "Linear", 100., "NULL", true, false)[0];
  BOOST_REQUIRE_CLOSE(sigma2R_CAMB, 0.29581015, 1.e-3);

  double sigma2R_CLASS = SG.sigma2R({radius}, "CLASS", redshift, true, "test", "Linear", 100., "NULL", true, false)[0];
  BOOST_REQUIRE_CLOSE(sigma2R_CLASS, 0.29595004, 1.e-3);

  cosmology->set_unit(false);
  cbl::cosmology::Sigma SG2(cosmology);
  
  sigma2R_EH = SG2.sigma2R({radius}, "EisensteinHu", redshift, true, "test", "Linear", 100., "NULL", true, false)[0];
  BOOST_REQUIRE_CLOSE(sigma2R_EH, 0.52166260, 1.e-3);

  sigma2R_CAMB = SG2.sigma2R({radius}, "CAMB", redshift, true, "test", "Linear", 100., "NULL", true, false)[0];
  BOOST_REQUIRE_CLOSE(sigma2R_CAMB, 0.51360216, 1.e-3);

  sigma2R_CLASS = SG2.sigma2R({radius}, "CLASS", redshift, true, "test", "Linear", 100., "NULL", true, false)[0];
  BOOST_REQUIRE_CLOSE(sigma2R_CLASS, 0.51370480, 1.e-3);
}

BOOST_AUTO_TEST_SUITE_END()
