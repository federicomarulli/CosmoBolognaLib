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

#include "Catalogue.h"
#include "Cosmologies.h"
#include "MapCalculator.h"
#include "Measure.h"

namespace std{
  ostream& operator<<(std::ostream& os, const std::vector<std::vector<double>>& vec) {
      for (const auto& inner_vec : vec) {
          for (const auto& value : inner_vec) {
              os << value << " ";
          }
          os << std::endl;
      }
      return os;
  }
}

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MapCalculator
#include <boost/test/included/unit_test.hpp>

using namespace std;
using namespace cbl::measure;


// =========================================================

BOOST_AUTO_TEST_SUITE(MapCalculator)

BOOST_AUTO_TEST_CASE(constructor_1)
{
    cbl::measure::MapCalculator map_test;

    std::vector<double> v1 (1, 1.);
    std::vector<double> v2 (1, 0.3);
    map_test.setNside(1);
    map_test.setNsigma(5);
    map_test.setWidths(v1);
    map_test.setZobs(v2);
    map_test.setScheme(RING);
    map_test.setSeed(0);

    const auto cosmo = make_shared<cbl::cosmology::LCDM>();
    std::string file_catalogue = "input_test.dat";
    cbl::catalogue::Catalogue cat (cbl::catalogue::ObjectType::_Mock_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmo, cbl::CoordinateUnits::_degrees_ );
    map_test.setCatalogue(cat);

    BOOST_CHECK_EQUAL(map_test.getNside(), 1); 
    BOOST_CHECK_EQUAL(map_test.getScheme(), RING);  
    BOOST_CHECK_EQUAL(map_test.getNsigma(), 5);
    
    for (size_t i=0; i<v1.size(); ++i) {
      BOOST_CHECK_EQUAL(map_test.getWidths()[i], v1[i]); 
      BOOST_CHECK_EQUAL(map_test.getZobs()[i], v2[i]);
    }
}

// =========================================================

BOOST_AUTO_TEST_CASE(constructor_2)
{
    const auto cosmo = make_shared<cbl::cosmology::LCDM>("Planck18");
    std::string file_catalogue = "input_test.dat";
    cbl::catalogue::Catalogue cat (cbl::catalogue::ObjectType::_Mock_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmo, cbl::CoordinateUnits::_degrees_ );
    
    std::vector<double> v1 (1, 1.);
    std::vector<double> v2 (1, 0.3);
    cbl::measure::MapCalculator map_test (1, 5, v1 , v2, cat);

    BOOST_CHECK_EQUAL(map_test.getNside(), 1); 
    BOOST_CHECK_EQUAL(map_test.getScheme(), RING);  
    BOOST_CHECK_EQUAL(map_test.getNsigma(), 5);

    for (size_t i=0; i<v1.size(); ++i) {
      BOOST_CHECK_EQUAL(map_test.getWidths()[i], v2[i]); 
      BOOST_CHECK_EQUAL(map_test.getZobs()[i], v1[i]);
    }
}

// =========================================================

BOOST_AUTO_TEST_CASE(constructor_3)
{
    const auto cosmo = make_shared<cbl::cosmology::LCDM>("Planck18");
    std::string file_catalogue = "input_test.dat";
    cbl::catalogue::Catalogue cat (cbl::catalogue::ObjectType::_Mock_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmo, cbl::CoordinateUnits::_degrees_ );
    
    std::vector<double> v1 (1, 1.);
    std::vector<double> v2 (1, 0.3);
    cbl::measure::MapCalculator map_test (1, RING, 5, v1 , v2, 0, cat);

    BOOST_CHECK_EQUAL(map_test.getNside(), 1); 
    BOOST_CHECK_EQUAL(map_test.getSeed(), 0);  
    BOOST_CHECK_EQUAL(map_test.getNsigma(), 5);
    for (size_t i=0; i<v1.size(); ++i) {
      BOOST_CHECK_EQUAL(map_test.getWidths()[i], v2[i]); 
      BOOST_CHECK_EQUAL(map_test.getZobs()[i], v1[i]);
    }
}

// =========================================================


BOOST_AUTO_TEST_CASE(HealpixMapping)
{
    const auto cosmo = make_shared<cbl::cosmology::LCDM>("Planck18");
    std::string file_catalogue = "input_test.dat";
    cbl::catalogue::Catalogue cat (cbl::catalogue::ObjectType::_Mock_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmo, cbl::CoordinateUnits::_degrees_ );
    
    std::vector<double> v1 (1, 1.);
    std::vector<double> v2 (1, 0.3);
    cbl::measure::MapCalculator map_test (1, 5, v1 , v2, cat);

    auto result = map_test.HealpixMapping();
    BOOST_TEST(result.size() == (map_test.getCatalogue())->nObjects());
}


// =========================================================


BOOST_AUTO_TEST_CASE(ADF_ARF)
{
    const auto cosmo = make_shared<cbl::cosmology::LCDM>("Planck18");
    string file_catalogue = "input_test.dat";
    cbl::catalogue::Catalogue cat (cbl::catalogue::ObjectType::_Mock_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmo, cbl::CoordinateUnits::_degrees_ );
    
    vector<double> v1 (1, 1.);
    vector<double> v2 (1, 0.3);
    cbl::measure::MapCalculator map_test (1, 5, v1 , v2, cat);

    vector<double> mskin (12 * map_test.getNside() * map_test.getNside(), 0.97);
    auto result = map_test.ADF(mskin);
    auto result2 = map_test.ADF_ARF(mskin);
    
    BOOST_TEST(result.size() == get<1>(result2).size());
    for (size_t i=0; i<result.size(); ++i) {
      BOOST_CHECK_EQUAL(result[i], get<0>(result2)[i]);
    }
}


// =========================================================


BOOST_AUTO_TEST_CASE(ADFr_ARFr)
{
    const auto cosmo = make_shared<cbl::cosmology::LCDM>("Planck18");
    string file_catalogue = "input_test.dat";
    cbl::catalogue::Catalogue cat (cbl::catalogue::ObjectType::_Mock_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmo, cbl::CoordinateUnits::_degrees_ );
    
    vector<double> v1 (1, 1.);
    vector<double> v2 (1, 0.3);
    cbl::measure::MapCalculator map_test (1, 5, v1 , v2, cat);

    vector<double> mskin (12 * map_test.getNside() * map_test.getNside(), 0.97);
    vector<double> zH (12 * map_test.getNside() * map_test.getNside(), 0.02);

    auto result = map_test.ADFr(mskin, zH);
    auto result2 = map_test.ADFr_ARFr(mskin, zH);
    BOOST_TEST(result.size() == get<1>(result2).size());
    
    for (size_t i=0; i<result.size(); ++i) {
      BOOST_CHECK_EQUAL(result[i], get<0>(result2)[i]);
    }
}


// =========================================================


BOOST_AUTO_TEST_CASE(SN)
{
    const auto cosmo = make_shared<cbl::cosmology::LCDM>("Planck18");
    string file_catalogue = "input_test.dat";
    cbl::catalogue::Catalogue cat (cbl::catalogue::ObjectType::_Mock_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmo, cbl::CoordinateUnits::_degrees_ );
    
    vector<double> v1 (1, 1.);
    vector<double> v2 (1, 0.3);
    cbl::measure::MapCalculator map_test (1, RING, 5, v1 , v2, 0, cat);

    vector<double> mskin (12 * map_test.getNside() * map_test.getNside(), 0.97);

    auto result = map_test.SN(mskin, 2);

    double elements = 0;
    for (size_t i=0; i<result.size(); ++i) {
      for (size_t j=0; j<result[i].size(); ++j) {
        for (size_t k=0; k<result[i][j].size(); ++k) {
          for (size_t l=0; l<result[i][j][k].size(); ++l) {
            for (size_t m=0; m<result[i][j][k][l].size(); ++m) {
              ++elements;
            }
          }
        }
      }
    }

    auto expected = (12 * map_test.getNside() * map_test.getNside() * map_test.getZobs().size() * map_test.getWidths().size() * 2 * 2);
    BOOST_CHECK_EQUAL(elements, expected);
}

BOOST_AUTO_TEST_SUITE_END()


