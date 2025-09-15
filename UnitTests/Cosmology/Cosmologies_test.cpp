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
#define BOOST_TEST_MODULE Cosmology
#include <boost/test/included/unit_test.hpp>

#include "Cosmologies.h"

using namespace std;
using namespace cbl::cosmology;


// =========================================================


BOOST_AUTO_TEST_SUITE(Cosmologies)

BOOST_AUTO_TEST_CASE(default_constructor_0)
{
  const LCDM cosmology_test1;
  BOOST_CHECK_EQUAL(cosmology_test1.model(), "LCDM");
  
  const LnuCDM cosmology_test2;
  BOOST_CHECK_EQUAL(cosmology_test2.model(), "LnuCDM");
  
  const w0CDM cosmology_test3;
  BOOST_CHECK_EQUAL(cosmology_test3.model(), "w0CDM");
  
  const kLCDM cosmology_test4;
  BOOST_CHECK_EQUAL(cosmology_test4.model(), "kLCDM");
  
  const w0nuCDM cosmology_test5;
  BOOST_CHECK_EQUAL(cosmology_test5.model(), "w0nuCDM");
  
  const w0waCDM cosmology_test6;
  BOOST_CHECK_EQUAL(cosmology_test6.model(), "w0waCDM");
  
  const kLnuCDM cosmology_test7;
  BOOST_CHECK_EQUAL(cosmology_test7.model(), "kLnuCDM");
  
  const kw0CDM cosmology_test8;
  BOOST_CHECK_EQUAL(cosmology_test8.model(), "kw0CDM");
  
  const w0wanuCDM cosmology_test9;
  BOOST_CHECK_EQUAL(cosmology_test9.model(), "w0wanuCDM");
  
  const kw0nuCDM cosmology_test10;
  BOOST_CHECK_EQUAL(cosmology_test10.model(), "kw0nuCDM");
  
  const kw0waCDM cosmology_test11;
  BOOST_CHECK_EQUAL(cosmology_test11.model(), "kw0waCDM");
  
  const kw0wanuCDM cosmology_test12;
  BOOST_CHECK_EQUAL(cosmology_test12.model(), "kw0wanuCDM");
}

// =========================================================


BOOST_AUTO_TEST_CASE(default_constructor_1)
{
  const auto LCDM = make_shared<cbl::cosmology::LCDM>();
  const auto LnuCDM = make_shared<cbl::cosmology::LnuCDM>();
  const auto w0CDM = make_shared<cbl::cosmology::w0CDM>();
  const auto kLCDM = make_shared<cbl::cosmology::kLCDM>();
  const auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>();
  const auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>();
  const auto kLnuCDM = make_shared<cbl::cosmology::kLnuCDM>();
  const auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>();
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>();
  const auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>();
  const auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>();
  const auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>();

  const vector<string> name = {"LCDM", "LnuCDM", "w0CDM", "kLCDM", "w0nuCDM", "w0waCDM", "kLnuCDM", "kw0CDM", "w0wanuCDM", "kw0nuCDM", "kw0waCDM", "kw0wanuCDM"};

  const vector<shared_ptr<cbl::cosmology::Cosmology>> cosmology = {LCDM, LnuCDM, w0CDM, kLCDM, w0nuCDM, w0waCDM, kLnuCDM, kw0CDM, w0wanuCDM, kw0nuCDM, kw0waCDM, kw0wanuCDM};
 
  for (size_t i=0; i<name.size(); ++i)
    BOOST_CHECK_EQUAL(cosmology[i]->model(), name[i]);
}


// =========================================================


BOOST_AUTO_TEST_CASE(constructor_WMAP5)
{
  const auto LCDM = make_shared<cbl::cosmology::LCDM>("WMAP5");
  const auto LnuCDM = make_shared<cbl::cosmology::LnuCDM>("WMAP5");
  const auto w0CDM = make_shared<cbl::cosmology::w0CDM>("WMAP5");
  const auto kLCDM = make_shared<cbl::cosmology::kLCDM>("WMAP5");
  const auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>("WMAP5");
  const auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>("WMAP5");
  const auto kLnuCDM = make_shared<cbl::cosmology::kLnuCDM>("WMAP5");
  const auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>("WMAP5");
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("WMAP5");
  const auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>("WMAP5");
  const auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>("WMAP5");
  const auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>("WMAP5");

  const vector<shared_ptr<cbl::cosmology::Cosmology>> cosmology = {LCDM, LnuCDM, w0CDM, kLCDM, w0nuCDM, w0waCDM, kLnuCDM, kw0CDM, w0wanuCDM, kw0nuCDM, kw0waCDM, kw0wanuCDM};

  for (auto &&cm : cosmology) {
    BOOST_REQUIRE_CLOSE(cm->Omega_radiation(), 4.17356e-05, 1.e-4);
    BOOST_REQUIRE_CLOSE(cm->Omega_matter(), 0.257958264, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->Omega_baryon(), 0.0441, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->little_h(), 0.719, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->scalar_pivot(), 0.002, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->scalar_amp(), 2.41e-09, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->n_spec(), 0.963, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->tau(), 0.087, 1.e-6);

    if (cm->model()=="LnuCDM" or cm->model()=="kLnuCDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0nuCDM" or cm->model()=="Lw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_neutrinos(), 0., 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->massless_neutrinos(), 3.04, 1.e-6);
      BOOST_CHECK_EQUAL(cm->massive_neutrinos(), 0);
    }

    if (cm->model()=="w0CDM" or cm->model()=="w0nuCDM" or cm->model()=="w0waCDM" or cm->model()=="kw0CDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->w0(), -1., 1.e-6);
    }

    if (cm->model()=="kLCDM" or cm->model()=="kLnuCDM" or cm->model()=="kw0CDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_DE(), 0.741999999, 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->Omega_k(), 0., 1.e-6);
    }

    if (cm->model()=="w0waCDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->wa(), 0., 1.e-6);
    }
  }
}


// =========================================================


BOOST_AUTO_TEST_CASE(constructor_WMAP7)
{
  const auto LCDM = make_shared<cbl::cosmology::LCDM>("WMAP7");
  const auto LnuCDM = make_shared<cbl::cosmology::LnuCDM>("WMAP7");
  const auto w0CDM = make_shared<cbl::cosmology::w0CDM>("WMAP7");
  const auto kLCDM = make_shared<cbl::cosmology::kLCDM>("WMAP7");
  const auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>("WMAP7");
  const auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>("WMAP7");
  const auto kLnuCDM = make_shared<cbl::cosmology::kLnuCDM>("WMAP7");
  const auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>("WMAP7");
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("WMAP7");
  const auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>("WMAP7");
  const auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>("WMAP7");
  const auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>("WMAP7");

  const vector<shared_ptr<cbl::cosmology::Cosmology>> cosmology = {LCDM, LnuCDM, w0CDM, kLCDM, w0nuCDM, w0waCDM, kLnuCDM, kw0CDM, w0wanuCDM, kw0nuCDM, kw0waCDM, kw0wanuCDM};

  for (auto &&cm : cosmology) {
    BOOST_REQUIRE_CLOSE(cm->Omega_radiation(), 4.17356e-05, 1.e-4);
    BOOST_REQUIRE_CLOSE(cm->Omega_matter(), 0.272958264, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->Omega_baryon(), 0.0455, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->little_h(), 0.704, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->scalar_pivot(), 0.002, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->scalar_amp(), 2.43e-09, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->n_spec(), 0.968, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->tau(), 0.088, 1.e-6);
    
    if (cm->model()=="LnuCDM" or cm->model()=="kLnuCDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0nuCDM" or cm->model()=="Lw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_neutrinos(), 0., 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->massless_neutrinos(), 3.04, 1.e-6);
      BOOST_CHECK_EQUAL(cm->massive_neutrinos(), 0);
    }

    if (cm->model()=="w0CDM" or cm->model()=="w0nuCDM" or cm->model()=="w0waCDM" or cm->model()=="kw0CDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->w0(), -1., 1.e-6);
    }

    if (cm->model()=="kLCDM" or cm->model()=="kLnuCDM" or cm->model()=="kw0CDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_DE(), 0.726999999, 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->Omega_k(), 0., 1.e-6);
    }

    if (cm->model()=="w0waCDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->wa(), 0., 1.e-6);
    }
  }
}


// =========================================================


BOOST_AUTO_TEST_CASE(constructor_WMAP9)
{
  const auto LCDM = make_shared<cbl::cosmology::LCDM>("WMAP9");
  const auto LnuCDM = make_shared<cbl::cosmology::LnuCDM>("WMAP9");
  const auto w0CDM = make_shared<cbl::cosmology::w0CDM>("WMAP9");
  const auto kLCDM = make_shared<cbl::cosmology::kLCDM>("WMAP9");
  const auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>("WMAP9");
  const auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>("WMAP9");
  const auto kLnuCDM = make_shared<cbl::cosmology::kLnuCDM>("WMAP9");
  const auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>("WMAP9");
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("WMAP9");
  const auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>("WMAP9");
  const auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>("WMAP9");
  const auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>("WMAP9");

  const vector<shared_ptr<cbl::cosmology::Cosmology>> cosmology = {LCDM, LnuCDM, w0CDM, kLCDM, w0nuCDM, w0waCDM, kLnuCDM, kw0CDM, w0wanuCDM, kw0nuCDM, kw0waCDM, kw0wanuCDM};

  for (auto &&cm : cosmology) {
    BOOST_REQUIRE_CLOSE(cm->Omega_radiation(), 8.53085e-05, 1.e-4);
    BOOST_REQUIRE_CLOSE(cm->Omega_matter(), 0.2793, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->Omega_baryon(), 0.0463, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->little_h(), 0.7, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->scalar_pivot(), 0.002, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->scalar_amp(), 2.41e-09, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->n_spec(), 0.972, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->tau(), 0.089, 1.e-6);
    
    if (cm->model()=="LnuCDM" or cm->model()=="kLnuCDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0nuCDM" or cm->model()=="Lw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_neutrinos(), 0., 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->massless_neutrinos(), 3.04, 1.e-6);
      BOOST_CHECK_EQUAL(cm->massive_neutrinos(), 0);
    }

    if (cm->model()=="w0CDM" or cm->model()=="w0nuCDM" or cm->model()=="w0waCDM" or cm->model()=="kw0CDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->w0(), -1., 1.e-6);
    }

    if (cm->model()=="kLCDM" or cm->model()=="kLnuCDM" or cm->model()=="kw0CDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_DE(), 0.720999999, 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->Omega_k(), -0.00038530849, 1.e-6);
    }

    if (cm->model()=="w0waCDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->wa(), 0., 1.e-6);
    }
  }
}


// =========================================================


BOOST_AUTO_TEST_CASE(constructor_Planck13)
{
  const auto LCDM = make_shared<cbl::cosmology::LCDM>("Planck13");
  const auto LnuCDM = make_shared<cbl::cosmology::LnuCDM>("Planck13");
  const auto w0CDM = make_shared<cbl::cosmology::w0CDM>("Planck13");
  const auto kLCDM = make_shared<cbl::cosmology::kLCDM>("Planck13");
  const auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>("Planck13");
  const auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>("Planck13");
  const auto kLnuCDM = make_shared<cbl::cosmology::kLnuCDM>("Planck13");
  const auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>("Planck13");
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck13");
  const auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>("Planck13");
  const auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>("Planck13");
  const auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>("Planck13");

  const vector<shared_ptr<cbl::cosmology::Cosmology>> cosmology = {LCDM, LnuCDM, w0CDM, kLCDM, w0nuCDM, w0waCDM, kLnuCDM, kw0CDM, w0wanuCDM, kw0nuCDM, kw0waCDM, kw0wanuCDM};

  for (auto &&cm : cosmology) {
    BOOST_REQUIRE_CLOSE(cm->Omega_radiation(), 9.28656e-05, 1.e-4);
    BOOST_REQUIRE_CLOSE(cm->Omega_matter(), 0.315, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->Omega_baryon(), 0.0487, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->little_h(), 0.673, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->scalar_pivot(), 0.05, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->scalar_amp(), 2.196e-09, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->n_spec(), 0.9603, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->tau(), 0.089, 1.e-6);

    if (cm->model()=="LnuCDM" or cm->model()=="kLnuCDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0nuCDM" or cm->model()=="Lw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_neutrinos(), 0.00142227929, 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->massless_neutrinos(), 2.04, 1.e-6);
      BOOST_CHECK_EQUAL(cm->massive_neutrinos(), 1);
    }

    if (cm->model()=="w0CDM" or cm->model()=="w0nuCDM" or cm->model()=="w0waCDM" or cm->model()=="kw0CDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->w0(), -1., 1.e-6);
    }

    if (cm->model()=="kLCDM" or cm->model()=="kLnuCDM" or cm->model()=="kw0CDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_DE(), 0.684907134, 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->Omega_k(), 0., 1.e-6);
    }

    if (cm->model()=="w0waCDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->wa(), 0., 1.e-6);
    }
  }
}


// =========================================================


BOOST_AUTO_TEST_CASE(constructor_Planck15)
{
  const auto LCDM = make_shared<cbl::cosmology::LCDM>("Planck15");
  const auto LnuCDM = make_shared<cbl::cosmology::LnuCDM>("Planck15");
  const auto w0CDM = make_shared<cbl::cosmology::w0CDM>("Planck15");
  const auto kLCDM = make_shared<cbl::cosmology::kLCDM>("Planck15");
  const auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>("Planck15");
  const auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>("Planck15");
  const auto kLnuCDM = make_shared<cbl::cosmology::kLnuCDM>("Planck15");
  const auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>("Planck15");
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck15");
  const auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>("Planck15");
  const auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>("Planck15");
  const auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>("Planck15");

  const vector<shared_ptr<cbl::cosmology::Cosmology>> cosmology = {LCDM, LnuCDM, w0CDM, kLCDM, w0nuCDM, w0waCDM, kLnuCDM, kw0CDM, w0wanuCDM, kw0nuCDM, kw0waCDM, kw0wanuCDM};

  for (auto &&cm : cosmology) {
    BOOST_REQUIRE_CLOSE(cm->Omega_radiation(), 9.22554e-05, 1.e-4);
    BOOST_REQUIRE_CLOSE(cm->Omega_matter(), 0.3121, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->Omega_baryon(), 0.0488, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->little_h(), 0.6751, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->scalar_pivot(), 0.05, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->scalar_amp(), 2.13e-09, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->n_spec(), 0.9653, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->tau(), 0.063, 1.e-6);

    if (cm->model()=="LnuCDM" or cm->model()=="kLnuCDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0nuCDM" or cm->model()=="Lw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_neutrinos(), 0.00141344463, 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->massless_neutrinos(), 2.04, 1.e-6);
      BOOST_CHECK_EQUAL(cm->massive_neutrinos(), 1);
    }

    if (cm->model()=="w0CDM" or cm->model()=="w0nuCDM" or cm->model()=="w0waCDM" or cm->model()=="kw0CDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->w0(), -1., 1.e-6);
    }

    if (cm->model()=="kLCDM" or cm->model()=="kLnuCDM" or cm->model()=="kw0CDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_DE(), 0.6878077446, 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->Omega_k(), 0., 1.e-6);
    }

    if (cm->model()=="w0waCDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->wa(), 0., 1.e-6);
    }
  }
}


// =========================================================


BOOST_AUTO_TEST_CASE(constructor_Planck15_TT)
{
  const auto LCDM = make_shared<cbl::cosmology::LCDM>("Planck15_TT");
  const auto LnuCDM = make_shared<cbl::cosmology::LnuCDM>("Planck15_TT");
  const auto w0CDM = make_shared<cbl::cosmology::w0CDM>("Planck15_TT");
  const auto kLCDM = make_shared<cbl::cosmology::kLCDM>("Planck15_TT");
  const auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>("Planck15_TT");
  const auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>("Planck15_TT");
  const auto kLnuCDM = make_shared<cbl::cosmology::kLnuCDM>("Planck15_TT");
  const auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>("Planck15_TT");
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck15_TT");
  const auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>("Planck15_TT");
  const auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>("Planck15_TT");
  const auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>("Planck15_TT");

  const vector<shared_ptr<cbl::cosmology::Cosmology>> cosmology = {LCDM, LnuCDM, w0CDM, kLCDM, w0nuCDM, w0waCDM, kLnuCDM, kw0CDM, w0wanuCDM, kw0nuCDM, kw0waCDM, kw0wanuCDM};

  for (auto &&cm : cosmology) {
    BOOST_REQUIRE_CLOSE(cm->Omega_radiation(), 9.15033e-05, 1.e-4);
    BOOST_REQUIRE_CLOSE(cm->Omega_matter(), 0.308, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->Omega_baryon(), 0.0484, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->little_h(), 0.6781, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->scalar_pivot(), 0.05, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->scalar_amp(), 2.139e-09, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->n_spec(), 0.9677, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->tau(), 0.066, 1.e-6);

    if (cm->model()=="LnuCDM" or cm->model()=="kLnuCDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0nuCDM" or cm->model()=="Lw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_neutrinos(), 0.00140096578, 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->massless_neutrinos(), 2.04, 1.e-6);
      BOOST_CHECK_EQUAL(cm->massive_neutrinos(), 1);
    }

    if (cm->model()=="w0CDM" or cm->model()=="w0nuCDM" or cm->model()=="w0waCDM" or cm->model()=="kw0CDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->w0(), -1., 1.e-6);
    }

    if (cm->model()=="kLCDM" or cm->model()=="kLnuCDM" or cm->model()=="kw0CDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_DE(), 0.6919084967, 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->Omega_k(), 0., 1.e-6);
    }

    if (cm->model()=="w0waCDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->wa(), 0., 1.e-6);
    }
  }
}


// =========================================================


BOOST_AUTO_TEST_CASE(constructor_Planck18)
{
  const auto LCDM = make_shared<cbl::cosmology::LCDM>("Planck18");
  const auto LnuCDM = make_shared<cbl::cosmology::LnuCDM>("Planck18");
  const auto w0CDM = make_shared<cbl::cosmology::w0CDM>("Planck18");
  const auto kLCDM = make_shared<cbl::cosmology::kLCDM>("Planck18");
  const auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>("Planck18");
  const auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>("Planck18");
  const auto kLnuCDM = make_shared<cbl::cosmology::kLnuCDM>("Planck18");
  const auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>("Planck18");
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");
  const auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>("Planck18");
  const auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>("Planck18");
  const auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>("Planck18");

  const vector<shared_ptr<cbl::cosmology::Cosmology>> cosmology = {LCDM, LnuCDM, w0CDM, kLCDM, w0nuCDM, w0waCDM, kLnuCDM, kw0CDM, w0wanuCDM, kw0nuCDM, kw0waCDM, kw0wanuCDM};

  for (auto &&cm : cosmology) {
    BOOST_REQUIRE_CLOSE(cm->Omega_radiation(), 9.26535e-05, 1.e-4);
    BOOST_REQUIRE_CLOSE(cm->Omega_matter(), 0.3153, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->Omega_baryon(), 0.0486, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->little_h(), 0.6736, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->scalar_pivot(), 0.05, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->scalar_amp(), 2.1e-09, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->n_spec(), 0.9649, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->tau(), 0.0544, 1.e-6);

    if (cm->model()=="LnuCDM" or cm->model()=="kLnuCDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0nuCDM" or cm->model()=="Lw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_neutrinos(), 0.001419746673, 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->massless_neutrinos(), 2.04, 1.e-6);
      BOOST_CHECK_EQUAL(cm->massive_neutrinos(), 1);
    }

    if (cm->model()=="w0CDM" or cm->model()=="w0nuCDM" or cm->model()=="w0waCDM" or cm->model()=="kw0CDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->w0(), -1., 1.e-6);
    }

    if (cm->model()=="kLCDM" or cm->model()=="kLnuCDM" or cm->model()=="kw0CDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_DE(), 0.6846073464, 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->Omega_k(), 0., 1.e-6);
    }

    if (cm->model()=="w0waCDM" or cm->model()=="w0wanuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->wa(), 0., 1.e-6);
    }
  }
}


// =========================================================


BOOST_AUTO_TEST_CASE(set_Omegas)
{
  const auto LCDM = make_shared<cbl::cosmology::LCDM>("Planck18");
  const auto LnuCDM = make_shared<cbl::cosmology::LnuCDM>("Planck18");
  const auto w0CDM = make_shared<cbl::cosmology::w0CDM>("Planck18");
  const auto kLCDM = make_shared<cbl::cosmology::kLCDM>("Planck18");
  const auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>("Planck18");
  const auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>("Planck18");
  const auto kLnuCDM = make_shared<cbl::cosmology::kLnuCDM>("Planck18");
  const auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>("Planck18");
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");
  const auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>("Planck18");
  const auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>("Planck18");
  const auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>("Planck18");

  const vector<shared_ptr<cbl::cosmology::Cosmology>> cosmology = {LCDM, LnuCDM, w0CDM, kLCDM, w0nuCDM, w0waCDM, kLnuCDM, kw0CDM, w0wanuCDM, kw0nuCDM, kw0waCDM, kw0wanuCDM};

  const std::vector<double> Omega1 = {0.0486, 0.233, 2.e-5};
  const std::vector<double> Omega2 = {0.0486, 0.233, 2.e-5, 1.e-5};
  const std::vector<double> Omega3 = {0.0486, 0.233, 2.e-5, -0.03162};
  const std::vector<double> Omega4 = {0.0486, 0.233, 2.e-5, -0.03163, 1.e-5};
  
  for (auto &&cm : cosmology) {
    if (cm->model()=="LCDM" or cm->model()=="w0CDM" or cm->model()=="w0waCDM") cm->set_parameters({"Omega_baryon","Omega_CDM","Omega_radiation"}, Omega1);
    else if (cm->model()=="LnuCDM" or cm->model()=="w0nuCDM" or cm->model()=="w0wanuCDM") cm->set_parameters({"Omega_baryon","Omega_CDM","Omega_radiation","Omega_neutrinos"}, Omega2);
    else if (cm->model()=="kLCDM" or cm->model()=="kw0CDM" or cm->model()=="kw0waCDM") cm->set_parameters({"Omega_baryon","Omega_CDM","Omega_radiation","Omega_k"}, Omega3);
    else if (cm->model()=="kLnuCDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0wanuCDM") cm->set_parameters({"Omega_baryon","Omega_CDM","Omega_radiation","Omega_k","Omega_neutrinos"}, Omega4);
    
    BOOST_REQUIRE_CLOSE(cm->Omega_baryon(), 0.0486, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->Omega_CDM(), 0.233, 1.e-6);
    BOOST_REQUIRE_CLOSE(cm->Omega_radiation(), 2.e-5, 1.e-6);

    if (cm->model()=="LCDM" or cm->model()=="w0CDM" or cm->model()=="w0waCDM")
      BOOST_REQUIRE_CLOSE(cm->Omega_matter(), cm->Omega_CDM()+cm->Omega_baryon(), 1.e-6);

    if (cm->model()=="LnuCDM" or cm->model()=="w0nuCDM" or cm->model()=="w0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_neutrinos(), 1.e-5, 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->Omega_matter(), cm->Omega_CDM()+cm->Omega_baryon()+cm->Omega_neutrinos(), 1.e-6);
    }
    
    else if (cm->model()=="kLCDM" or cm->model()=="kw0CDM" or cm->model()=="kw0waCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_DE(), 0.75, 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->Omega_k(), 1.-cm->Omega_matter()-cm->Omega_radiation()-cm->Omega_DE(), 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->Omega_matter(), cm->Omega_CDM()+cm->Omega_baryon(), 1.e-6);
    }
    
    else if (cm->model()=="kLnuCDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0wanuCDM") {
      BOOST_REQUIRE_CLOSE(cm->Omega_neutrinos(), 1.e-5, 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->Omega_DE(), 0.75, 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->Omega_matter(), cm->Omega_CDM()+cm->Omega_baryon()+cm->Omega_neutrinos(), 1.e-6);
      BOOST_REQUIRE_CLOSE(cm->Omega_k(), 1.-cm->Omega_matter()-cm->Omega_radiation()-cm->Omega_DE(), 1.e-6);
    }
    
  }
}


// =========================================================


BOOST_AUTO_TEST_CASE(set_Omega_baryon)
{
  const auto LCDM = make_shared<cbl::cosmology::LCDM>("Planck18");
  const auto LnuCDM = make_shared<cbl::cosmology::LnuCDM>("Planck18");
  const auto w0CDM = make_shared<cbl::cosmology::w0CDM>("Planck18");
  const auto kLCDM = make_shared<cbl::cosmology::kLCDM>("Planck18");
  const auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>("Planck18");
  const auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>("Planck18");
  const auto kLnuCDM = make_shared<cbl::cosmology::kLnuCDM>("Planck18");
  const auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>("Planck18");
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");
  const auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>("Planck18");
  const auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>("Planck18");
  const auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>("Planck18");

  const vector<shared_ptr<cbl::cosmology::Cosmology>> cosmology = {LCDM, LnuCDM, w0CDM, kLCDM, w0nuCDM, w0waCDM, kLnuCDM, kw0CDM, w0wanuCDM, kw0nuCDM, kw0waCDM, kw0wanuCDM};

  for (auto &&cm : cosmology) {
    cm->set_parameter("Omega_baryon",0.0486);
    BOOST_REQUIRE_CLOSE(cm->Omega_baryon(), 0.0486, 1.e-6);
    
    if (cm->model()!="LnuCDM" and cm->model()!="w0nuCDM" and cm->model()!="kLnuCDM" and cm->model()!="w0wanuCDM" and cm->model()!="kw0nuCDM" and cm->model()!="kw0wanuCDM") BOOST_REQUIRE_CLOSE(cm->Omega_CDM(), 0.2667, 1.e-6);
    else BOOST_REQUIRE_CLOSE(cm->Omega_CDM(), 0.2652802533, 1.e-6);
  }
}


// =========================================================


BOOST_AUTO_TEST_CASE(set_Omega_matter)
{
  const auto LCDM = make_shared<cbl::cosmology::LCDM>("Planck18");
  const auto LnuCDM = make_shared<cbl::cosmology::LnuCDM>("Planck18");
  const auto w0CDM = make_shared<cbl::cosmology::w0CDM>("Planck18");
  const auto kLCDM = make_shared<cbl::cosmology::kLCDM>("Planck18");
  const auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>("Planck18");
  const auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>("Planck18");
  const auto kLnuCDM = make_shared<cbl::cosmology::kLnuCDM>("Planck18");
  const auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>("Planck18");
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");
  const auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>("Planck18");
  const auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>("Planck18");
  const auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>("Planck18");

  const vector<shared_ptr<cbl::cosmology::Cosmology>> cosmology = {LCDM, LnuCDM, w0CDM, kLCDM, w0nuCDM, w0waCDM, kLnuCDM, kw0CDM, w0wanuCDM, kw0nuCDM, kw0waCDM, kw0wanuCDM};

  for (auto &&cm : cosmology) {
    cm->set_parameter("Omega_CDM",0.26);
    BOOST_REQUIRE_CLOSE(cm->Omega_CDM(), 0.26, 1.e-6);

    if (cm->model().find('k') != std::string::npos)
      cm->set_parameter("Omega_k",0.0553);
      
    if (cm->model()=="kLCDM" or cm->model()=="kLnuCDM" or cm->model()=="kw0CDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") BOOST_REQUIRE_CLOSE(cm->Omega_k(), 0.0553, 1.e-6);
  }
}


// =========================================================


BOOST_AUTO_TEST_CASE(set_Omega_radiation)
{
  const auto LCDM = make_shared<cbl::cosmology::LCDM>("Planck18");
  const auto LnuCDM = make_shared<cbl::cosmology::LnuCDM>("Planck18");
  const auto w0CDM = make_shared<cbl::cosmology::w0CDM>("Planck18");
  const auto kLCDM = make_shared<cbl::cosmology::kLCDM>("Planck18");
  const auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>("Planck18");
  const auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>("Planck18");
  const auto kLnuCDM = make_shared<cbl::cosmology::kLnuCDM>("Planck18");
  const auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>("Planck18");
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");
  const auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>("Planck18");
  const auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>("Planck18");
  const auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>("Planck18");

  const vector<shared_ptr<cbl::cosmology::Cosmology>> cosmology = {LCDM, LnuCDM, w0CDM, kLCDM, w0nuCDM, w0waCDM, kLnuCDM, kw0CDM, w0wanuCDM, kw0nuCDM, kw0waCDM, kw0wanuCDM};

  for (auto &&cm : cosmology) {
    cm->set_parameter("Omega_radiation",9.2e-05);
    BOOST_REQUIRE_CLOSE(cm->Omega_radiation(), 9.2e-05, 1.e-6);

    if (cm->model().find('k') != std::string::npos)
      cm->set_parameter("Omega_k",6.5354099e-7);
    
    if (cm->model()=="kLCDM" or cm->model()=="kLnuCDM" or cm->model()=="kw0CDM" or cm->model()=="kw0nuCDM" or cm->model()=="kw0waCDM" or cm->model()=="kw0wanuCDM") BOOST_REQUIRE_CLOSE(cm->Omega_k(), 6.5354099e-7, 1.e-6);
  }
}


// =========================================================


BOOST_AUTO_TEST_CASE(Omega_matter)
{
  w0CDM Planck18("Planck18");
  
  BOOST_REQUIRE_CLOSE(Planck18.Omega_matter(1.5), 0.87742609561413765, 1.e-6);
  BOOST_REQUIRE_CLOSE(Planck18.Omega_matter(1.5), Planck18.Omega_w(Planck18.Omega_matter(0.), 0., 1.5), 1.e-6);
}


// =========================================================


BOOST_AUTO_TEST_CASE(Omega_radiation)
{
  w0CDM Planck18("Planck18");
  
  BOOST_REQUIRE_CLOSE(Planck18.Omega_radiation(1.5), 0.00064459748428896387, 1.e-6);
  BOOST_REQUIRE_CLOSE(Planck18.Omega_radiation(1.5), Planck18.Omega_w(Planck18.Omega_radiation(0.), 1./3., 1.5), 1.e-6);
}


// =========================================================


BOOST_AUTO_TEST_CASE(f_DE)
{
  const auto w0CDM = make_shared<cbl::cosmology::w0CDM>("Planck18");
  const auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>("Planck18");
  const auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>("Planck18");
  const auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>("Planck18");
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");
  const auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>("Planck18");
  const auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>("Planck18");
  const auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>("Planck18");

  const vector<shared_ptr<cbl::cosmology::Cosmology>> cosmology_1 = {w0CDM, w0nuCDM, kw0CDM, kw0nuCDM};
  const vector<shared_ptr<cbl::cosmology::Cosmology>> cosmology_2 = {w0waCDM, w0wanuCDM, kw0waCDM, kw0wanuCDM};

  for (auto &&cm : cosmology_1) {
    cm->set_parameter("w0",-1.1);
    BOOST_REQUIRE_CLOSE(cm->f_DE(1.), 0.81225239635623536, 1.e-6);
  }

  for (auto &&cm : cosmology_2) {
    cm->set_parameter("w0",-1.1);
    cm->set_parameter("wa",0.1);
    BOOST_REQUIRE_CLOSE(cm->f_DE(1.), 0.86070797642505759, 1.e-6);
  }
  
}


// =========================================================


BOOST_AUTO_TEST_CASE(Omega_DE)
{
  const auto LCDM = make_shared<cbl::cosmology::LCDM>("Planck18");
  const auto LnuCDM = make_shared<cbl::cosmology::LnuCDM>("Planck18");
  const auto w0CDM = make_shared<cbl::cosmology::w0CDM>("Planck18");
  const auto kLCDM = make_shared<cbl::cosmology::kLCDM>("Planck18");
  const auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>("Planck18");
  const auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>("Planck18");
  const auto kLnuCDM = make_shared<cbl::cosmology::kLnuCDM>("Planck18");
  const auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>("Planck18");
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");
  const auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>("Planck18");
  const auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>("Planck18");
  const auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(LCDM->Omega_DE(1.1), 0.18983109370981133, 1.e-6);
  
  BOOST_REQUIRE_CLOSE(LnuCDM->Omega_DE(1.1), 0.18983109370981133, 1.e-6);

  w0CDM->set_parameter("w0",-1.1);
  BOOST_REQUIRE_CLOSE(w0CDM->Omega_DE(1.1), 0.15793292442721263, 1.e-6);
  
  BOOST_REQUIRE_CLOSE(kLCDM->Omega_DE(1.1), 0.18983109370981133, 1.e-6);

  w0nuCDM->set_parameter("w0",-1.1);
  BOOST_REQUIRE_CLOSE(w0nuCDM->Omega_DE(1.1), 0.15793292442721263, 1.e-6);

  w0waCDM->set_parameter("w0",-1.1);
  w0waCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0waCDM->Omega_DE(1.1), 0.16683159292913258, 1.e-6);

  BOOST_REQUIRE_CLOSE(kLnuCDM->Omega_DE(1.1), 0.18983109370981133, 1.e-6);

  kw0CDM->set_parameter("w0",-1.1);
  BOOST_REQUIRE_CLOSE(kw0CDM->Omega_DE(1.1), 0.15793292442721263, 1.e-6);

  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_DE(1.1), 0.16683159292913258, 1.e-6);

  kw0nuCDM->set_parameter("w0",-1.1);
  BOOST_REQUIRE_CLOSE(kw0nuCDM->Omega_DE(1.1), 0.15793292442721263, 1.e-6);

  kw0waCDM->set_parameter("w0",-1.1);
  kw0waCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(kw0waCDM->Omega_DE(1.1), 0.1668315929291326, 1.e-6);

  kw0wanuCDM->set_parameter("w0",-1.1);
  kw0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(kw0wanuCDM->Omega_DE(1.1), 0.1668315929291326, 1.e-6); 
}


// =========================================================


BOOST_AUTO_TEST_CASE(EE)
{
  const auto LCDM = make_shared<cbl::cosmology::LCDM>("Planck18");
  const auto LnuCDM = make_shared<cbl::cosmology::LnuCDM>("Planck18");
  const auto w0CDM = make_shared<cbl::cosmology::w0CDM>("Planck18");
  const auto kLCDM = make_shared<cbl::cosmology::kLCDM>("Planck18");
  const auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>("Planck18");
  const auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>("Planck18");
  const auto kLnuCDM = make_shared<cbl::cosmology::kLnuCDM>("Planck18");
  const auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>("Planck18");
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");
  const auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>("Planck18");
  const auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>("Planck18");
  const auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(LCDM->EE(1.1), 1.8990530750322905, 1.e-6);
  
  BOOST_REQUIRE_CLOSE(LnuCDM->EE(1.1), 1.8990530750322905, 1.e-6);

  w0CDM->set_parameter("w0",-1.1);
  BOOST_REQUIRE_CLOSE(w0CDM->EE(1.1), 1.8627370155534999, 1.e-6);
  
  BOOST_REQUIRE_CLOSE(kLCDM->EE(1.1), 1.8990530750322905, 1.e-6);

  w0nuCDM->set_parameter("w0",-1.1);
  BOOST_REQUIRE_CLOSE(w0nuCDM->EE(1.1), 1.8627370155534999, 1.e-6);

  w0waCDM->set_parameter("w0",-1.1);
  w0waCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0waCDM->EE(1.1), 1.8726580916072342, 1.e-6);

  BOOST_REQUIRE_CLOSE(kLnuCDM->EE(1.1), 1.8990530750322905, 1.e-6);

  kw0CDM->set_parameter("w0",-1.1);
  BOOST_REQUIRE_CLOSE(kw0CDM->EE(1.1), 1.8627370155534999, 1.e-6);

  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->EE(1.1), 1.8726580916072342, 1.e-6);

  kw0nuCDM->set_parameter("w0",-1.1);
  BOOST_REQUIRE_CLOSE(kw0nuCDM->EE(1.1), 1.8627370155534999, 1.e-6);

  kw0waCDM->set_parameter("w0",-1.1);
  kw0waCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(kw0waCDM->EE(1.1), 1.8726580916072342, 1.e-6);

  kw0wanuCDM->set_parameter("w0",-1.1);
  kw0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(kw0wanuCDM->EE(1.1), 1.8726580916072342, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(Hubble)
{
  const auto LCDM = make_shared<cbl::cosmology::LCDM>("Planck18");
  const auto LnuCDM = make_shared<cbl::cosmology::LnuCDM>("Planck18");
  const auto w0CDM = make_shared<cbl::cosmology::w0CDM>("Planck18");
  const auto kLCDM = make_shared<cbl::cosmology::kLCDM>("Planck18");
  const auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>("Planck18");
  const auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>("Planck18");
  const auto kLnuCDM = make_shared<cbl::cosmology::kLnuCDM>("Planck18");
  const auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>("Planck18");
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");
  const auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>("Planck18");
  const auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>("Planck18");
  const auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>("Planck18");

  LCDM->set_unit(false);
  BOOST_REQUIRE_CLOSE(LCDM->Hubble(1.1), 127.92021513417508, 1.e-6);

  LnuCDM->set_unit(false);
  BOOST_REQUIRE_CLOSE(LnuCDM->Hubble(1.1), 127.92021513417508, 1.e-6);

  w0CDM->set_unit(false);
  w0CDM->set_parameter("w0",-1.1);
  BOOST_REQUIRE_CLOSE(w0CDM->Hubble(1.1), 125.47396536768375, 1.e-6);

  kLCDM->set_unit(false);
  BOOST_REQUIRE_CLOSE(kLCDM->Hubble(1.1), 127.92021513417508, 1.e-6);

  w0nuCDM->set_unit(false);
  w0nuCDM->set_parameter("w0",-1.1);
  BOOST_REQUIRE_CLOSE(w0nuCDM->Hubble(1.1), 125.47396536768375, 1.e-6);

  w0waCDM->set_unit(false);
  w0waCDM->set_parameter("w0",-1.1);
  w0waCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0waCDM->Hubble(1.1), 126.1422490506633, 1.e-6);

  kLnuCDM->set_unit(false);
  BOOST_REQUIRE_CLOSE(kLnuCDM->Hubble(1.1), 127.92021513417508, 1.e-6);

  kw0CDM->set_unit(false);
  kw0CDM->set_parameter("w0",-1.1);
  BOOST_REQUIRE_CLOSE(kw0CDM->Hubble(1.1), 125.47396536768375, 1.e-6);

  w0wanuCDM->set_unit(false);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Hubble(1.1), 126.1422490506633, 1.e-6);

  kw0nuCDM->set_unit(false);
  kw0nuCDM->set_parameter("w0",-1.1);
  BOOST_REQUIRE_CLOSE(kw0nuCDM->Hubble(1.1), 125.47396536768375, 1.e-6);

  kw0waCDM->set_unit(false);
  kw0waCDM->set_parameter("w0",-1.1);
  kw0waCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(kw0waCDM->Hubble(1.1), 126.1422490506633, 1.e-6);

  kw0wanuCDM->set_unit(false);
  kw0wanuCDM->set_parameter("w0",-1.1);
  kw0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(kw0wanuCDM->Hubble(1.1), 126.1422490506633, 1.e-6);

}


BOOST_AUTO_TEST_SUITE_END()


