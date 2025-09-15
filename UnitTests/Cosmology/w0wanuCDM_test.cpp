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

#include "w0wanuCDM.h"

using namespace std;
using namespace cbl::cosmology;


BOOST_AUTO_TEST_SUITE(w0wanuCDM)

// =========================================================

BOOST_AUTO_TEST_CASE(default_constructor)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>();
  BOOST_CHECK_EQUAL(w0wanuCDM->model(), "w0wanuCDM");
}


// =========================================================


BOOST_AUTO_TEST_CASE(constructor_WMAP5)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("WMAP5");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_radiation(), 4.17356e-05, 1.e-4);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_matter(), 0.257958264, 1.e-6);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_baryon(), 0.0441, 1.e-6);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->little_h(), 0.719, 1.e-6);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->scalar_pivot(), 0.002, 1.e-6);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->scalar_amp(), 2.41e-09, 1.e-6);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->n_spec(), 0.963, 1.e-6);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->tau(), 0.087, 1.e-6);  
  w0wanuCDM->Omega_matter();
  w0wanuCDM->Omega_k();
  w0wanuCDM->Omega_DE();
  w0wanuCDM->Omega_neutrinos();
    
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_neutrinos(), 0., 1.e-6);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->massless_neutrinos(), 3.04, 1.e-6);
  BOOST_CHECK_EQUAL(w0wanuCDM->massive_neutrinos(), 0);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->w0(), -1., 1.e-6);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->wa(), 0., 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(constructor_WMAP7)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("WMAP7");

    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_radiation(), 4.17356e-05, 1.e-4);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_matter(), 0.272958264, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_baryon(), 0.0455, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->little_h(), 0.704, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->scalar_pivot(), 0.002, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->scalar_amp(), 2.43e-09, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->n_spec(), 0.968, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->tau(), 0.088, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_neutrinos(), 0., 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->massless_neutrinos(), 3.04, 1.e-6);
    BOOST_CHECK_EQUAL(w0wanuCDM->massive_neutrinos(), 0);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->w0(), -1., 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->wa(), 0., 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(constructor_WMAP9)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("WMAP9");

    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_radiation(), 8.53085e-05, 1.e-4);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_matter(), 0.2793, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_baryon(), 0.0463, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->little_h(), 0.7, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->scalar_pivot(), 0.002, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->scalar_amp(), 2.41e-09, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->n_spec(), 0.972, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->tau(), 0.089, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_neutrinos(), 0., 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->massless_neutrinos(), 3.04, 1.e-6);
    BOOST_CHECK_EQUAL(w0wanuCDM->massive_neutrinos(), 0);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->w0(), -1., 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->wa(), 0., 1.e-6);
}


// =========================================================


BOOST_AUTO_TEST_CASE(constructor_Planck13)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck13");

    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_radiation(), 9.28656e-05, 1.e-4);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_matter(), 0.315, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_baryon(), 0.0487, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->little_h(), 0.673, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->scalar_pivot(), 0.05, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->scalar_amp(), 2.196e-09, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->n_spec(), 0.9603, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->tau(), 0.089, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_neutrinos(), 0.00142227929, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->massless_neutrinos(), 2.04, 1.e-6);
    BOOST_CHECK_EQUAL(w0wanuCDM->massive_neutrinos(), 1);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->w0(), -1., 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->wa(), 0., 1.e-6);
}


// =========================================================


BOOST_AUTO_TEST_CASE(constructor_Planck15)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck15");

    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_radiation(), 9.22554e-05, 1.e-4);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_matter(), 0.3121, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_baryon(), 0.0488, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->little_h(), 0.6751, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->scalar_pivot(), 0.05, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->scalar_amp(), 2.13e-09, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->n_spec(), 0.9653, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->tau(), 0.063, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_neutrinos(), 0.00141344463, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->massless_neutrinos(), 2.04, 1.e-6);
    BOOST_CHECK_EQUAL(w0wanuCDM->massive_neutrinos(), 1);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->w0(), -1., 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->wa(), 0., 1.e-6);
}


// =========================================================


BOOST_AUTO_TEST_CASE(constructor_Planck15_TT)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck15_TT");

    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_radiation(), 9.15033e-05, 1.e-4);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_matter(), 0.308, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_baryon(), 0.0484, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->little_h(), 0.6781, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->scalar_pivot(), 0.05, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->scalar_amp(), 2.139e-09, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->n_spec(), 0.9677, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->tau(), 0.066, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_neutrinos(), 0.00140096578, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->massless_neutrinos(), 2.04, 1.e-6);
    BOOST_CHECK_EQUAL(w0wanuCDM->massive_neutrinos(), 1);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->w0(), -1., 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->wa(), 0., 1.e-6);
}


// =========================================================


BOOST_AUTO_TEST_CASE(constructor_Planck18)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_radiation(), 9.26535e-05, 1.e-4);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_matter(), 0.3153, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_baryon(), 0.0486, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->little_h(), 0.6736, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->scalar_pivot(), 0.05, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->scalar_amp(), 2.1e-09, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->n_spec(), 0.9649, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->tau(), 0.0544, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_neutrinos(), 0.001419746673, 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->massless_neutrinos(), 2.04, 1.e-6);
    BOOST_CHECK_EQUAL(w0wanuCDM->massive_neutrinos(), 1);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->w0(), -1., 1.e-6);
    BOOST_REQUIRE_CLOSE(w0wanuCDM->wa(), 0., 1.e-6);
}


// =========================================================


BOOST_AUTO_TEST_CASE(set_Omegas)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  const std::vector<double> Omegas = {0.0486, 0.233, 2.e-5, 1.e-5};
  
  w0wanuCDM->set_parameters({"Omega_baryon","Omega_CDM","Omega_radiation","Omega_neutrinos"}, Omegas);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_baryon(), 0.0486, 1.e-6);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_CDM(), 0.233, 1.e-6);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_radiation(), 2.e-5, 1.e-6);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_neutrinos(), 1.e-5, 1.e-6);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_matter(), w0wanuCDM->Omega_CDM()+w0wanuCDM->Omega_baryon()+w0wanuCDM->Omega_neutrinos(), 1.e-6);
 
}

// =========================================================


BOOST_AUTO_TEST_CASE(set_Omega_baryon)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  w0wanuCDM->set_parameter("Omega_baryon",0.0486);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_baryon(), 0.0486, 1.e-6);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_CDM(), 0.2652802533, 1.e-6);
  
}


// =========================================================


BOOST_AUTO_TEST_CASE(set_Omega_matter)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  w0wanuCDM->set_parameter("Omega_CDM",0.26);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_matter(), 0.31001974667353033, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(set_Omega_radiation)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  w0wanuCDM->set_parameter("Omega_radiation",9.2e-05);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_radiation(), 9.2e-05, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(set_Omega_neutrinos)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");
  const double Omega_Nu{0.001419746673}, massless{2.04};
  const int massive{1};

  w0wanuCDM->set_parameter("massless_neutrinos", massless);
  w0wanuCDM->set_parameter("massive_neutrinos",massive);
  w0wanuCDM->set_parameter("Omega_neutrinos",Omega_Nu);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_neutrinos(), 0.001419746673, 1.e-6);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->massless_neutrinos(), 2.04, 1.e-6);
  BOOST_CHECK_EQUAL(w0wanuCDM->massive_neutrinos(), 1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_CDM(), w0wanuCDM->Omega_matter()-w0wanuCDM->Omega_baryon()-w0wanuCDM->Omega_neutrinos(), 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(set_w0)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  w0wanuCDM->set_parameter("w0",-1.);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->w0(), -1., 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(set_wa)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  w0wanuCDM->set_parameter("wa",0.);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->wa(), 0., 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(Omega_matter)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");
  
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_matter(1.), 0.7861642563274404, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(Omega_radiation)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");
  
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_radiation(1.), 0.000462041878535081, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(f_DE)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->f_DE(1.), 1., 1.e-6);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->f_DE(1.), 0.86070797642505759, 1.e-6);

  
}


// =========================================================


BOOST_AUTO_TEST_CASE(Omega_DE)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_DE(1.), 0.21337370179402448, 1.e-6);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_DE(1.), 0.1892780275597836, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(neutrino_mass)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->neutrino_mass(), 0.06, 1.e-4);

}


// =========================================================


BOOST_AUTO_TEST_CASE(Omega_Nu)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_Nu(1.), 0.003539974905706805, 1.e-6);
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Omega_Nu(1.), 0.0034908892411867616, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(EE)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->EE(1.), 1.7912257822828754, 1.e-6);
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->EE(1.), 1.7643663120802111, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(Hubble)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  w0wanuCDM->set_unit(false);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Hubble(1.1), 126.1422490506633, 1.e-6);

}

// =========================================================


BOOST_AUTO_TEST_CASE(linear_growth_rate)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->linear_growth_rate(1.), 0.8764079469816577, 1.e-6);
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->linear_growth_rate(1.), 0.89123107436855409, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(D_C)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->D_C(1.), 2291.83192879128, 1.e-6);
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->D_C(1.), 2326.6433916619712, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(D_M)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_CHECK_EQUAL(w0wanuCDM->D_M(1.), w0wanuCDM->D_C(1.));

}


// =========================================================


BOOST_AUTO_TEST_CASE(D_A)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->D_A(1.,0.), 1145.91596439564, 1.e-6);
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->D_A(1.,0.), 1163.3216958309856, 1.e-6);

}

// =========================================================


BOOST_AUTO_TEST_CASE(qq)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->qq(1.), 0.1801704682482308, 1.e-6);
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->qq(1.), 0.20210311519792779, 1.e-6);

}

// =========================================================


BOOST_AUTO_TEST_CASE(z_acc)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->z_acc(), 0.6315564992697738, 1.e-6);
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->z_acc(), 0.5732264832525702, 1.e-6);

}

// =========================================================


BOOST_AUTO_TEST_CASE(z_eq)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->z_eq(), 0.29496725164319604, 1.e-6);
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->z_eq(), 0.26792877096434342, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(z_eq_rad)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->z_eq_rad(), 3444.5972266223116, 1.e-6);
  w0wanuCDM->set_parameter("Omega_CDM",0.32);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->z_eq_rad(), 4042.4008664355629, 1.e-6);

}

// =========================================================


BOOST_AUTO_TEST_CASE(Mag_Volume_limited)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->Mag_Volume_limited(), -63.30606380597709, 1.e-6);
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Mag_Volume_limited(), -63.33879909522279, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(Lum_bol)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->Lum_bol(1.), 264019124.3972557, 1.e-6);
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Lum_bol(1.), 272100601.68040013, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(Redshift)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->Redshift(), 0.0003335904227546314, 1.e-6);
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Redshift(), 0.0003335847112521807, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(Redshift_time)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->Redshift_time(1.,1.,2.), 1.9998779296875, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(Volume)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->Volume(1.,2.,1.), 3432036.2946982444, 1.e-6); 
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Volume(1.), 50423954925.428955, 1.e-6);
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Volume(1.,2.,1.), 3568522.084401865, 1.e-6);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Volume(1.), 52756755563.912437, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(max_redshift)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->max_redshift(5.e5,1.,1.), 1.174406486683231, 1.e-6);
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->max_redshift(5.e5,1.,1.), 1.1672729134141602, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(rho_m)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->rho_m(1.), 696904221361.8322, 1.e-6); 
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->rho_m(1.), 696904221361.83228, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(Delta_c)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->Delta_c(1.), 158.38063054442432, 1.e-6); 
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->Delta_c(1.), 160.73350251164788, 1.e-6);

}


// ========================================================


BOOST_AUTO_TEST_CASE(sound_horizon_at_drag_epoch)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->sound_horizon_at_drag_epoch(), 99.13265602367865, 1.e-6); 
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->sound_horizon_at_drag_epoch(), 99.134230002413133, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_CASE(As)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->As(0.81), 2.045147767767742e-09, 1.e-10); 
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->As(0.81), 1.8977153401027489e-09, 1.e-10);

}


// =========================================================


BOOST_AUTO_TEST_CASE(sigma8_interpolated)
{
  const auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

  BOOST_REQUIRE_CLOSE(w0wanuCDM->sigma8_interpolated(0.), 0.8128679966040663, 1.e-6); 
  w0wanuCDM->set_parameter("massless_neutrinos",3.04);
  w0wanuCDM->set_parameter("massive_neutrinos",0);
  w0wanuCDM->set_parameter("Omega_neutrinos",0.0014);
  w0wanuCDM->set_parameter("w0",-1.1);
  w0wanuCDM->set_parameter("wa",0.1);
  BOOST_REQUIRE_CLOSE(w0wanuCDM->sigma8_interpolated(0.), 0.82612849303404834, 1.e-6);

}


// =========================================================


BOOST_AUTO_TEST_SUITE_END()


