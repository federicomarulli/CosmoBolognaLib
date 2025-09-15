/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  Modelling/ShearShearLensing/ModelFunction_PowerSpectrum_Shear.cpp
 *
 *  @brief Functions to model the shear-shear angular power spectrum
 *
 *  This file contains the implementation of the functions used to
 *  model the shear-shear angular power spectrum

 */

#include "ModelFunction_PowerSpectrum_Angular.h"
#include "ModelFunction_PowerSpectrum_Shear.h"
#include "Modelling_PowerSpectrum_Shear.h"

using namespace std;
using namespace cbl;


// ===============================================================================================

   double cbl::modelling::sslensing::window_function_IA (const double zz, const std::shared_ptr<cbl::cosmology::Cosmology> cosmo, const std::shared_ptr<cbl::glob::FuncGrid> nz)
   { 
     double W_IA_val; 
     W_IA_val = cosmo->Hubble(zz)*nz->operator()(zz)*pow(cbl::par::cc,-1);     

     return W_IA_val;
   }
   
// ============================================================================================


   double cbl::modelling::sslensing::window_function_shear (const double zz, const double z_max, const std::shared_ptr<cbl::cosmology::Cosmology> cosmo, const std::shared_ptr<cbl::glob::FuncGrid> nz, cbl::glob::FuncGrid DC )
   {
     auto integrand = [&zz,&nz,&DC] (const double z)
     {
       double N_z = nz->operator()(z);
       double D_ = DC(zz)/DC(z);
       return N_z*(1.-D_);
     };
     double integral = wrapper::gsl::GSL_integrate_cquad(integrand, zz, z_max);
     double W_shear_val = 1.5*pow(cosmo->Hubble(0.)*pow(cbl::par::cc,-1),2)*cosmo->Omega_matter(0.)*(1.0 + zz)*cosmo->D_C(zz)*integral;

   return W_shear_val;
   }

// ============================================================================================

   double cbl::modelling::sslensing::Intrinsic_Alignment (const double zz,const double Om_0,const double z_p, const std::vector<double> pars, cbl::glob::FuncGrid DN, cbl::glob::FuncGrid L_ratio)
   {  
     const double precfactor = -1.*pars[0]*pars[1]*Om_0/DN(zz);
     const double z_ratio = pow((1+zz)*(1+z_p),pars[2]);
     const double LL_ratio = pow(L_ratio(zz),pars[3]);
     double f_IA = precfactor*z_ratio*LL_ratio;
      
     return f_IA;
   }
   
// ============================================================================================

    std::vector<double> cbl::modelling::sslensing::C_ij_LL_eNLA (const std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
    { 
      // structure contaning the required input data
      shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

      // redefine the cosmology
      auto cosmo = pp->cosmology->clone();
  
      // set the cosmological parameters 
      for (size_t i=0; i<pp->Cpar.size(); ++i)
        cosmo->set_parameter(pp->Cpar[i], parameter[i]);
    
      vector<double> z_vector = linear_bin_vector(150,pp->z_min,pp->z_max);
      vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);

      cosmology::PkXi PX(cosmo);
      vector<vector<double>> Pk = PX.Pk_matter(kk, pp->method_Pk,pp->NL, z_vector, false, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec,pp->file_par, pp->unit);
      vector<double> Pk0 = PX.Pk_matter(kk, pp->method_Pk,pp->NL, 0, false,pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec,pp->file_par, pp->unit); 
      glob::FuncGrid2D pk_interp(z_vector,kk, Pk, "Cubic");  
  
      vector<double> comoving_distance, hubble,DN_z,f_IA_z;
      for (size_t i=0; i<z_vector.size(); ++i)
      {
        hubble.emplace_back(cosmo->Hubble(z_vector[i]));
        comoving_distance.emplace_back(cosmo->D_C(z_vector[i]));  
        DN_z.emplace_back(sqrt(Pk[i][0]/Pk0[0]));
      }

      glob::FuncGrid HH(z_vector, hubble, "Spline");
      glob::FuncGrid DC(z_vector, comoving_distance, "Spline");
      glob::FuncGrid DN(z_vector, DN_z, "Spline");  
  
      for (size_t i=0; i<z_vector.size(); ++i) 
        f_IA_z.emplace_back(Intrinsic_Alignment(z_vector[i],cosmo->Omega_matter(0.),0.,pp->pars_IA,DN,pp->L_ratio));

      glob::FuncGrid f_IA(z_vector, f_IA_z, "Spline");
  
      int n_size = 0.5*pp->Nz_all.size()*(pp->Nz_all.size()-1)+pp->Nz_all.size();
      int l_size = l.size()/n_size;
      vector<double> ells(l_size,0.), pre_factor(l_size,0.);
      vector<double> Cl;
      for (int i=0; i<l_size;++i)
      {
        ells[i] = l[i*n_size];
        pre_factor[i] = (ells[i] + 2.0) * (ells[i] + 1.0) * ells[i] * (ells[i] - 1.0)/pow((ells[i] + 0.5),4);
      }
 
      vector<vector<double>> WW_IA (pp->Nz_all.size());
      vector<vector<double>> WW_shear (pp->Nz_all.size());
      for (size_t j=0; j<WW_IA.size(); ++j)
      {
        for (size_t i=0; i<z_vector.size(); ++i)
        {
          WW_IA[j].emplace_back(window_function_IA(z_vector[i],cosmo,pp->Nz_all[j]));
          WW_shear[j].emplace_back(window_function_shear(z_vector[i],pp->z_max,cosmo,pp->Nz_all[j],DC));
        }
      }
  
      vector<glob::FuncGrid> Wij_IA (pp->Nz_all.size());
      vector<glob::FuncGrid> Wij_shear (pp->Nz_all.size());
      for (size_t i=0; i<Wij_IA.size(); ++i)
      {
        Wij_IA[i]=glob::FuncGrid(z_vector, WW_IA[i], "Spline");
        Wij_shear[i]=glob::FuncGrid(z_vector, WW_shear[i], "Spline");
      }
  
      for (size_t i=0; i<ells.size(); ++i)
      {
        for (size_t j=0; j<pp->Nz_all.size(); ++j)
        {
          for (size_t k=j; k<pp->Nz_all.size();++k)
            Cl.push_back(pre_factor[i]*cbl::modelling::angularpk::C_ij_xy(ells[i],pp->z_min,pp->z_max,Wij_IA[j],Wij_IA[k],Wij_shear[j],Wij_shear[k],HH, DC, pk_interp, f_IA));
        }
      }

      return Cl;
    }


