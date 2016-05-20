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
 *  @file Cosmology/Lib/PkXi.cpp
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
 *  @author federico.marulli3@unbo.it
 */

#include "Cosmology.h"
using namespace cosmobl;


// =====================================================================================


double cosmobl::Cosmology::As (const double sigma8) const
{ 
  double zero = 0.;
  return pow(sigma8/1.79e4*pow(m_Omega_baryon*m_hh*m_hh/0.024,1./3.)*pow(m_Omega_matter*m_hh*m_hh/0.14,-0.563)*pow(7.808*m_hh,0.5*(1.-m_n_spec))*pow(m_hh/0.72,-0.693)*0.76/gg(zero),2);
}


// =====================================================================================


string cosmobl::Cosmology::Pk_output_file (const string code, const bool NL, const double redshift, const bool run, const string output_root, const double k_max, const string file_par)
{
  string dir_loc = fullpath(par::DirLoc);
  string dir_cosmo = fullpath(par::DirCosmo);
   
  string sdir;
  if (NL==0) sdir = "output_linear/";
  else if (NL==1) sdir= "output_nonlinear/";
  else ErrorMsg("Error in cosmobl::Cosmology::Table_PkCodes of PkXi.cpp!");
  sdir += "h"+conv(m_hh,par::fDP3)+"_OmB"+conv(m_Omega_baryon,par::fDP3)+"_OmCDM"+conv(m_Omega_CDM,par::fDP3)+"_OmL"+conv(m_Omega_DE,par::fDP3)+"_OmN"+conv(m_Omega_neutrinos,par::fDP3)+"_Z"+conv(redshift,par::fDP3)+"_scalar_amp"+conv(m_scalar_amp,par::ee3)+"_n"+conv(m_n_spec,par::fDP3)+"/";

  string dir = dir_cosmo+"External/"+code+"/";
  string dirC = dir_cosmo+"External/CAMB/";
  if (chdir (dirC.c_str())) {};

  string dir_output = dir+sdir+"pk.dat";

  if (run) {
    vector<double> lgkk, lgPk;
    Table_PkCodes(code, NL, lgkk, lgPk, redshift, output_root, k_max, file_par);   
  }

  return dir_output;
}


// =====================================================================================


double cosmobl::Cosmology::Pk_UnNorm (const double kk, const double redshift, const string method_Pk) const
{ 
  double Pk = -1.;
  if (method_Pk=="EisensteinHu") { // Eisenstein & Hu  
    TFmdm_set_cosm(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massive_neutrinos, m_Omega_DE, m_hh, redshift); 
    double func = (m_unit) ? TFmdm_onek_hmpc(kk) : TFmdm_onek_mpc(kk); // transfer function
    Pk = pow(kk,m_n_spec)*pow(func,2);
  }
  
  if (Pk<0) ErrorMsg("Error in cosmobl::Cosmology::Pk_UnNorm of PkXi.cpp!");

  return Pk;
}


// =====================================================================================


void cosmobl::Cosmology::run_CAMB (const bool NL, const double redshift, const string output_root, const double k_max, const string file_par) const 
{
  string dir = par::DirCosmo+"External/CAMB/";

  string File_par = file_par;

  if (File_par==par::defaultString) {
    
    // --------------------------------------------------------------------------
    // --------- set the cosmological parameters in the file params.ini ---------
    // --------------------------------------------------------------------------

    string file_par_default = dir+"params.ini";
    File_par = dir+"params_"+output_root+".ini";
  
    string CP = "cp "+file_par_default+" "+File_par; if (system (CP.c_str())) {};

    double HH0 = m_hh*100.;

    string sed;
    sed = "sed '/test/s//"+output_root+"/g' "+File_par+" > temp; mv temp "+File_par; if (system (sed.c_str())) {};
    sed = "sed '/do_nonlinear = 0/s//do_nonlinear = "+conv(NL,par::fINT)+"/g' "+File_par+" > temp; mv temp "+File_par+""; if (system (sed.c_str())) {};
    sed = "sed '/hubble = 70/s//hubble = "+conv(HH0,par::fDP6)+"/g' "+File_par+" > temp; mv temp "+File_par+""; if (system (sed.c_str())) {};
    sed = "sed '/omega_baryon = 0.0462/s//omega_baryon = "+conv(m_Omega_baryon,par::fDP6)+"/g' "+File_par+" > temp; mv temp "+File_par+""; if (system (sed.c_str())) {};
    sed = "sed '/omega_cdm = 0.2538/s//omega_cdm = "+conv(m_Omega_CDM,par::fDP6)+"/g' "+File_par+" > temp; mv temp "+File_par+""; if (system (sed.c_str())) {};
    sed = "sed '/omega_lambda = 0.7/s//omega_lambda = "+conv(m_Omega_DE,par::fDP6)+"/g' "+File_par+" > temp; mv temp "+File_par+""; if (system (sed.c_str())) {};
    sed = "sed '/omega_neutrino = 0/s//omega_neutrino = "+conv(m_Omega_neutrinos,par::fDP6)+"/g' "+File_par+" > temp; mv temp "+File_par+""; if (system (sed.c_str())) {};
    sed = "sed '/transfer_redshift(1) = 0/s//transfer_redshift(1) = "+conv(redshift,par::fDP6)+"/g' "+File_par+" > temp; mv temp "+File_par; if (system (sed.c_str())) {};
    sed = "sed '/massless_neutrinos = 3.046/s//massless_neutrinos = "+conv(m_massless_neutrinos,par::fDP6)+"/g' "+File_par+" > temp; mv temp "+File_par; if (system (sed.c_str())) {};
    sed = "sed '/massive_neutrinos = 0/s//massive_neutrinos = "+conv(m_massive_neutrinos,par::fINT)+"/g' "+File_par+" > temp; mv temp "+File_par; if (system (sed.c_str())) {};
    sed = "sed '/scalar_spectral_index(1) = 0.96/s//scalar_spectral_index(1) = "+conv(m_n_spec,par::fDP6)+"/g' "+File_par+" > temp; mv temp "+File_par; if (system (sed.c_str())) {};
    sed = "sed '/w = -1/s//w = "+conv(m_w0,par::fDP2)+"/g' "+File_par+" > temp; mv temp "+File_par; if (system (sed.c_str())) {};
    sed = "sed '/wa = 0/s//wa = "+conv(m_wa,par::fDP2)+"/g' "+File_par+" > temp; mv temp "+File_par; if (system (sed.c_str())) {};
    if (m_scalar_amp>0) {sed = "sed '/scalar_amp(1) = 2.1e-9/s//scalar_amp(1) = "+conv(m_scalar_amp,par::ee3)+"/g' "+File_par+" > temp; mv temp "+File_par; if (system (sed.c_str())) {};}
    sed = "sed '/transfer_kmax = 2/s//transfer_kmax = "+conv(k_max,par::fDP2)+"/g' "+File_par+" > temp; mv temp "+File_par; if (system (sed.c_str())) {};
  
  }


  // --------------------------------------------------------------------------
  

  if (chdir (dir.c_str())) {};

  string Camb = "./camb "+File_par;
  if (system (Camb.c_str())) {};

  string RM = "rm -f "+File_par+" "+output_root+"_params.ini "+output_root+"_transfer_out.dat "+output_root+"_matterpower.dat"; 
  if (system (RM.c_str())) {};

  if (chdir(par::DirLoc.c_str())) {};
}


// =====================================================================================


void cosmobl::Cosmology::Table_PkCodes (const string code, const bool NL, vector<double> &lgkk, vector<double> &lgPk, const double redshift, const string output_root, const double k_max, string file_par) const 
{ 
  if (code=="MPTbreeze-v1" && m_sigma8<0) 
    ErrorMsg("Error in cosmobl::Cosmology::Table_PkCodes of PkXi.cpp: sigma8 must be >0 if MPTbreeze-v1 is used!"); 
  
  lgkk.erase(lgkk.begin(), lgkk.end());
  lgPk.erase(lgPk.begin(), lgPk.end());
  
  string dir_loc = fullpath(par::DirLoc);
  string dir_cosmo = fullpath(par::DirCosmo);
  string sdir;
  if (NL==0) sdir = "output_linear/";
  else if (NL==1) sdir= "output_nonlinear/";
  else ErrorMsg("Error in cosmobl::Cosmology::Table_PkCodes of PkXi.cpp!");
  sdir += "h"+conv(m_hh,par::fDP3)+"_OmB"+conv(m_Omega_baryon,par::fDP3)+"_OmCDM"+conv(m_Omega_CDM,par::fDP3)+"_OmL"+conv(m_Omega_DE,par::fDP3)+"_OmN"+conv(m_Omega_neutrinos,par::fDP3)+"_Z"+conv(redshift,par::fDP3)+"_scalar_amp"+conv(m_scalar_amp,par::ee3)+"_n"+conv(m_n_spec,par::fDP3)+"/";

  string dir = dir_cosmo+"External/"+code+"/";
  string dirC = dir_cosmo+"External/CAMB/";
  if (chdir (dirC.c_str())) {};

  string dir_output = dir+sdir;
  
  string file_in = dir_output+"Pk.dat";
  ifstream fin;
  fin.open(file_in.c_str());

  string file_par0;

  if (!fin) {
   
    string dirCr = (code!="MPTbreeze-v1") ? dir : dirC;
    if (chdir (dirCr.c_str())) {};
    
    if (file_par==par::defaultString) {
      
      // --------------------------------------------------------------------------
      // --------- set the cosmological parameters in the file params.ini ---------
      // --------------------------------------------------------------------------

      string file_par_default = dirCr+"params.ini";
      file_par0 = "params_"+output_root+".ini";
      file_par = dirCr+"params_"+output_root+".ini";
      string CP = "cp "+file_par_default+" "+file_par; if (system (CP.c_str())) {};

      double HH0 = m_hh*100.;

      string sed;
    
      if (code=="CAMB" || code=="MPTbreeze-v1") {
      
	sed = "sed '/test/s//"+output_root+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};
	sed = "sed '/do_nonlinear = 0/s//do_nonlinear = "+conv(NL,par::fINT)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system (sed.c_str())) {};
	sed = "sed '/hubble = 70/s//hubble = "+conv(HH0,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system (sed.c_str())) {};
	sed = "sed '/omega_baryon = 0.0462/s//omega_baryon = "+conv(m_Omega_baryon,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system (sed.c_str())) {};
	sed = "sed '/omega_cdm = 0.2538/s//omega_cdm = "+conv(m_Omega_CDM,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system (sed.c_str())) {};
	sed = "sed '/omega_lambda = 0.7/s//omega_lambda = "+conv(m_Omega_DE,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system (sed.c_str())) {};
	sed = "sed '/omega_neutrino = 0/s//omega_neutrino = "+conv(m_Omega_neutrinos,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system (sed.c_str())) {};
	sed = "sed '/transfer_redshift(1) = 0/s//transfer_redshift(1) = "+conv(redshift,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};
	sed = "sed '/massless_neutrinos = 3.046/s//massless_neutrinos = "+conv(m_massless_neutrinos,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};
	sed = "sed '/massive_neutrinos = 0/s//massive_neutrinos = "+conv(m_massive_neutrinos,par::fINT)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};
	sed = "sed '/scalar_spectral_index(1) = 0.96/s//scalar_spectral_index(1) = "+conv(m_n_spec,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};
	sed = "sed '/w = -1/s//w = "+conv(m_w0,par::fDP2)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};
	sed = "sed '/wa = 0/s//wa = "+conv(m_wa,par::fDP2)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};
	if (m_scalar_amp>0) {sed = "sed '/scalar_amp(1) = 2.1e-9/s//scalar_amp(1) = "+conv(m_scalar_amp,par::ee3)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};}
	sed = "sed '/transfer_kmax = 2/s//transfer_kmax = "+conv(k_max,par::fDP2)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};
      
      }
    
      if (code=="classgal_v1") {
      
	sed = "sed '/output\\/test_/s//"+output_root+"_/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};
	sed = "sed '/non linear = /s//non linear = halofit/g' "+file_par+" > temp; mv temp "+file_par+""; if (system (sed.c_str())) {};
	sed = "sed '/h = 0.7/s//h = "+conv(m_hh,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system (sed.c_str())) {};
	sed = "sed '/Omega_b = 0.05/s//Omega_b = "+conv(m_Omega_baryon,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system (sed.c_str())) {};
	sed = "sed '/Omega_cdm = 0.25/s//Omega_cdm = "+conv(m_Omega_CDM,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system (sed.c_str())) {};
	sed = "sed '/Omega_Lambda = 0.7/s//Omega_Lambda = "+conv(m_Omega_DE,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system (sed.c_str())) {};
	sed = "sed '/Omega_k = 0./s//Omega_k = "+conv(m_Omega_k,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system (sed.c_str())) {};
	sed = "sed '/Omega_ncdm = /s//Omega_ncdm = "+conv(m_Omega_neutrinos,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system (sed.c_str())) {};
	sed = "sed '/z_pk = 0/s//z_pk = "+conv(redshift,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};
	sed = "sed '/N_eff = 3.04/s//N_eff = "+conv(m_massless_neutrinos,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};
	sed = "sed '/N_ncdm = 0/s//N_ncdm = "+conv(m_massive_neutrinos,par::fINT)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};
	sed = "sed '/n_s = 1./s//n_s = "+conv(m_n_spec,par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};

	double w00 = max(-0.999,m_w0); // check!!!
	sed = "sed '/w0_fld = -0.9/s//w0_fld = "+conv(w00,par::fDP3)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};

	sed = "sed '/wa_fld = 0./s//wa_fld = "+conv(m_wa,par::fDP2)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};
	if (m_scalar_amp>0) {sed = "sed '/A_s = 2.3e-9/s//A_s = "+conv(m_scalar_amp,par::ee3)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};}
	sed = "sed '/P_k_max_h\\/Mpc = 1/s//P_k_max_h\\/Mpc = "+conv(k_max,par::fDP2)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system (sed.c_str())) {};
      }
    }
    
    else {
      string CP = "cp "+file_par+" file_par_temp"; if (system (CP.c_str())) {};
      file_par = "file_par_temp";
    }
    
  
    // -----------------------------------------------------------------------------
    // --------- run the code to get P(k): either CAMB, MPTbreeze or CLASS ---------
    // -----------------------------------------------------------------------------


    string MK = "mkdir -p "+dir_output; 
    if (code=="MPTbreeze-v1") MK += " "+dirCr+sdir; 
    if (system (MK.c_str())) {};

    
    if (code=="CAMB" || code=="MPTbreeze-v1") {
    
      string GO = "./camb "+file_par;
      cout <<"---> "<<GO<<endl<<endl;
      if (system (GO.c_str())) {};
      
      if (code=="MPTbreeze-v1") {
	if (chdir (dir.c_str())) {};
	string GOM = "./mptbreeze -noverbose -camb ../CAMB/"+file_par0+" -fileTF ../CAMB/"+output_root+"_transfer_out.dat -sigma8 "+conv(m_sigma8,par::fDP5)+" -omegam "+conv(m_Omega_matter,par::fDP5)+" -ns "+conv(m_n_spec,par::fDP5)+" -w "+conv(m_w0,par::fDP5)+" -redshift "+conv(redshift,par::fDP5)+" -filePk "+output_root+"_matterpower.dat";	
	if (system (GOM.c_str())) {};
      }

      string MV = "mv "+output_root+"_matterpower*dat "+file_in+"; mv "+output_root+"_transfer*dat "+dir_output+"transfer_out.dat; mv "+file_par+" "+dir_output+"param.dat";
      if (system (MV.c_str())) {};

    }

    if (code=="classgal_v1") {
      if (chdir (dir.c_str())) {};
      string GO = "./class "+file_par;
      if (system (GO.c_str())) {};
      string MV = (NL==0) ? "mv "+output_root+"_pk.dat "+file_in : "mv "+output_root+"_pk_nl.dat "+file_in;
      if (system (MV.c_str())) {};   
    } 
  
    string RM = "rm -f "+file_par+" *_params.ini"; if (system (RM.c_str())) {};
  }
  fin.clear(); fin.close();


  // ---------------------------------------
  // --------- get the output P(k) ---------
  // ---------------------------------------

  if (chdir(dir_loc.c_str())) {};

  double KK, PK, PK0, PK1, PK2;
  
  if (code=="CAMB" || code=="MPTbreeze-v1") {

    string file_inCAMB = dirC+sdir+"Pk.dat";
    fin.open(file_inCAMB.c_str()); checkIO (file_inCAMB,1); 
    
    while (fin >>KK>>PK) 
      if (KK>0 && PK>0) {
	lgkk.push_back(log10(KK));
	lgPk.push_back(log10(PK));
      } 
    fin.clear(); fin.close();

    
    if (code=="MPTbreeze-v1") {

      vector<double> lgkkM, lgPkM;
      fin.open(file_in.c_str()); checkIO (file_in,1); 

      while (fin >>KK>>PK0>>PK1>>PK2) {
	if (KK>0 && PK0>0) {
	  lgkkM.push_back(log10(KK));
	  lgPkM.push_back(log10(pow(2.*par::pi,3)*(PK0+PK1+PK2)));
	} 
      }
      fin.clear(); fin.close();

      double lgm = Min(lgkkM), lgM = Max(lgkkM);
      
      for (unsigned int i=0; i<lgkk.size(); i++)
	if (lgm<lgkk[i] && lgkk[i]<lgM) 
	  lgPk[i] = interpolated(lgkk[i], lgkkM, lgPkM, "Linear");
    }
  }

  if (code=="classgal_v1") {
    fin.open(file_in.c_str()); checkIO (file_in,1); 
    string line; for (int i=0; i<4; i++) getline(fin,line);
    while (fin >>KK>>PK) 
      if (KK>0 && PK>0) {
	lgkk.push_back(log10(KK));
	lgPk.push_back(log10(PK));
      } 
    fin.clear(); fin.close();
  }
}


// =====================================================================================


void cosmobl::Cosmology::Table_XiCodes (const string code, const bool NL, vector<double> &rr, vector<double> &xi, const double redshift, const string output_root, const double k_max, string file_par) const
{
  vector<double> lgkk, lgPk;
  Table_PkCodes(code, NL, lgkk, lgPk, redshift, output_root, k_max, file_par);

  string dir_loc = fullpath(par::DirLoc);
  string dir_cosmo = fullpath(par::DirCosmo);

  string sdir;
  if (NL==0) sdir = "output_linear/";
  else if (NL==1) sdir= "output_nonlinear/";
  else ErrorMsg("Error in cosmobl::Cosmology::Table_XiCodes of PkXi.cpp!");
  sdir += "h"+conv(m_hh,par::fDP3)+"_OmB"+conv(m_Omega_baryon,par::fDP3)+"_OmCDM"+conv(m_Omega_CDM,par::fDP3)+"_OmL"+conv(m_Omega_DE,par::fDP3)+"_OmN"+conv(m_Omega_neutrinos,par::fDP3)+"_Z"+conv(redshift,par::fDP3)+"_scalar_amp"+conv(m_scalar_amp,par::ee3)+"_n"+conv(m_n_spec,par::fDP3)+"/";

  string dir = dir_cosmo+"External/"+code+"/";
  string dirFFT = dir_cosmo+"External/fftlog-f90-master/";
  if (chdir(dirFFT.c_str())) {};
  string dir_output = dir+sdir;

  string file_in = dir_output+"Pk.dat";
  string file_out = dir_output+"Xi.dat";

  ifstream fin;
  fin.open(file_out.c_str());
  
  if (!fin) {
    string cmd = "./fftlog-f90 "+file_in+" "+file_out+" 600";
    if(system(cmd.c_str())) {}
  }
  
  if (chdir(dir_loc.c_str())) {};

  fin.clear(); fin.close();

  rr.erase(rr.begin(), rr.end());
  xi.erase(xi.begin(), xi.end());

  fin.open(file_out.c_str()); 

  double RR, XXII;
  while (fin >>RR >> XXII) {
    rr.push_back(RR);
    xi.push_back(XXII);
  }

  fin.clear(); fin.close();
}


// =====================================================================================


void cosmobl::Cosmology::Pk_0 (const string method_Pk, const double redshift, const string output_root, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par)
{  
  if (m_sigma8<0) ErrorMsg("Error in cosmobl::Cosmology::Pk_0 of PkXi.cpp, sigma8<0!");

  double RR = 8.; // sigma_8 = sigma(8Mpc/h)
  double RHO = Rho(m_Omega_matter, m_Omega_neutrinos); 
  double MM = Mass(RR, RHO);

  bool NL = 0;
  double Int = -1.;

  if (GSL) {
    int limit_size = 1000;
    gsl_integration_workspace *ww = gsl_integration_workspace_alloc (limit_size);
    gsl_function Func;

    double error = -1.; 
 
    if (method_Pk=="EisensteinHu") {
      cosmobl::glob::STR_SSM_EH str;
      str.Omega_matter = m_Omega_matter; 
      str.Omega_baryon = m_Omega_baryon; 
      str.Omega_neutrinos = m_Omega_neutrinos; 
      str.massless_neutrinos = m_massless_neutrinos; 
      str.massive_neutrinos = m_massive_neutrinos; 
      str.Omega_DE = m_Omega_DE; 
      str.Omega_radiation = m_Omega_radiation; 
      str.hh = m_hh; 
      str.scalar_amp = m_scalar_amp; 
      str.n_spec = m_n_spec;
      str.w0 = m_w0; 
      str.wa = m_wa; 
      str.fNL = m_fNL;
      str.type_NG = m_type_NG;
      str.model = m_model;
      str.unit = m_unit;
      str.method_Pk = method_Pk;  
      str.redshift = redshift;    
      str.mass = MM;
      str.rho = RHO;

      Func.function = &glob::func_SSM_EH_GSL;
      Func.params = &str;
      gsl_integration_qag (&Func, k_min, k_max, 0., prec, limit_size, 6, ww, &Int, &error); 
    }

    else if (method_Pk=="CAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="classgal_v1") {
      vector<double> lgkk, lgPk;
      Table_PkCodes (method_Pk, NL, lgkk, lgPk, redshift, output_root, k_max, file_par);
   
      cosmobl::glob::STR_SSM str;
      str.unit = m_unit;
      str.hh = m_hh;
      str.n_spec = m_n_spec;
      str.mass = MM;
      str.rho = RHO;
      str.lgkk = lgkk;
      str.lgPk = lgPk;

      Func.function = &glob::func_SSM_GSL;
      Func.params = &str;
      gsl_integration_qag (&Func, k_min, k_max, 0., prec, limit_size, 6, ww, &Int, &error); 
    }
    
    else ErrorMsg("Error in cosmobl::Cosmology::sigma8_Pk of PkXi.cpp: method_Pk is wrong!");

    gsl_integration_workspace_free (ww);
  
    if (method_Pk=="EisensteinHu") m_Pk0_EH = 2.*pow(par::pi*m_sigma8,2)/Int*pow(DD(redshift),2);
    if (method_Pk=="CAMB") m_Pk0_CAMB = 2.*pow(par::pi*m_sigma8,2)/Int*pow(DD(redshift),2);
    if (method_Pk=="MPTbreeze-v1") m_Pk0_MPTbreeze = 2.*pow(par::pi*m_sigma8,2)/Int*pow(DD(redshift),2);
    if (method_Pk=="classgal_v1") m_Pk0_CLASS = 2.*pow(par::pi*m_sigma8,2)/Int*pow(DD(redshift),2);

  } 
  
  else { // using Numerical libraries

    if (method_Pk=="EisensteinHu") {
      cosmobl::classfunc::func_SSM func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, method_Pk, MM, redshift);
    
      Midpnt<cosmobl::classfunc::func_SSM> q1(func,0.,1.); 
      Midinf<cosmobl::classfunc::func_SSM> q2(func,1.,1.e30);
      Int = qromo(q1)+qromo(q2);
      m_Pk0_EH = 2.*pow(par::pi*m_sigma8,2)/Int*pow(DD(redshift),2);
    } 

    else if (method_Pk=="CAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="classgal_v1") {
      vector<double> lgkk, lgPk;
      Table_PkCodes (method_Pk, NL, lgkk, lgPk, redshift, output_root, k_max, file_par);
       
      cosmobl::classfunc::func_SSM_Table func (m_hh, m_n_spec, RHO, m_unit, lgkk, lgPk, MM); 

      Midpnt<cosmobl::classfunc::func_SSM_Table> q1(func,0.,1.); 
      Midinf<cosmobl::classfunc::func_SSM_Table> q2(func,1.,1.e30);
      Int = qromo(q1)+qromo(q2);

      if (method_Pk=="CAMB") m_Pk0_CAMB = 2.*pow(par::pi*m_sigma8,2)/Int*pow(DD(redshift),2);
      if (method_Pk=="MPTbreeze-v1") m_Pk0_MPTbreeze = 2.*pow(par::pi*m_sigma8,2)/Int*pow(DD(redshift),2);
      if (method_Pk=="classgal_v1") m_Pk0_CLASS = 2.*pow(par::pi*m_sigma8,2)/Int*pow(DD(redshift),2);
    }
  
    else ErrorMsg("Error in cosmobl::Cosmology::Pk_0 of Cosmology.cpp: method_Pk is wrong!");
  }

}


// =====================================================================================


double cosmobl::Cosmology::Pk (const double kk, const string method_Pk, const bool NL, const double redshift, const string output_root, const int norm, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par)
{ 
  int Norm = norm;
  if (Norm==-1) Norm = (m_sigma8>0) ? 1 : 0;
 
  if (method_Pk=="MPTbreeze-v1") Norm = 0; // check!!!

  if (Norm==1) Pk_0(method_Pk, redshift, output_root, k_min, k_max, GSL, prec, file_par);
  else { m_Pk0_EH = 1.; m_Pk0_CAMB = 1.; m_Pk0_MPTbreeze = 1.; m_Pk0_CLASS = 1.; }

  if (method_Pk=="EisensteinHu")  // NL is not used!!!
    return m_Pk0_EH*Pk_UnNorm(kk, redshift, method_Pk);
  
  if (method_Pk=="CAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="classgal_v1") {
    double fact = (m_unit) ? 1. : m_hh;
    double lgk = log10(kk/fact);
    vector<double> lgkk, lgPk;
    Table_PkCodes (method_Pk, NL, lgkk, lgPk, redshift, output_root, k_max, file_par);

    double lgPK = interpolated(lgk, lgkk, lgPk, "Linear");

    double PP0 = -1.;
    if (method_Pk=="CAMB") PP0 = m_Pk0_CAMB;
    if (method_Pk=="MPTbreeze-v1") PP0 = m_Pk0_MPTbreeze;
    if (method_Pk=="classgal_v1") PP0 = m_Pk0_CLASS;

    return PP0*pow(10., lgPK)/pow(fact, m_n_spec);
  }

  else { ErrorMsg("Error in cosmobl::Cosmology::Pk of PkXi.cpp: method_Pk is wrong!"); return 0; }

}


// =====================================================================================


void cosmobl::Cosmology::Pk_Kaiser_multipoles (vector<double> &Pk0, vector<double> &Pk2, vector<double> &Pk4, const vector<double> kk, const string method_Pk, const bool NL, const double redshift, const double bias, const double sigma_NL, const string output_root, const int norm, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par)
{ 
  vector<double> Pk_arr;
  size_t nbin_k = kk.size();
  double f = linear_growth_rate(redshift);

  if(sigma_NL==0){
    for(size_t i=0;i<nbin_k;i++)
      Pk_arr.push_back(Pk(kk[i],method_Pk,NL,redshift,output_root,norm,k_min,k_max,GSL,prec,file_par));
  }
  else{
    for(size_t i=0;i<nbin_k;i++)
      Pk_arr.push_back(Pk_DeWiggle(kk[i], redshift, sigma_NL , output_root, norm, k_min, k_max, prec));
  }

  vector<vector<double> > Pk_multipoles = cosmobl::Pkl_Kaiser({0,2,4}, kk, Pk_arr, bias, f);

  Pk0.erase(Pk0.begin(),Pk0.end());
  Pk2.erase(Pk2.begin(),Pk2.end());
  Pk4.erase(Pk4.begin(),Pk4.end());
  for(size_t i=0;i<nbin_k;i++){
    Pk0.push_back(Pk_multipoles[0][i]);
    Pk2.push_back(Pk_multipoles[1][i]);
    Pk4.push_back(Pk_multipoles[2][i]);
  }

}


// =====================================================================================

/// @cond glob

double cosmobl::glob::func_xi_EH_GSL (double kk, void *params) 
{
  struct cosmobl::glob::STR_xi_EH *pp = (struct cosmobl::glob::STR_xi_EH *) params;

  Cosmology cosm (pp->Omega_matter, pp->Omega_baryon, pp->Omega_neutrinos, pp->massless_neutrinos, pp->massive_neutrinos, pp->Omega_DE, pp->Omega_radiation, pp->hh, pp->scalar_amp, pp->n_spec, pp->w0, pp->wa, pp->fNL, pp->type_NG, pp->model, pp->unit);

  double Int = cosm.Pk_UnNorm(kk,pp->redshift,pp->method_Pk)*sin(kk*pp->rr)*kk/pp->rr;
 
  return Int * exp(-kk*kk*pp->aa*pp->aa); // eq. 24 of Anderson et al. 2012  
}

/// @endcond

// =====================================================================================


double cosmobl::Cosmology::xi_DM (const double rr, const string method_Pk, const double redshift, const string output_root, const bool NL, const int norm, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const string file_par) 
{
  int Norm = norm;
  if (Norm==-1) Norm = (m_sigma8>0) ? 1 : 0;

  if (method_Pk=="MPTbreeze-v1") Norm = 0; // check!!!
  if (method_Pk=="MPTbreeze-v1" && NL==0) ErrorMsg("Error in cosmobl::Cosmology::xi_DM of PkXi.cpp: MPTbreeze is non-linear!");  

  double Int = -1.; 
  double fact = (GSL) ? 1./(2.*pow(par::pi,2)) : 1.;

  if (GSL) {
    int limit_size = 1000;
    gsl_integration_workspace *ww = gsl_integration_workspace_alloc (limit_size);
    gsl_function Func;

    double error = -1.; 
 
    if (method_Pk=="EisensteinHu") {
   
      if (m_sigma8<0) ErrorMsg("Error in cosmobl::Cosmology::xi_DM of PkXi.cpp: sigma8<0!");
      if (NL==1) WarningMsg("Attention: the correlation function by Eisenstein&Hu is linear (see xi_DM of PkXi.cpp)!");
  
      cosmobl::glob::STR_xi_EH str;
      str.Omega_matter = m_Omega_matter; 
      str.Omega_baryon = m_Omega_baryon; 
      str.Omega_neutrinos = m_Omega_neutrinos; 
      str.massless_neutrinos = m_massless_neutrinos; 
      str.massive_neutrinos = m_massive_neutrinos; 
      str.Omega_DE = m_Omega_DE; 
      str.Omega_radiation = m_Omega_radiation; 
      str.hh = m_hh; 
      str.scalar_amp = m_scalar_amp; 
      str.n_spec = m_n_spec;
      str.w0 = m_w0; 
      str.wa = m_wa; 
      str.fNL = m_fNL;
      str.type_NG = m_type_NG;
      str.model = m_model;
      str.unit = m_unit;
      str.rr = rr;
      str.aa = aa;
      str.redshift = redshift;
      str.method_Pk = method_Pk;

      Func.function = &glob::func_xi_EH_GSL;
      Func.params = &str;
      gsl_integration_qag (&Func, k_min, k_max, 0., prec, limit_size, 6, ww, &Int, &error); 
    }

    else if (method_Pk=="CAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="classgal_v1") {
      vector<double> lgkk, lgPk;
      Table_PkCodes (method_Pk, NL, lgkk, lgPk, redshift, output_root, k_max, file_par);

      cosmobl::glob::STR_xi str;
      str.rr = rr;
      str.aa = aa;
      str.lgkk = lgkk;
      str.lgPk = lgPk;

      Func.function = &glob::func_xi_GSL;
      Func.params = &str;
      gsl_integration_qag (&Func, k_min, k_max, 0., prec, limit_size, 5, ww, &Int, &error); 
    }

    else ErrorMsg("Error in cosmobl::Cosmology::xi_DM of PkXi.cpp: method_Pk is wrong!");

    gsl_integration_workspace_free (ww);
  }


  else { // using FFTLOG

    if (method_Pk=="EisensteinHu") 
	ErrorMsg("Error in xi_DM of Cosmology, EisensteinHu method only works with GSL integration");

    else if (method_Pk=="CAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="classgal_v1") {
      vector<double> r,xi;
      Table_XiCodes (method_Pk, NL, r, xi, redshift, output_root, k_max, file_par);
      Int = interpolated(rr, r, xi, "Spline");
    }

    else ErrorMsg("Error in cosmobl::Cosmology::xi_DM of PkXi.cpp: method_Pk is wrong!");

  }


  if (Norm==1) Pk_0(method_Pk, redshift, output_root, k_min, k_max, GSL, prec, file_par); 

  double PP0 = -1.;
  if (method_Pk=="EisensteinHu") PP0 = m_Pk0_EH;
  if (method_Pk=="CAMB") PP0 = m_Pk0_CAMB;
  if (method_Pk=="MPTbreeze-v1") PP0 = m_Pk0_MPTbreeze;
  if (method_Pk=="classgal_v1") PP0 = m_Pk0_CLASS;

  return PP0*fact*Int;
}


// =====================================================================================


double cosmobl::Cosmology::wp_DM (const double rp, const string method_Pk, const double redshift, const string output_root, const bool NL, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const string file_par)
{
  int Norm = norm;
  if (Norm==-1) Norm = (m_sigma8>0) ? 1 : 0;

  if (method_Pk=="MPTbreeze-v1") Norm = 0; // check!!!


  // check if the table with lg(r)-lg(xi) already exists

  string mDir = (GSL==0) ? "Numerical" : "GSL";

  string dir = par::DirCosmo+"Cosmology/Tables/"+mDir+"/"+method_Pk+"/h"+conv(m_hh,par::fDP3)+"_OmB"+conv(m_Omega_baryon,par::fDP3)+"_OmCDM"+conv(m_Omega_CDM,par::fDP3)+"_OmL"+conv(m_Omega_DE,par::fDP3)+"_OmN"+conv(m_Omega_neutrinos,par::fDP3)+"_Z"+conv(redshift,par::fDP3)+"_scalar_amp"+conv(m_scalar_amp,par::ee3)+"_n"+conv(m_n_spec,par::fDP3)+"/";
  
  string file_table = (NL) ? dir+"xiDM_NL.dat" : dir+"xiDM_Lin.dat";
  ifstream fin;
  fin.open(file_table.c_str());
  
  double RR, XI;
  vector<double> rr, Xi;

  if (fin)  // read the table
    while (fin>>RR>>XI) {
      rr.push_back(RR);
      Xi.push_back(XI);
    }

  else { // create the table
    cout <<"I'm writing the file: "<<file_table<<endl;
    
    string MK = "mkdir -p "+dir; if (system (MK.c_str())) {};
    ofstream fout (file_table.c_str()); checkIO (file_table,0); 
    
    int step = 200;
    vector<double> rad = logarithmic_bin_vector(step, r_min, r_max);

    int index = 0;
    while (index<int(rad.size())) {
      RR = rad[index];
      XI = xi_DM(rad[index], method_Pk, redshift, output_root, NL, Norm, k_min, k_max, aa, GSL, prec, file_par);
      fout <<RR<<"   "<<XI<<endl;
      cout <<"xi("<<RR<<") = "<<XI<<endl;
      rr.push_back(RR);
      Xi.push_back(XI);
      index ++;
    }
    fout.clear(); fout.close(); cout <<"I wrote the file: "<<file_table<<endl;
  }
    
  fin.clear(); fin.close();


  return wp(rp, rr, Xi, r_max);
}


// =====================================================================================

                            
double cosmobl::Cosmology::sigmaR_DM (const double RR, const int corrType, const string method_Pk, const double redshift, const string output_root, const bool NL, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const string file_par)
{
  // check if the table with lg(r)-lg(xi) already exists

  string mDir = (GSL==0) ? "Numerical" : "GSL";
  
  string dir = par::DirCosmo+"Cosmology/Tables/"+mDir+"/"+method_Pk+"/h"+conv(m_hh,par::fDP3)+"_OmB"+conv(m_Omega_baryon,par::fDP3)+"_OmCDM"+conv(m_Omega_CDM,par::fDP3)+"_OmL"+conv(m_Omega_DE,par::fDP3)+"_OmN"+conv(m_Omega_neutrinos,par::fDP3)+"_Z"+conv(redshift,par::fDP3)+"_scalar_amp"+conv(m_scalar_amp,par::ee3)+"_n"+conv(m_n_spec,par::fDP3)+"/";
  
  string file_table = (corrType==1) ? dir+"xi_DM.dat" : dir+"wp_DM.dat";
  ifstream fin;
  fin.open(file_table.c_str());

  double RRR, XI;
  vector<double> rr, Xi;

  if (fin) { // read the table
    while (fin>>RRR>>XI) {
      rr.push_back(RRR);
      Xi.push_back(XI);
    }
  }

  else { // create the table
    cout <<"I'm writing the file: "<<file_table<<endl;

    string MK = "mkdir -p "+dir; if (system (MK.c_str())) {};
    ofstream fout (file_table.c_str()); checkIO (file_table,0); 

    int step = 200; 
    vector<double> rad = logarithmic_bin_vector(step, r_min, r_max);

    int index = 0;
    while (index<int(rad.size())) {
      RRR = rad[index];
      XI = (corrType==1) ? xi_DM(rad[index], method_Pk, redshift, output_root, NL, norm, k_min, k_max, aa, GSL, prec, file_par) : wp_DM(rad[index], method_Pk, redshift, output_root, NL, norm, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
      fout <<RRR<<"   "<<XI<<endl;
      cout <<"xi("<<RRR<<") = "<<XI<<endl;
      rr.push_back(RRR);
      Xi.push_back(XI);
      index ++;
    }
    fout.clear(); fout.close(); cout <<"I wrote the file: "<<file_table<<endl;
  }
    
  fin.clear(); fin.close();

  return sigmaR(RR, corrType, rr, Xi);
}


// =====================================================================================

/// @cond glob

double cosmobl::glob::func_SSM_EH_GSL (double kk, void *params) 
{
  struct cosmobl::glob::STR_SSM_EH *pp = (struct cosmobl::glob::STR_SSM_EH *) params;

  Cosmology cosm (pp->Omega_matter, pp->Omega_baryon, pp->Omega_neutrinos, pp->massless_neutrinos, pp->massive_neutrinos, pp->Omega_DE, pp->Omega_radiation, pp->hh, pp->scalar_amp, pp->n_spec, pp->w0, pp->wa, pp->fNL, pp->type_NG, pp->model, pp->unit);

  double RHO = cosm.RhoZero();
  if (pp->unit==0) {cosm.set_unit(1); RHO = cosm.Rho(pp->Omega_matter,pp->Omega_neutrinos); cosm.set_unit(0);}
  double rr = Radius(pp->mass,RHO);

  return cosm.Pk_UnNorm(kk,pp->redshift,pp->method_Pk)*pow(TopHat_WF(kk*rr)*kk,2);  
}

/// @endcond

// =====================================================================================


double cosmobl::Cosmology::sigma8_Pk (const string method_Pk, const double redshift, const string output_root, const bool NL, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par) const 
{  
  if (NL==1) {cout <<"Attention: sigma8 is defined for the linear P(k)!"<<endl;}

  double RR = 8.; // sigma_8 = sigma(8Mpc/h)
  double RHO = Rho(m_Omega_matter,m_Omega_neutrinos); 
  double MM = Mass(RR,RHO);
  double Int = -1.;

  if (GSL) {
    int limit_size = 1000;
    gsl_integration_workspace *ww = gsl_integration_workspace_alloc (limit_size);
    gsl_function Func;

    double error = -1.; 
 
    if (method_Pk=="EisensteinHu") {
      cosmobl::glob::STR_SSM_EH str;
      str.Omega_matter = m_Omega_matter; 
      str.Omega_baryon = m_Omega_baryon; 
      str.Omega_neutrinos = m_Omega_neutrinos; 
      str.massless_neutrinos = m_massless_neutrinos; 
      str.massive_neutrinos = m_massive_neutrinos; 
      str.Omega_DE = m_Omega_DE; 
      str.Omega_radiation = m_Omega_radiation; 
      str.hh = m_hh; 
      str.scalar_amp = m_scalar_amp; 
      str.n_spec = m_n_spec;
      str.w0 = m_w0; 
      str.wa = m_wa; 
      str.fNL = m_fNL;
      str.type_NG = m_type_NG;
      str.model = m_model;
      str.unit = m_unit;
      str.method_Pk = method_Pk;  
      str.redshift = redshift;    
      str.mass = MM;
      str.rho = RHO;

      Func.function = &glob::func_SSM_EH_GSL;
      Func.params = &str;
      gsl_integration_qag (&Func, k_min, k_max, 0., prec, limit_size, 6, ww, &Int, &error); 
    }

    else if (method_Pk=="CAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="classgal_v1") {
      vector<double> lgkk, lgPk;
      Table_PkCodes (method_Pk, NL, lgkk, lgPk, redshift, output_root, k_max, file_par);
   
      cosmobl::glob::STR_SSM str;
      str.unit = m_unit;
      str.hh = m_hh;
      str.n_spec = m_n_spec;
      str.mass = MM;
      str.rho = RHO;
      str.lgkk = lgkk;
      str.lgPk = lgPk;

      Func.function = &glob::func_SSM_GSL;
      Func.params = &str;
      gsl_integration_qag (&Func, k_min, k_max, 0., prec, limit_size, 6, ww, &Int, &error); 
    }
    
    else ErrorMsg("Error in cosmobl::Cosmology::sigma8_Pk of PkXi.cpp: method_Pk is wrong!");

    gsl_integration_workspace_free (ww);
  } 
  
  else { // using Numerical libraries

    if (method_Pk=="EisensteinHu") {
      cosmobl::classfunc::func_SSM func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, method_Pk, MM, redshift);
    
      Midpnt<cosmobl::classfunc::func_SSM> q1(func,0.,1.); 
      Midinf<cosmobl::classfunc::func_SSM> q2(func,1.,1.e30);
      Int = qromo(q1)+qromo(q2);
    } 

    else if (method_Pk=="CAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="classgal_v1") {
      vector<double> lgkk, lgPk;
      Table_PkCodes (method_Pk, NL, lgkk, lgPk, redshift, output_root, k_max, file_par);
    
      cosmobl::classfunc::func_SSM_Table func (m_hh, m_n_spec, RHO, m_unit, lgkk, lgPk, MM); 
      Midpnt<cosmobl::classfunc::func_SSM_Table> q1(func,0.,1.); 
      Midinf<cosmobl::classfunc::func_SSM_Table> q2(func,1.,1.e30);
      Int = qromo(q1)+qromo(q2);
    }

    else ErrorMsg("Error in cosmobl::Cosmology::sigma8_Pk of PkXi.cpp: method_Pk is wrong!");
  }

  return sqrt(1./(2.*par::pi*par::pi)*Int);
}


// =====================================================================================


double cosmobl::Cosmology::Sn_PT (const int nn, const double RR, const string method_SS, const string output_root, const string interpType, const int Num, const double stepsize, const double k_max, const string file_par) const 
{
  if (3>nn || nn>5) { string Err = "Error in cosmobl::Cosmology::Sn_PT of PkXi.cpp: nn = " + conv(nn,par::fINT); ErrorMsg(Err); }
  
  double redshift = 0.; // (the hierarchical moments predicted by the PT do not depend on the redshift)

  double gamma1 = 1., gamma2 = -1., gamma3 = -1., d2S = -1., d3S = -1., Sn = 1.;

  double RHO = Rho(m_Omega_matter,m_Omega_neutrinos); 
  double MASS = Mass(RR,RHO);
  double SSS = SSM_norm(MASS,method_SS,redshift,output_root,k_max,file_par);

  //double FF = 4.*par::pi*RHO*pow(RR,3)/SSS;
  //gamma1 = FF*dnSM(1,MASS,method_SS,redshift,output_root,interpType,Num,stepsize,k_max,file_par); 
  
  gamma1 = RR/SSS*dnSR(1,RR,method_SS,redshift,output_root,interpType,Num,stepsize,k_max,file_par);

  if (nn>3) {
    d2S = dnSR(2,RR,method_SS,redshift,output_root,interpType,Num,stepsize,k_max,file_par);
    gamma2 = gamma1+pow(RR,2)/SSS*d2S;
  }

  if (nn>4) {
    d3S = dnSR(3,RR,method_SS,redshift,output_root,interpType,Num,stepsize,k_max,file_par);
    gamma3 = gamma2+pow(RR,2)/SSS*(2.*d2S+RR*d3S);
  }

  if (nn==3) Sn = 34./7.+gamma1;
  if (nn==4) Sn = 60712./1323.+62./3.*gamma1+7./3.*pow(gamma1,2)+2./3.*gamma2;
  if (nn==5) Sn = 200575880./305613.+1847200./3969.*gamma1+6940./63.*pow(gamma1,2)+235./27.*pow(gamma1,3)+1490./63.*gamma2+50./9.*gamma1*gamma2+10./27.*gamma3;
  
  return Sn;
}


// =====================================================================================


double cosmobl::Cosmology::Sigman_PT (const int nn, const double RR, const string method_SS, const string output_root, const string interpType, const int Num, const double stepsize, const double k_max, const string file_par) const 
{
  if (3>nn || nn>5) { string Err = "Error in cosmobl::Cosmology::Sigma_PT of PkXi.cpp: nn = " + conv(nn,par::fINT); ErrorMsg(Err); }
 
  double redshift = 0.; // (the hierarchical moments predicted by the PT do not depend on the redshift)
  
  double RHO = Rho(m_Omega_matter, m_Omega_neutrinos); 
  double MASS = Mass(RR, RHO);
  double SSS = SSM_norm(MASS, method_SS, redshift, output_root, k_max, file_par);

  double gamma1 = RR/SSS*dnSR(1, RR, method_SS, redshift, output_root, interpType, Num, stepsize, k_max, file_par);

  double Sn = -1.;

  if (nn==3) Sn = 36./7.+3./2.*(gamma1+1.);
  if (nn==4) Sn = 2540./49.+33.*(gamma1+1.)+21./4.*pow(gamma1+1.,2);
  if (nn==5) Sn = 793.+794.*(gamma1+1.)+265.*pow(gamma1+1.,2)+29.4*pow(gamma1+1.,3);
  
  return Sn;
}


// =====================================================================================


double cosmobl::Cosmology::k_star (const string method_Pk, const double redshift, const string output_root, const double k_max, const string file_par) const 
{  
  if (method_Pk=="EisensteinHu") ErrorMsg("Work in progress... (in k_star of PkXi.cpp)");

  vector<double> lgkk, lgPk;
  bool do_nonlinear = 0;
  Table_PkCodes (method_Pk, do_nonlinear, lgkk, lgPk, redshift, output_root, k_max, file_par);

  cosmobl::classfunc::func_kstar func (m_hh, m_unit, lgkk, lgPk);
  
  Midpnt<cosmobl::classfunc::func_kstar> q1(func,0.,1.); 
  Midinf<cosmobl::classfunc::func_kstar> q2(func,1.,1.e30);
  double Int = qromo(q1)+qromo(q2);
  
  return pow(1./(3.*par::pi*par::pi)*Int,-0.5);
}


// =====================================================================================


void cosmobl::Cosmology::get_xi (vector<double> &rr, vector<double> &Xi, const string method_Pk, const double redshift, const string output_root, const bool xiType, const double k_star, const bool xiNL, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const string file_par)
{
  int Norm = norm;
  if (Norm==-1) Norm = (m_sigma8>0) ? 1 : 0;

  if (method_Pk=="MPTbreeze-v1") Norm = 0; // check!!!

  bool XiNL = xiNL;
  if (XiNL && method_Pk=="EisensteinHu") { WarningMsg("The P(k) of EisensteinHu is linear! --> XiNL = 0"); XiNL = 0; }
  

  // ----- compute the real space DM xi(r) ----- 

  // check if the table with lg(r)-lg(xi) already exists

  string mDir = (GSL==0) ? "Numerical" : "GSL";
  string nDir = (xiType==0) ? method_Pk : "CWmodel";

  string dir = par::DirCosmo+"Cosmology/Tables/"+mDir+"/"+nDir+"/h"+conv(m_hh,par::fDP3)+"_OmB"+conv(m_Omega_baryon,par::fDP3)+"_OmCDM"+conv(m_Omega_CDM,par::fDP3)+"_OmL"+conv(m_Omega_DE,par::fDP3)+"_OmN"+conv(m_Omega_neutrinos,par::fDP3)+"_Z"+conv(redshift,par::fDP3)+"_scalar_amp"+conv(m_scalar_amp,par::ee3)+"_n"+conv(m_n_spec,par::fDP3)+"/";
  
  string file_table = (XiNL) ? dir+"xi_DM.dat": dir+"xi_DM_lin.dat";
  
  //cout <<endl<<"file with tabulated values of xi(r): "<<file_table<<endl<<endl;
  
  ifstream fin (file_table.c_str()); 
  
  double RR, XI;

  if (fin) { // read the table
    while (fin >>RR>>XI) {
      rr.push_back(RR);
      Xi.push_back(XI);
    }
    fin.clear(); fin.close();
  }

  else { // create the table

    string MK = "mkdir -p "+dir; if (system (MK.c_str())) {};
    ofstream fout (file_table.c_str()); checkIO (file_table,0); 

    int step = 1000;
    vector<double> rad = linear_bin_vector(step, r_min, r_max);

    int index = 0;
    while (index<int(rad.size())) {
      RR = rad[index];
      XI = (xiType==0) ? xi_DM(rad[index], method_Pk, redshift, output_root, XiNL, Norm, k_min, k_max, aa, GSL, prec, file_par) : 
	xi_star(rad[index], redshift, output_root, k_star, k_max, k_max, GSL, prec, file_par);

      fout <<RR<<"   "<<XI<<endl;
      cout <<"xi("<<RR<<") = "<<XI<<endl;

      rr.push_back(RR);
      Xi.push_back(XI);
      index ++;
    }
    fout.clear(); fout.close(); cout <<"I wrote the file: "<<file_table<<endl;
  }

}


// =====================================================================================


void cosmobl::Cosmology::get_barred_xi (vector<double> rr, vector<double> Xi, vector<double> &Xi_, vector<double> &Xi__, const string method_Pk, const double redshift, const bool xiType, const double k_star, const bool xiNL, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const string file_par) const
{
  int Norm = norm;
  if (Norm==-1) Norm = (m_sigma8>0) ? 1 : 0;

  if (method_Pk=="MPTbreeze-v1") Norm = 0; // check!!!

  bool XiNL = xiNL;
  if (XiNL && method_Pk=="EisensteinHu") { WarningMsg("The P(k) of EisensteinHu is linear! --> XiNL = 0"); XiNL = 0;}
  
  
  // ----- compute the barred functions: xi_(r) and xi__(r) ----- 
  
  // check if the table with lg(r)-lg(xi) already exists

  string mDir = (GSL==0) ? "Numerical" : "GSL";
  string nDir = (xiType==0) ? method_Pk : "CWmodel";

  string dir = par::DirCosmo+"Cosmology/Tables/"+mDir+"/"+nDir+"/h"+conv(m_hh,par::fDP3)+"_OmB"+conv(m_Omega_baryon,par::fDP3)+"_OmCDM"+conv(m_Omega_CDM,par::fDP3)+"_OmL"+conv(m_Omega_DE,par::fDP3)+"_OmN"+conv(m_Omega_neutrinos,par::fDP3)+"_Z"+conv(redshift,par::fDP3)+"_scalar_amp"+conv(m_scalar_amp,par::ee3)+"_n"+conv(m_n_spec,par::fDP3)+"/";
  
  string file_table = (XiNL) ? dir+"xi_DM.dat": dir+"xi_DM_lin.dat";

  ifstream fin(file_table.c_str()); 

  string file_tableb = (XiNL) ? dir+"xibarred_DM.dat": dir+"xibarred_DM_lin.dat";

  ifstream finb (file_tableb.c_str()); 

  double RR, XI_, XI__;

  bool DO = (rr.size()==0) ? 1 : 0;

  if (finb) { // read the table
    while (finb >>RR>>XI_>>XI__) {
      if (DO) rr.push_back(RR);
      Xi_.push_back(XI_);
      Xi__.push_back(XI__);
    }
    finb.clear(); finb.close();
  }

  else { // create the table
    cout <<"I'm writing the file: "<<file_tableb<<endl;
    
    string MK = "mkdir -p "+dir; if (system (MK.c_str())) {};
    ofstream foutb (file_tableb.c_str()); checkIO (file_tableb,0); 

    vector<double> rad;

    if (DO) {
      int step = 1000;
      rad = linear_bin_vector(step, r_min, r_max);
    } 
    else rad = rr;

    int index = 0;
    while (index<int(rad.size())) {
      XI_ = barred_xi_direct(rad[index],rr,Xi,0,-1,-1);
      XI__ = barred_xi__direct(rad[index],rr,Xi,0,-1,-1);

      foutb <<rad[index]<<"   "<<XI_<<"   "<<XI__<<endl;
      cout <<"r = "<<rad[index]<<" --> xi_ = "<<XI_<<", xi__ = "<<XI__<<endl;

      Xi_.push_back(XI_);
      Xi__.push_back(XI__);

      index ++;
    }
    foutb.clear(); foutb.close(); cout <<"I wrote the file: "<<file_tableb<<endl;
  }
}


// =====================================================================================


double cosmobl::Cosmology::Pk_DeWiggle (const double kk, const double redshift, const double sigma_NL, const string output_root, const bool norm, const double k_min, const double k_max, const double aa, const double prec)
{
  bool NL = 0;

  string author1 = "CAMB";
  string author2 = "EisensteinHu";

  double PkCamb = Pk(kk, author1, NL, redshift, output_root, norm, k_min, k_max, 1, prec);
  double PkEH =  Pk(kk, author2, NL, redshift, output_root, norm, k_min, k_max, 1, prec);

  double PkDEW = PkEH*(1+(PkCamb/PkEH-1)*exp(-0.5*pow(kk*sigma_NL, 2)));
  return PkDEW;
}


// =====================================================================================


double cosmobl::Cosmology::xi_DM_DeWiggle (const double rr, const double redshift, const double sigma_NL, const string output_root, const bool norm, const double k_min, const double k_max, const double aa, const double prec)
{
  bool NL = 0;

  string author1 = "CAMB";
  string author2 = "EisensteinHu";

  vector<double> kk,PkCamb,PkM;
  Table_PkCodes (author1, NL, kk, PkCamb, redshift, output_root, k_max);

  for (size_t i = 0; i<kk.size(); i++) {
    kk[i] = pow(10,kk[i]);
    PkCamb[i] = Pk(kk[i], author1, NL, redshift, output_root, norm, k_min, k_max, 1, prec);
    double PkEH = Pk(kk[i], author2, NL, redshift, output_root, norm, k_min, k_max, 1, prec);
    PkM.push_back(PkEH*(1+(PkCamb[i]/PkEH-1)*exp(-0.5*pow(kk[i]*sigma_NL, 2))));
  }

  return xi_from_Pk(rr, kk, PkM, k_min, k_max, aa, 1, prec);
}


// =====================================================================================


vector<vector<double> > cosmobl::Cosmology::get_XiMonopole_covariance(const int nbins, const double rMin, const double rMax, const double nn, const double Volume, const vector<double> kk, const vector<double> Pk0, const int IntegrationMethod)
{
  int nbins_k = kk.size();
  vector<double> r = linear_bin_vector(nbins,rMin,rMax);
  vector<vector<double>> covariance(nbins,vector<double>(nbins,0));
  double dr=r[1]-r[0];

  vector<vector<double>> sigma2 = cosmobl::sigma2_k(nn, Volume,kk,{Pk0},{0});
  vector<vector<double>> jr(nbins,vector<double>(nbins_k,0));

  for(size_t j=0;j<r.size();j++){
     for(size_t i=0;i<kk.size();i++){
        jr[j][i] = jl_distance_average(kk[i], 0, r[j]-dr, r[j]+dr);
     }
  }

  if (IntegrationMethod==0) //Perform trapezoid integration
  { 
     vector<double> integrand(kk.size(),0);

     for (int i=0;i<nbins;i++){
        for (int j=i;j<nbins;j++){
           for(int ii=0;ii<nbins_k;ii++){
              integrand[ii] =kk[ii]*kk[ii]*sigma2[0][ii]*jr[i][ii]*jr[j][ii]; 
           }

           double Int = trapezoid_integration(kk,integrand);

           Int = Int/(2.*par::pi*par::pi);
           covariance[i][j] = Int;
           covariance[j][i] = Int;
        }
     }
  }
  else if (IntegrationMethod==1) //Perform integration with GSL
  {
     cosmobl::glob::STR_covariance_XiMultipoles_integrand params;
     int limit_size = 1000;

     gsl_function Func;
     Func.function = &covariance_XiMultipoles_integrand;
     Func.params = &params;

     double k_min=1.e-4;
     double k_max = 1.e0; 
     double prec = 1.e-2;

     classfunc::func_grid_GSL s2(kk,sigma2[0],"Spline");
     params.s2 = &s2;

     for (int i=0;i<nbins;i++){
        classfunc::func_grid_GSL jl1r1(kk,jr[i],"Spline");
        for (int j=i;j<nbins;j++){
           classfunc::func_grid_GSL jl2r2(kk,jr[j],"Spline");
           params.jl1r1 = &jl1r1;
           params.jl2r2 = &jl2r2;

           double Int = GSL_integrate_qag(Func,k_min,k_max,prec,limit_size,6);
           Int = Int/(2.*par::pi*par::pi);
           covariance[i][j] = Int;
           covariance[j][i] = Int;
           jl2r2.free();
        } 
        jl1r1.free();
     }
     s2.free();
  }  

  return covariance;

}


// =====================================================================================


vector<vector<double> > cosmobl::Cosmology::get_XiMultipoles_covariance (const int nbins, const double rMin, const double rMax, const double nn, const double Volume, const vector<double> kk, const vector<double> Pk0, const vector<double> Pk2, const vector<double> Pk4, const int IntegrationMethod)
{
  int n_leg = 3;
  int nbins_k = kk.size();
  vector<double> r = linear_bin_vector(nbins, rMin, rMax);
  vector<vector<double>> covariance(n_leg*nbins, vector<double>(n_leg*nbins, 0));

  vector<vector<double>> sigma2 = cosmobl::sigma2_k(nn, Volume, kk, {Pk0, Pk2, Pk4}, {0,2,4});
  double dr = r[1]-r[0];

  vector<vector<vector<double> >> jr(n_leg, vector<vector<double>>(nbins, vector<double>(nbins_k, 0)));

  for (int ll=0; ll<n_leg; ll++)
    for (size_t j=0; j<r.size(); j++)
      for (size_t i=0; i<kk.size(); i++)
	jr[ll][j][i] = jl_distance_average(kk[i], ll*2, r[j]-dr, r[j]+dr);

  if (IntegrationMethod==0) { // perform trapezoid integration
  
    vector<double> integrand(kk.size(), 0);

    for (int l1=0; l1<n_leg; l1++) {
      for (int l2=l1; l2<n_leg; l2++) {
	int index = l2+n_leg*l1;
	for (int i=0; i<nbins; i++) {
	  for (int j=i; j<nbins; j++) {
	    for(int ii=0; ii<nbins_k; ii++) {
	      integrand[ii] = kk[ii]*kk[ii]*sigma2[index][ii]*jr[l1][i][ii]*jr[l2][j][ii]; 
	    }

	    double Int = trapezoid_integration(kk, integrand);
 
	    Int = Int/(2.*par::pi*par::pi);
	    covariance[i+nbins*l1][j+nbins*l2] = Int;
	    covariance[j+nbins*l1][i+nbins*l2] = Int;
	    covariance[i+nbins*l2][j+nbins*l1] = Int;
	    covariance[j+nbins*l2][i+nbins*l1] = Int;

	  }
	}
      }
    }
  }
  
  else if (IntegrationMethod==1) // perform integration with GSL
  {
    cosmobl::glob::STR_covariance_XiMultipoles_integrand params;
    int limit_size = 1000;

    gsl_function Func;
    Func.function = &covariance_XiMultipoles_integrand;
    Func.params = &params;

    double k_min = 1.e-4;
    double k_max = 1.e0; 
    double prec = 1.e-2;

    for (int l1=0; l1<n_leg; l1++) {
      for (int l2=l1; l2<n_leg; l2++) {
	int index = l2+n_leg*l1;
	classfunc::func_grid_GSL s2(kk, sigma2[index], "Spline");
	params.s2 = &s2;

	for (int i=0; i<nbins; i++) {
	  classfunc::func_grid_GSL jl1r1(kk, jr[l1][i], "Spline");
	  for (int j=i;j<nbins;j++){
	    classfunc::func_grid_GSL jl2r2(kk, jr[l2][j], "Spline");
	    params.jl1r1 = &jl1r1;
	    params.jl2r2 = &jl2r2;

	    double Int = GSL_integrate_qag(Func, k_min, k_max, prec, limit_size, 6);
	    Int = Int/(2.*par::pi*par::pi);
	    covariance[i+nbins*l1][j+nbins*l2] = Int;
	    covariance[j+nbins*l1][i+nbins*l2] = Int;
	    covariance[i+nbins*l2][j+nbins*l1] = Int;
	    covariance[j+nbins*l2][i+nbins*l1] = Int;
	    //cout << l1 << " " <<l2 << " " << i << " " << j << " " << Int << endl;
	    jl2r2.free();
	  }
	  jl1r1.free();
	}
	s2.free();
      }  
    }
  }

  return covariance;
}


// =====================================================================================


vector<vector<double> > cosmobl::Cosmology::get_XiMultipoles (const int nbins, const double rMin, const double rMax, const vector<double> kk, const vector<double> Pk0, const vector<double> Pk2, const vector<double> Pk4, const int IntegrationMethod)
{
  int nbins_k = kk.size();
  vector<double> r = linear_bin_vector(nbins,rMin,rMax);
  double f0 = 1./(2.*par::pi*par::pi), f2 = -1./(2.*par::pi*par::pi), f4 =1./(2.*par::pi*par::pi); 

  vector<double> xi0(nbins,0), xi2(nbins,0), xi4(nbins,0);

  vector<vector<double> > j0_vec(nbins,vector<double>(nbins_k,0)), j2_vec(nbins,vector<double>(nbins_k,0)), j4_vec(nbins,vector<double>(nbins_k,0));

  for(size_t i=0;i<r.size();i++){
     for(size_t j=0;j<kk.size();j++){
        double xx = kk[j]*r[i];
        j0_vec[i][j] = j0(xx); 
        j2_vec[i][j] = j2(xx); 
        j4_vec[i][j] = j4(xx); 
     }
  }

  if (IntegrationMethod==0) //Perform trapezoid integration
  { 

     for (int i=0;i<nbins;i++){
        
        vector<double> i0(kk.size(),0), i2(kk.size(),0), i4(kk.size(),0);

        for(int j=0;j<nbins_k;j++){
           double xx = kk[j]*r[i];
           i0[j] =kk[j]*kk[j]*Pk0[j]*j0(xx); 
           i2[j] =kk[j]*kk[j]*Pk2[j]*j2(xx); 
           i4[j] =kk[j]*kk[j]*Pk4[j]*j4(xx); 
        }

        xi0[i] = trapezoid_integration(kk,i0)*f0;
        xi2[i] = trapezoid_integration(kk,i2)*f2;
        xi4[i] = trapezoid_integration(kk,i4)*f4;
     }

  }
  else if (IntegrationMethod==1) //Perform integration with GSL
  {
     cosmobl::glob::STR_XiMultipoles_integrand params;
     int limit_size = 1000;

     gsl_function Func;
     Func.function = &XiMultipoles_integrand;
     Func.params = &params;

     double k_min=1.e-4;
     double k_max = 1.e2; 
     double prec = 1.e-3;

     classfunc::func_grid_GSL Pk0_interp(kk,Pk0,"Spline");
     classfunc::func_grid_GSL Pk2_interp(kk,Pk2,"Spline");
     classfunc::func_grid_GSL Pk4_interp(kk,Pk4,"Spline");
     for (int i=0;i<nbins;i++){
        params.r = r[i];
        
        params.Pkl = &Pk0_interp;
        params.l = 0;
        xi0[i] = GSL_integrate_qag(Func,k_min,k_max,prec,limit_size,6)*f0;

        params.Pkl = &Pk2_interp;
        params.l = 2;
        xi2[i] = GSL_integrate_qag(Func,k_min,k_max,prec,limit_size,6)*f2;

        params.Pkl = &Pk4_interp;
        params.l = 4;
        xi4[i] = GSL_integrate_qag(Func,k_min,k_max,prec,limit_size,6)*f4;

     }

     Pk0_interp.free();
     Pk2_interp.free();
     Pk4_interp.free();  

  }

  return {xi0,xi2,xi4};
}

