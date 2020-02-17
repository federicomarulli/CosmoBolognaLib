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

using namespace std;

using namespace cbl;
using namespace cosmology;


// =====================================================================================


double cbl::cosmology::Cosmology::As (const double sigma8) const
{ 
  return pow(sigma8/1.79e4*pow(m_Omega_baryon*m_hh*m_hh/0.024, 1./3.)*pow(m_Omega_matter*m_hh*m_hh/0.14, -0.563)*pow(7.808*m_hh, 0.5*(1.-m_n_spec))*pow(m_hh/0.72, -0.693)*0.76/gg(0.), 2);
}


// =====================================================================================


double cbl::cosmology::Cosmology::sigma8_interpolated (const double redshift) const
{
  double wm = m_Omega_matter*m_hh*m_hh;
  double wb = m_Omega_baryon*m_hh*m_hh;
  double sigma8 = 0.1058*pow(m_scalar_amp/2.196e-9, 0.5)*pow(wm/0.1426, 0.520)*pow(wb/0.02205, -0.294)*pow(m_hh/0.673, 0.683)*pow((m_massless_neutrinos+m_massive_neutrinos)/3.046, -0.24)*exp(0.3727*(m_n_spec-0.96))*pow(1-m_Omega_k, 0.175)*DD(redshift)/DD(9); 

  return ((m_Omega_neutrinos>0) ? 0.995*sigma8 : sigma8);
}


// =====================================================================================


std::string cbl::cosmology::Cosmology::Pk_output_file (const string code, const bool NL, const double redshift, const bool run, const bool store_output, const string output_root, const double k_max, const string file_par)
{
  string dir_loc = fullpath(par::DirLoc);
  string dir_cosmo = fullpath(par::DirCosmo);
   
  string dir_grid;
  if (NL==0) dir_grid = "output_linear/";
  else if (NL==1) dir_grid= "output_nonlinear/";
  else ErrorCBL("", "Pk_output_file", "PkXi.cpp");
  dir_grid += "h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6)+"/";

  string dir = dir_cosmo+"External/"+code+"/";
  string dirC = dir_cosmo+"External/CAMB/";
  if (chdir (dirC.c_str())) {}

  string dir_output = dir+dir_grid+"pk.dat";

  if (run) {
    vector<double> lgkk, lgPk;
    Table_PkCodes(code, NL, lgkk, lgPk, redshift, store_output, output_root, k_max, file_par);   
  }

  return dir_output;
}


// =====================================================================================


void cbl::cosmology::Cosmology::run_CAMB (const bool NL, const double redshift, const std::string output_root, const std::string output_dir, const double k_max, const std::string file_par) const 
{
  string dir = fullpath(par::DirCosmo)+"External/CAMB/";

  string File_par = file_par;

  bool delete_output = (output_dir==par::defaultString) ? true : false;

  string mkdir = "mkdir -p "+output_dir; if (system(mkdir.c_str())) {} 

  string OutputRoot = (output_dir==par::defaultString) ? dir+output_root : output_dir+"/"+output_root;

  OutputRoot = (omp_get_max_threads()>1) ? OutputRoot+"_t"+conv(omp_get_thread_num(), par::fINT) : OutputRoot;

  if (File_par==par::defaultString) {

    // --------------------------------------------------------------------------
    // --------- set the cosmological parameters in the file params.ini ---------
    // --------------------------------------------------------------------------

    string file_par_default = dir+"params_cut.ini";
    File_par = OutputRoot+"params.ini";

    string CP = "cp "+file_par_default+" "+File_par; if (system(CP.c_str())) {}
    ofstream fout(File_par.c_str(), std::ios_base::app | std::ios_base::out);
    double HH0 = m_hh*100.;

    fout << "output_root = " << OutputRoot << endl;
    fout << "do_nonlinear = " << conv(NL, par::fINT) << endl;
    fout << "hubble = " << conv(HH0, par::fDP6) << endl;
    fout << "omega_baryon = " << conv(m_Omega_baryon, par::fDP6) << endl;
    fout << "omega_cdm = " << conv(m_Omega_CDM, par::fDP6) << endl;
    fout << "omega_lambda = " << conv(m_Omega_DE, par::fDP6) << endl;
    fout << "omega_neutrino = " << conv(m_Omega_neutrinos, par::fDP6) << endl;
    fout << "transfer_redshift(1) = " << conv(redshift, par::fDP6) << endl;
    fout << "massless_neutrinos = " << conv(m_massless_neutrinos, par::fDP6) << endl;
    fout << "massive_neutrinos = " << conv(m_massive_neutrinos, par::fINT) << endl;
    fout << "scalar_spectral_index(1) = " << conv(m_n_spec, par::fDP6) << endl;
    fout << "w = " << conv(m_w0, par::fDP6) << endl;
    fout << "wa = " << conv(m_wa, par::fDP6) << endl;
    if (m_scalar_amp>0) {
      fout << "scalar_amp(1) = " << conv(m_scalar_amp, par::ee3) << endl;
      fout << "pivot_scalar = " << conv(m_scalar_pivot, par::fDP6) << endl;
    }
    fout << "transfer_kmax = "+conv(k_max, par::fDP6) << endl;
    fout << "re_optical_depth = "+conv(m_tau, par::fDP6) << endl;
    fout << endl;

    fout.clear(); fout.close();

  }

  // --------------------------------------------------------------------------

  string Camb = dir+"camb "+File_par;
  if (system(Camb.c_str())) {}

  if (delete_output) {
    string RM = (output_dir == "./") ? "rm -rf "+OutputRoot+"*" : "rm -rf "+output_dir; 
    if (system(RM.c_str())) {}
  }
}


// =====================================================================================


void cbl::cosmology::Cosmology::run_CAMB (std::vector<double> &kk, std::vector<double> &Pk, const bool NL, const double redshift, const std::string output_root, const std::string output_dir, const double k_max, const std::string file_par) const 
{
  string dir = fullpath(par::DirCosmo)+"External/CAMB/";

  string File_par = file_par;

  bool delete_output = (output_dir==par::defaultString) ? true : false;

  string mkdir = "mkdir -p "+output_dir; if (system(mkdir.c_str())) {} 

  string OutputRoot = (output_dir==par::defaultString) ? dir+output_root : output_dir+"/"+output_root;


  OutputRoot = (omp_get_max_threads() > 1) ? OutputRoot+"_t"+conv(omp_get_thread_num(), par::fINT) : OutputRoot;


  if (File_par==par::defaultString) {

    // --------------------------------------------------------------------------
    // --------- set the cosmological parameters in the file params.ini ---------
    // --------------------------------------------------------------------------

    string file_par_default = dir+"params_cut.ini";
    File_par = OutputRoot+"params.ini";

    string CP = "cp "+file_par_default+" "+File_par; if (system(CP.c_str())) {}
    ofstream fout(File_par.c_str(), std::ios_base::app | std::ios_base::out);
    double HH0 = m_hh*100.;

    fout << "output_root = " << OutputRoot << endl;
    fout << "do_nonlinear = " << conv(NL, par::fINT) << endl;
    fout << "hubble = " << conv(HH0, par::fDP6) << endl;
    fout << "omega_baryon = " << conv(m_Omega_baryon, par::fDP6) << endl;
    fout << "omega_cdm = " << conv(m_Omega_CDM, par::fDP6) << endl;
    fout << "omega_lambda = " << conv(m_Omega_DE, par::fDP6) << endl;
    fout << "omega_neutrino = " << conv(m_Omega_neutrinos, par::fDP6) << endl;
    fout << "transfer_redshift(1) = " << conv(redshift, par::fDP6) << endl;
    fout << "massless_neutrinos = " << conv(m_massless_neutrinos, par::fDP6) << endl;
    fout << "massive_neutrinos = " << conv(m_massive_neutrinos, par::fINT) << endl;
    fout << "scalar_spectral_index(1) = " << conv(m_n_spec, par::fDP6) << endl;
    fout << "w = " << conv(m_w0, par::fDP6) << endl;
    fout << "wa = " << conv(m_wa, par::fDP6) << endl;
    if (m_scalar_amp>0) {
      fout << "scalar_amp(1) = " << conv(m_scalar_amp, par::ee3) << endl;
      fout << "pivot_scalar = " << conv(m_scalar_pivot, par::fDP6) << endl;
    }
    fout << "transfer_kmax = "+conv(k_max, par::fDP6) << endl;
    fout << "re_optical_depth = "+conv(m_tau, par::fDP6) << endl;
    fout << "feedback_level = -1" << endl;
    fout << endl;

    fout.clear(); fout.close();
  }


  // --------------------------------------------------------------------------
  
  string Camb = dir+"camb "+File_par;
  if (system(Camb.c_str())) {}

  // read the power spectrum
  string file_inCAMB = OutputRoot+"_matterpower.dat";
  ifstream fin(file_inCAMB.c_str()); checkIO(fin, file_inCAMB); 

  kk.erase(kk.begin(), kk.end());
  Pk.erase(Pk.begin(), Pk.end());

  double KK, PK;
  while (fin >> KK >> PK) 
    if (KK>0 && PK>0) {
      kk.push_back(KK);
      Pk.push_back(PK);
    } 
  fin.clear(); fin.close();

  if (delete_output) {
    string RM = (output_dir == "./") ? "rm -rf "+OutputRoot+"*" : "rm -rf "+output_dir; 
    if (system(RM.c_str())) {}
  }
}


// =====================================================================================


void cbl::cosmology::Cosmology::Table_PkCodes (const std::string code, const bool NL, std::vector<double> &lgkk, std::vector<double> &lgPk, const double redshift, const bool store_output, const std::string output_root, const double k_max, const string file_par) const 
{

  if (file_par==par::defaultString) {

    if (code=="CAMB" || code=="MPTbreeze-v1") 
      m_Table_Pk_CAMB_MPTbreeze(code, NL, lgkk, lgPk, redshift, store_output, output_root, k_max);
    
    else if (code=="CLASS")
      m_Table_Pk_CLASS(NL, lgkk, lgPk, redshift, store_output, output_root, k_max);

    else
      ErrorCBL("the choosen code is not allowed!", "Table_PkCodes", "PkXi.cpp");
  }

  else {
    WarningMsgCBL("the input k_max parameters will not be used", "m_Table_Pk_parameterFile", "PkXi.cpp");
    m_Table_Pk_parameterFile(code, file_par, NL, lgkk, lgPk, redshift, output_root);   
  }
  
}


// =====================================================================================


void cbl::cosmology::Cosmology::m_Table_Pk_CAMB_MPTbreeze (const string code, const bool NL, std::vector<double> &lgkk, std::vector<double> &lgPk, const double redshift, const bool store_output, const std::string output_root, const double k_max) const
{
  lgkk.erase(lgkk.begin(), lgkk.end());
  lgPk.erase(lgPk.begin(), lgPk.end());
  
  string dir_grid = (NL) ? "output_nonlinear/" : "output_linear/";

  dir_grid += "h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6);

  const string dirC = fullpath(par::DirCosmo)+"External/CAMB/";
  if (chdir(dirC.c_str())) {}

  const string dir_output = dirC+dir_grid+"/";
  
  const string file_in = dir_output+"Pk.dat"; 
  ifstream fin;
  fin.open(file_in.c_str());

  string file_par;
  
  if (!fin) {
      
    // --------------------------------------------------------------------------
    // --------- set the cosmological parameters in the file params.ini ---------
    // --------------------------------------------------------------------------

    file_par = dirC+"params_"+output_root+"_t"+conv(omp_get_thread_num(), par::fINT)+".ini";
    if (system(("cp params.ini "+file_par).c_str())) {}
    
    string sed;
   
    sed = "sed '/test/s//"+output_root+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/do_nonlinear = 0/s//do_nonlinear = "+conv(NL, par::fINT)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/hubble = 70/s//hubble = "+conv(m_hh*100., par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/omega_baryon = 0.0462/s//omega_baryon = "+conv(m_Omega_baryon, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/omega_cdm = 0.2538/s//omega_cdm = "+conv(m_Omega_CDM, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/omega_lambda = 0.7/s//omega_lambda = "+conv(m_Omega_DE, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/omega_neutrino = 0/s//omega_neutrino = "+conv(m_Omega_neutrinos, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/transfer_redshift(1) = 0/s//transfer_redshift(1) = "+conv(redshift, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/massless_neutrinos = 3.046/s//massless_neutrinos = "+conv(m_massless_neutrinos, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/massive_neutrinos = 0/s//massive_neutrinos = "+conv(m_massive_neutrinos, par::fINT)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/scalar_spectral_index(1) = 0.96/s//scalar_spectral_index(1) = "+conv(m_n_spec, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/w = -1/s//w = "+conv(m_w0, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/wa = 0/s//wa = "+conv(m_wa, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
    if (m_scalar_amp>0) {
      sed = "sed '/scalar_amp(1) = 2.1e-9/s//scalar_amp(1) = "+conv(m_scalar_amp, par::ee3)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
      sed = "sed '/pivot_scalar = 0.05/s//pivot_scalar = "+conv(m_scalar_pivot, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
    }
    sed = "sed '/transfer_kmax = 2/s//transfer_kmax = "+conv(k_max, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/re_optical_depth = 0.09/s//re_optical_depth = "+conv(m_tau, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
      
  
    // ----------------------------
    // --------- run CAMB ---------
    // ----------------------------

    if (system(("mkdir -p "+dir_output).c_str())) {}
    
    if (system(("./camb "+file_par).c_str())) {}
    
    if (system(("mv "+output_root+"_matterpower*dat "+dir_output+"Pk.dat").c_str())) {}
    if (system(("mv "+output_root+"_transfer*dat "+dir_output+"transfer_out.dat").c_str())) {}
  
    if (system(("rm -f "+file_par+" *_params.ini").c_str())) {}


    // ---------------------------------
    // --------- run MPTbreeze ---------
    // ---------------------------------
    
    if (code=="MPTbreeze-v1") {
      WarningMsgCBL("check the consistency of the input redshift with the one set in the parameter file", "m_Table_Pk_parameterFile", "PkXi.cpp");

      if (chdir((fullpath(par::DirCosmo)+"External/MPTbreeze-v1").c_str())) {}
  
      if (system(("./mptbreeze -noverbose -camb ../CAMB/"+file_par+" -fileTF ../CAMB/"+output_root+"_transfer_out.dat -sigma8 "+conv(m_sigma8, par::fDP6)+" -omegam "+conv(m_Omega_matter, par::fDP6)+" -ns "+conv(m_n_spec, par::fDP6)+" -w "+conv(m_w0, par::fDP6)+" -redshift "+conv(redshift, par::fDP6)+" -filePk ../CAMB/"+output_root+"_matterpower.dat").c_str())) {}
      if (system(("mv ../CAMB/"+output_root+"_matterpower*dat ../CAMB/Pk.dat").c_str())) {}
    }
    
  }
  
  fin.clear(); fin.close();


  // ---------------------------------------
  // --------- get the output P(k) ---------
  // ---------------------------------------

  if (chdir(fullpath(par::DirLoc).c_str())) {}
  
  fin.open(file_in.c_str()); checkIO(fin, file_in); 
  string line;

  while (getline(fin, line)) {
    stringstream ss(line);
    vector<double> num;
    double aa;
    while (ss>>aa) num.push_back(aa);
    
    if (num[0]>0 && num[1]>0) {
      lgkk.push_back(log10(num[0]));
      if (code=="MPTbreeze-v1") 
	lgPk.push_back(log10(pow(2.*par::pi, 3)*(num[1]+num[2]+num[3])));
      else
	lgPk.push_back(log10(num[1]));
    }
  }
  
  fin.clear(); fin.close();
  
  if (lgkk.size()==0 || lgPk.size()==0)
    ErrorCBL("lgkk.size()="+conv(lgkk.size(), par::fINT)+", lgPk.size()="+conv(lgPk.size(), par::fINT), "m_Table_Pk_CAMB_MPTbreeze", "PkXi.cpp");
  
  if (!store_output) 
    if (system(("rm -rf "+dir_output).c_str())) {}
  
  (void)code; (void)store_output; (void)output_root; (void)k_max;
}


// =====================================================================================


void cbl::cosmology::Cosmology::m_Table_Pk_CLASS (const bool NL, std::vector<double> &lgkk, std::vector<double> &lgPk, const double redshift, const bool store_output, const std::string output_root, const double k_max) const
{
  lgkk.erase(lgkk.begin(), lgkk.end());
  lgPk.erase(lgPk.begin(), lgPk.end());

  string dir_grid = (NL) ? "output_nonlinear/" : "output_linear/";
  
  dir_grid += "h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6)+"/";
  
  const string dirC = fullpath(par::DirCosmo)+"External/CLASS/";
  if (chdir(dirC.c_str())) {}

  const string dir_output = dirC+dir_grid;
  
  const string file_in = dir_output+"Pk.dat"; 
  ifstream fin;
  fin.open(file_in.c_str());

  string file_par;
  
  if (!fin) {
    
    // --------------------------------------------------------------------------
    // --------- set the cosmological parameters in the file params.ini ---------
    // --------------------------------------------------------------------------

    file_par = dirC+"params_"+output_root+"_t"+conv(omp_get_thread_num(), par::fINT)+".ini";
    if (system(("cp params.ini "+file_par).c_str())) {}
    
    string sed;

    sed = "sed '/output\\/test_/s//"+output_root+"_/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/non linear = /s//non linear = halofit/g' "+file_par+" > temp; mv temp "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/h = 0.7/s//h = "+conv(m_hh, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/Omega_b = 0.05/s//Omega_b = "+conv(m_Omega_baryon, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/Omega_cdm = 0.25/s//Omega_cdm = "+conv(m_Omega_CDM, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/N_ncdm = 0/s//N_ncdm = "+conv(m_massive_neutrinos, par::fINT)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/Omega_ncdm = /s//Omega_ncdm = "+conv(m_Omega_neutrinos, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system(sed.c_str())) {}

    if (m_Omega_neutrinos>0) { // check!!!
      sed = "sed '/m_ncdm = 0.04, 0.04, 0.04/s//#m_ncdm/g' "+file_par+" > temp; mv temp "+file_par+""; 
      if (system(sed.c_str())) {}
    }

    sed = "sed '/Omega_Lambda = 0.7/s//Omega_Lambda = "+conv(m_Omega_DE, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/Omega_k = 0./s//Omega_k = "+conv(m_Omega_k, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/z_pk = 0/s//z_pk = "+conv(redshift, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/N_eff = 3.04/s//N_eff = "+conv(m_massless_neutrinos, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/N_ncdm = 0/s//N_ncdm = "+conv(m_massive_neutrinos, par::fINT)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/n_s = 1./s//n_s = "+conv(m_n_spec, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}

    double w00 = max(-0.999,m_w0); // check!!!
    sed = "sed '/w0_fld = -0.9/s//w0_fld = "+conv(w00, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}

    sed = "sed '/wa_fld = 0./s//wa_fld = "+conv(m_wa, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
    if (m_scalar_amp>0) {
      sed = "sed '/A_s = 2.3e-9/s//A_s = "+conv(m_scalar_amp, par::ee3)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
      sed = "sed '/k_pivot = 0.05/s//k_pivot = "+conv(m_scalar_pivot, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
    }
    sed = "sed '/P_k_max_h\\/Mpc = 1/s//P_k_max_h\\/Mpc = "+conv(k_max, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/tau_reio = 0.09/s//tau_reio = "+conv(m_tau, par::fDP6)+"/g' "+file_par+" > temp; mv temp "+file_par; if (system(sed.c_str())) {}

    
    // -----------------------------
    // --------- run CLASS ---------
    // -----------------------------

    if (system(("mkdir -p "+dir_output).c_str())) {}
 
    if (system(("./class "+file_par).c_str())) {}
    
    if (system(("mv "+output_root+"_pk"+((NL) ? ".dat" : "_nl.dat")+" "+dir_output+"Pk.dat").c_str())) {}   
    if (system(("rm -f "+output_root+"_pk* "+file_par).c_str())) {}
  }

  fin.clear(); fin.close();


  // ---------------------------------------
  // --------- get the output P(k) ---------
  // ---------------------------------------

  if (chdir(fullpath(par::DirLoc).c_str())) {}
  
  fin.open(file_in.c_str()); checkIO(fin, file_in); 
  string line;

  for (short i=0; i<4; ++i)
    getline(fin, line);

  while (getline(fin, line)) {
    stringstream ss(line);
    vector<double> num;
    double aa;
    while (ss>>aa) num.push_back(aa);
    
    if (num[0]>0 && num[1]>0) {
      lgkk.push_back(log10(num[0]));
      lgPk.push_back(log10(num[1]));
    }
  }
  
  fin.clear(); fin.close();
  
  if (lgkk.size()==0 || lgPk.size()==0)
    ErrorCBL("lgkk.size()="+conv(lgkk.size(), par::fINT)+", lgPk.size()="+conv(lgPk.size(), par::fINT), "m_Table_Pk_CAMB_MPTbreeze", "PkXi.cpp");
  
  if (!store_output) 
    if (system(("rm -rf "+dir_output).c_str())) {}
}


// =====================================================================================


void cbl::cosmology::Cosmology::m_Table_Pk_parameterFile (const std::string code, const std::string file_par, const bool NL, std::vector<double> &lgkk, std::vector<double> &lgPk, const double redshift, const std::string output_root) const
{
  WarningMsgCBL("Check the consistency between the output_root parameter given in input and the one set in the parameter file", "m_Table_Pk_parameterFile", "PkXi.cpp");
  
  lgkk.erase(lgkk.begin(), lgkk.end());
  lgPk.erase(lgPk.begin(), lgPk.end());
  
  const string dirB = (code=="CAMB" || code=="MPTbreeze-v1") ? fullpath(par::DirCosmo)+"External/CAMB/"
    : fullpath(par::DirCosmo)+"External/"+code+"/";
  
  const string dir_output = dirB+"output_ParameterFiles/"+file_par+"/";

  if (system(("mkdir -p "+dir_output).c_str())) {}
  
  const string file_out = dir_output+"Pk.dat";

  ifstream fin;
  fin.open(file_out.c_str());

  if (!fin) {
  
    const string File_par = file_par+"_t"+conv(omp_get_thread_num(), par::fINT)+".ini";
    if (system(("cp "+file_par+" "+dirB+File_par).c_str())) {}

    if (chdir(dirB.c_str())) {}
  
  
    // --------------------------------------------
    // --------- run the Boltzmann solver ---------
    // --------------------------------------------

    if (code=="CAMB" || code=="MPTbreeze-v1") {
      if (system(("./camb "+File_par).c_str())) {}
      if (system(("mv "+output_root+"_matterpower*dat "+dir_output+"Pk.dat").c_str())) {}
    
      if (code=="MPTbreeze-v1") {
	WarningMsgCBL("check the consistency of the input redshift with the one set in the parameter file", "m_Table_Pk_parameterFile", "PkXi.cpp");
	if (chdir((fullpath(par::DirCosmo)+"External/"+code).c_str())) {}

	if (system(("./mptbreeze -noverbose -camb ../CAMB/"+File_par+" -fileTF ../CAMB/"+output_root+"_transfer_out.dat -sigma8 "+conv(m_sigma8, par::fDP6)+" -omegam "+conv(m_Omega_matter, par::fDP6)+" -ns "+conv(m_n_spec, par::fDP6)+" -w "+conv(m_w0, par::fDP6)+" -redshift "+conv(redshift, par::fDP6)+" -filePk ../CAMB/"+output_root+"_matterpower.dat").c_str())) {}
	if (system(("mv ../CAMB/"+output_root+"_matterpower*dat "+dir_output+"Pk.dat").c_str())) {}
      }

    }
  
    else if (code=="CLASS") {
      WarningMsgCBL("check the consistency of the input NL value with the one set in the parameter file", "m_Table_Pk_parameterFile", "PkXi.cpp");
      if (system(("./class "+File_par).c_str())) {}
      if (system(("mv "+output_root+"pk"+((NL) ? ".dat" : "_nl.dat")+" "+dir_output+"Pk.dat").c_str())) {}   
    }
  
    else ErrorCBL("the chosen code is not allowed!", "m_Table_parameterFile", "PkXi.cpp");
  
    if (system(("rm -f "+File_par+" *_params.ini").c_str())) {}
    
    if (chdir(fullpath(par::DirLoc).c_str())) {}
    
    fin.clear(); fin.close();

  }
  
  
  // ---------------------------------------
  // --------- get the output P(k) ---------
  // ---------------------------------------

  string line;

  if (code=="CLASS")
    for (short i=0; i<4; ++i)
      getline(fin, line);

  while (getline(fin, line)) {
    stringstream ss(line);
    vector<double> num;
    double aa;
    while (ss>>aa) num.push_back(aa);
    
    if (num[0]>0 && num[1]>0) {
      lgkk.push_back(log10(num[0]));
      if (code=="MPTbreeze-v1") 
	lgPk.push_back(log10(pow(2.*par::pi, 3)*(num[1]+num[2]+num[3])));
      else
	lgPk.push_back(log10(num[1]));
    }
    
  }

  fin.clear(); fin.close();
  
  if (lgkk.size()==0 || lgPk.size()==0)
    ErrorCBL("lgkk.size()="+conv(lgkk.size(), par::fINT)+", lgPk.size()="+conv(lgPk.size(), par::fINT), "Table_PkCodes", "PkXi.cpp");
    
  if (system(("rm -f "+dirB+output_root+"* ").c_str())) {}
  
}


// =====================================================================================


void cbl::cosmology::Cosmology::Table_XiCodes (const std::string code, const bool NL, std::vector<double> &rr, std::vector<double> &xi, const double redshift, const bool store_output, const std::string output_root, const double k_max, const std::string file_par) const
{
  vector<double> lgkk, lgPk;
  Table_PkCodes(code, NL, lgkk, lgPk, redshift, store_output, output_root, k_max, file_par);

  if (lgkk.size()==0 || lgPk.size()==0)
    ErrorCBL("lgkk.size()="+conv(lgkk.size(), par::fINT)+", lgPk.size()="+conv(lgPk.size(), par::fINT), "Table_XiCodes", "PkXi.cpp");
  
  string dir_loc = fullpath(par::DirLoc);
  string dir_cosmo = fullpath(par::DirCosmo);
  
  string dir_grid;
  if (NL==false) dir_grid = "output_linear/";
  else dir_grid= "output_nonlinear/";
  
  dir_grid += "h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6);

  string dir = dir_cosmo+"External/"+code+"/";
  string dirFFT = dir_cosmo+"External/fftlog-f90-master/";
  if (chdir(dirFFT.c_str())) {}
  string dir_output = dir+dir_grid+"/";

  string file_in = (code=="MPTbreeze-v1") ? dir_output+"Pk2.dat" : dir_output+"Pk.dat";
  string file_out = dir_output+"Xi.dat";

  const string mkdir = "mkdir -p "+dir_output; if (system(mkdir.c_str())) {} 
  
  ifstream fin;
  fin.open(file_out.c_str());
  
  if (!fin) {
    vector<double> kk = lgkk, Pk = lgPk;
    for_each( kk.begin(), kk.end(), [] (double &vv) { vv = pow(10., vv); } );
    for_each( Pk.begin(), Pk.end(), [] (double &vv) { vv = pow(10., vv); } );
    
    cbl::wrapper::fftlog::transform_FFTlog(rr, xi, 1, kk, Pk);

    ofstream fout(file_out.c_str()); checkIO(fout, file_out);

    for (size_t i=0; i<rr.size(); ++i) 
      fout << rr[i] << "   " << xi[i] << endl;
    
    fout.clear(); fout.close();
  }
  
  if (chdir(dir_loc.c_str())) {}

  fin.clear(); fin.close();

  rr.erase(rr.begin(), rr.end());
  xi.erase(xi.begin(), xi.end());

  fin.open(file_out.c_str()); 

  double RR, XI;
  while (fin >> RR >> XI) {
    rr.push_back(RR);
    xi.push_back(XI);
  }

  fin.clear(); fin.close();
  
  if (!store_output) {
    string RM = "rm -rf "+dir_grid;
    if (system(RM.c_str())) {}
  }
}


// =====================================================================================


void cbl::cosmology::Cosmology::Pk_0 (const std::string method_Pk, const double redshift, const bool store_output, const std::string output_root, const double k_min, const double k_max, const double prec, const std::string file_par)
{  
  if (m_sigma8<0) ErrorCBL("sigma8<0!", "Pk_0", "PkXi.cpp");

  double RR = 8.; // sigma_8 = sigma(8Mpc/h)
  double RHO = rho_m(0., true); 
  double MM = Mass(RR, RHO);

  bool NL = false;
  double Int = -1.;
  double error = -1.; 

  if (method_Pk=="EisensteinHu") {

    EisensteinHu eh;

    eh.TFmdm_set_cosm(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massive_neutrinos, m_Omega_DE, m_hh, redshift, m_scalar_amp, m_scalar_pivot, m_n_spec); 

    auto func = [&] (const double kk)
      {
	return pow(TopHat_WF(kk*RR)*kk, 2)*eh.Pk(kk); 
      };
    
    Int = wrapper::gsl::GSL_integrate_qag(func, k_min, k_max, prec);
  }

  else if (method_Pk=="CAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="CLASS") {
    vector<double> lgkk, lgPk;
    Table_PkCodes(method_Pk, NL, lgkk, lgPk, redshift, store_output, output_root, k_max, file_par);

    int limit_size = 1000;
    gsl_integration_workspace *ww = gsl_integration_workspace_alloc(limit_size);
    gsl_function Func;

    cbl::glob::STR_SSM str;
    str.unit = true;
    str.hh = m_hh;
    str.n_spec = m_n_spec;
    str.mass = MM;
    str.rho = RHO;
    str.lgkk = lgkk;
    str.lgPk = lgPk;

    Func.function = &glob::func_SSM_GSL;
    Func.params = &str;
    gsl_integration_qag(&Func, k_min, k_max, 0., prec, limit_size, 6, ww, &Int, &error); 
    gsl_integration_workspace_free (ww);
  }

  else ErrorCBL("method_Pk is wrong!", "Pk_0", "PkXi.cpp");


  if (method_Pk=="EisensteinHu") m_Pk0_EH = 2.*pow(par::pi*m_sigma8,2)/Int*pow(DD(redshift)/DD(0),2);
  if (method_Pk=="CAMB") m_Pk0_CAMB = 2.*pow(par::pi*m_sigma8,2)/Int*pow(DD(redshift)/DD(0),2);
  if (method_Pk=="MPTbreeze-v1") m_Pk0_MPTbreeze = 2.*pow(par::pi*m_sigma8,2)/Int*pow(DD(redshift)/DD(0),2);
  if (method_Pk=="CLASS") m_Pk0_CLASS = 2.*pow(par::pi*m_sigma8,2)/Int*pow(DD(redshift)/DD(0),2);

}


// =====================================================================================


double cbl::cosmology::Cosmology::Pk (const double kk, const std::string method_Pk, const bool NL, const double redshift, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec, const std::string file_par, const bool unit1)
{ 
  int Norm = norm;

  double fact1 = (m_unit || unit1) ? 1. : m_hh;
  double fact2 = pow(fact1, 3);

  if (Norm==-1) Norm = (m_sigma8>0) ? 1 : 0;
 
  if (method_Pk=="MPTbreeze-v1") Norm = 0; // check!!!

  if (Norm==1) Pk_0(method_Pk, redshift, store_output, output_root, k_min, k_max, prec, file_par);
  else { m_Pk0_EH = 1.; m_Pk0_CAMB = 1.; m_Pk0_MPTbreeze = 1.; m_Pk0_CLASS = 1.; }

  if (method_Pk=="EisensteinHu") {  // NL is not used!!!
    EisensteinHu eh;

    eh.TFmdm_set_cosm(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massive_neutrinos, m_Omega_DE, m_hh, redshift, m_scalar_amp, m_scalar_pivot, m_n_spec); 

    return m_Pk0_EH*eh.Pk(kk/fact1)/fact2;
  }
  
  if (method_Pk=="CAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="CLASS") {
    double lgk = log10(kk/fact1);
    vector<double> lgkk, lgPk;

    Table_PkCodes(method_Pk, NL, lgkk, lgPk, redshift, store_output, output_root, k_max, file_par);

    double lgPK = interpolated(lgk, lgkk, lgPk, "Linear");
    double PP0 = -1.;
    if (method_Pk=="CAMB") PP0 = m_Pk0_CAMB;
    if (method_Pk=="MPTbreeze-v1") PP0 = m_Pk0_MPTbreeze;
    if (method_Pk=="CLASS") PP0 = m_Pk0_CLASS;

    return PP0*pow(10., lgPK)/fact2;
  }

  else { ErrorCBL("method_Pk is wrong!", "Pk", "PkXi.cpp"); return 0; }
  
}

// =====================================================================================


std::vector<double> cbl::cosmology::Cosmology::Pk (const std::vector<double> kk, const std::string method_Pk, const bool NL, const double redshift, const std::string output_dir, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec, const std::string file_par, const bool unit1) 
{

  double fact1 = (m_unit || unit1) ? 1. : 1./m_hh;
  double fact2 = pow(fact1, 3);

  vector<double> newk = kk;
  if (fact1!=1)
    for(size_t i=0; i<newk.size(); i++)
      newk[i] *= fact1;

  // define the normalization
  int Norm = norm;
  if (Norm==-1) Norm = (m_sigma8>0) ? 1 : 0;

  vector<double> Pk(kk.size()); 

  if (method_Pk=="EisensteinHu") { // NL is not used!!!

    if (Norm==1) Pk_0(method_Pk, redshift, store_output, output_root, k_min, k_max, prec, file_par);
    else { m_Pk0_EH = 1.;}

    EisensteinHu eh;

    eh.TFmdm_set_cosm(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massive_neutrinos, m_Omega_DE, m_hh, redshift, m_scalar_amp, m_scalar_pivot, m_n_spec); 

    for (size_t i=0; i<kk.size(); i++) {
      Pk[i] = m_Pk0_EH*eh.Pk(newk[i])*fact2;
    }
  }

  else if (method_Pk=="CAMB") {

    vector<double> _kk, _pk;

    run_CAMB(_kk, _pk, NL, redshift, output_root, output_dir, Max(newk), file_par);
    glob::FuncGrid interp_Pk(_kk, _pk, "Spline");
    Pk = interp_Pk.eval_func(newk);

    if (Norm==1) {

      double sigma8;

      if (NL==true)
	sigma8 = sigma8_Pk("CAMB", 0., store_output, "test", false, k_min, k_max, prec);
      
      else {
	const double RR = 8.;
	glob::FuncGrid interpPk(newk, Pk, "Spline");
	auto func_sigma = [&] (double _k){
	  return pow(TopHat_WF(_k*RR)*_k, 2)*interpPk(_k); 
	};
	sigma8 = sqrt(1./(2.*pow(par::pi, 2))*wrapper::gsl::GSL_integrate_qag (func_sigma, k_min, k_max, 1.e-5))*DD(0.)/DD(redshift);
      }

      m_Pk0_CAMB = pow(m_sigma8/sigma8,2);
    }
    else { m_Pk0_CAMB = 1.;}

    for (size_t i=0; i<kk.size(); i++)
      Pk[i] *= m_Pk0_CAMB*fact2;
  
  }

  else { ErrorCBL("method_Pk is wrong!", "Pk", "PkXi.cpp");  vector<double> vv; return vv; }

  return Pk;
}

// =====================================================================================


void cbl::cosmology::Cosmology::Pk_Kaiser_multipoles (std::vector<double> &Pk0, std::vector<double> &Pk2, std::vector<double> &Pk4, const std::vector<double> kk, const std::string method_Pk, const bool NL, const double redshift, const double bias, const double sigma_NL, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec, const std::string file_par)
{ 
  vector<double> Pk_arr;
  size_t nbin_k = kk.size();
  double f = linear_growth_rate(redshift, 1.);

  if (sigma_NL==0) {
    for (size_t i=0; i<nbin_k; i++)
      Pk_arr.push_back(Pk(kk[i], method_Pk, NL, redshift, store_output, output_root, norm, k_min, k_max, prec,file_par));
  }
  else {
    for (size_t i=0; i<nbin_k; i++)
      Pk_arr.push_back(Pk_DeWiggle(kk[i], redshift, sigma_NL, store_output, output_root, norm, k_min, k_max, prec));
  }

  vector<vector<double>> Pk_multipoles = cbl::Pkl_Kaiser({0,2,4}, kk, Pk_arr, bias, f);
  
  Pk0.erase(Pk0.begin(), Pk0.end());
  Pk2.erase(Pk2.begin(), Pk2.end());
  Pk4.erase(Pk4.begin(), Pk4.end());

  for (size_t i=0; i<nbin_k; i++) {
    Pk0.push_back(Pk_multipoles[0][i]);
    Pk2.push_back(Pk_multipoles[1][i]);
    Pk4.push_back(Pk_multipoles[2][i]);
  }

}


// =====================================================================================

/// @cond glob

double cbl::glob::func_xi_EH_GSL (double kk, void *params) 
{
  struct cbl::glob::STR_xi_EH *pp = (struct cbl::glob::STR_xi_EH *) params;

  EisensteinHu eh;

  eh.TFmdm_set_cosm(pp->Omega_matter, pp->Omega_baryon, pp->Omega_neutrinos, pp->massive_neutrinos, pp->Omega_DE, pp->hh, pp->redshift, pp->scalar_amp, pp->scalar_pivot, pp->n_spec); 
  
  double Int = eh.Pk(kk)*sin(kk*pp->rr)*kk/pp->rr;
 
  return Int*exp(-kk*kk*pp->aa*pp->aa); // eq. 24 of Anderson et al. 2012  
}

/// @endcond

// =====================================================================================


double cbl::cosmology::Cosmology::xi_DM (const double rr, const std::string method_Pk, const double redshift, const bool store_output, const std::string output_root, const bool NL, const int norm, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const std::string file_par) 
{
  bool gsl = GSL;
  if (gsl==false && method_Pk=="EisensteinHu") {
    //WarningMsgCBL("EisensteinHu method only works with GSL integration", "xi_DM", "PkXi.cpp");
    gsl = true;
  }
  
  int Norm = norm;
  if (Norm==-1) Norm = (m_sigma8>0) ? 1 : 0;

  if (method_Pk=="MPTbreeze-v1") Norm = 0; // check!!!
  if (method_Pk=="MPTbreeze-v1" && NL==0) ErrorCBL("MPTbreeze is non-linear!", "xi_DM", "PkXi.cpp");  

  double Int = -1.; 
  double fact = (gsl) ? 1./(2.*pow(par::pi,2)) : 1.;
 
  if (gsl) {
    int limit_size = 1000;
    gsl_integration_workspace *ww = gsl_integration_workspace_alloc(limit_size);
    gsl_function Func;

    double error = -1.; 
 
    if (method_Pk=="EisensteinHu") {
   
      if (m_sigma8<0) ErrorCBL("sigma8<0!", "xi_DM", "PkXi.cpp");
      if (NL==1) WarningMsgCBL("the correlation function by Eisenstein&Hu is linear (see xi_DM of PkXi.cpp)!", "xi_DM", "PkXi.cpp");
  
      cbl::glob::STR_xi_EH str;
      str.Omega_matter = m_Omega_matter; 
      str.Omega_baryon = m_Omega_baryon; 
      str.Omega_neutrinos = m_Omega_neutrinos; 
      str.massless_neutrinos = m_massless_neutrinos; 
      str.massive_neutrinos = m_massive_neutrinos; 
      str.Omega_DE = m_Omega_DE; 
      str.Omega_radiation = m_Omega_radiation; 
      str.hh = m_hh; 
      str.scalar_amp = m_scalar_amp;
      str.scalar_pivot = m_scalar_pivot;
      str.n_spec = m_n_spec;
      str.w0 = m_w0; 
      str.wa = m_wa; 
      str.fNL = m_fNL;
      str.type_NG = m_type_NG;
      str.tau = m_tau;
      str.model = m_model;
      str.unit = m_unit;
      str.rr = rr;
      str.aa = aa;
      str.redshift = redshift;
      str.method_Pk = method_Pk;

      Func.function = &glob::func_xi_EH_GSL;
      Func.params = &str;
      gsl_integration_qag(&Func, k_min, k_max, 0., prec, limit_size, 6, ww, &Int, &error); 
    }

    else if (method_Pk=="CAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="CLASS") {
      vector<double> lgkk, lgPk;
      Table_PkCodes(method_Pk, NL, lgkk, lgPk, redshift, store_output, output_root, k_max, file_par);
      
      glob::STR_xi str;
      str.rr = rr;
      str.aa = aa;
      str.lgkk = lgkk;
      str.lgPk = lgPk;

      Func.function = &glob::func_xi_GSL;
      Func.params = &str;
      gsl_integration_qag(&Func, k_min, k_max, 0., prec, limit_size, 5, ww, &Int, &error); 
    }

    else ErrorCBL("method_Pk is wrong!", "xi_DM", "PkXi.cpp");

    gsl_integration_workspace_free(ww);
  }


  else { // using FFTLOG
    if (method_Pk=="CAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="CLASS") {
      vector<double> r, xi;
      Table_XiCodes(method_Pk, NL, r, xi, redshift, store_output, output_root, k_max, file_par);
      Int = interpolated(rr, r, xi, "Spline");
    }
    
    else ErrorCBL("method_Pk is wrong!", "xi_DM", "PkXi.cpp");
    
  }
  

  if (Norm==1) Pk_0(method_Pk, redshift, store_output, output_root, k_min, k_max, prec, file_par); 

  double PP0 = -1.;
  if (method_Pk=="EisensteinHu") PP0 = m_Pk0_EH;
  if (method_Pk=="CAMB") PP0 = m_Pk0_CAMB;
  if (method_Pk=="MPTbreeze-v1") PP0 = m_Pk0_MPTbreeze;
  if (method_Pk=="CLASS") PP0 = m_Pk0_CLASS;

  return PP0*fact*Int;
}


// =====================================================================================


double cbl::cosmology::Cosmology::wp_DM (const double rp, const std::string method_Pk, const double redshift, const double pimax, const bool store_output, const std::string output_root, const bool NL, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const std::string file_par)
{
  int Norm = norm;
  if (Norm==-1) Norm = (m_sigma8>0) ? 1 : 0;

  if (method_Pk=="MPTbreeze-v1") Norm = 0; // check!!!
  if (method_Pk=="MPTbreeze-v1" && NL==0) ErrorCBL("MPTbreeze is non-linear!", "wp_DM", "PkXi.cpp");  

  // check if the table with lg(r)-lg(xi) already exists

  string mDir = (GSL==0) ? "fftlog" : "GSL";

  string dir_grid = fullpath(par::DirCosmo)+"Cosmology/Tables/"+mDir+"/"+method_Pk+"/h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6)+"/";
  
  string file_table = (NL) ? dir_grid+"xiDM_NL.dat" : dir_grid+"xiDM_Lin.dat";
  ifstream fin;
  fin.open(file_table.c_str());
  
  double RR, XI;
  vector<double> rr, Xi;

  if (fin) // read the table
    while (fin >> RR >> XI) {
      rr.push_back(RR);
      Xi.push_back(XI);
    }

  else { // create the table
    coutCBL <<"I'm writing the file: "<<file_table<<endl;
    
    string MK = "mkdir -p "+dir_grid; if (system(MK.c_str())) {}
    ofstream fout(file_table.c_str()); checkIO(fout, file_table); 
    
    int step = 200;
    vector<double> rad = logarithmic_bin_vector(step, r_min, r_max);

    int index = 0;
    while (index<int(rad.size())) {
      RR = rad[index];
      XI = xi_DM(rad[index], method_Pk, redshift, store_output, output_root, NL, 0, k_min, k_max, aa, GSL, prec, file_par);
      fout << RR << "   " << XI << endl;
      coutCBL << "xi(" << RR << ") = " << XI << endl;
      rr.push_back(RR);
      Xi.push_back(XI);
      index ++;
    }
    fout.clear(); fout.close(); coutCBL <<"I wrote the file: "<<file_table<<endl;
  }
    
  fin.clear(); fin.close();

  double rmax_integral = sqrt(rp*rp+pimax*pimax);
  double Int = wp(rp, rr, Xi, rmax_integral);
  
  if (Norm==1) Pk_0(method_Pk, redshift, store_output, output_root, k_min, k_max, prec, file_par); 

  double PP0 = -1.;
  if (method_Pk=="EisensteinHu") PP0 = m_Pk0_EH;
  if (method_Pk=="CAMB") PP0 = m_Pk0_CAMB;
  if (method_Pk=="MPTbreeze-v1") PP0 = m_Pk0_MPTbreeze;
  if (method_Pk=="CLASS") PP0 = m_Pk0_CLASS;
  return PP0*Int;
}


// =====================================================================================

                            
double cbl::cosmology::Cosmology::sigmaR_DM (const double RR, const int corrType, const std::string method_Pk, const double redshift, const double pimax, const bool store_output, const std::string output_root, const bool NL, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const std::string file_par)
{
  // check if the table with lg(r)-lg(xi) already exists

  string mDir = "GSL";
  
  string dir_grid = fullpath(par::DirCosmo)+"Cosmology/Tables/"+mDir+"/"+method_Pk+"/h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6)+"/";
  
  string file_table = (corrType==1) ? dir_grid+"xi_DM.dat" : dir_grid+"wp_DM.dat";
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
    coutCBL <<"I'm writing the file: "<<file_table<<endl;

    string MK = "mkdir -p "+dir_grid; if (system(MK.c_str())) {}
    ofstream fout(file_table.c_str()); checkIO(fout, file_table); 

    int step = 200; 
    vector<double> rad = logarithmic_bin_vector(step, r_min, r_max);

    int index = 0;
    while (index<int(rad.size())) {
      RRR = rad[index];
      XI = (corrType==1) ? xi_DM(rad[index], method_Pk, redshift, store_output, output_root, NL, norm, k_min, k_max, aa, GSL, prec, file_par) : wp_DM(rad[index], method_Pk, redshift, pimax, store_output, output_root, NL, norm, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
      fout <<RRR<<"   "<<XI<<endl;
      coutCBL <<"xi("<<RRR<<") = "<<XI<<endl;
      rr.push_back(RRR);
      Xi.push_back(XI);
      index ++;
    }
    fout.clear(); fout.close(); coutCBL <<"I wrote the file: "<<file_table<<endl;
  }
    
  fin.clear(); fin.close();

  return sigmaR(RR, corrType, rr, Xi);
}


// =====================================================================================

/// @cond glob

double cbl::glob::func_sigma2M_EH_GSL (double kk, void *params) 
{
  struct cbl::glob::STR_sigma2M_EH *pp = (struct cbl::glob::STR_sigma2M_EH *) params;

  Cosmology cosm(pp->Omega_matter, pp->Omega_baryon, pp->Omega_neutrinos, pp->massless_neutrinos, pp->massive_neutrinos, pp->Omega_DE, pp->Omega_radiation, pp->hh, pp->scalar_amp, pp->scalar_pivot, pp->n_spec, pp->w0, pp->wa, pp->fNL, pp->type_NG, pp->tau, pp->model, pp->unit);

  double RHO = cosm.rho_m(0., true);
  double rr = Radius(pp->mass,RHO);

  EisensteinHu eh;

  eh.TFmdm_set_cosm(pp->Omega_matter, pp->Omega_baryon, pp->Omega_neutrinos, pp->massive_neutrinos, pp->Omega_DE, pp->hh, pp->redshift, pp->scalar_amp, pp->scalar_pivot, pp->n_spec); 
  
  return eh.Pk(kk)*pow(TopHat_WF(kk*rr)*kk, 2);  
}

/// @endcond

// =====================================================================================


double cbl::cosmology::Cosmology::sigma8_Pk (const std::string method_Pk, const double redshift, const bool store_output, const std::string output_root, const bool NL, const double k_min, const double k_max, const double prec, const std::string file_par) const 
{ 
  if (NL) WarningMsgCBL("sigma8 is defined for the linear P(k)!", "sigma8_Pk", "PkXi.cpp");

  if (m_sigma8>0) return m_sigma8*DD(redshift)/DD(0.);

  else {
  
    const double RR = 8.; // sigma_8 = sigma(8Mpc/h)
    const double RHO = rho_m(0., true); 
    const double MM = Mass(RR, RHO);
    double Int = -1., error = -1.;

    const int limit_size = 1000;
    gsl_integration_workspace *ww = gsl_integration_workspace_alloc(limit_size);
    gsl_function Func;

    if (method_Pk=="EisensteinHu") {

      EisensteinHu eh;

      eh.TFmdm_set_cosm(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massive_neutrinos, m_Omega_DE, m_hh, redshift, m_scalar_amp, m_scalar_pivot, m_n_spec); 

      auto func = [&] (double kk)
	{
	  return pow(TopHat_WF(kk*RR)*kk, 2)*eh.Pk(kk); 
	};
    
      Int = wrapper::gsl::GSL_integrate_qag(func, k_min, k_max, prec, limit_size);
    
    }

    else if (method_Pk=="CAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="CLASS") {
      vector<double> lgkk, lgPk;
      Table_PkCodes(method_Pk, NL, lgkk, lgPk, redshift, store_output, output_root, k_max, file_par);

      cbl::glob::STR_SSM str;
      str.unit = true;
      str.hh = m_hh;
      str.n_spec = m_n_spec;
      str.mass = MM;
      str.rho = RHO;
      str.lgkk = lgkk;
      str.lgPk = lgPk;

      Func.function = &glob::func_SSM_GSL;
      Func.params = &str;
      gsl_integration_qag(&Func, k_min, k_max, 0., prec, limit_size, 6, ww, &Int, &error); 
    }

    else ErrorCBL("method_Pk is wrong!", "sigma8_Pk", "PkXi.cpp");

    gsl_integration_workspace_free(ww);

    return sqrt(1./(2.*par::pi*par::pi)*Int);
  }
}


// =====================================================================================


double cbl::cosmology::Cosmology::Sn_PT (const int nn, const double RR, const std::string method_SS, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file) const 
{
  if (3>nn || nn>5) ErrorCBL("nn = " + conv(nn, par::fINT), "Sn_PT", "PkXi.cpp"); 
  
  double redshift = 0.; // (the hierarchical moments predicted by the PT do not depend on the redshift)

  double gamma1 = 1., gamma2 = -1., gamma3 = -1., d2S = -1., d3S = -1., Sn = 1.;

  double RHO = rho_m(0., true); 
  double MASS = Mass(RR,RHO);
  double SSS = sigma2M(MASS, method_SS, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  gamma1 = RR/SSS*dnsigma2R(1, RR, method_SS, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  if (nn>3) {
    d2S = dnsigma2R(2, RR, method_SS, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file);
    gamma2 = gamma1+pow(RR, 2)/SSS*d2S;
  }

  if (nn>4) {
    d3S = dnsigma2R(3, RR, method_SS, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file);
    gamma3 = gamma2+pow(RR, 2)/SSS*(2.*d2S+RR*d3S);
  }

  if (nn==3) Sn = 34./7.+gamma1;
  if (nn==4) Sn = 60712./1323.+62./3.*gamma1+7./3.*pow(gamma1, 2)+2./3.*gamma2;
  if (nn==5) Sn = 200575880./305613.+1847200./3969.*gamma1+6940./63.*pow(gamma1, 2)+235./27.*pow(gamma1, 3)+1490./63.*gamma2+50./9.*gamma1*gamma2+10./27.*gamma3;
  
  return Sn;
}


// =====================================================================================


double cbl::cosmology::Cosmology::Sigman_PT (const int nn, const double RR, const std::string method_SS, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file) const 
{
  if (3>nn || nn>5) ErrorCBL("nn = " + conv(nn, par::fINT), "Sigman_PT", "PkXi.cpp");
 
  double redshift = 0.; // (the hierarchical moments predicted by the PT do not depend on the redshift)
  
  double RHO = rho_m(redshift, true); 
  double MASS = Mass(RR, RHO);
  double SSS = sigma2M(MASS, method_SS, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  double gamma1 = RR/SSS*dnsigma2R(1, RR, method_SS, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  double Sn = -1.;

  if (nn==3) Sn = 36./7.+3./2.*(gamma1+1.);
  if (nn==4) Sn = 2540./49.+33.*(gamma1+1.)+21./4.*pow(gamma1+1.,2);
  if (nn==5) Sn = 793.+794.*(gamma1+1.)+265.*pow(gamma1+1.,2)+29.4*pow(gamma1+1.,3);
  
  return Sn;
}


// =====================================================================================


double cbl::cosmology::Cosmology::k_star (const std::string method_Pk, const double redshift, const bool store_output, const std::string output_root, const double k_max, const std::string file_par) const 
{  
  if (method_Pk=="EisensteinHu") ErrorCBL("", "k_star", "PkXi.cpp", glob::ExitCode::_workInProgress_);

  vector<double> lgkk, lgPk;
  bool do_nonlinear = 0;
  Table_PkCodes(method_Pk, do_nonlinear, lgkk, lgPk, redshift, store_output, output_root, k_max, file_par);

  cbl::classfunc::func_kstar func (m_hh, m_unit, lgkk, lgPk);

  function<double(double)> ff = bind(&cbl::classfunc::func_kstar::operator(), func, std::placeholders::_1);
  double Int1 = wrapper::gsl::GSL_integrate_qag(ff, 0., 1., 1.e-4);
  double Int2 = wrapper::gsl::GSL_integrate_qag(ff, 1., 1.e30, 1.e-4);
  
  double Int = Int1+Int2;
  
  return pow(1./(3.*par::pi*par::pi)*Int,-0.5);
}


// =====================================================================================


void cbl::cosmology::Cosmology::get_xi (std::vector<double> &rr, std::vector<double> &Xi, const std::string method_Pk, const double redshift, const bool store_output, const std::string output_root, const bool xiType, const double k_star, const bool xiNL, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const std::string file_par)
{
  int Norm = norm;
  if (Norm==-1) Norm = (m_sigma8>0) ? 1 : 0;

  if (method_Pk=="MPTbreeze-v1") Norm = 0; // check!!!

  bool XiNL = xiNL;
  if (XiNL && method_Pk=="EisensteinHu")
    { WarningMsgCBL("The P(k) of EisensteinHu is linear! --> XiNL = 0", "get_xi", "PkXi.cpp"); XiNL = 0; }
  

  // ----- compute the real space DM xi(r) ----- 

  // check if the table with lg(r)-lg(xi) already exists

  string mDir = "GSL";
  string nDir = (xiType==0) ? method_Pk : "CWmodel";

  string dir_cosmo = fullpath(par::DirCosmo);
  
  string dir_grid = dir_cosmo+"Cosmology/Tables/"+mDir+"/"+nDir+"/h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6)+"/";
  
  string file_table = (XiNL) ? dir_grid+"xi_DM.dat": dir_grid+"xi_DM_lin.dat";
  
  //coutCBL <<endl<<"file with tabulated values of xi(r): "<<file_table<<endl<<endl;
  
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

    string MK = "mkdir -p "+dir_grid; if (system(MK.c_str())) {}
    ofstream fout(file_table.c_str()); checkIO(fout, file_table); 

    int step = 1000;
    vector<double> rad = linear_bin_vector(step, r_min, r_max);

    int index = 0;
    while (index<int(rad.size())) {
      RR = rad[index];
      XI = (xiType==0) ? xi_DM(rad[index], method_Pk, redshift, store_output, output_root, XiNL, Norm, k_min, k_max, aa, GSL, prec, file_par) : 
	xi_star(rad[index], redshift, store_output, output_root, k_star, k_max, k_max, prec, file_par);

      fout <<RR<<"   "<<XI<<endl;
      coutCBL <<"xi("<<RR<<") = "<<XI<<endl;

      rr.push_back(RR);
      Xi.push_back(XI);
      index ++;
    }
    fout.clear(); fout.close(); coutCBL <<"I wrote the file: "<<file_table<<endl;
  }

}


// =====================================================================================


void cbl::cosmology::Cosmology::get_barred_xi (std::vector<double> rr, std::vector<double> Xi, std::vector<double> &Xi_, std::vector<double> &Xi__, const std::string method_Pk, const double redshift, const bool xiType, const double k_star, const bool xiNL, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const double prec, const std::string file_par) const
{
  (void)k_star; (void)k_min; (void)k_max; (void)aa; (void)prec; (void)file_par;
  
  int Norm = norm;
  if (Norm==-1) Norm = (m_sigma8>0) ? 1 : 0;

  if (method_Pk=="MPTbreeze-v1") Norm = 0; // check!!!

  bool XiNL = xiNL;
  if (XiNL && method_Pk=="EisensteinHu")
    { WarningMsgCBL("The P(k) of EisensteinHu is linear! --> XiNL = 0", "get_barred_xi", "PkXi.cpp"); XiNL = 0;}
  
  
  // ----- compute the barred functions: xi_(r) and xi__(r) ----- 
  
  // check if the table with lg(r)-lg(xi) already exists

  string mDir = "GSL";
  string nDir = (xiType==0) ? method_Pk : "CWmodel";

  string dir_cosmo = fullpath(par::DirCosmo);

  string dir_grid = dir_cosmo+"Cosmology/Tables/"+mDir+"/"+nDir+"/h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6)+"/";
  
  string file_table = (XiNL) ? dir_grid+"xi_DM.dat": dir_grid+"xi_DM_lin.dat";

  ifstream fin(file_table.c_str()); 

  string file_tableb = (XiNL) ? dir_grid+"xibarred_DM.dat": dir_grid+"xibarred_DM_lin.dat";

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
    coutCBL <<"I'm writing the file: "<<file_tableb<<endl;
    
    string MK = "mkdir -p "+dir_grid; if (system(MK.c_str())) {}
    ofstream foutb(file_tableb.c_str()); checkIO(foutb, file_tableb); 

    vector<double> rad;

    if (DO) {
      int step = 1000;
      rad = linear_bin_vector(step, r_min, r_max);
    } 
    else rad = rr;

    int index = 0;
    while (index<int(rad.size())) {
      XI_ = barred_xi_direct(rad[index], rr, Xi, 0, -1, -1);
      XI__ = barred_xi__direct(rad[index], rr, Xi, 0, -1, -1);

      foutb <<rad[index]<<"   "<<XI_<<"   "<<XI__<<endl;
      coutCBL <<"r = "<<rad[index]<<" --> xi_ = "<<XI_<<", xi__ = "<<XI__<<endl;

      Xi_.push_back(XI_);
      Xi__.push_back(XI__);

      index ++;
    }
    foutb.clear(); foutb.close(); coutCBL <<"I wrote the file: "<<file_tableb<<endl;
  }
}


// =====================================================================================


double cbl::cosmology::Cosmology::Pk_DeWiggle (const double kk, const double redshift, const double sigma_NL, const bool store_output, const std::string output_root, const bool norm, const double k_min, const double k_max, const double aa, const double prec)
{
  (void) aa;
  
  bool NL = 0;

  string author1 = "CAMB";
  string author2 = "EisensteinHu";

  double PkCamb = Pk(kk, author1, NL, redshift, store_output, output_root, norm, k_min, k_max, prec);
  double PkEH =  Pk(kk, author2, NL, redshift, store_output, output_root, norm, k_min, k_max, prec);

  double PkDEW = PkEH*(1+(PkCamb/PkEH-1)*exp(-0.5*pow(kk*sigma_NL, 2)));
  return PkDEW;
}


// =====================================================================================


double cbl::cosmology::Cosmology::xi_DM_DeWiggle (const double rr, const double redshift, const double sigma_NL, const bool store_output, const std::string output_root, const bool norm, const double k_min, const double k_max, const double aa, const double prec)
{
  bool NL = false;

  string author1 = "CAMB";
  string author2 = "EisensteinHu";

  vector<double> kk, PkCamb, PkM;
  Table_PkCodes(author1, NL, kk, PkCamb, redshift, store_output, output_root, k_max);

  for (size_t i = 0; i<kk.size(); i++) {
    kk[i] = pow(10,kk[i]);
    PkCamb[i] = Pk(kk[i], author1, NL, redshift, store_output, output_root, norm, k_min, k_max, prec);
    double PkEH = Pk(kk[i], author2, NL, redshift, store_output, output_root, norm, k_min, k_max, prec);
    PkM.push_back(PkEH*(1+(PkCamb[i]/PkEH-1)*exp(-0.5*pow(kk[i]*sigma_NL, 2))));
  }

  return xi_from_Pk(rr, kk, PkM, k_min, k_max, aa, prec);
}


// =====================================================================================


std::vector<std::vector<double> > cbl::cosmology::Cosmology::XiMonopole_covariance (const int nbins, const double rMin, const double rMax, const double nn, const double Volume, const std::vector<double> kk, const std::vector<double> Pk0, const int IntegrationMethod)
{
  int nbins_k = kk.size();
  vector<double> r = linear_bin_vector(nbins,rMin,rMax);
  vector<vector<double>> covariance(nbins,vector<double>(nbins,0));
  double dr=r[1]-r[0];

  vector<vector<double>> sigma2 = cbl::sigma2_k(nn, Volume,kk,{Pk0},{0});
  vector<vector<double>> jr(nbins,vector<double>(nbins_k,0));

  for (size_t j=0; j<r.size(); j++) {
    for (size_t i=0; i<kk.size(); i++) {
      jr[j][i] = jl_distance_average(kk[i], 0, r[j]-dr, r[j]+dr);
    }
  }

  if (IntegrationMethod==0) //Perform trapezoid integration
    { 
      vector<double> integrand(kk.size(),0);

      for (int i=0; i<nbins; i++) {
        for (int j=i; j<nbins; j++) {
	  for (int ii=0; ii<nbins_k; ii++) {
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
      cbl::glob::STR_covariance_XiMultipoles_integrand params;
      int limit_size = 1000;

      gsl_function Func;
      Func.function = &covariance_XiMultipoles_integrand;
      Func.params = &params;

      double k_min=1.e-4;
      double k_max = 1.e0; 
      double prec = 1.e-2;

      glob::FuncGrid s2(kk,sigma2[0], "Spline");
      params.s2 = &s2;

      for (int i=0; i<nbins; i++) {
	glob::FuncGrid jl1r1(kk,jr[i], "Spline");
        for (int j=i; j<nbins; j++) {
	  glob::FuncGrid jl2r2(kk,jr[j], "Spline");
	  params.jl1r1 = &jl1r1;
	  params.jl2r2 = &jl2r2;

	  double Int = wrapper::gsl::GSL_integrate_qag(Func,k_min,k_max, prec,limit_size,6);
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


std::vector<std::vector<double> > cbl::cosmology::Cosmology::XiMultipoles_covariance (const int nbins, const double rMin, const double rMax, const double nn, const double Volume, const std::vector<double> kk, const std::vector<double> Pk0, const std::vector<double> Pk2, const std::vector<double> Pk4, const int IntegrationMethod)
{
  int n_leg = 3;
  int nbins_k = kk.size();
  vector<double> r = linear_bin_vector(nbins, rMin, rMax);
  vector<vector<double>> covariance(n_leg*nbins, vector<double>(n_leg*nbins, 0));

  vector<vector<double>> sigma2 = cbl::sigma2_k(nn, Volume, kk, {Pk0, Pk2, Pk4}, {0,2,4});
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
	    for (int ii=0; ii<nbins_k; ii++) {
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
      cbl::glob::STR_covariance_XiMultipoles_integrand params;
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
	  glob::FuncGrid s2(kk, sigma2[index], "Spline");
	  params.s2 = &s2;

	  for (int i=0; i<nbins; i++) {
	    glob::FuncGrid jl1r1(kk, jr[l1][i], "Spline");
	    for (int j=i; j<nbins; j++) {
	      glob::FuncGrid jl2r2(kk, jr[l2][j], "Spline");
	      params.jl1r1 = &jl1r1;
	      params.jl2r2 = &jl2r2;

	      double Int = wrapper::gsl::GSL_integrate_qag(Func, k_min, k_max, prec, limit_size, 6);
	      Int = Int/(2.*par::pi*par::pi);
	      covariance[i+nbins*l1][j+nbins*l2] = Int;
	      covariance[j+nbins*l1][i+nbins*l2] = Int;
	      covariance[i+nbins*l2][j+nbins*l1] = Int;
	      covariance[j+nbins*l2][i+nbins*l1] = Int;
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


std::vector<std::vector<double> > cbl::cosmology::Cosmology::XiMultipoles (const int nbins, const double rMin, const double rMax, const std::vector<double> kk, const std::vector<double> Pk0, const std::vector<double> Pk2, const std::vector<double> Pk4, const int IntegrationMethod)
{
  int nbins_k = kk.size();
  vector<double> r = linear_bin_vector(nbins,rMin,rMax);
  double f0 = 1./(2.*par::pi*par::pi), f2 = -1./(2.*par::pi*par::pi), f4 =1./(2.*par::pi*par::pi); 

  vector<double> xi0(nbins,0), xi2(nbins,0), xi4(nbins,0);

  vector<vector<double> > j0_vec(nbins,vector<double>(nbins_k,0)), j2_vec(nbins,vector<double>(nbins_k,0)), j4_vec(nbins,vector<double>(nbins_k,0));

  for (size_t i=0; i<r.size(); i++) {
    for (size_t j=0; j<kk.size(); j++) {
      double xx = kk[j]*r[i];
      j0_vec[i][j] = j0(xx); 
      j2_vec[i][j] = j2(xx); 
      j4_vec[i][j] = j4(xx); 
    }
  }

  if (IntegrationMethod==0) { // perform trapezoid integration
    
    for (int i=0; i<nbins; i++) {
      
      vector<double> i0(kk.size(),0), i2(kk.size(),0), i4(kk.size(),0);
      
      for (int j=0; j<nbins_k; j++) {
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
  
  else if (IntegrationMethod==1) { // perform integration with GSL
    
    cbl::glob::STR_XiMultipoles_integrand params;
    int limit_size = 1000;
    
    gsl_function Func;
    Func.function = &XiMultipoles_integrand;
    Func.params = &params;
    
    double k_min = 1.e-4;
    double k_max = 1.e2; 
    double prec = 1.e-3;
    
    glob::FuncGrid Pk0_interp(kk, Pk0, "Spline");
    glob::FuncGrid Pk2_interp(kk, Pk2, "Spline");
    glob::FuncGrid Pk4_interp(kk, Pk4, "Spline");
    
    for (int i=0; i<nbins; i++) {
      params.r = r[i];
      
      params.Pkl = &Pk0_interp;
      params.l = 0;
      params.k_cut = 0.7;
      params.cut_pow = 2.;

      xi0[i] = wrapper::gsl::GSL_integrate_qag(Func,k_min,k_max, prec,limit_size,6)*f0;
      
      params.Pkl = &Pk2_interp;
      params.l = 2;
      params.k_cut = 0.58;
      params.cut_pow = 4.;
      xi2[i] = wrapper::gsl::GSL_integrate_qag(Func, k_min, k_max, prec, limit_size, 6)*f2;
      
      params.Pkl = &Pk4_interp;
      params.l = 4;
      params.k_cut = 0.6;
      params.cut_pow = 2.;
      xi4[i] = wrapper::gsl::GSL_integrate_qag(Func, k_min, k_max, prec, limit_size, 6)*f4;
    }
    
    Pk0_interp.free();
    Pk2_interp.free();
    Pk4_interp.free();  
    
  }

  return {xi0, xi2, xi4};
}


// =====================================================================================


double cbl::cosmology::Cosmology::wtheta_DM (const double theta, const std::vector<double> zz, const std::vector<double> phiz, const std::string interpolationType, const CoordinateUnits coordUnits, const bool GSL, const std::string method_Pk, const bool NL, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec, const std::string file_par)
{
  if (NL)
    ErrorCBL("non linearities in angular correlation function not yet implemented!", "wtheta_DM", "PkXi.cpp", glob::ExitCode::_workInProgress_);

  double theta_rad = cbl::converted_angle (theta, coordUnits, cbl::CoordinateUnits::_radians_);

  vector<double> kk = cbl::logarithmic_bin_vector(200, k_min, k_max);
  vector<double> Pk;

  vector<double> dc;
  vector<double> _zz = cbl::linear_bin_vector(1000, 0., 2*Max(zz));

  for(size_t i=0; i<_zz.size(); i++)
    dc.push_back(this->D_C(_zz[i]));
  cbl::glob::FuncGrid DC_interp(_zz, dc, interpolationType);

  cbl::glob::Distribution phi(cbl::glob::DistributionType::_Interpolated_, zz, phiz, 0, interpolationType);

  double zmin, zmax, zmean;
  zmin = cbl::Min(zz);
  zmax = cbl::Max(zz);
  zmean = phi.mean();

  for (size_t i=0; i<kk.size(); i++)
    Pk.push_back( this->Pk(kk[i], method_Pk, NL, zmean, store_output, output_root, norm, k_min, k_max, prec, file_par));

  vector<double> r, xi;
  cbl::wrapper::fftlog::transform_FFTlog(r, xi, 1, kk, Pk);
  cbl::glob::FuncGrid xi_interp(r, xi, interpolationType);

  if (GSL) {
    auto integrand = [&] (double z1)
      {
	double r1 = DC_interp(z1);

	auto integrand_z2 = [&] (double z2) {
	  double r2 = DC_interp(z2);
	  double ss = sqrt(pow(r1,2)+pow(r2,2)-2*r1*r2*cos(theta_rad));
	  return xi_interp(ss)*phi(z2);
	};

	return cbl::wrapper::gsl::GSL_integrate_qag(integrand_z2, zmin, zmax)*phi(z1);
      };
    return cbl::wrapper::gsl::GSL_integrate_qag(integrand, zmin, zmax);
  }
  else{
    auto integrand = [&] (vector<double> zz)
      {
	double r1 = DC_interp(zz[0]);
	double r2 = DC_interp(zz[1]);
	double ss = sqrt(pow(r1,2)+pow(r2,2)-2*r1*r2*cos(theta_rad));
	return xi_interp(ss)*phi(zz[0])*phi(zz[1]);
      };

    cbl::wrapper::cuba::CUBAwrapper integrator(integrand, 2);

    return integrator.IntegrateCuhre( {{zmin, zmax}, {zmin, zmax}});
  }
}


// =====================================================================================


double cbl::cosmology::Cosmology::wtheta_DM (const double theta, const std::vector<double> kk, const std::vector<double> Pk, const std::vector<double> zz, const std::vector<double> nz, const std::vector<double> phiz, const std::string interpolationType, const CoordinateUnits coordUnits, const bool GSL, const double redshift_Pk)
{
  const double theta_rad = cbl::converted_angle (theta, coordUnits, cbl::CoordinateUnits::_radians_);

  const double zmin = cbl::Min(zz);
  const double zmax = cbl::Max(zz);


  // set the distribution function

  vector<double> nphi(zz.size(), 0);

  for(size_t i=0; i<nphi.size(); i++)
    nphi[i] = (phiz.size()!=0) ? phiz[i]*nz[i] : nz[i];

  cbl::glob::Distribution phi(cbl::glob::DistributionType::_Interpolated_, zz, nphi, 0, interpolationType);
  //const double zmean = phi.mean();

  
  // set the distirbution integrand and normalization

  auto normalization_integrand = [&] (const double redshift)
    {
      return phi(redshift)*this->dV_dZdOmega(redshift, 1);
    };
  
  double normalization = wrapper::gsl::GSL_integrate_qag(normalization_integrand, zmin, zmax);

  
  // compute the 2PCF model

  vector<double> r, xi;
  cbl::wrapper::fftlog::transform_FFTlog(r, xi, 1, kk, Pk);
  cbl::glob::FuncGrid xi_interp(r, xi, interpolationType);
  double integral;

  if (GSL) {
    auto integrand = [&] (double z1)
      {
	double r1 = this->D_C(z1);

	auto integrand_z2 = [&] (double z2) {
	  double DD = this->DD((zz[0]+zz[1])*0.5)/this->DD(redshift_Pk);
	  double r2 = this->D_C(z2);
	  double ss = sqrt(pow(r1,2)+pow(r2,2)-2*r1*r2*cos(theta_rad));
	  return DD*DD*xi_interp(ss)*normalization_integrand(z2);
	};

	return cbl::wrapper::gsl::GSL_integrate_qag(integrand_z2, zmin, zmax)*normalization_integrand(z1);
      };

    integral = cbl::wrapper::gsl::GSL_integrate_qag(integrand, zmin, zmax);

  }
  
  else {
    auto integrand = [&] (vector<double> zz)
      {
	double DD = this->DD((zz[0]+zz[1])*0.5)/this->DD(redshift_Pk);
	double r1 = this->D_C(zz[0]);
	double r2 = this->D_C(zz[1]);
	double ss = sqrt(pow(r1,2)+pow(r2,2)-2*r1*r2*cos(theta_rad));
	return DD*DD*xi_interp(ss)*normalization_integrand(zz[0])*normalization_integrand(zz[1]);
      };

    cbl::wrapper::cuba::CUBAwrapper integrator(integrand, 2);

    integral = integrator.IntegrateCuhre( {{zmin, zmax}, {zmin, zmax}});
  }

  return integral/normalization;
}


// =====================================================================================


std::vector<double> cbl::cosmology::Cosmology::C_l_DM (const int lmax, const std::vector<double> zz, const std::vector<double> phiz, const std::string interpolationMethod, const std::string method_Pk, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec, const std::string file_par)
{
  const double zmin = Min(zz);
  const double zmax = Max(zz);

  vector<double> kk = cbl::logarithmic_bin_vector(200, k_min, k_max);

  vector<double> Pk;
  for(size_t i=0; i<kk.size(); i++)
    Pk.push_back( this->Pk(kk[i], method_Pk, false, 0., store_output, output_root, norm, k_min, k_max, prec, file_par));
  cbl::glob::FuncGrid Pk_interp(kk, Pk, interpolationMethod);

  auto integrand_sbao = [&] (const double kk) 
  {
    return Pk_interp(kk);
  };

  double sbao = sqrt(4.*par::pi*cbl::wrapper::gsl::GSL_integrate_qag(integrand_sbao, 1.e-4, 1)/3./pow(2*par::pi,3));

  cbl::glob::FuncGrid phi(zz, phiz, interpolationMethod);

  vector<double> C_l;

  for (int l=0; l<lmax+1; l++) {
    double integral;
    
    if (l<60) {
      auto integrand = [&] ( const double kk)
      {
	auto integrand_z = [&] (const double zz)
	{
	  return DD(zz)*jl(kk*D_C(zz), l);
	};

	double integral_z = cbl::wrapper::gsl::GSL_integrate_qag(integrand_z, zmin, zmax);
	return kk*kk*Pk_interp(kk)*pow(integral_z, 2)*exp(-kk*kk*sbao*sbao);
      };

      integral = 2*cbl::wrapper::gsl::GSL_integrate_qag(integrand, 1.e-4, 10)/par::pi;
    }
    else {
    
      auto integrand = [&] (const double zz)
      {
	double dc = D_C(zz);
	double kk = (l+0.5)/dc;
	return pow(DD(zz)*phi(zz), 2)*Pk_interp(kk)*HH(zz)/(par::cc*dc*dc);
      };

      integral = cbl::wrapper::gsl::GSL_integrate_qag(integrand, zmin, zmax);
    }   
    C_l.push_back(integral);
  }

  return C_l;
}

/*
std::vector<double> cbl::cosmology::Cosmology::C_l_DM (const int lmax, const std::vector<double> zz, const std::vector<double> phiz, const std::string interpolationMethod, const std::string method_Pk, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec, const std::string file_par)
{
  const double zmin = Min(zz);
  const double zmax = Max(zz);

  vector<double> kk = cbl::logarithmic_bin_vector(200, k_min, k_max);
  vector<double> Pk;
  for(size_t i=0; i<kk.size(); i++)
    Pk.push_back( this->Pk(kk[i], method_Pk, false, 0., output_root, norm, k_min, k_max, prec, file_par));

  cbl::glob::FuncGrid Pk_interp(kk, Pk, interpolationMethod);
  cbl::glob::FuncGrid phi(zz, phiz, interpolationMethod);

  vector<double> C_l;

  for (int l=0; l<lmax+1; l++) {
    auto integrand = [&] ( const double redshift)
    {
      double dc = D_C(redshift);
      double _kk = double(l)/dc;
      return pow(phi(redshift), 2)*HH(redshift)/pow(dc, 2)*pow(DD(redshift),2)*Pk_interp(_kk);
    };
    C_l.push_back(cbl::wrapper::gsl::GSL_integrate_qag(integrand, zmin, zmax)/par::cc);
  }

  return C_l;
 }
*/



// =====================================================================================


std::vector<std::vector<double>> cbl::cosmology::Cosmology::Pk_TNS_AB_multipoles (std::vector<double> kk, const std::string method, const double redshift, const std::string output_dir, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec)
{
  string dir = fullpath(par::DirCosmo)+"External/CPT_Library/";
  string output_tmpCPT = dir+"tmpCPT/";
  string MKout = "mkdir -p " + output_tmpCPT; if (system(MKout.c_str())) {}
  double sigma8_z0 = sigma8_Pk(method, 0., store_output, output_root);

  // Pklin_z0
  const vector<double> Pklin = Pk(kk, method, false, 0., output_dir, store_output, output_root, norm, k_min, k_max, prec);
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
  vector<vector<double>> Pk_AB = {pkA0_new, pkA2_new, pkA4_new, pkB0_new, pkB2_new, pkB4_new};

  string RM = "rm -rf " + output_tmpCPT; if (system (RM.c_str())) {}

  return Pk_AB;

}


//====================================================================
std::vector<std::vector<double>> cbl::cosmology::Cosmology::Pk_TNS_AB_1loop (std::vector<double> kk, const double mu, const std::string method, const double redshift, const std::string output_dir, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec)
{
  vector<vector<double>> Pk_AB_multipoles = Pk_TNS_AB_multipoles(kk, method, redshift, output_dir, store_output, output_root, norm, k_min, k_max, prec);
  vector<double> Pk_A, Pk_B;
  for(size_t ii=0; ii < Pk_AB_multipoles[0].size(); ++ii){
    Pk_A.emplace_back(Pk_AB_multipoles[0][ii] + Pk_AB_multipoles[1][ii]*cbl::legendre_polynomial(mu, 2) + Pk_AB_multipoles[2][ii]*cbl::legendre_polynomial(mu, 4));
    Pk_B.emplace_back(Pk_AB_multipoles[3][ii] + Pk_AB_multipoles[4][ii]*cbl::legendre_polynomial(mu, 2) + Pk_AB_multipoles[5][ii]*cbl::legendre_polynomial(mu, 4));
  }

  vector<vector<double>> Pk_AB = {Pk_A, Pk_B};

  return Pk_AB;
}


//====================================================================


std::vector<std::vector<double>> cbl::cosmology::Cosmology::Pk_TNS_AB_terms_1loop (std::vector<double> kk, const std::string method, const double redshift, const std::string output_dir, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec)
{
  string dir = fullpath(par::DirCosmo)+"External/CPT_Library/";
  string output_tmpCPT = dir+"tmpCPT/";
  string MKout = "mkdir -p " + output_tmpCPT; if (system(MKout.c_str())) {}
  double sigma8_z0 = sigma8_Pk(method, 0., store_output, output_root);

  // input Pklin_z0
  const vector<double> Pklin = Pk(kk, method, false, 0., output_dir, store_output, output_root, norm, k_min, k_max, prec);
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
    auto func_sigma = [&] (double _k){
      return pow(TopHat_WF(_k*RR)*_k, 2)*interpPk(_k);
    };

    sigma8 = sqrt(1./(2.*pow(par::pi, 2))*wrapper::gsl::GSL_integrate_qag (func_sigma, k_min, k_max, 1.e-5))*DD(0.)/DD(redshift);
    m_Pk0_CAMB = pow(m_sigma8/sigma8,2);
    }
    else { m_Pk0_CAMB = 1.;}

   for (size_t i=0; i<kk.size(); i++){
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

  vector<vector<double>> Pk_AB = {pk_A11_new, pk_A12_new, pk_A22_new, pk_A23_new, pk_A33_new, pk_B12_new, pk_B13_new, pk_B14_new, pk_B22_new, pk_B23_new, pk_B24_new, pk_B33_new, pk_B34_new, pk_B44_new};

  string RM = "rm -rf " + output_tmpCPT; if (system (RM.c_str())) {}

  return Pk_AB;
}


//====================================================================


std::vector<std::vector<double>> cbl::cosmology::Cosmology::Pk_TNS_AB_1loop (std::vector<double> kk, const double mu, const double linear_growth_rate, const double bias, const std::string method, const double redshift, const std::string output_dir, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec)
{
  double beta = linear_growth_rate/bias;

  vector<vector<double>> Pk_AB = Pk_TNS_AB_terms_1loop(kk, method, redshift, output_dir, store_output, output_root, norm, k_min, k_max, prec);
  vector<double> Pk_A, Pk_B;

  for(size_t ii=0; ii < kk.size(); ++ii){
    Pk_A.emplace_back(beta*pow(mu,2)*Pk_AB[0][ii] + pow(beta,2)*(pow(mu,2)*Pk_AB[1][ii] + pow(mu,4)*Pk_AB[2][ii]) + pow(beta,3)*(pow(mu,4)*Pk_AB[3][ii] + pow(mu,6)*Pk_AB[4][ii]));
    Pk_B.emplace_back(pow(mu,2)*(pow(beta,2)*Pk_AB[5][ii] + pow(beta,3)*Pk_AB[6][ii] + pow(beta,4)*Pk_AB[7][ii]) + pow(mu,4)*(pow(beta,2)*Pk_AB[8][ii] + pow(beta,3)*Pk_AB[9][ii] + pow(beta,4)*Pk_AB[10][ii]) + pow(mu,6)*(pow(beta,3)*Pk_AB[11][ii] + pow(beta,4)*Pk_AB[12][ii]) + pow(mu,8)*pow(beta,4)*Pk_AB[13][ii]);
  }

  vector<vector<double>> Pk_AB_TNS = {Pk_A, Pk_B};

  return Pk_AB_TNS;
}

//====================================================================


std::vector<std::vector<double>> cbl::cosmology::Cosmology::Pk_TNS_dd_dt_tt (std::vector<double> kk, const std::string method, const double redshift, const std::string output_dir, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec)
{
  string dir = fullpath(par::DirCosmo)+"External/CPT_Library/";
  string output_tmpCPT = dir+"tmpCPT/";
  string MKout = "mkdir -p " + output_tmpCPT; if (system(MKout.c_str())) {}
  double sigma8_z0 = sigma8_Pk(method, 0., store_output, output_root);

  // input Pklin_z0
  const vector<double> Pklin = Pk(kk, method, false, 0., output_dir, store_output, output_root, norm, k_min, k_max, prec);
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
    auto func_sigma = [&] (double _k){
      return pow(TopHat_WF(_k*RR)*_k, 2)*interpPk(_k);
    };

    sigma8 = sqrt(1./(2.*pow(par::pi, 2))*wrapper::gsl::GSL_integrate_qag (func_sigma, k_min, k_max, 1.e-5))*DD(0.)/DD(redshift);
    m_Pk0_CAMB = pow(m_sigma8/sigma8,2);
    }
    else { m_Pk0_CAMB = 1.;}

    for (size_t i=0; i<kk.size(); i++){
   Pkdd_new[i] *= m_Pk0_CAMB;
   Pkdt_new[i] *= m_Pk0_CAMB;
   Pktt_new[i] *= m_Pk0_CAMB;
 }

  vector<vector<double>> Pk_TNS_terms = {Pkdd_new, Pkdt_new, Pktt_new};

  string RM = "rm -rf " + output_tmpCPT; if (system (RM.c_str())) {}

  return Pk_TNS_terms;
}

