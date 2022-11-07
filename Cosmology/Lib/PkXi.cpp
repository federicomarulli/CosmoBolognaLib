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
 *  @author federico.marulli3@unibo.it
 */

#include "FuncGrid_Bspline.h"
#include "Cosmology.h"
#include <regex>

using namespace std;

using namespace cbl;
using namespace cosmology;


// =====================================================================================


double cbl::cosmology::Cosmology::As (const double sigma8) const //check
{
  return pow(sigma8/1.79e4*pow(m_Omega_baryon*m_hh*m_hh/0.024, 1./3.)*pow(m_Omega_matter*m_hh*m_hh/0.14, -0.563)*pow(7.808*m_hh, 0.5*(1.-m_n_spec))*pow(m_hh/0.72, -0.693)*0.76/gg(0.), 2);
}


// =====================================================================================


double cbl::cosmology::Cosmology::sigma8_interpolated (const double redshift) const
{
  const double wm = m_Omega_matter*m_hh*m_hh;
  const double wb = m_Omega_baryon*m_hh*m_hh;
  const double sigma8 = 0.1058*pow(m_scalar_amp/2.196e-9, 0.5)*pow(wm/0.1426, 0.520)*pow(wb/0.02205, -0.294)*pow(m_hh/0.673, 0.683)*pow((m_massless_neutrinos+m_massive_neutrinos)/3.046, -0.24)*exp(0.3727*(m_n_spec-0.96))*pow(1-m_Omega_k, 0.175)*DN(redshift, 9.);

  return ((m_Omega_neutrinos>0) ? 0.995*sigma8 : sigma8);
}


// =====================================================================================


std::string cbl::cosmology::Cosmology::Pk_output_file (const string code, const bool NL, const double redshift, const bool run, const bool store_output, const string output_root, const double k_max, const string file_par)
{
  cbl::Path path;
  string dir_loc = path.fullpath(path.DirLoc());
  string dir_cosmo = path.DirCosmo();

  string dir_grid;
  if (NL==0) dir_grid = "output_linear/";
  else if (NL==1) dir_grid = "output_nonlinear/";
  else ErrorCBL("", "Pk_output_file", "PkXi.cpp");
  
  dir_grid += "h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6)+"/";

  //dir_grid += (output_root=="test") ? "/" : "_"+output_root;

  string dir = dir_cosmo+"/External/"+code+"/";
  string dirC = dir_cosmo+"/External/CAMB/";
  if (chdir(dirC.c_str())) {}

  string dir_output = dir+dir_grid+"pk.dat";

  if (run) {
    vector<double> lgkk, lgPk;
    Table_PkCodes(code, NL, lgkk, lgPk, redshift, store_output, output_root, k_max, file_par);
  }

  if (chdir(dir_loc.c_str())) {}

  return dir_output;
}


// =====================================================================================


void cbl::cosmology::Cosmology::run_CAMB (const bool NL, const double redshift, const std::string output_root, const std::string output_dir, const double k_max, const std::string file_par) const
{
  cbl::Path path;
  string dir_CAMB = path.DirCosmo()+"/External/CAMB/fortran/";

  string File_par = file_par;

  bool delete_output = (output_dir==par::defaultString) ? true : false;

  if (output_dir!=par::defaultString) if (system(("mkdir -p "+output_dir).c_str())) {}

  string OutputRoot = (output_dir==par::defaultString) ? dir_CAMB+"../inifiles/"+output_root : path.DirLoc()+output_dir+"/"+output_root;

  OutputRoot = (omp_get_max_threads()>1) ? OutputRoot+"_t"+conv(omp_get_thread_num(), par::fINT) : OutputRoot;

  if (File_par==par::defaultString) {

    // --------------------------------------------------------------------------
    // --------- set the cosmological parameters in the file params.ini ---------
    // --------------------------------------------------------------------------

    string file_par_default = dir_CAMB+"../inifiles/params_cut.ini";
    File_par = OutputRoot+"_params.ini";

    if (system(("cp "+file_par_default+" "+File_par).c_str())) {}
    ofstream fout(File_par.c_str(), std::ios_base::app | std::ios_base::out);
    double HH0 = m_hh*100.;

    fout << "output_root = " << OutputRoot << endl;
    fout << "do_nonlinear = " << conv(NL, par::fINT) << endl;
    fout << "hubble = " << conv(HH0, par::fDP6) << endl;
    fout << "ombh2 =" << conv(m_Omega_baryon*m_hh*m_hh, par::fDP6) << endl;
    fout << "omch2 = " << conv(m_Omega_CDM*m_hh*m_hh, par::fDP6) << endl;
    fout << "omk = " << conv(m_Omega_k, par::fDP6) << endl;
    fout << "omnuh2 = " << conv(m_Omega_neutrinos*m_hh*m_hh, par::fDP6) << endl;
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
    fout << "feedback_level = 1" << endl;
    fout << "print_sigma8 = T" << endl;
    fout << endl;

    fout.clear(); fout.close();

  }

  // --------------------------------------------------------------------------

  if (system((dir_CAMB+"camb "+File_par).c_str())) {}

  if (delete_output) {
    string RM = "rm -f "+OutputRoot+"*";
    if (system(RM.c_str())) {}
  }

}


// =====================================================================================


void cbl::cosmology::Cosmology::run_CAMB (std::vector<double> &kk, std::vector<double> &Pk, const bool NL, const double redshift, const std::string output_root, const std::string output_dir, const double k_max, const std::string file_par) const
{
  cbl::Path path;
  string dir_CAMB = path.DirCosmo()+"/External/CAMB/fortran/";
  string File_par = file_par;
  bool delete_output = (output_dir==par::defaultString) ? true : false;

  if (output_dir!=par::defaultString) if (system(("mkdir -p "+output_dir).c_str())) {}

  string OutputRoot = (output_dir==par::defaultString) ? dir_CAMB+"../inifiles/"+output_root : path.DirLoc()+output_dir+"/"+output_root;

  OutputRoot = (omp_get_max_threads() > 1) ? OutputRoot+"_t"+conv(omp_get_thread_num(), par::fINT) : OutputRoot;

  if (File_par==par::defaultString) {

    // --------------------------------------------------------------------------
    // --------- set the cosmological parameters in the file params.ini ---------
    // --------------------------------------------------------------------------

    string file_par_default = dir_CAMB+"../inifiles/params_cut.ini";
    File_par = OutputRoot+"_params.ini";

    if (system(("cp "+file_par_default+" "+File_par).c_str())) {}
    ofstream fout(File_par.c_str(), std::ios_base::app | std::ios_base::out);
    double HH0 = m_hh*100.;

    fout << "output_root = " << OutputRoot << endl;
    fout << "do_nonlinear = " << conv(NL, par::fINT) << endl;
    fout << "hubble = " << conv(HH0, par::fDP6) << endl;
    fout << "ombh2 =" << conv(m_Omega_baryon*m_hh*m_hh, par::fDP6) << endl;
    fout << "omch2 = " << conv(m_Omega_CDM*m_hh*m_hh, par::fDP6) << endl;
    fout << "omk = " << conv(m_Omega_k, par::fDP6) << endl;
    fout << "omnuh2 = " << conv(m_Omega_neutrinos*m_hh*m_hh, par::fDP6) << endl;
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
    fout << "print_sigma8 = F" << endl;
    fout << endl;

    fout.clear(); fout.close();
  }


  // --------------------------------------------------------------------------

  if (system((dir_CAMB+"camb "+File_par).c_str())) {}

  // read the power spectrum
  string file_inCAMB = OutputRoot+"_matterpower.dat";
  ifstream fin(file_inCAMB.c_str()); checkIO(fin, file_inCAMB);

  kk.erase(kk.begin(), kk.end());
  Pk.erase(Pk.begin(), Pk.end());

  string line;
  getline(fin, line);

  double KK, PK;
  while (fin >> KK >> PK)
    if (KK>0 && PK>0) {
      kk.push_back(KK);
      Pk.push_back(PK);
    }
  fin.clear(); fin.close();

  if (delete_output) if (system(("rm -f "+OutputRoot+"*").c_str())) {}
}


// =====================================================================================


void cbl::cosmology::Cosmology::Table_PkCodes (const std::string code, const bool NL, std::vector<double> &lgkk, std::vector<double> &lgPk, const double redshift, const bool store_output, const std::string output_root, const double k_max, const std::string file_par) const
{
  if (code=="MPTbreeze-v1") {
    if (m_sigma8<0)
      ErrorCBL("sigma8<0! The function set_sigma8() can be used to set the value of sigma8!", "Table_PkCodes", "PkXi.cpp");
    if (NL)
      WarningMsgCBL("NL is ignored by MPTbreeze-v1, that provides in output the non-linear power spectrum", "Table_PkCodes", "PkXi.cpp");
  }

  if (file_par==par::defaultString) {

    if (code=="CAMB" || code=="MGCAMB" || code=="MPTbreeze-v1") 
      m_Table_Pk_CAMB_MPTbreeze(code, NL, lgkk, lgPk, redshift, store_output, output_root, k_max);

    else if (code=="CLASS")
      m_Table_Pk_CLASS(NL, lgkk, lgPk, redshift, store_output, output_root, k_max);

    else
      ErrorCBL("the choosen code is not allowed!", "Table_PkCodes", "PkXi.cpp");
  }

  else {
    WarningMsgCBL("the input k_max parameters will not be used", "Table_PkCodes", "PkXi.cpp");
    m_Table_Pk_parameterFile(code, file_par, NL, lgkk, lgPk, redshift, output_root);
  }

}


// =====================================================================================


void cbl::cosmology::Cosmology::Table_PkCodes (const std::string code, const bool NL, std::vector<std::vector<double>> &lgkk, std::vector<std::vector<double>> &lgPk, const std::vector<double> redshift, const bool store_output, const std::string output_root, const double k_max, const std::string file_par) const
{
  if (code=="MPTbreeze-v1") {
    if (m_sigma8<0)
      ErrorCBL("sigma8<0! The function set_sigma8() can be used to set the value of sigma8!", "Table_PkCodes", "PkXi.cpp");
    if (NL)
      WarningMsgCBL("NL is ignored by MPTbreeze-v1, that provides in output the non-linear power spectrum", "Table_PkCodes", "PkXi.cpp");
  }

  lgkk.resize(redshift.size(), vector<double>(1, 0.));
  lgPk.resize(redshift.size(), vector<double>(1, 0.));

  if (file_par==par::defaultString) {

    if (code=="CAMB" || code=="MGCAMB" || code=="MPTbreeze-v1") 
      m_Table_Pk_CAMB_MPTbreeze(code, NL, lgkk, lgPk, redshift, store_output, output_root, k_max);

    else if (code=="CLASS")
      m_Table_Pk_CLASS(NL, lgkk, lgPk, redshift, store_output, output_root, k_max);
	
    else
      ErrorCBL("the choosen code is not allowed!", "Table_PkCodes", "PkXi.cpp");
 
  }

  else {
    WarningMsgCBL("the input k_max parameters will not be used", "Table_PkCodes", "PkXi.cpp");
    for (size_t zz=0; redshift.size(); zz++)
      m_Table_Pk_parameterFile(code, file_par, NL, lgkk[zz], lgPk[zz], redshift[zz], output_root);
  }

}


// =====================================================================================


void cbl::cosmology::Cosmology::m_Table_Pk_CAMB_MPTbreeze (const std::string code, const bool NL, std::vector<double> &lgkk, std::vector<double> &lgPk, const double redshift, const bool store_output, const std::string output_root, const double k_max) const
{
  lgkk.erase(lgkk.begin(), lgkk.end());
  lgPk.erase(lgPk.begin(), lgPk.end());

  string dir_grid;
  if (code=="CAMB" || code=="MGCAMB")
    dir_grid = (NL) ? "output_nonlinear/" : "output_linear/";
  else
    dir_grid = "output/";

  string filename = "h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6);

  filename += (output_root=="test") ? "/" : "_"+output_root+"/";

  string file_par;

  cbl::Path path;
  const string dirBZ = path.DirCosmo()+"/External/"+code+"/";
  const string dirBZ_output = dirBZ+dir_grid+filename;
  string fileBZ_in = dirBZ_output+"Pk.dat";
  ifstream fin_BZ;
  fin_BZ.open(fileBZ_in.c_str());

  const string nn = output_root+"_t"+conv(omp_get_thread_num(), par::fINT);
  const string tot_output_root = dirBZ+nn;
  string dir_output_root = dirBZ+nn;
  dir_output_root = regex_replace(dir_output_root, std::regex("/"), "\\/");


  // ------------------------------------------------------------------------------------------------------------
  // --------- set the CAMB parameter file if it does not exist, or if the MPTbreeze one does not exist ---------
  // ------------------------------------------------------------------------------------------------------------

  if (!fin_BZ && (code=="CAMB" || code=="MPTbreeze-v1")) {

    const string dirCAMB = path.DirCosmo()+"/External/CAMB/fortran/";
    
    file_par = dirCAMB+"../inifiles/params_"+nn+".ini";
    if (system(("cp "+dirCAMB+"../inifiles/params.ini "+file_par).c_str())) {}

    string sed;

    sed = "sed '/test/s//"+dir_output_root+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par;
    if (system(sed.c_str())) {}
    sed = "sed '/do_nonlinear = 0/s//do_nonlinear = "+conv(NL, par::fINT)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/hubble = 70/s//hubble = "+conv(m_hh*100., par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/ombh2 = 0.0226/s//ombh2 = "+conv(m_Omega_baryon*m_hh*m_hh, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/omch2 = 0.112/s//omch2 = "+conv(m_Omega_CDM*m_hh*m_hh, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/omk = 0/s//omk = "+conv(m_Omega_k, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/omnuh2 = 0.00064/s//omnuh2 = "+conv(m_Omega_neutrinos*m_hh*m_hh, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/transfer_redshift(1) = 0/s//transfer_redshift(1) = "+conv(redshift, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/massless_neutrinos = 2.046/s//massless_neutrinos = "+conv(m_massless_neutrinos, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/massive_neutrinos = 1/s//massive_neutrinos = "+conv(m_massive_neutrinos, par::fINT)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/scalar_spectral_index(1) = 0.96/s//scalar_spectral_index(1) = "+conv(m_n_spec, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/w = -1/s//w = "+conv(m_w0, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/wa = 0/s//wa = "+conv(m_wa, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    if (m_scalar_amp>0) {
      sed = "sed '/scalar_amp(1) = 2.1e-9/s//scalar_amp(1) = "+conv(m_scalar_amp, par::ee3)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
      sed = "sed '/pivot_scalar = 0.05/s//pivot_scalar = "+conv(m_scalar_pivot, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    }
    sed = "sed '/transfer_kmax = 2/s//transfer_kmax = "+conv(k_max, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/re_optical_depth = 0.09/s//re_optical_depth = "+conv(m_tau, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}

    // ----------------------------
    // --------- run CAMB ---------
    // ----------------------------
    if (system((dirCAMB+"camb "+file_par).c_str())) {}
  }

  else if (!fin_BZ && code=="MGCAMB") {

    coutCBL << "Set the parameters for modified gravity in the file External/MGCAMB/params_MG.ini" << endl;
 
    file_par = dirBZ+"params_"+nn+".ini";
    if (system(("cp "+dirBZ+"/params.ini "+file_par).c_str())) {}
    
    string sed;
   
    sed = "sed '/test/s//"+dir_output_root+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/do_nonlinear = 0/s//do_nonlinear = "+conv(NL, par::fINT)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/hubble = 70/s//hubble = "+conv(m_hh*100., par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/ombh2 = 0.0226/s//ombh2 = "+conv(m_Omega_baryon*m_hh*m_hh, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/omch2 = 0.112/s//omch2 = "+conv(m_Omega_CDM*m_hh*m_hh, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/omk = 0/s//omk = "+conv(m_Omega_k, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/omnuh2 = 0.00064/s//omnuh2 = "+conv(m_Omega_neutrinos*m_hh*m_hh, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/transfer_redshift(1) = 0/s//transfer_redshift(1) = "+conv(redshift, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}    
    sed = "sed '/massless_neutrinos = 2.046/s//massless_neutrinos = "+conv(m_massless_neutrinos, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/massive_neutrinos = 1/s//massive_neutrinos = "+conv(m_massive_neutrinos, par::fINT)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/scalar_spectral_index(1) = 0.96/s//scalar_spectral_index(1) = "+conv(m_n_spec, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/w = -1/s//w = "+conv(m_w0, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/wa = 0/s//wa = "+conv(m_wa, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    if (m_scalar_amp>0) {
      sed = "sed '/scalar_amp(1) = 2.1e-9/s//scalar_amp(1) = "+conv(m_scalar_amp, par::ee3)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
      sed = "sed '/pivot_scalar = 0.05/s//pivot_scalar = "+conv(m_scalar_pivot, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    }
    sed = "sed '/transfer_kmax = 2/s//transfer_kmax = "+conv(k_max, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/re_optical_depth = 0.09/s//re_optical_depth = "+conv(m_tau, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}

    // ----------------------------
    // -------- run MGCAMB --------
    // ----------------------------
    if (system((dirBZ+"camb "+file_par).c_str())) {}
    
  }

  // -----------------------------------------------------------------------------
  // --------- check if the MPTbreeze table exists, if not run MPTbreeze ---------
  // -----------------------------------------------------------------------------

  if (!fin_BZ && code=="MPTbreeze-v1") {

    // ---------------------------------
    // --------- run MPTbreeze ---------
    // ---------------------------------

    if (system((dirBZ+"mptbreeze -noverbose -camb "+file_par+" -fileTF "+tot_output_root+"_transfer_out.dat -sigma8 "+conv(m_sigma8, par::fDP6)+" -omegam "+conv(m_Omega_matter, par::fDP6)+" -ns "+conv(m_n_spec, par::fDP6)+" -w "+conv(m_w0, par::fDP6)+" -redshift "+conv(redshift, par::fDP6)+" -filePk "+tot_output_root+"_matterpower.dat").c_str())) {}
    
  }

  bool remove = false;
  
  if (!fin_BZ) {
    if (store_output) {
      if (system(("mkdir -p "+dirBZ_output).c_str())) {}
      if (system(("mv "+tot_output_root+"_matterpower.dat "+dirBZ_output+"Pk.dat").c_str())) {}
      if (system(("mv "+tot_output_root+"_transfer_out.dat "+dirBZ_output+"transfer_out.dat").c_str())) {}
    }
    else {
      fileBZ_in = tot_output_root+"_matterpower.dat";
      remove = true;
    }
  }

  if (system(("rm -f "+file_par+" "+tot_output_root+"_params.ini").c_str())) {}
  fin_BZ.clear(); fin_BZ.close();


  // ---------------------------------------
  // --------- get the output P(k) ---------
  // ---------------------------------------

  ifstream fin(fileBZ_in.c_str()); checkIO(fin, fileBZ_in);
  
  string line;

  while (getline(fin, line)) {
    if (line.find("#")!=0) {
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
  }
  
  if (remove) if (system(("rm "+fileBZ_in+" "+tot_output_root+"_transfer_out.dat").c_str())) {}
  fin.clear(); fin.close();

  if (lgkk.size()==0 || lgPk.size()==0)
    ErrorCBL("lgkk.size()="+conv(lgkk.size(), par::fINT)+", lgPk.size()="+conv(lgPk.size(), par::fINT), "m_Table_Pk_CAMB_MPTbreeze", "PkXi.cpp");

}


// =====================================================================================


void cbl::cosmology::Cosmology::m_Table_Pk_CAMB_MPTbreeze (const std::string code, const bool NL, std::vector<std::vector<double>> &lgkk, std::vector<std::vector<double>> &lgPk, const std::vector<double> redshift, const bool store_output, const std::string output_root, const double k_max) const
{
  string dir_grid;
  if (code=="CAMB" || code=="MGCAMB")
    dir_grid = (NL) ? "output_nonlinear/" : "output_linear/";
  else
    dir_grid = "output/";

  string add_string;
  add_string = (output_root=="test") ? "/" : "_"+output_root+"/";

  vector<double> zz = redshift;
  std::sort(zz.begin(), zz.end(), greater<double>());
  vector<string> filename_z(zz.size());
  vector<string> filename_z_original(redshift.size());
  
  for (size_t ii=0; ii<zz.size(); ii++) {
    filename_z[ii] = "h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(zz[ii], par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6);
    filename_z[ii]+=add_string;
  }

  for (size_t ii=0; ii<redshift.size(); ii++) {
    filename_z_original[ii] = "h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift[ii], par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6);
    filename_z_original[ii]+=add_string;
  }

  string file_par;

  cbl::Path path;
  const string dirBZ = path.DirCosmo()+"/External/"+code+"/";
  const string dirBZ_output = dirBZ+dir_grid;

  // check if the power spectrum has been already stored for that cosmology and redshift
  for (size_t ii=zz.size(); ii-->0;) {
    string fileBZ_in = dirBZ+dir_grid+filename_z[ii]+"Pk.dat";
    ifstream fin_BZ;
    fin_BZ.open(fileBZ_in.c_str());
    if (fin_BZ) { // remove the redshift from the vector if already computed
      zz.erase(zz.begin()+ii); 
      filename_z.erase(filename_z.begin()+ii);
    }
    fin_BZ.clear(); fin_BZ.close();
  }

  const string nn = output_root+"_t"+conv(omp_get_thread_num(), par::fINT);
  const string tot_output_root = dirBZ+nn;
  string dir_output_root = dirBZ+nn;
  dir_output_root = regex_replace(dir_output_root, std::regex("/"), "\\/");

  // ------------------------------------------------------------------------------------------------------------
  // --------- set the CAMB parameter file if it does not exist, or if the MPTbreeze one does not exist ---------
  // ------------------------------------------------------------------------------------------------------------

  if (zz.size()>0 && (code=="CAMB" || code=="MPTbreeze-v1")) {

    cbl::Path path;
    const string dirCAMB = path.DirCosmo()+"/External/CAMB/fortran/";
    
    file_par = dirCAMB+"../inifiles/params_"+nn+".ini";
    if (system(("cp "+dirCAMB+"../inifiles/params.ini "+file_par).c_str())) {}

    string sed;

    sed = "sed '/test/s//"+dir_output_root+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par;
    if (system(sed.c_str())) {}
    sed = "sed '/do_nonlinear = 0/s//do_nonlinear = "+conv(NL, par::fINT)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/hubble = 70/s//hubble = "+conv(m_hh*100., par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/ombh2 = 0.0226/s//ombh2 = "+conv(m_Omega_baryon*m_hh*m_hh, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/omch2 = 0.112/s//omch2 = "+conv(m_Omega_CDM*m_hh*m_hh, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/omk = 0/s//omk = "+conv(m_Omega_k, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/omnuh2 = 0.00064/s//omnuh2 = "+conv(m_Omega_neutrinos*m_hh*m_hh, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/massless_neutrinos = 2.046/s//massless_neutrinos = "+conv(m_massless_neutrinos, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/massive_neutrinos = 1/s//massive_neutrinos = "+conv(m_massive_neutrinos, par::fINT)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/scalar_spectral_index(1) = 0.96/s//scalar_spectral_index(1) = "+conv(m_n_spec, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}

    string transfer_z = "//";
    string transfer_file = "//";
    string transfer_matter = "//";
    for (size_t ii=0; ii<zz.size(); ii++) {
      transfer_z+="transfer_redshift("+std::to_string(ii+1)+") = "+conv(zz[ii], par::fDP6)+" \\\n";
      transfer_file+="transfer_filename("+std::to_string(ii+1)+") = transfer_out_"+std::to_string(ii+1)+".dat \\\n";
      transfer_matter+="transfer_matterpower("+std::to_string(ii+1)+") = matterpower_"+std::to_string(ii+1)+".dat \\\n";
    }

    sed = "sed '/transfer_num_redshifts  = 1/s//transfer_num_redshifts  = "+conv(zz.size(), par::fINT)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/transfer_redshift(1) = 0/s"+transfer_z+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/transfer_filename(1) = transfer_out.dat/s"+transfer_file+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/transfer_matterpower(1) = matterpower.dat/s"+transfer_matter+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    
    sed = "sed '/w = -1/s//w = "+conv(m_w0, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/wa = 0/s//wa = "+conv(m_wa, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    if (m_scalar_amp>0) {
      sed = "sed '/scalar_amp(1) = 2.1e-9/s//scalar_amp(1) = "+conv(m_scalar_amp, par::ee3)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
      sed = "sed '/pivot_scalar = 0.05/s//pivot_scalar = "+conv(m_scalar_pivot, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    }
    sed = "sed '/transfer_kmax = 2/s//transfer_kmax = "+conv(k_max, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/re_optical_depth = 0.09/s//re_optical_depth = "+conv(m_tau, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}

    // ----------------------------
    // --------- run CAMB ---------
    // ----------------------------

    if (system((dirCAMB+"camb "+file_par).c_str())) {}
  }

  else if (zz.size()>0 && code=="MGCAMB") {

    coutCBL << "Set the parameters for modified gravity in the file External/MGCAMB/params_MG.ini" << endl;
 
    file_par = dirBZ+"params_"+nn+".ini";
    if (system(("cp "+dirBZ+"/params.ini "+file_par).c_str())) {}
    
    string sed;
   
    sed = "sed '/test/s//"+dir_output_root+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/do_nonlinear = 0/s//do_nonlinear = "+conv(NL, par::fINT)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/hubble = 70/s//hubble = "+conv(m_hh*100., par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/ombh2 = 0.0226/s//ombh2 = "+conv(m_Omega_baryon*m_hh*m_hh, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/omch2 = 0.112/s//omch2 = "+conv(m_Omega_CDM*m_hh*m_hh, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/omk = 0/s//omk = "+conv(m_Omega_k, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/omnuh2 = 0.00064/s//omnuh2 = "+conv(m_Omega_neutrinos*m_hh*m_hh, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}

    string transfer_z = "//";
    string transfer_file = "//";
    string transfer_matter = "//";
    for (size_t ii=0; ii<zz.size(); ii++) {
      transfer_z+="transfer_redshift("+std::to_string(ii+1)+") = "+conv(zz[ii], par::fDP6)+" \\\n";
      transfer_file+="transfer_filename("+std::to_string(ii+1)+") = transfer_out_"+std::to_string(ii+1)+".dat \\\n";
      transfer_matter+="transfer_matterpower("+std::to_string(ii+1)+") = matterpower_"+std::to_string(ii+1)+".dat \\\n";
    }
    sed = "sed '/transfer_num_redshifts = 1/s//transfer_num_redshifts  = "+conv(zz.size(), par::fINT)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/transfer_redshift(1) = 0/s"+transfer_z+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/transfer_filename(1) = transfer_out.dat/s"+transfer_file+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/transfer_matterpower(1) = matterpower.dat/s"+transfer_matter+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}

    sed = "sed '/massless_neutrinos = 2.046/s//massless_neutrinos = "+conv(m_massless_neutrinos, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/massive_neutrinos = 1/s//massive_neutrinos = "+conv(m_massive_neutrinos, par::fINT)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/scalar_spectral_index(1) = 0.96/s//scalar_spectral_index(1) = "+conv(m_n_spec, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/w = -1/s//w = "+conv(m_w0, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/wa = 0/s//wa = "+conv(m_wa, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    if (m_scalar_amp>0) {
      sed = "sed '/scalar_amp(1) = 2.1e-9/s//scalar_amp(1) = "+conv(m_scalar_amp, par::ee3)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
      sed = "sed '/pivot_scalar = 0.05/s//pivot_scalar = "+conv(m_scalar_pivot, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    }
    sed = "sed '/transfer_kmax = 2/s//transfer_kmax = "+conv(k_max, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/re_optical_depth = 0.09/s//re_optical_depth = "+conv(m_tau, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    
    // ----------------------------
    // -------- run MGCAMB --------
    // ----------------------------

    if (system((dirBZ+"camb "+file_par).c_str())) {}
    
  }
    
  // -----------------------------------------------------------------------------
  // --------- check if the MPTbreeze table exists, if not run MPTbreeze ---------
  // -----------------------------------------------------------------------------

  if (zz.size()>0 && code=="MPTbreeze-v1") {

    // ---------------------------------
    // --------- run MPTbreeze ---------
    // ---------------------------------
    
    for (size_t ii=0; ii<zz.size(); ii++) {
      if (system((dirBZ+"mptbreeze -noverbose -camb "+file_par+" -fileTF "+tot_output_root+"_transfer_out_"+std::to_string(ii+1)+".dat -sigma8 "+conv(m_sigma8, par::fDP6)+" -omegam "+conv(m_Omega_matter, par::fDP6)+" -ns "+conv(m_n_spec, par::fDP6)+" -w "+conv(m_w0, par::fDP6)+" -redshift "+conv(zz[ii], par::fDP6)+" -filePk "+tot_output_root+"_matterpower_"+std::to_string(ii+1)+".dat").c_str())) {}
    }

  }

  if (store_output)
    for (size_t ii=0; ii<zz.size(); ii++) {
      if (system(("mkdir -p "+dirBZ_output+filename_z[ii]).c_str())) {}
      if (system(("mv "+tot_output_root+"_matterpower_"+std::to_string(ii+1)+".dat "+dirBZ_output+filename_z[ii]+"Pk.dat").c_str())) {}
      if (system(("mv "+tot_output_root+"_transfer_out_"+std::to_string(ii+1)+".dat "+dirBZ_output+filename_z[ii]+"transfer_out.dat").c_str())) {}
    }
  
  if (system(("rm -f "+file_par+" "+tot_output_root+"_params.ini").c_str())) {}

  // ---------------------------------------
  // --------- get the output P(k) ---------
  // ---------------------------------------

  int counter = 0;
  
  for (size_t ii=redshift.size(); ii-->0;) {

    lgkk[ii].erase(lgkk[ii].begin(), lgkk[ii].end());
    lgPk[ii].erase(lgPk[ii].begin(), lgPk[ii].end());

    string fileBZ_in = dirBZ+dir_grid+filename_z_original[ii]+"Pk.dat";
    ifstream fin_check(fileBZ_in.c_str());
    
    if (fin_check and !store_output) counter++;
    if (!fin_check) fileBZ_in = tot_output_root+"_matterpower_"+std::to_string(redshift.size()-ii-counter)+".dat";

    ifstream fin(fileBZ_in.c_str());
    string line;

    while (getline(fin, line)) {
      
      if (line.find("#")!=0) {
	stringstream ss(line);
	vector<double> num;
	double aa;
	while (ss>>aa) num.push_back(aa);

	if (num[0]>0 && num[1]>0) {

	  lgkk[ii].push_back(log10(num[0]));
	  if (code=="MPTbreeze-v1")
	    lgPk[ii].push_back(log10(pow(2.*par::pi, 3)*(num[1]+num[2]+num[3])));
	  else
	    lgPk[ii].push_back(log10(num[1]));

	}
      }
    }

    if (!fin_check) if (system(("rm "+fileBZ_in+" "+tot_output_root+"_transfer_out_"+std::to_string(redshift.size()-ii-counter)+".dat").c_str())) {}

    fin_check.clear(); fin_check.close();
    fin.clear(); fin.close();
    
    if (lgkk[ii].size()==0 || lgPk[ii].size()==0)
      ErrorCBL("lgkk.size()="+conv(lgkk[ii].size(), par::fINT)+", lgPk.size()="+conv(lgPk[ii].size(), par::fINT), "m_Table_Pk_CAMB_MPTbreeze", "PkXi.cpp");

  }
  
}


// =====================================================================================


void cbl::cosmology::Cosmology::m_Table_Pk_CLASS (const bool NL, std::vector<double> &lgkk, std::vector<double> &lgPk, const double redshift, const bool store_output, const std::string output_root, const double k_max) const
{
  lgkk.erase(lgkk.begin(), lgkk.end());
  lgPk.erase(lgPk.begin(), lgPk.end());

  cbl::Path path;
  const string dirC = path.DirCosmo()+"/External/CLASS/";
  
  string dir_grid = (NL) ? dirC+"output_nonlinear/" : dirC+"output_linear/";

  dir_grid += "h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6);

  dir_grid += (output_root=="test") ? "/" : "_"+output_root+"/";

  const string nn = output_root+"_t"+conv(omp_get_thread_num(), par::fINT);
  const string dir_output = dirC+nn+"_00";

  string file_in = dir_grid+"Pk.dat";
  ifstream fin;
  fin.open(file_in.c_str());

  string file_par;

  if (!fin) {

    // --------------------------------------------------------------------------
    // --------- set the cosmological parameters in the file params.ini ---------
    // --------------------------------------------------------------------------

    string dir_output_root = dirC+nn;
    dir_output_root = regex_replace(dir_output_root, std::regex("/"), "\\/");

    file_par = dirC+"params_"+nn+".ini";
    if (system(("cp "+dirC+"params.ini "+file_par).c_str())) {}

    string sed;

    sed = "sed '/output\\/test_/s//"+dir_output_root+"_/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}

    if (NL) {
      sed = "sed '/non linear =/s//non linear = halofit/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+"";
      if (system(sed.c_str())) {}
    }

    sed = "sed '/h = 0.67556/s//h = "+conv(m_hh, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/Omega_b = 0.022032/s//Omega_b = "+conv(m_Omega_baryon, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/Omega_cdm = 0.12038/s//Omega_cdm = "+conv(m_Omega_CDM, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/N_ncdm = 0/s//N_ncdm = "+conv(m_massive_neutrinos, par::fINT)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/Omega_ncdm = 0/s//Omega_ncdm = "+conv(m_Omega_neutrinos, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}

    if (m_Omega_neutrinos>0) { // check!!!
      sed = "sed '/m_ncdm = 0.04, 0.04, 0.04/s//#m_ncdm/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+"";
      if (system(sed.c_str())) {}
    }

    sed = "sed '/Omega_Lambda = 0.7/s//Omega_Lambda = "+conv(m_Omega_DE, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/Omega_k = 0./s//Omega_k = "+conv(m_Omega_k, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/z_pk = 0/s//z_pk = "+conv(redshift, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/N_ur = 3.046/s//N_ur = "+conv(m_massless_neutrinos, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/N_ncdm = 0/s//N_ncdm = "+conv(m_massive_neutrinos, par::fINT)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/n_s = 0.9619/s//n_s = "+conv(m_n_spec, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}

    double w00 = max(-0.999,m_w0); // check!!!
    sed = "sed '/w0_fld = -0.9/s//w0_fld = "+conv(w00, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}

    sed = "sed '/wa_fld = 0./s//wa_fld = "+conv(m_wa, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    if (m_scalar_amp>0) {
      sed = "sed '/A_s = 2.215e-9/s//A_s = "+conv(m_scalar_amp, par::ee3)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
      sed = "sed '/k_pivot = 0.05/s//k_pivot = "+conv(m_scalar_pivot, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    }
    sed = "sed '/P_k_max_h\\/Mpc = 1/s//P_k_max_h\\/Mpc = "+conv(k_max, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/tau_reio = 0.0925/s//tau_reio = "+conv(m_tau, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}

    // -----------------------------
    // --------- run CLASS ---------
    // -----------------------------

    if (system((dirC+"/class "+file_par).c_str())) {}
    if (system(("rm -f "+file_par).c_str())) {}
  }

  fin.clear(); fin.close();

  bool remove = false;
  
  if (store_output) {
    if (system(("mkdir -p "+dir_grid).c_str())) {}
    if (system(("mv "+dir_output+"_pk"+((NL) ? "_nl.dat" : ".dat")+" "+dir_grid+"Pk.dat").c_str())) {}
  }

  else {
    file_in = dir_output+"_pk"+((NL) ? "_nl.dat" : ".dat");
    remove = true;
  }
  
  // ---------------------------------------
  // --------- get the output P(k) ---------
  // ---------------------------------------

  fin.open(file_in.c_str()); checkIO(fin, file_in);
  string line;

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

  if (remove) if (system(("rm -f "+dir_output+"_pk*.dat").c_str())) {}
  fin.clear(); fin.close();

  if (lgkk.size()==0 || lgPk.size()==0)
    ErrorCBL("lgkk.size()="+conv(lgkk.size(), par::fINT)+", lgPk.size()="+conv(lgPk.size(), par::fINT), "m_Table_Pk_CLASS", "PkXi.cpp");

}


// =====================================================================================


void cbl::cosmology::Cosmology::m_Table_Pk_CLASS (const bool NL, std::vector<std::vector<double>> &lgkk, std::vector<std::vector<double>> &lgPk, const std::vector<double> redshift, const bool store_output, const std::string output_root, const double k_max) const
{
  cbl::Path path;
  const string dirC = path.DirCosmo()+"/External/CLASS/";
  string dir_grid = (NL) ? dirC+"output_nonlinear/" : dirC+"output_linear/";

  vector<string> filename_z_original(redshift.size());
  vector<string> dir_output_z_original(redshift.size());


  for (size_t ii=0; ii<redshift.size(); ii++) {
    filename_z_original[ii] = dir_grid+"h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift[ii], par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6);
    filename_z_original[ii] += (output_root=="test") ? "/" : "_"+output_root+"/";
    dir_output_z_original[ii] = filename_z_original[ii];
  }

  vector<double> zz = redshift;

  for (size_t ii=zz.size(); ii-->0;) {
    string file_in = dir_output_z_original[ii]+"Pk.dat";
    ifstream fin;
    fin.open(file_in.c_str());
    if (fin) zz.erase(zz.begin()+ii); // remove the redshift from the vector if already computed
    fin.clear(); fin.close();
  }

  vector<string> filename_z(zz.size());
  vector<string> dir_output_z(zz.size());

  for (size_t ii=0; ii<zz.size(); ii++) {
    filename_z[ii] = dir_grid+"h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(zz[ii], par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6);
    filename_z[ii] += (output_root=="test") ? "/" : "_"+output_root+"/";
    dir_output_z[ii] = filename_z[ii];
  }

  const string nn = output_root+"_t"+conv(omp_get_thread_num(), par::fINT);

  if (zz.size()>0) {

    // --------------------------------------------------------------------------
    // --------- set the cosmological parameters in the file params.ini ---------
    // --------------------------------------------------------------------------

    string dir_output_root = dirC+nn;
    dir_output_root = regex_replace(dir_output_root, std::regex("/"), "\\/");

    string file_par = dirC+"params_"+output_root+"_t"+conv(omp_get_thread_num(), par::fINT)+".ini";
    if (system(("cp "+dirC+"params.ini "+file_par).c_str())) {}

    string sed;

    sed = "sed '/output\\/test_/s//"+dir_output_root+"_/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}

    if (NL) {
      sed = "sed '/non linear =/s//non linear = halofit/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+"";
      if (system(sed.c_str())) {}
    }

    sed = "sed '/h = 0.67556/s//h = "+conv(m_hh, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/Omega_b = 0.022032/s//Omega_b = "+conv(m_Omega_baryon, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/Omega_cdm = 0.12038/s//Omega_cdm = "+conv(m_Omega_CDM, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/N_ncdm = 0/s//N_ncdm = "+conv(m_massive_neutrinos, par::fINT)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/Omega_ncdm = 0/s//Omega_ncdm = "+conv(m_Omega_neutrinos, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}

    if (m_Omega_neutrinos>0) { // check!!!
      sed = "sed '/m_ncdm = 0.04, 0.04, 0.04/s//#m_ncdm/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+"";
      if (system(sed.c_str())) {}
    }

    sed = "sed '/Omega_Lambda = 0.7/s//Omega_Lambda = "+conv(m_Omega_DE, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}
    sed = "sed '/Omega_k = 0./s//Omega_k = "+conv(m_Omega_k, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par+""; if (system(sed.c_str())) {}

    string z_pk = "z_pk = ";
    for (size_t ii=0; ii<zz.size(); ii++) {
      z_pk+=conv(zz[ii], par::fDP6);
      if (ii<zz.size()-1) z_pk += ",";
    }
    
    sed = "sed '/z_pk = 0/s//"+z_pk+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    
    sed = "sed '/N_ur = 3.046/s//N_ur = "+conv(m_massless_neutrinos, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/N_ncdm = 0/s//N_ncdm = "+conv(m_massive_neutrinos, par::fINT)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/n_s = 0.9619/s//n_s = "+conv(m_n_spec, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}

    double w00 = max(-0.999,m_w0); // check!!!
    sed = "sed '/w0_fld = -0.9/s//w0_fld = "+conv(w00, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}

    sed = "sed '/wa_fld = 0./s//wa_fld = "+conv(m_wa, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    if (m_scalar_amp>0) {
      sed = "sed '/A_s = 2.215e-9/s//A_s = "+conv(m_scalar_amp, par::ee3)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
      sed = "sed '/k_pivot = 0.05/s//k_pivot = "+conv(m_scalar_pivot, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    }
    sed = "sed '/P_k_max_h\\/Mpc = 1/s//P_k_max_h\\/Mpc = "+conv(k_max, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}
    sed = "sed '/tau_reio = 0.0925/s//tau_reio = "+conv(m_tau, par::fDP6)+"/g' "+file_par+" > temp_"+nn+"; mv temp_"+nn+" "+file_par; if (system(sed.c_str())) {}


    // -----------------------------
    // --------- run CLASS ---------
    // -----------------------------

    if (system((dirC+"/class "+file_par).c_str())) {}
    if (system(("rm -f "+file_par).c_str())) {}
    
    if (store_output) {
      for (size_t ii=0; ii<zz.size(); ii++) {
	if (system(("mkdir -p "+dir_output_z[ii]).c_str())) {}

	if (zz.size()==1) {if (system(("mv "+dirC+nn+"_00_pk"+((NL) ? "_00_nl.dat" : ".dat")+" "+dir_output_z[ii]+"Pk.dat").c_str())) {}}
	else {if (system(("mv "+dirC+nn+"_00_z"+std::to_string(ii+1)+((NL) ? "_pk_nl.dat" : "_pk.dat")+" "+dir_output_z[ii]+"Pk.dat").c_str())) {}}
      }
    }
 
  }
  
  // ---------------------------------------
  // --------- get the output P(k) ---------
  // ---------------------------------------

  int counter = 0;
  
  for (size_t ii=0; ii<redshift.size(); ii++) {

    lgkk[ii].erase(lgkk[ii].begin(), lgkk[ii].end());
    lgPk[ii].erase(lgPk[ii].begin(), lgPk[ii].end());

    string file_in = dir_output_z_original[ii]+"Pk.dat";
    ifstream fin_check(file_in.c_str());

    if (fin_check and !store_output) counter++;
    if (!fin_check) {
      if (zz.size()==1) file_in = dirC+nn+"_00_pk"+((NL) ? "_00_nl.dat" : ".dat");
      else file_in = dirC+nn+"_00_z"+std::to_string(ii+1-counter)+((NL) ? "_pk_nl.dat" : "_pk.dat");
    }

    ifstream fin(file_in.c_str());
    string line;

    while (getline(fin, line)) {
      stringstream ss(line);
      vector<double> num;
      double aa;
      while (ss>>aa) num.push_back(aa);

      if (num[0]>0 && num[1]>0) {
	lgkk[ii].push_back(log10(num[0]));
	lgPk[ii].push_back(log10(num[1]));
      }
    }
    
    fin_check.clear(); fin_check.close();
    fin.clear(); fin.close();

    if (lgkk[ii].size()==0 || lgPk[ii].size()==0)
      ErrorCBL("lgkk.size()="+conv(lgkk[ii].size(), par::fINT)+", lgPk.size()="+conv(lgPk[ii].size(), par::fINT), "m_Table_Pk_CLASS", "PkXi.cpp");
  }

  if (system(("rm -f "+dirC+nn+"_00_*.dat").c_str())) {} // to remove also *_cb.dat files

}


// =====================================================================================


void cbl::cosmology::Cosmology::remove_output_Pk_tables (const string code, const bool NL, const double redshift, const string output_root) const
{ 
  string dir_grid;
  if (code=="CAMB" || code=="MGCAMB" || code=="CLASS")
    dir_grid = (NL) ? "output_nonlinear/" : "output_linear/";
  else
    dir_grid = "output/";

  dir_grid += "h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6);

  if (output_root!="test") dir_grid += "_"+output_root;

  string dir_output;
  
  cbl::Path path;
  if (code=="CAMB") dir_output = path.DirCosmo()+"/External/CAMB/"+dir_grid;
  if (code=="MGCAMB") dir_output = path.DirCosmo()+"/External/MGCAMB/"+dir_grid;
  else if (code=="CLASS") dir_output = path.DirCosmo()+"/External/CLASS/"+dir_grid;
  else dir_output = path.DirCosmo()+"/External/MPTbreeze-v1/"+dir_grid;

  if (system(("rm -rf "+dir_output+" > /dev/null 2>&1").c_str())) {}

}


// =====================================================================================


void cbl::cosmology::Cosmology::m_Table_Pk_parameterFile (const std::string code, const std::string file_par, const bool NL, std::vector<double> &lgkk, std::vector<double> &lgPk, const double redshift, const std::string output_root) const
{
  WarningMsgCBL("Check the consistency between the output_root/root parameter given in input and the one set in the parameter file", "m_Table_Pk_parameterFile", "PkXi.cpp");

  lgkk.erase(lgkk.begin(), lgkk.end());
  lgPk.erase(lgPk.begin(), lgPk.end());

  cbl::Path path;
  const string dir_loc = path.fullpath(path.DirLoc());
  const string dirB = (code=="CAMB" || code=="MPTbreeze-v1") ? path.DirCosmo()+"/External/CAMB/fortran/"
    : path.DirCosmo()+"/External/"+code+"/";

  const string dir_output = (code=="CAMB" || code=="MPTbreeze-v1") ? dirB+"../output_ParameterFiles/"+file_par+"/"
    : dirB+"output_ParameterFiles/"+file_par+"/";

  if (system(("mkdir -p "+dir_output).c_str())) {}

  const string file_out = dir_output+"Pk.dat";
  
  const string dirPar = (code=="CAMB" || code=="MPTbreeze-v1") ? "../inifiles/" : "./";
  const string File_par = dirB+dirPar+file_par;

  ifstream check;
  check.open(File_par.c_str());
  if (!check) ErrorCBL("Parameter file not found! Make sure that this file is located in CosmoBolognaLib/External/CAMB/inifiles/ for CAMB and MPTbreeze-v1, and in CosmoBolognaLib/External/_NameBoltzmannSolver_/ for the remaining cases", "m_Table_Pk_parameterFile", "PkXi.cpp");
  check.clear(); check.close();
  
  ifstream fin;
  fin.open(file_out.c_str());

  if (!fin) {

    if (chdir(dirB.c_str())) {}
    
    // --------------------------------------------
    // --------- run the Boltzmann solver ---------
    // --------------------------------------------

    if (code=="CAMB" || code=="MGCAMB" || code=="MPTbreeze-v1") {
      if (system(("./camb "+File_par).c_str())) {}
      if (system(("mv "+output_root+"_matterpower*dat "+dir_output+"Pk.dat").c_str())) {}

      if (code=="MPTbreeze-v1") {
	WarningMsgCBL("Check the consistency of the input redshift with the one set in the parameter file", "m_Table_Pk_parameterFile", "PkXi.cpp");
	if (chdir((path.DirCosmo()+"/External/"+code).c_str())) {}

	if (system(("./mptbreeze -noverbose -camb ../../CAMB/"+File_par+" -fileTF ../../CAMB/"+output_root+"_transfer_out.dat -sigma8 "+conv(m_sigma8, par::fDP6)+" -omegam "+conv(m_Omega_matter, par::fDP6)+" -ns "+conv(m_n_spec, par::fDP6)+" -w "+conv(m_w0, par::fDP6)+" -redshift "+conv(redshift, par::fDP6)+" -filePk ../../CAMB/"+output_root+"_matterpower.dat").c_str())) {}
	if (system(("mv ../../CAMB/"+output_root+"_matterpower*dat "+dir_output+"Pk.dat").c_str())) {}
      }

    }

    else if (code=="CLASS") {
      WarningMsgCBL("Check the consistency of the input NL value with the one set in the parameter file", "m_Table_Pk_parameterFile", "PkXi.cpp");
      if (system(("./class "+File_par).c_str())) {}
      if (system(("mv "+output_root+"pk"+((NL) ? ".dat" : "_nl.dat")+" "+dir_output+"Pk.dat").c_str())) {}
    }

    else ErrorCBL("the chosen code is not allowed!", "m_Table_parameterFile", "PkXi.cpp");

    if (system(("rm -f "+output_root+"*_params.ini "+output_root+"*_transfer*" ).c_str())) {}
    if (chdir(dir_loc.c_str())) {}

  }
  
  fin.clear(); fin.close();

  // ---------------------------------------
  // --------- get the output P(k) ---------
  // ---------------------------------------

  fin.open(file_out.c_str()); checkIO(fin, file_out);
  string line;

  while (getline(fin, line)) {
    if (line.find("#")!=0) {
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
  }

  fin.clear(); fin.close();

  if (lgkk.size()==0 || lgPk.size()==0)
    ErrorCBL("lgkk.size()="+conv(lgkk.size(), par::fINT)+", lgPk.size()="+conv(lgPk.size(), par::fINT), "Table_PkCodes", "PkXi.cpp");

  if (chdir(dir_loc.c_str())) {}
}



// =====================================================================================


void cbl::cosmology::Cosmology::Table_XiCodes (const std::string code, const bool NL, std::vector<double> &rr, std::vector<double> &xi, const double redshift, const bool store_output, const std::string output_root, const double k_max, const std::string file_par) const
{
  vector<double> _lgkk, _lgPk;
  Table_PkCodes(code, NL, _lgkk, _lgPk, redshift, store_output, output_root, k_max, file_par);

  if (_lgkk.size()==0 || _lgPk.size()==0)
    ErrorCBL("_lgkk.size()="+conv(_lgkk.size(), par::fINT)+", _lgPk.size()="+conv(_lgPk.size(), par::fINT), "Table_XiCodes", "PkXi.cpp");

  // if the size of the input dataset is less than 1000, interpolate the values in a denser grid
  //WarningMsgCBL("the number of k-P(k) grid points will be increased to 1000, by interpolation", "Table_XiCodes", "PkXi.cpp");
  vector<double> lgkk, lgPk;
  if (_lgkk.size()<1000) {
    lgkk = linear_bin_vector(1000, min(1.e-4, Min(_lgkk)), max(1., Max(_lgkk)));
    for (size_t i=0; i<lgkk.size(); ++i)
      lgPk.emplace_back(interpolated(lgkk[i], _lgkk, _lgPk, "Linear"));
  }
  else {
    lgkk = _lgkk;
    lgPk = _lgPk;
  }
  cbl::Path path;
  const string dir_loc = path.DirLoc();
  const string dir_cosmo = path.DirCosmo();

  string dir_grid;
  if (code=="MPTbreeze-v1") dir_grid = "output/";
  else if (NL==false) dir_grid = "output_linear/";
  else dir_grid = "output_nonlinear/";

  dir_grid += "h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6);

  if (output_root!="test") dir_grid += "_"+output_root;

  string dir = dir_cosmo+"/External/"+code+"/";
  string dirFFT = dir_cosmo+"/External/fftlog-f90-master/";
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

    wrapper::fftlog::transform_FFTlog(rr, xi, 1, kk, Pk);

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

  else if (method_Pk == "CAMB" || method_Pk=="CAMB_wrapper" || method_Pk=="MGCAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="CLASS") {

    vector<double> lgkk, lgPk;

    if (method_Pk == "CAMB_wrapper") {
      
      const int npoints = 500;
      vector<double> Pk = wrapper::camb::Pk_CAMB (false, redshift, k_min, k_max, npoints, m_Omega_baryon*m_hh*m_hh, m_Omega_CDM*m_hh*m_hh, m_Omega_neutrinos*m_hh*m_hh, m_massless_neutrinos, m_massive_neutrinos, m_Omega_k, m_hh*100., m_n_spec, m_scalar_amp, m_scalar_pivot, m_w0, m_wa, m_tau);      
      vector<double> kk = logarithmic_bin_vector(npoints, k_min, k_max);
      
      for (size_t i=0; i<Pk.size(); i++) {
	lgkk.emplace_back(log10(kk[i]/m_hh));
	lgPk.emplace_back(log10(Pk[i]));
      }
      
    }
    
    else
      Table_PkCodes(method_Pk, NL, lgkk, lgPk, redshift, store_output, output_root, k_max, file_par);

    int limit_size = 1000;
    gsl_integration_workspace *ww = gsl_integration_workspace_alloc(limit_size);
    gsl_function Func;

    glob::STR_SSM str;
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

  const double D_N = DN(redshift);
  
  if (method_Pk=="EisensteinHu") m_Pk0_EH = 2.*pow(par::pi*m_sigma8,2)/Int*D_N*D_N;
  if (method_Pk=="CAMB" || method_Pk=="CAMB_wrapper" || method_Pk=="MGCAMB") m_Pk0_CAMB = 2.*pow(par::pi*m_sigma8,2)/Int*D_N*D_N;
  if (method_Pk=="MPTbreeze-v1") m_Pk0_MPTbreeze = 2.*pow(par::pi*m_sigma8,2)/Int*D_N*D_N;
  if (method_Pk=="CLASS") m_Pk0_CLASS = 2.*pow(par::pi*m_sigma8,2)/Int*D_N*D_N;

}


// =====================================================================================


std::vector<double> cbl::cosmology::Cosmology::Pk_matter (const std::vector<double> kk, const std::string method_Pk, const bool NL, const double redshift, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec, const std::string file_par, const bool unit1)
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

  if (Norm==1) Pk_0(method_Pk, redshift, store_output, output_root, k_min, k_max, prec, file_par);
  else { m_Pk0_EH = 1.; m_Pk0_CAMB = 1.; m_Pk0_MPTbreeze = 1.; m_Pk0_CLASS = 1.; }

  vector<double> Pk(kk.size());

  if (method_Pk=="EisensteinHu") { // NL is not used!!!

    EisensteinHu eh;

    eh.TFmdm_set_cosm(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massive_neutrinos, m_Omega_DE, m_hh, redshift, m_scalar_amp, m_scalar_pivot, m_n_spec);

    for (size_t i=0; i<kk.size(); i++) {
      Pk[i] = m_Pk0_EH*eh.Pk(newk[i])*fact2;
    }
  }

  else if (method_Pk=="CAMB" || method_Pk=="MGCAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="CLASS" || (method_Pk=="CAMB_wrapper" && kk.size()==1)) {

    vector<double> _kk, _pk;
    Table_PkCodes(method_Pk, NL, _kk, _pk, redshift, store_output, output_root, k_max, file_par);

    for(size_t i=0; i<_kk.size(); i++) {
      _kk[i] = pow(10., _kk[i]);
      _pk[i] = pow(10., _pk[i]);
    }
    
    glob::FuncGrid interp_Pk(_kk, _pk, "Spline");

    Pk = interp_Pk.eval_func(newk);
    
    double Pk0 = 1.;
    if (method_Pk=="CAMB") Pk0 = m_Pk0_CAMB;
    else if (method_Pk=="MGCAMB") Pk0 = m_Pk0_CAMB;
    else if (method_Pk=="MPTbreeze-v1") Pk0 = m_Pk0_MPTbreeze;
    else if (method_Pk=="CLASS") Pk0 = m_Pk0_CLASS;

    for (size_t i=0; i<kk.size(); i++)
      Pk[i] *= Pk0*fact2;
  }

  else if (method_Pk=="CAMB_wrapper" && kk.size()>1) {

    double fact3 = (m_unit || unit1) ? m_hh : 1.;

    Pk = wrapper::camb::Pk_CAMB (NL, redshift, Min(kk)*fact3, Max(kk)*fact3, (int)(kk.size()), m_Omega_baryon*m_hh*m_hh, m_Omega_CDM*m_hh*m_hh, m_Omega_neutrinos*m_hh*m_hh, m_massless_neutrinos, m_massive_neutrinos, m_Omega_k, m_hh*100., m_n_spec, m_scalar_amp, m_scalar_pivot, m_w0, m_wa, m_tau);

    for (size_t i=0; i<kk.size(); i++)
      Pk[i] *= m_Pk0_CAMB*fact2;
    
  }

  else { ErrorCBL("method_Pk is wrong!", "Pk", "PkXi.cpp");  vector<double> vv; return vv; }

  return Pk;
}


// =====================================================================================


std::vector<std::vector<double>> cbl::cosmology::Cosmology::Pk_matter (const std::vector<double> kk, const std::string method_Pk, const bool NL, const std::vector<double> redshift, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec, const std::string file_par, const bool unit1)
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

  vector<vector<double>> Pk(redshift.size(), vector<double>(kk.size()));
  double Pk0 = 1.;

  if (method_Pk=="EisensteinHu") { // NL is not used!!!

    for (size_t zz=0; zz<redshift.size(); zz++) {

      if (Norm==1) Pk_0(method_Pk, redshift[zz], store_output, output_root, k_min, k_max, prec, file_par);
      else { m_Pk0_EH = 1.;}

      EisensteinHu eh;

      eh.TFmdm_set_cosm(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massive_neutrinos, m_Omega_DE, m_hh, redshift[zz], m_scalar_amp, m_scalar_pivot, m_n_spec);

      for (size_t i=0; i<kk.size(); i++) {
	Pk[zz][i] = m_Pk0_EH*eh.Pk(newk[i])*fact2;
      }
    }
  }

  else if (method_Pk=="CAMB" || method_Pk=="MGCAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="CLASS") {

    vector<vector<double>> _kk, _pk;
    Table_PkCodes(method_Pk, NL, _kk, _pk, redshift, store_output, output_root, k_max, file_par);

    for (size_t zz=0; zz<redshift.size(); zz++) {

      for(size_t i=0; i<_kk[zz].size(); i++) {
	_kk[zz][i] = pow(10., _kk[zz][i]);
	_pk[zz][i] = pow(10., _pk[zz][i]);
      }

      glob::FuncGrid interp_Pk(_kk[zz], _pk[zz], "Spline");
    
      if (Norm==1) {

	double sigma8;

	if (NL==true)
	  sigma8 = sigma8_Pk("CAMB", 0., store_output, "test", false, k_min, k_max, prec);

	else {
	  const double RR = 8.;
	  auto func_sigma = [&] (double _k){
			      return pow(TopHat_WF(_k*RR)*_k, 2)*interp_Pk(_k);
			    };
	  sigma8 = sqrt(1./(2.*pow(par::pi, 2))*wrapper::gsl::GSL_integrate_qag (func_sigma, k_min, k_max, 1.e-5))/DN(redshift[zz], 0.);
	}

	Pk0 = pow(m_sigma8/sigma8,2);
      }

      Pk[zz] = interp_Pk.eval_func(newk);

      for (size_t i=0; i<kk.size(); i++)
	Pk[zz][i] *= Pk0*fact2;

    }
  }

  else { ErrorCBL("method_Pk is wrong!", "Pk", "PkXi.cpp");  vector<vector<double>> vv; return vv; }

  return Pk;
}


// =====================================================================================

/// @cond glob

double cbl::glob::func_xi_EH_GSL (double kk, void *params)
{
  struct glob::STR_xi_EH *pp = (struct glob::STR_xi_EH *) params;

  EisensteinHu eh;

  eh.TFmdm_set_cosm(pp->Omega_matter, pp->Omega_baryon, pp->Omega_neutrinos, pp->massive_neutrinos, pp->Omega_DE, pp->hh, pp->redshift, pp->scalar_amp, pp->scalar_pivot, pp->n_spec);

  double Int = eh.Pk(kk)*sin(kk*pp->rr)*kk/pp->rr;

  return Int*exp(-kk*kk*pp->aa*pp->aa); // eq. 24 of Anderson et al. 2012
}

/// @endcond


// =====================================================================================


double cbl::cosmology::Cosmology::xi_matter (const double rr, const std::string method_Pk, const bool NL, const double redshift, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const std::string file_par)
{
  bool gsl = GSL;
  if (gsl==false && method_Pk=="EisensteinHu") {
    //WarningMsgCBL("EisensteinHu method only works with GSL integration", "xi_matter", "PkXi.cpp");
    gsl = true;
  }

  int Norm = norm;
  if (Norm==-1) Norm = (m_sigma8>0) ? 1 : 0;

  if (method_Pk=="MPTbreeze-v1") Norm = 0; // check!!!

  double Int = -1.;
  double fact = (gsl) ? 1./(2.*pow(par::pi,2)) : 1.;

  if (gsl) {
    int limit_size = 1000;
    gsl_integration_workspace *ww = gsl_integration_workspace_alloc(limit_size);
    gsl_function Func;

    double error = -1.;

    if (method_Pk=="EisensteinHu") {

      if (m_sigma8<0) ErrorCBL("sigma8<0!", "xi_matter", "PkXi.cpp");
      if (NL==1) WarningMsgCBL("the correlation function by Eisenstein&Hu is linear (see xi_matter of PkXi.cpp)!", "xi_matter", "PkXi.cpp");

      glob::STR_xi_EH str;
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

    else if (method_Pk=="CAMB" || method_Pk=="MGCAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="CLASS") {
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

    else ErrorCBL("method_Pk is wrong!", "xi_matter", "PkXi.cpp");

    gsl_integration_workspace_free(ww);
  }


  else { // using FFTLOG
    if (method_Pk=="CAMB" || method_Pk=="MGCAMB" || method_Pk=="MPTbreeze-v1" || method_Pk=="CLASS") {
      vector<double> r, xi;
      Table_XiCodes(method_Pk, NL, r, xi, redshift, store_output, output_root, k_max, file_par);
      Int = interpolated(rr, r, xi, "Spline");
    }

    else ErrorCBL("method_Pk is wrong!", "xi_matter", "PkXi.cpp");

  }


  if (Norm==1) Pk_0(method_Pk, redshift, store_output, output_root, k_min, k_max, prec, file_par);

  double PP0 = -1.;
  if (method_Pk=="EisensteinHu") PP0 = m_Pk0_EH;
  if (method_Pk=="CAMB" || method_Pk=="MGCAMB") PP0 = m_Pk0_CAMB;
  if (method_Pk=="MPTbreeze-v1") PP0 = m_Pk0_MPTbreeze;
  if (method_Pk=="CLASS") PP0 = m_Pk0_CLASS;

  return PP0*fact*Int;
}


// =====================================================================================


double cbl::cosmology::Cosmology::wp_DM (const double rp, const std::string method_Pk, const bool NL, const double redshift, const double pimax, const bool store_output, const std::string output_root, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const std::string file_par)
{
  int Norm = norm;
  if (Norm==-1) Norm = (m_sigma8>0) ? 1 : 0;

  if (method_Pk=="MPTbreeze-v1") Norm = 0; // check!!!
  if (method_Pk=="MPTbreeze-v1" && NL==0) ErrorCBL("MPTbreeze is non-linear!", "wp_DM", "PkXi.cpp");

  // check if the table with lg(r)-lg(xi) already exists

  string mDir = (GSL==0) ? "fftlog" : "GSL";
  
  cbl::Path path;
  string dir_grid = path.DirCosmo()+"/Cosmology/Tables/"+mDir+"/"+method_Pk+"/h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6)+"/";

  //dir_grid += (output_root=="test") ? "/" : "_"+output_root;

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
      XI = xi_matter(rad[index], method_Pk, NL, redshift, store_output, output_root, 0, k_min, k_max, aa, GSL, prec, file_par);
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
  if (method_Pk=="CAMB" || method_Pk=="MGCAMB") PP0 = m_Pk0_CAMB;
  if (method_Pk=="MPTbreeze-v1") PP0 = m_Pk0_MPTbreeze;
  if (method_Pk=="CLASS") PP0 = m_Pk0_CLASS;

  return PP0*Int;
}


// =====================================================================================


double cbl::cosmology::Cosmology::sigmaR_DM (const double RR, const int corrType, const std::string method_Pk, const double redshift, const double pimax, const bool store_output, const std::string output_root, const bool NL, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const std::string file_par)
{
  // check if the table with lg(r)-lg(xi) already exists

  string mDir = "GSL";

  cbl::Path path;
  string dir_grid = path.DirCosmo()+"/Cosmology/Tables/"+mDir+"/"+method_Pk+"/h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6)+"/";

  //dir_grid += (output_root=="test") ? "/" : "_"+output_root;

  string file_table = (corrType==1) ? dir_grid+"xi_matter.dat" : dir_grid+"wp_DM.dat";
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
      XI = (corrType==1) ? xi_matter(rad[index], method_Pk, NL, redshift, store_output, output_root, norm, k_min, k_max, aa, GSL, prec, file_par) : wp_DM(rad[index], method_Pk, NL, redshift, pimax, store_output, output_root, norm, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
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
  struct glob::STR_sigma2M_EH *pp = (struct glob::STR_sigma2M_EH *) params;

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
  
  if (m_sigma8>0) return m_sigma8*DN(redshift);

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

      glob::STR_SSM str;
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

  classfunc::func_kstar func (m_hh, m_unit, lgkk, lgPk);

  function<double(double)> ff = bind(&classfunc::func_kstar::operator(), func, std::placeholders::_1);
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

  cbl::Path path;
  string dir_cosmo = path.DirCosmo();

  string dir_grid = dir_cosmo+"/Cosmology/Tables/"+mDir+"/"+nDir+"/h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6)+"/";

  //dir_grid += (output_root=="test") ? "/" : "_"+output_root;

  string file_table = (XiNL) ? dir_grid+"xi_matter.dat": dir_grid+"xi_matter_lin.dat";

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
      XI = (xiType==0) ? xi_matter(rad[index], method_Pk, XiNL, redshift, store_output, output_root, Norm, k_min, k_max, aa, GSL, prec, file_par) :
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

  cbl::Path path;
  string dir_cosmo = path.DirCosmo();

  string dir_grid = dir_cosmo+"/Cosmology/Tables/"+mDir+"/"+nDir+"/h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6);

  string file_table = (XiNL) ? dir_grid+"xi_matter.dat": dir_grid+"xi_matter_lin.dat";

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


vector<double> cbl::cosmology::Cosmology::Pk_matter_NoWiggles_gaussian (const vector<double> kk, const vector<double> PkLin, const vector<double> PkApprox, const double lambda, const string kind)
{
  vector<double> PkNW(kk.size());
  vector<double> OF(kk.size());

  for (size_t i=0; i<kk.size(); i++)
    OF[i] = PkLin[i]/PkApprox[i];

  // Get Pk normalization
  glob::FuncGrid interp_OF(kk, OF, "Spline");

  // Smooth the oscillatory part
  if (kind=="gaussian_3d")
    {
      double norm = sqrt(2/par::pi)/lambda*log(10.);
      for (size_t i=0; i<kk.size(); i++)
	{
	  auto integrand = [&] (const double log_q) {
			     double qq = pow(10, log_q);
			     double x = qq*kk[i]/(lambda*lambda);
			     double fact = -(qq*qq+kk[i]*kk[i])/(2*lambda*lambda)+gsl_sf_lnsinh(x)-log(kk[i]*qq)+3*log(qq)+log(interp_OF(qq));
			     return gsl_sf_exp(fact);
			   };

	  PkNW[i] = wrapper::gsl::GSL_integrate_cquad(integrand, -5, 3)*PkApprox[i]*norm;
	}
    }
  else if (kind=="gaussian_1d")
    {
      double norm = 1./(sqrt(2*par::pi)*lambda);
      for (size_t i=0; i<kk.size(); i++)
	{
	  double log_k = log10(kk[i]);

	  double log_qmin = log_k-4*lambda;
	  double log_qmax = log_k+4*lambda;

	  auto integrand = [&] (const double log_q) {
			     return interp_OF(pow(10., log_q))*exp(-pow(log_k-log_q, 2)/(2*lambda*lambda));
			   };

	  PkNW[i] = wrapper::gsl::GSL_integrate_cquad(integrand, log_qmin, log_qmax)*norm*PkApprox[i];
	}
    }
  else
    ErrorCBL("wrong name for gaussian wiggles smoothing!", "Pk_matter_NoWiggles_gaussian", "PkXi.cpp", glob::ExitCode::_error_);

  return PkNW;
}


// =====================================================================================


vector<double> cbl::cosmology::Cosmology::Pk_matter_NoWiggles_bspline (const vector<double> kk, const vector<double> PkLin, const vector<double> PkApprox, const int order, const int nknots)
{
  vector<double> log_kk(kk.size());
  vector<double> PkNW(kk.size());
  vector<double> OF(kk.size());

  for (size_t i=0; i<kk.size(); i++) {
    log_kk[i] = log10(kk[i]);
    OF[i] = PkLin[i]/PkApprox[i];
  }

  // Create b-spline
  glob::FuncGrid_Bspline bspline(log_kk, OF, nknots, order);
  vector<double> interp = bspline.eval_func(log_kk);
  for (size_t i=0; i<kk.size(); i++)
    PkNW[i] = interp[i]*PkApprox[i];

  return PkNW;
}


// =====================================================================================


vector<double> cbl::cosmology::Cosmology::Pk_matter_NoWiggles (const string method, const vector<double> kk, const double redshift, const string linear_method, const int order, const int nknots, const double lambda, const bool store_output, const std::string output_root, const bool norm, const double prec)
{
  vector<double> PkNW;
  if (method == "EisensteinHu") {
    PkNW = Pk_matter(kk, "EisensteinHu", false, redshift, store_output, output_root, norm, 1.e-4, 100., prec);
  }
  else if (method == "bspline" or method=="gaussian_1d" or method=="gaussian_3d"){
    vector<double> PkLin = Pk_matter(kk, linear_method, false, redshift, store_output, output_root, norm, 1.e-4, 100., prec);
    vector<double> PkApprox = Pk_matter(kk, "EisensteinHu", false, redshift, store_output, output_root, norm, 1.e-4, 100., prec);

    if (method == "bspline") {
      PkNW = Pk_matter_NoWiggles_bspline(kk, PkLin, PkApprox, order, nknots);
    }
    else if (method == "gaussian_3d" or method == "gaussian_1d") {
      PkNW = Pk_matter_NoWiggles_gaussian(kk, PkLin, PkApprox, lambda, method);
    }

    // Get Pk normalization
    glob::FuncGrid interp_camb(kk, PkLin, "Spline");
    double integral = interp_camb.integrate_qag(1.e-5, 1.e3);

    glob::FuncGrid interp_nw(kk, PkNW, "Spline");
    double interp_integral = interp_nw.integrate_qag(1.e-5, 1.e3);

    for (size_t i=0; i<kk.size(); i++)
      PkNW[i] = PkNW[i]*integral/interp_integral;
  }
  else
    ErrorCBL("wrong name for pk no wiggles!", "Pk_NoWiggles", "PkXi.cpp", glob::ExitCode::_error_);

  return PkNW;
}


// =====================================================================================


vector<double> cbl::cosmology::Cosmology::Pk_matter_Linear (const string method, const vector<double> kk, const double redshift, const bool store_output, const std::string output_root, const bool norm, const double prec)
{
  vector<double> pk;
  if (method=="CAMB" or method=="CLASS")
    pk = Pk_matter(kk, method, false, redshift, store_output, output_root, norm, 1.e-4, 100., prec);
  else
    ErrorCBL("wrong name for linear!", "Pk_matter_Linear", "PkXi.cpp", glob::ExitCode::_error_);

  return pk;
}


// =====================================================================================


vector<double> cbl::cosmology::Cosmology::Pk_matter_DeWiggled (const string linear_method, const string nowiggles_method, const vector<double> kk, const double redshift, const double sigma_NL, const int order, const int nknots, const double lambda, const bool store_output, const std::string output_root, const bool norm, const double prec)
{
  vector<double> PkLin = Pk_matter_Linear(linear_method, kk, redshift, store_output, output_root, norm, prec);
  vector<double> PkNW = Pk_matter_NoWiggles(nowiggles_method, kk, redshift, linear_method, order, nknots, lambda, store_output, output_root, norm, prec);
  vector<double> PkDW(kk.size(), 0);

  for (size_t i=0; i<kk.size(); i++)
    PkDW[i] = PkNW[i]+exp(0.5*pow(kk[i]*sigma_NL, 2))*(PkLin[i]-PkNW[i]);

  return PkDW;
}


// =====================================================================================


double cbl::cosmology::Cosmology::xi_matter_DeWiggle (const double rr, const double redshift, const double sigma_NL, const bool store_output, const std::string output_root, const bool norm, const double k_min, const double k_max, const double aa, const double prec)
{
  bool NL = false;

  string author1 = "CAMB";
  string author2 = "EisensteinHu";

  vector<double> kk, PkCamb, PkM, PkEH;
  Table_PkCodes(author1, NL, kk, PkCamb, redshift, store_output, output_root, k_max);

  for (size_t i = 0; i<kk.size(); i++)
    kk[i] = pow(10,kk[i]);

  PkCamb = Pk_matter(kk, author1, NL, redshift, store_output, output_root, norm, k_min, k_max, prec);
  PkEH = Pk_matter(kk, author2, NL, redshift, store_output, output_root, norm, k_min, k_max, prec);
  
  for (size_t i = 0; i<kk.size(); i++)
    PkM.push_back(PkEH[i]*(1+(PkCamb[i]/PkEH[i]-1)*exp(-0.5*pow(kk[i]*sigma_NL, 2))));

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
      glob::STR_covariance_XiMultipoles_integrand params;
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

  vector<vector<double>> sigma2 = sigma2_k(nn, Volume, kk, {Pk0, Pk2, Pk4}, {0,2,4});
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
      glob::STR_covariance_XiMultipoles_integrand params;
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

    glob::STR_XiMultipoles_integrand params;
    int limit_size = 1000;

    gsl_function Func;
    Func.function = &XiMultipoles_integrand;
    Func.params = &params;

    double k_min = 1.e-4;
    double k_max = 1.e0;
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

  double theta_rad = converted_angle (theta, coordUnits, CoordinateUnits::_radians_);

  vector<double> kk = logarithmic_bin_vector(200, k_min, k_max);
  vector<double> Pk;

  vector<double> dc;
  vector<double> _zz = linear_bin_vector(1000, 0., 2*Max(zz));

  for(size_t i=0; i<_zz.size(); i++)
    dc.push_back(this->D_C(_zz[i]));
  glob::FuncGrid DC_interp(_zz, dc, interpolationType);

  glob::Distribution phi(glob::DistributionType::_Interpolated_, zz, phiz, 0, interpolationType);

  double zmin, zmax, zmean;
  zmin = Min(zz);
  zmax = Max(zz);
  zmean = phi.mean();

  Pk = this->Pk_matter(kk, method_Pk, NL, zmean, store_output, output_root, norm, k_min, k_max, prec, file_par);

  vector<double> r, xi;
  wrapper::fftlog::transform_FFTlog(r, xi, 1, kk, Pk);
  glob::FuncGrid xi_interp(r, xi, interpolationType);

  if (GSL) {
    auto integrand = [&] (double z1)
		     {
		       double r1 = DC_interp(z1);

		       auto integrand_z2 = [&] (double z2) {
					     double r2 = DC_interp(z2);
					     double ss = sqrt(pow(r1,2)+pow(r2,2)-2*r1*r2*cos(theta_rad));
					     return xi_interp(ss)*phi(z2);
					   };

		       return wrapper::gsl::GSL_integrate_qag(integrand_z2, zmin, zmax)*phi(z1);
		     };
    return wrapper::gsl::GSL_integrate_qag(integrand, zmin, zmax);
  }
  else{
    auto integrand = [&] (vector<double> zz)
		     {
		       double r1 = DC_interp(zz[0]);
		       double r2 = DC_interp(zz[1]);
		       double ss = sqrt(pow(r1,2)+pow(r2,2)-2*r1*r2*cos(theta_rad));
		       return xi_interp(ss)*phi(zz[0])*phi(zz[1]);
		     };

    wrapper::cuba::CUBAwrapper integrator(integrand, 2);

    return integrator.IntegrateCuhre( {{zmin, zmax}, {zmin, zmax}});
  }
}


// =====================================================================================


double cbl::cosmology::Cosmology::wtheta_DM (const double theta, const std::vector<double> kk, const std::vector<double> Pk, const std::vector<double> zz, const std::vector<double> nz, const std::vector<double> phiz, const std::string interpolationType, const CoordinateUnits coordUnits, const bool GSL, const double redshift_Pk)
{
  const double theta_rad = converted_angle (theta, coordUnits, CoordinateUnits::_radians_);

  const double zmin = Min(zz);
  const double zmax = Max(zz);


  // set the distribution function

  vector<double> nphi(zz.size(), 0);

  for(size_t i=0; i<nphi.size(); i++)
    nphi[i] = (phiz.size()!=0) ? phiz[i]*nz[i] : nz[i];

  glob::Distribution phi(glob::DistributionType::_Interpolated_, zz, nphi, 0, interpolationType);
  //const double zmean = phi.mean();


  // set the distirbution integrand and normalization

  auto normalization_integrand = [&] (const double redshift)
				 {
				   return phi(redshift)*this->dV_dZdOmega(redshift, 1);
				 };

  double normalization = wrapper::gsl::GSL_integrate_qag(normalization_integrand, zmin, zmax);


  // compute the 2PCF model

  vector<double> r, xi;
  wrapper::fftlog::transform_FFTlog(r, xi, 1, kk, Pk);
  glob::FuncGrid xi_interp(r, xi, interpolationType);
  double integral;

  if (GSL) {
    auto integrand = [&] (double z1)
		     {
		       double r1 = this->D_C(z1);

		       auto integrand_z2 = [&] (double z2) {
					     double DD = this->DN((zz[0]+zz[1])*0.5, redshift_Pk);
					     double r2 = this->D_C(z2);
					     double ss = sqrt(pow(r1,2)+pow(r2,2)-2*r1*r2*cos(theta_rad));
					     return DD*DD*xi_interp(ss)*normalization_integrand(z2);
					   };

		       return wrapper::gsl::GSL_integrate_qag(integrand_z2, zmin, zmax)*normalization_integrand(z1);
		     };

    integral = wrapper::gsl::GSL_integrate_qag(integrand, zmin, zmax);

  }

  else {
    auto integrand = [&] (vector<double> zz)
		     {
		       double DD = this->DN((zz[0]+zz[1])*0.5, redshift_Pk);
		       double r1 = this->D_C(zz[0]);
		       double r2 = this->D_C(zz[1]);
		       double ss = sqrt(pow(r1,2)+pow(r2,2)-2*r1*r2*cos(theta_rad));
		       return DD*DD*xi_interp(ss)*normalization_integrand(zz[0])*normalization_integrand(zz[1]);
		     };

    wrapper::cuba::CUBAwrapper integrator(integrand, 2);

    integral = integrator.IntegrateCuhre( {{zmin, zmax}, {zmin, zmax}});
  }

  return integral/normalization;
}


// =====================================================================================


std::vector<double> cbl::cosmology::Cosmology::C_l_DM (const int lmax, const std::vector<double> zz, const std::vector<double> phiz, const std::string interpolationMethod, const std::string method_Pk, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const double prec, const std::string file_par)
{
  const double zmin = Min(zz);
  const double zmax = Max(zz);

  vector<double> kk = logarithmic_bin_vector(200, k_min, k_max);

  vector<double> Pk = this->Pk_matter(kk, method_Pk, false, 0., store_output, output_root, norm, k_min, k_max, prec, file_par);
  glob::FuncGrid Pk_interp(kk, Pk, interpolationMethod);

  auto integrand_sbao = [&] (const double kk)
			{
			  return Pk_interp(kk);
			};

  double sbao = sqrt(4.*par::pi*wrapper::gsl::GSL_integrate_qag(integrand_sbao, 1.e-4, 1)/3./pow(2*par::pi,3));

  glob::FuncGrid phi(zz, phiz, interpolationMethod);

  vector<double> C_l;

  ErrorCBL("Check the normalization of D(z)", "C_l_DM", "PkXi.cpp", glob::ExitCode::_workInProgress_);
  
  for (int l=0; l<lmax+1; l++) {
    double integral;

    if (l<60) {
      auto integrand = [&] ( const double kk)
		       {
			 auto integrand_z = [&] (const double zz)
					    {
					      return DN(zz)*jl(kk*D_C(zz), l);
					    };

			 double integral_z = wrapper::gsl::GSL_integrate_qag(integrand_z, zmin, zmax);
			 return kk*kk*Pk_interp(kk)*pow(integral_z, 2)*exp(-kk*kk*sbao*sbao);
		       };

      integral = 2*wrapper::gsl::GSL_integrate_qag(integrand, 1.e-4, 10)/par::pi;
    }
    else {

      auto integrand = [&] (const double zz)
		       {
			 double dc = D_C(zz);
			 double kk = (l+0.5)/dc;
			 return pow(DN(zz)*phi(zz), 2)*Pk_interp(kk)*HH(zz)/(par::cc*dc*dc);
		       };

      integral = wrapper::gsl::GSL_integrate_qag(integrand, zmin, zmax);
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

  vector<double> kk = logarithmic_bin_vector(200, k_min, k_max);
  vector<double> Pk;
  for(size_t i=0; i<kk.size(); i++)
  Pk.push_back( this->Pk_matter(kk[i], method_Pk, false, 0., output_root, norm, k_min, k_max, prec, file_par));

  glob::FuncGrid Pk_interp(kk, Pk, interpolationMethod);
  glob::FuncGrid phi(zz, phiz, interpolationMethod);

  vector<double> C_l;

  for (int l=0; l<lmax+1; l++) {
  auto integrand = [&] ( const double redshift)
  {
  double dc = D_C(redshift);
  double _kk = double(l)/dc;
  return pow(phi(redshift), 2)*HH(redshift)/pow(dc, 2)*pow(DN(redshift),2)*Pk_interp(_kk);
  };
  C_l.push_back(wrapper::gsl::GSL_integrate_qag(integrand, zmin, zmax)/par::cc);
  }

  return C_l;
  }
*/



