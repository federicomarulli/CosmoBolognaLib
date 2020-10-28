/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli                         *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file CosmoBolognaLib/Kernel/Kernel.cpp
 *
 *  @brief Useful generic functions
 *
 *  This file contains the implementation of a large set of useful
 *  functions of wide use
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unibo.it
 */

#include "Kernel.h"

using namespace std;

using namespace cbl;
using namespace glob;


// ============================================================================


string cbl::fullpath (std::string path, const bool isDir)
{ 
  const string find = "~";
  const string replace = getenv("HOME");
  char buff[PATH_MAX];

  size_t pos = 0;
  while ((pos=path.find(find, pos))!=string::npos) {
    path.replace(pos, find.length(), replace);
    pos += replace.length();
  }

  return string(realpath(path.c_str(),buff))+((isDir) ? "/" : "");
}


// ============================================================================================


short cbl::ShortSwap (const short s)
{
  unsigned char b1, b2;
  b1 = s & 255;
  b2 = (s>>8) & 255;
  return (b1<<8) + b2;
}


// ============================================================================================


int cbl::IntSwap (const int i)
{
  unsigned char b1, b2, b3, b4;
  b1 = i & 255;
  b2 = (i>>8) & 255;
  b3 = (i>>16) & 255;
  b4 = (i>>24) & 255;
  return ((int)b1<<24) + ((int)b2<<16) + ((int)b3<<8) + b4;
}


// ============================================================================================


long long cbl::LongSwap (const long long i)
{
  unsigned char b1, b2, b3, b4, b5, b6, b7, b8;
  b1 = i & 255;
  b2 = (i>>8) & 255;
  b3 = (i>>16) & 255;
  b4 = (i>>24) & 255;
  b5 = (i>>32) & 255;
  b6 = (i>>40) & 255;
  b7 = (i>>48) & 255;
  b8 = (i>>56) & 255;
  return ((long long)b1<<56) + ((long long)b2<<48) + ((long long)b3<<40) + ((long long)b4<<32) + ((long long)b5<<24) + ((long long)b6<<16) + ((long long)b7<<8) + b8;
}


// ============================================================================================


float cbl::FloatSwap (const float f)
{
  union {
    float f;
    unsigned char b[4];
  } dat1, dat2;
  
  dat1.f = f;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
  
  return dat2.f;
}


// ============================================================================================


double cbl::DoubleSwap (const double d)
{
  union {
    double d;
    unsigned char b[8];
  } dat1, dat2;
  
  dat1.d = d;
  dat2.b[0] = dat1.b[7];
  dat2.b[1] = dat1.b[6];
  dat2.b[2] = dat1.b[5];
  dat2.b[3] = dat1.b[4];
  dat2.b[4] = dat1.b[3];
  dat2.b[5] = dat1.b[2];
  dat2.b[6] = dat1.b[1];
  dat2.b[7] = dat1.b[0];
  
  return dat2.d;
}


// ============================================================================================


double cbl::round_to_digits (const double num, const int ndigits)
{
  int exp_base10 = round(log10(num));
  double man_base10 = num*pow(10., -exp_base10);
  double factor = pow(10., -ndigits+1);  
  double truncated_man_base10 = man_base10-fmod(man_base10, factor);
  double rounded_remainder = fmod(man_base10, factor)/factor;

  rounded_remainder = rounded_remainder > 0.5 ? 1.0*factor : 0.0;

  return (truncated_man_base10+rounded_remainder)*pow(10.0, exp_base10);
}


// ============================================================================================


double cbl::round_to_precision (const double num, const int ndigits)
{
  const double tenth = pow(10., ndigits);
  return floor(num*tenth)/tenth;
}


// ============================================================================


void cbl::checkIO (const std::ifstream &fin, const std::string file)
{
  if (fin.fail()) {    
    string Err = "Error in opening the input file";
    if (file!="NULL") Err += ": " + file;
    ErrorCBL(Err, "checkIO", "Kernel.cpp", ExitCode::_IO_); 
  }
}


// ============================================================================


void cbl::checkIO (const std::ofstream &fout, const std::string file)
{
  if (fout.fail()) {    
    string Err = "Error in opening the output file";
    if (file!="NULL") Err += ": " + file;
    ErrorCBL(Err, "checkIO", "Kernel.cpp", ExitCode::_IO_); 
  }
}


// ============================================================================================


void cbl::set_EnvVar (std::vector<std::string> Var) 
{
  for (size_t vv=0; vv<Var.size(); vv++) 
    putenv(&Var[0][0]);
}


// ============================================================================================


void cbl::check_EnvVar (const std::string Var) 
{
  string COM = "if [ $"+Var+" ]; then touch tmp; fi";
  if (system (COM.c_str())) {};
  ifstream fin_check("tmp");
  if (!fin_check) 
    WarningMsgCBL("the variable " + Var + " has not been defined!", "check_EnvVar", "Kenrel.cpp");
  fin_check.clear(); fin_check.close();
  if (system("rm -f tmp")) {};
}


// ============================================================================================


int cbl::used_memory (const int type)
{
#ifdef LINUX
  int memory = -1;

  string mem;
  if (type==1) mem = "VmRSS:";
  else if (type==2) mem = "VmSize:";
  else ErrorCBL("the input value of type is not allowed!", "used_memory", "Kernel.cpp");
  
  string file = "/proc/self/status";
  ifstream fin(file.c_str()); checkIO(fin, file);
  
  string line, aa;
  while (getline(fin, line)) {
    stringstream ss(line);
    vector<string> val;
    while (ss>>aa) val.push_back(aa);
    if (val.size()==3 && val[0]==mem) {
      memory = atoi(val[1].c_str());
      break;
    }
  }
  
  fin.clear(); fin.close();
  
  return memory;

#else 
  (void)type;
  //WarningMsgCBL("this function works only on Linux systems", "used_memory", "Kernel.cpp");
  return 1;

#endif
  
}


// ============================================================================


int cbl::check_memory (const double frac, const bool exit, const std::string func, const int type)
{
#ifdef LINUX
  struct sysinfo memInfo;
  sysinfo (&memInfo);
  
  long long freePhysMem = memInfo.freeram;
  freePhysMem *= memInfo.mem_unit;
  int used = used_memory(type);
  
  if (used > freePhysMem*0.001*frac) { // 0.001 is to convert kbytes in bytes
    string Err = "possible memory problem";
    Err += (func.empty()) ? "!\n" : " in "+func+" !\n";
    Err += "freePhysMem = "+conv((double)freePhysMem*1.e-9, par::fDP3)+" GB\n";
    Err += "memory used by the process: = "+conv((double)used*1.e-6, par::fDP3)+" GB\n";
    if (exit) ErrorCBL(Err, "check_memory", "Kernel.cpp");
    else { WarningMsgCBL(Err, "check_memory", "Kernel.cpp"); return 0; }
  }
  return 1;
  
#else
  (void)frac; (void)exit; (void)func; (void)type;
  //WarningMsgCBL("this function works only on Linux systems", "checked_memory", "Kernel.cpp");
  return 1;
  
#endif
}


// ============================================================================


void cbl::unique_unsorted (std::vector<int> &vv) // erase all equal elements
{
  sort(vv.begin(),vv.end());
  vector<int>::iterator it;
  it = unique(vv.begin(),vv.end());
  vv.resize(it-vv.begin()); 
}


// ============================================================================


void cbl::unique_unsorted (std::vector<double> &vv) // erase all equal elements
{
  sort(vv.begin(),vv.end());
  vector<double>::iterator it;
  it = unique(vv.begin(),vv.end());
  vv.resize(it-vv.begin()); 
}


// ============================================================================


void cbl::unique_unsorted (std::vector<std::string> &vv) // erase all equal elements
{
  sort(vv.begin(),vv.end());
  vector<string>::iterator it;
  it = unique(vv.begin(),vv.end());
  vv.resize(it-vv.begin()); 
}


// ============================================================================

/// @cond glob
bool cbl::glob::operator<(const cbl::glob::CL &c1, const cbl::glob::CL &c2) {return c1.VV[0] < c2.VV[0];}
/// @endcond

void cbl::sort_2vectors (std::vector<double>::iterator p1, std::vector<double>::iterator p2, const int dim) 
{
  int temp = 0;
  vector<cbl::glob::CL> ccc; 
  for (int i=0; i<dim; i++) { 
    vector<double> vect = {*p1, *p2}; 
    cbl::glob::CL cl(vect); 
    ccc.push_back(cl);
    if (i+1<dim) {*p1 ++; *p2 ++; temp ++;}
  }
  sort (ccc.begin(),ccc.end());
  p1 -= temp; p2 -= temp;
  for (int i=0; i<dim; i++) { 
    *p1 = ccc[i].VV[0]; *p2 = ccc[i].VV[1];
    if (i+1<dim) {*p1 ++; *p2 ++;} 
  }
}

void cbl::sort_3vectors (std::vector<double>::iterator p1, std::vector<double>::iterator p2, std::vector<double>::iterator p3, const int dim) 
{
  int temp = 0;
  vector<cbl::glob::CL> ccc;
  for (int i=0; i<dim; i++) { 
    vector<double> vect = {*p1, *p2, *p3};
    cbl::glob::CL cl(vect); 
    ccc.push_back(cl);
    if (i+1<dim) {*p1 ++; *p2 ++; *p3 ++; temp ++;}
  }
  sort (ccc.begin(),ccc.end());
  p1 -= temp; p2 -= temp; p3 -= temp;
  for (int i=0; i<dim; i++) { 
    *p1 = ccc[i].VV[0]; *p2 = ccc[i].VV[1]; *p3 = ccc[i].VV[2];
    if (i+1<dim) {*p1 ++; *p2 ++; *p3 ++;} 
  }
}

void cbl::sort_4vectors (std::vector<double>::iterator p1, std::vector<double>::iterator p2, std::vector<double>::iterator p3, std::vector<double>::iterator p4, const int dim) 
{
  int temp = 0;
  vector<cbl::glob::CL> ccc;
  for (int i=0; i<dim; i++) { 
    vector<double> vect = {*p1, *p2, *p3, *p4};
    cbl::glob::CL cl(vect); 
    ccc.push_back(cl);
    if (i+1<dim) {*p1 ++; *p2 ++; *p3 ++; *p4 ++; temp ++;}
  }
  sort (ccc.begin(),ccc.end());
  p1 -= temp; p2 -= temp; p3 -= temp; p4 -= temp;
  for (int i=0; i<dim; i++) { 
    *p1 = ccc[i].VV[0]; *p2 = ccc[i].VV[1]; *p3 = ccc[i].VV[2]; *p4 = ccc[i].VV[3];
    if (i+1<dim) {*p1 ++; *p2 ++; *p3 ++; *p4 ++;} 
  }
}


// ============================================================================


int cbl::makeDir (std::string path, const std::string rootPath, const mode_t mode, const bool verbose)
{
  struct stat st;

  for (string::iterator iter=path.begin(); iter!=path.end();) {
    string::iterator newIter = find(iter, path.end(), '/');
    string newPath = rootPath+"/"+string(path.begin(), newIter);
    
    if (stat(newPath.c_str(), &st)!=0) {           
      if (mkdir(newPath.c_str(), mode)!=0 && errno!=EEXIST) 
	return ErrorCBL("cannot create the directory "+newPath+strerror(errno), "makeDir", "Kernel.cpp");
    }
    else
      if (!S_ISDIR(st.st_mode)) {
	errno = ENOTDIR;
	return ErrorCBL(newPath+" is not a directory", "makeDir", "Kernel.cpp");
      }
      else if (verbose)
	WarningMsgCBL(newPath+" already exists", "makeDir", "Kernel.cpp");

    iter = newIter;
    if (newIter!=path.end()) ++iter;
  }
  
  return 0;
}
