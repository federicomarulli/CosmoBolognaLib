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
 *  @file CosmoBolognaLib/GlobalFunc/Func.cpp
 *
 *  @brief Generic functions that use one or more classes of the
 *  CosmoBolognaLib
 *
 *  This file contains the implementation of a set of generic
 *  functions that use one or more classes of the CosmoBolognaLib
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#include "GlobalFunc.h"

using namespace std;

using namespace cbl;


// ============================================================================


void cbl::redshift_range (const double mean_redshift, const double boxSide, cosmology::Cosmology &real_cosm, double &redshift_min, double &redshift_max) 
{
  coutCBL <<"I'm computing the redshift range..."<<endl; 

  double z_min = 0.;
  double lll = real_cosm.D_C(mean_redshift)+boxSide;
  double zf1 = mean_redshift, zf2 = zf1+12.;
  double z_max = real_cosm.Redshift(lll, zf1, zf2);
  int step = 50000;
  double delta_z = (z_max-z_min)/step;
  double zz1 = z_min, zz2, L1, L2, LL, dist, dist_min = 1.e20;
  
  for (int i=0; i<step; i++) {
    L1 = real_cosm.D_C(zz1);
    zz2 = 2.*mean_redshift-zz1;
    if (zz2<zz1) break;
    L2 = real_cosm.D_C(zz2);
    LL = L2-L1;
    dist = fabs(LL-boxSide);
    if (dist<dist_min) {
      dist_min = dist;
      redshift_min = zz1;
      redshift_max = zz2; 
    }			
    zz1 += delta_z;
  }	
  coutCBL <<"z1 = "<<redshift_min<<"; z2 = "<<redshift_max<<" (L_subBox = "<<real_cosm.D_C(redshift_max)-real_cosm.D_C(redshift_min)<<" ~ "<<boxSide<<")"<<endl;  
}


// ============================================================================


double cbl::volume (const double boxSize, const int frac, const double Bord, const double mean_redshift, cosmology::Cosmology &real_cosm)
{
  double redshift_min, redshift_max;
  double boxSide = boxSize/double(frac);
  redshift_range(mean_redshift, boxSide, real_cosm, redshift_min, redshift_max);
  redshift_min += Bord;
  redshift_max -= Bord;
  double Lmin = real_cosm.D_C(redshift_min);
  double Lmax = real_cosm.D_C(redshift_max);
  return pow(Lmax-Lmin, 3.);
}


// ============================================================================


void cbl::coord_zSpace (std::vector<double> &ra, std::vector<double> &dec, std::vector<double> &redshift, std::vector<double> &xx, std::vector<double> &yy, std::vector<double> &zz, const std::vector<double> vx, const std::vector<double> vy, const std::vector<double> vz, const double sigmaV, cosmology::Cosmology &real_cosm, const double mean_redshift, const double redshift_min, const double redshift_max, const int seed) 
{
  if (ra.size()==0 && xx.size()==0)
    ErrorCBL("both ra.size() and xx.size() are equal 0!", "coord_zSpace", "GlobalFunc/Func.cpp"); 
  
  vector<double> redshift_bin = linear_bin_vector(10000, redshift_min, redshift_max);
  vector<double> DC_bin(10000);

  for (int i=0; i<10000; i++)
    DC_bin[i] = real_cosm.D_C(redshift_bin[i]);

  glob::FuncGrid interp_DC(redshift_bin, DC_bin, "Spline");
  glob::FuncGrid interp_Z(DC_bin, redshift_bin, "Spline");
  
  if (ra.size()==0 && xx.size()>0) { 
    ra.resize(xx.size(), 0);
    dec.resize(xx.size(), 0);
    redshift.resize(xx.size(), 0);
    vector<double> dc(xx.size(), 0);

    cbl::polar_coord(xx, yy, zz, ra, dec, dc); 
    for (size_t i=0; i<dc.size(); ++i) 
      redshift[i] = interp_Z(dc[i]);
  }
  

  
  // ----- real-space --> redshift-space -----

  vector<double> dc;

  double catastrophic_error = sigmaV-int(sigmaV);
  double SigmaV = sigmaV-catastrophic_error;
  
  coutCBL <<"sigmaV = "<<SigmaV<<", catastrophic_error = "<<catastrophic_error<<endl;

  random::NormalRandomNumbers ran(0., SigmaV, seed);
  random::UniformRandomNumbers ran2(0., 1., seed);
  
  vector<double>::iterator pos;
  pos = min_element(xx.begin(),xx.end()); double xm = *pos;
  pos = max_element(xx.begin(),xx.end()); double xM = *pos; 
  pos = min_element(yy.begin(),yy.end()); double ym = *pos; 
  pos = max_element(yy.begin(),yy.end()); double yM = *pos; 
  pos = min_element(zz.begin(),zz.end()); double zm = *pos;
  pos = max_element(zz.begin(),zz.end()); double zM = *pos; 
 
  
  for (size_t i=0; i<ra.size(); i++) {
    
    if (ran2()<catastrophic_error) {  // catastrophic error     
      
      /*
      // test
      redshift[i] = ran2.doub()*(redshift_max-redshift_min)+redshift_min;
      dc.push_back(real_cosm.D_C(redshift[i]));    
      */
      
      // default
      double XX = ran2()*(xM-xm)+xm;
      double YY = ran2()*(yM-ym)+ym;
      double ZZ = ran2()*(zM-zm)+zm;
      double Dc = sqrt(XX*XX+YY*YY+ZZ*ZZ);
      dc.push_back(Dc);
      ra[i] = atan(XX/YY);
      dec[i] = asin(ZZ/Dc);
      double Zguess_min = 0.8*0.5, Zguess_max = 1.2*2.; 
      redshift[i] = real_cosm.Redshift(Dc, Zguess_min, Zguess_max);
      
    }

    else {

      double vrad = vx[i]*cos(dec[i])*sin(ra[i])+vy[i]*cos(dec[i])*cos(ra[i])+vz[i]*sin(dec[i]);
   
      double gerr = (SigmaV>0) ? ran()/par::cc : 0.;
      redshift[i] += vrad/par::cc*(1.+mean_redshift)+gerr; // peculiar velocities + gaussian error

      // dc.push_back(real_cosm.D_C(redshift[i]));          
      dc.push_back(interp_DC(redshift[i]));          
    
    }
  }

  cartesian_coord(ra, dec, dc, xx, yy, zz);

  coutCBL <<"done!"<<endl;
}


// ============================================================================


void cbl::create_mocks (const std::vector<double> xx, const std::vector<double> yy, const std::vector<double> zz,  const std::vector<double> vx, const std::vector<double> vy, const std::vector<double> vz, const std::vector<double> var1, const std::vector<double> var2, const std::vector<double> var3, const std::string output_dir, const double boxSize, const int frac, const double Bord, const double mean_redshift, cosmology::Cosmology &real_cosm, const int REAL, const double sigmaV, const int idum, double &Volume) 
{   
  coutCBL <<endl<<"I'm creating the mock files..."<<endl;

  
  //  ------- compute the redshift range -------  

  double redshift_min = -1., redshift_max = -1.;
  double boxSide = boxSize/double(frac);
  redshift_range(mean_redshift, boxSide, real_cosm, redshift_min, redshift_max);
  
  double Lmin = real_cosm.D_C(redshift_min);
  double Lmax = real_cosm.D_C(redshift_max);

  vector<double> shift;
  double SH = 0.;
  for (int i=0; i<frac; i++) {
    shift.push_back(SH);
    SH += boxSize/double(frac);
  }
  
  double fact1 = 1./(boxSize/double(frac)), fact2 = boxSize/double(frac)*0.5;
  vector<double> ra, dec, red, xx_temp, yy_temp, zz_temp;
  vector<int> subCube_temp, subCube;

  
  for (size_t i=0; i<xx.size(); i++) {
    
    double XX = xx[i];
    double YY = yy[i];
    double ZZ = zz[i];


    // ------- divide the box in sub-boxes ------- 

    int subx = min(int(XX*fact1), int(shift.size()-1));
    int suby = min(int(YY*fact1), int(shift.size()-1));
    int subz = min(int(ZZ*fact1), int(shift.size()-1));


    // ------- rescale the coordinates ------- 
    
    XX -= (shift[subx]+fact2);
    YY -= (shift[suby]-Lmin);
    ZZ -= (shift[subz]+fact2); 


    // ------- store the temporary vectors ------- 

    xx_temp.push_back(XX);
    yy_temp.push_back(YY);
    zz_temp.push_back(ZZ);
    subCube_temp.push_back(subx*frac*frac+suby*frac+subz);	

  }


  // ------- get polar coordinates ------- 
  
  vector<double> ra_temp(xx_temp.size()), dec_temp(xx_temp.size()), dd_temp(xx_temp.size()), red_temp(xx_temp.size());
  polar_coord (xx_temp, yy_temp, zz_temp, ra_temp, dec_temp, dd_temp);
 
  double Zguess_min = redshift_min*0.5, Zguess_max = redshift_max*2.;
  
  for (size_t i=0; i<dd_temp.size(); i++) red_temp[i] = real_cosm.Redshift(dd_temp[i], Zguess_min, Zguess_max);
  
  vector<double> avx, avy, avz;
  for (size_t i=0; i<vx.size(); i++) {avx.push_back(fabs(vx[i])); avy.push_back(fabs(vy[i])); avz.push_back(fabs(vz[i]));}
  coutCBL <<"<|v_x|> = "<<Average(avx)<<", <|v_y|> = "<<Average(avy)<<", <|v_z|> = "<<Average(avz)<<endl;


  // ------- add redshift-space distortions ------- 
  
  if (REAL==0) coord_zSpace(ra_temp, dec_temp, red_temp, xx_temp, yy_temp, zz_temp, vx, vy, vz, sigmaV, real_cosm, mean_redshift, redshift_min, redshift_max, idum);


  // ------- cut the borders of the box ------- 

  redshift_min += Bord;
  redshift_max -= Bord;
  Lmin = real_cosm.D_C(redshift_min);
  Lmax = real_cosm.D_C(redshift_max);
  double Lnew = (Lmax-Lmin)*0.5;
  Volume = pow(Lmax-Lmin, 3.);
  
  if (Lnew<0) ErrorCBL("", "create_mocks", "GlobalFunc/Func.cpp");
  
  coutCBL <<redshift_min<<" < z < "<<redshift_max<<" --> L = "<<Lnew*2.<<" --> Volume = "<<Volume<<endl;

  vector<double> vx_new, vy_new, vz_new, var1_new, var2_new, var3_new;
  
  for (size_t i=0; i<xx_temp.size(); i++)
    if (-Lnew<xx_temp[i] && xx_temp[i]<Lnew && Lmin<yy_temp[i] && yy_temp[i]<Lmax && -Lnew<zz_temp[i] && zz_temp[i]<Lnew) {
      ra.push_back(ra_temp[i]);
      dec.push_back(dec_temp[i]);
      red.push_back(red_temp[i]);
      vx_new.push_back(vx[i]);
      vy_new.push_back(vy[i]);
      vz_new.push_back(vz[i]);
      var1_new.push_back(var1[i]);
      var2_new.push_back(var2[i]);
      var3_new.push_back(var3[i]);
      subCube.push_back(subCube_temp[i]);	
    }

 
  // ------- store the data ------- 

  string MK = "mkdir -p "+output_dir;
  if (system(MK.c_str())) {};

  int nTOT = 0;
  
  for (int cub=0; cub<pow(frac,3.); cub++) {
	    
    string file_out = output_dir+"mock_"+conv(cub,par::fINT)+".dat";
    ofstream fout(file_out.c_str()); checkIO(fout, file_out);

    for (size_t i=0; i<ra.size(); i++) 
      if (subCube[i]==cub) {
	fout <<ra[i]<<"   "<<dec[i]<<"   "<<red[i]<<"   "<<vx_new[i]<<"   "<<vy_new[i]<<"   "<<vz_new[i]<<"   "<<var1_new[i]<<"   "<<var2_new[i]<<"   "<<var3_new[i]<<endl;
	nTOT ++;
      }
    fout.clear(); fout.close(); coutCBL <<"I wrote the file: "<<file_out<<endl;
  }

  coutCBL <<"Total number of haloes in the subcubes: "<<nTOT<<endl;
}
