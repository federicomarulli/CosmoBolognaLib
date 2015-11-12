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
 *  @file GlobalFunc/Func.cpp
 *
 *  @brief Generic functions that use one or more classes of the
 *  CosmoBolognaLib
 *
 *  This file contains the implementation of a set of generic
 *  functions that use one or more classes of the CosmoBolognaLib
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#include "GlobalFunc.h"
using namespace cosmobl;


// ============================================================================


void cosmobl::redshift_range (double &mean_redshift, double &boxSide, Cosmology &real_cosm, double *redshift_min, double *redshift_max) 
{
  cout <<"I'm computing the redshift range..."<<endl; 

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
      *redshift_min = zz1;
      *redshift_max = zz2; 
    }			
    zz1 += delta_z;
  }	
  cout <<"z1 = "<<*redshift_min<<"; z2 = "<<*redshift_max<<" (L_subBox = "<<real_cosm.D_C(*redshift_max)-real_cosm.D_C(*redshift_min)<<" ~ "<<boxSide<<")"<<endl;  
}


// ============================================================================


double cosmobl::volume (double &boxSize, int &frac, double &Bord, double &mean_redshift, Cosmology &real_cosm)
{
  double redshift_min, redshift_max;
  double boxSide = boxSize/double(frac);
  redshift_range(mean_redshift, boxSide, real_cosm, &redshift_min, &redshift_max);
  redshift_min += Bord;
  redshift_max -= Bord;
  double Lmin = real_cosm.D_C(redshift_min);
  double Lmax = real_cosm.D_C(redshift_max);
  return pow(Lmax-Lmin, 3.);
}


// ============================================================================


void cosmobl::coord_zSpace (vector<double> &ra, vector<double> &dec, vector<double> &redshift, vector<double> &xx, vector<double> &yy, vector<double> &zz, vector<double> vx, vector<double> vy, vector<double> vz, double &sigmaV, Cosmology &real_cosm, double &mean_redshift, __attribute__((unused)) double &redshift_min, __attribute__((unused)) double &redshift_max, int &idum) 
{
  if (ra.size()==0) ErrorMsg("Error in coord_zSpace of GlobalFunc.cpp!");
  

  // ----- real-space --> redshift-space -----

  vector<double> dc;

  double catastrophic_error = sigmaV-int(sigmaV);
  double SigmaV = sigmaV-catastrophic_error;
  
  cout <<"sigmaV = "<<SigmaV<<", catastrophic_error = "<<catastrophic_error<<endl;

  Normaldev ran(0, SigmaV, idum);
  Ran ran2(idum);

  
  vector<double>::iterator pos;
  pos = min_element(xx.begin(),xx.end()); double xm = *pos;
  pos = max_element(xx.begin(),xx.end()); double xM = *pos; 
  pos = min_element(yy.begin(),yy.end()); double ym = *pos; 
  pos = max_element(yy.begin(),yy.end()); double yM = *pos; 
  pos = min_element(zz.begin(),zz.end()); double zm = *pos;
  pos = max_element(zz.begin(),zz.end()); double zM = *pos; 
 
  
  for (unsigned int i=0; i<ra.size(); i++) {
    
    if (ran2.doub()<catastrophic_error) {  // catastrophic error     
      
      /*
      // test
      redshift[i] = ran2.doub()*(redshift_max-redshift_min)+redshift_min;
      dc.push_back(real_cosm.D_C(redshift[i]));    
      */
      
      // default
      double XX = ran2.doub()*(xM-xm)+xm;
      double YY = ran2.doub()*(yM-ym)+ym;
      double ZZ = ran2.doub()*(zM-zm)+zm;
      double Dc = sqrt(XX*XX+YY*YY+ZZ*ZZ);
      dc.push_back(Dc);
      ra[i] = atan(XX/YY);
      dec[i] = asin(ZZ/Dc);
      double Zguess_min = 0.8*0.5, Zguess_max = 1.2*2.; 
      redshift[i] = real_cosm.Redshift(Dc, Zguess_min, Zguess_max);
      
    }

    else {

      double vrad = vx[i]*cos(dec[i])*sin(ra[i])+vy[i]*cos(dec[i])*cos(ra[i])+vz[i]*sin(dec[i]);
   
      redshift[i] += vrad/par::cc*(1.+mean_redshift) + ran.dev()/par::cc; // peculiar velocities + gaussian error

      dc.push_back(real_cosm.D_C(redshift[i]));          
    
    }
  }

  cartesian_coord(ra, dec, dc, xx, yy, zz);

  cout <<"done!"<<endl;
}


// ============================================================================


void cosmobl::create_mocks (vector<double> xx, vector<double> yy, vector<double> zz, vector<double> vx, vector<double> vy, vector<double> vz, vector<double> var1, vector<double> var2, vector<double> var3, string &output_dir, double &boxSize, int &frac, double &Bord, double &mean_redshift, Cosmology &real_cosm, int &REAL, double &sigmaV, int &idum, double *Volume) 
{   
  cout <<endl<<"I'm creating the mock files..."<<endl;

  
  //  ------- compute the redshift range -------  

  double redshift_min = -1., redshift_max = -1.;
  double boxSide = boxSize/double(frac);
  redshift_range(mean_redshift, boxSide, real_cosm, &redshift_min, &redshift_max);
  
  double Lmin = real_cosm.D_C(redshift_min);
  double Lmax = real_cosm.D_C(redshift_max);

  vector<double> shift;
  double SH = 0.;
  for (int i=0; i<frac; i++) {
    shift.push_back(SH);
    SH += boxSize/double(frac);
  }
  
  double fact1 = 1./(boxSize/double(frac)), fact2 = boxSize/double(frac)*0.5;
  vector<double> ra, dec, red, 
    xx_temp, yy_temp, zz_temp, vx_temp, vy_temp, vz_temp, vr_temp, var1_temp, var2_temp, var3_temp;
  vector<int> subCube_temp, subCube;

  
  for (unsigned int i=0; i<xx.size(); i++) {
    
    double XX = xx[i];
    double YY = yy[i];
    double ZZ = zz[i];
    double Vx = vx[i];
    double Vy = vy[i];
    double Vz = vz[i];


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
    vx_temp.push_back(Vx);
    vy_temp.push_back(Vy);
    vz_temp.push_back(Vz);
    var1_temp.push_back(var1[i]);
    var2_temp.push_back(var2[i]);
    var3_temp.push_back(var3[i]);
    subCube_temp.push_back(subx*frac*frac+suby*frac+subz);	

  }


  // ------- get polar coordinates ------- 
  
  vector<double> ra_temp(xx_temp.size()), dec_temp(xx_temp.size()), dd_temp(xx_temp.size()), red_temp(xx_temp.size());
  polar_coord (xx_temp, yy_temp, zz_temp, ra_temp, dec_temp, dd_temp);
 
  double Zguess_min = redshift_min*0.5, Zguess_max = redshift_max*2.;
  
  for (unsigned int i=0; i<dd_temp.size(); i++) red_temp[i] = real_cosm.Redshift(dd_temp[i], Zguess_min, Zguess_max);
  
  vector<double> avx, avy, avz;
  for (unsigned int i=0; i<vx.size(); i++) {avx.push_back(fabs(vx[i])); avy.push_back(fabs(vy[i])); avz.push_back(fabs(vz[i]));}
  cout <<"<|v_x|> = "<<Average(avx)<<", <|v_y|> = "<<Average(avy)<<", <|v_z|> = "<<Average(avz)<<endl;


  // ------- add redshift-space distortions ------- 

  if (REAL==0) coord_zSpace(ra_temp, dec_temp, red_temp, xx_temp, yy_temp, zz_temp, vx_temp, vy_temp, vz_temp, sigmaV, real_cosm, mean_redshift, redshift_min, redshift_max, idum);


  // ------- cut the borders of the box ------- 

  redshift_min += Bord;
  redshift_max -= Bord;
  Lmin = real_cosm.D_C(redshift_min);
  Lmax = real_cosm.D_C(redshift_max);
  double Lnew = (Lmax-Lmin)*0.5;
  *Volume = pow(Lmax-Lmin,3.);
  if (Lnew<0) ErrorMsg("Error in create_mocks of GlobalFunc.h!");
  cout <<redshift_min<<" < z < "<<redshift_max<<" --> L = "<<Lnew*2.<<" --> Volume = "<<*Volume<<endl;

  for (unsigned int i=0; i<xx_temp.size(); i++)
    if (-Lnew<xx_temp[i] && xx_temp[i]<Lnew && Lmin<yy_temp[i] && yy_temp[i]<Lmax && -Lnew<zz_temp[i] && zz_temp[i]<Lnew) {
      ra.push_back(ra_temp[i]);
      dec.push_back(dec_temp[i]);
      red.push_back(red_temp[i]);
      vx.push_back(vx_temp[i]);
      vy.push_back(vy_temp[i]);
      vz.push_back(vz_temp[i]);
      var1.push_back(var1_temp[i]);
      var2.push_back(var2_temp[i]);
      var3.push_back(var3_temp[i]);
      subCube.push_back(subCube_temp[i]);	
    }

 
  // ------- store the data ------- 

  string MK = "mkdir -p "+output_dir;
  if (system(MK.c_str())) {};

  int nTOT = 0;
  
  for (int cub=0; cub<pow(frac,3.); cub++) {
	    
    string file_out = output_dir+"mock_"+conv(cub,par::fINT)+".dat";
    ofstream fout (file_out.c_str()); checkIO (file_out,0);

    for (unsigned int i=0; i<ra.size(); i++) 
      if (subCube[i]==cub) {
	fout <<ra[i]<<"   "<<dec[i]<<"   "<<red[i]<<"   "<<vx[i]<<"   "<<vy[i]<<"   "<<vz[i]<<"   "<<var1[i]<<"   "<<var2[i]<<"   "<<var3[i]<<endl;
	nTOT ++;
      }
    fout.clear(); fout.close(); cout <<"I wrote the file: "<<file_out<<endl;
  }

  cout <<"Total number of haloes in the subcubes: "<<nTOT<<endl;
}
