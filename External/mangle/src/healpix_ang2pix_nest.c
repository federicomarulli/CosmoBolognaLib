/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2005 Krzysztof M. Gorski, Eric Hivon, 
 *                          Benjamin D. Wandelt, Anthony J. Banday, 
 *                          Matthias Bartelmann, 
 *                          Reza Ansari & Kenneth M. Ganga 
 *
 *
 *  This file is part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.jpl.nasa.gov
 *
 *----------------------------------------------------------------------------- */
/* healpix_ang2pix_nest.c */

/* Standard Includes */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "pi.h"
#include "manglefn.h"

void healpix_ang2pix_nest( const int nside, long double theta, long double phi, int *ipix) {

  /* =======================================================================
   * subroutine ang2pix_nest(nside, theta, phi, ipix)
   * =======================================================================
   * gives the pixel number ipix (NESTED) corresponding to angles theta and phi
   *
   * the computation is made to the highest resolution available (nside=8192)
   * and then degraded to that required (by integer division)
   * this doesn't cost more, and it makes sure that the treatement of round-off 
   * will be consistent for every resolution
   * =======================================================================
   */
  
  long double z, za, z0, tt, tp, tmp;
  int    face_num,jp,jm;
  int    ifp, ifm;
  int    ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf, ntt;
  int    i, K, IP, I, J, id;
  int    ns_max = 8192;
  static int x2pix[128], y2pix[128];
  static char setup_done = 0;
  
  if( nside < 1 || nside > ns_max ) {
    fprintf(stderr, "healpix_ang2pix_nest: nside out of range: %d\n", nside);
    exit(0);
  }
  if( theta < 0. || theta > PI ) {
    fprintf(stderr, "healpix_ang2pix_nest: theta out of range: %Lf\n", theta);
    exit(0);
  }
  if( !setup_done ) {
    /* set the array giving the number of the pixel lying in (x,y), where x and y
       are in {1,128} and the pixel number is in {0,128^2-1} */
    for (i = 0; i < 127; i++) x2pix[i] = 0;
    for (I = 1; I <= 128; I++) {
        J = I - 1;       // !pixel numbers
	K = 0;
	IP = 1;
	truc : if (J == 0) {
	  x2pix[I-1] = K;
	  y2pix[I-1] = 2*K;
	}
	else {
	  id = (int)fmodl(J,2);
	  J  = J/2;
	  K  = IP*id + K;
	  IP = IP*4;
	  goto truc;
	}
    }
    setup_done = 1;
  }
  
  z  = cosl(theta);
  za = fabsl(z);
  z0 = 2./3.;
  if( phi >= TWOPI ) phi = phi - TWOPI;
  if( phi < 0. )     phi = phi + TWOPI;
  tt = phi / PIBYTWO; /* in [0,4] */
  
  if( za <= z0 ) { /* equatorial region */
    
    /* (the index of edge lines increase when the longitude=phi goes up) */
    jp = (int)floorl(ns_max*(0.5 + tt - z*0.75)); /* ascending edge line index */
    jm = (int)floorl(ns_max*(0.5 + tt + z*0.75)); /* descending edge line index */
    
    /* finds the face */
    ifp = jp / ns_max; /* in {0,4} */
    ifm = jm / ns_max;
    
    if( ifp == ifm ) face_num = (int)fmodl(ifp,4) + 4; /* faces 4 to 7 */
    else if( ifp < ifm ) face_num = (int)fmodl(ifp,4); /* (half-)faces 0 to 3 */
    else face_num = (int)fmodl(ifm,4) + 8;             /* (half-)faces 8 to 11 */
    
    ix = (int)fmodl(jm, ns_max);
    iy = ns_max - (int)fmodl(jp, ns_max) - 1;
  }
  else { /* polar region, za > 2/3 */
    
    ntt = (int)floorl(tt);
    if( ntt >= 4 ) ntt = 3;
    tp = tt - ntt;
    tmp = sqrtl( 3.*(1. - za) ); /* in ]0,1] */
    
    /* (the index of edge lines increase when distance from the closest pole
     * goes up)
     */
    /* line going toward the pole as phi increases */
    jp = (int)floorl( ns_max * tp          * tmp ); 

    /* that one goes away of the closest pole */
    jm = (int)floorl( ns_max * (1. - tp) * tmp );
    jp = (int)(jp < ns_max-1 ? jp : ns_max-1);
    jm = (int)(jm < ns_max-1 ? jm : ns_max-1);
    
    /* finds the face and pixel's (x,y) */
    if( z >= 0 ) {
      face_num = ntt; /* in {0,3} */
      ix = ns_max - jm - 1;
      iy = ns_max - jp - 1;
    }
    else {
      face_num = ntt + 8; /* in {8,11} */
      ix =  jp;
      iy =  jm;
    }
  }
  
  ix_low = (int)fmodl(ix,128);
  ix_hi  =     ix/128;
  iy_low = (int)fmodl(iy,128);
  iy_hi  =     iy/128;

  ipf = (x2pix[ix_hi] + y2pix[iy_hi]) * (128 * 128) + (x2pix[ix_low] + y2pix[iy_low]);
  ipf = (int)(ipf / powl(ns_max/nside,2));     /* in {0, nside**2 - 1} */
  *ipix =(int)( ipf + face_num*powl(nside,2)); /* in {0, 12*nside**2 - 1} */
}
