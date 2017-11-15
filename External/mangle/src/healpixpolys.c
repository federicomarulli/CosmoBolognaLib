/*-------------------------------------------------------------
(C) J C Hill 2006
-------------------------------------------------------------*/
#include <stdlib.h>
#include <math.h>
#include "pi.h"
#include "manglefn.h"

/*-----------------------------------------------------------
  get_healpix_poly: uses the HEALPix Fortran subroutine
                     pix2vec_nest to construct the HEALPix
                     pixels at the desired resolution

  Input: nside = HEALPix parameter describing resolution (for
                 res >= 1, simply defined by 2^(res-1); for
		 res = 0, define nside = 0)
         hpix = id number of HEALPix pixel that you wish to construct
  Return value: pointer to polygon if successful
                0x0 if error occurred
*/

polygon *get_healpix_poly(int nside, int hpix)
{
  int nv, nvmax, i, pix_n, pix_e, pix_s, pix_w, ev[1];
  vertices *vert;
  polygon *pixel, *pixelbetter;
  long double verts_vec[12], verts_vec_n[12], verts_vec_e[12], verts_vec_s[12], verts_vec_w[12], dist_n, dist_w, dist_s, dist_e;
  vec center, center_n, center_e, center_s, center_w, vertices_vec[4], vertices_vec_n[4], vertices_vec_e[4], vertices_vec_s[4], vertices_vec_w[4];
  azel *vertices_azel[8], vertices[8];

  for(i=0;i<=7;i++) vertices_azel[i] = &(vertices[i]);

  if (nside == 0) {
     pixel=new_poly(0);
     pixel->weight=1;
     pixel->pixel=0;
     
     return(pixel);
  }

  else {
    healpix_verts(nside, hpix, center, verts_vec);
    
    /* north vertex */
    for(i=0;i<=2;i++) (vertices_vec[0])[i] = verts_vec[i];

    /* east vertex */
    for(i=0;i<=2;i++) (vertices_vec[1])[i] = verts_vec[i+3];

    /* south vertex */
    for(i=0;i<=2;i++) (vertices_vec[2])[i] = verts_vec[i+6];

    /* west vertex */
    for(i=0;i<=2;i++) (vertices_vec[3])[i] = verts_vec[i+9];

    rp_to_azel(vertices_vec[0], vertices_azel[0]);
    rp_to_azel(vertices_vec[3], vertices_azel[2]);
    rp_to_azel(vertices_vec[2], vertices_azel[4]);
    rp_to_azel(vertices_vec[1], vertices_azel[6]);

    pix_n = 4*hpix + 3;
    pix_e = 4*hpix + 1;
    pix_s = 4*hpix;
    pix_w = 4*hpix + 2;

    healpix_verts(nside*2, pix_n, center_n, verts_vec_n);
    healpix_verts(nside*2, pix_e, center_e, verts_vec_e);
    healpix_verts(nside*2, pix_s, center_s, verts_vec_s);
    healpix_verts(nside*2, pix_w, center_w, verts_vec_w);

    /* north vertex of each child pixel */
    for(i=0;i<=2;i++){
      (vertices_vec_n[0])[i] = verts_vec_n[i];
      (vertices_vec_e[0])[i] = verts_vec_e[i];
      (vertices_vec_s[0])[i] = verts_vec_s[i];
      (vertices_vec_w[0])[i] = verts_vec_w[i];
    }

    /* east vertex of each child pixel */
    for(i=0;i<=2;i++){
      (vertices_vec_n[1])[i] = verts_vec_n[i+3];
      (vertices_vec_e[1])[i] = verts_vec_e[i+3];
      (vertices_vec_s[1])[i] = verts_vec_s[i+3];
      (vertices_vec_w[1])[i] = verts_vec_w[i+3];
    }

    /* south vertex of each child pixel */
    for(i=0;i<=2;i++){
      (vertices_vec_n[2])[i] = verts_vec_n[i+6];
      (vertices_vec_e[2])[i] = verts_vec_e[i+6];
      (vertices_vec_s[2])[i] = verts_vec_s[i+6];
      (vertices_vec_w[2])[i] = verts_vec_w[i+6];
    }

    /* west vertex of each child pixel */
    for(i=0;i<=2;i++){
      (vertices_vec_n[3])[i] = verts_vec_n[i+9];
      (vertices_vec_e[3])[i] = verts_vec_e[i+9];
      (vertices_vec_s[3])[i] = verts_vec_s[i+9];
      (vertices_vec_w[3])[i] = verts_vec_w[i+9];
    }
    
    rp_to_azel(vertices_vec_n[3], vertices_azel[1]);
    rp_to_azel(vertices_vec_w[2], vertices_azel[3]);
    rp_to_azel(vertices_vec_s[1], vertices_azel[5]);
    rp_to_azel(vertices_vec_e[0], vertices_azel[7]);

    for(i=0; i<8; i++){
      if(vertices[i].az < 0.) vertices[i].az = vertices[i].az + TWOPI;
      else {};
    }

    for(i=0; i<8; i++){
      if(vertices[i].az >= TWOPI) vertices[i].az = vertices[i].az - TWOPI;
      else {};
    } 

    nv=8; nvmax=8;
    vert=new_vert(nvmax);
    if(!vert){
      fprintf(stderr, "error in get_healpix_poly: failed to allocate memory for 8 vertices\n");
      return(0x0);
    }
    vert->nv=nv; vert->v=&vertices[0];

    pixel=new_poly(4);
    
    if(!pixel){
      fprintf(stderr, "error in get_healpix_poly: failed to allocate memory for polygon of 4 caps\n");
      return(0x0);
    }

    ev[0] = 8;

    edge_to_poly(vert, 2, &ev[0], pixel);
    pixel->id = hpix;

    pixelbetter=new_poly(5);
    
    if(!pixelbetter){
      fprintf(stderr, "error in get_healpix_poly: failed to allocate memory for polygon of 5 caps\n");
      return(0x0);
    }

    pixelbetter->np = 5;
    pixelbetter->npmax = 5;

    for(i=0; i<=3; i++) {
      pixelbetter->rp[i][0] = pixel->rp[i][0];
      pixelbetter->rp[i][1] = pixel->rp[i][1];
      pixelbetter->rp[i][2] = pixel->rp[i][2];
      pixelbetter->cm[i] = pixel->cm[i];
    }

    pixelbetter->rp[4][0] = center[0]; pixelbetter->rp[4][1] = center[1]; pixelbetter->rp[4][2] = center[2];
    pixelbetter->id = hpix;

    dist_n = cmrpirpj(center, vertices_vec[0]);
    dist_w = cmrpirpj(center, vertices_vec[3]);
    dist_s = cmrpirpj(center, vertices_vec[2]);
    dist_e = cmrpirpj(center, vertices_vec[1]);

    if(dist_n>=dist_w && dist_n>=dist_s && dist_n>=dist_e){
      pixelbetter->cm[4] = dist_n+0.000001;
    }
    else if(dist_w>=dist_n && dist_w>=dist_s && dist_w>=dist_e){
      pixelbetter->cm[4] = dist_w+0.000001;
    }
    else if(dist_s>=dist_n && dist_s>=dist_w && dist_s>=dist_e){
      pixelbetter->cm[4] = dist_s+0.000001;
    }
    else if(dist_e>=dist_n && dist_e>=dist_s && dist_e>=dist_w){
      pixelbetter->cm[4] = dist_e+0.000001;
    }
    else{
      fprintf(stderr, "error in get_healpix_poly: cannot find correct fifth cap\n");
      return(0x0);
    }


    if(!pixelbetter){
      fprintf(stderr, "error in get_healpix_poly: polygon is NULL.\n");
      return(0x0);
    }

    return(pixelbetter);
  }
}

/*------------------------------------------------------------
  get_nside: determines the HEALPix nside parameter (related to
             the resolution) from the number of weights (i.e.,
	     polygons) listed in the input file

  Input: nweights = number of weights in HEALPix_weight input file
  Return value: nside parameter
*/

int get_nside(int nweights)
{
  int res, nside;
  long double res_d;

  if (nweights == 1) {
     nside = 0;
     return(nside);
  }

  else {
    res_d = (long double)(logl(((long double)nweights)/3.0)/logl(4.0));

    /* res_d is often slightly under the correct res, so we add 0.1 to make it correct upon truncation */
    res = (int)(res_d+0.1);
    nside = (int)(powl(2,(int)(res-1)));
    return(nside);
  }
}

/*------------------------------------------------------------
  healpix_verts: interface to Fortran subroutine that calculates 
                the exact vertices of any HEALPix pixel at a
                given nside

  Input: nside = HEALPix parameter describing resolution
         pix = pixel number, where the first pixel number
               for any resolution is always 0
  Output: center = (x,y,z) coordinates of the pixel center
          verts = array of (x,y,z) positions of the four
	          pixel vertices in the order N,E,S,W
*/

void healpix_verts(int nside, int pix, vec center, long double verts[12])
{
  pix2vec_nest__(&nside, &pix, &(center[0]), &(center[1]), &(center[2]), &(verts[0]), &(verts[1]), &(verts[2]), &(verts[3]), &(verts[4]), &(verts[5]), &(verts[6]), &(verts[7]), &(verts[8]), &(verts[9]), &(verts[10]), &(verts[11]));
}

/*-------------------------------------------------------------
  cmrpirpj: C version of Fortran subroutine that calculates the value
         of (1-cosl(th(ij))), where th(ij) is the angle between the
	 unit vectors rpi and rpj

  Input: rpi, rpj = unit vectors
  Return value: 1- cosl(th(ij))
*/

long double cmrpirpj(vec rpi, vec rpj)
{
  long double cmij;
  
  cmij = (powl((rpi[0]-rpj[0]),2)+powl((rpi[1]-rpj[1]),2)+powl((rpi[2]-rpj[2]),2))/2.;
  return(cmij);
}
