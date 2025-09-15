/*
 *  This file is part of libsharp.
 *
 *  libsharp is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libsharp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libsharp; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/* libsharp is being developed at the Max-Planck-Institut fuer Astrophysik */

/*! \file sharp_core_inc.c
 *  Computational core
 *
 *  Copyright (C) 2012-2019 Max-Planck-Society
 *  \author Martin Reinecke
 */

// FIXME: special ugly workaround for problems on OSX
#if (!defined(__APPLE__)) || (!defined(__AVX512F__))

#if (defined(MULTIARCH) || defined(GENERIC_ARCH))

#define XCONCATX(a,b) a##_##b
#define XCONCATX2(a,b) XCONCATX(a,b)
#define XARCH(a) XCONCATX2(a,ARCH)

#include <complex.h>
#include <math.h>
#include <string.h>
#include "libsharp/sharp_vecsupport.h"
#include "libsharp/sharp.h"
#include "libsharp/sharp_internal.h"
#include "c_utils/c_utils.h"

typedef complex double dcmplx;

#define nv0 (128/VLEN)
#define nvx (64/VLEN)

typedef Tv Tbv0[nv0];
typedef double Tbs0[nv0*VLEN];

typedef struct
  {
  Tbv0 sth, corfac, scale, lam1, lam2, csq, p1r, p1i, p2r, p2i;
  } s0data_v;

typedef struct
  {
  Tbs0 sth, corfac, scale, lam1, lam2, csq, p1r, p1i, p2r, p2i;
  } s0data_s;

typedef union
  {
  s0data_v v;
  s0data_s s;
  } s0data_u;

typedef Tv Tbvx[nvx];
typedef double Tbsx[nvx*VLEN];

typedef struct
  {
  Tbvx sth, cfp, cfm, scp, scm, l1p, l2p, l1m, l2m, cth,
       p1pr, p1pi, p2pr, p2pi, p1mr, p1mi, p2mr, p2mi;
  } sxdata_v;

typedef struct
  {
  Tbsx sth, cfp, cfm, scp, scm, l1p, l2p, l1m, l2m, cth,
       p1pr, p1pi, p2pr, p2pi, p1mr, p1mi, p2mr, p2mi;
  } sxdata_s;

typedef union
  {
  sxdata_v v;
  sxdata_s s;
  } sxdata_u;

static inline void Tvnormalize (Tv * restrict val, Tv * restrict scale,
  double maxval)
  {
  const Tv vfmin=vload(sharp_fsmall*maxval), vfmax=vload(maxval);
  const Tv vfsmall=vload(sharp_fsmall), vfbig=vload(sharp_fbig);
  Tm mask = vgt(vabs(*val),vfmax);
  while (vanyTrue(mask))
    {
    vmuleq_mask(mask,*val,vfsmall);
    vaddeq_mask(mask,*scale,vone);
    mask = vgt(vabs(*val),vfmax);
    }
  mask = vand_mask(vlt(vabs(*val),vfmin),vne(*val,vzero));
  while (vanyTrue(mask))
    {
    vmuleq_mask(mask,*val,vfbig);
    vsubeq_mask(mask,*scale,vone);
    mask = vand_mask(vlt(vabs(*val),vfmin),vne(*val,vzero));
    }
  }

static void mypow(Tv val, int npow, const double * restrict powlimit,
  Tv * restrict resd, Tv * restrict ress)
  {
  Tv vminv=vload(powlimit[npow]);
  Tm mask = vlt(vabs(val),vminv);
  if (!vanyTrue(mask)) // no underflows possible, use quick algoritm
    {
    Tv res=vone;
    do
      {
      if (npow&1)
        res*=val;
      val*=val;
      }
    while(npow>>=1);
    *resd=res;
    *ress=vzero;
    }
  else
    {
    Tv scale=vzero, scaleint=vzero, res=vone;
    Tvnormalize(&val,&scaleint,sharp_fbighalf);
    do
      {
      if (npow&1)
        {
        res*=val;
        scale+=scaleint;
        Tvnormalize(&res,&scale,sharp_fbighalf);
        }
      val*=val;
      scaleint+=scaleint;
      Tvnormalize(&val,&scaleint,sharp_fbighalf);
      }
    while(npow>>=1);
    *resd=res;
    *ress=scale;
    }
  }

static inline void getCorfac(Tv scale, Tv * restrict corfac,
  const double * restrict cf)
  {
  typedef union
    { Tv v; double s[VLEN]; } Tvu;

  Tvu sc, corf;
  sc.v=scale;
  for (int i=0; i<VLEN; ++i)
    corf.s[i] = (sc.s[i]<sharp_minscale) ?
      0. : cf[(int)(sc.s[i])-sharp_minscale];
  *corfac=corf.v;
  }

static inline int rescale(Tv * restrict v1, Tv * restrict v2, Tv * restrict s, Tv eps)
  {
  Tm mask = vgt(vabs(*v2),eps);
  if (vanyTrue(mask))
    {
    vmuleq_mask(mask,*v1,vload(sharp_fsmall));
    vmuleq_mask(mask,*v2,vload(sharp_fsmall));
    vaddeq_mask(mask,*s,vone);
    return 1;
    }
  return 0;
  }

NOINLINE static void iter_to_ieee(const sharp_Ylmgen_C * restrict gen,
  s0data_v * restrict d, int * restrict l_, int * restrict il_, int nv2)
  {
  int l=gen->m, il=0;
  Tv mfac = vload((gen->m&1) ? -gen->mfac[gen->m]:gen->mfac[gen->m]);
  Tv limscale=vload(sharp_limscale);
  int below_limit = 1;
  for (int i=0; i<nv2; ++i)
    {
    d->lam1[i]=vzero;
    mypow(d->sth[i],gen->m,gen->powlimit,&d->lam2[i],&d->scale[i]);
    d->lam2[i] *= mfac;
    Tvnormalize(&d->lam2[i],&d->scale[i],sharp_ftol);
    below_limit &= vallTrue(vlt(d->scale[i],limscale));
    }

  while (below_limit)
    {
    if (l+4>gen->lmax) {*l_=gen->lmax+1;return;}
    below_limit=1;
    Tv a1=vload(gen->coef[il  ].a), b1=vload(gen->coef[il  ].b);
    Tv a2=vload(gen->coef[il+1].a), b2=vload(gen->coef[il+1].b);
    for (int i=0; i<nv2; ++i)
      {
      d->lam1[i] = (a1*d->csq[i] + b1)*d->lam2[i] + d->lam1[i];
      d->lam2[i] = (a2*d->csq[i] + b2)*d->lam1[i] + d->lam2[i];
      if (rescale(&d->lam1[i], &d->lam2[i], &d->scale[i], vload(sharp_ftol)))
        below_limit &= vallTrue(vlt(d->scale[i],vload(sharp_limscale)));
      }
    l+=4; il+=2;
    }
  *l_=l; *il_=il;
  }

NOINLINE static void alm2map_kernel(s0data_v * restrict d,
  const sharp_ylmgen_dbl2 * restrict coef, const dcmplx * restrict alm,
  int l, int il, int lmax, int nv2)
  {
  if (nv2==nv0)
    {
    for (; l<=lmax-2; il+=2, l+=4)
      {
      Tv ar1=vload(creal(alm[l  ])), ai1=vload(cimag(alm[l  ]));
      Tv ar2=vload(creal(alm[l+1])), ai2=vload(cimag(alm[l+1]));
      Tv ar3=vload(creal(alm[l+2])), ai3=vload(cimag(alm[l+2]));
      Tv ar4=vload(creal(alm[l+3])), ai4=vload(cimag(alm[l+3]));
      Tv a1=vload(coef[il  ].a), b1=vload(coef[il  ].b);
      Tv a2=vload(coef[il+1].a), b2=vload(coef[il+1].b);
      for (int i=0; i<nv0; ++i)
        {
        d->p1r[i] += d->lam2[i]*ar1;
        d->p1i[i] += d->lam2[i]*ai1;
        d->p2r[i] += d->lam2[i]*ar2;
        d->p2i[i] += d->lam2[i]*ai2;
        d->lam1[i] = (a1*d->csq[i] + b1)*d->lam2[i] + d->lam1[i];
        d->p1r[i] += d->lam1[i]*ar3;
        d->p1i[i] += d->lam1[i]*ai3;
        d->p2r[i] += d->lam1[i]*ar4;
        d->p2i[i] += d->lam1[i]*ai4;
        d->lam2[i] = (a2*d->csq[i] + b2)*d->lam1[i] + d->lam2[i];
        }
      }
    }
  else
    {
    for (; l<=lmax-2; il+=2, l+=4)
      {
      Tv ar1=vload(creal(alm[l  ])), ai1=vload(cimag(alm[l  ]));
      Tv ar2=vload(creal(alm[l+1])), ai2=vload(cimag(alm[l+1]));
      Tv ar3=vload(creal(alm[l+2])), ai3=vload(cimag(alm[l+2]));
      Tv ar4=vload(creal(alm[l+3])), ai4=vload(cimag(alm[l+3]));
      Tv a1=vload(coef[il  ].a), b1=vload(coef[il  ].b);
      Tv a2=vload(coef[il+1].a), b2=vload(coef[il+1].b);
      for (int i=0; i<nv2; ++i)
        {
        d->p1r[i] += d->lam2[i]*ar1;
        d->p1i[i] += d->lam2[i]*ai1;
        d->p2r[i] += d->lam2[i]*ar2;
        d->p2i[i] += d->lam2[i]*ai2;
        d->lam1[i] = (a1*d->csq[i] + b1)*d->lam2[i] + d->lam1[i];
        d->p1r[i] += d->lam1[i]*ar3;
        d->p1i[i] += d->lam1[i]*ai3;
        d->p2r[i] += d->lam1[i]*ar4;
        d->p2i[i] += d->lam1[i]*ai4;
        d->lam2[i] = (a2*d->csq[i] + b2)*d->lam1[i] + d->lam2[i];
        }
      }
    }
  for (; l<=lmax; ++il, l+=2)
    {
    Tv ar1=vload(creal(alm[l  ])), ai1=vload(cimag(alm[l  ]));
    Tv ar2=vload(creal(alm[l+1])), ai2=vload(cimag(alm[l+1]));
    Tv a=vload(coef[il].a), b=vload(coef[il].b);
    for (int i=0; i<nv2; ++i)
      {
      d->p1r[i] += d->lam2[i]*ar1;
      d->p1i[i] += d->lam2[i]*ai1;
      d->p2r[i] += d->lam2[i]*ar2;
      d->p2i[i] += d->lam2[i]*ai2;
      Tv tmp = (a*d->csq[i] + b)*d->lam2[i] + d->lam1[i];
      d->lam1[i] = d->lam2[i];
      d->lam2[i] = tmp;
      }
    }
  }

NOINLINE static void calc_alm2map (sharp_job * restrict job,
  const sharp_Ylmgen_C * restrict gen, s0data_v * restrict d, int nth)
  {
  int l,il,lmax=gen->lmax;
  int nv2 = (nth+VLEN-1)/VLEN;
  iter_to_ieee(gen, d, &l, &il, nv2);
  job->opcnt += il * 4*nth;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 6*nth;

  const sharp_ylmgen_dbl2 * restrict coef = gen->coef;
  const dcmplx * restrict alm=job->almtmp;
  int full_ieee=1;
  for (int i=0; i<nv2; ++i)
    {
    getCorfac(d->scale[i], &d->corfac[i], gen->cf);
    full_ieee &= vallTrue(vge(d->scale[i],vload(sharp_minscale)));
    }

  while((!full_ieee) && (l<=lmax))
    {
    Tv ar1=vload(creal(alm[l  ])), ai1=vload(cimag(alm[l  ]));
    Tv ar2=vload(creal(alm[l+1])), ai2=vload(cimag(alm[l+1]));
    Tv a=vload(coef[il].a), b=vload(coef[il].b);
    full_ieee=1;
    for (int i=0; i<nv2; ++i)
      {
      d->p1r[i] += d->lam2[i]*d->corfac[i]*ar1;
      d->p1i[i] += d->lam2[i]*d->corfac[i]*ai1;
      d->p2r[i] += d->lam2[i]*d->corfac[i]*ar2;
      d->p2i[i] += d->lam2[i]*d->corfac[i]*ai2;
      Tv tmp = (a*d->csq[i] + b)*d->lam2[i] + d->lam1[i];
      d->lam1[i] = d->lam2[i];
      d->lam2[i] = tmp;
      if (rescale(&d->lam1[i], &d->lam2[i], &d->scale[i], vload(sharp_ftol)))
        getCorfac(d->scale[i], &d->corfac[i], gen->cf);
      full_ieee &= vallTrue(vge(d->scale[i],vload(sharp_minscale)));
      }
    l+=2; ++il;
    }
  if (l>lmax) return;

  for (int i=0; i<nv2; ++i)
    {
    d->lam1[i] *= d->corfac[i];
    d->lam2[i] *= d->corfac[i];
    }
  alm2map_kernel(d, coef, alm, l, il, lmax, nv2);
  }

NOINLINE static void map2alm_kernel(s0data_v * restrict d,
  const sharp_ylmgen_dbl2 * restrict coef, dcmplx * restrict alm, int l,
  int il, int lmax, int nv2)
  {
  for (; l<=lmax-2; il+=2, l+=4)
    {
    Tv a1=vload(coef[il  ].a), b1=vload(coef[il  ].b);
    Tv a2=vload(coef[il+1].a), b2=vload(coef[il+1].b);
    Tv atmp1[4] = {vzero, vzero, vzero, vzero};
    Tv atmp2[4] = {vzero, vzero, vzero, vzero};
    for (int i=0; i<nv2; ++i)
      {
      atmp1[0] += d->lam2[i]*d->p1r[i];
      atmp1[1] += d->lam2[i]*d->p1i[i];
      atmp1[2] += d->lam2[i]*d->p2r[i];
      atmp1[3] += d->lam2[i]*d->p2i[i];
      d->lam1[i] = (a1*d->csq[i] + b1)*d->lam2[i] + d->lam1[i];
      atmp2[0] += d->lam1[i]*d->p1r[i];
      atmp2[1] += d->lam1[i]*d->p1i[i];
      atmp2[2] += d->lam1[i]*d->p2r[i];
      atmp2[3] += d->lam1[i]*d->p2i[i];
      d->lam2[i] = (a2*d->csq[i] + b2)*d->lam1[i] + d->lam2[i];
      }
    vhsum_cmplx_special (atmp1[0], atmp1[1], atmp1[2], atmp1[3], &alm[l  ]);
    vhsum_cmplx_special (atmp2[0], atmp2[1], atmp2[2], atmp2[3], &alm[l+2]);
    }
  for (; l<=lmax; ++il, l+=2)
    {
    Tv a=vload(coef[il].a), b=vload(coef[il].b);
    Tv atmp[4] = {vzero, vzero, vzero, vzero};
    for (int i=0; i<nv2; ++i)
      {
      atmp[0] += d->lam2[i]*d->p1r[i];
      atmp[1] += d->lam2[i]*d->p1i[i];
      atmp[2] += d->lam2[i]*d->p2r[i];
      atmp[3] += d->lam2[i]*d->p2i[i];
      Tv tmp = (a*d->csq[i] + b)*d->lam2[i] + d->lam1[i];
      d->lam1[i] = d->lam2[i];
      d->lam2[i] = tmp;
      }
    vhsum_cmplx_special (atmp[0], atmp[1], atmp[2], atmp[3], &alm[l]);
    }
  }

NOINLINE static void calc_map2alm (sharp_job * restrict job,
  const sharp_Ylmgen_C * restrict gen, s0data_v * restrict d, int nth)
  {
  int l,il,lmax=gen->lmax;
  int nv2 = (nth+VLEN-1)/VLEN;
  iter_to_ieee(gen, d, &l, &il, nv2);
  job->opcnt += il * 4*nth;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 6*nth;

  const sharp_ylmgen_dbl2 * restrict coef = gen->coef;
  dcmplx * restrict alm=job->almtmp;
  int full_ieee=1;
  for (int i=0; i<nv2; ++i)
    {
    getCorfac(d->scale[i], &d->corfac[i], gen->cf);
    full_ieee &= vallTrue(vge(d->scale[i],vload(sharp_minscale)));
    }

  while((!full_ieee) && (l<=lmax))
    {
    Tv a=vload(coef[il].a), b=vload(coef[il].b);
    Tv atmp[4] = {vzero, vzero, vzero, vzero};
    full_ieee=1;
    for (int i=0; i<nv2; ++i)
      {
      atmp[0] += d->lam2[i]*d->corfac[i]*d->p1r[i];
      atmp[1] += d->lam2[i]*d->corfac[i]*d->p1i[i];
      atmp[2] += d->lam2[i]*d->corfac[i]*d->p2r[i];
      atmp[3] += d->lam2[i]*d->corfac[i]*d->p2i[i];
      Tv tmp = (a*d->csq[i] + b)*d->lam2[i] + d->lam1[i];
      d->lam1[i] = d->lam2[i];
      d->lam2[i] = tmp;
      if (rescale(&d->lam1[i], &d->lam2[i], &d->scale[i], vload(sharp_ftol)))
        getCorfac(d->scale[i], &d->corfac[i], gen->cf);
      full_ieee &= vallTrue(vge(d->scale[i],vload(sharp_minscale)));
      }
    vhsum_cmplx_special (atmp[0], atmp[1], atmp[2], atmp[3], &alm[l]);
    l+=2; ++il;
    }
  if (l>lmax) return;

  for (int i=0; i<nv2; ++i)
    {
    d->lam1[i] *= d->corfac[i];
    d->lam2[i] *= d->corfac[i];
    }
  map2alm_kernel(d, coef, alm, l, il, lmax, nv2);
  }

NOINLINE static void iter_to_ieee_spin (const sharp_Ylmgen_C * restrict gen,
  sxdata_v * restrict d, int * restrict l_, int nv2)
  {
  const sharp_ylmgen_dbl2 * restrict fx = gen->coef;
  Tv prefac=vload(gen->prefac[gen->m]),
     prescale=vload(gen->fscale[gen->m]);
  Tv limscale=vload(sharp_limscale);
  int below_limit=1;
  for (int i=0; i<nv2; ++i)
    {
    Tv cth2=vmax(vload(1e-15),vsqrt((vone+d->cth[i])*vload(0.5)));
    Tv sth2=vmax(vload(1e-15),vsqrt((vone-d->cth[i])*vload(0.5)));
    Tm mask=vlt(d->sth[i],vzero);
    vmuleq_mask(vand_mask(mask,vlt(d->cth[i],vzero)),cth2,vload(-1.));
    vmuleq_mask(vand_mask(mask,vgt(d->cth[i],vzero)),sth2,vload(-1.));

    Tv ccp, ccps, ssp, ssps, csp, csps, scp, scps;
    mypow(cth2,gen->cosPow,gen->powlimit,&ccp,&ccps);
    mypow(sth2,gen->sinPow,gen->powlimit,&ssp,&ssps);
    mypow(cth2,gen->sinPow,gen->powlimit,&csp,&csps);
    mypow(sth2,gen->cosPow,gen->powlimit,&scp,&scps);

    d->l1p[i] = vzero;
    d->l1m[i] = vzero;
    d->l2p[i] = prefac*ccp;
    d->scp[i] = prescale+ccps;
    d->l2m[i] = prefac*csp;
    d->scm[i] = prescale+csps;
    Tvnormalize(&d->l2m[i],&d->scm[i],sharp_fbighalf);
    Tvnormalize(&d->l2p[i],&d->scp[i],sharp_fbighalf);
    d->l2p[i] *= ssp;
    d->scp[i] += ssps;
    d->l2m[i] *= scp;
    d->scm[i] += scps;
    if (gen->preMinus_p)
      d->l2p[i] = vneg(d->l2p[i]);
    if (gen->preMinus_m)
      d->l2m[i] = vneg(d->l2m[i]);
    if (gen->s&1)
      d->l2p[i] = vneg(d->l2p[i]);

    Tvnormalize(&d->l2m[i],&d->scm[i],sharp_ftol);
    Tvnormalize(&d->l2p[i],&d->scp[i],sharp_ftol);

    below_limit &= vallTrue(vlt(d->scm[i],limscale)) &&
                   vallTrue(vlt(d->scp[i],limscale));
    }

  int l=gen->mhi;

  while (below_limit)
    {
    if (l+2>gen->lmax) {*l_=gen->lmax+1;return;}
    below_limit=1;
    Tv fx10=vload(fx[l+1].a),fx11=vload(fx[l+1].b);
    Tv fx20=vload(fx[l+2].a),fx21=vload(fx[l+2].b);
    for (int i=0; i<nv2; ++i)
      {
      d->l1p[i] = (d->cth[i]*fx10 - fx11)*d->l2p[i] - d->l1p[i];
      d->l1m[i] = (d->cth[i]*fx10 + fx11)*d->l2m[i] - d->l1m[i];
      d->l2p[i] = (d->cth[i]*fx20 - fx21)*d->l1p[i] - d->l2p[i];
      d->l2m[i] = (d->cth[i]*fx20 + fx21)*d->l1m[i] - d->l2m[i];
      if (rescale(&d->l1p[i],&d->l2p[i],&d->scp[i],vload(sharp_ftol)) ||
          rescale(&d->l1m[i],&d->l2m[i],&d->scm[i],vload(sharp_ftol)))
        below_limit &= vallTrue(vlt(d->scp[i],limscale)) &&
                       vallTrue(vlt(d->scm[i],limscale));
      }
    l+=2;
    }

  *l_=l;
  }

NOINLINE static void alm2map_spin_kernel(sxdata_v * restrict d,
  const sharp_ylmgen_dbl2 * restrict fx, const dcmplx * restrict alm,
  int l, int lmax, int nv2)
  {
  int lsave = l;
  while (l<=lmax)
    {
    Tv fx10=vload(fx[l+1].a),fx11=vload(fx[l+1].b);
    Tv fx20=vload(fx[l+2].a),fx21=vload(fx[l+2].b);
    Tv agr1=vload(creal(alm[2*l  ])), agi1=vload(cimag(alm[2*l  ])),
       acr1=vload(creal(alm[2*l+1])), aci1=vload(cimag(alm[2*l+1]));
    Tv agr2=vload(creal(alm[2*l+2])), agi2=vload(cimag(alm[2*l+2])),
       acr2=vload(creal(alm[2*l+3])), aci2=vload(cimag(alm[2*l+3]));
    for (int i=0; i<nv2; ++i)
      {
      d->l1p[i] = (d->cth[i]*fx10 - fx11)*d->l2p[i] - d->l1p[i];
      d->p1pr[i] += agr1*d->l2p[i];
      d->p1pi[i] += agi1*d->l2p[i];
      d->p1mr[i] += acr1*d->l2p[i];
      d->p1mi[i] += aci1*d->l2p[i];

      d->p1pr[i] += aci2*d->l1p[i];
      d->p1pi[i] -= acr2*d->l1p[i];
      d->p1mr[i] -= agi2*d->l1p[i];
      d->p1mi[i] += agr2*d->l1p[i];
      d->l2p[i] = (d->cth[i]*fx20 - fx21)*d->l1p[i] - d->l2p[i];
      }
    l+=2;
    }
  l=lsave;
  while (l<=lmax)
    {
    Tv fx10=vload(fx[l+1].a),fx11=vload(fx[l+1].b);
    Tv fx20=vload(fx[l+2].a),fx21=vload(fx[l+2].b);
    Tv agr1=vload(creal(alm[2*l  ])), agi1=vload(cimag(alm[2*l  ])),
       acr1=vload(creal(alm[2*l+1])), aci1=vload(cimag(alm[2*l+1]));
    Tv agr2=vload(creal(alm[2*l+2])), agi2=vload(cimag(alm[2*l+2])),
       acr2=vload(creal(alm[2*l+3])), aci2=vload(cimag(alm[2*l+3]));
    for (int i=0; i<nv2; ++i)
      {
      d->l1m[i] = (d->cth[i]*fx10 + fx11)*d->l2m[i] - d->l1m[i];
      d->p2pr[i] -= aci1*d->l2m[i];
      d->p2pi[i] += acr1*d->l2m[i];
      d->p2mr[i] += agi1*d->l2m[i];
      d->p2mi[i] -= agr1*d->l2m[i];

      d->p2pr[i] += agr2*d->l1m[i];
      d->p2pi[i] += agi2*d->l1m[i];
      d->p2mr[i] += acr2*d->l1m[i];
      d->p2mi[i] += aci2*d->l1m[i];
      d->l2m[i] = (d->cth[i]*fx20 + fx21)*d->l1m[i] - d->l2m[i];
      }
    l+=2;
    }
  }

NOINLINE static void calc_alm2map_spin (sharp_job * restrict job,
  const sharp_Ylmgen_C * restrict gen, sxdata_v * restrict d, int nth)
  {
  int l,lmax=gen->lmax;
  int nv2 = (nth+VLEN-1)/VLEN;
  iter_to_ieee_spin(gen, d, &l, nv2);
  job->opcnt += (l-gen->mhi) * 7*nth;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 23*nth;

  const sharp_ylmgen_dbl2 * restrict fx = gen->coef;
  const dcmplx * restrict alm=job->almtmp;
  int full_ieee=1;
  for (int i=0; i<nv2; ++i)
    {
    getCorfac(d->scp[i], &d->cfp[i], gen->cf);
    getCorfac(d->scm[i], &d->cfm[i], gen->cf);
    full_ieee &= vallTrue(vge(d->scp[i],vload(sharp_minscale))) &&
                 vallTrue(vge(d->scm[i],vload(sharp_minscale)));
    }

  while((!full_ieee) && (l<=lmax))
    {
    Tv fx10=vload(fx[l+1].a),fx11=vload(fx[l+1].b);
    Tv fx20=vload(fx[l+2].a),fx21=vload(fx[l+2].b);
    Tv agr1=vload(creal(alm[2*l  ])), agi1=vload(cimag(alm[2*l  ])),
       acr1=vload(creal(alm[2*l+1])), aci1=vload(cimag(alm[2*l+1]));
    Tv agr2=vload(creal(alm[2*l+2])), agi2=vload(cimag(alm[2*l+2])),
       acr2=vload(creal(alm[2*l+3])), aci2=vload(cimag(alm[2*l+3]));
    full_ieee=1;
    for (int i=0; i<nv2; ++i)
      {
      d->l1p[i] = (d->cth[i]*fx10 - fx11)*d->l2p[i] - d->l1p[i];
      d->l1m[i] = (d->cth[i]*fx10 + fx11)*d->l2m[i] - d->l1m[i];

      Tv l2p=d->l2p[i]*d->cfp[i], l2m=d->l2m[i]*d->cfm[i];
      Tv l1m=d->l1m[i]*d->cfm[i], l1p=d->l1p[i]*d->cfp[i];

      d->p1pr[i] += agr1*l2p + aci2*l1p;
      d->p1pi[i] += agi1*l2p - acr2*l1p;
      d->p1mr[i] += acr1*l2p - agi2*l1p;
      d->p1mi[i] += aci1*l2p + agr2*l1p;

      d->p2pr[i] += agr2*l1m - aci1*l2m;
      d->p2pi[i] += agi2*l1m + acr1*l2m;
      d->p2mr[i] += acr2*l1m + agi1*l2m;
      d->p2mi[i] += aci2*l1m - agr1*l2m;

      d->l2p[i] = (d->cth[i]*fx20 - fx21)*d->l1p[i] - d->l2p[i];
      d->l2m[i] = (d->cth[i]*fx20 + fx21)*d->l1m[i] - d->l2m[i];
      if (rescale(&d->l1p[i], &d->l2p[i], &d->scp[i], vload(sharp_ftol)))
        getCorfac(d->scp[i], &d->cfp[i], gen->cf);
      full_ieee &= vallTrue(vge(d->scp[i],vload(sharp_minscale)));
      if (rescale(&d->l1m[i], &d->l2m[i], &d->scm[i], vload(sharp_ftol)))
        getCorfac(d->scm[i], &d->cfm[i], gen->cf);
      full_ieee &= vallTrue(vge(d->scm[i],vload(sharp_minscale)));
      }
    l+=2;
    }
//  if (l>lmax) return;

  for (int i=0; i<nv2; ++i)
    {
    d->l1p[i] *= d->cfp[i];
    d->l2p[i] *= d->cfp[i];
    d->l1m[i] *= d->cfm[i];
    d->l2m[i] *= d->cfm[i];
    }
  alm2map_spin_kernel(d, fx, alm, l, lmax, nv2);

  for (int i=0; i<nv2; ++i)
    {
    Tv tmp;
    tmp = d->p1pr[i]; d->p1pr[i] -= d->p2mi[i]; d->p2mi[i] += tmp;
    tmp = d->p1pi[i]; d->p1pi[i] += d->p2mr[i]; d->p2mr[i] -= tmp;
    tmp = d->p1mr[i]; d->p1mr[i] += d->p2pi[i]; d->p2pi[i] -= tmp;
    tmp = d->p1mi[i]; d->p1mi[i] -= d->p2pr[i]; d->p2pr[i] += tmp;
    }
  }

NOINLINE static void map2alm_spin_kernel(sxdata_v * restrict d,
  const sharp_ylmgen_dbl2 * restrict fx, dcmplx * restrict alm,
  int l, int lmax, int nv2)
  {
  int lsave=l;
  while (l<=lmax)
    {
    Tv fx10=vload(fx[l+1].a),fx11=vload(fx[l+1].b);
    Tv fx20=vload(fx[l+2].a),fx21=vload(fx[l+2].b);
    Tv agr1=vzero, agi1=vzero, acr1=vzero, aci1=vzero;
    Tv agr2=vzero, agi2=vzero, acr2=vzero, aci2=vzero;
    for (int i=0; i<nv2; ++i)
      {
      d->l1p[i] = (d->cth[i]*fx10 - fx11)*d->l2p[i] - d->l1p[i];
      agr1 += d->p2mi[i]*d->l2p[i];
      agi1 -= d->p2mr[i]*d->l2p[i];
      acr1 -= d->p2pi[i]*d->l2p[i];
      aci1 += d->p2pr[i]*d->l2p[i];
      agr2 += d->p2pr[i]*d->l1p[i];
      agi2 += d->p2pi[i]*d->l1p[i];
      acr2 += d->p2mr[i]*d->l1p[i];
      aci2 += d->p2mi[i]*d->l1p[i];
      d->l2p[i] = (d->cth[i]*fx20 - fx21)*d->l1p[i] - d->l2p[i];
      }
    vhsum_cmplx_special (agr1,agi1,acr1,aci1,&alm[2*l]);
    vhsum_cmplx_special (agr2,agi2,acr2,aci2,&alm[2*l+2]);
    l+=2;
    }
  l=lsave;
  while (l<=lmax)
    {
    Tv fx10=vload(fx[l+1].a),fx11=vload(fx[l+1].b);
    Tv fx20=vload(fx[l+2].a),fx21=vload(fx[l+2].b);
    Tv agr1=vzero, agi1=vzero, acr1=vzero, aci1=vzero;
    Tv agr2=vzero, agi2=vzero, acr2=vzero, aci2=vzero;
    for (int i=0; i<nv2; ++i)
      {
      d->l1m[i] = (d->cth[i]*fx10 + fx11)*d->l2m[i] - d->l1m[i];
      agr1 += d->p1pr[i]*d->l2m[i];
      agi1 += d->p1pi[i]*d->l2m[i];
      acr1 += d->p1mr[i]*d->l2m[i];
      aci1 += d->p1mi[i]*d->l2m[i];
      agr2 -= d->p1mi[i]*d->l1m[i];
      agi2 += d->p1mr[i]*d->l1m[i];
      acr2 += d->p1pi[i]*d->l1m[i];
      aci2 -= d->p1pr[i]*d->l1m[i];
      d->l2m[i] = (d->cth[i]*fx20 + fx21)*d->l1m[i] - d->l2m[i];
      }
    vhsum_cmplx_special (agr1,agi1,acr1,aci1,&alm[2*l]);
    vhsum_cmplx_special (agr2,agi2,acr2,aci2,&alm[2*l+2]);
    l+=2;
    }
  }

NOINLINE static void calc_map2alm_spin (sharp_job * restrict job,
  const sharp_Ylmgen_C * restrict gen, sxdata_v * restrict d, int nth)
  {
  int l,lmax=gen->lmax;
  int nv2 = (nth+VLEN-1)/VLEN;
  iter_to_ieee_spin(gen, d, &l, nv2);
  job->opcnt += (l-gen->mhi) * 7*nth;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 23*nth;

  const sharp_ylmgen_dbl2 * restrict fx = gen->coef;
  dcmplx * restrict alm=job->almtmp;
  int full_ieee=1;
  for (int i=0; i<nv2; ++i)
    {
    getCorfac(d->scp[i], &d->cfp[i], gen->cf);
    getCorfac(d->scm[i], &d->cfm[i], gen->cf);
    full_ieee &= vallTrue(vge(d->scp[i],vload(sharp_minscale))) &&
                 vallTrue(vge(d->scm[i],vload(sharp_minscale)));
    }
  for (int i=0; i<nv2; ++i)
    {
    Tv tmp;
    tmp = d->p1pr[i]; d->p1pr[i] -= d->p2mi[i]; d->p2mi[i] += tmp;
    tmp = d->p1pi[i]; d->p1pi[i] += d->p2mr[i]; d->p2mr[i] -= tmp;
    tmp = d->p1mr[i]; d->p1mr[i] += d->p2pi[i]; d->p2pi[i] -= tmp;
    tmp = d->p1mi[i]; d->p1mi[i] -= d->p2pr[i]; d->p2pr[i] += tmp;
    }

  while((!full_ieee) && (l<=lmax))
    {
    Tv fx10=vload(fx[l+1].a),fx11=vload(fx[l+1].b);
    Tv fx20=vload(fx[l+2].a),fx21=vload(fx[l+2].b);
    Tv agr1=vzero, agi1=vzero, acr1=vzero, aci1=vzero;
    Tv agr2=vzero, agi2=vzero, acr2=vzero, aci2=vzero;
    full_ieee=1;
    for (int i=0; i<nv2; ++i)
      {
      d->l1p[i] = (d->cth[i]*fx10 - fx11)*d->l2p[i] - d->l1p[i];
      d->l1m[i] = (d->cth[i]*fx10 + fx11)*d->l2m[i] - d->l1m[i];
      Tv l2p = d->l2p[i]*d->cfp[i], l2m = d->l2m[i]*d->cfm[i];
      Tv l1p = d->l1p[i]*d->cfp[i], l1m = d->l1m[i]*d->cfm[i];
      agr1 += d->p1pr[i]*l2m + d->p2mi[i]*l2p;
      agi1 += d->p1pi[i]*l2m - d->p2mr[i]*l2p;
      acr1 += d->p1mr[i]*l2m - d->p2pi[i]*l2p;
      aci1 += d->p1mi[i]*l2m + d->p2pr[i]*l2p;
      agr2 += d->p2pr[i]*l1p - d->p1mi[i]*l1m;
      agi2 += d->p2pi[i]*l1p + d->p1mr[i]*l1m;
      acr2 += d->p2mr[i]*l1p + d->p1pi[i]*l1m;
      aci2 += d->p2mi[i]*l1p - d->p1pr[i]*l1m;

      d->l2p[i] = (d->cth[i]*fx20 - fx21)*d->l1p[i] - d->l2p[i];
      d->l2m[i] = (d->cth[i]*fx20 + fx21)*d->l1m[i] - d->l2m[i];
      if (rescale(&d->l1p[i], &d->l2p[i], &d->scp[i], vload(sharp_ftol)))
        getCorfac(d->scp[i], &d->cfp[i], gen->cf);
      full_ieee &= vallTrue(vge(d->scp[i],vload(sharp_minscale)));
      if (rescale(&d->l1m[i], &d->l2m[i], &d->scm[i], vload(sharp_ftol)))
        getCorfac(d->scm[i], &d->cfm[i], gen->cf);
      full_ieee &= vallTrue(vge(d->scm[i],vload(sharp_minscale)));
      }
    vhsum_cmplx_special (agr1,agi1,acr1,aci1,&alm[2*l]);
    vhsum_cmplx_special (agr2,agi2,acr2,aci2,&alm[2*l+2]);
    l+=2;
    }
  if (l>lmax) return;

  for (int i=0; i<nv2; ++i)
    {
    d->l1p[i] *= d->cfp[i];
    d->l2p[i] *= d->cfp[i];
    d->l1m[i] *= d->cfm[i];
    d->l2m[i] *= d->cfm[i];
    }
  map2alm_spin_kernel(d, fx, alm, l, lmax, nv2);
  }


NOINLINE static void alm2map_deriv1_kernel(sxdata_v * restrict d,
  const sharp_ylmgen_dbl2 * restrict fx, const dcmplx * restrict alm,
  int l, int lmax, int nv2)
  {
  int lsave=l;
  while (l<=lmax)
    {
    Tv fx10=vload(fx[l+1].a),fx11=vload(fx[l+1].b);
    Tv fx20=vload(fx[l+2].a),fx21=vload(fx[l+2].b);
    Tv ar1=vload(creal(alm[l  ])), ai1=vload(cimag(alm[l  ])),
       ar2=vload(creal(alm[l+1])), ai2=vload(cimag(alm[l+1]));
    for (int i=0; i<nv2; ++i)
      {
      d->l1p[i] = (d->cth[i]*fx10 - fx11)*d->l2p[i] - d->l1p[i];
      d->p1pr[i] += ar1*d->l2p[i];
      d->p1pi[i] += ai1*d->l2p[i];

      d->p1mr[i] -= ai2*d->l1p[i];
      d->p1mi[i] += ar2*d->l1p[i];
      d->l2p[i] = (d->cth[i]*fx20 - fx21)*d->l1p[i] - d->l2p[i];
      }
    l+=2;
    }
  l=lsave;
  while (l<=lmax)
    {
    Tv fx10=vload(fx[l+1].a),fx11=vload(fx[l+1].b);
    Tv fx20=vload(fx[l+2].a),fx21=vload(fx[l+2].b);
    Tv ar1=vload(creal(alm[l  ])), ai1=vload(cimag(alm[l  ])),
       ar2=vload(creal(alm[l+1])), ai2=vload(cimag(alm[l+1]));
    for (int i=0; i<nv2; ++i)
      {
      d->l1m[i] = (d->cth[i]*fx10 + fx11)*d->l2m[i] - d->l1m[i];
      d->p2mr[i] += ai1*d->l2m[i];
      d->p2mi[i] -= ar1*d->l2m[i];

      d->p2pr[i] += ar2*d->l1m[i];
      d->p2pi[i] += ai2*d->l1m[i];
      d->l2m[i] = (d->cth[i]*fx20 + fx21)*d->l1m[i] - d->l2m[i];
      }
    l+=2;
    }
  }

NOINLINE static void calc_alm2map_deriv1(sharp_job * restrict job,
  const sharp_Ylmgen_C * restrict gen, sxdata_v * restrict d, int nth)
  {
  int l,lmax=gen->lmax;
  int nv2 = (nth+VLEN-1)/VLEN;
  iter_to_ieee_spin(gen, d, &l, nv2);
  job->opcnt += (l-gen->mhi) * 7*nth;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 15*nth;

  const sharp_ylmgen_dbl2 * restrict fx = gen->coef;
  const dcmplx * restrict alm=job->almtmp;
  int full_ieee=1;
  for (int i=0; i<nv2; ++i)
    {
    getCorfac(d->scp[i], &d->cfp[i], gen->cf);
    getCorfac(d->scm[i], &d->cfm[i], gen->cf);
    full_ieee &= vallTrue(vge(d->scp[i],vload(sharp_minscale))) &&
                 vallTrue(vge(d->scm[i],vload(sharp_minscale)));
    }

  while((!full_ieee) && (l<=lmax))
    {
    Tv fx10=vload(fx[l+1].a),fx11=vload(fx[l+1].b);
    Tv fx20=vload(fx[l+2].a),fx21=vload(fx[l+2].b);
    Tv ar1=vload(creal(alm[l  ])), ai1=vload(cimag(alm[l  ])),
       ar2=vload(creal(alm[l+1])), ai2=vload(cimag(alm[l+1]));
    full_ieee=1;
    for (int i=0; i<nv2; ++i)
      {
      d->l1p[i] = (d->cth[i]*fx10 - fx11)*d->l2p[i] - d->l1p[i];
      d->l1m[i] = (d->cth[i]*fx10 + fx11)*d->l2m[i] - d->l1m[i];

      Tv l2p=d->l2p[i]*d->cfp[i], l2m=d->l2m[i]*d->cfm[i];
      Tv l1m=d->l1m[i]*d->cfm[i], l1p=d->l1p[i]*d->cfp[i];

      d->p1pr[i] += ar1*l2p;
      d->p1pi[i] += ai1*l2p;
      d->p1mr[i] -= ai2*l1p;
      d->p1mi[i] += ar2*l1p;

      d->p2pr[i] += ar2*l1m;
      d->p2pi[i] += ai2*l1m;
      d->p2mr[i] += ai1*l2m;
      d->p2mi[i] -= ar1*l2m;

      d->l2p[i] = (d->cth[i]*fx20 - fx21)*d->l1p[i] - d->l2p[i];
      d->l2m[i] = (d->cth[i]*fx20 + fx21)*d->l1m[i] - d->l2m[i];
      if (rescale(&d->l1p[i], &d->l2p[i], &d->scp[i], vload(sharp_ftol)))
        getCorfac(d->scp[i], &d->cfp[i], gen->cf);
      full_ieee &= vallTrue(vge(d->scp[i],vload(sharp_minscale)));
      if (rescale(&d->l1m[i], &d->l2m[i], &d->scm[i], vload(sharp_ftol)))
        getCorfac(d->scm[i], &d->cfm[i], gen->cf);
      full_ieee &= vallTrue(vge(d->scm[i],vload(sharp_minscale)));
      }
    l+=2;
    }
//  if (l>lmax) return;

  for (int i=0; i<nv2; ++i)
    {
    d->l1p[i] *= d->cfp[i];
    d->l2p[i] *= d->cfp[i];
    d->l1m[i] *= d->cfm[i];
    d->l2m[i] *= d->cfm[i];
    }
  alm2map_deriv1_kernel(d, fx, alm, l, lmax, nv2);

  for (int i=0; i<nv2; ++i)
    {
    Tv tmp;
    tmp = d->p1pr[i]; d->p1pr[i] -= d->p2mi[i]; d->p2mi[i] += tmp;
    tmp = d->p1pi[i]; d->p1pi[i] += d->p2mr[i]; d->p2mr[i] -= tmp;
    tmp = d->p1mr[i]; d->p1mr[i] += d->p2pi[i]; d->p2pi[i] -= tmp;
    tmp = d->p1mi[i]; d->p1mi[i] -= d->p2pr[i]; d->p2pr[i] += tmp;
    }
  }


#define VZERO(var) do { memset(&(var),0,sizeof(var)); } while(0)

NOINLINE static void inner_loop_a2m(sharp_job *job, const int *ispair,
  const double *cth_, const double *sth_, int llim, int ulim,
  sharp_Ylmgen_C *gen, int mi, const int *mlim)
  {
  const int m = job->ainfo->mval[mi];
  sharp_Ylmgen_prepare (gen, m);

  switch (job->type)
    {
    case SHARP_ALM2MAP:
    case SHARP_ALM2MAP_DERIV1:
      {
      if (job->spin==0)
        {
        //adjust the a_lm for the new algorithm
        dcmplx * restrict alm=job->almtmp;
        for (int il=0, l=gen->m; l<=gen->lmax; ++il,l+=2)
          {
          dcmplx al = alm[l];
          dcmplx al1 = (l+1>gen->lmax) ? 0. : alm[l+1];
          dcmplx al2 = (l+2>gen->lmax) ? 0. : alm[l+2];
          alm[l  ] = gen->alpha[il]*(gen->eps[l+1]*al + gen->eps[l+2]*al2);
          alm[l+1] = gen->alpha[il]*al1;
          }

        const int nval=nv0*VLEN;
        int ith=0;
        int itgt[nval];
        while (ith<ulim-llim)
          {
          s0data_u d;
          VZERO(d.s.p1r); VZERO(d.s.p1i); VZERO(d.s.p2r); VZERO(d.s.p2i);
          int nth=0;
          while ((nth<nval)&&(ith<ulim-llim))
            {
            if (mlim[ith]>=m)
              {
              itgt[nth] = ith;
              d.s.csq[nth]=cth_[ith]*cth_[ith];
              d.s.sth[nth]=sth_[ith];
              ++nth;
              }
            else
              {
              int phas_idx = ith*job->s_th + mi*job->s_m;
              job->phase[phas_idx] = job->phase[phas_idx+1] = 0;
              }
            ++ith;
            }
          if (nth>0)
            {
            int i2=((nth+VLEN-1)/VLEN)*VLEN;
            for (int i=nth; i<i2; ++i)
              {
              d.s.csq[i]=d.s.csq[nth-1];
              d.s.sth[i]=d.s.sth[nth-1];
              d.s.p1r[i]=d.s.p1i[i]=d.s.p2r[i]=d.s.p2i[i]=0.;
              }
            calc_alm2map (job, gen, &d.v, nth);
            for (int i=0; i<nth; ++i)
              {
              int tgt=itgt[i];
              //adjust for new algorithm
              d.s.p2r[i]*=cth_[tgt];
              d.s.p2i[i]*=cth_[tgt];
              int phas_idx = tgt*job->s_th + mi*job->s_m;
              complex double r1 = d.s.p1r[i] + d.s.p1i[i]*_Complex_I,
                             r2 = d.s.p2r[i] + d.s.p2i[i]*_Complex_I;
              job->phase[phas_idx] = r1+r2;
              if (ispair[tgt])
                job->phase[phas_idx+1] = r1-r2;
              }
            }
          }
        }
      else
        {
        //adjust the a_lm for the new algorithm
        if (job->nalm==2)
          for (int l=gen->mhi; l<=gen->lmax+1; ++l)
            {
            job->almtmp[2*l  ]*=gen->alpha[l];
            job->almtmp[2*l+1]*=gen->alpha[l];
            }
        else
          for (int l=gen->mhi; l<=gen->lmax+1; ++l)
            job->almtmp[l]*=gen->alpha[l];

        const int nval=nvx*VLEN;
        int ith=0;
        int itgt[nval];
        while (ith<ulim-llim)
          {
          sxdata_u d;
          VZERO(d.s.p1pr); VZERO(d.s.p1pi); VZERO(d.s.p2pr); VZERO(d.s.p2pi);
          VZERO(d.s.p1mr); VZERO(d.s.p1mi); VZERO(d.s.p2mr); VZERO(d.s.p2mi);
          int nth=0;
          while ((nth<nval)&&(ith<ulim-llim))
            {
            if (mlim[ith]>=m)
              {
              itgt[nth] = ith;
              d.s.cth[nth]=cth_[ith]; d.s.sth[nth]=sth_[ith];
              ++nth;
              }
            else
              {
              int phas_idx = ith*job->s_th + mi*job->s_m;
              job->phase[phas_idx  ] = job->phase[phas_idx+1] = 0;
              job->phase[phas_idx+2] = job->phase[phas_idx+3] = 0;
              }
            ++ith;
            }
          if (nth>0)
            {
            int i2=((nth+VLEN-1)/VLEN)*VLEN;
            for (int i=nth; i<i2; ++i)
              {
              d.s.cth[i]=d.s.cth[nth-1];
              d.s.sth[i]=d.s.sth[nth-1];
              d.s.p1pr[i]=d.s.p1pi[i]=d.s.p2pr[i]=d.s.p2pi[i]=0.;
              d.s.p1mr[i]=d.s.p1mi[i]=d.s.p2mr[i]=d.s.p2mi[i]=0.;
              }
            (job->type==SHARP_ALM2MAP) ?
              calc_alm2map_spin  (job, gen, &d.v, nth) :
              calc_alm2map_deriv1(job, gen, &d.v, nth);
            for (int i=0; i<nth; ++i)
              {
              int tgt=itgt[i];
              int phas_idx = tgt*job->s_th + mi*job->s_m;
              complex double q1 = d.s.p1pr[i] + d.s.p1pi[i]*_Complex_I,
                             q2 = d.s.p2pr[i] + d.s.p2pi[i]*_Complex_I,
                             u1 = d.s.p1mr[i] + d.s.p1mi[i]*_Complex_I,
                             u2 = d.s.p2mr[i] + d.s.p2mi[i]*_Complex_I;
              job->phase[phas_idx  ] = q1+q2;
              job->phase[phas_idx+2] = u1+u2;
              if (ispair[tgt])
                {
                dcmplx *phQ = &(job->phase[phas_idx+1]),
                       *phU = &(job->phase[phas_idx+3]);
                *phQ = q1-q2;
                *phU = u1-u2;
                if ((gen->mhi-gen->m+gen->s)&1)
                  { *phQ=-(*phQ); *phU=-(*phU); }
                }
              }
            }
          }
        }
      break;
      }
    default:
      {
      UTIL_FAIL("must not happen");
      break;
      }
    }
  }

NOINLINE static void inner_loop_m2a(sharp_job *job, const int *ispair,
  const double *cth_, const double *sth_, int llim, int ulim,
  sharp_Ylmgen_C *gen, int mi, const int *mlim)
  {
  const int m = job->ainfo->mval[mi];
  sharp_Ylmgen_prepare (gen, m);

  switch (job->type)
    {
    case SHARP_MAP2ALM:
      {
      if (job->spin==0)
        {
        const int nval=nv0*VLEN;
        int ith=0;
        while (ith<ulim-llim)
          {
          s0data_u d;
          int nth=0;
          while ((nth<nval)&&(ith<ulim-llim))
            {
            if (mlim[ith]>=m)
              {
              d.s.csq[nth]=cth_[ith]*cth_[ith]; d.s.sth[nth]=sth_[ith];
              int phas_idx = ith*job->s_th + mi*job->s_m;
              dcmplx ph1=job->phase[phas_idx];
              dcmplx ph2=ispair[ith] ? job->phase[phas_idx+1] : 0.;
              d.s.p1r[nth]=creal(ph1+ph2); d.s.p1i[nth]=cimag(ph1+ph2);
              d.s.p2r[nth]=creal(ph1-ph2); d.s.p2i[nth]=cimag(ph1-ph2);
              //adjust for new algorithm
              d.s.p2r[nth]*=cth_[ith];
              d.s.p2i[nth]*=cth_[ith];
              ++nth;
              }
            ++ith;
            }
          if (nth>0)
            {
            int i2=((nth+VLEN-1)/VLEN)*VLEN;
            for (int i=nth; i<i2; ++i)
              {
              d.s.csq[i]=d.s.csq[nth-1];
              d.s.sth[i]=d.s.sth[nth-1];
              d.s.p1r[i]=d.s.p1i[i]=d.s.p2r[i]=d.s.p2i[i]=0.;
              }
            calc_map2alm (job, gen, &d.v, nth);
            }
          }
        //adjust the a_lm for the new algorithm
        dcmplx * restrict alm=job->almtmp;
        dcmplx alm2 = 0.;
        double alold=0;
        for (int il=0, l=gen->m; l<=gen->lmax; ++il,l+=2)
          {
          dcmplx al = alm[l];
          dcmplx al1 = (l+1>gen->lmax) ? 0. : alm[l+1];
          alm[l  ] = gen->alpha[il]*gen->eps[l+1]*al + alold*gen->eps[l]*alm2;
          alm[l+1] = gen->alpha[il]*al1;
          alm2=al;
          alold=gen->alpha[il];
          }
        }
      else
        {
        const int nval=nvx*VLEN;
        int ith=0;
        while (ith<ulim-llim)
          {
          sxdata_u d;
          int nth=0;
          while ((nth<nval)&&(ith<ulim-llim))
            {
            if (mlim[ith]>=m)
              {
              d.s.cth[nth]=cth_[ith]; d.s.sth[nth]=sth_[ith];
              int phas_idx = ith*job->s_th + mi*job->s_m;
              dcmplx p1Q=job->phase[phas_idx],
                     p1U=job->phase[phas_idx+2],
                     p2Q=ispair[ith] ? job->phase[phas_idx+1]:0.,
                     p2U=ispair[ith] ? job->phase[phas_idx+3]:0.;
              if ((gen->mhi-gen->m+gen->s)&1)
                { p2Q=-p2Q; p2U=-p2U; }
              d.s.p1pr[nth]=creal(p1Q+p2Q); d.s.p1pi[nth]=cimag(p1Q+p2Q);
              d.s.p1mr[nth]=creal(p1U+p2U); d.s.p1mi[nth]=cimag(p1U+p2U);
              d.s.p2pr[nth]=creal(p1Q-p2Q); d.s.p2pi[nth]=cimag(p1Q-p2Q);
              d.s.p2mr[nth]=creal(p1U-p2U); d.s.p2mi[nth]=cimag(p1U-p2U);
              ++nth;
              }
            ++ith;
            }
          if (nth>0)
            {
            int i2=((nth+VLEN-1)/VLEN)*VLEN;
            for (int i=nth; i<i2; ++i)
              {
              d.s.cth[i]=d.s.cth[nth-1];
              d.s.sth[i]=d.s.sth[nth-1];
              d.s.p1pr[i]=d.s.p1pi[i]=d.s.p2pr[i]=d.s.p2pi[i]=0.;
              d.s.p1mr[i]=d.s.p1mi[i]=d.s.p2mr[i]=d.s.p2mi[i]=0.;
              }
            calc_map2alm_spin(job, gen, &d.v, nth);
            }
          }
        //adjust the a_lm for the new algorithm
        for (int l=gen->mhi; l<=gen->lmax; ++l)
          {
          job->almtmp[2*l  ]*=gen->alpha[l];
          job->almtmp[2*l+1]*=gen->alpha[l];
          }
        }
      break;
      }
    default:
      {
      UTIL_FAIL("must not happen");
      break;
      }
    }
  }

void XARCH(inner_loop) (sharp_job *job, const int *ispair,
  const double *cth_, const double *sth_, int llim, int ulim,
  sharp_Ylmgen_C *gen, int mi, const int *mlim);
void XARCH(inner_loop) (sharp_job *job, const int *ispair,
  const double *cth_, const double *sth_, int llim, int ulim,
  sharp_Ylmgen_C *gen, int mi, const int *mlim)
  {
  (job->type==SHARP_MAP2ALM) ?
    inner_loop_m2a(job,ispair,cth_,sth_,llim,ulim,gen,mi,mlim) :
    inner_loop_a2m(job,ispair,cth_,sth_,llim,ulim,gen,mi,mlim);
  }

#undef VZERO

int XARCH(sharp_veclen)(void);
int XARCH(sharp_veclen)(void)
  {
  return VLEN;
  }

int XARCH(sharp_max_nvec)(int spin);
int XARCH(sharp_max_nvec)(int spin)
  {
  return (spin==0) ? nv0 : nvx;
  }

#define xstr(a) str(a)
#define str(a) #a
const char *XARCH(sharp_architecture)(void);
const char *XARCH(sharp_architecture)(void)
  {
  return xstr(ARCH);
  }

#endif

#endif
