/*
 *  This file is part of Healpix Java.
 *
 *  This code is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This code is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this code; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.sourceforge.net
 */
package healpix.essentials.test;

import junit.framework.TestCase;
import java.util.Random;
import java.text.DecimalFormat;

import healpix.essentials.*;

/** @author Martin Reinecke */
public class HealpixBaseTest extends TestCase {

  static private final int nsamples=10000; // influences number of correctness tests
  static private final DecimalFormat form = new DecimalFormat("##0.00");

  private Zphi random_zphi(Random rng)
    { return new Zphi(2*rng.nextDouble()-1,2*Math.PI*rng.nextDouble()); }

  private Pointing random_dir(Random rng)
    { return new Pointing(Math.acos(2*rng.nextDouble()-1),2*Math.PI*rng.nextDouble()); }

  public void test_boundaries()throws Exception
    {
    System.out.println("Testing whether the boundary shapes look sane");

    for (int nside=1; nside<=5; ++nside)
      {
      HealpixBase base=new HealpixBase(nside,Scheme.RING);
      for (int pix=0; pix<base.getNpix(); ++pix)
        {
        for (int res=1; res<=50; res+=7)
          {
          Vec3[] points = base.boundaries(pix,res);
          double dmin=1e30, dmax=-1e30;
          for (int i=0; i<points.length; ++i)
            {
            double dv=(points[i].sub(points[(i+1)%points.length])).length();
            assertTrue("error in boundaries",dv!=0);
            dmin = Math.min(dv,dmin);
            dmax = Math.max(dv,dmax);
            }
          assertTrue("error in boundaries",dmax/dmin<=2);
          }
        }
      }
    }

  public void test_accuracy()throws Exception
    {
    System.out.println("Testing accuracy near the poles");

    for (int order=0; order<=HealpixBase.order_max; ++order)
      {
      HealpixBase base = new HealpixBase (1L<<order,Scheme.RING);
      assertTrue("incorrect pix2ang",base.pix2ang(1).theta>0.0);
      }
    }

  public void test_ringnestring()throws Exception
    {
    System.out.println("Testing identity ring2nest(nest2ring(i))==i");
    Random rng = new Random(5);

    for (int order=0; order<=HealpixBase.order_max; ++order)
      {
      HealpixBase base = new HealpixBase (1L<<order,Scheme.RING);
      for (int m=0; m<nsamples; ++m)
        {
        long pix = (long)(rng.nextDouble()*base.getNpix());
        assertEquals ("ringnestring problem",
          base.ring2nest(base.nest2ring(pix)),pix);
        }
      }
    }

  public void test_pixzphipix()throws Exception
    {
    System.out.println("Testing identity pix2zphi(zphi2pix(i))==i");
    int omax=HealpixBase.order_max;
    Random rng = new Random(5);
    for (int order=0; order<=omax; ++order)
      {
      HealpixBase base1 = new HealpixBase (1L<<order,Scheme.RING),
                  base2 = new HealpixBase (1L<<order,Scheme.NESTED);
      for (int m=0; m<nsamples; ++m)
        {
        long pix = (long)(rng.nextDouble()*base1.getNpix());
        assertEquals ("pixzphipix problem",
          pix, base1.zphi2pix(base1.pix2zphi(pix)));
        assertEquals ("pixzphipix problem",
          pix, base2.zphi2pix(base2.pix2zphi(pix)));
        }
      }
    for (long nside=3; nside<(1L<<omax); nside+=nside/2+1)
      {
      HealpixBase base = new HealpixBase (nside,Scheme.RING);
      for (int m=0; m<nsamples; ++m)
        {
        long pix = (long)(rng.nextDouble()*base.getNpix());
        assertEquals ("pixzphipix problem",
          pix, base.zphi2pix(base.pix2zphi(pix)));
        }
      }
    }

  public void test_zphipixzphi() throws Exception
    {
    System.out.println
      ("Testing near-identity zphi2pix(pix2zphi(i))approx i");
    int omax=HealpixBase.order_max;
    Random rng = new Random(5);
    for (int order=0; order<=omax; ++order)
      {
      HealpixBase base1 = new HealpixBase (1L<<order,Scheme.RING),
                  base2 = new HealpixBase (1L<<order,Scheme.NESTED);
      double mincos = Math.min (Math.cos(base1.maxPixrad()),0.999999999999999);
      for (int m=0; m<nsamples; ++m)
        {
        Zphi zp1 = random_zphi(rng);
        Zphi zp2 = base1.pix2zphi(base1.zphi2pix(zp1));
        double cd1 = HealpixUtils.cosdist_zphi(zp1.z,zp1.phi,zp2.z,zp2.phi);
        assertTrue("zphipixzphi problem "+cd1+" "+mincos,cd1>=mincos);
        zp2 = base2.pix2zphi(base2.zphi2pix(zp1));
        cd1 = HealpixUtils.cosdist_zphi(zp1.z,zp1.phi,zp2.z,zp2.phi);
        assertTrue("zphipixzphi problem",cd1>=mincos);
        }
      }
    for (long nside=3; nside<(1L<<omax); nside+=nside/2+1)
      {
      HealpixBase base = new HealpixBase (nside,Scheme.RING);
      double mincos = Math.min (Math.cos(base.maxPixrad()),0.999999999999999);
      for (int m=0; m<nsamples; ++m)
        {
        Zphi zp1 = random_zphi(rng);
        Zphi zp2 = base.pix2zphi(base.zphi2pix(zp1));
        double cd1 = HealpixUtils.cosdist_zphi(zp1.z,zp1.phi,zp2.z,zp2.phi);
        assertTrue("zphipixzphi problem",cd1>=mincos);
        }
      }
    }

  public void test_neighbours() throws Exception
    {
    System.out.println("Testing neighbour function");
    int omax=HealpixBase.order_max;
    Random rng = new Random(5);
    for (int order=0; order<=omax; ++order)
      {
      HealpixBase base1 = new HealpixBase (1L<<order,Scheme.RING),
                  base2 = new HealpixBase (1L<<order,Scheme.NESTED);
      double maxang = 2.01*base1.maxPixrad();
      for (int m=0; m<nsamples; ++m)
        {
        long pix = (long)(rng.nextDouble()*base1.getNpix());
        Vec3 v = base1.pix2vec(pix);
        long [] nb = base1.neighbours(pix);
        long [] nb2 = base2.neighbours(base2.ring2nest(pix));
        int nnb=0;
        for (int n=0; n<8; ++n)
          {
          if (nb[n]<0)
            assertTrue("neighbour problem 1",nb2[n]<0);
          else
            {
            ++nnb;
            assertEquals("neighbour problem 2",base1.ring2nest(nb[n]),nb2[n]);
            assertTrue("neighbour problem 3",base1.pix2vec(nb[n]).angle(v)<maxang);
            }
          }
        assertTrue("neighbour problem 4 "+order+" "+nnb,(nnb>=7)||((order==0)&&(nnb>=6)));
        }
      }
    for (long nside=3; nside<(1L<<omax); nside+=nside/2+1)
      {
      HealpixBase base = new HealpixBase (nside,Scheme.RING);
      double maxang = 2.01*base.maxPixrad();
      for (int m=0; m<nsamples; ++m)
        {
        long pix = (long)(rng.nextDouble()*base.getNpix());
        Vec3 v = base.pix2vec(pix);
        long[] nb = base.neighbours(pix);
        int nnb=0;
        for (int n=0; n<8; ++n)
          if (nb[n]>=0)
            {
            ++nnb;
            assertTrue("neighbour problem 5",base.pix2vec(nb[n]).angle(v)<maxang);
            }
        assertTrue("neighbour problem 6",nnb>=7);
        }
      }
    }
  public void test_query_disc_strict() throws Exception
    {
    System.out.println("Testing non-inclusive queryDisc()");
    Random rng = new Random(5);
    for (int order=0; order<=5; ++order)
      {
      HealpixBase base = new HealpixBase (1L<<order,Scheme.NESTED);
      boolean[] map = new boolean[(int)base.getNpix()];
      Vec3[] vmap = new Vec3[(int)base.getNpix()];
      for (int m=0; m<base.getNpix(); ++m)
        {
        map[m]=false;
        vmap[m]=base.pix2vec(m);
        }
      for (int m=0; m<nsamples; ++m)
        {
        Pointing ptg = random_dir (rng);
        double rad = Math.PI * rng.nextDouble();
        RangeSet rs = base.queryDisc(ptg,rad);
        Vec3 vptg = new Vec3(ptg);
        double cosrad=Math.cos(rad);
        for (int i=0; i<rs.nranges(); ++i)
          for (long j=rs.ivbegin(i); j<rs.ivend(i); ++j)
            map[(int)j]=true;
        for (int i=0; i<base.getNpix(); ++i)
          {
          boolean inside = vmap[i].dot(vptg)>cosrad;
          assertFalse ("query_disc_strict problem",inside^map[i]);
          }
        for (int i=0; i<rs.nranges(); ++i)
          for (long j=rs.ivbegin(i); j<rs.ivend(i); ++j)
            map[(int)j]=false;
        }
      }
    }

  public void test_query_disc() throws Exception
    {
    System.out.println("Testing queryDisc() empirically");
    int omax=17;
    Random rng = new Random(5);
    for (int order=0; order<=omax; ++order)
      {
      HealpixBase rbase = new HealpixBase (1L<<order,Scheme.RING),
                  nbase = new HealpixBase (1L<<order,Scheme.NESTED);
      int niter=Math.max(1,Math.min(nsamples/1000,100000>>order));
      for (int m=0; m<niter; ++m)
        {
        Pointing ptg = random_dir (rng);
        double rad = Math.PI * rng.nextDouble();
        RangeSet rs = rbase.queryDisc(ptg,rad);
        long nval = rs.nval();
        rs = nbase.queryDisc(ptg,rad);
        assertEquals("queryDisc problem 1", nval,rs.nval());
        rs = rbase.queryDiscInclusive(ptg,rad,4);
        long nv1 = rs.nval();
        rs = nbase.queryDiscInclusive(ptg,rad,4);
        long nv2 = rs.nval();
        assertTrue("queryDisc problem 2", nv1>=nv2);
        assertTrue("queryDisc problem 3", nv2>=nval);
        }
      }
    }

  public void testQueryPolygon() throws Exception
    {
    System.out.println("Testing queryPolygon()");
    HealpixBase base = new HealpixBase(1024,Scheme.NESTED);
    Pointing[] corner = new Pointing[4];
    corner[0]=new Pointing(new Vec3(1,0.01,0.01));
    corner[1]=new Pointing(new Vec3(1,1,-0.3));
    corner[2]=new Pointing(new Vec3(0.01,1,0.01));
    corner[3]=new Pointing(new Vec3(0.01,0.01,1));
    RangeSet lrs=base.queryPolygon(corner);
    assertEquals("QueryPolygon problem",lrs.nval(),1696714);
    lrs=base.queryPolygonInclusive(corner,4);
    assertEquals("QueryPolygon problem",lrs.nval(),1700206);
    base = new HealpixBase(1024,Scheme.RING);
    lrs=base.queryPolygon(corner);
    assertEquals("QueryPolygon problem",lrs.nval(),1696714);
    lrs=base.queryPolygonInclusive(corner,4);
    assertEquals("QueryPolygon problem",lrs.nval(),1700206);
    }

  public void testQueryPolygon2() throws Exception
    {
    System.out.println("Testing queryPolygon() empirically");
    int omax=17;
    Random rng = new Random(5);
    for (int order=0; order<=omax; ++order)
      {
      HealpixBase rbase = new HealpixBase (1L<<order,Scheme.RING),
                  nbase = new HealpixBase (1L<<order,Scheme.NESTED);
      int niter=Math.max(1,Math.min(nsamples/1000,100000>>order));
      for (int m=0; m<niter; ++m)
        {
        Pointing[] corner = new Pointing[3];
        corner[0]=random_dir(rng);
        corner[1]=random_dir(rng);
        corner[2]=random_dir(rng);
        RangeSet rs = rbase.queryPolygon(corner);
        long nval = rs.nval();
        rs = nbase.queryPolygon(corner);
        assertEquals("queryPolygon problem 1", nval,rs.nval());
        rs = rbase.queryPolygonInclusive(corner,4);
        long nv1 = rs.nval();
        rs = nbase.queryPolygonInclusive(corner,4);
        long nv2 = rs.nval();
        assertTrue("queryPolygon problem 2", nv1>=nv2);
        assertTrue("queryPolygon problem 3", nv2>=nval);
        }
      }
    }

  public void testSmallPoly()  throws Exception
    {
    System.out.println("Testing small polygon");
    HealpixBase base = new HealpixBase(1024,Scheme.NESTED);
    Pointing[] corner = new Pointing[4];
    corner[0]=new Pointing(2.7190983373838,1.50840508184223);
    corner[1]=new Pointing(2.7190980835706067,1.5084057923700163);
    corner[2]=new Pointing(2.719101654815075,1.5084144044115078);
    corner[3]=new Pointing(2.719101908630556,1.5084136938829422);
    RangeSet lrs=base.queryPolygon(corner);
    }

  public void test() throws Exception
    {
    int nside=256;
    HealpixBase base  = new HealpixBase(nside,Scheme.NESTED);
    HealpixBase base2 = new HealpixBase(nside,Scheme.RING);
    for (int i=0; i<12*nside*nside; ++i)
      {
      assertEquals ("pixel mismatch_nest",i,base.ang2pix(base.pix2ang(i)));
      assertEquals ("pixel mismatch_nest",i,base.vec2pix(base.pix2vec(i)));
      assertEquals ("pixel mismatch_ring",i,base2.ang2pix(base2.pix2ang(i)));
      assertEquals ("pixel mismatch_ring",i,base2.vec2pix(base2.pix2vec(i)));
      assertEquals ("pixel mismatch",i,base.ring2nest(base2.ang2pix(base.pix2ang(i))));
      assertEquals ("pixel mismatch_ringnestring",i,base.ring2nest(base.nest2ring(i)));
      }
    }

  }
