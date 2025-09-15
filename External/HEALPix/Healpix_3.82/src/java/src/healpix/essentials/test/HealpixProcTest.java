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

/** Tests for the HealpixProc class

    @copyright 2014 Max-Planck-Society
    @author Martin Reinecke */
public class HealpixProcTest extends TestCase {

  static private final int nsamples=10000; // influences number of correctness tests
  static private final DecimalFormat form = new DecimalFormat("##0.00");

  private Zphi random_zphi(Random rng)
    { return new Zphi(2*rng.nextDouble()-1,2*Math.PI*rng.nextDouble()); }

  private Pointing random_dir(Random rng)
    { return new Pointing(Math.acos(2*rng.nextDouble()-1),2*Math.PI*rng.nextDouble()); }

  public void test_accuracy()throws Exception
    {
    System.out.println("Testing accuracy near the poles");

    for (int order=0; order<=HealpixBase.order_max; ++order)
      assertTrue("incorrect pix2ang",
        HealpixProc.pix2angRing(order,1).theta>0.0);
    }

  public void test_ringnestring()throws Exception
    {
    System.out.println("Testing identity ring2nest(nest2ring(i))==i");
    Random rng = new Random(5);

    for (int o=0; o<=HealpixBase.order_max; ++o)
      {
      for (int m=0; m<nsamples; ++m)
        {
        long pix = (long)(rng.nextDouble()*HealpixProc.order2Npix(o));
        assertEquals ("ringnestring problem",
          HealpixProc.ring2nest(o,HealpixProc.nest2ring(o,pix)),pix);
        }
      }
    }

  public void test_pixangpix()throws Exception
    {
    System.out.println("Testing identity pix2ang(ang2pix(i))==i");
    int omax=HealpixBase.order_max;
    Random rng = new Random(5);
    for (int o=0; o<=omax; ++o)
      {
      for (int m=0; m<nsamples; ++m)
        {
        long pix = (long)(rng.nextDouble()*HealpixProc.order2Npix(o));
        assertEquals ("pixangpix problem",
          pix, HealpixProc.ang2pixNest(o,HealpixProc.pix2angNest(o,pix)));
        assertEquals ("pixangpix problem",
          pix, HealpixProc.ang2pixRing(o,HealpixProc.pix2angRing(o,pix)));
        }
      }
    }

  public void test_neighbours() throws Exception
    {
    System.out.println("Testing neighbour function");
    int omax=HealpixBase.order_max;
    Random rng = new Random(5);
    for (int o=0; o<=omax; ++o)
      {
      double maxang = 2.01*HealpixProc.maxPixrad(o);
      for (int m=0; m<nsamples; ++m)
        {
        long pix = (long)(rng.nextDouble()*HealpixProc.order2Npix(o));
        Vec3 v = HealpixProc.pix2vecRing(o,pix);
        long [] nb = HealpixProc.neighboursRing(o,pix);
        long [] nb2 = HealpixProc.neighboursNest(o,HealpixProc.ring2nest(o,pix));
        int nnb=0;
        for (int n=0; n<8; ++n)
          {
          if (nb[n]<0)
            assertTrue("neighbour problem 1",nb2[n]<0);
          else
            {
            ++nnb;
            assertEquals("neighbour problem 2",HealpixProc.ring2nest(o,nb[n]),nb2[n]);
            assertTrue("neighbour problem 3",HealpixProc.pix2vecRing(o,nb[n]).angle(v)<maxang);
            }
          }
        assertTrue("neighbour problem 4 "+o+" "+nnb,(nnb>=7)||((o==0)&&(nnb>=6)));
        }
      }
    }

  public void test_query_disc_strict() throws Exception
    {
    System.out.println("Testing non-inclusive queryDisc()");
    Random rng = new Random(5);
    for (int o=0; o<=5; ++o)
      {
      int npix=(int)HealpixProc.order2Npix(o);
      boolean[] map = new boolean[npix];
      Vec3[] vmap = new Vec3[npix];
      for (int m=0; m<npix; ++m)
        {
        map[m]=false;
        vmap[m]=HealpixProc.pix2vecRing(o,m);
        }
      for (int m=0; m<nsamples; ++m)
        {
        Pointing ptg = random_dir (rng);
        double rad = Math.PI * rng.nextDouble();
        RangeSet rs = HealpixProc.queryDiscRing(o,ptg,rad);
        Vec3 vptg = new Vec3(ptg);
        double cosrad=Math.cos(rad);
        for (int i=0; i<rs.nranges(); ++i)
          for (long j=rs.ivbegin(i); j<rs.ivend(i); ++j)
            map[(int)j]=true;
        for (int i=0; i<npix; ++i)
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
    for (int o=0; o<=omax; ++o)
      {
      int niter=Math.max(1,Math.min(nsamples/1000,100000>>o));
      for (int m=0; m<niter; ++m)
        {
        Pointing ptg = random_dir (rng);
        double rad = Math.PI * rng.nextDouble();
        RangeSet rs = HealpixProc.queryDiscRing(o,ptg,rad);
        long nval = rs.nval();
        rs = HealpixProc.queryDiscNest(o,ptg,rad);
        assertEquals("queryDisc problem 1", nval,rs.nval());
        rs = HealpixProc.queryDiscInclusiveRing(o,ptg,rad,4);
        long nv1 = rs.nval();
        rs = HealpixProc.queryDiscInclusiveNest(o,ptg,rad,4);
        long nv2 = rs.nval();
        assertTrue("queryDisc problem 2", nv1>=nv2);
        assertTrue("queryDisc problem 3", nv2>=nval);
        }
      }
    }

  public void testQueryPolygon() throws Exception
    {
    System.out.println("Testing queryPolygon()");
    Pointing[] corner = new Pointing[4];
    corner[0]=new Pointing(new Vec3(1,0.01,0.01));
    corner[1]=new Pointing(new Vec3(1,1,-0.3));
    corner[2]=new Pointing(new Vec3(0.01,1,0.01));
    corner[3]=new Pointing(new Vec3(0.01,0.01,1));
    RangeSet lrs=HealpixProc.queryPolygonNest(10,corner);
    assertEquals("QueryPolygon problem",lrs.nval(),1696714);
    lrs=HealpixProc.queryPolygonInclusiveNest(10,corner,4);
    assertEquals("QueryPolygon problem",lrs.nval(),1700206);
    lrs=HealpixProc.queryPolygonRing(10,corner);
    assertEquals("QueryPolygon problem",lrs.nval(),1696714);
    lrs=HealpixProc.queryPolygonInclusiveRing(10,corner,4);
    assertEquals("QueryPolygon problem",lrs.nval(),1700206);
    }

  public void testQueryPolygon2() throws Exception
    {
    System.out.println("Testing queryPolygon() empirically");
    int omax=17;
    Random rng = new Random(5);
    for (int o=0; o<=omax; ++o)
      {
      int niter=Math.max(1,Math.min(nsamples/1000,100000>>o));
      for (int m=0; m<niter; ++m)
        {
        Pointing[] corner = new Pointing[3];
        corner[0]=random_dir(rng);
        corner[1]=random_dir(rng);
        corner[2]=random_dir(rng);
        RangeSet rs = HealpixProc.queryPolygonRing(o,corner);
        long nval = rs.nval();
        rs = HealpixProc.queryPolygonNest(o,corner);
        assertEquals("queryPolygon problem 1", nval,rs.nval());
        rs = HealpixProc.queryPolygonInclusiveRing(o,corner,4);
        long nv1 = rs.nval();
        rs = HealpixProc.queryPolygonInclusiveNest(o,corner,4);
        long nv2 = rs.nval();
        assertTrue("queryPolygon problem 2", nv1>=nv2);
        assertTrue("queryPolygon problem 3", nv2>=nval);
        }
      }
    }

  public void test() throws Exception
    {
    int o=8;
    for (int i=0; i<HealpixProc.order2Npix(o); ++i)
      {
      assertEquals ("pixel mismatch_nest",i,HealpixProc.ang2pixNest(o,HealpixProc.pix2angNest(o,i)));
      assertEquals ("pixel mismatch_nest",i,HealpixProc.vec2pixNest(o,HealpixProc.pix2vecNest(o,i)));
      assertEquals ("pixel mismatch_ring",i,HealpixProc.ang2pixRing(o,HealpixProc.pix2angRing(o,i)));
      assertEquals ("pixel mismatch_ring",i,HealpixProc.vec2pixRing(o,HealpixProc.pix2vecRing(o,i)));
      assertEquals ("pixel mismatch",i,HealpixProc.ring2nest(o,HealpixProc.ang2pixRing(o,HealpixProc.pix2angNest(o,i))));
      assertEquals ("pixel mismatch_ringnestring",i,HealpixProc.ring2nest(o,HealpixProc.nest2ring(o,i)));
      }
    }
  }
