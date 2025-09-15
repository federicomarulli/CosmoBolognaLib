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
import java.text.DecimalFormat;

import healpix.essentials.*;

/** Benchmarks for the HealpixBase class

    @copyright 2014 Max-Planck-Society
    @author Martin Reinecke */
public class HealpixBasePerfTest extends TestCase
  {
  static private final DecimalFormat form = new DecimalFormat("##0.00");
  static private final int nsteps=100000; // influences number of performance tests

  private double subtest_perf_neighbours(String name, Scheme scheme)
    throws Exception
    {
    long cnt=0;
    double dummy=0;
    long tstart = System.nanoTime();
    int omax=HealpixBase.order_max;
    for (int order=0; order<=omax; ++order)
      {
      HealpixBase base = new HealpixBase (1L<<order,scheme);
      long dpix=Math.max(base.getNpix()/nsteps,1L);
      for (long pix=0; pix<base.getNpix(); pix+=dpix)
        {
        long nres[] = base.neighbours(pix);
        dummy+=nres[0];
        ++cnt;
        }
      }
    double time = 1e-9*(System.nanoTime()-tstart);
    System.out.println(name + ": " + form.format(cnt/time*1e-6) + "MOps/s");
    return dummy;
    }

  private double subtest_perf_pix2ang(String name, Scheme scheme)
    throws Exception
    {
    long cnt=0;
    double dummy=0;
    long tstart = System.nanoTime();
    int omax=HealpixBase.order_max;
    for (int order=0; order<=omax; ++order)
      {
      HealpixBase base = new HealpixBase (1L<<order,scheme);
      long dpix=Math.max(base.getNpix()/nsteps,1L);
      for (long pix=0; pix<base.getNpix(); pix+=dpix)
        {
        Pointing p=base.pix2ang(pix);
        dummy+=p.theta+p.phi;
        ++cnt;
        }
      }
    double time = 1e-9*(System.nanoTime()-tstart);
    System.out.println(name + ": " + form.format(cnt/time*1e-6) + "MOps/s");
    return dummy;
    }

  private double subtest_perf_pix2vec(String name, Scheme scheme)
    throws Exception
    {
    long cnt=0;
    double dummy=0;
    long tstart = System.nanoTime();
    int omax=HealpixBase.order_max;
    for (int order=0; order<=omax; ++order)
      {
      HealpixBase base = new HealpixBase (1L<<order,scheme);
      long dpix=Math.max(base.getNpix()/nsteps,1L);
      for (long pix=0; pix<base.getNpix(); pix+=dpix)
        {
        Vec3 v=base.pix2vec(pix);
        dummy+=v.x+v.y+v.z;
        ++cnt;
        }
      }
    double time = 1e-9*(System.nanoTime()-tstart);
    System.out.println(name + ": " + form.format(cnt/time*1e-6) + "MOps/s");
    return dummy;
    }

  private double subtest_perf_pix2zphi(String name, Scheme scheme)
    throws Exception
    {
    long cnt=0;
    double dummy=0;
    long tstart = System.nanoTime();
    int omax=HealpixBase.order_max;
    for (int order=0; order<=omax; ++order)
      {
      HealpixBase base = new HealpixBase (1L<<order,scheme);
      long dpix=Math.max(base.getNpix()/nsteps,1L);
      for (long pix=0; pix<base.getNpix(); pix+=dpix)
        {
        Zphi blah=base.pix2zphi(pix);
        dummy+=blah.z;
        ++cnt;
        }
      }
    double time = 1e-9*(System.nanoTime()-tstart);
    System.out.println(name + ": " + form.format(cnt/time*1e-6) + "MOps/s");
    return dummy;
    }

  private double subtest_perf_ring2nest() throws Exception
    {
    long cnt=0;
    double dummy=0;
    long tstart = System.nanoTime();
    int omax=HealpixBase.order_max;
    for (int order=0; order<=omax; ++order)
      {
      HealpixBase base = new HealpixBase (1L<<order,Scheme.RING);
      long dpix=Math.max(base.getNpix()/nsteps,1L);
      for (long pix=0; pix<base.getNpix(); pix+=dpix)
        {
        dummy+=base.ring2nest(pix);
        ++cnt;
        }
      }
    double time = 1e-9*(System.nanoTime()-tstart);
    System.out.println("ring2nest       : " + form.format(cnt/time*1e-6) + "MOps/s");
    return dummy;
    }

  private double subtest_perf_nest2ring() throws Exception
    {
    long cnt=0;
    double dummy=0;
    long tstart = System.nanoTime();
    int omax=HealpixBase.order_max;
    for (int order=0; order<=omax; ++order)
      {
      HealpixBase base = new HealpixBase (1L<<order,Scheme.RING);
      long dpix=Math.max(base.getNpix()/nsteps,1L);
      for (long pix=0; pix<base.getNpix(); pix+=dpix)
        {
        dummy+=base.nest2ring(pix);
        ++cnt;
        }
      }
    double time = 1e-9*(System.nanoTime()-tstart);
    System.out.println("nest2ring       : " + form.format(cnt/time*1e-6) + "MOps/s");
    return dummy;
    }

  private double subtest_perf_ang2pix(String name, Scheme scheme)
    throws Exception
    {
    long cnt=0;
    double dummy=0;
    long tstart = System.nanoTime();
    int omax=HealpixBase.order_max;
    double dth=Math.PI/Math.sqrt(nsteps);
    double dph=Constants.twopi/Math.sqrt(nsteps);
    for (int order=0; order<=omax; ++order)
      {
      HealpixBase base = new HealpixBase (1L<<order,scheme);
      for (double theta=0; theta<Math.PI; theta+=dth)
        for (double phi=0; phi<2*Math.PI; phi+=dph)
          {
          dummy+=base.ang2pix(new Pointing(theta+1e-15*phi,phi));
          ++cnt;
          }
      }
    double time = 1e-9*(System.nanoTime()-tstart);
    System.out.println(name+": " + form.format(cnt/time*1e-6) + "MOps/s");
    return dummy;
    }

  private double subtest_perf_zphi2pix(String name, Scheme scheme)
    throws Exception
    {
    long cnt=0;
    double dummy=0;
    long tstart = System.nanoTime();
    int omax=HealpixBase.order_max;
    double dz=2./Math.sqrt(nsteps);
    double dph=Constants.twopi/Math.sqrt(nsteps);
    for (int order=0; order<=omax; ++order)
      {
      HealpixBase base = new HealpixBase (1L<<order,scheme);
      for (double z=-1; z<1; z+=dz)
        for (double phi=0; phi<Constants.twopi; phi+=dph)
          {
          dummy+=base.zphi2pix(new Zphi(z,phi));
          ++cnt;
          }
      }
    double time = 1e-9*(System.nanoTime()-tstart);
    System.out.println(name+": " + form.format(cnt/time*1e-6) + "MOps/s");
    return dummy;
    }

  private double subtest_perf_query_disc(String name,
    Scheme scheme) throws Exception
    {
    long cnt=0;
    double dummy=0;
    HealpixBase base = new HealpixBase (1024,scheme);
    long tstart = System.nanoTime();
    for (int m=0; m<1000; ++m)
      {
      RangeSet lrs=base.queryDisc(new Pointing(new Vec3(1,0,0)),Constants.halfpi/9.);
      dummy+=lrs.nranges();
      ++cnt;
      }
    double time = 1e-9*(System.nanoTime()-tstart);
    System.out.println(name+": " + form.format(cnt/time*1e-3) + "kOps/s");
    return dummy;
    }

  private double subtest_perf_query_polygon(String name,
    Scheme scheme) throws Exception
    {
    long cnt=0;
    double dummy=0;
    HealpixBase base = new HealpixBase (1024,scheme);
    Pointing[] corner = new Pointing[4];
    corner[0]=new Pointing(new Vec3(1,0.01,0.01));
    corner[1]=new Pointing(new Vec3(1,1,-0.3));
    corner[2]=new Pointing(new Vec3(0.01,1,0.01));
    corner[3]=new Pointing(new Vec3(0.01,0.01,1));
    long tstart = System.nanoTime();
    for (int m=0; m<1000; ++m)
      {
      RangeSet lrs=base.queryPolygon(corner);
      dummy+=lrs.nranges();
      ++cnt;
      }
    double time = 1e-9*(System.nanoTime()-tstart);
    System.out.println(name+": " + form.format(cnt/time*1e-3) + "kOps/s");
    return dummy;
    }

  public void test_perf() throws Exception
    {
    System.out.println("Performance tests of HealpixBase methods");
    double d=0;
    d+=subtest_perf_neighbours("neighbours(NEST)",Scheme.NESTED);
    d+=subtest_perf_neighbours("neighbours(RING)",Scheme.RING);
    d+=subtest_perf_pix2ang   ("pix2ang   (NEST)",Scheme.NESTED);
    d+=subtest_perf_pix2ang   ("pix2ang   (RING)",Scheme.RING);
    d+=subtest_perf_ang2pix   ("ang2pix   (NEST)",Scheme.NESTED);
    d+=subtest_perf_ang2pix   ("ang2pix   (RING)",Scheme.RING);
    d+=subtest_perf_pix2vec   ("pix2vec   (NEST)",Scheme.NESTED);
    d+=subtest_perf_pix2vec   ("pix2vec   (RING)",Scheme.RING);
    d+=subtest_perf_pix2zphi  ("pix2zphi  (NEST)",Scheme.NESTED);
    d+=subtest_perf_pix2zphi  ("pix2zphi  (RING)",Scheme.RING);
    d+=subtest_perf_zphi2pix  ("zphi2pix  (NEST)",Scheme.NESTED);
    d+=subtest_perf_zphi2pix  ("zphi2pix  (RING)",Scheme.RING);
    d+=subtest_perf_ring2nest();
    d+=subtest_perf_nest2ring();
    d+=subtest_perf_query_disc("disc      (NEST)",Scheme.NESTED);
    d+=subtest_perf_query_disc("disc      (RING)",Scheme.RING);
    d+=subtest_perf_query_polygon("polygon   (NEST)",Scheme.NESTED);
    d+=subtest_perf_query_polygon("polygon   (RING)",Scheme.RING);
    }
  }
