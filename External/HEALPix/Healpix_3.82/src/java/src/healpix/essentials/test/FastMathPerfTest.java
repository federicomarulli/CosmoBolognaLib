package healpix.essentials.test;

import junit.framework.TestCase;
import java.text.DecimalFormat;
import static healpix.essentials.FastMath.*;

/** FastMath benchmarking code

    @copyright 2014 Max-Planck-Society
    @author Martin Reinecke */
public class FastMathPerfTest extends TestCase
  {
  static private final DecimalFormat form = new DecimalFormat("##0.00");
  final int ncomp=10000000;

  public void test_perf_math()
    {
    System.out.println("Rough performance test of math functions");
    long cnt, tstart;
    double dummy=0, di, time;
    di=1./ncomp;
    tstart = System.nanoTime();
    cnt=0;
    for (double i=0; i<1; i+=di)
       { dummy+=Math.sqrt(i); ++cnt; }
    time = 1e-9*(System.nanoTime()-tstart);
    System.out.println("sqrt: " + form.format(cnt/time*1e-6) + "MOps/s");
    di=4*Math.PI/ncomp;
    tstart = System.nanoTime();
    cnt=0;
    for (double i=-2*Math.PI; i<2*Math.PI; i+=di)
      { dummy+=Math.sin(i); ++cnt; }
    time = 1e-9*(System.nanoTime()-tstart);
    System.out.println("sin : " + form.format(cnt/time*1e-6) + "MOps/s");
    tstart = System.nanoTime();
    cnt=0;
    for (double i=-2*Math.PI; i<2*Math.PI; i+=di)
      { dummy+=Math.cos(i); ++cnt; }
    time = 1e-9*(System.nanoTime()-tstart);
    System.out.println("cos : " + form.format(cnt/time*1e-6) + "MOps/s");
    di=2./ncomp;
    tstart = System.nanoTime();
    cnt=0;
    for (double i=-1; i<1; i+=di)
      { dummy+=Math.acos(i); ++cnt; }
    time = 1e-9*(System.nanoTime()-tstart);
    System.out.println("acos: " + form.format(cnt/time*1e-6) + "MOps/s");
    di=1000./ncomp;
    tstart = System.nanoTime();
    cnt=0;
    for (double i=-500; i<500; i+=di)
      { dummy+=Math.atan(i); ++cnt; }
    time = 1e-9*(System.nanoTime()-tstart);
    System.out.println("atan: " + form.format(cnt/time*1e-6) + "MOps/s");
    di=1000./ncomp;
    tstart = System.nanoTime();
    cnt=0;
    for (double i=-500; i<500; i+=di)
      { dummy+=Math.atan2(i,2.9); ++cnt; }
    time = 1e-9*(System.nanoTime()-tstart);
    System.out.println("atan2: " + form.format(cnt/time*1e-6) + "MOps/s");
    di=1000./ncomp;
    }
  public void test_perf_fastmath()
    {
    System.out.println("Rough performance test of fast math functions");
    long cnt,tstart;
    double dummy=0, di, time;
    di=4*Math.PI/ncomp;
    tstart = System.nanoTime();
    cnt=0;
    for (double i=-2*Math.PI; i<2*Math.PI; i+=di)
      { dummy+=sin(i); ++cnt; }
    time = 1e-9*(System.nanoTime()-tstart);
    System.out.println("sin : " + form.format(cnt/time*1e-6) + "MOps/s");
    tstart = System.nanoTime();
    cnt=0;
    for (double i=-2*Math.PI; i<2*Math.PI; i+=di)
      { dummy+=cos(i); ++cnt; }
    time = 1e-9*(System.nanoTime()-tstart);
    System.out.println("cos : " + form.format(cnt/time*1e-6) + "MOps/s");
    di=2./ncomp;
    tstart = System.nanoTime();
    cnt=0;
    for (double i=-1; i<1; i+=di)
      { dummy+=acos(i); ++cnt; }
    time = 1e-9*(System.nanoTime()-tstart);
    System.out.println("acos: " + form.format(cnt/time*1e-6) + "MOps/s");
    di=1000./ncomp;
    tstart = System.nanoTime();
    cnt=0;
    for (double i=-500; i<500; i+=di)
      { dummy+=atan(i); ++cnt; }
    time = 1e-9*(System.nanoTime()-tstart);
    System.out.println("atan: " + form.format(cnt/time*1e-6) + "MOps/s");
    di=1000./ncomp;
    tstart = System.nanoTime();
    cnt=0;
    for (double i=-500; i<500; i+=di)
      { dummy+=atan2(i,2.9); ++cnt; }
    time = 1e-9*(System.nanoTime()-tstart);
    System.out.println("atan2: " + form.format(cnt/time*1e-6) + "MOps/s");
    }
  }
