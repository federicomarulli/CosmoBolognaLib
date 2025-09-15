package healpix.essentials.test;

import junit.framework.TestCase;
import java.text.DecimalFormat;
import static healpix.essentials.FastMath.*;

/** A class to perform correctness and speed tests for FastMath class
 *
 * @author Naoki Shibata
 */
public class FastMathTest extends TestCase
  {
  private static final double p0=+0.0, m0=-0.0, pinf=Double.POSITIVE_INFINITY,
    minf=Double.NEGATIVE_INFINITY, nan=Double.NaN;

  private static boolean isPlusZero(double x)
    { return x == 0 && Math.copySign(1, x) == 1; }
  private static boolean isMinusZero(double x)
    { return x == 0 && Math.copySign(1, x) == -1; }

  private static boolean cmpDenorm(double x, double y) {
    if (Double.isNaN(x) && Double.isNaN(y)) return true;
    if (x == pinf && y == pinf) return true;
    if (x == minf && y == minf) return true;
    if (!Double.isNaN(x) && !Double.isNaN(y) &&
        !Double.isInfinite(x) && !Double.isInfinite(y))
      return true;
    return false;
    }

  private static void ulptest (double v1, double v2, double ulp)
    { assertTrue ("inconsistency",Math.abs((v1-v2)/Math.ulp(v2))<=ulp); }

  private static void check1 (double d)
    {
    ulptest (sin(d),Math.sin(d),2.);
    ulptest (cos(d),Math.cos(d),2.);
    ulptest (atan(d),Math.atan(d),2.);
    }
  private static void check2 (double d)
    {
    ulptest (asin(d),Math.asin(d),3.);
    ulptest (acos(d),Math.acos(d),2.);
    }

  public void testmain() throws Exception
    {
    assertEquals("inconsistency", Math.PI,(atan2(p0, m0)));
    assertEquals("inconsistency",-Math.PI,(atan2(m0, m0)));
    assertTrue("inconsistency",isPlusZero(atan2(p0, p0)));
    assertTrue("inconsistency",isMinusZero(atan2(m0, p0)));
    assertEquals("inconsistency",3*Math.PI/4,atan2(pinf, minf));
    assertEquals("inconsistency",-3*Math.PI/4,atan2(minf, minf));
    assertEquals("inconsistency",Math.PI/4,atan2(pinf, pinf));
    assertEquals("inconsistency",-Math.PI/4,atan2(minf, pinf));
    assertTrue("inconsistency",Double.isNaN(atan2(nan,nan)));

    double[] a = { 100000.5, 100000, 3, 2.5, 2, 1.5, 1.0, 0.5 };
    double[] b = { nan, pinf, minf };
    double[] c = { nan, pinf, minf, 2, -2 };

    for(int i=0;i<a.length;i++)
      {
      assertEquals("inconsistency",Math.PI,atan2(p0, -a[i]));
      assertEquals("inconsistency",-Math.PI,atan2(m0, -a[i]));
      assertEquals("inconsistency",-Math.PI/2,atan2(-a[i], p0));
      assertEquals("inconsistency",-Math.PI/2,atan2(-a[i], m0));
      assertEquals("inconsistency",Math.PI/2,atan2(a[i], p0));
      assertEquals("inconsistency",Math.PI/2,atan2(a[i], m0));
      assertEquals("inconsistency",Math.PI/2,atan2(pinf,a[i]));
      assertEquals("inconsistency",Math.PI/2,atan2(pinf,-a[i]));
      assertEquals("inconsistency",-Math.PI/2,atan2(minf,a[i]));
      assertEquals("inconsistency",-Math.PI/2,atan2(minf,-a[i]));
      assertEquals("inconsistency",Math.PI,atan2(a[i],minf));
      assertEquals("inconsistency",-Math.PI,atan2(-a[i],minf));
      assertEquals("inconsistency",p0,atan2(a[i],pinf));
      assertEquals("inconsistency",m0,atan2(-a[i],pinf));
      assertTrue("inconsistency",Double.isNaN(atan2(a[i],nan)));
      assertTrue("inconsistency",Double.isNaN(atan2(-a[i],nan)));
      assertTrue("inconsistency",Double.isNaN(atan2(nan,a[i])));
      assertTrue("inconsistency",Double.isNaN(atan2(nan,-a[i])));
      }

    for(int i=0;i<b.length;i++)
      {
      assertTrue("inconsistency",cmpDenorm(sin(b[i]), Math.sin(b[i])));
      assertTrue("inconsistency",cmpDenorm(cos(b[i]), Math.cos(b[i])));
      assertTrue("inconsistency",cmpDenorm(atan(b[i]), Math.atan(b[i])));
      }
    for(int i=0;i<c.length;i++)
      {
      assertTrue("inconsistency",cmpDenorm(asin(c[i]), Math.asin(c[i])));
      assertTrue("inconsistency",cmpDenorm(acos(c[i]), Math.acos(c[i])));
      }

    for(double d = -10;d < 10;d += 0.000001)
      check1(d);
    for(double d = -10000;d < 10000;d += 0.001)
      check1(d);

    for(double d = -1;d < 1;d += 0.0000001)
      check2(d);

    for(double y = -10;y < 10;y += 0.01)
        for(double x = -10;x < 10;x += 0.01)
          ulptest (atan2(y,x),Math.atan2(y,x),2.);
    for(double y = -1000;y < 1000;y += 1.01)
        for(double x = -1000;x < 1000;x += 1.01)
          ulptest (atan2(y,x),Math.atan2(y,x),2.);
    }
  }
