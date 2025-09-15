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

package healpix.essentials;

/** Procedural interface to the {@link HealpixBase} functionality.
    This class is intended for users who prefer a procedural (instead of
    object-oriented) interface to the HEALPix functionality. It should make
    transition from packages like PixTools easier.

    @copyright 2012 Max-Planck-Society
    @author Martin Reinecke */
public abstract class HealpixProc extends HealpixBase
  {
  static final HealpixBase[] bn=new HealpixBase[HealpixBase.order_max+1],
                             br=new HealpixBase[HealpixBase.order_max+1];
  static final double[] mpr =new double[HealpixBase.order_max+1],
                        cmpr=new double[HealpixBase.order_max+1],
                        smpr=new double[HealpixBase.order_max+1];

  static
    {
    try {
      for (int i=0; i<=HealpixBase.order_max; ++i)
        {
        bn[i]=new HealpixBase(1L<<i, Scheme.NESTED);
        br[i]=new HealpixBase(1L<<i, Scheme.RING);
        mpr[i]=bn[i].maxPixrad();
        cmpr[i]=FastMath.cos(mpr[i]);
        smpr[i]=FastMath.sin(mpr[i]);
        }
      }
    catch (Exception Ex) {/*doesn't happen*/}
    }

  private HealpixProc() {} // this class should not be instantiated

  public static double maxPixrad (int order)
    { return bn[order].maxPixrad(); }

  public static long ang2pixNest (int order, Pointing ptg) throws Exception
    { return bn[order].ang2pix(ptg); }
  public static long ang2pixRing (int order, Pointing ptg) throws Exception
    { return br[order].ang2pix(ptg); }

  public static Pointing pix2angNest (int order, long pix) throws Exception
    { return bn[order].pix2ang(pix); }
  public static Pointing pix2angRing (int order, long pix) throws Exception
    { return br[order].pix2ang(pix); }

  public static long vec2pixNest (int order, Vec3 vec) throws Exception
    { return bn[order].vec2pix(vec); }
  public static long vec2pixRing (int order, Vec3 vec) throws Exception
    { return br[order].vec2pix(vec); }

  public static Vec3 pix2vecNest (int order, long pix) throws Exception
    { return bn[order].pix2vec(pix); }
  public static Vec3 pix2vecRing (int order, long pix) throws Exception
    { return br[order].pix2vec(pix); }

  public static long ring2nest (int order, long pix) throws Exception
    { return bn[order].ring2nest(pix); }
  public static long nest2ring (int order, long pix) throws Exception
    { return bn[order].nest2ring(pix); }

  public static long[] neighboursNest (int order, long pix) throws Exception
    { return bn[order].neighbours(pix); }
  public static long[] neighboursRing (int order, long pix) throws Exception
    { return br[order].neighbours(pix); }

  public static Vec3[] boundariesNest(int order, long pix, int step)
    throws Exception
    { return bn[order].boundaries(pix,step); }
  public static Vec3[] boundariesRing(int order, long pix, int step)
    throws Exception
    { return br[order].boundaries(pix,step); }

  public static RangeSet queryDiscNest(int order, Pointing ptg, double radius)
    throws Exception
    { return bn[order].queryDisc(ptg,radius); }
  public static RangeSet queryDiscRing(int order, Pointing ptg, double radius)
    throws Exception
    { return br[order].queryDisc(ptg,radius); }
  public static RangeSet queryDiscInclusiveNest(int order, Pointing ptg,
    double radius, int fact) throws Exception
    { return bn[order].queryDiscInclusive(ptg,radius,fact); }
  public static RangeSet queryDiscInclusiveRing(int order, Pointing ptg,
    double radius, int fact) throws Exception
    { return br[order].queryDiscInclusive(ptg,radius,fact); }

  public static RangeSet queryPolygonNest(int order, Pointing[] vertex)
    throws Exception
    { return bn[order].queryPolygon(vertex); }
  public static RangeSet queryPolygonRing(int order, Pointing[] vertex)
    throws Exception
    { return br[order].queryPolygon(vertex); }
  public static RangeSet queryPolygonInclusiveNest(int order, Pointing[] vertex,
    int fact) throws Exception
    { return bn[order].queryPolygonInclusive(vertex,fact); }
  public static RangeSet queryPolygonInclusiveRing(int order, Pointing[] vertex,
    int fact) throws Exception
    { return br[order].queryPolygonInclusive(vertex,fact); }

  public static RangeSet queryStripNest(int order, double theta1, double theta2,
    boolean inclusive) throws Exception
    { return bn[order].queryStrip(theta1,theta2,inclusive); }
  public static RangeSet queryStripRing(int order, double theta1, double theta2,
    boolean inclusive) throws Exception
    { return br[order].queryStrip(theta1,theta2,inclusive); }
  };
