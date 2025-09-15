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

/** Collection of utility functions.

    @copyright 2011-2014 Max-Planck-Society
    @author Martin Reinecke */
public class HealpixUtils
  {
  /** Evaluates cond; if it is false, throws an Exception containing errtxt. */
  static public void check(boolean cond, String errtxt) throws Exception
    { if (!cond) throw new Exception(errtxt); }

  /** Integer base 2 logarithm.
      @param arg
      @return the largest integer {@code n} that fulfills {@code 2^n<=arg}.
      For negative arguments and zero, 0 is returned. */
  static public int ilog2(long arg)
    {
    return 63-Long.numberOfLeadingZeros(Math.max(arg,1L));
    }

  /** Integer square root.
      @param arg
      @return the integer {@code n} which fulfills {@code n^2<=arg<(n+1)^2} */
  static public int isqrt(long arg)
    {
    long res = (long)Math.sqrt(((double)arg)+0.5);
    if (arg<(1L<<50)) return (int)res;
    if (res*res>arg)
      --res;
    else if ((res+1)*(res+1)<=arg)
      ++res;
    return (int)res;
    }

  /** Computes the cosine of the angular distance between two z, phi positions
      on the unit sphere. */
  static public double cosdist_zphi (double z1, double phi1,
    double z2, double phi2)
    {
    return z1*z2+ FastMath.cos(phi1-phi2)* Math.sqrt((1.0-z1*z1)*(1.0-z2*z2));
    }
  /** Computes the cosine of the angular distance between two z, phi positions
      on the unit sphere. */
  static public double cosdist_zphi (Zphi zp1, Zphi zp2)
    { return cosdist_zphi(zp1.z,zp1.phi,zp2.z,zp2.phi); }

  /** Returns the remainder of the division {@code v1/v2}.
      The result is non-negative.
      @param v1 dividend; can be positive or negative
      @param v2 divisor; must be positive
      @return remainder of the division; positive and smaller than {@code v2} */
  static public double fmodulo (double v1, double v2)
    {
    if (v1>=0)
      return (v1<v2) ? v1 : v1%v2;
    double tmp=v1%v2+v2;
    return (tmp==v2) ? 0. : tmp;
    }

  static public boolean approx (float a, float b, float epsilon)
    { return Math.abs(a-b) < (epsilon*Math.abs(b)); }
  static public boolean approx (double a, double b, double epsilon)
    { return Math.abs(a-b) < (epsilon*Math.abs(b)); }

  static public int uniq2order (long uniq)
    { return ilog2(uniq>>>2)>>>1; }

  static private long nest_peano_helper (long pix, int order, int dir)
    {
    final short arr2[] =
      {   0, 35, 65, 66, 68,  5,103,  6,110,109, 15, 44, 72,  9,107, 10,
         31,126, 60,125, 81, 16, 82, 51,123, 88, 26, 25,119, 84, 22, 21,
         42, 75, 41,104, 12, 47, 77, 78, 38, 71, 37,100, 98, 97,  3, 32,
         53, 54,116, 87, 57, 58,120, 91, 19,114, 48,113, 93, 28, 94, 63,
         64,  1, 99,  2, 46, 79, 45,108,  4, 39, 69, 70,  8, 43, 73, 74,
         85, 20, 86, 55,115, 80, 18, 17, 89, 24, 90, 59, 61, 62,124, 95,
        106,105, 11, 40,102,101,  7, 36, 76, 13,111, 14, 34, 67, 33, 96,
        127, 92, 30, 29, 27,122, 56,121, 49, 50,112, 83, 23,118, 52,117,
        128,194,195,161,196,133,135,230,204,141,143,238,171,233,232,138,
        149,212,214,183,221,159,158,252,217,155,154,248,178,243,241,144,
        175,237,236,142,235,170,168,201,227,162,160,193,132,198,199,165,
        186,251,249,152,242,176,177,211,246,180,181,215,157,220,222,191,
        192,129,131,226,136,202,203,169,140,206,207,173,231,166,164,197,
        213,151,150,244,145,208,210,179,153,216,218,187,254,188,189,223,
        239,174,172,205,167,229,228,134,163,225,224,130,200,137,139,234,
        250,184,185,219,190,255,253,156,182,247,245,148,209,147,146,240 };
    final byte arr[] = { 16, 1,27, 2,31,20, 6, 5,10,19, 9,24,13,14,28,23,
                          0,11,17,18,21, 4,22,15,26,25, 3, 8, 7,30,12,29,
                         48,33,35,58,53,39,38,60,59,42,40,49,62,44,45,55,
                         32,50,51,41,37,52,54,47,43,57,56,34,46,63,61,36 };
    final byte peano_face2path[][] =
      { { 2,5,2,5,3,6,3,6,2,3,2,3 }, { 2,6,2,3,3,5,2,6,2,3,3,5 } };
    final byte peano_face2face[][] =
      { { 0,5,6,11,10,1,4,7,2,3,8,9 }, { 0,5,8,9,6,1,2,7,10,11,4,3 } };

    int face = (int)(pix>>>(2*order));
    long result = 0L;
    int state = ((peano_face2path[dir][face]<<4))|(dir<<7);
    int shift=2*order-4;
    for (; shift>=0; shift-=4)
      {
      state=arr2[(state&0xF0) | (((int)(pix>>>shift))&0xF)];
      result = (result<<4) | (state&0xF);
      }
    if (shift==-2)
      {
      state=arr[((state>>>2)&0xFC) | ((int)(pix)&0x3)];
      result = (result<<2) | (state&0x3);
      }

    return result + (((long)peano_face2face[dir][face])<<(2*order));
    }

  static public long nest2peano(long pix, int order)
    { return nest_peano_helper(pix,order,0); }
  static public long peano2nest(long pix, int order)
    { return nest_peano_helper(pix,order,1); }

  }
