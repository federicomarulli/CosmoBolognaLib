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

import healpix.essentials.*;
import junit.framework.TestCase;
import java.io.*;
import java.util.Random;
import java.util.Arrays;

/** Tests for the Moc class

    @copyright 2015 Max-Planck-Society
    @author Martin Reinecke */
public class MocTest extends TestCase
  {
  private static Moc randomMoc(int num, long start, int dist)
    {
    Random rng = new Random();
    Moc moc = new Moc();
    long curval=start+(1L<<58);
    for (int i=0; i<num; ++i)
      {
      long v1=curval+1+rng.nextInt(dist);
      long v2=v1+1+rng.nextInt(dist);
      moc.addPixelRange(29,v1,v2);
      curval=v2;
      }
    return moc;
    }

  public void testSimple() throws Exception
    {
    Moc moc=new Moc();
    moc.addPixelRange(0,4,5);
    moc.addPixelRange(0,6,7);
    moc.addPixelRange(2,4,17);
    moc.addPixelRange(10,3000000,3000001);

    assertEquals("inconsistency",moc,moc.complement().complement());
    assertEquals("inconsistency",moc,MocStringIO.mocFromString(" 0/4, 6 2/ \t 4 -16 10/3000000 \t\n "));
    assertEquals("inconsistency",moc,MocStringIO.mocFromString("0/6 2/ 5 2/4 2/6- 16 0/4  10/3000000"));
    assertEquals("inconsistency",moc,MocStringIO.mocFromString
      ("{\"0\":[6] , \"2\": [5 ], \"2\":[  4,6,7,8,9,10,11,12,13,14,15,16], \"0\":[4],  \"10\":[3000000]}"));
    assertEquals("inconsistency",moc,MocStringIO.mocFromString(MocStringIO.mocToStringASCII(moc)));
    assertEquals("inconsistency",moc,MocStringIO.mocFromString(MocStringIO.mocToStringJSON(moc)));
    assertEquals("inconsistency",moc,Moc.fromUniq(moc.toUniq()));
    assertEquals("inconsistency",moc.maxOrder(),10);
    Moc xtmp = moc.degradedToOrder(8,false);
    assertTrue("inconsistency",moc.contains(xtmp));
    assertFalse("inconsistency",xtmp.contains(moc));
    assertTrue("inconsistency",xtmp.overlaps(moc));
    xtmp=moc.degradedToOrder(8,true);
    assertFalse("inconsistency",moc.contains(xtmp));
    assertTrue("inconsistency",xtmp.contains(moc));
    assertEquals("inconsistency",moc,Moc.fromCompressed(moc.toCompressed()));
    }

  public void testPeano() throws Exception
    {
    HealpixBase base = new HealpixBase (8192,Scheme.NESTED);
    RangeSet lrs=base.queryDisc(new Pointing(new Vec3(1,0,0)),Constants.halfpi/9.);
    Moc moc=new Moc(lrs,13);

     long[] arr=moc.toUniq();
    for (int i=0; i<arr.length; ++i)
      {
      int order=HealpixUtils.uniq2order(arr[i]);
      long shift=1L<<(2*order+2);
      arr[i]=HealpixUtils.nest2peano(arr[i]-shift,order)+shift;
      }
    Arrays.sort(arr);
    Moc pmoc=Moc.fromUniq(arr);
    arr=pmoc.toUniq();
    for (int i=0; i<arr.length; ++i)
      {
      int order=HealpixUtils.uniq2order(arr[i]);
      long shift=1L<<(2*order+2);
      arr[i]=HealpixUtils.peano2nest(arr[i]-shift,order)+shift;
      }
    Arrays.sort(arr);
    Moc moc2=Moc.fromUniq(arr);
    assertEquals(moc,moc2);
    }

  public void testOps() throws Exception
    {
    int niter = 100;
    Moc full = MocStringIO.mocFromString("0/0-11");
    Moc empty = MocStringIO.mocFromString("");
    for (int iter=0; iter<niter; ++iter)
      {
      Moc a = randomMoc(1000, 0, 100);
      assertEquals(a.complement().complement(),a);
      assertFalse(a.overlaps(a.complement()));
      assertEquals(a.union(a.complement()),full);
      assertEquals(a.intersection(a.complement()),empty);
      }
    }
  }
