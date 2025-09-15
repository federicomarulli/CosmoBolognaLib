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
import java.util.UUID;
import java.io.File;

/** Tests for the FitsUtil class

    @copyright 2014 Max-Planck-Society
    @author Martin Reinecke */
public class FitsUtilTest extends TestCase
  {
  public void testHPmapf() throws Exception
    {
    String name = UUID.randomUUID().toString()+".fits";
    HealpixMapFloat hpmf = new HealpixMapFloat(128,Scheme.NESTED);
    for (long i=0; i<hpmf.getNpix(); ++i)
      hpmf.setPixel(i,(float)i);
    FitsUtil.writeFloatMap(hpmf,name);
    hpmf.setNsideAndScheme (64, Scheme.RING);
    hpmf=FitsUtil.getFloatMap(name,2,1);
    new File(name).delete();
    assertEquals("Scheme problem",Scheme.NESTED,hpmf.getScheme());
    assertEquals("Nside problem",128,hpmf.getNside());
    for (long i=0; i<hpmf.getNpix(); ++i)
      assertEquals("Value problem",(float)i,hpmf.getPixel(i));
    }
  public void testHPmapd() throws Exception
    {
    String name = UUID.randomUUID().toString()+".fits";
    HealpixMapDouble hpmd = new HealpixMapDouble(128,Scheme.NESTED);
    for (long i=0; i<hpmd.getNpix(); ++i)
      hpmd.setPixel(i,(double)i);
    FitsUtil.writeDoubleMap(hpmd,name);
    hpmd.setNsideAndScheme (64, Scheme.RING);
    hpmd=FitsUtil.getDoubleMap(name,2,1);
    new File(name).delete();
    assertEquals("Scheme problem",Scheme.NESTED,hpmd.getScheme());
    assertEquals("Nside problem",128,hpmd.getNside());
    for (long i=0; i<hpmd.getNpix(); ++i)
      assertEquals("Value problem",(double)i,hpmd.getPixel(i));
    }
  }
