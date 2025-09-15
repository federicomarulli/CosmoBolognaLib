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

import healpix.essentials.Compressor;
import junit.framework.TestCase;
import java.io.*;
import java.util.Random;

/** Tests for the linear interpolative coding compression

    @copyright 2015 Max-Planck-Society
    @author Martin Reinecke */
public class CompressorTest extends TestCase {

  public void testCompressDecompress() throws Exception
    {
    int num=100000;
    Random randomGenerator = new Random();
    long[] input=new long[num];
    input[0]=0;
    for (int i=1; i<num; ++i) input[i]=input[i-1]+1+randomGenerator.nextInt(100);
    byte[] compressed=Compressor.interpol_encode(input,0,input.length);
    long[] output=Compressor.interpol_decode(compressed);
    assertEquals("inconsistency",input.length,output.length);
    for (int i=0; i<output.length; ++i)
      assertEquals("inconsistency",input[i],output[i]);
    }
}
