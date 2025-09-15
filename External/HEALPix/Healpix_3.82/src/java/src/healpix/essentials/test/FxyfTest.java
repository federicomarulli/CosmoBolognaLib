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

import healpix.essentials.*;

/** Tests for the Fxyf class

    @copyright 2015 Max-Planck-Society
    @author Martin Reinecke */
public class FxyfTest extends TestCase {

  public void testFxyf() throws Exception
    {
    System.out.println("Testing Fxyf");
    for (int i=0; i<12; ++i)
      for (int j=0; j<=100; ++j)
        for (int k=0; k<=100; ++k)
          {
          double fx=(0.01*j)*(1-1e-14)+.5e-14,
                 fy=(0.01*k)*(1-1e-14)+.5e-14;
          Fxyf res=new Fxyf(new Fxyf(fx,fy,i).toVec3());
          assertEquals (res.face,i);
          assertEquals (res.fx,fx,5e-15);
          assertEquals (res.fy,fy,5e-15);
          }
    }

  }
