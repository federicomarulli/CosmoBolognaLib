/*
 *  This file is part of Healpix Java.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */
package healpix.essentials;

import java.io.FileOutputStream;

import nom.tam.fits.BasicHDU;
import nom.tam.fits.BinaryTable;
import nom.tam.fits.BinaryTableHDU;
import nom.tam.fits.Fits;
import nom.tam.fits.FitsFactory;
import nom.tam.fits.Header;
import nom.tam.util.BufferedDataOutputStream;

/** Basic FITS I/O support for Healpix maps.

    @copyright 2011 Max-Planck-Society
    @author Martin Reinecke */
public class FitsUtil
  {
  private static String getKey (Header head, String key) throws Exception
    {
    String upkey = key.toUpperCase();
    if (head.containsKey(upkey))
      return head.findCard(upkey).getValue().trim();
    else if (head.containsKey("HIERARCH."+key))
      return head.findCard("HIERARCH."+key).getValue().trim();
    else
      HealpixUtils.check(false,"key not found: "+key);
    return null;
    }

  private static Scheme getScheme(Header head) throws Exception
    {
    String fscheme=getKey(head,"Ordering");
    HealpixUtils.check((fscheme.equals("RING")) || (fscheme.equals("NESTED")),
      "unsupported ordering scheme");
    return (fscheme.equals("RING")) ? Scheme.RING : Scheme.NESTED;
    }

  public static HealpixMapFloat getFloatMap(String filename, int hdu, int col)
    throws Exception
    {
    FitsFactory.setUseHierarch(true);
    FitsFactory.setUseAsciiTables(false);
    FitsFactory.setLongStringsEnabled(true);
    BasicHDU<?> bhdu = (new Fits(filename)).getHDU(hdu-1);
    Header head = bhdu.getHeader();

    Scheme scheme = getScheme(head);
    int nside = Integer.parseInt(getKey(head,"Nside"));

    Object tmp = ((BinaryTable)bhdu.getData()).getFlattenedColumn(col-1);
    float[] data=null;
    if (tmp instanceof float[])
      data = (float[])tmp;
    else if (tmp instanceof double[])
      {
      double[] tmp2 = (double[])tmp;
      data = new float[tmp2.length];
      for (int i=0; i<tmp2.length; ++i)
        data[i] = (float)tmp2[i];
      }
    else
      HealpixUtils.check(false,"unsupported data format");

    HealpixUtils.check(nside*nside*12==data.length,
      "inconsistent Nside and column length");
    return new HealpixMapFloat(data, scheme);
    }
  public static HealpixMapDouble getDoubleMap(String filename, int hdu, int col)
    throws Exception
    {
    FitsFactory.setUseHierarch(true);
    FitsFactory.setUseAsciiTables(false);
    FitsFactory.setLongStringsEnabled(true);
    BasicHDU<?> bhdu = (new Fits(filename)).getHDU(hdu-1);
    Header head = bhdu.getHeader();

    Scheme scheme = getScheme(head);
    int nside = Integer.parseInt(getKey(head,"Nside"));

    Object tmp = ((BinaryTable)bhdu.getData()).getFlattenedColumn(col-1);
    double[] data=null;
    if (tmp instanceof double[])
      data = (double[])tmp;
    else if (tmp instanceof float[])
      {
      float[] tmp2 = (float[])tmp;
      data = new double[tmp2.length];
      for (int i=0; i<tmp2.length; ++i)
        data[i] = tmp2[i];
      }
    else
      HealpixUtils.check(false,"unsupported data format");

    HealpixUtils.check(nside*nside*12==data.length,
      "inconsistent Nside and column length");
    return new HealpixMapDouble(data, scheme);
    }

  public static void writeFloatMap(HealpixMapFloat map, String filename)
    throws Exception
    {
    FitsFactory.setUseHierarch(true);
    FitsFactory.setUseAsciiTables(false);
    Fits f = new Fits();
    Object[] table = new Object[1];
    table[0] = map.getData();

    f.addHDU(Fits.makeHDU(table));
    BinaryTableHDU bhdu = (BinaryTableHDU) f.getHDU(1);

    bhdu.setColumnName(0, "data", "values");

    bhdu.addValue("PIXTYPE", "HEALPIX", "This is a HEALPix map");
    bhdu.addValue("NSIDE", map.getNside(), "HEALPix NSIDE parameter");
    bhdu.addValue("ORDERING", map.getScheme().toString().toUpperCase(),
      "HEALPix ordering scheme");

    FileOutputStream fos = new FileOutputStream(filename);
    BufferedDataOutputStream s = new BufferedDataOutputStream(fos);

    f.write(s);
    s.flush();
    s.close();
    }
  public static void writeDoubleMap(HealpixMapDouble map, String filename)
    throws Exception
    {
    FitsFactory.setUseHierarch(true);
    FitsFactory.setUseAsciiTables(false);
    Fits f = new Fits();
    Object[] table = new Object[1];
    table[0] = map.getData();

    f.addHDU(Fits.makeHDU(table));
    BinaryTableHDU bhdu = (BinaryTableHDU) f.getHDU(1);

    bhdu.setColumnName(0, "data", "values");

    bhdu.addValue("PIXTYPE", "HEALPIX", "This is a HEALPix map");
    bhdu.addValue("NSIDE", map.getNside(), "HEALPix NSIDE parameter");
    bhdu.addValue("ORDERING", map.getScheme().toString().toUpperCase(),
      "HEALPix ordering scheme");

    FileOutputStream fos = new FileOutputStream(filename);
    BufferedDataOutputStream s = new BufferedDataOutputStream(fos);

    f.write(s);
    s.flush();
    s.close();
    }
  }
