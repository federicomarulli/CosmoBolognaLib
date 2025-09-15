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

import java.io.*;
import nom.tam.fits.*;

/** Moc FITS I/O routines.
    @copyright 2014-2015 Max-Planck-Society
    @author Martin Reinecke */
public class MocFitsIO
  {
  /** Converts the contents of a FITS input stream to a MOC. */
  public static Moc mocFromFits(InputStream inp) throws Exception
    {
    FitsFactory.setUseHierarch(true);
    FitsFactory.setUseAsciiTables(false);
    FitsFactory.setLongStringsEnabled(true);
    BasicHDU<?> bhdu = (new Fits(inp)).getHDU(1);
    Header head = bhdu.getHeader();

    Object tmp = ((BinaryTable)bhdu.getData()).getFlattenedColumn(0);
    long[] data=null;
    if (tmp instanceof long[])
      data = (long[])tmp;
    else if (tmp instanceof int[])
      {
      int[] tmp2 = (int[])tmp;
      data = new long[tmp2.length];
      for (int i=0; i<tmp2.length; ++i)
        data[i] = (long)tmp2[i];
      }
    else if (tmp instanceof short[])
      {
      short[] tmp2 = (short[])tmp;
      data = new long[tmp2.length];
      for (int i=0; i<tmp2.length; ++i)
        data[i] = (long)tmp2[i];
      }
    else
      HealpixUtils.check(false,"unsupported data format");

    RangeSet ru = new RangeSet();
    for (int i=0; i<data.length; ++i)
      ru.append(data[i]);
    return Moc.fromUniqRS(ru);
    }
  /** Converts the contents of a FITS file to a MOC. */
  public static Moc mocFromFits(String filename) throws Exception
    {
    FileInputStream inp = new FileInputStream(filename);
    Moc moc = mocFromFits(inp);
    inp.close();
    return moc;
    }

  /** Writes the provided Moc to the stream in FITS format. */
  public static void mocToFits(Moc moc, OutputStream out) throws Exception
    {
    FitsFactory.setUseHierarch(true);
    FitsFactory.setUseAsciiTables(false);
    FitsFactory.setLongStringsEnabled(true);
    Fits f = new Fits();
    Object[] table = new Object[1];
    long[] data=moc.toUniq();
    long maxval=0;
    if (data.length>0) maxval=data[data.length-1];
    if (maxval<=0x7fff)
      {
      short[] dtmp = new short[data.length];
      for (int i=0; i<data.length; ++i)
        dtmp[i]=(short)data[i];
      table[0]=dtmp;
      }
    else if (maxval<=0x7FFFFFFF)
      {
      int[] dtmp = new int[data.length];
      for (int i=0; i<data.length; ++i)
        dtmp[i]=(int)data[i];
      table[0]=dtmp;
      }
    else
      table[0] = data;

    f.addHDU(Fits.makeHDU(table));
    BinaryTableHDU bhdu = (BinaryTableHDU) f.getHDU(1);

    bhdu.setColumnName(0, "PIXEL", "");
    bhdu.addValue("PIXTYPE", "HEALPIX", "HEALPix magic value");
    bhdu.addValue("ORDERING", "NUNIQ", "NUNIQ coding method");
    bhdu.addValue("COORDSYS", "C", "mandated by MOC standard");
    bhdu.addValue("MOCORDER", moc.maxOrder(), "MOC resolution (best order)");

    f.write(new DataOutputStream(out));
    out.flush();
    }
  /** Writes the provided Moc to the specified file in FITS format. */
  public static void mocToFits(Moc moc, String filename) throws Exception
    {
    FileOutputStream out = new FileOutputStream(filename);
    mocToFits(moc,out);
    out.close();
    }
  }
