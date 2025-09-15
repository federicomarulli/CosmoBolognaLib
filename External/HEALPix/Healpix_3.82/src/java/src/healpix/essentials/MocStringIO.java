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

/** Moc string I/O routines.
    @copyright 2014-2015 Max-Planck-Society
    @author Martin Reinecke */
public class MocStringIO
  {
  /** Parses a string following either the basic ASCII or JSON syntax given in
      the MOC standard document, and converts it into a MOC.
      The string need not obey the rules for a well-formed MOC given in the
      document. */
  public static Moc mocFromString(String in)
    {
    in=in.replaceAll(",+"," "); // replace commas with spaces
    in=in.replaceAll("[\\[\\]{}\"]+",""); // get rid of useless characters
    in=in.replaceAll("\\s*[:/]\\s*","/ "); // stick order indicator to order
    in=in.replaceAll("\\s*-\\s*","-"); // fuse ranges into one token
    String[] tok = in.split("[\\s]+"); // split at whitespace
    RangeSet ru = new RangeSet();
    int order=0;
    long ofs=0;
    for (int i=0; i<tok.length; ++i)
      {
      if (tok[i].isEmpty()) continue;
      if (tok[i].contains("/")) // new order
        {
        tok[i]=tok[i].split("/")[0];
        order=Integer.parseInt(tok[i]);
        ofs = 4*(1L<<(2*order));
        }
      else if (tok[i].contains("-")) // range of values
        {
        String [] lim=tok[i].split("-");
        ru.add(Long.parseLong(lim[0])+ofs,
               Long.parseLong(lim[1])+ofs+1);
        }
      else // single value
        ru.add(Long.parseLong(tok[i])+ofs);
      }
    return Moc.fromUniqRS(ru);
    }
  private static String mocToStringGeneral(Moc moc, boolean json)
    {
    RangeSet ru = moc.toUniqRS();
    StringBuilder s = new StringBuilder();
    if (json) s.append("{");
    boolean firstOrder=true;
    int omax = moc.maxOrder();
    for (int o=0; o<=omax; ++o)
      {
      RangeSet rt = new RangeSet();
      long offset=4*(1L<<(2*o));
      rt.append(offset, 4*offset);
      rt=rt.intersection(ru);
      boolean prefix=false;
      if (!rt.isEmpty())
        {
        for (int iv=0; iv<rt.nranges(); ++iv)
          {
          long a=rt.ivbegin(iv)-offset,
               b=rt.ivend(iv)-offset;
          if (!prefix)
            {
            if (!firstOrder) s.append(json ? ", " : " ");
            firstOrder=false;
            if (json)
              s.append("\"").append(o).append("\":[");
            else
              s.append(o).append("/");
            prefix=true;
            }
          else
            s.append(",");
          if (json)
            {
            for (long i=a;i<b-1;++i)
              s.append(i).append(",");
            s.append(b-1);
            }
          else
            {
            s.append(a);
            if (b>a+1) s.append("-").append(b-1);
            }
          }
        }
      if (json&&prefix)
        s.append("]");
      }
    if (json) s.append("}");
    return s.toString();
    }
  /** Converts the Moc to its basic ASCII representation as described in the MOC
      standard document. The result is well-formed. */
  public static String mocToStringASCII (Moc moc)
    { return mocToStringGeneral(moc,false); }
  /** Converts the Moc to its JSON representation as described in the MOC
      standard document. The result is well-formed. */
  public static String mocToStringJSON(Moc moc)
    { return mocToStringGeneral(moc,true); }
  }
