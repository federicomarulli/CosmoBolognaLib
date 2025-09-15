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

/** Support for MOC algorithms.
    @copyright 2014-2015 Max-Planck-Society
    @author Martin Reinecke */
public class Moc
  {
  private final static int maxorder=29;
  private RangeSet rs;

  /** Creates a new, empty Moc. */
  public Moc()
    { rs=new RangeSet(); }
  /** Creates a new Moc, which is identical to "other". */
  public Moc(Moc other)
    { rs=new RangeSet(other.rs); }
  /** Creates a new Moc from the range set of NESTED pixels at the given order.
      */
  public Moc(RangeSet rs2, int order)
    {
    rs=new RangeSet(rs2.nranges());
    int shift=2*(maxorder-order);
    for (int i=0; i<rs2.nranges(); ++i)
      rs.append(rs2.ivbegin(i)<<shift,rs2.ivend(i)<<shift);
    }

  static private Moc fromNewRangeSet(RangeSet rngset)
    {
    Moc res = new Moc();
    res.rs = rngset;
    return res;
    }

  /** Returns the internal range set. */
  public RangeSet getRangeSet()
    { return new RangeSet(rs); }

  /** Returns the maximum HEALPix order necessary to describe the Moc exactly.
      */
  public int maxOrder()
    {
    long combo=0;
    for (int i=0; i<rs.nranges(); ++i)
      combo|=rs.ivbegin(i)|rs.ivend(i);
    return maxorder-(Long.numberOfTrailingZeros(combo)>>>1);
    }

  /** Returns a Moc with degraded resolution.
      @param order the maximum HEALPix order used for the result
      @param keepPartialCells if true, partially filled cells will be included
        in the output Moc; if false, they will be ignored. */
  public Moc degradedToOrder (int order, boolean keepPartialCells)
    {
    int shift=2*(maxorder-order);
    long ofs=(1L<<shift)-1;
    long mask = ~ofs;
    long adda = keepPartialCells ? 0L : ofs,
         addb = keepPartialCells ? ofs : 0L;
    RangeSet rs2=new RangeSet();
    for (int i=0; i<rs.nranges(); ++i)
      {
      long a=(rs.ivbegin(i)+adda)&mask;
      long b=(rs.ivend  (i)+addb)&mask;
      if (b>a) rs2.append(a,b);
      }
    return fromNewRangeSet(rs2);
    }

  /** Adds a range of pixels at a given HEALPix order to the Moc.
      @param order the HEALPix order of the added pixels
      @param p1 the first pixel of the range
      @param p2 the one-after-last pixel of the range */
  public void addPixelRange (int order, long p1, long p2)
    {
    int shift=2*(maxorder-order);
    rs.add(p1<<shift,p2<<shift);
    }
  public void addPixel (int order, long p)
    { addPixelRange(order, p, p+1); }
  /** Returns a new Moc that contains the union of this Moc and "other". */
  public Moc union (Moc other)
    { return fromNewRangeSet(rs.union(other.rs)); }
  /** Returns a new Moc that contains the intersection of this Moc and "other".
      */
  public Moc intersection (Moc other)
    { return fromNewRangeSet(rs.intersection(other.rs)); }
  /** Returns a new Moc that contains all parts of this Moc that are not
      contained in "other". */
  public Moc subtraction (Moc other)
    { return fromNewRangeSet(rs.difference(other.rs)); }
  /** Returns the complement of this Moc. */
  public Moc complement()
    {
    RangeSet full = new RangeSet(new long[]{0L,12L*(1L<<(2*maxorder))});
    return fromNewRangeSet(full.difference(rs));
    }
  /** @return true, if "other" is a subset of this Moc, else false. */
  public boolean contains(Moc other) // FIXME: needs optimization!
    { return rs.contains(other.rs); }
  /** @return true, if the intersection of this Moc and "other" is not empty. */
  public boolean overlaps(Moc other) // FIXME: needs optimization!
    { return rs.overlaps(other.rs); }
  /** @return A RangeSet containing all HEALPix pixels (in NUNIQ order) covered
      by this Moc. The result is well-formed in the sense that every pixel is
      given at its lowest possible HEALPix order. */
  public RangeSet toUniqRS() // should be tuned!
    {
    RangeSet r2 = new RangeSet(rs);
    RangeSet r3 = new RangeSet();
    RangeSet res= new RangeSet();
    for (int o=0; o<=maxorder; ++o)
      {
      if (r2.isEmpty()) return res;

      int shift = 2*(maxorder-o);
      long ofs=(1L<<shift)-1;
      long ofs2 = 1L<<(2*o+2);
      r3.clear();
      for (int iv=0; iv<r2.nranges(); ++iv)
        {
        long a=(r2.ivbegin(iv)+ofs)>>>shift,
             b=r2.ivend(iv)>>>shift;
        r3.append(a<<shift, b<<shift);
        res.append(a+ofs2,b+ofs2);
        }
      if (!r3.isEmpty())
        r2 = r2.difference(r3);
      }
    return res;
    }
  public long[] toUniq() // should be tuned!
    { return toUniqRS().toArray(); }

  /** @return A Moc built from the RangeSet of NUNIQ HEALPix pixels given in
      "ru". "ru" need not be well-formed. */
  public static Moc fromUniqRS (RangeSet ru) // should be tuned!
    {
    RangeSet r= new RangeSet();
    RangeSet rtmp = new RangeSet();
    int lastorder=0;
    int shift=2*maxorder;
    for (int i=0; i<ru.nranges(); ++i)
      for (long j=ru.ivbegin(i); j<ru.ivend(i); ++j)
        {
        int order = HealpixUtils.uniq2order(j);
        if (order!=lastorder)
          {
          r=r.union(rtmp);
          rtmp.clear();
          lastorder=order;
          shift=2*(maxorder-order);
          }
        long pix = j-(1L<<(2*order+2));
        rtmp.append (pix<<shift,(pix+1)<<shift);
        }
    r=r.union(rtmp);
    return fromNewRangeSet(r);
    }

  public static Moc fromUniq (long []u) // should be tuned!
    {
    return fromUniqRS(RangeSet.fromArray(u));
    }

  /** @return A compressed representation of the Moc obtained by interpolative
      coding. */
  public byte[] toCompressed() throws Exception
    { return rs.toCompressed(); }
  /** @return A Moc built from the compressed representation given in "data". */
  public static Moc fromCompressed(byte[] data) throws Exception
    { return fromNewRangeSet(RangeSet.fromCompressed(data)); }

  public boolean equals(Object obj)
    {
    if (this == obj)
      return true;
    if ((obj==null) || (!(obj instanceof Moc)))
      return false;
    Moc other = (Moc) obj;
    return rs.equals(other.rs);
    }
  public int hashCode()
    { return rs.hashCode(); }

  public int nranges()
    { return rs.nranges(); }
  }
