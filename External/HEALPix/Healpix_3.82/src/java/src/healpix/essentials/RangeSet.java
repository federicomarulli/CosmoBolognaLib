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
import java.util.NoSuchElementException;

/** Class for dealing with sets of integer ranges.
    Ranges are described by the first element and the one-past-last element.
    This code was inspired by Jan Kotek's "LongRangeSet" class, but has been
    completely reimplemented.

    @copyright 2011-2015 Max-Planck-Society
    @author Martin Reinecke */
public class RangeSet
  {
  /** Interface describing an iterator for going through all values in
      a RangeSet object. */
  public interface ValueIterator {
    public boolean hasNext();
    public long next();
    }

  private static final ValueIterator EMPTY_ITER = new ValueIterator()
    {
    public boolean hasNext() { return false; }
    public long next() { throw new NoSuchElementException(); }
    };

  /** Sorted list of interval boundaries. */
  protected long[] r;
  /** Current number of active entries. */
  protected int sz;

  /** Construct new object with initial space for 4 ranges. */
  public RangeSet() { this(4); }
  /** Construct new object with initial capacity for a given number of ranges.
      @param cap number of initially reserved ranges. */
  public RangeSet(int cap)
    {
    if (cap<0) throw new IllegalArgumentException("capacity must be positive");
    r = new long[cap<<1];
    sz=0;
    }
  /** Construct new object from an array of longs.
      @param data */
   public RangeSet(long[] data)
     {
     sz=data.length;
     r = new long[sz];
     System.arraycopy(data,0,r,0,sz);
     checkConsistency();
     }
  /** Construct new object from another RangeSet
      @param other */
  public RangeSet(RangeSet other)
    {
    sz=other.sz;
    r = new long[sz];
    System.arraycopy(other.r,0,r,0,sz);
    }

  /** Checks the object for internal consistency. If a problem is detected,
      an IllegalArgumentException is thrown. */
  public void checkConsistency()
    {
    if ((sz&1)!=0)
      throw new IllegalArgumentException("invalid number of entries");
    for (int i=1; i<sz; ++i)
      if (r[i]<=r[i-1])
        throw new IllegalArgumentException("inconsistent entries");
    }

  private void resize(int newsize)
    {
    if (newsize<sz)
      throw new IllegalArgumentException("requested array size too small");
    if (newsize==r.length) return;
    long[] rnew = new long[newsize];
    System.arraycopy(r,0,rnew,0,sz);
    r = rnew;
    }

  /** Make sure the object can hold at least the given number of entries. */
  public void ensureCapacity(int cap)
    { if (r.length<cap) resize (Math.max(2*r.length,cap)); }
  /** Shrinks the array for the entries to minimum size. */
  public void trimSize()
    { resize(sz); }
  /** Shrinks the array for the entries to minimum size, if it is more than
      twice the minimum size */
  public void trimIfTooLarge()
    { if (r.length-sz>=sz) resize(sz); }

  /** Returns an internal representation of the interval a number belongs to.
      @param val number whose interval is requested
      @return interval number, starting with -1 (smaller than all numbers in the
      RangeSet), 0 (first "on" interval), 1 (first "off" interval etc.), up to
      (and including) sz-1 (larger than all numbers in the RangeSet). */
  private int iiv (long val)
    {
    int count=sz, first=0;
    while (count>0)
      {
      int step=count>>>1, it = first+step;
      if (r[it]<=val)
        { first=++it; count-=step+1; }
      else
        count=step;
      }
    return first-1;
    }

  /** Append a single-value range to the object.
      @param val value to append */
  public void append(long val)
    { append(val,val+1); }

  /** Append a range to the object.
      @param a first long in range
      @param b one-after-last long in range */
  public void append (long a, long b)
    {
    if (a>=b) return;
    if ((sz>0) && (a<=r[sz-1]))
      {
      if (a<r[sz-2]) throw new IllegalArgumentException("bad append operation");
      if (b>r[sz-1]) r[sz-1]=b;
      return;
      }
    ensureCapacity(sz+2);

    r[sz] = a;
    r[sz+1] = b;
    sz+=2;
    }

  /** Append an entire range set to the object. */
  public void append (RangeSet other)
    {
    for (int i=0; i<other.sz; i+=2)
      append(other.r[i],other.r[i+1]);
    }

  /** @return number of ranges in the set. */
  public int nranges()
    { return sz>>>1; }

  /** @return true if no entries are stored, else false. */
  public boolean isEmpty()
    { return sz==0; }

  /** @return first number in range iv. */
  public long ivbegin(int iv)
    { return r[2*iv]; }
  /** @return one-past-last number in range iv. */
  public long ivend(int iv)
    { return r[2*iv+1]; }

  /** Remove all entries in the set. */
  public void clear()
    { sz=0; }

  /** Push a single entry at the end of the entry vector. */
  private void pushv(long v)
    { ensureCapacity(sz+1); r[sz++]=v; }

  /** Estimate a good strategy for set operations involving two RangeSets. */
  private static int strategy (int sza, int szb)
    {
    final double fct1 = 1.;
    final double fct2 = 1.;
    int slo = sza<szb ? sza : szb,
        shi = sza<szb ? szb : sza;
    double cost1 = fct1 * (sza+szb);
    double cost2 = fct2 * slo * Math.max(1.,HealpixUtils.ilog2(shi));
    return (cost1<=cost2) ? 1 : (slo==sza) ? 2 : 3;
    }

  private static boolean generalAllOrNothing1 (RangeSet a, RangeSet b,
    boolean flip_a, boolean flip_b)
    {
    boolean state_a=flip_a, state_b=flip_b, state_res=state_a||state_b;
    int ia=0, ea=a.sz, ib=0, eb=b.sz;
    boolean runa = ia!=ea, runb = ib!=eb;
    while(runa||runb)
      {
      long va = runa ? a.r[ia] : 0L,
           vb = runb ? b.r[ib] : 0L;
      boolean adv_a = runa && (!runb || (va<=vb)),
              adv_b = runb && (!runa || (vb<=va));
      if (adv_a) { state_a=!state_a; ++ia; runa = ia!=ea; }
      if (adv_b) { state_b=!state_b; ++ib; runb = ib!=eb; }
      if ((state_a||state_b)!=state_res)
        return false;
      }
    return true;
    }

  private static boolean generalAllOrNothing2 (RangeSet a, RangeSet b,
    boolean flip_a, boolean flip_b)
    {
    int iva = flip_a ? 0 : -1;
    while (iva<a.sz)
      {
      if (iva==-1) // implies that flip_a==false
        { if ((!flip_b)||(b.r[0]<a.r[0])) return false; }
      else if (iva==a.sz-1) // implies that flip_a==false
        { if ((!flip_b)||(b.r[b.sz-1]>a.r[a.sz-1])) return false; }
      else
        {
        int ivb=b.iiv(a.r[iva]);
        if ((ivb!=b.sz-1)&&(b.r[ivb+1]<a.r[iva+1])) return false;
        if (flip_b==((ivb&1)==0)) return false;
        }
      iva+=2;
      }
    return true;
    }

  private static boolean generalAllOrNothing (RangeSet a, RangeSet b,
    boolean flip_a, boolean flip_b)
    {
    if (a.isEmpty())
      return flip_a ? true : b.isEmpty();
    if (b.isEmpty())
      return flip_b ? true : a.isEmpty();
    int strat = strategy (a.nranges(), b.nranges());
    return (strat==1) ? generalAllOrNothing1(a,b,flip_a,flip_b) :
             ((strat==2) ? generalAllOrNothing2(a,b,flip_a,flip_b)
                         : generalAllOrNothing2(b,a,flip_b,flip_a));
    }

  /** Internal helper function for constructing unions, intersections
      and differences of two RangeSets. */
  private static RangeSet generalUnion1 (RangeSet a, RangeSet b,
    boolean flip_a, boolean flip_b)
    {
    RangeSet res = new RangeSet();

    boolean state_a=flip_a, state_b=flip_b, state_res=state_a||state_b;
    int ia=0, ea=a.sz, ib=0, eb=b.sz;
    boolean runa = ia!=ea, runb = ib!=eb;
    while(runa||runb)
      {
      long va = runa ? a.r[ia] : 0L,
           vb = runb ? b.r[ib] : 0L;
      boolean adv_a = runa && (!runb || (va<=vb)),
              adv_b = runb && (!runa || (vb<=va));
      if (adv_a) { state_a=!state_a; ++ia; runa = ia!=ea; }
      if (adv_b) { state_b=!state_b; ++ib; runb = ib!=eb; }
      if ((state_a||state_b)!=state_res)
        { res.pushv(adv_a ? va : vb); state_res = !state_res; }
      }
    return res;
    }
  /** Internal helper function for constructing unions, intersections
      and differences of two RangeSets. */
  private static RangeSet generalUnion2 (RangeSet a, RangeSet b,
    boolean flip_a, boolean flip_b)
    {
    RangeSet res = new RangeSet();
    int iva = flip_a ? 0 : -1;
    while (iva<a.sz)
      {
      int ivb = (iva==-1) ? -1 : b.iiv(a.r[iva]);
      boolean state_b = flip_b^((ivb&1)==0);
      if ((iva>-1) && (!state_b)) res.pushv(a.r[iva]);
      while((ivb<b.sz-1)&&((iva==a.sz-1)||(b.r[ivb+1]<a.r[iva+1])))
        { ++ivb; state_b=!state_b; res.pushv(b.r[ivb]); }
      if ((iva<a.sz-1)&&(!state_b)) res.pushv(a.r[iva+1]);
      iva+=2;
      }
    return res;
    }
  private static RangeSet generalUnion (RangeSet a, RangeSet b,
    boolean flip_a, boolean flip_b)
    {
    if (a.isEmpty())
      return flip_a ? new RangeSet() : new RangeSet(b);
    if (b.isEmpty())
      return flip_b ? new RangeSet() : new RangeSet(a);
    int strat = strategy (a.nranges(), b.nranges());
    return (strat==1) ? generalUnion1(a,b,flip_a,flip_b) :
             ((strat==2) ? generalUnion2(a,b,flip_a,flip_b)
                         : generalUnion2(b,a,flip_b,flip_a));
    }

  /** Return the union of this RangeSet and other. */
  public RangeSet union (RangeSet other)
    { return generalUnion (this,other,false,false); }
  /** Return the intersection of this RangeSet and other. */
  public RangeSet intersection (RangeSet other)
    { return generalUnion (this,other,true,true); }
  /** Return the difference of this RangeSet and other. */
  public RangeSet difference (RangeSet other)
    { return generalUnion (this,other,true,false); }

  /** Returns true if a is contained in the set, else false. */
  public boolean contains (long a)
    { return ((iiv(a)&1)==0); }
  /** Returns true if all numbers [a;b[ are contained in the set, else false. */
  public boolean contains (long a,long b)
    {
    int res=iiv(a);
    if ((res&1)!=0) return false;
    return (b<=r[res+1]);
    }
  @Deprecated
  public boolean containsAll (long a,long b)
    { return contains(a,b); }
  /** Returns true if any of the numbers [a;b[ are contained in the set,
      else false. */
  public boolean overlaps (long a,long b)
    {
    int res=iiv(a);
    if ((res&1)==0) return true;
    if (res==sz-1) return false; // beyond the end of the set
    return (r[res+1]<b);
    }
  @Deprecated
  public boolean containsAny (long a,long b)
    { return overlaps(a,b); }
  /** Returns true if the set completely contains "other", else false. */
  public boolean contains (RangeSet other)
    { return generalAllOrNothing(this,other,false,true); }
  @Deprecated
  public boolean containsAll (RangeSet other)
    { return contains(other); }
  /** Returns true if there is overlap between the set and "other",
      else false. */
  public boolean overlaps (RangeSet other)
    { return !generalAllOrNothing(this,other,true,true); }
  @Deprecated
  public boolean containsAny (RangeSet other)
    { return overlaps(other); }
  /** Returns true the object represents an identical set of ranges as obj. */
  public boolean equals(Object obj) {
    if (this == obj)
      return true;
    if ((obj==null) || (!(obj instanceof RangeSet)))
      return false;
    RangeSet other = (RangeSet) obj;
    if (other.sz!=sz) return false;
    for (int i=0; i<sz; ++i)
      if (other.r[i]!=r[i]) return false;
    return true;
    }
  public int hashCode()
    {
    int result = Integer.valueOf(sz).hashCode();
    for (int i=0; i<sz; ++i)
      result = 31 * result + Long.valueOf(r[sz]).hashCode();
    return result;
    }

  /** @return total number of values (not ranges) in the set. */
  public long nval()
    {
    long res = 0;
    for (int i=0; i<sz; i+=2)
      res+=r[i+1]-r[i];
    return res;
    }

  /** Internal helper function for building unions and differences of the
      RangeSet with a single range. */
  private void addRemove (long a, long b, int v)
    {
    int pos1=iiv(a), pos2=iiv(b);
    if ((pos1>=0) && (r[pos1]==a)) --pos1;
    // first to delete is at pos1+1; last is at pos2
    boolean insert_a = (pos1&1)==v;
    boolean insert_b = (pos2&1)==v;
    int rmstart=pos1+1+(insert_a ? 1 : 0);
    int rmend  =pos2-(insert_b?1:0);

    if (((rmend-rmstart)&1)==0)
      throw new IllegalArgumentException("cannot happen: "+rmstart+" "+rmend);

    if (insert_a && insert_b && (pos1+1>pos2)) // insert
      {
      ensureCapacity(sz+2);
      System.arraycopy(r,pos1+1,r,pos1+3,sz-pos1-1); // move to right
      r[pos1+1]=a;
      r[pos1+2]=b;
      sz+=2;
      }
    else
      {
      if (insert_a) r[pos1+1]=a;
      if (insert_b) r[pos2]=b;
      if (rmstart!=rmend+1)
        System.arraycopy(r,rmend+1,r,rmstart,sz-rmend-1); // move to left
      sz-=rmend-rmstart+1;
      }
    }

  /** After this operation, the RangeSet contains the intersection of itself and
      [a;b[. */
  public void intersect (long a, long b)
    {
    int pos1=iiv(a), pos2=iiv(b);
    if ((pos2>=0) && (r[pos2]==b)) --pos2;
    // delete all up to pos1 (inclusive); and starting from pos2+1
    boolean insert_a = (pos1&1)==0;
    boolean insert_b = (pos2&1)==0;

    // cut off end
    sz=pos2+1;
    if (insert_b) r[sz++]=b;

    // erase start
    if (insert_a) r[pos1--]=a;
    if (pos1>=0)
      System.arraycopy(r,pos1+1,r,0,sz-pos1-1); // move to left

    sz-=pos1+1;
    if ((sz&1)!=0)
      throw new IllegalArgumentException("cannot happen");
    }

  /** After this operation, the RangeSet contains the union of itself and
      [a;b[. */
  public void add (long a, long b)
    { if ((sz==0)||(a>=r[sz-1])) append(a,b); else addRemove(a,b,1); }
  /** After this operation, the RangeSet contains the union of itself and
      [a;a+1[. */
  public void add (long a)
    { if ((sz==0)||(a>=r[sz-1])) append(a,a+1); else addRemove(a,a+1,1); }
  /** After this operation, the RangeSet contains the difference of itself and
      [a;b[. */
  public void remove (long a, long b)
    { addRemove(a,b,0); }
  /** After this operation, the RangeSet contains the difference of itself and
      [a;a+1[. */
  public void remove (long a)
    { addRemove(a,a+1,0); }

  /** Creates an array containing all the numbers in the RangeSet.
      Not recommended, because the arrays can become prohibitively large.
      It is preferrable to use a ValueIterator or explicit loops. */
  public long[] toArray(){
    long[] res = new long[(int)nval()];
    int ofs=0;
    for (int i=0; i<sz; i+=2)
      for (long j=r[i]; j<r[i+1]; ++j)
        res[ofs++]=j;
    return res;
    }
  public static RangeSet fromArray(long[]v){
    RangeSet res = new RangeSet();
    for (int i=0; i<v.length; i++)
      res.append(v[i]);
    return res;
    }

  public String toString()
    {
    StringBuilder s = new StringBuilder();
    s.append("{ ");
    for (int i=0; i<sz; i+=2)
      {
      s.append("[").append(r[i]).append(";").append(r[i + 1]).append("[");
      if (i<sz-2) s.append(",");
      }
    s.append(" }");
    return s.toString();
    }

  /** Returns a ValueIterator, which iterates over all individual numbers
      in the RangeSet. */
  public ValueIterator valueIterator()
    {
    if(sz == 0) return EMPTY_ITER;
    return new ValueIterator()
      {
      int pos = 0;
      long value = (sz>0) ? r[0] : 0;

      public boolean hasNext()
        { return (pos<sz); }

      public long next() {
        if (pos>sz)
          throw new NoSuchElementException();
        long ret = value;
        if (++value==r[pos+1])
          {
          pos+=2;
          if (pos<sz)
            value = r[pos];
          }
        return ret;
        }
      };
    }

  /** Returns a compressed representation of the RangeSet, using interpolative
      coding. */
  public byte[] toCompressed() throws Exception
    { return Compressor.interpol_encode (r, 0, sz); }
  /** Returns a RangeSet obtained by decompressing a byte array which was
      originally generated by toCompressed(). */
  public static RangeSet fromCompressed (byte[] data) throws Exception
    { return new RangeSet(Compressor.interpol_decode(data)); }
  }
