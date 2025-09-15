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
import java.util.ArrayList;

/** Support for MOC queries.
    @copyright 2015 Max-Planck-Society
    @author Martin Reinecke */
public class MocQuery
  {
  static public enum MocQueryOp { AND,OR,XOR,NOT,NONE }

  static public class MocQueryComponent
    {
    public MocQueryOp op;
    public Vec3 center;
    public double radius;
    public int nops;
    public MocQueryComponent(MocQueryOp op_) throws Exception
      {
      op= op_;
      HealpixUtils.check(op_!=MocQueryOp.NONE,"bad operator");
      switch (op)
        {
        case AND:
        case OR:
        case XOR:
          nops=2;
          break;
        case NOT:
          nops=1;
          break;
        case NONE:
          nops=0;
          break;
        }
      }
    public MocQueryComponent(MocQueryOp op_, int nops_) throws Exception
      {
      op= op_;
      nops=nops_;
      switch (op)
        {
        case AND:
        case OR:
          HealpixUtils.check(nops>=2,"bad nops");
          break;
        case XOR:
          HealpixUtils.check(nops==2,"bad nops");
          break;
        case NOT:
          HealpixUtils.check(nops==1,"bad nops");
          break;
        case NONE:
          HealpixUtils.check(false,"bad operator");
          break;
        }
      }
    public MocQueryComponent(Vec3 cnt, double rad)
      {
      op = MocQueryOp.NONE;
      center = new Vec3(cnt);
      center.normalize();
      radius = rad;
      nops=0;
      }
    }

  static private class queryHelper
    {
    private final class pstack
      {
      private long[] p;
      private int[] o;
      private int s, m;

      public pstack (int sz)
        { p=new long[sz]; o=new int[sz]; s=m=0; }
      public void push (long p_, int o_)
        { p[s]=p_; o[s]=o_; ++s; }
      public void pop ()
        { --s; }
      public void popToMark ()
        { s=m; }
      public int size ()
        { return s; }
      public void mark ()
        { m=s; }
      public int otop ()
        { return o[s-1]; }
      public long ptop ()
        { return p[s-1]; }
      public void clear ()
        { s=m=0; }
      public boolean empty ()
        { return s==0; }
      }

    private int order, omax, ncomp;
    private boolean inclusive;
    private HealpixBase base[];
    private double cr[], crmin[][], crmax[][];
    private int shortcut[];
    private MocQueryOp op[];
    private int nops[];
    private Vec3 center[];

    private pstack stk; // stack for pixel numbers and their orders
    private long pix;
    private int o;
    private Vec3 pv;
    private int loc;

    private void check_pixel (int zone, Moc pixset)
      {
      if (zone==0) return;
      if (o<order)
        {
        if (zone>=3)
          pixset.addPixel(o,pix); // output all subpixels
        else // (zone>=1)
          for (int i=0; i<4; ++i)
            stk.push(4*pix+3-i,o+1); // add children
        }
      else if (o>order) // this implies that inclusive==true
        {
        if (zone>=2) // pixel center in shape
          {
          pixset.addPixel(order,pix>>>(2*(o-order))); // output parent pixel
          stk.popToMark(); // unwind the stack
          }
        else // (zone>=1): pixel center in safety range
          {
          if (o<omax) // check sublevels
            for (int i=0; i<4; ++i) // add children in reverse order
              stk.push(4*pix+3-i,o+1);
          else // at resolution limit
            {
            pixset.addPixel(order,pix>>>(2*(o-order))); // output parent pixel
            stk.popToMark(); // unwind the stack
            }
          }
        }
      else // o==order
        {
        if (zone>=2)
          pixset.addPixel(order,pix);
        else if (inclusive) // and (zone>=1)
          {
          if (order<omax) // check sublevels
            {
            stk.mark(); // remember current stack position
            for (int i=0; i<4; ++i) // add children in reverse order
              stk.push(4*pix+3-i,o+1);
            }
          else // at resolution limit
            pixset.addPixel(order,pix); // output the pixel
          }
        }
      }

    void correctLoc() throws Exception
      {
      int myloc=loc--;
      HealpixUtils.check((myloc>=0)&&(myloc<ncomp),"inconsistency");
      for (int i=0; i<nops[myloc]; ++i)
        correctLoc();
      }
    int getZone (int zmin, int zmax)
      {
      if (zmin==zmax) { loc=shortcut[loc]; return zmin; } // short-circuit
      int myloc=loc--;
      switch (op[myloc])
        {
        case AND:
          {
          int z1=zmax;
          for (int i=0; i<nops[myloc]; ++i)
            z1 = getZone(zmin,z1);
          return z1;
          }
        case OR:
          {
          int z1=zmin;
          for (int i=0; i<nops[myloc]; ++i)
            z1 = getZone(z1,zmax);
          return z1;
          }
        case XOR:
          {
          int z1=getZone(0,3);
          int z2=getZone(0,3);
          return Math.max(zmin,Math.min(zmax,
            Math.max(Math.min(z1,3-z2),Math.min(3-z1,z2))));
          }
        case NOT:
          return 3-getZone(3-zmax,3-zmin);
        case NONE:
          {
          int res=zmax;
          double crad=pv.dot(center[myloc]);
          if (crad<=crmax[o][myloc]) res=0;
          else if (crad<=cr[myloc]) res=1;
          else if (crad<=crmin[o][myloc]) res=2;
          return Math.max(zmin,Math.min(zmax,res));
          }
        }
      return -1; // unreachable
      }

    public queryHelper (int order_, int omax_, boolean inclusive_,
      ArrayList<MocQueryComponent> comp) throws Exception
      {
      order=order_;
      omax=omax_;
      ncomp=comp.size();
      inclusive=inclusive_;
      base=new HealpixBase[omax+1];
      cr=new double[ncomp];
      crmin=new double [omax+1][ncomp];
      crmax=new double [omax+1][ncomp];
      HealpixUtils.check(ncomp>=1,"bad query component ArrayList");
      HealpixUtils.check(order<=omax,"order>omax");
      if (!inclusive) HealpixUtils.check(order==omax,"inconsistency");
      HealpixUtils.check(omax<=HealpixBase.order_max,"omax too high");

      op=new MocQueryOp[ncomp];
      nops=new int[ncomp];
      center=new Vec3[ncomp];
      for (int i=0; i<ncomp; ++i)
        {
        op[i]=comp.get(i).op;
        nops[i]=comp.get(i).nops;
        center[i]=comp.get(i).center;
        if (op[i]==MocQueryOp.NONE) // it's a cap
          cr[i]=Math.cos(comp.get(i).radius);
        }
      for (o=0; o<=omax; ++o) // prepare data at the required orders
        {
        base[o] = new HealpixBase(1<<o,Scheme.NESTED);
        double dr=base[o].maxPixrad(); // safety distance
        for (int i=0; i<ncomp; ++i)
          if (op[i]==MocQueryOp.NONE) // it's a cap
            {
            double r=comp.get(i).radius;
            crmax[o][i] = (r+dr>=Math.PI) ?
              -1.01 : Math.cos(r+dr);
            crmin[o][i] = (r-dr<=0.) ?
               1.01 : Math.cos(r-dr);
            }
        }
      stk=new pstack(12+3*omax); // reserve maximum size to avoid reallocation

      shortcut=new int[ncomp];
      for (int i=0; i<ncomp; ++i)
        {
        loc=i;
        correctLoc();
        shortcut[i]=loc;
        }
      }
    Moc result() throws Exception
      {
      Moc pixset = new Moc();
      stk.clear();
      stk.mark();
      for (int i=0; i<12; ++i) // insert the 12 base pixels in reverse order
        stk.push(11-i,0);

      while (!stk.empty()) // as long as there are pixels on the stack
        {
        // pop current pixel number and order from the stack
        pix=stk.ptop();
        o=stk.otop();
        stk.pop();
        pv = base[o].pix2vec(pix);

        loc=ncomp-1;
        int zone=getZone(0,3);
        check_pixel (zone, pixset);
        HealpixUtils.check(loc==-1,"stack not used up");
        }
      return pixset;
      }
    }

  static public Moc doMocQuery (int order, ArrayList<MocQueryComponent> comp)
    throws Exception
    { return (new queryHelper(order,order,false,comp)).result(); }

  static public Moc doMocQueryInclusive (int order, int omax,
    ArrayList<MocQueryComponent> comp)
    throws Exception
    { return (new queryHelper(order,omax,true,comp)).result(); }

  static private double isLeft (Vec3 a, Vec3 b, Vec3 c)
    {
    return (a.cross(b)).dot(c);
    }

  // adapted from code available at http://geomalgorithms.com/a12-_hull-3.html
  // Original copyright notice follows:
  // Copyright 2001 softSurfer, 2012 Dan Sunday
  // This code may be freely used and modified for any purpose
  // providing that this copyright notice is included with it.
  // SoftSurfer makes no warranty for this code, and cannot be held
  // liable for any real or imagined damage resulting from its use.
  // Users of this code must verify correctness for their application.
  static private int[] getHull (Vec3 vert[], int P[])
    throws Exception
    {
    // initialize a deque D[] from bottom to top so that the
    // 1st three vertices of P[] are a ccw triangle
    int n = P.length;
    int D[] = new int[2*n+1];
    int bot = n-2, top = bot+3;    // initial bottom and top deque indices
    D[bot] = D[top] = P[2];      // 3rd vertex is at both bot and top
    if (isLeft(vert[P[0]], vert[P[1]], vert[P[2]]) > 0)
      {
      D[bot+1] = P[0];
      D[bot+2] = P[1];           // ccw vertices are: 2,0,1,2
      }
    else
      {
      D[bot+1] = P[1];
      D[bot+2] = P[0];           // ccw vertices are: 2,1,0,2
      }

    // compute the hull on the deque D[]
    for (int i=3; i<n; i++)
      {   // process the rest of vertices
      // test if next vertex is inside the deque hull
      if ((isLeft(vert[D[bot]], vert[D[bot+1]], vert[P[i]]) > 0) &&
          (isLeft(vert[D[top-1]], vert[D[top]], vert[P[i]]) > 0) )
        continue;         // skip an interior vertex

      // incrementally add an exterior vertex to the deque hull
      // get the rightmost tangent at the deque bot
      while (isLeft(vert[D[bot]], vert[D[bot+1]], vert[P[i]]) <= 0)
        ++bot;                 // remove bot of deque
      D[--bot] = P[i];         // insert P[i] at bot of deque

      // get the leftmost tangent at the deque top
      while (isLeft(vert[D[top-1]], vert[D[top]], vert[P[i]]) <= 0)
        --top;                 // pop top of deque
      D[++top] = P[i];         // push P[i] onto top of deque
      }

    // transcribe deque D[] to the output hull array H[]
    int nout = top-bot;
    int res[] = new int[nout];
    for (int h=0; h<nout; h++)
      res[h] = D[bot + h +1];

    return res;
    }

  static public ArrayList<MocQueryComponent> prepPolyHelper (Vec3 vv[],
    int P[], ArrayList<MocQueryComponent> comp, boolean doLast)
    throws Exception
    {
    int hull[]=getHull(vv,P);
    boolean addHull[]=new boolean[hull.length];

    // sync both sequences at the first point of the convex hull
    int ihull=0, ipoly=0, nhull=hull.length, npoly=P.length;
    while (hull[ihull]!=P[ipoly]) ++ipoly;

    // iterate over the pockets between the polygon and its convex hull
    int npockets=0;
    if (P.length==3)
      for (int i=0; i<3; i++) addHull[i]=true;
    else
      {
      do
        {
        int ihull_next = (ihull+1)%nhull,
            ipoly_next = (ipoly+1)%npoly;
        if (hull[ihull_next]==P[ipoly_next]) // no pocket found
          { addHull[ihull]=true; ihull=ihull_next; ipoly=ipoly_next; }
        else // query pocket
          {
          int nvpocket=2; // number of vertices for this pocket
          while (P[ipoly_next]!=hull[ihull_next])
            {
            ipoly_next = (ipoly_next+1)%npoly;
            ++nvpocket;
            }
          int ppocket[] = new int[nvpocket];
          int idx=0;
          int ipoly_bw=ipoly_next;
          while (P[ipoly_bw]!=hull[ihull])
            {
            ppocket[idx++]=P[ipoly_bw];
            ipoly_bw=(ipoly_bw+npoly-1)%npoly;
            }
          ppocket[idx]=hull[ihull];
          // process pocket recursively
          ++npockets;
          comp=prepPolyHelper (vv, ppocket, comp, false);
          ihull=ihull_next;
          ipoly=ipoly_next;
          }
        } while (ihull!=0);
      }
    if (npockets>1) 
      comp.add(new MocQueryComponent(MocQueryOp.OR,npockets));
    if (npockets>0) 
      comp.add(new MocQueryComponent(MocQueryOp.NOT));

    if (!doLast)
      addHull[hull.length-1]=false;

    // add convex hull
    for (int i=0; i<hull.length; ++i)
      if (addHull[i])
        comp.add(new MocQueryComponent
          (vv[hull[i]].cross(vv[hull[(i+1)%hull.length]]).norm(),0.5*Math.PI));

    int num_and = 0;
    for (int i=0; i<hull.length; ++i)
      if (addHull[i]) ++num_and;
    if (npockets>0) ++num_and;
    if (num_and>1) 
      comp.add(new MocQueryComponent(MocQueryOp.AND,num_and));

    return comp;
    }

  static public ArrayList<MocQueryComponent> prepPolygon(ArrayList<Vec3> vertex)
    throws Exception
    {
    HealpixUtils.check(vertex.size()>=3,"not enough vertices in polygon");
    Vec3 vv[] = new Vec3[vertex.size()];
    for (int i=0; i<vertex.size(); ++i)
      vv[i]=vertex.get(i).norm();

    int[] P=new int[vv.length];
    for (int i=0; i<P.length; ++i)
      P[i]=i;
    ArrayList<MocQueryComponent> comp = new ArrayList<MocQueryComponent>();
    return prepPolyHelper(vv,P,comp,true);
    }

  static public Moc queryGeneralPolygon (ArrayList<Vec3> vertex, int order)
    throws Exception
    { return doMocQuery (order, MocQuery.prepPolygon(vertex)); }
  static public Moc queryGeneralPolygonInclusive (ArrayList<Vec3> vertex,
    int order, int omax)
    throws Exception
    { return doMocQueryInclusive (order, omax, MocQuery.prepPolygon(vertex)); }
  }
