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

/** Finds the smallest enclosing cone for a point set on the sphere according to
    Barequet & Elber: Information Processing Letters 93(2005), p.83.
    All points are expected to be passed as unit vectors.
    The enclosing cone must have an opening angle smaller than pi/2.

    @copyright 2011 Max-Planck-Society
    @author Martin Reinecke */
final class CircleFinder
  {
  private Vec3 center;
  private double cosrad;

  private void getCircle (Vec3[] point, int q1, int q2)
    {
    center = point[q1].add(point[q2]); center.normalize();
    cosrad = point[q1].dot(center);
    for (int i=0; i<q1; ++i)
      if (point[i].dot(center)<cosrad) // point outside the current circle
        {
        center=(point[q1].sub(point[i])).cross(point[q2].sub(point[i]));
        center.normalize();
        cosrad=point[i].dot(center);
        if (cosrad<0)
          { center.flip(); cosrad=-cosrad; }
        }
    }
  private void getCircle (Vec3[] point, int q)
    {
    center = point[0].add(point[q]); center.normalize();
    cosrad = point[0].dot(center);
    for (int i=1; i<q; ++i)
      if (point[i].dot(center)<cosrad) // point outside the current circle
        getCircle(point,i,q);
    }

  public CircleFinder (Vec3[] point) throws Exception
    {
    int np=point.length;
    HealpixUtils.check(np>=2,"too few points");
    center = point[0].add(point[1]); center.normalize();
    cosrad = point[0].dot(center);
    for (int i=2; i<np; ++i)
      if (point[i].dot(center)<cosrad) // point outside the current circle
        getCircle(point,i);
    }

  public Vec3 getCenter() { return new Vec3(center); }
  public double getCosrad() { return cosrad; }
  }
