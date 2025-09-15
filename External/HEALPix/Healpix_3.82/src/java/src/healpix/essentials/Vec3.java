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

/** Cartesian 3-vector.
    Instead of using the javax.vecmath.Vector3d class, this separate class was
    implemented for two reasons: first, to avoid the external dependency from
    vecmath.jar, and also because the function Vector3d.angle(Vector3d v1) is
    too inaccurate for usage in Healpix for very small angles.

    @copyright (C) 2011-2015 Max-Planck-Society
    @author Martin Reinecke */
public class Vec3
  {
  public double x, y, z;

  /** Default constructor. */
  public Vec3() {}

  public Vec3 (Vec3 v)
    { x=v.x; y=v.y; z=v.z; }
  /** Creation from individual components */
  public Vec3 (double x1, double y1, double z1)
    { x=x1; y=y1; z=z1; }
  /** Conversion from {@link Pointing} */
  public Vec3 (Pointing ptg)
    {
    double sth = FastMath.sin(ptg.theta);
    x=sth*FastMath.cos(ptg.phi);
    y=sth*FastMath.sin(ptg.phi);
    z=FastMath.cos(ptg.theta);
    }
  /** Conversion from {@link Zphi} */
  public Vec3 (Zphi zphi)
    {
    double sth = Math.sqrt((1.0-zphi.z)*(1.0+zphi.z));
    x=sth*FastMath.cos(zphi.phi);
    y=sth*FastMath.sin(zphi.phi);
    z=zphi.z;
    }

  public Vec3 (double[] arr)
    {
    if(arr.length!=3) throw new IllegalArgumentException("Wrong array size");
    x = arr[0]; y = arr[1]; z = arr[2];
    }

  /** Vector length
      @return the length of the vector. */
  public final double length()
    { return Math.sqrt(lengthSquared()); }

  /** Squared vector length
      @return the squared length of the vector. */
  public final double lengthSquared()
    { return x*x + y*y + z*z; }

  /** Normalize the vector */
  public void normalize()
    {
    double d = 1./length();
    x *= d; y *= d; z *= d;
    }
  /** Return normalized vector */
  public Vec3 norm()
    {
    double d = 1./length();
    return new Vec3(x*d,y*d,z*d);
    }

  /** Angle between two vectors.
      @param v1 another vector
      @return the angle in radians between this vector and {@code v1};
        constrained to the range [0,PI]. */
  public final double angle(Vec3 v1)
    { return FastMath.atan2(cross(v1).length(), dot(v1)); }

  /** Vector cross product.
      @param v another vector
      @return the vector cross product between this vector and {@code v} */
  public Vec3 cross(Vec3 v)
    { return new Vec3(y*v.z - v.y*z, z*v.x - v.z*x, x*v.y - v.x*y); }

  /** Vector scaling.
    * @param n the scale number to be multiply to the coordinates {@code x,y,z}
    * @return the vector with coordinates multiplied by {@code n}
    */
  public Vec3 mul(double n)
    { return new Vec3(n*x, n*y, n*z); }

  /** Invert the signs of all components */
  public void flip()
    { x=-x; y=-y; z=-z; }
  public Vec3 flipped()
    { return new Vec3(-x,-y,-z); }

  /** Scale the vector by a given factor
      @param n the scale factor */
  public void scale(double n)
    { x*=n; y*=n; z*=n; }

  /** Computes the dot product of the this vector and {@code v1}.
    * @param v1 another vector
    * @return dot product */
  public final double dot(Vec3 v1)
    { return x*v1.x + y*v1.y + z*v1.z; }

  /** Vector addition
    * @param v the vector to be added
    * @return addition result */
  public Vec3 add(Vec3 v)
    { return new Vec3(x+v.x, y+v.y, z+v.z); }

  /** Vector subtraction
    * @param v the vector to be subtracted
    * @return subtraction result */
  public Vec3 sub(Vec3 v)
    { return new Vec3(x-v.x, y-v.y, z-v.z); }

  public String toString()
    {
    StringBuilder s = new StringBuilder();
    s.append("vec3(");s.append(x);
    s.append(",");s.append(y);
    s.append(",");s.append(z);
    s.append(")");
    return s.toString();
    }

  public double[] toArray()
    { return new double[]{x,y,z}; }

  public void toArray(double[] arr)
    {
    if(arr.length!=3) throw new IllegalArgumentException("wrong array size");
    arr[0] = x; arr[1] = y; arr[2] = z;
    }

  public boolean equals(Object o)
    {
    if (this == o) return true;
    if (o == null || getClass() != o.getClass()) return false;
    Vec3 vec3 = (Vec3) o;
    if (Double.compare(vec3.x, x) != 0) return false;
    if (Double.compare(vec3.y, y) != 0) return false;
    if (Double.compare(vec3.z, z) != 0) return false;
    return true;
    }

  public int hashCode()
    {
    int result = Double.valueOf(x).hashCode();
    result = 31 * result + Double.valueOf(y).hashCode();
    result = 31 * result + Double.valueOf(z).hashCode();
    return result;
    }
  }
