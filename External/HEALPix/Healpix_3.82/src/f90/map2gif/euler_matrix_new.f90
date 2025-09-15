!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2013 Krzysztof M. Gorski, Eric Hivon,
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.sourceforge.net
!
!-----------------------------------------------------------------------------

subroutine euler_matrix_new(       &
                            a1,    &
                            a2,    &
                            a3,    &
                            matrix,&
                            itype   &
                           )

  ! --------------------------------------------------------------
  !
  ! this subroutine computes the Euler matrix for different prescription
  !          X: +0, Y: +1, Z: +2
  !
  !itype=1 :    rotation a1 around original Z
  !             rotation a2 around interm   Y   (+1)
  !             rotation a3 around final    X   (+0)
  !     aeronautics convention                   =1
  !
  !itype=2 :    rotation a1 around original Z
  !             rotation a2 around interm   X   (+0)
  !             rotation a3 around final    Z   (+2)
  !     classical mechanics convention           =2
  !
  !
  !itype=3 :    rotation a1 around original Z
  !             rotation a2 around interm   Y   (+1)
  !             rotation a3 around final    Z   (+2)
  !     quantum mechanics convention             =3
  !
  ! see H. Goldestein, Classical Mechanics (2nd Ed.) p. 147 for
  ! discussion
  ! INPUTS:
  !        a1, a2, a3 = Euler angles, in radians
  !         all the angles are measured counterclockwise
  !
  ! MODIFICATION HISTORY:
  !    October 1997, Eric Hivon, TAC (original IDL code)
  !    January 1999, A.J. Banday MPA (conversion to F90)
  !    March 22 1999, E.H. Caltech   (modification to match IDL version,
  !    correction of a bug on the matrix product)
  !    May 2002, E. H. Caltech (transformation of euler_matrix into euler_matrix_new)
  !       M_new(a,b,c)  =  M_old(-a, b,-c)
  ! --------------------------------------------------------------

   USE healpix_types

   IMPLICIT NONE

   REAL(DP)                , INTENT(IN)  :: a1,a2,a3
   REAL(DP), DIMENSION(3,3), INTENT(OUT) :: matrix
   INTEGER(I4B)            , INTENT(IN)  :: itype

   REAL(DP)                 :: c1,c2,c3, s1,s2,s3
   REAL(DP)                 :: ze,un
   REAL(DP), DIMENSION(3,3) :: m1,m2,m3
   character(len=*), parameter :: code ="Euler_Matrix_New"

   ! --------------------------------------------------------------

   !-- convert to sin and cosine values --
   c1 = COS(a1)
   s1 = SIN(a1)
   c2 = COS(a2)
   s2 = SIN(a2)
   c3 = COS(a3)
   s3 = SIN(a3)

   ze = 0.0_dp
   un = 1.0_dp

   !-- form full rotation matrix --
   !-- FOR REFERENCE:

   ! itype : 1
   ! m1 = [[ c1,-s1,  0],[ s1, c1,  0],[  0,  0,  1]] ; around   z
   ! m2 = [[ c2,  0, s2],[  0,  1,  0],[-s2,  0, c2]] ; around   y
   ! m3 = [[  1,  0,  0],[  0, c3,-s3],[  0, s3, c3]] ; around   x

   ! itype : 2
   ! m1 = [[ c1,-s1,  0],[ s1, c1,  0],[  0,  0,  1]] ; around   z
   ! m2 = [[  1,  0,  0],[  0, c2,-s2],[  0, s2, c2]] ; around   x
   ! m3 = [[ c3,-s3,  0],[ s3, c3,  0],[  0,  0,  1]] ; around   z

   ! itype : 3
   ! m1 = [[ c1,-s1,  0],[ s1, c1,  0],[  0,  0,  1]] ; around   z
   ! m2 = [[ c2,  0, s2],[  0,  1,  0],[-s2,  0, c2]] ; around   y
   ! m3 = [[ c3,-s3,  0],[ s3, c3,  0],[  0,  0,  1]] ; around   z


   ! matrix = m1 ( m2 m3 )

   m1(:,1) = (/ c1,-s1, ze /)  ! 1st column
   m1(:,2) = (/ s1, c1, ze /)  ! 2nd column
   m1(:,3) = (/ ze, ze, un /)  ! 3rd column

   select case (itype)

   case(1) ! ZYX

      m2(:,1) = (/ c2, ze, s2 /)
      m2(:,2) = (/ ze, un, ze /)
      m2(:,3) = (/-s2, ze, c2 /)

      m3(:,1) = (/ un, ze, ze /)
      m3(:,2) = (/ ze, c3,-s3 /)
      m3(:,3) = (/ ze, s3, c3 /)

   case(2) ! ZXZ

      m2(:,1) = (/ un, ze, ze /)
      m2(:,2) = (/ ze, c2,-s2 /)
      m2(:,3) = (/ ze, s2, c2 /)

      m3(:,1) = (/ c3,-s3, ze /)
      m3(:,2) = (/ s3, c3, ze /)
      m3(:,3) = (/ ze, ze, un /)

   case(3) ! ZYZ

      m2(:,1) = (/ c2, ze, s2 /)
      m2(:,2) = (/ ze, un, ze /)
      m2(:,3) = (/-s2, ze, c2 /)

      m3(:,1) = (/ c3,-s3, ze /)
      m3(:,2) = (/ s3, c3, ze /)
      m3(:,3) = (/ ze, ze, un /)

   case default

      print*,'Unknow rotation prescription in '//code
      print*,'Must be in [1,2,3]'
      print*,'is ',itype
      stop

   end select

   matrix = matmul(m1,matmul(m2,m3))

 end subroutine euler_matrix_new
