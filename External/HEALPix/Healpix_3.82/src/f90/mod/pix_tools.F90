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
module pix_tools
!==================================================================
!    Various subroutines for converting angles and unit vectors
!    to pixel numbers in both the NESTED and RING schemes,
!    as well as looking for neighbours and making pixel queries
!
!    last update : Oct 4, 2002, EH, computation of pixel vertex in
!        pix2vec_ring and pix2vec_nest
!    2004-05-28, EH: edition of next_in_line_nest for nside=1
!                    addition of optional 'mask' in remove_dipole
!    2004-06-01      correction of bug at low phi in in_ring
!    2004-08-25      added multi-D convert_inplace_* and convert_nest2ring, convert_ring2nest
!    2004-10-07      added template_pixel_*, same_shape_pixels_*
!    2004-12-15      generic forms for convert_inplace (I,R,D with 1 or N dimensions), 
!                     remove_dipole (R,D)
!    2005-01-25      ensure backward compatibility of remove_dipole
!    2005-08-03      edit same_shape_pixels_nest to comply with GFORTRAN
!    2006-06-28      edit to neighbours_nest to correct error in Nside=1 case
!    2008-02-03:     addition of weights in remove dipole
!    2008-02-04:     check ordering in [1,2] in remove dipole ; check nest in [0,1] in query_*
!    2008-03-28:     edit ring_num and query_strip to improve inclusive query_strip
!                    bug fix on query_disc
!    2009-03-04:     starts supporting 64-bit npix
!              largest supported Nside is now 2^28 instead of 2^13
!             Note: the theoretical limit Nside=2^29 is not implemented to avoid
!               dealing with too many 64 bit integer variables that slow down execution
!    2011-08-25 : 2011-08-*: improves accuracy of pixel routines close to Poles
!    2011-10-13: improved query_triangle, introduced process_intervals
!    2012-07-17: Parallel OpenMP implementation of medfiltmap
!    2012-08-27: correction of a bug affecting neighbours_nest and next_in_line_nest at Nside>8192
!    2013-04-02: bug correction in query_disc in inclusive mode
!    2013-05-07: G95-compatible
!    2015-09-02: added nest2uniq, uniq2nest
!    2018-05-22: added nside2npweights
!    2018-10-22: added apply_mask
!==================================================================
  ! subroutine query_strip                          Done (To be Tested) depends on in_ring
  ! subroutine query_polygon                        Done (To be Tested) depends on isort
  ! subroutine query_triangle                       Done (To be Tested)
  ! subroutine getdisc_ring                         Done depends on query_disc
  ! subroutine query_disc                           Done (To be Tested)
  ! function ring_num                                     OK
  ! function ring2z                                 Done  
  ! subroutine in_ring                              Done (To be Tested) depends in next_in_line_nest
  ! subroutine intrs_intrv                                OK
  ! subroutine process_intervals               ok
  ! function fudge_query_radius                ok
  !!!!!! subroutine correct_ring_phi                           OK
  !                                         
  ! subroutine pix2ang_ring                          Done  
  ! subroutine pix2vec_ring                          Done  
  ! subroutine ang2pix_ring                          Done  
  ! subroutine vec2pix_ring                          Done  
  ! subroutine ang2pix_nest                          Done  
  ! subroutine vec2pix_nest                          Done  
  ! subroutine nest2ring                             Done  
  ! subroutine ring2nest                             Done  
  ! function npix2nside                              Done  
  !                                         
  ! subroutine warning_oldbounds                           OK
  ! subroutine remove_dipole_real                    Done depends on query_strip
  ! subroutine remove_dipole_double                  Done
  ! subroutine remove_dipole_real_old                Done
  ! subroutine remove_dipole_double_old              Done
  ! subroutine remove_dipole_real_v12                Done
  ! subroutine remove_dipole_double_v12              Done
  ! subroutine add_dipole_real                      Done  
  ! subroutine add_dipole_double                    Done  
  ! subroutine medfiltmap_S                         Done (To be tested)
  ! subroutine medfiltmap_D                         Done (To be tested)

  ! subroutine convert_nest2ring_int_1d             Done
  ! subroutine convert_nest2ring_int8_1d            Done
  ! subroutine convert_nest2ring_real_1d            Done
  ! subroutine convert_nest2ring_double_1d          Done      
  ! subroutine convert_ring2nest_int_1d             Done
  ! subroutine convert_ring2nest_int8_1d            Done
  ! subroutine convert_ring2nest_real_1d            Done
  ! subroutine convert_ring2nest_double_1d          Done
  ! subroutine convert_nest2ring_int_nd             Done
  ! subroutine convert_nest2ring_int8_nd            Done
  ! subroutine convert_nest2ring_real_nd            Done
  ! subroutine convert_nest2ring_double_nd          Done      
  ! subroutine convert_ring2nest_int_nd             Done
  ! subroutine convert_ring2nest_int8_nd            Done
  ! subroutine convert_ring2nest_real_nd            Done
  ! subroutine convert_ring2nest_double_nd          Done
  !
  ! subroutine convert_inplace_int_1d               Done (to be tested)
  ! subroutine convert_inplace_real_1d              Done (to be tested)
  ! subroutine convert_inplace_double_1d            Done (to be tested)
  ! subroutine convert_inplace_int_nd               Done (to be tested)
  ! subroutine convert_inplace_real_nd              Done (to be tested)
  ! subroutine convert_inplace_double_nd            Done (to be tested)
  !
  ! subroutine xy2pix_nest                          Done (to be tested)               
  ! subroutine pix2xy_nest                          Done (to be tested)
  ! subroutine mk_pix2xy                                  OK
  ! subroutine mk_xy2pix                                  OK
  ! subroutine mk_xy2pix1                                 OK
  ! subroutine neighbours_nest                      Done (to be tested)               
  ! subroutine next_in_line_nest                          Done  depends on xy2pix_nest pix2xy_nest             
  ! subroutine ang2vec                                    OK
  ! subroutine vec2ang                                    OK
  ! function nside2npix                            Done  
  ! subroutine surface_triangle                           OK
  ! subroutine angdist                                    OK
  ! subroutine vect_prod                                  OK
  ! function nside2ntemplates                       Done
  ! subroutine template_pixel_ring                  Done (to be tested)               
  ! subroutine same_shape_pixels_ring               Done (to be tested)               
  ! subroutine template_pixel_nest                  Done (to be tested)               
  ! subroutine same_shape_pixels_nest               Done (to be tested)
! 
! 2011-08-*: improves accuracy close to Poles
! ang2pix_nest  small
! ang2pix_ring  small
! ang2vec     OK      
! angdist     Done  
! pix2ang_nest  small, kshift
! pix2ang_ring  isqrt, small 
! pix2vec_nest  small,  [kshift TODO]
! pix2vec_ring  isqrt, small
! cheap_isqrt Added     
! vec2ang     Done
! vec2pix_nest  small
! vec2pix_ring  small
! ring2nest     isqrt

  USE healpix_types
  USE misc_utils
  USE long_intrinsic, only: long_size
  IMPLICIT none

  INTEGER(KIND=i4b), private, PARAMETER :: ns_max4=8192     ! 2^13
  INTEGER(KIND=i4b), private, PARAMETER :: ns_max8=268435456! 2^28
#ifdef NO64BITS
  INTEGER(KIND=i4b), private, PARAMETER :: ns_max=ns_max4 ! largest nside available
#else
  INTEGER(KIND=i4b), private, PARAMETER :: ns_max=ns_max8 ! largest nside available
#endif

  !initialise array x2pix, y2pix and pix2x, pix2y used in several routines
!  integer(KIND=i4b), private, save, dimension(128) :: x2pix=0,y2pix=0
  integer(KIND=i4b), private, save, dimension(0:127) :: x2pix1=-1,y2pix1=-1

  integer(KIND=i4b), private, save, dimension(0:1023) :: pix2x=-1, pix2y=-1

  ! coordinate of the lowest corner of each face
  INTEGER(KIND=I4B), private, save, dimension(0:11) :: jrll1 = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
  INTEGER(KIND=I4B), private, save, dimension(0:11) :: jpll1 = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2

  ! obsolete
  interface convert_inplace_real
     module procedure convert_inplace_real_1d, convert_inplace_real_nd
  end interface 
  interface convert_inplace_int
     module procedure convert_inplace_int_1d, convert_inplace_int_nd
  end interface 

  ! generic form
  interface convert_inplace
     module procedure convert_inplace_real_1d, convert_inplace_real_nd, &
          &           convert_inplace_int_1d, convert_inplace_int_nd, &
          &           convert_inplace_double_1d, convert_inplace_double_nd
  end interface 

  interface remove_dipole
     module procedure remove_dipole_real, remove_dipole_double, &
          &           remove_dipole_real_old, remove_dipole_double_old, &
          &           remove_dipole_real_v12, remove_dipole_double_v12
  end interface

  interface add_dipole
     module procedure add_dipole_real, add_dipole_double
  end interface

  interface apply_mask
     module procedure apply_mask_real, apply_mask_double
  end interface

  interface convert_nest2ring
     module procedure convert_nest2ring_int_1d, &
          &           convert_nest2ring_real_1d, &
          &           convert_nest2ring_double_1d, &
          &           convert_nest2ring_int_nd, &
          &           convert_nest2ring_real_nd, &
          &           convert_nest2ring_double_nd
#ifndef NO64BITS
     module procedure convert_nest2ring_int8_1d, &
          convert_nest2ring_int8_nd
#endif
  end interface 

  interface convert_ring2nest
     module procedure convert_ring2nest_int_1d, &
          &           convert_ring2nest_real_1d, &
          &           convert_ring2nest_double_1d, &
          &           convert_ring2nest_int_nd, &
          &           convert_ring2nest_real_nd, &
          &           convert_ring2nest_double_nd
#ifndef NO64BITS
     module procedure convert_ring2nest_int8_1d, &
          convert_ring2nest_int8_nd
#endif
  end interface 

  interface medfiltmap
     module procedure medfiltmap_s, medfiltmap_d
  end interface 

#ifdef NO64BITS
  interface npix2nside
     module procedure npix2nside
  end interface
  interface ang2pix_nest
     module procedure ang2pix_nest
  end interface
  interface vec2pix_nest
     module procedure vec2pix_nest
  end interface
  interface ang2pix_ring
     module procedure ang2pix_ring
  end interface
  interface vec2pix_ring
     module procedure vec2pix_ring
  end interface
  interface ring2nest
     module procedure ring2nest
  end interface
  interface nest2ring
     module procedure nest2ring
  end interface
  interface pix2ang_ring
     module procedure pix2ang_ring
  end interface
  interface pix2vec_ring
     module procedure pix2vec_ring
  end interface
  interface pix2ang_nest
     module procedure pix2ang_nest
  end interface
  interface pix2vec_nest
     module procedure pix2vec_nest
  end interface
  interface in_ring
     module procedure in_ring
  end interface
  interface getdisc_ring
     module procedure getdisc_ring
  end interface
  interface query_strip
     module procedure query_strip
  end interface
  interface query_disc_old ! <<<<<<<<<<<<<<<<<<<<<<
     module procedure query_disc_old
  end interface
  interface query_disc
     module procedure query_disc
  end interface
  interface query_polygon
     module procedure query_polygon
  end interface
  interface query_triangle
     module procedure query_triangle
  end interface
  interface xy2pix_nest
     module procedure xy2pix_nest
  end interface
  interface pix2xy_nest
     module procedure pix2xy_nest
  end interface
  interface neighbours_nest
     module procedure neighbours_nest
  end interface
  interface next_in_line_nest
     module procedure next_in_line_nest
  end interface
  interface template_pixel_ring
     module procedure template_pixel_ring
  end interface
  interface same_shape_pixels_ring
     module procedure same_shape_pixels_ring
  end interface
  interface template_pixel_nest
     module procedure template_pixel_nest
  end interface
  interface same_shape_pixels_nest
     module procedure same_shape_pixels_nest
  end interface
  interface cheap_isqrt
     module procedure cheap_isqrt_4
  end interface
  interface discedge2fulldisc
     module procedure discedge2fulldisc
  end interface
  interface nest2uniq
     module procedure nest2uniq
  end interface
  interface uniq2nest
     module procedure uniq2nest
  end interface
#else
  interface npix2nside
     module procedure npix2nside, npix2nside_8
  end interface
  interface ang2pix_nest
     module procedure ang2pix_nest, ang2pix_nest_8
  end interface
  interface vec2pix_nest
     module procedure vec2pix_nest, vec2pix_nest_8
  end interface
  interface ang2pix_ring
     module procedure ang2pix_ring, ang2pix_ring_8
  end interface
  interface vec2pix_ring
     module procedure vec2pix_ring, vec2pix_ring_8
  end interface
  interface ring2nest
     module procedure ring2nest, ring2nest_8
  end interface
  interface nest2ring
     module procedure nest2ring, nest2ring_8
  end interface
  interface pix2ang_ring
     module procedure pix2ang_ring, pix2ang_ring_8
  end interface
  interface pix2vec_ring
     module procedure pix2vec_ring, pix2vec_ring_8
  end interface
  interface pix2ang_nest
     module procedure pix2ang_nest, pix2ang_nest_8
  end interface
  interface pix2vec_nest
     module procedure pix2vec_nest, pix2vec_nest_8
  end interface
  interface in_ring
     module procedure in_ring, in_ring_8
  end interface
  interface getdisc_ring
     module procedure getdisc_ring, getdisc_ring_8
  end interface
  interface query_strip
     module procedure query_strip, query_strip_8
  end interface
  interface query_disc_old ! <<<<<<<<<<<<<<<<<<<<<<
     module procedure query_disc_old, query_disc_old_8
  end interface
  interface query_disc
     module procedure query_disc, query_disc_8
  end interface
  interface query_polygon
     module procedure query_polygon, query_polygon_8
  end interface
  interface query_triangle
     module procedure query_triangle, query_triangle_8
  end interface
  interface xy2pix_nest
     module procedure xy2pix_nest, xy2pix_nest_8
  end interface
  interface pix2xy_nest
     module procedure pix2xy_nest, pix2xy_nest_8
  end interface
  interface neighbours_nest
     module procedure neighbours_nest, neighbours_nest_8
  end interface
  interface next_in_line_nest
     module procedure next_in_line_nest, next_in_line_nest_8
  end interface
  interface template_pixel_ring
     module procedure template_pixel_ring, template_pixel_ring_8
  end interface
  interface same_shape_pixels_ring
     module procedure same_shape_pixels_ring, same_shape_pixels_ring_8
  end interface
  interface template_pixel_nest
     module procedure template_pixel_nest, template_pixel_nest_8
  end interface
  interface same_shape_pixels_nest
     module procedure same_shape_pixels_nest, same_shape_pixels_nest_8
  end interface
  interface cheap_isqrt
     module procedure cheap_isqrt_4, cheap_isqrt_8
  end interface
  interface discedge2fulldisc
     module procedure discedge2fulldisc, discedge2fulldisc_8
  end interface
  interface nest2uniq
     module procedure nest2uniq, nest2uniq_8
  end interface
  interface uniq2nest
     module procedure uniq2nest, uniq2nest_8
  end interface
#endif

  private

  public :: remove_dipole, add_dipole, &
       & query_strip, &
       & query_polygon, &
       & query_triangle, &
       & query_disc, &
       & pix2ang_ring, pix2vec_ring, ang2pix_ring, vec2pix_ring, &
       & pix2ang_nest, pix2vec_nest, ang2pix_nest, vec2pix_nest, &
       & convert_nest2ring, convert_ring2nest, &
       & convert_inplace, convert_inplace_real, convert_inplace_int, &
       & nest2ring, ring2nest, xy2pix_nest, pix2xy_nest, &
!!!       & mk_pix2xy, mk_xy2pix, &
       & neighbours_nest, &
       & next_in_line_nest, &
       & ang2vec, vec2ang, &
       & npix2nside, nside2npix, &
       & surface_triangle, angdist, vect_prod

!   public :: ang2pix_nest_old, vec2pix_nest_old, &
!        ang2pix_ring_old, vec2pix_ring_old, &
!        ring2nest_old, nest2ring_old, &
!        pix2ang_ring_old, pix2vec_ring_old, &
!        pix2ang_nest_old, pix2vec_nest_old, &
!        template_pixel_ring_old, same_shape_pixels_ring_old, &
!        template_pixel_nest_old, same_shape_pixels_nest_old

  public :: query_disc_old

  public :: nside2ntemplates, &
       & template_pixel_ring, same_shape_pixels_ring, &
       & template_pixel_nest, same_shape_pixels_nest

  public :: medfiltmap

  public :: intrs_intrv, in_ring, ring_num, ring2z ! arcane usage
  public :: cheap_isqrt
  public :: process_intervals, fudge_query_radius, discphirange_at_z

  public :: getdisc_ring  ! obsolete

  public :: nest2uniq, uniq2nest

  public :: nside2npweights

  public :: apply_mask


contains


!   !=======================================================================
!   function discphirange_at_z (vcenter, radius, z, phi0) result(dphi)
!     !=======================================================================
!     !  for a disc centered on  vcenter and given radius,
!     !  and for location z=cos(theta)
!     !  returns dphi such that the disc has   abs(phi-phi0) <= dphi
!     !
!     ! solving disc equation on sphere:
!     !  sin(theta) sin(theta0) cos(phi-phi0) + cos(theta) cos(theta0) >= cos(radius)
!     !
!     ! 2011-06-21: adapted from IDL routine of same name
!     !=======================================================================
!     real(DP), dimension(1:), intent(in)      :: vcenter
!     real(DP),                intent(in)      :: radius, z
!     real(DP),                intent(out)     :: phi0
!     real(DP) :: dphi 
!     !
!     real(DP) :: cosang, cosdphi, norm, x0, y0, z0, st0, diff, st
!     real(DP), parameter :: zero=0.0_dp, one=1.0_dp

!     cosang = cos(radius)
!     norm = sqrt(sum(vcenter(1:3)**2))
!     x0 = vcenter(1) / norm
!     y0 = vcenter(2) / norm
!     z0 = vcenter(3) / norm

!     phi0 = zero
!     if ((x0 /= zero) .or. (y0 /= zero)) phi0 = atan2(y0, x0)

!     st0 = x0*x0 + y0*y0 ! sin(theta0)^2
!     diff = cosang - z*z0 ! cos(rad) - cos(theta)cos(theta0)
    
!     dphi = -1000.0_dp ! out of disc
!     if (st0 == zero) then ! polar cap
!        if (diff <= zero) dphi = PI ! fully in cap
!     else
!        st = max(one - z*z, 1.0e-12_dp)! sin(theta)^2
!        cosdphi = diff/sqrt(st0*st)
!        if (cosdphi < -one) dphi = PI ! fully in disc
!        if (abs(cosdphi) <= one) dphi = acos(cosdphi)
!     endif

!   end function discphirange_at_z
  !=======================================================================
  function ring_num (nside, z, shift) result(ring_num_result)
    !=======================================================================
    ! ring = ring_num(nside, z [, shift=])
    !     returns the ring number in {1, 4*nside-1}
    !     from the z coordinate
    ! usually returns the ring closest to the z provided
    ! if shift < 0, returns the ring immediatly north (of smaller index) of z
    ! if shift > 0, returns the ring immediatly south (of smaller index) of z
    !
    !=======================================================================
    INTEGER(KIND=I4B)             :: ring_num_result
    REAL(KIND=DP),     INTENT(IN) :: z
    INTEGER(KIND=I4B), INTENT(IN) :: nside
    integer(i4b),      intent(in), optional :: shift

    INTEGER(KIND=I4B) :: iring
    real(DP) :: my_shift
    !=======================================================================


    my_shift = 0.0_dp
    if (present(shift)) my_shift = shift * 0.5_dp

    !     ----- equatorial regime ---------
    iring = NINT( nside*(2.0_dp-1.500_dp*z)   + my_shift )

    !     ----- north cap ------
    if (z > twothird) then
       iring = NINT( nside* SQRT(3.0_dp*(1.0_dp-z))  + my_shift )
       if (iring == 0) iring = 1
    endif

    !     ----- south cap -----
    if (z < -twothird   ) then
       ! beware that we do a -shift in the south cap
       iring = NINT( nside* SQRT(3.0_dp*(1.0_dp+z))   - my_shift )
       if (iring == 0) iring = 1
       iring = 4*nside - iring
    endif

    ring_num_result = iring

    return
  end function ring_num
  !=======================================================================
  function ring2z (nside, ir) result(z)
    !=======================================================================
    !     returns the z coordinate of ring ir for Nside
    ! 2009-03-25: accepts Nside > 8192
    !=======================================================================
    integer(kind=I4B), intent(in) :: ir
    integer(kind=I4B), intent(in) :: nside
    real(kind=DP)                 :: z

    real(DP)     :: fn, tmp
    !=======================================================================

    fn = real(nside,kind=DP)

    if (ir < nside) then  ! polar cap (north)
       tmp = real(ir, DP)
       z = 1.0_dp - (tmp * tmp) / (3.0_dp * fn * fn)
    else if (ir < 3*nside) then ! tropical band
       z = real( 2*nside-ir, kind=DP) * 2.0_dp / (3.0_dp * fn)
    else                  ! polar cap (south)
       tmp = real(4*nside - ir, DP)
       z = -1.0_dp + (tmp * tmp) / (3.0_dp * fn * fn)
    endif
    
    return
  end function ring2z
  !=======================================================================
  subroutine intrs_intrv( d1, d2, di, ni)
    !=======================================================================
    ! computes the intersection di
    ! of 2 intervals d1 (= [a1,b1]) and d2 (= [a2,b2])
    ! on the periodic domain ( = [A,B], where A and B are arbitrary)
    ! ni is the resulting number of intervals (0,1, or 2)
    !
    ! if a1<b1 then d1 = {x | a1 <= x <= b1}
    ! if a1>b1 then d1 = {x | a1 <= x <= B  U  A <= x <= b1}
    !=======================================================================
    real(kind=DP), dimension(1:), INTENT(IN)  :: d1, d2
    real(kind=DP), dimension(1:), INTENT(OUT) :: di
    integer(kind=I4B), INTENT(OUT) :: ni

    real(kind=DP), dimension(1:4) :: dk
    integer(kind=I4B) :: ik
    logical(kind=LGT) :: tr12, tr21, tr34, tr43, tr13, tr31, tr24, tr42, tr14, tr32
    !=======================================================================

    tr12 = (d1(1) < d1(2))
    tr21 = .NOT. tr12
    tr34 = (d2(1) < d2(2))
    tr43 = .NOT. tr34
    tr13 = (d1(1) < d2(1))
    tr31 = .NOT. tr13
    tr24 = (d1(2) < d2(2))
    tr42 = .NOT. tr24
    tr14 = (d1(1) < d2(2))
    tr32 = (d2(1) < d1(2))

    ik = 0
    dk(1:4) = -1.0e9_dp


    if ((tr31.AND.tr14) .OR. (tr43.AND.(tr31.OR.tr14))) then
       ik = ik + 1
       dk(ik) = d1(1)  ! a1
    endif
    if ((tr13.AND.tr32) .OR. (tr21.AND.(tr13.OR.tr32))) then
       ik = ik + 1
       dk(ik) = d2(1)  ! a2
    endif
    if ((tr32.AND.tr24) .OR. (tr43.AND.(tr32.OR.tr24))) then
       ik = ik + 1
       dk(ik) = d1(2)  ! b1
    endif
    if ((tr14.AND.tr42) .OR. (tr21.AND.(tr14.OR.tr42))) then
       ik = ik + 1
       dk(ik) =  d2(2)  ! b2
    endif

    di(1:4) = 0.0_dp
    select case (ik)
    case (0)
       ni = 0
    case (2)
       ni = 1
       di(1:2) = (/ dk(1), dk(2) /) ! [a1,b1] or [a1,b2] or [a2,b1] or [a2,b2]
    case (4)
       ni = 2
       di(1:4) = (/ dk(1), dk(4), dk(2), dk(3) /) ! [a1,b2] U [a2,b1]
    case default
       print*,"error in intrs_intrv", ik
       print*,dk
       print*,d1,d2
    end select

    return
  end subroutine intrs_intrv

  !=======================================================================
  subroutine process_intervals(interval1, interval_list, interval_out, n_out)
  !=======================================================================
    ! intersection of 1 interval (defined by its 2 boundaries)
    ! with an arbitrary list of intervals (defined as n*2 boundaries)
    !=======================================================================
    real(dp), dimension(1:), intent(in) :: interval1
    real(dp), dimension(1:), intent(in) :: interval_list
    real(dp), dimension(1:), intent(out):: interval_out
    integer(i4b),          intent(out) :: n_out
    !
    integer(i4b) :: n_in, i_in, n_tmp_out
    real(dp), dimension(1:4) :: w4
    real(dp), dimension(1:30) :: work
    !=======================================================================

    n_out = 0
    interval_out = 0.d0
    n_in = size(interval_list)/2
    work(1:2*n_in) = interval_list(1:2*n_in) ! copy input array to avoid overwritting it

    if (n_in > 0) then
       do i_in=0, n_in-1
          ! intersection of 2 intervals -> 0,1,2 intervals
          call intrs_intrv(interval1, work(1+2*i_in:2+2*i_in), w4, n_tmp_out)
          if (n_tmp_out > 0) then ! n_tmp_out = 1 or 2
             interval_out(2*n_out+1: 2*(n_out+n_tmp_out)) = w4(1:2*n_tmp_out)
             n_out = n_out + n_tmp_out
          endif
       enddo
    endif
    
    return
  end subroutine process_intervals
  !=======================================================================
  function fudge_query_radius(nside, radius_in, quadratic) result(radius_out)
  !=======================================================================
    !  radius_out = fudge_query_radius( nside, radius_in, QUADRATIC=)
    ! 
    !  with 
    !    radius_out = radius_in + fudge            (default)
    !  or
    !    radius_out = sqrt(radius_in^2 + fudge^2)  (quadratic)
    ! 
    !  if absent, radius_in = 0
    ! where
    !    fudge = factor(nside) * Pi / (4 *Nside)
    !  with factor = sqrt( 5/9 + (8/Pi)^2 / 5 * (1-1/(2*Nside)) )
    ! 
    !  an upper bound of the actual largest radius 
    !  (determined by pixel located at z=2/3 and its North corner).
    !  
    ! 
    !  2011-10-14, EH, v1.0
    !=======================================================================
    integer(i4b), intent(in)           :: nside
    real(dp),     intent(in), optional :: radius_in
    logical(lgt), intent(in), optional :: quadratic
    real(dp)                           :: radius_out
    logical(LGT)                       :: do_quad
    real(dp)                           :: factor, fudge
    !=======================================================================
    radius_out = 0.0_dp
    if (present(radius_in)) radius_out = radius_in
    do_quad = .false.
    if (present(quadratic)) do_quad = quadratic

    factor = sqrt( 0.55555555555555555555_dp + 1.29691115062192347448165712908_dp*(1.0_dp-0.5_dp/nside) )
    fudge  = factor * PI / (4.0_dp * nside) ! increase effective radius

    if (do_quad) then 
       radius_out = sqrt(radius_out**2 + fudge**2)
    else
       radius_out = radius_out + fudge
    endif
    radius_out = min(radius_out, PI)

    return
  end function fudge_query_radius
  !=======================================================================
  subroutine discphirange_at_z(vector0, radius, z, nz, dphi, phi0)
    !=======================================================================
    !  for a disc centered on  vcenter and given radius,
    !  and for location z=cos(theta)
    !  returns dphi such that the disc has   abs(phi-phi0) <= dphi
    !
    ! solving disc equation on sphere:
    !  sin(theta) sin(theta0) cos(phi-phi0) + cos(theta) cos(theta0) >= cos(radius)
    !
    !=======================================================================
    real(dp), dimension(1:3), intent(in) :: vector0
    real(dp),                 intent(in) :: radius
    real(dp), dimension(1:),  intent(in) :: z
    integer(i4b),             intent(in) :: nz
    real(dp), dimension(1:),  intent(out) :: dphi
    real(dp),                 intent(out) :: phi0
    !
    real(dp) :: cosang, cosphi0, x0, y0, z0, a, b, c, cosdphi, norm
    integer(i4b) ::  i
    real(DP), parameter :: zero=0.0_dp, one=1.0_dp
    
    cosang = cos(radius)
    norm =  sqrt(sum(vector0(1:3)**2))
    x0 = vector0(1) / norm
    y0 = vector0(2) / norm
    z0 = vector0(3) / norm

    phi0 = zero
    if ((x0 /= zero).or.(y0 /= zero)) phi0 = atan2(y0, x0)  ! in ]-Pi, Pi]
    !cosphi0 = cos(phi0)
    a = x0*x0 + y0*y0   ! sin(theta0)^2  

    do i=1, nz
       dphi(i) = -1000.0_dp ! default value for out of disc

       b = cosang - z(i)*z0   ! cos(rad) - cos(theta)cos(theta0)
       if (a == zero) then ! poles
          if (b <= zero) dphi(i) = PI
       else
          c = max(one - z(i)*z(i), 1.0e-12_dp) ! sin(theta)^2
          cosdphi = b / sqrt(a*c)
          if (cosdphi      < -one) dphi(i) = PI ! all the pixels at this elevation are in the disc
          if (abs(cosdphi) <= one) dphi(i) = acos(cosdphi) ! in [0,Pi]
       endif

    enddo

    return
  end subroutine discphirange_at_z
  !=======================================================================
  subroutine pixels_on_edge(nside, irmin, irmax, phi0, dphi, ringphi, ngr)
  !=======================================================================
  !=======================================================================
    integer(i4b),                intent(in) :: nside, irmin, irmax
    real(dp),                    intent(in) :: phi0
    real(dp),     dimension(1:), intent(in) :: dphi
    integer(i4b), dimension(1:,1:), intent(out) :: ringphi
    integer(i4b),                   intent(out) :: ngr
    integer(i4b), parameter :: badvalue = -9999
    integer(i4b) :: nrings, ir, k, thisring, npr
    real(dp)    , parameter :: zero = 0.0_dp
    integer(i4b) :: iphi_low, iphi_hi, kshift
    real(dp) :: shift

    nrings = irmax - irmin + 1

    ngr = 0
    do thisring = irmin, irmax
       ir = thisring - irmin + 1
       call pixels_per_ring(nside, thisring, npr, kshift)
       if (dphi(ir) >= PI) then ! full ring
          ngr = ngr + 1
          ringphi(1, ngr) = thisring
          ringphi(2, ngr) = 0
          ringphi(3, ngr) = npr-1
       elseif (dphi(ir) >= zero) then ! partial ring
          shift = kshift * 0.5_dp
          iphi_low = ceiling (npr * (phi0 - dphi(ir)) / TWOPI - shift)
          iphi_hi  = floor   (npr * (phi0 + dphi(ir)) / TWOPI - shift)
          if (iphi_hi >= iphi_low) then ! pixel center in range
             ngr = ngr + 1
             ringphi(1, ngr) = thisring
             ringphi(2, ngr) = modulo(iphi_low, npr)
             ringphi(3, ngr) = modulo(iphi_hi,  npr)
          endif
       endif
    enddo

  end subroutine pixels_on_edge
  !=======================================================================
  subroutine pixels_per_ring(nside, ring, npr, kshift, npnorth)
  !=======================================================================
    ! for a given Nside and ring index in [1,4*Nside-1], 
    ! returns the number of pixels in ring, their shift (0 or 1) in azimuth
    ! and the number of pixels in current ring and above (=North)
    !
    ! NB: 'rings' 0 and 4*Nside respectively are North and South Poles
    !=======================================================================
    integer(i4b), intent(in) :: nside, ring
    integer(i4b), intent(out) :: npr, kshift
    integer(i8b), intent(out), optional :: npnorth
    integer(i8b) :: ncap, npix, ir
    
    ! number of pixels in current ring
    npr = min(nside, ring, 4*nside-ring) * 4
    ! shift
    kshift = mod(ring + 1, 2) ! 1 for even, 0 for odd
    if (nside == 1) kshift = 1 - kshift ! except for Nside=1
    if (npr < 4*nside) kshift = 1 ! 1 on polar cap
    ! Number of pixels in current ring and above
    if (present(npnorth)) then
       if (ring <= nside) then ! in North cap
          npnorth = ring*(ring+1_i8b)*2_i8b
       elseif (ring <= 3*nside) then ! in Equatorial region
          ncap = nside*(nside+1_i8b)*2_i8b
          ir = ring-nside
          npnorth = ncap + 4_i8b*nside*ir
       else ! in South cap
          npix = nside2npix(nside)
          ir = 4_i8b*nside-ring - 1 ! count ring from south
          npnorth = npix - ir*(ir+1_i8b)*2_i8b
       endif
    endif
    return
  end subroutine pixels_per_ring
  !=======================================================================
  subroutine check_edge_pixels(nside, nsboost, irmin, irmax, phi0, dphi, ringphi, ngr)
  !=======================================================================
    integer(i4b), intent(in) :: nside, nsboost, irmin, irmax
    real(dp),                       intent(in) :: phi0
    real(dp),     dimension(1:),    intent(in) :: dphi
    integer(i4b), dimension(1:,1:), intent(inout) :: ringphi
    integer(i4b),                   intent(inout) :: ngr

    !integer(i4b) :: i, j, k, kk, ngr_out, diff, iphi, i0, nrh
    integer(i4b) :: i, j, k, kk, ngr_out, iphi, i0, nrh, nshrink, np
    real(dp), dimension(1:2*nsboost+1) :: phiw, phie
    real(dp) :: dd, dph, phic
    ! 2020-08-31: bug correction

  !=======================================================================
    if (nsboost <= 1) return

    nrh = irmax-irmin+1
    do i=1, ngr ! loop on low-res rings
       i0 = ringphi(1,i) * nsboost - nsboost - irmin
       do k=-1,1,2 ! West and East side of disc
          nshrink = 0
          kk = (k+5)/2 ! 2 or 3
222       continue
          iphi = ringphi(kk, i)
          !if (ringphi(2,i) <= ringphi(3,i) .and. iphi >= 0) then
          if (iphi >= 0) then ! corrected 2020-08-31
             call find_pixel_bounds(nside, nsboost, ringphi(1,i), iphi, phiw, phie, np=np)
             do j=1, 2*nsboost+1
                if (i0+j >= 1 .and. i0+j <= nrh) then
                   phic = (phie(j)+phiw(j))*0.5_dp ! pixel center
                   dph  = (phie(j)-phiw(j))*0.5_dp + dphi(i0+j) ! pixel size + circle radius
                   dd = abs(phi0 - phic)    ! distance from disc center to pixel border sample
                   dd = min(dd, twopi - dd) ! in [0,Pi]
                   if (dd <= dph) goto 1000 ! pixel touched by disc, move to next one
                endif
             enddo
             !ringphi(kk, i)= iphi - k ! pixel not in disc, move edge pixel inwards
             ringphi(kk, i) = mod(iphi + np - k, np) ! pixel not in disc, move edge pixel inwards in [0,np-1]
             nshrink = nshrink + 1
             if (nshrink <= 2) then
                goto 222 ! try next pixel inward
             else
                ringphi(kk,i) = -999
             endif
1000         continue
          endif
       enddo ! loop on k East-West sides
    enddo ! loop on low-res rings

    ! remove empty rings
    ngr_out = 0
    do i=1,ngr
!        diff = ringphi(3,i) - ringphi(2,i)
!        if (ringphi(2,i) >=0 .and. ringphi(3,i) >=0 .and. diff /= -2 .and. diff /= -1) then
       if (ringphi(2,i) >=0 .and. ringphi(3,i) >=0) then ! 2020-08-31
          ngr_out = ngr_out + 1
          ringphi(1:3, ngr_out) = ringphi(1:3, i)
       endif
    enddo
    ! set empty rings to -1
    do i=ngr_out+1, ngr
       ringphi(2:3, i) = -1
    enddo

    ngr = ngr_out
    return
  end subroutine check_edge_pixels
  !=======================================================================
  subroutine find_pixel_bounds (nside, nsboost, iring, iphi, phiw, phie, np)
    !=======================================================================
    integer(i4b),               intent(in)  :: nside, nsboost, iring, iphi
    real(dp),     dimension(1:2*nsboost+1), intent(out) :: phiw, phie
    integer(i4b), optional,                 intent(out) :: np
    
    real(dp),     dimension(1:2*nsboost+1) :: f, f1, phiw_t, phie_t
    real(dp) :: c0, quad, phie1, phie2, phiw1, phiw2, cv
    integer(i4b) :: npr, kshift, nq, ip, i
    logical(lgt) :: transition
  !=======================================================================

    call pixels_per_ring(nside, iring, npr, kshift)
    !f = ((/ (i,i=0,2*nsboost) /) - nsboost) / nsboost
    f = ((/ (i*1.d0,i=0,2*nsboost) /) - nsboost*1.d0) / nsboost

    if (present(np)) np = npr
    nq = npr/4 ! number of pixels on current ring in [0,Pi/2] (quadrant)
    transition = (iring == nside .or. iring == nside*3)

    if (nq == nside .or. transition) then ! equatorial region (and transition rings)

       f1 = (1.0_dp-abs(f))*0.5_dp    ! triangle of height 1/2
       f1 = halfpi * f1 / nq
       c0 = halfpi * (iphi + kshift*0.5_dp) / nq
       phiw = c0 - f1
       phie = c0 + f1
       if (transition) then ! store for future use
          phiw_t = phiw
          phie_t = phie
       endif
    endif

    if (nq < nside .or. transition) then ! polar regions and transition rings
       ip = mod(iphi,nq) ! in [0,nq-1]
       quad = iphi / nq ! quadrant in [0,3]
       if (iring <= nside*2) then
          f1 = halfpi / (nq + f) 
       else
          f1 = halfpi / (nq - f)! swap sign for South pole
       endif
       do i=1, 2*nsboost+1
          cv = f1(i)
          phiw1 = min(cv *     ip,    halfpi)
          phie1 = min(cv * (   ip+1), halfpi)
          phiw2 = min(cv * (nq-ip-1), halfpi)
          phie2 = min(cv * (nq-ip),   halfpi)
          phiw(i) = max(phiw1, halfpi - phie2) + (quad * halfpi)
          phie(i) = min(phie1, halfpi - phiw2) + (quad * halfpi)
       enddo
    endif

    if (transition) then 
       if (iring == nside) then ! transition in N hemisphere
          phiw(nsboost+2:2*nsboost+1) = phiw_t(nsboost+2:2*nsboost+1) 
          phie(nsboost+2:2*nsboost+1) = phie_t(nsboost+2:2*nsboost+1)
       else ! transition in S hemisphere
          phiw(1:nsboost+1) = phiw_t(1:nsboost+1)
          phie(1:nsboost+1) = phie_t(1:nsboost+1)
       endif
    endif

    return
  end subroutine find_pixel_bounds
  
!   !=======================================================================
!   subroutine correct_ring_phi(location, iring, iphi)
!     !=======================================================================
!     ! returns ring and phi indexes corrected from round-off errors
!     ! appearing at large Nside
!     ! if phi <0 :      move 1 ring North
!     ! if phi >4*iring: move 1 ring South
!     ! rings are counted from the closest pole starting at 1.
!     integer(i4b), intent(in)    :: location !+1:North, -1:South
!     integer(i4b), intent(inout) :: iring, iphi
!     integer(i4b) :: delta
!     !-----------------------------------------------------------------------
!     delta = 0
!     if (iphi < 0)        delta =  -1
!     if (iphi >= 4*iring) delta =  +1
!     if (delta /= 0) then
!        if (abs(location) /= 1) then 
!           stop 'wrong location'
!        endif
!        if (delta > 0) then ! too large iphi
!           iphi  = iphi  - 4*iring  ! use old iring
!           iring = iring + location ! move south
!        else ! too small iphi
!           iring = iring - location ! move north
!           iphi  = iphi  +  4*iring ! use new iring
!        endif
!     endif
!     return
!   end subroutine correct_ring_phi
  !=======================================================================
  ! CHEAP_ISQRT
  !       Returns exact Floor(sqrt(x)) where x is a (64 bit) integer.
  !             y^2 <= x < (y+1)^2         (1)
  !       The double precision floating point operation is not accurate enough
  !       when dealing with 64 bit integers, especially in the vicinity of 
  !       perfect squares. 
  !=======================================================================
#ifndef NO64BITS
  function cheap_isqrt_8(lin) result (lout)
    integer(i8b), intent(in) :: lin
    integer(i8b) :: lout, diff
    real(DP) :: dout, din
    lout = floor(sqrt(dble(lin)), kind=I8B) ! round-off error may offset result
    diff = lin - lout*lout ! test Eq (1)
    if (diff <0)      lout = lout - 1
    if (diff >2*lout) lout = lout + 1
    return
  end function cheap_isqrt_8
#endif
  function cheap_isqrt_4(lin) result (lout)
    integer(i4b), intent(in) :: lin
    integer(i4b) :: lout
    lout = floor(sqrt(dble(lin)), kind=I4B)
    return
  end function cheap_isqrt_4
  !=======================================================================

  ! perform 2 compilations of the same source file
#undef DOI8Bi
#include "pixel_routines.F90"
#ifndef NO64BITS
#define DOI8B
#include "pixel_routines.F90"
#endif
!
! pix2ang_ring
! pix2vec_ring
! ang2pix_ring
! vec2pix_ring
! pix2ang_nest
! pix2vec_nest
! ang2pix_nest
! vec2pix_nest
! nest2ring
! npix2nside

  !*************************************************************
  !
  !                    MAP Manipulations
  !
  !*************************************************************
  subroutine warning_oldbounds(code, cos_theta_cut, zbounds)
    character(len=*), intent(in)          :: code
    real(DP), intent(in)                  :: cos_theta_cut
    real(DP), intent(out), dimension(1:2) :: zbounds

    if (cos_theta_cut <= 0.0_dp) then ! no cut
       zbounds(1) =  -1.0_dp
       zbounds(2) =   1.0_dp
    else
       zbounds(1) =   cos_theta_cut
       zbounds(2) = - cos_theta_cut
    endif
    print*,' -------------------------------------'
    print*,'WARNING: obsolete interface to '//code
    print*,'    cos_theta_cut currently a DP scalar with value'
    write(*,9000) '    ',cos_theta_cut
    print*,'    shoud now be replaced with a 2-element vector with values:'
    write(*,9001) '    ',zbounds(1),zbounds(2)
    print*,'    See documentation for details.'
    print*,' -------------------------------------'
9000 format (a,g13.6)
9001 format (a,g13.6,g13.6)

    return
  end subroutine warning_oldbounds
  !==============================================================
  ! REMOVE_DIPOLE( nside, map, ordering, degree, multipoles, zbounds, fmissval, mask, weights)
  !
  ! removes monopole (and dipole) from a map
  !
  ! Nside:     I4,       IN   : Healpix resolution parameter
  ! map:       KMAP, array,INOUT: Heapix map (see Notes below)
  ! ordering:  I4,       IN:   Healpix scheme 1:RING, 2: NESTED
  ! degree:    I4,       IN:   multipole to remove, 1: monopole, 2: monopole and dipole
  ! multipoles:R8, array,OUT:  value of monopole and dipole
  ! zbounds   :R8,      IN:  2-el vector
  ! range of z in [-1,1] on which to estimate monopole and dipole
  ! fmissval:  KMAP, Option, IN: value used to flag bad pixel on input, default=-1.6375e30
  !                            Pixels with map = fmissval are not used for fit
  ! mask    :  KMAP, Option, IN: Pixels with |mask|<1.e-10 are not used for fit
  !                              others are kept as is
  !                               note : the mask in NOT applied to the map
  ! weights  : KMAP, Option, IN: pixel weight to be applied to the map
  !
  ! KMAP: either R4 or R8
  !
  ! note : if degree= 1, or 2, the map is modified on output
  !     * the monopole (and dipole) is/are removed
  !     * pixels within the symmetric cut parameterized
  !       by cos_theta_cut are set to fmissval (or its default value)
  !  if degree = 0, nothing is done
  !  all other values of degree are invalid
  !
  ! v1.0, EH, Caltech, Jan-2002, based on homonyme IDL routine
  ! 2008-02-03: addition of weights
  !
  !==============================================================
  subroutine remove_dipole_real( nside, map, ordering, degree, multipoles, zbounds, &
       & fmissval, mask, weights)
    !============================================================
    use num_rec, only : dsvdcmp, dsvbksb
    ! parameters
    character(len=*),      parameter :: code = "REMOVE_DIPOLE_REAL"
    integer(kind=I4B),     parameter :: KMAP = SP
    real   (kind=KMAP),    parameter :: fbad_value = -1.6375e30_KMAP
    !
    include 'remove_dipole_inc.f90'
  end subroutine remove_dipole_real
  !==============================================================
  subroutine remove_dipole_double( nside, map, ordering, degree, multipoles, zbounds, &
       & fmissval, mask, weights)
    !============================================================
    use num_rec, only : dsvdcmp, dsvbksb
    ! parameters
    character(len=*),      parameter :: code = "REMOVE_DIPOLE_DOUBLE"
    integer(kind=I4B),     parameter :: KMAP = DP
    real   (kind=KMAP),    parameter :: fbad_value = -1.6375e30_KMAP
    !
    include 'remove_dipole_inc.f90'
  end subroutine remove_dipole_double

  !==============================================================
  subroutine remove_dipole_real_old( nside, map, ordering, degree, multipoles, cos_theta_cut, fmissval, mask)
    !============================================================
    ! parameters
    character(len=*),      parameter :: code = "REMOVE_DIPOLE_REAL"
    integer(kind=I4B),     parameter :: KMAP = SP
    ! dummy
    integer(kind=i4b),                  intent(in)    :: nside
    integer(kind=i4b),                  intent(in)    :: ordering, degree
    real   (kind=DP),   dimension(0:),  intent(out)   :: multipoles
    real   (kind=KMAP), dimension(0:),  intent(inout) :: map
    real   (kind=DP),                   intent(in)    :: cos_theta_cut
    real   (kind=KMAP),                 intent(in), optional :: fmissval
    real   (kind=KMAP), dimension(0:),  intent(in), optional :: mask
    ! local
    real(DP),                          dimension(1:2) :: zbounds

    call warning_oldbounds(code, cos_theta_cut, zbounds)
    call remove_dipole(nside, map, ordering, degree, multipoles, zbounds, fmissval, mask)
  end subroutine remove_dipole_real_old
  !==============================================================
  subroutine remove_dipole_double_old( nside, map, ordering, degree, multipoles, cos_theta_cut, fmissval, mask)
    !============================================================
     ! parameters
    character(len=*),      parameter :: code = "REMOVE_DIPOLE_DOUBLE"
    integer(kind=I4B),     parameter :: KMAP = DP
    ! dummy
    integer(kind=i4b),                  intent(in)    :: nside
    integer(kind=i4b),                  intent(in)    :: ordering, degree
    real   (kind=DP),   dimension(0:),  intent(out)   :: multipoles
    real   (kind=KMAP), dimension(0:),  intent(inout) :: map
    real   (kind=DP),                   intent(in)    :: cos_theta_cut
    real   (kind=KMAP),                 intent(in), optional :: fmissval
    real   (kind=KMAP), dimension(0:),  intent(in), optional :: mask
    ! local
    real(DP),                          dimension(1:2) :: zbounds

    call warning_oldbounds(code, cos_theta_cut, zbounds)
    call remove_dipole(nside, map, ordering, degree, multipoles, zbounds, fmissval, mask) 
  end subroutine remove_dipole_double_old

  !==============================================================
  subroutine remove_dipole_real_v12( nside, map, nmaps, ordering, degree, multipoles, cos_theta_cut, fmissval, mask)
    !============================================================
    ! parameters
    character(len=*),      parameter :: code = "REMOVE_DIPOLE_REAL"
    integer(kind=I4B),     parameter :: KMAP = SP
    ! dummy
    integer(kind=i4b),                  intent(in)    :: nside
    integer(kind=i4b),                  intent(in)    :: ordering, degree, nmaps
    real   (kind=DP),   dimension(0:),  intent(out)   :: multipoles
    real   (kind=KMAP), dimension(0:),  intent(inout) :: map
    real   (kind=DP),                   intent(in)    :: cos_theta_cut
    real   (kind=KMAP),                 intent(in), optional :: fmissval
    real   (kind=KMAP), dimension(0:),  intent(in), optional :: mask
    ! local
    real(DP),                          dimension(1:2) :: zbounds

    print*,'=========================================================='
    print*,'WARNING: Interface to remove_dipole has changed'
    print*,' from'
    print*,'call remove_dipole(nside, map, NMAPS, ordering, degree, multipoles, COS_THETA_CUT, fmissval, mask)'
    print*,' to'
    print*,'call remove_dipole(nside, map,        ordering, degree, multipoles, ZBOUNDS,       fmissval, mask)'
    print*,'=========================================================='
    
    call warning_oldbounds(code, cos_theta_cut, zbounds)
    call remove_dipole(nside, map, ordering, degree, multipoles, zbounds, fmissval, mask)
  end subroutine remove_dipole_real_v12
  !==============================================================
  subroutine remove_dipole_double_v12( nside, map, nmaps, ordering, degree, multipoles, cos_theta_cut, fmissval, mask)
    !============================================================
     ! parameters
    character(len=*),      parameter :: code = "REMOVE_DIPOLE_DOUBLE"
    integer(kind=I4B),     parameter :: KMAP = DP
    ! dummy
    integer(kind=i4b),                  intent(in)    :: nside
    integer(kind=i4b),                  intent(in)    :: ordering, degree, nmaps
    real   (kind=DP),   dimension(0:),  intent(out)   :: multipoles
    real   (kind=KMAP), dimension(0:),  intent(inout) :: map
    real   (kind=DP),                   intent(in)    :: cos_theta_cut
    real   (kind=KMAP),                 intent(in), optional :: fmissval
    real   (kind=KMAP), dimension(0:),  intent(in), optional :: mask
    ! local
    real(DP),                          dimension(1:2) :: zbounds

    print*,'=========================================================='
    print*,'WARNING: Interface to remove_dipole has changed'
    print*,' from'
    print*,'call remove_dipole(nside, map, NMAPS, ordering, degree, multipoles, COS_THETA_CUT, fmissval, mask)'
    print*,' to'
    print*,'call remove_dipole(nside, map,        ordering, degree, multipoles, ZBOUNDS,       fmissval, mask)'
    print*,'=========================================================='

    call warning_oldbounds(code, cos_theta_cut, zbounds)
    call remove_dipole(nside, map, ordering, degree, multipoles, zbounds, fmissval, mask) 
  end subroutine remove_dipole_double_v12

  !=======================================================================
  ! ADD_DIPOLE( nside, map, ordering, degree, multipoles, fmissval)
  !
  ! removes monopole (and dipole) from a map
  !
  ! Nside:     I4,       IN   : Healpix resolution parameter
  ! map:       KMAP, array,INOUT: Heapix map (see Notes below)
  ! ordering:  I4,       IN:   Healpix scheme 1:RING, 2: NESTED
  ! degree:    I4,       IN:   multipole to remove, 1: monopole, 2: monopole and dipole
  ! multipoles:R8, array,IN:  value of monopole and dipole
  ! fmissval:  KMAP, Option, IN: value used to flag bad pixel on input, default=-1.6375e30
  !                            Pixels with map = fmissval are left unchanged
  ! 2009-03-25: accepts Nside > 8192
  !=======================================================================
  subroutine add_dipole_real(nside, map, ordering, degree, multipoles, fmissval)
    !=======================================================================
    ! single precision
    !=======================================================================
    character(len=*),      parameter :: code = "ADD_DIPOLE_REAL"
    integer(kind=I4B),     parameter :: KMAP = SP
    real   (kind=KMAP),    parameter :: fbad_value = -1.6375e30_KMAP
    ! dummy
    integer(kind=i4b),                  intent(in)    :: nside
    integer(kind=i4b),                  intent(in)    :: ordering, degree
    real   (kind=DP),   dimension(0:),  intent(in)    :: multipoles
    real   (kind=KMAP), dimension(0:),  intent(inout) :: map
    real   (kind=KMAP),                 intent(in), optional :: fmissval

    real   (kind=KMAP)                :: fmiss_effct
    integer(kind=i8b)                 :: ipix, npix
    logical(lgt)                      :: dodipole !, do_mask
    real(kind=dp), dimension(1:3)     :: vec
    !=======================================================================
    
    npix = nside2npix(nside)
    fmiss_effct = fbad_value
    if (present(fmissval)) fmiss_effct = fmissval

    if (degree == 0) then
       print*," No monopole nor dipole to add"
       return
    elseif (degree == 1) then
       dodipole = .false.
    else if (degree == 2) then
       dodipole = .true.
    else
       print*,code//"> degree can only be "
       print*,"      1: monopole (l=0) addition or "
       print*,"      2: monopole and dipole (l=0,1) addition"
       print*,code//"> ABORT ! "
       call fatal_error
    endif

    do ipix = 0, npix-1
       if ( abs(map(ipix) - fmiss_effct) <= abs(1.e-5*fmiss_effct) ) goto 20
!        if (do_mask) then
!           if (abs(mask(ipix)) <= 1.e-10) goto 20
!        endif
       map(ipix) = map(ipix) + multipoles(0)
       if (dodipole) then
          ! computes dipole basis functions
          ! pixel -> vector
          if (ordering == 1) call pix2vec_ring( nside, ipix, vec)
          if (ordering == 2) call pix2vec_nest( nside, ipix, vec)
          map(ipix) = map(ipix) + sum(multipoles(1:3) * vec(1:3))
       endif

20     continue
    enddo

    return
  end subroutine add_dipole_real
  !=======================================================================
  subroutine add_dipole_double(nside, map, ordering, degree, multipoles, fmissval)
    !=======================================================================
    ! single precision
    !=======================================================================
    character(len=*),      parameter :: code = "ADD_DIPOLE_DOUBLE"
    integer(kind=I4B),     parameter :: KMAP = DP
    real   (kind=KMAP),    parameter :: fbad_value = -1.6375e30_KMAP
    ! dummy
    integer(kind=i4b),                  intent(in)    :: nside
    integer(kind=i4b),                  intent(in)    :: ordering, degree
    real   (kind=DP),   dimension(0:),  intent(in)    :: multipoles
    real   (kind=KMAP), dimension(0:),  intent(inout) :: map
    real   (kind=KMAP),                 intent(in), optional :: fmissval

    real   (kind=KMAP)                :: fmiss_effct
    integer(kind=i8b)                 :: ipix, npix
    logical(lgt)                      :: dodipole !, do_mask
    real(kind=dp), dimension(1:3)     :: vec
    !=======================================================================
    
    npix = nside2npix(nside)
    fmiss_effct = fbad_value
    if (present(fmissval)) fmiss_effct = fmissval

    if (degree == 0) then
       print*," No monopole nor dipole to add"
       return
    elseif (degree == 1) then
       dodipole = .false.
    else if (degree == 2) then
       dodipole = .true.
    else
       print*,code//"> degree can only be "
       print*,"      1: monopole (l=0) addition or "
       print*,"      2: monopole and dipole (l=0,1) addition"
       print*,code//"> ABORT ! "
       call fatal_error
    endif

    do ipix = 0, npix-1
       if ( abs(map(ipix) - fmiss_effct) <= abs(1.e-5*fmiss_effct) ) goto 20
!        if (do_mask) then
!           if (abs(mask(ipix)) <= 1.e-10) goto 20
!        endif
       map(ipix) = map(ipix) + multipoles(0)
       if (dodipole) then
          ! computes dipole basis functions
          ! pixel -> vector
          if (ordering == 1) call pix2vec_ring( nside, ipix, vec)
          if (ordering == 2) call pix2vec_nest( nside, ipix, vec)
          map(ipix) = map(ipix) + sum(multipoles(1:3) * vec(1:3))
       endif

20     continue
    enddo

    return
  end subroutine add_dipole_double

  !====================================================
  ! apply_mask (map, order, mask=mask, zbounds=zbounds)
  !   applies a mask and/or a zbounds to a TQU maps,
  !   map and mask should have the same ordering
  !  map :  SP/DP, INOUT,   (0:,1:)
  !  order: I4B,   IN,  1=RING, 2=NESTED
  !  mask:  SP/DP  IN, (0:,1:) optional
  !  zbounds: DP   IN, (1:2)   optional
  !=====================================================================
  subroutine apply_mask_real(map, order, mask, zbounds)
    !=====================================================================
    integer(I4B), parameter :: KMAP = SP
    include 'apply_mask_inc.f90'
  end subroutine apply_mask_real
  !=====================================================================
  subroutine apply_mask_double(map, order, mask, zbounds)
    !=====================================================================
    integer(I4B), parameter :: KMAP = DP
    include 'apply_mask_inc.f90'
  end subroutine apply_mask_double
  !====================================================
  ! medfiltmap
  !   compute the median filtered map of a given Healpix map
  !   in_map: SP/DP  input Healpix full sky map
  !   radius: DP     radius in Radians
  !   med_map: SP/DP output Healpix full sky map
  !   nest:    I4B, OPT   either 0 (ring scheme) or 1 (nested scheme)
  !   fmissval:SP/DP, OPT sentinel value given to missing pixels
  !   fill_holes: LGT, OPT 
  ! 2012-07-17: Parallel OpenMP implementation
  !=================================================================
  subroutine medfiltmap_S( in_map, radius, med_map, nest, fmissval, fill_holes)
  !=================================================================
    use statistics, only: median
    integer(I4B), parameter :: KMAP = SP
    !
    real(KMAP), dimension(0:),intent(in), target       :: in_map
    real(DP),                 intent(in)               :: radius
    real(KMAP), dimension(0:),intent(out)              :: med_map
    integer(I4B),             intent(in),  optional    :: nest
    real(KMAP),               intent(in),  optional    :: fmissval
    logical(LGT),             intent(in),  optional    :: fill_holes
    !
    integer(I8B) :: npix, np
    integer(I4B) :: nside, p, nlist, status
    logical(LGT) :: do_nest, do_fill
    real(DP)     :: fraction
    real(KMAP)   :: fmissval_in
    integer(I4B), dimension(:),  allocatable :: listpix
    real(DP),     dimension(1:3)             :: vector
    character(len=*), parameter :: code = 'medfiltmap'
    !-----------------------------------------------
    npix = long_size(in_map)
    nside = npix2nside(npix)
    call assert(nside > 0, code//": invalid map size")

    fraction = 0.5_DP * (1.0_dp - cos(radius))
    np = npix * fraction * 1.2 + 50
    call assert(np < MAX_I4B, code//": too many pixels to compute median")

    do_nest = .false.
    if (present(nest)) then
       call assert(nest>=0 .and. nest <=1,code//': invalid NEST flag')
       do_nest = (nest == 1)
    endif

    do_fill = .false.
    if (present(fill_holes)) do_fill = fill_holes

    fmissval_in = hpx_Sbadval
    if (present(fmissval)) fmissval_in = fmissval

    ! make sure common arrays are initiated
    call mk_pix2xy()
!!!    print*,'************* Parallel Median **************'
!$OMP parallel default(none) &
!$OMP   shared(in_map, med_map, pix2x, pix2y, &
!$OMP          nside, npix, radius, np, do_nest, do_fill, nest, fmissval_in) &
!$OMP  private(listpix, vector, p, nlist, status)
    allocate(listpix(0:np-1),stat=status)
    call assert_alloc(status,code,'listpix')
!$OMP do schedule(dynamic, 64)
    do p = 0, npix-1
       ! find pixel location
       if (do_nest) then
          call pix2vec_nest( nside, p, vector)
       else
          call pix2vec_ring( nside, p, vector)
       endif

       ! find disc centered on pixel
       call query_disc(nside, vector, radius, listpix, nlist, nest=nest)

       if (do_fill .or. abs(in_map(p)-fmissval_in) > abs(fmissval_in*1.e-7)) then
          med_map(p) = median(in_map(listpix(0:nlist-1)), badval = fmissval_in, even=(nlist<100))
       else
          med_map(p) = in_map(p)
       endif
    enddo
!$OMP end do
    deallocate(listpix)
!$OMP end parallel
    return
  end subroutine medfiltmap_S
  !=================================================================
  subroutine medfiltmap_D( in_map, radius, med_map, nest, fmissval, fill_holes)
  !=================================================================
    use statistics, only: median
    integer(I4B), parameter :: KMAP = DP
    !
    real(KMAP), dimension(0:),intent(in), target       :: in_map
    real(DP),                 intent(in)               :: radius
    real(KMAP), dimension(0:),intent(out)              :: med_map
    integer(I4B),             intent(in),  optional    :: nest
    real(KMAP),               intent(in),  optional    :: fmissval
    logical(LGT),             intent(in),  optional    :: fill_holes
    !
    integer(I8B) :: npix, np
    integer(I4B) :: nside, p, nlist, status
    logical(LGT) :: do_nest, do_fill
    real(DP)     :: fraction
    real(KMAP)   :: fmissval_in
    integer(I4B), dimension(:),  allocatable :: listpix
    real(DP),     dimension(1:3)             :: vector
    character(len=*), parameter :: code = 'medfiltmap'
    !-----------------------------------------------
    npix = long_size(in_map)
    nside = npix2nside(npix)
    call assert(nside > 0, code//": invalid map size")

    fraction = 0.5_DP * (1.0_dp - cos(radius))
    np = npix * fraction * 1.2 + 50
    call assert(np < MAX_I4B, code//": too many pixels to compute median")

    do_nest = .false.
    if (present(nest)) then
       call assert(nest>=0 .and. nest <=1,code//': invalid NEST flag')
       do_nest = (nest == 1)
    endif

    do_fill = .false.
    if (present(fill_holes)) do_fill = fill_holes

    fmissval_in = hpx_Dbadval
    if (present(fmissval)) fmissval_in = fmissval

    ! make sure common arrays are initiated
    call mk_pix2xy()
!!!    print*,'************* Parallel Median **************'
!$OMP parallel default(none) &
!$OMP   shared(in_map, med_map, pix2x, pix2y, &
!$OMP          nside, npix, radius, np, do_nest, do_fill, nest, fmissval_in) &
!$OMP  private(listpix, vector, p, nlist, status)
    allocate(listpix(0:np-1),stat=status)
    call assert_alloc(status,code,'listpix')
!$OMP do schedule(dynamic, 64)
    do p = 0, npix-1
       ! find pixel location
       if (do_nest) then
          call pix2vec_nest( nside, p, vector)
       else
          call pix2vec_ring( nside, p, vector)
       endif

       ! find disc centered on pixel
       call query_disc(nside, vector, radius, listpix, nlist, nest=nest)

       if (do_fill .or. abs(in_map(p)-fmissval_in) > abs(fmissval_in*1.e-7)) then
          med_map(p) = median(in_map(listpix(0:nlist-1)), badval = fmissval_in, even=(nlist<100))
       else
          med_map(p) = in_map(p)
       endif
    enddo
!$OMP end do
    deallocate(listpix)
!$OMP end parallel
    return
  end subroutine medfiltmap_D
  !**************************************************************
  ! following 2x3 routines: out of place conversions for 1D maps
  ! if size(map) = Npix, peak memory = 2*Npix
  ! These routines are parallelized for shared memory architecture 
  !   (using OpenMP directives)
  !**************************************************************
  !=======================================================================
  !     makes the conversion NEST to RING
  !=======================================================================
  subroutine convert_nest2ring_int_1d(nside, map)
    !=======================================================================
    character(len=*), parameter :: code = 'convert_nest2ring_int_1d'
    integer(kind=I4B),   parameter :: KMAP = I4B
    integer(kind=KMAP),  dimension(0:), intent(inout) ::  map
    integer(kind=KMAP),  dimension(:),  allocatable :: map_tmp
    include 'convert_nest2ring_1d_inc.f90'
  end subroutine convert_nest2ring_int_1d
#ifndef NO64BITS
  !=======================================================================
  subroutine convert_nest2ring_int8_1d(nside, map)
    !=======================================================================
    character(len=*), parameter :: code = 'convert_nest2ring_int8_1d'
    integer(kind=I4B),   parameter :: KMAP = I8B
    integer(kind=KMAP),  dimension(0:), intent(inout) ::  map
    integer(kind=KMAP),  dimension(:),  allocatable :: map_tmp
    include 'convert_nest2ring_1d_inc.f90'
  end subroutine convert_nest2ring_int8_1d
#endif
  !=======================================================================
  subroutine convert_nest2ring_real_1d(nside, map)
    !=======================================================================
    character(len=*), parameter :: code = 'convert_nest2ring_real_1d'
    integer(kind=I4B),   parameter :: KMAP = SP
    real   (kind=KMAP),  dimension(0:), intent(inout) ::  map
    real   (kind=KMAP),  dimension(:),  allocatable :: map_tmp
    include 'convert_nest2ring_1d_inc.f90'
  end subroutine convert_nest2ring_real_1d
  !=======================================================================
  subroutine convert_nest2ring_double_1d(nside, map)
    !=======================================================================
    character(len=*), parameter :: code = 'convert_nest2ring_double_1d'
    integer(kind=I4B),   parameter :: KMAP = DP
    real   (kind=KMAP),  dimension(0:), intent(inout) ::  map
    real   (kind=KMAP),  dimension(:),  allocatable :: map_tmp
    include 'convert_nest2ring_1d_inc.f90'
  end subroutine convert_nest2ring_double_1d
  !=======================================================================
  !     makes the conversion RING to NEST
  !=======================================================================
  subroutine convert_ring2nest_int_1d(nside, map)
    !=======================================================================
    character(len=*), parameter :: code = 'convert_ring2nest_int_1d'
    integer(kind=I4B),   parameter :: KMAP = I4B
    integer(kind=KMAP),  dimension(0:), intent(inout) ::  map
    integer(kind=KMAP),  dimension(:), allocatable :: map_tmp
    include 'convert_ring2nest_1d_inc.f90'
  end subroutine convert_ring2nest_int_1d
#ifndef NO64BITS
  !=======================================================================
  subroutine convert_ring2nest_int8_1d(nside, map)
    !=======================================================================
    character(len=*), parameter :: code = 'convert_ring2nest_int8_1d'
    integer(kind=I4B),   parameter :: KMAP = I8B
    integer(kind=KMAP),  dimension(0:), intent(inout) ::  map
    integer(kind=KMAP),  dimension(:), allocatable :: map_tmp
    include 'convert_ring2nest_1d_inc.f90'
  end subroutine convert_ring2nest_int8_1d
#endif
  !=======================================================================
  subroutine convert_ring2nest_real_1d(nside, map)
    !=======================================================================
    character(len=*), parameter :: code = 'convert_ring2nest_real_1d'
    integer(kind=I4B),   parameter :: KMAP = SP
    real   (kind=KMAP),  dimension(0:), intent(inout) ::  map
    real   (kind=KMAP),  dimension(:), allocatable :: map_tmp
    include 'convert_ring2nest_1d_inc.f90'
  end subroutine convert_ring2nest_real_1d
  !=======================================================================
  subroutine convert_ring2nest_double_1d(nside, map)
    !=======================================================================
    character(len=*), parameter :: code = 'convert_ring2nest_double_1d'
    integer(kind=I4B),   parameter :: KMAP = DP
    real   (kind=KMAP),  dimension(0:), intent(inout) ::  map
    real   (kind=KMAP),  dimension(:), allocatable :: map_tmp
    include 'convert_ring2nest_1d_inc.f90'
  end subroutine convert_ring2nest_double_1d

  !**************************************************************
  ! following 6 routines: out of place conversions for N-Dim maps
  ! if size(map) = Npix*Nd, peak memory = (Nd+2)*Npix
  ! 2004-08-25, EH
  !**************************************************************
  !=======================================================================
  subroutine convert_nest2ring_int_nd(nside, map)
    !=======================================================================
    !   NEST to RING conversion: 1D, integer
    !=======================================================================
    character(len=*),  parameter :: code = "convert_nest2ring_int_nd"
    integer(kind=I4B), parameter :: KMAP = I4B
    integer(kind=KMAP), dimension(0:,1:), intent(inout), target ::  map
    integer(kind=KMAP), dimension(:),               allocatable :: map_tmp
    integer(kind=KMAP), dimension(:),                   pointer :: map1
    include 'convert_nest2ring_nd_inc.f90'
  end subroutine convert_nest2ring_int_nd
#ifndef NO64BITS
  !=======================================================================
  subroutine convert_nest2ring_int8_nd(nside, map)
    !=======================================================================
    !   NEST to RING conversion: 1D, integer*8
    !=======================================================================
    character(len=*),  parameter :: code = "convert_nest2ring_int_nd"
    integer(kind=I4B), parameter :: KMAP = I8B
    integer(kind=KMAP), dimension(0:,1:), intent(inout), target ::  map
    integer(kind=KMAP), dimension(:),               allocatable :: map_tmp
    integer(kind=KMAP), dimension(:),                   pointer :: map1
    include 'convert_nest2ring_nd_inc.f90'
  end subroutine convert_nest2ring_int8_nd
#endif
  !=======================================================================
  subroutine convert_nest2ring_real_nd(nside, map)
    !=======================================================================
    !   NEST to RING conversion: 1D, real
    !=======================================================================
    character(len=*),  parameter :: code = "convert_nest2ring_real_nd"
    integer(kind=I4B), parameter :: KMAP = SP
    real(kind=KMAP), dimension(0:,1:), intent(inout), target ::  map
    real(kind=KMAP), dimension(:),               allocatable :: map_tmp
    real(kind=KMAP), dimension(:),                   pointer :: map1
    include 'convert_nest2ring_nd_inc.f90'
  end subroutine convert_nest2ring_real_nd
  !=======================================================================
  subroutine convert_nest2ring_double_nd(nside, map)
    !=======================================================================
    !   NEST to RING conversion: 1D, double
    !=======================================================================
    character(len=*),  parameter :: code = "convert_nest2ring_double_nd"
    integer(kind=I4B), parameter :: KMAP = DP
    real(kind=KMAP), dimension(0:,1:), intent(inout), target ::  map
    real(kind=KMAP), dimension(:),               allocatable :: map_tmp
    real(kind=KMAP), dimension(:),                   pointer :: map1
    include 'convert_nest2ring_nd_inc.f90'
  end subroutine convert_nest2ring_double_nd
  !=======================================================================
  subroutine convert_ring2nest_int_nd(nside, map)
    !=======================================================================
    !   RING to NEST conversion: 1D, integer
    !=======================================================================
    character(len=*),  parameter :: code = "convert_ring2nest_int_nd"
    integer(kind=I4B), parameter :: KMAP = I4B
    integer(kind=KMAP), dimension(0:,1:), intent(inout), target ::  map
    integer(kind=KMAP), dimension(:),               allocatable :: map_tmp
    integer(kind=KMAP), dimension(:),                   pointer :: map1
    include 'convert_ring2nest_nd_inc.f90'
  end subroutine convert_ring2nest_int_nd
#ifndef NO64BITS
  !=======================================================================
  subroutine convert_ring2nest_int8_nd(nside, map)
    !=======================================================================
    !   RING to NEST conversion: 1D, integer
    !=======================================================================
    character(len=*),  parameter :: code = "convert_ring2nest_int8_nd"
    integer(kind=I4B), parameter :: KMAP = I8B
    integer(kind=KMAP), dimension(0:,1:), intent(inout), target ::  map
    integer(kind=KMAP), dimension(:),               allocatable :: map_tmp
    integer(kind=KMAP), dimension(:),                   pointer :: map1
    include 'convert_ring2nest_nd_inc.f90'
  end subroutine convert_ring2nest_int8_nd
#endif
  !=======================================================================
  subroutine convert_ring2nest_real_nd(nside, map)
    !=======================================================================
    !   RING to NEST conversion: 1D, real
    !=======================================================================
    character(len=*),  parameter :: code = "convert_ring2nest_real_nd"
    integer(kind=I4B), parameter :: KMAP = SP
    real(kind=KMAP), dimension(0:,1:), intent(inout), target ::  map
    real(kind=KMAP), dimension(:),               allocatable :: map_tmp
    real(kind=KMAP), dimension(:),                   pointer :: map1
    include 'convert_ring2nest_nd_inc.f90'
  end subroutine convert_ring2nest_real_nd
  !=======================================================================
  subroutine convert_ring2nest_double_nd(nside, map)
    !=======================================================================
    !   RING to NEST conversion: 1D, double
    !=======================================================================
    character(len=*),  parameter :: code = "convert_ring2nest_double_nd"
    integer(kind=I4B), parameter :: KMAP = DP
    real(kind=KMAP), dimension(0:,1:), intent(inout), target ::  map
    real(kind=KMAP), dimension(:),               allocatable :: map_tmp
    real(kind=KMAP), dimension(:),                   pointer :: map1
    include 'convert_ring2nest_nd_inc.f90'
  end subroutine convert_ring2nest_double_nd
  !====================================================================
  ! The following 6 routines convert in place integer, real, and double
  ! arrays between the NESTED and RING schemes.
  !
  ! in place: without allocating a temporary map. This routine is more general,
  ! but slower than convert_nest2ring.
  !
  !     This is a wrapper for the toolbox functions "ring2nest" and
  !     "nest2ring". Their names are supplied in the "subcall"
  !     argument.
  !
  ! Author: Benjamin D. Wandelt October 1997
  ! Added to pix_tools for version 1.00 in March 1999.
  ! 2004-08-25: EH, added N-Dim facility and double precision IO
  !====================================================================
  !====================================================================
  subroutine convert_inplace_int_1d(subcall,map)
    !==================================================================
    ! 1D, integer implementation
    !==================================================================
    character(len=*),  parameter :: code = "convert_inplace_int_1d"
    integer(kind=I4B), parameter :: KMAP = I4B
    integer(kind=KMAP), dimension(0:)    :: map
    integer(kind=KMAP)                   :: pixbuf1,pixbuf2
    include 'convert_inplace_1d_inc.f90'
  end subroutine convert_inplace_int_1d
  !====================================================================
  subroutine convert_inplace_real_1d(subcall,map)
    !====================================================================
    ! 1D, real implementation
    !==================================================================
    character(len=*),  parameter :: code = "convert_inplace_real_1d"
    integer(kind=I4B), parameter :: KMAP = SP
    real   (kind=KMAP), dimension(0:)    :: map
    real   (kind=KMAP)                   :: pixbuf1,pixbuf2
    include 'convert_inplace_1d_inc.f90'
  end subroutine convert_inplace_real_1d
  !====================================================================
  subroutine convert_inplace_double_1d(subcall,map)
    !====================================================================
    ! 1D, double precision implementation
    !==================================================================
    character(len=*),  parameter :: code = "convert_inplace_double_1d"
    integer(kind=I4B), parameter :: KMAP = DP
    real   (kind=KMAP), dimension(0:)    :: map
    real   (kind=KMAP)                   :: pixbuf1,pixbuf2
    include 'convert_inplace_1d_inc.f90'
  end subroutine convert_inplace_double_1d
  !====================================================================
  subroutine convert_inplace_int_nd(subcall,map)
    !==================================================================
    ! ND, integer implementation
    !==================================================================
    character(len=*),  parameter :: code = "convert_inplace_int_nd"
    integer(kind=i4b), parameter :: ND_MAX = 10
    integer(kind=I4B), parameter :: KMAP = I4B
    integer(kind=KMAP), dimension(0:,1:)    :: map
    integer(kind=KMAP), dimension(1:ND_MAX) :: pixbuf1,pixbuf2
    include 'convert_inplace_nd_inc.f90'
  end subroutine convert_inplace_int_nd
  !====================================================================
  subroutine convert_inplace_real_nd(subcall,map)
    !====================================================================
    ! ND, real implementation
    !==================================================================
    character(len=*),  parameter :: code = "convert_inplace_real_nd"
    integer(kind=i4b), parameter :: ND_MAX = 10
    integer(kind=I4B), parameter :: KMAP = SP
    real   (kind=KMAP), dimension(0:,1:)    :: map
    real   (kind=KMAP), dimension(1:ND_MAX) :: pixbuf1,pixbuf2
    include 'convert_inplace_nd_inc.f90'
  end subroutine convert_inplace_real_nd
  !====================================================================
  subroutine convert_inplace_double_nd(subcall,map)
    !====================================================================
    ! ND, double precision implementation
    !==================================================================
    character(len=*),  parameter :: code = "convert_inplace_double_nd"
    integer(kind=i4b), parameter :: ND_MAX = 10
    integer(kind=I4B), parameter :: KMAP = DP
    real   (kind=KMAP), dimension(0:,1:)    :: map
    real   (kind=KMAP), dimension(1:ND_MAX) :: pixbuf1,pixbuf2
    include 'convert_inplace_nd_inc.f90'
  end subroutine convert_inplace_double_nd
  !=======================================================================
  subroutine mk_pix2xy()
    !=======================================================================
    !     constructs the array giving x and y in the face from pixel number
    !     for the nested (quad-cube like) ordering of pixels
    !
    !     the bits corresponding to x and y are interleaved in the pixel number
    !     one breaks up the pixel number by even and odd bits
    !=======================================================================
    INTEGER(KIND=I4B) ::  kpix, jpix, ix, iy, ip, id

    !cc cf block data      data      pix2x(1023) /0/
    !-----------------------------------------------------------------------
    !      print *, 'initiate pix2xy'
    do kpix=0,1023          ! pixel number
       jpix = kpix
       IX = 0
       IY = 0
       IP = 1               ! bit position (in x and y)
!        do while (jpix/=0) ! go through all the bits
       do
          if (jpix == 0) exit ! go through all the bits
          ID = MODULO(jpix,2)  ! bit value (in kpix), goes in ix
          jpix = jpix/2
          IX = ID*IP+IX

          ID = MODULO(jpix,2)  ! bit value (in kpix), goes in iy
          jpix = jpix/2
          IY = ID*IP+IY

          IP = 2*IP         ! next bit (in x and y)
       enddo
       pix2x(kpix) = IX     ! in 0,31
       pix2y(kpix) = IY     ! in 0,31
    enddo

    return
  end subroutine mk_pix2xy
!   !=======================================================================
!   subroutine mk_xy2pix()
!     !=======================================================================
!     !     sets the array giving the number of the pixel lying in (x,y)
!     !     x and y are in {1,128}
!     !     the pixel number is in {0,128**2-1}
!     !
!     !     if  i-1 = sum_p=0  b_p * 2^p
!     !     then ix = sum_p=0  b_p * 4^p
!     !          iy = 2*ix
!     !     ix + iy in {0, 128**2 -1}
!     !=======================================================================
!     INTEGER(KIND=I4B):: k,ip,i,j,id
!     !=======================================================================

!     do i = 1,128           !for converting x,y into
!        j  = i-1            !pixel numbers
!        k  = 0
!        ip = 1

!        do
!           if (j==0) then
!              x2pix(i) = k
!              y2pix(i) = 2*k
!              exit
!           else
!              id = MODULO(J,2)
!              j  = j/2
!              k  = ip*id+k
!              ip = ip*4
!           endif
!        enddo

!     enddo

!     return
!   end subroutine mk_xy2pix
  !=======================================================================
  subroutine mk_xy2pix1()
    !=======================================================================
    !     sets the array giving the number of the pixel lying in (x,y)
    !     x and y are in {1,128}
    !     the pixel number is in {0,128**2-1}
    !
    !     if  i-1 = sum_p=0  b_p * 2^p
    !     then ix = sum_p=0  b_p * 4^p
    !          iy = 2*ix
    !     ix + iy in {0, 128**2 -1}
    !=======================================================================
    INTEGER(KIND=I4B):: k,ip,i,j,id
    !=======================================================================

    do i = 0,127           !for converting x,y into
       j  = i           !pixel numbers
       k  = 0
       ip = 1

       do
          if (j==0) then
             x2pix1(i) = k
             y2pix1(i) = 2*k
             exit
          else
             id = MODULO(J,2)
             j  = j/2
             k  = ip*id+k
             ip = ip*4
          endif
       enddo

    enddo

    return
  end subroutine mk_xy2pix1
  !=======================================================================

  !=======================================================================
  subroutine ang2vec(theta, phi, vector)
    !=======================================================================
    !     renders the vector (x,y,z) corresponding to angles
    !     theta (co-latitude measured from North pole, in [0,Pi] radians)
    !     and phi (longitude measured eastward, in radians)
    !     North pole is (x,y,z)=(0,0,1)
    !     added by EH, Feb 2000
    !=======================================================================
    REAL(KIND=DP), INTENT(IN) :: theta, phi
    REAL(KIND=DP), INTENT(OUT), dimension(1:) :: vector

    REAL(KIND=DP) :: sintheta
    !=======================================================================

    if (theta<0.0_dp .or. theta>pi)  then
       print*,"ANG2VEC: theta : ",theta," is out of range [0, Pi]"
       call fatal_error
    endif
    sintheta = SIN(theta)

    vector(1) = sintheta * COS(phi)
    vector(2) = sintheta * SIN(phi)
    vector(3) = COS(theta)

    return
  end subroutine ang2vec
  !=======================================================================
  subroutine vec2ang(vector, theta, phi)
    !=======================================================================
    !     renders the angles theta, phi corresponding to vector (x,y,z)
    !     theta (co-latitude measured from North pole, in [0,Pi] radians)
    !     and phi (longitude measured eastward, in [0,2Pi[ radians)
    !     North pole is (x,y,z)=(0,0,1)
    !     added by EH, Feb 2000
    ! 2011-08: replaced ACOS(z) by more accurate ATAN2(r,z)
    !=======================================================================
    REAL(KIND=DP), INTENT(IN), dimension(1:) :: vector
    REAL(KIND=DP), INTENT(OUT) :: theta, phi

    !=======================================================================

    theta = atan2(sqrt(vector(1)**2 + vector(2)**2), vector(3))

    phi = 0.0_dp
    if (vector(1) /= 0.0_dp .or. vector(2) /= 0.0_dp) &
         &     phi = ATAN2(vector(2),vector(1)) ! phi in ]-pi,pi]
    if (phi < 0.0)     phi = phi + twopi ! phi in [0,2pi[

    return
  end subroutine vec2ang
  !=======================================================================
  function nside2npix(nside) result(npix_result)
    !=======================================================================
    ! given nside, returns npix such that npix = 12*nside^2
    !  nside should be a power of 2 smaller than ns_max
    !  if not, -1 is returned
    ! EH, Feb-2000
    ! 2009-03-04: returns i8b result, faster
    !=======================================================================
    INTEGER(KIND=I8B)             :: npix_result
    INTEGER(KIND=I4B), INTENT(IN) :: nside

    INTEGER(KIND=I8B) :: npix
    CHARACTER(LEN=*), PARAMETER :: code = "nside2npix"
    !=======================================================================

    npix = (12_i8b*nside)*nside
    if (nside < 1 .or. nside > ns_max .or. iand(nside-1,nside) /= 0) then
       print*,code//": Nside="//trim(string(nside))//" is not a power of 2."
       npix = -1
    endif
    npix_result = npix

    return
  end function nside2npix
  !=======================================================================
  subroutine surface_triangle(vec1, vec2, vec3, surface)
    !=======================================================================
    ! returns the surface in steradians
    !  of the spherical triangle with vertices vec1, vec2, vec3
    !
    ! algorithm : finds triangle sides and uses l'Huilier formula to compute
    ! "spherical excess" = surface area of triangle on a sphere of radius one
    ! see, eg Bronshtein, Semendyayev Eq 2.86
    !=======================================================================
    real(kind=dp), dimension(1:), intent(in) :: vec1, vec2, vec3
    real(kind=dp), intent(out) :: surface

    !   real(kind=dp), dimension(1:3) :: v1, v2, v3
    real(kind=dp), dimension(1:3) :: side
    !  real(kind=dp) :: hp
    real(kind=dp) :: x0, x1, x2, x3
    !=======================================================================

    ! half perimeter
!   hp = 0.5_dp * (side1 + side2 + side3)

!   ! l'Huilier formula
!   x0 = tan( hp          * 0.5_dp)
!   x1 = tan((hp - side1) * 0.5_dp)
!   x2 = tan((hp - side2) * 0.5_dp)
!   x3 = tan((hp - side3) * 0.5_dp)

    ! find triangle sides
    call angdist(vec2, vec3, side(1))
    call angdist(vec3, vec1, side(2))
    call angdist(vec1, vec2, side(3))
    ! divide by 4
    side(1:3) = side(1:3) * 0.25_dp

    ! l'Huilier formula
    x0 = tan( side(1) + side(2) + side(3) )
    x1 = tan(-side(1) + side(2) + side(3) )
    x2 = tan( side(1) - side(2) + side(3) )
    x3 = tan( side(1) + side(2) - side(3) )
    surface = 4.0_dp * atan( sqrt(x0 * x1 * x2 * x3) )

    return
  end subroutine surface_triangle
  !=======================================================================
  subroutine angdist(v1, v2, dist)
    !=======================================================================
    ! call angdist(v1, v2, dist)
    ! computes the angular distance dist (in rad) between 2 vectors v1 and v2
    ! dist = atan2 ( |v1 x v2| / (v1 . v2) )
    ! (more accurate than acos(v1. v2) when dist close to 0 or Pi.
    ! 2011-08-25: replaced ACOS with ATAN2
    !=======================================================================
    real(kind=DP), intent(IN), dimension(1:) :: v1, v2
    real(kind=DP), intent(OUT) :: dist

    real(kind=DP), dimension(1:3) :: v3
    real(kind=DP) :: sprod, vprod
    !=======================================================================

    ! scalar product     s = A. cos(theta)
    sprod = dot_product(v1, v2)
    ! vectorial product  v = A. sin(theta)
    call vect_prod(v1, v2, v3)
    vprod = sqrt(dot_product(v3,v3))
    ! theta = atan( |v|/s) in [0,Pi]
    dist = atan2( vprod, sprod)

    return
  end subroutine angdist
  !=======================================================================
  subroutine vect_prod(v1, v2, v3) !
    !=======================================================================
    !     returns in v3 the vectorial product of the 2 3-dimensional
    !     vectors v1 and v2
    !=======================================================================
    real(kind=DP), dimension(1:), INTENT(IN)  :: v1, v2
    real(kind=DP), dimension(1:), INTENT(OUT) :: v3
    !=======================================================================

    v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
    v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
    v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

    return
  end subroutine vect_prod

  !****************************************************************************

  !=======================================================================
  function nside2ntemplates(nside) result(ntemplates)
    !=======================================================================
    ! returns the number of template pixels
    ! 2009-03-24: accepts Nside > 8192, and returns I8B result
    !

    integer(kind=i8b) :: ntemplates
    integer(kind=i4b), intent(IN) :: nside

    ntemplates = 1_I8B + nside * (nside + 6_I8B)
    ntemplates = ntemplates / 4_I8B

    return
  end function nside2ntemplates

!******************************************************************************
  function nside2npweights(nside) result(npweights)
!---------------------------------------
!  returns the number of full sky weights (in compressed form) for a given Nside.
! Given the symmetries of the Healpix lay out, the number of non-redundant weights
! is  nf = ((3*Nside+1)*(Nside+1))/4 ~ Npix/16
!
! adapted from nside2npweights.pro
! 2018-05-18
!---------------------------------------
  
    integer(i8b) :: npweights
    integer(i4b), intent(in) :: nside
    integer(i8b), parameter  :: one  = 1_i8b
    integer(i8b), parameter  :: four = 4_i8b

    npweights = ((3*Nside + one)*(Nside + one)) / four

    return
  end function nside2npweights


end module pix_tools



